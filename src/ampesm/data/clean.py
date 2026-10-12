"""Limpieza: filtros, deduplicado y emparejado de longitudes.

Cada paso cuenta cuanto descarta y ese conteo va al reporte. La razon es la
leccion de P-01: un conteo que no cuadra con lo esperado es la senal mas barata
de que algo se rompio en silencio.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..config import (
    AMINO_ACIDS,
    CLEAN_REPORT_PATH,
    PEPTIDES_PATH,
    RAW_DIR,
    ensure_dirs,
    load_params,
)
from .parse import PeptideRecord, parse_source
from .sources import get_sources

STANDARD = set(AMINO_ACIDS)

# Prioridad al deduplicar: ante la misma secuencia en varias fuentes, se
# conserva el registro de la fuente mas curada, porque trae mejores metadatos.
SOURCE_PRIORITY = {"apd3": 0, "dramp": 1, "uniprot_amp": 2, "uniprot_neg": 3}


def seq_sha1(sequence: str) -> str:
    """Hash de la secuencia: clave de deduplicado y del cache de embeddings.

    Se usa la secuencia y no el id porque el mismo peptido aparece con ids
    distintos en fuentes distintas, y queremos calcular su embedding una vez.
    """
    return hashlib.sha1(sequence.encode("utf-8")).hexdigest()


@dataclass
class CleanStats:
    """Contadores de un paso de filtrado, por fuente."""

    counts: dict[str, Counter]

    @classmethod
    def empty(cls) -> CleanStats:
        return cls(counts=defaultdict(Counter))

    def add(self, source: str, key: str, n: int = 1) -> None:
        self.counts[source][key] += n

    def to_dict(self) -> dict:
        return {src: dict(c) for src, c in self.counts.items()}


# --- Filtros ---------------------------------------------------------------


def is_standard(sequence: str) -> bool:
    """True si la secuencia solo contiene los 20 aminoacidos estandar.

    Se descartan las que tienen X (desconocido), B/J/Z (ambiguos), U/O (no
    estandar), '*' (parada) o huecos de alineamiento. Motivo: ESM-2 los mapea a
    un token desconocido y los descriptores fisicoquimicos no tienen valor
    definido para ellos, asi que incluirlos mete ruido en ambas ramas.
    """
    return bool(sequence) and set(sequence) <= STANDARD


def clean_records(
    records: list[PeptideRecord], min_length: int, max_length: int
) -> tuple[list[PeptideRecord], CleanStats]:
    """Aplica los filtros de alfabeto y longitud."""
    stats = CleanStats.empty()
    kept: list[PeptideRecord] = []

    for rec in records:
        src = rec.source
        stats.add(src, "input")
        seq = rec.sequence.strip().upper().replace(" ", "")

        if not seq:
            stats.add(src, "dropped_empty")
            continue
        if not is_standard(seq):
            stats.add(src, "dropped_non_standard")
            continue
        if len(seq) < min_length:
            stats.add(src, "dropped_too_short")
            continue
        if len(seq) > max_length:
            stats.add(src, "dropped_too_long")
            continue

        rec.sequence = seq
        kept.append(rec)
        stats.add(src, "kept")

    return kept, stats


# --- Deduplicado -----------------------------------------------------------


def deduplicate(records: list[PeptideRecord]) -> tuple[list[PeptideRecord], dict]:
    """Deduplica por secuencia exacta y resuelve conflictos de etiqueta.

    Dos situaciones distintas:

    1. La misma secuencia en varias fuentes de la misma clase. Normal: hay mucho
       solapamiento entre APD3 y DRAMP. Se conserva la de mayor prioridad.

    2. La misma secuencia como positiva y como negativa. Esto es ruido de
       etiquetado: si esta en APD3 tiene actividad verificada, asi que la
       anotacion ausente en UniProt es el error. Se conserva la positiva y se
       cuenta el caso, porque el numero es interesante por si mismo: mide
       directamente la tasa de falsos negativos del conjunto de negativos.
    """
    by_hash: dict[str, list[PeptideRecord]] = defaultdict(list)
    for rec in records:
        by_hash[seq_sha1(rec.sequence)].append(rec)

    kept: list[PeptideRecord] = []
    label_conflicts = 0
    duplicate_pairs = 0
    conflict_examples: list[str] = []

    for group in by_hash.values():
        duplicate_pairs += len(group) - 1
        labels = {r.label for r in group}

        if len(labels) > 1:
            label_conflicts += 1
            if len(conflict_examples) < 10:
                conflict_examples.append(group[0].sequence)
            group = [r for r in group if r.label == 1]

        group.sort(key=lambda r: SOURCE_PRIORITY.get(r.source, 99))
        winner = group[0]
        winner.extra["n_sources"] = len({r.source for r in group})
        winner.extra["all_sources"] = ",".join(sorted({r.source for r in group}))
        kept.append(winner)

    report = {
        "unique_sequences": len(kept),
        "duplicate_records_removed": duplicate_pairs,
        "label_conflicts": label_conflicts,
        "label_conflict_note": (
            "Secuencias presentes como positivas y como negativas. Se conservo la "
            "positiva. Este numero es una estimacion por lo bajo de la tasa de "
            "falsos negativos en el conjunto de negativos (ver docs/01-objetivos.md)."
        ),
        "label_conflict_examples": conflict_examples,
    }
    return kept, report


# --- Emparejado de longitudes ----------------------------------------------


def match_lengths(df: pd.DataFrame, ratio: float | None, seed: int) -> tuple[pd.DataFrame, dict]:
    """Submuestrea negativos para imitar la distribucion de longitud positiva.

    Sin esto, el modelo aprende la longitud como atajo: si los negativos de
    UniProt son sistematicamente mas largos o mas cortos, basta contar residuos
    para separar las clases, y el ROC-AUC sale alto sin haber aprendido nada de
    biologia.

    Se trabaja longitud a longitud: para cada valor L, se toman hasta
    `ratio * n_positivos(L)` negativos de esa misma longitud. Si no hay
    suficientes negativos de longitud L, se toman todos los que haya y el
    deficit se registra, pero no se compensa con otras longitudes: hacerlo
    volveria a desbalancear la distribucion, que es justo lo que se intenta
    evitar.
    """
    rng = np.random.default_rng(seed)
    pos = df[df.label == 1]
    neg = df[df.label == 0]

    pos_by_len = pos.length.value_counts().to_dict()
    chosen_idx: list = []
    deficit = 0

    for length, neg_group in neg.groupby("length"):
        n_pos = pos_by_len.get(length, 0)
        if n_pos == 0:
            continue  # ningun positivo de esta longitud: no se necesitan negativos
        want = len(neg_group) if ratio is None else int(round(ratio * n_pos))
        take = min(want, len(neg_group))
        deficit += max(0, want - len(neg_group))
        picked = rng.choice(neg_group.index.to_numpy(), size=take, replace=False)
        chosen_idx.extend(picked.tolist())

    out = pd.concat([pos, neg.loc[chosen_idx]]).sort_values("id").reset_index(drop=True)

    report = {
        "applied": True,
        "ratio_requested": ratio,
        "negatives_before": int(len(neg)),
        "negatives_after": int(len(chosen_idx)),
        "positives": int(len(pos)),
        "deficit": int(deficit),
        "deficit_note": (
            "Negativos que se habrian querido para igualar la proporcion en alguna "
            "longitud y no existian en la fuente."
        ),
    }
    return out, report


# --- Orquestacion ----------------------------------------------------------


def build_dataset(keys: list[str] | None = None) -> pd.DataFrame:
    """Lee los crudos, limpia, deduplica, empareja longitudes y escribe el CSV."""
    ensure_dirs()
    params = load_params()
    dcfg = params.data
    min_len = int(dcfg.get("min_length", 5))
    max_len = int(dcfg.get("max_length", 50))

    sources = get_sources(keys)
    all_records: list[PeptideRecord] = []
    parse_reports = []

    for src in sources:
        path = RAW_DIR / src.filename
        if not path.exists():
            raise FileNotFoundError(f"Falta {path}. Corre 'make data' o 'ampesm download' primero.")
        recs, rep = parse_source(path, src)
        print(f"[limpieza] {src.key}: {rep['parsed']} registros parseados")
        if src.key == "apd3" and rep.get("dropped_empty") != 1:
            # P-01: se espera exactamente una cabecera sin secuencia.
            print(
                f"  AVISO: se esperaba 1 cabecera vacia en APD3 y hubo "
                f"{rep.get('dropped_empty')}. El formato pudo cambiar; revisa "
                f"docs/09-problemas-encontrados.md (P-01)."
            )
        all_records.extend(recs)
        parse_reports.append(rep)

    kept, filter_stats = clean_records(all_records, min_len, max_len)
    print(f"[limpieza] tras filtros de alfabeto y longitud: {len(kept)}")

    deduped, dedup_report = deduplicate(kept)
    print(
        f"[limpieza] tras deduplicado: {len(deduped)} "
        f"({dedup_report['duplicate_records_removed']} duplicados, "
        f"{dedup_report['label_conflicts']} conflictos de etiqueta)"
    )

    df = pd.DataFrame(
        {
            "id": [f"{r.source}:{r.id}" for r in deduped],
            "sequence": [r.sequence for r in deduped],
            "length": [len(r.sequence) for r in deduped],
            "label": [r.label for r in deduped],
            "source": [r.source for r in deduped],
            "seq_sha1": [seq_sha1(r.sequence) for r in deduped],
            "n_sources": [r.extra.get("n_sources", 1) for r in deduped],
            "all_sources": [r.extra.get("all_sources", r.source) for r in deduped],
        }
    )

    # P-10b: la deduplicacion es por secuencia, no por id, asi que dos
    # secuencias distintas bajo el mismo accession dejan ids repetidos. Mas
    # abajo en la tuberia eso pisaria una de las dos en silencio, asi que se
    # desambigua aqui con los primeros caracteres del hash de la secuencia.
    dup_mask = df.id.duplicated(keep=False)
    n_dups = int(dup_mask.sum())
    if n_dups:
        df.loc[dup_mask, "id"] = df.loc[dup_mask, "id"] + "#" + df.loc[dup_mask, "seq_sha1"].str[:8]
        print(f"[limpieza] {n_dups} ids duplicados desambiguados con el hash de la secuencia")
    duplicate_ids_fixed = n_dups

    if dcfg.get("match_negative_lengths", True):
        ratio = dcfg.get("negative_ratio")
        df, match_report = match_lengths(df, ratio, params.seed)
        print(
            f"[limpieza] emparejado de longitudes: "
            f"{match_report['negatives_before']} -> {match_report['negatives_after']} negativos"
        )
    else:
        match_report = {"applied": False}

    # P-10c: un conjunto vacio o de una sola clase llega a np.vstack([]) dos
    # etapas mas tarde y falla con "need at least one array to concatenate",
    # que no dice nada sobre la causa real.
    if df.empty:
        raise RuntimeError(
            "El conjunto quedo vacio tras los filtros. Revisa que las fuentes se "
            "descargaran bien y los limites de longitud en params.yaml."
        )
    if df.label.nunique() < 2:
        solo = "positivos" if df.label.iloc[0] == 1 else "negativos"
        raise RuntimeError(
            f"El conjunto solo tiene {solo}. Hacen falta ambas clases: revisa que "
            f"se incluyan fuentes positivas (apd3, dramp, uniprot_amp) y negativas "
            f"(uniprot_neg). Fuentes usadas: {sorted(df.source.unique())}."
        )

    df.to_csv(PEPTIDES_PATH, index=False)

    report = {
        "min_length": min_len,
        "max_length": max_len,
        "parse": parse_reports,
        "filters": filter_stats.to_dict(),
        "deduplicate": dedup_report,
        "duplicate_ids_disambiguated": duplicate_ids_fixed,
        "length_matching": match_report,
        "final": {
            "n_total": int(len(df)),
            "n_positive": int((df.label == 1).sum()),
            "n_negative": int((df.label == 0).sum()),
            "by_source": {k: int(v) for k, v in df.source.value_counts().items()},
            "length_mean_positive": float(df[df.label == 1].length.mean()),
            "length_mean_negative": float(df[df.label == 0].length.mean()),
        },
    }
    CLEAN_REPORT_PATH.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    f = report["final"]
    print(f"[limpieza] final: {f['n_total']} peptidos ({f['n_positive']}+ / {f['n_negative']}-)")
    print(
        f"[limpieza] longitud media: positivos {f['length_mean_positive']:.1f}, "
        f"negativos {f['length_mean_negative']:.1f}"
    )
    print(f"[limpieza] {PEPTIDES_PATH}")
    return df
