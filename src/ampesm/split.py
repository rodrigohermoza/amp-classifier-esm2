"""Clustering por identidad de secuencia y division train/val/test.

Es la etapa que decide si los numeros del proyecto significan algo. La
justificacion completa esta en docs/04-division-clusters.md; aqui solo el como.

Invariante que debe cumplirse siempre: ningun cluster aparece en mas de un
split. Se verifica en tests/test_split.py. Si ese test falla, todas las
metricas del proyecto son invalidas.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from .config import PEPTIDES_PATH, SPLIT_REPORT_PATH, SPLITS_PATH, ensure_dirs, load_params

# --- Utilidades ------------------------------------------------------------


def write_fasta(sequences: dict[str, str], path: Path) -> None:
    """Escribe un FASTA con los ids como cabecera."""
    with open(path, "w", encoding="utf-8") as fh:
        for seq_id, seq in sequences.items():
            fh.write(f">{seq_id}\n{seq}\n")


def tool_version(cmd: list[str]) -> str | None:
    """Version de una herramienta externa, o None si no esta disponible."""
    if shutil.which(cmd[0]) is None:
        return None
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        text = (out.stdout or out.stderr).strip()
        return text.splitlines()[0] if text else "desconocida"
    except (subprocess.SubprocessError, OSError):
        return None


# --- Metodo 1: MMseqs2 -----------------------------------------------------


def cluster_mmseqs2(
    sequences: dict[str, str], identity: float, coverage: float
) -> dict[str, int] | None:
    """Agrupa con MMseqs2. Devuelve id -> cluster, o None si no esta instalado.

    --cov-mode 1 mide la cobertura respecto a la secuencia mas corta, que es lo
    que interesa al mezclar peptidos de 5 y de 50 residuos: con el modo por
    defecto, un peptido corto nunca alcanzaria la cobertura exigida frente a uno
    largo y quedaria en su propio cluster, perdiendo la agrupacion.
    """
    if shutil.which("mmseqs") is None:
        return None

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        fasta = tmpdir / "input.fasta"
        prefix = tmpdir / "clust"
        write_fasta(sequences, fasta)

        cmd = [
            "mmseqs", "easy-cluster", str(fasta), str(prefix), str(tmpdir / "tmp"),
            "--min-seq-id", str(identity),
            "-c", str(coverage),
            "--cov-mode", "1",
            "-v", "1",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            stderr = getattr(exc, "stderr", "") or ""
            print(f"[split] MMseqs2 fallo: {stderr[:500]}")
            return None

        tsv = prefix.parent / f"{prefix.name}_cluster.tsv"
        if not tsv.exists():
            print("[split] MMseqs2 no produjo el archivo de clusteres esperado")
            return None

        # Formato: representante<TAB>miembro, una linea por miembro.
        rep_to_id: dict[str, int] = {}
        mapping: dict[str, int] = {}
        with open(tsv, encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 2:
                    continue
                rep, member = parts
                if rep not in rep_to_id:
                    rep_to_id[rep] = len(rep_to_id)
                mapping[member] = rep_to_id[rep]

    missing = set(sequences) - set(mapping)
    if missing:
        # No deberia pasar, pero si pasa es mejor darles cluster propio que
        # perderlos en silencio.
        next_id = max(mapping.values(), default=-1) + 1
        for i, seq_id in enumerate(sorted(missing)):
            mapping[seq_id] = next_id + i
        print(f"[split] AVISO: {len(missing)} secuencias sin cluster; se les asigno uno propio")

    return mapping


# --- Metodo 2: CD-HIT ------------------------------------------------------


def cluster_cdhit(
    sequences: dict[str, str], identity: float, coverage: float
) -> dict[str, int] | None:
    """Agrupa con CD-HIT. Devuelve id -> cluster, o None si no esta instalado.

    Dos ajustes obligados para peptidos cortos: el tamano de palabra tiene que
    bajar con el umbral de identidad (-n 2 para 0.4, que es el minimo que acepta
    CD-HIT), y -l 4 permite secuencias de 5 residuos, que el valor por defecto
    de 10 descartaria.
    """
    if shutil.which("cd-hit") is None:
        return None

    word_size = 5 if identity >= 0.7 else 4 if identity >= 0.6 else 3 if identity >= 0.5 else 2

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        fasta = tmpdir / "input.fasta"
        out = tmpdir / "out"
        write_fasta(sequences, fasta)

        cmd = [
            "cd-hit", "-i", str(fasta), "-o", str(out),
            "-c", str(identity), "-n", str(word_size),
            "-l", "4", "-aS", str(coverage), "-d", "0", "-M", "2000",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            stderr = getattr(exc, "stderr", "") or ""
            print(f"[split] CD-HIT fallo: {stderr[:500]}")
            return None

        clstr = Path(str(out) + ".clstr")
        if not clstr.exists():
            return None

        mapping: dict[str, int] = {}
        current = -1
        with open(clstr, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">Cluster"):
                    current = int(line.split()[1])
                elif ">" in line:
                    # '0  33aa, >AP00001... *'
                    name = line.split(">", 1)[1].split("...")[0].strip()
                    mapping[name] = current

    return mapping if mapping else None


# --- Metodo 3: respaldo propio --------------------------------------------


def kmer_set(sequence: str, k: int) -> frozenset[str]:
    """Conjunto de k-meros. Si la secuencia es mas corta que k, ella misma."""
    if len(sequence) < k:
        return frozenset({sequence})
    return frozenset(sequence[i : i + k] for i in range(len(sequence) - k + 1))


def _make_aligner():
    """Alineador global con parametros tipo BLASTP.

    Se importa Biopython aqui y no arriba para que el modulo se pueda importar
    sin Bio instalado cuando se usa MMseqs2.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def identity_and_coverage(seq_a: str, seq_b: str, aligner) -> tuple[float, float]:
    """Identidad y cobertura de dos secuencias, por alineamiento global.

    Identidad = residuos identicos / longitud de la mas corta.
    Cobertura  = columnas alineadas / longitud de la mas corta.

    Se normaliza por la mas corta y no por la longitud del alineamiento porque
    con peptidos de longitudes muy distintas lo segundo castigaria a un peptido
    corto que esta contenido entero en uno largo, que es precisamente el caso
    que queremos detectar como relacionado.
    """
    try:
        alignment = aligner.align(seq_a, seq_b)[0]
    except (ValueError, IndexError):
        return 0.0, 0.0

    shorter = min(len(seq_a), len(seq_b))
    matches = 0
    aligned_cols = 0
    for (a_start, a_end), (b_start, b_end) in zip(*alignment.aligned):
        block = a_end - a_start
        aligned_cols += block
        matches += sum(
            1 for i in range(block) if seq_a[a_start + i] == seq_b[b_start + i]
        )

    return matches / shorter, aligned_cols / shorter


def cluster_builtin(
    sequences: dict[str, str],
    identity: float,
    coverage: float,
    kmer_size: int = 3,
    jaccard_prefilter: float = 0.15,
) -> dict[str, int]:
    """Clustering voraz al estilo CD-HIT, sin herramientas externas.

    1. Ordenar de mas larga a mas corta (las largas son mejores representantes).
    2. Para cada secuencia, alinear solo contra los representantes que comparten
       k-meros por encima del umbral de Jaccard, usando un indice invertido de
       k-meros para no recorrer todos los representantes.
    3. Si alguno pasa el umbral de identidad y cobertura, unirse a ese cluster;
       si no, volverse representante de un cluster nuevo.

    El prefiltro es lo que hace esto viable: alinear cada secuencia contra miles
    de representantes en Python puro seria inviable, pero las operaciones de
    conjuntos descartan la gran mayoria de pares antes de alinear.
    """
    aligner = _make_aligner()
    order = sorted(sequences, key=lambda i: (-len(sequences[i]), i))

    kmer_index: dict[str, set[int]] = defaultdict(set)  # k-mero -> indices de repr.
    rep_ids: list[str] = []
    rep_kmers: list[frozenset[str]] = []
    mapping: dict[str, int] = {}
    alignments_done = 0

    for seq_id in order:
        seq = sequences[seq_id]
        kmers = kmer_set(seq, kmer_size)

        # Candidatos: representantes que comparten al menos un k-mero.
        shared: dict[int, int] = defaultdict(int)
        for km in kmers:
            for rep_idx in kmer_index.get(km, ()):
                shared[rep_idx] += 1

        # Se prueban primero los de mayor solapamiento: si uno cercano pasa el
        # umbral, se corta y se ahorran los alineamientos restantes.
        candidates = sorted(shared.items(), key=lambda kv: -kv[1])
        assigned: int | None = None

        for rep_idx, n_shared in candidates:
            union = len(kmers | rep_kmers[rep_idx])
            if union == 0 or (n_shared / union) < jaccard_prefilter:
                continue
            ident, cov = identity_and_coverage(sequences[rep_ids[rep_idx]], seq, aligner)
            alignments_done += 1
            if ident >= identity and cov >= coverage:
                assigned = rep_idx
                break

        if assigned is None:
            assigned = len(rep_ids)
            rep_ids.append(seq_id)
            rep_kmers.append(kmers)
            for km in kmers:
                kmer_index[km].add(assigned)

        mapping[seq_id] = assigned

    print(
        f"[split] respaldo propio: {len(rep_ids)} clusteres, "
        f"{alignments_done:,} alineamientos realizados"
    )
    return mapping


# --- Asignacion de clusteres a splits -------------------------------------


def assign_splits(
    df: pd.DataFrame, fractions: dict[str, float], seed: int
) -> dict[int, str]:
    """Reparte clusteres completos entre train/val/test.

    Heuristica voraz: se recorren los clusteres de mayor a menor tamano y cada
    uno va al split que esta relativamente mas lejos de cumplir su cuota; a
    igualdad, al que degrada menos el balance de clases.

    Los grandes van primero a proposito. Si se dejan para el final, un unico
    cluster grande puede desbordar la cuota de test y descuadrar las
    proporciones por completo.

    El reparto exacto es un problema de particion (NP-duro) y no hace falta
    resolverlo: lo que importa es que no haya fuga, no que los porcentajes
    salgan exactos. Los porcentajes reales se registran en el reporte.
    """
    rng = np.random.default_rng(seed)
    total = len(df)
    targets = {k: v * total for k, v in fractions.items()}

    sizes = df.groupby("cluster_id").size()
    pos_counts = df[df.label == 1].groupby("cluster_id").size()

    clusters = [
        {
            "id": int(cid),
            "n": int(sizes[cid]),
            "n_pos": int(pos_counts.get(cid, 0)),
        }
        for cid in sizes.index
    ]
    # El desempate aleatorio evita que el orden dependa del id del cluster, que
    # a su vez depende del orden del archivo de entrada.
    rng.shuffle(clusters)
    clusters.sort(key=lambda c: -c["n"])

    assigned: dict[int, str] = {}
    current = {k: 0 for k in fractions}
    current_pos = {k: 0 for k in fractions}
    global_pos_rate = float(df.label.mean())

    for cluster in clusters:
        best_split, best_score = None, None
        for split in fractions:
            # Deficit relativo: que fraccion de su cuota le falta a este split.
            deficit = (targets[split] - current[split]) / max(targets[split], 1.0)

            # Desviacion del balance de clases si se le anade este cluster.
            n_after = current[split] + cluster["n"]
            pos_after = current_pos[split] + cluster["n_pos"]
            balance_penalty = abs((pos_after / n_after) - global_pos_rate)

            score = deficit - 0.5 * balance_penalty
            if best_score is None or score > best_score:
                best_split, best_score = split, score

        assigned[cluster["id"]] = best_split
        current[best_split] += cluster["n"]
        current_pos[best_split] += cluster["n_pos"]

    return assigned


# --- Orquestacion ---------------------------------------------------------


def build_splits(force_method: str | None = None) -> pd.DataFrame:
    """Agrupa y divide. Escribe splits.csv y split_report.json."""
    ensure_dirs()
    params = load_params()
    scfg = params.split

    identity = float(scfg.get("identity_threshold", 0.4))
    coverage = float(scfg.get("min_coverage", 0.5))
    fractions = dict(scfg.get("fractions", {"train": 0.70, "val": 0.15, "test": 0.15}))
    methods = list(scfg.get("methods", ["mmseqs2", "cdhit", "builtin"]))
    bcfg = dict(scfg.get("builtin", {}))

    if not PEPTIDES_PATH.exists():
        raise FileNotFoundError(f"Falta {PEPTIDES_PATH}. Corre 'make data' primero.")
    df = pd.read_csv(PEPTIDES_PATH)
    sequences = dict(zip(df.id, df.sequence))
    print(f"[split] {len(sequences)} secuencias, umbral de identidad {identity}")

    if force_method:
        methods = [force_method]

    mapping: dict[str, int] | None = None
    used_method = None
    used_version = None

    for method in methods:
        if method == "mmseqs2":
            version = tool_version(["mmseqs", "version"])
            if version is None:
                print("[split] mmseqs no disponible")
                continue
            print(f"[split] usando MMseqs2 ({version})")
            mapping = cluster_mmseqs2(sequences, identity, coverage)
        elif method == "cdhit":
            version = tool_version(["cd-hit", "-h"])
            if version is None:
                print("[split] cd-hit no disponible")
                continue
            print(f"[split] usando CD-HIT ({version})")
            mapping = cluster_cdhit(sequences, identity, coverage)
        elif method == "builtin":
            version = "interno"
            print("[split] usando el respaldo propio en Python (puede tardar minutos)")
            mapping = cluster_builtin(
                sequences,
                identity,
                coverage,
                kmer_size=int(bcfg.get("kmer_size", 3)),
                jaccard_prefilter=float(bcfg.get("kmer_jaccard_prefilter", 0.15)),
            )
        else:
            raise ValueError(f"Metodo de clustering desconocido: {method!r}")

        if mapping:
            used_method, used_version = method, version
            break

    if not mapping:
        raise RuntimeError(f"Ningun metodo de clustering funciono. Probados: {methods}")

    df["cluster_id"] = df.id.map(mapping)
    if df.cluster_id.isna().any():
        missing = int(df.cluster_id.isna().sum())
        raise RuntimeError(f"{missing} secuencias quedaron sin cluster: hay un bug en el mapeo")
    df["cluster_id"] = df.cluster_id.astype(int)

    assignment = assign_splits(df, fractions, params.seed)
    df["split"] = df.cluster_id.map(assignment)

    # Verificacion del invariante, aqui y no solo en los tests: si esto falla,
    # es mejor no escribir el archivo que escribir uno invalido.
    per_cluster = df.groupby("cluster_id").split.nunique()
    if (per_cluster > 1).any():
        culprits = per_cluster[per_cluster > 1].index.tolist()[:10]
        raise RuntimeError(f"FUGA: clusteres repartidos entre splits: {culprits}")

    df.to_csv(SPLITS_PATH, index=False)

    cluster_sizes = df.groupby("cluster_id").size()
    mixed = df.groupby("cluster_id").label.nunique()
    report = {
        "method": used_method,
        "tool_version": used_version,
        "identity_threshold": identity,
        "min_coverage": coverage,
        "seed": params.seed,
        "n_sequences": int(len(df)),
        "n_clusters": int(df.cluster_id.nunique()),
        "cluster_size": {
            "max": int(cluster_sizes.max()),
            "mean": float(cluster_sizes.mean()),
            "median": float(cluster_sizes.median()),
            "singletons": int((cluster_sizes == 1).sum()),
        },
        "mixed_label_clusters": int((mixed > 1).sum()),
        "mixed_label_note": (
            "Clusteres con secuencias de ambas clases: secuencias parecidas con "
            "etiquetas distintas. Probablemente es donde se concentra el ruido de "
            "etiquetado. Ver docs/08-decisiones.md, D-03."
        ),
        "splits": {},
    }
    for split in fractions:
        sub = df[df.split == split]
        report["splits"][split] = {
            "n": int(len(sub)),
            "fraction_actual": round(len(sub) / len(df), 4),
            "fraction_target": fractions[split],
            "n_positive": int((sub.label == 1).sum()),
            "n_negative": int((sub.label == 0).sum()),
            "positive_rate": round(float(sub.label.mean()), 4) if len(sub) else 0.0,
            "n_clusters": int(sub.cluster_id.nunique()),
        }

    SPLIT_REPORT_PATH.write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(f"[split] {report['n_clusters']} clusteres ({report['cluster_size']['singletons']} de 1)")
    print(f"[split] {report['mixed_label_clusters']} clusteres con ambas clases")
    for split, info in report["splits"].items():
        print(
            f"[split] {split:5s}: {info['n']:6d} ({info['fraction_actual']:.1%}) "
            f"{info['n_positive']}+ / {info['n_negative']}-  "
            f"tasa positiva {info['positive_rate']:.3f}"
        )
    print(f"[split] {SPLITS_PATH}")
    return df
