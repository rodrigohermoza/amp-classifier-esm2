"""Parseo de archivos crudos a registros uniformes.

Ninguna funcion de aqui filtra por longitud ni por alfabeto: eso es trabajo de
clean.py. Aqui solo se convierte formato en registros, y se cuenta lo que se
descarta por estar malformado.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .sources import Source


@dataclass
class PeptideRecord:
    """Un peptido, antes de limpiar."""

    id: str
    sequence: str
    source: str
    label: int
    extra: dict[str, Any] = field(default_factory=dict)


# --- FASTA -----------------------------------------------------------------


def iter_fasta(path: Path) -> list[tuple[str, str]]:
    """Lector de FASTA minimo: devuelve (cabecera, secuencia).

    Se implementa a mano en vez de usar Biopython por dos razones: no tener que
    cargar Bio para la etapa de datos, y poder tolerar el encabezado espurio de
    APD3 sin pelearse con un parser estricto.

    La cabecera se devuelve sin el '>' y sin recortar nada mas; interpretarla es
    responsabilidad de quien llame.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            elif header is not None:
                chunks.append(line)
            # Una linea de secuencia antes de cualquier cabecera se ignora:
            # es un archivo malformado, no un registro.

    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def _id_from_header(header: str, source_key: str, index: int) -> str:
    """Extrae un identificador estable de una cabecera FASTA.

    APD3: '>AP00001' -> AP00001
    UniProt: '>sp|P12345|NAME_HUMAN Descripcion...' -> P12345
    Cualquier otra cosa: el primer token, o un indice si no hay nada usable.
    """
    header = header.strip()
    if not header:
        return f"{source_key}_{index}"

    # Formato UniProt: db|accession|entry_name
    parts = header.split("|")
    if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
        return parts[1]

    token = header.split()[0]
    return token or f"{source_key}_{index}"


def parse_fasta_source(path: Path, source: Source) -> tuple[list[PeptideRecord], dict]:
    """Parsea un FASTA y devuelve (registros, reporte).

    El reporte incluye `dropped_empty`, que es la senal de P-01: el FASTA de
    APD3 arranca con '>Your search led to 3306 peptides', una cabecera sin
    secuencia. Si ese contador no vale exactamente 1 para apd3, el formato del
    archivo cambio y hay que mirarlo antes de seguir.
    """
    raw = iter_fasta(path)
    records: list[PeptideRecord] = []
    dropped_empty: list[str] = []

    for i, (header, seq) in enumerate(raw):
        if not seq:
            dropped_empty.append(header[:80])
            continue
        records.append(
            PeptideRecord(
                id=_id_from_header(header, source.key, i),
                sequence=seq.upper(),
                source=source.key,
                label=source.label,
                extra={"header": header[:200]},
            )
        )

    report = {
        "source": source.key,
        "file": path.name,
        "fasta_entries": len(raw),
        "parsed": len(records),
        "dropped_empty": len(dropped_empty),
        "dropped_empty_headers": dropped_empty[:5],
    }
    return records, report


# --- XLSX ------------------------------------------------------------------

# La hoja de DRAMP tiene decenas de columnas y sus nombres han cambiado entre
# versiones. En vez de fijar un nombre, se busca por patron y se registra que
# columna se acabo usando.
SEQ_COL_PATTERNS = (r"^sequence$", r"sequence", r"^seq$")
ID_COL_PATTERNS = (r"^dramp[_ ]?id$", r"^id$", r"dramp", r"accession")

# Metadatos que vale la pena conservar para el analisis de errores (docs/07).
EXTRA_COL_PATTERNS = {
    "target_organism": (r"target.*organism", r"^target$"),
    "activity": (r"^activity$", r"antimicrobial.*activity"),
    "source_organism": (r"source.*organism", r"^source$"),
    "name": (r"^peptide[_ ]?name$", r"^name$"),
}


def _find_column(columns: list[str], patterns: tuple[str, ...]) -> str | None:
    """Primera columna que coincide, probando los patrones en orden.

    El orden importa: se prueba la coincidencia exacta antes que la parcial,
    para no elegir 'sequence_length' cuando existe 'sequence'.
    """
    lowered = {c: str(c).strip().lower() for c in columns}
    for pat in patterns:
        rx = re.compile(pat)
        for col, low in lowered.items():
            if rx.search(low):
                return col
    return None


def parse_xlsx_source(path: Path, source: Source) -> tuple[list[PeptideRecord], dict]:
    """Parsea una hoja de calculo (DRAMP) y devuelve (registros, reporte)."""
    import pandas as pd

    df = pd.read_excel(path, engine="openpyxl")
    columns = list(df.columns)

    seq_col = _find_column(columns, SEQ_COL_PATTERNS)
    if seq_col is None:
        raise ValueError(
            f"No se encontro columna de secuencia en {path.name}. "
            f"Columnas disponibles: {columns[:30]}"
        )
    id_col = _find_column(columns, ID_COL_PATTERNS)
    extra_cols = {
        name: col
        for name, pats in EXTRA_COL_PATTERNS.items()
        if (col := _find_column(columns, pats)) is not None and col != seq_col
    }

    records: list[PeptideRecord] = []
    dropped_empty = 0

    for i, row in df.iterrows():
        seq = row[seq_col]
        if not isinstance(seq, str) or not seq.strip():
            dropped_empty += 1
            continue
        rec_id = str(row[id_col]).strip() if id_col is not None else f"{source.key}_{i}"
        extra = {
            name: (str(row[col]) if row[col] is not None else "")
            for name, col in extra_cols.items()
        }
        records.append(
            PeptideRecord(
                id=rec_id or f"{source.key}_{i}",
                sequence=seq.strip().upper(),
                source=source.key,
                label=source.label,
                extra=extra,
            )
        )

    report = {
        "source": source.key,
        "file": path.name,
        "rows": int(len(df)),
        "parsed": len(records),
        "dropped_empty": dropped_empty,
        "sequence_column": str(seq_col),
        "id_column": str(id_col) if id_col else None,
        "extra_columns": {k: str(v) for k, v in extra_cols.items()},
        "n_columns": len(columns),
    }
    return records, report


# --- Despachador -----------------------------------------------------------

PARSERS = {"fasta": parse_fasta_source, "xlsx": parse_xlsx_source}


def parse_source(path: Path, source: Source) -> tuple[list[PeptideRecord], dict]:
    if source.parser not in PARSERS:
        raise ValueError(f"Parser no soportado: {source.parser!r} (fuente {source.key})")
    return PARSERS[source.parser](path, source)
