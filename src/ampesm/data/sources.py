"""Registro declarativo de las fuentes de datos.

Separado de la descarga a proposito: este archivo es la respuesta a "¿de donde
salieron estos datos?" y se puede leer sin leer codigo de red. Las licencias y
las citas estan aqui porque es donde alguien las va a buscar.

Detalle del esquema en docs/03-datos-fuentes.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from urllib.parse import urlencode

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# Palabra clave "Antimicrobial" de UniProt. Define a la vez los positivos de
# UniProt y lo que se excluye de los negativos.
KW_ANTIMICROBIAL = "KW-0929"


@dataclass(frozen=True)
class Source:
    """Una fuente descargable.

    key:      identificador corto, usado como prefijo de los ids y en la
              columna `source` del CSV.
    label:    1 = positivos (AMP), 0 = negativos.
    parser:   "fasta" | "xlsx".
    license:  resumen en una linea; el detalle esta en docs/03.
    """

    key: str
    url: str
    filename: str
    label: int
    parser: str
    description: str
    license: str
    citation: str
    homepage: str

    @property
    def is_positive(self) -> bool:
        return self.label == 1


def _uniprot_url(query: str, fmt: str = "fasta") -> str:
    return f"{UNIPROT_STREAM}?{urlencode({'query': query, 'format': fmt})}"


def uniprot_count_url(query: str) -> str:
    """URL de busqueda con size=1, solo para leer la cabecera x-total-results.

    Se usa para registrar en el manifiesto cuantos registros declaraba UniProt
    en el momento de la descarga, que es independiente de cuantos acabamos
    conservando tras los filtros.
    """
    return f"{UNIPROT_SEARCH}?{urlencode({'query': query, 'format': 'fasta', 'size': 1})}"


# Consultas de UniProt, como constantes nombradas para que aparezcan tal cual en
# la documentacion y en el manifiesto.
QUERY_UNIPROT_POS = f"(keyword:{KW_ANTIMICROBIAL}) AND (length:[5 TO 50])"

# Negativos: revisados (Swiss-Prot), sin la palabra clave antimicrobiana y sin
# la palabra "antimicrobial" en el campo de funcion. Lo segundo captura entradas
# donde la actividad se describe en texto libre pero nadie asigno la keyword.
QUERY_UNIPROT_NEG = (
    "(length:[5 TO 50]) AND (reviewed:true) "
    f"NOT (keyword:{KW_ANTIMICROBIAL}) NOT (cc_function:antimicrobial)"
)


SOURCES: tuple[Source, ...] = (
    Source(
        key="apd3",
        url="https://aps.unmc.edu/assets/sequences/naturalAMPs_APD2024a.fasta",
        filename="apd3_naturalAMPs_APD2024a.fasta",
        label=1,
        parser="fasta",
        description="AMPs naturales con actividad verificada experimentalmente (APD3).",
        license="Uso academico con citacion; sin licencia abierta declarada. No se redistribuye.",
        citation=(
            "Wang G., Li X., Wang Z. APD3: the antimicrobial peptide database as a tool "
            "for research and education. Nucleic Acids Research, 2016."
        ),
        homepage="https://aps.unmc.edu/",
    ),
    Source(
        key="dramp",
        url=(
            "https://dramp.cpu-bioinfor.org/downloads/download.php"
            "?filename=download_data/DRAMP3.0_new/general_amps.xlsx"
        ),
        filename="dramp_general_amps.xlsx",
        label=1,
        parser="xlsx",
        description="AMPs generales de DRAMP 3.0, con metadatos de organismo diana y actividad.",
        license="Uso academico y no comercial. No se redistribuye.",
        citation=(
            "Shi G. et al. DRAMP 3.0: an enhanced comprehensive data repository of "
            "antimicrobial peptides. Nucleic Acids Research, 2022."
        ),
        homepage="https://dramp.cpu-bioinfor.org/",
    ),
    Source(
        key="uniprot_amp",
        url=_uniprot_url(QUERY_UNIPROT_POS),
        filename="uniprot_antimicrobial.fasta",
        label=1,
        parser="fasta",
        description=(
            "Peptidos de UniProtKB con la palabra clave Antimicrobial (KW-0929). "
            "Sirven ademas como control de sesgo de procedencia."
        ),
        license="CC BY 4.0.",
        citation="The UniProt Consortium. UniProt: the Universal Protein Knowledgebase.",
        homepage="https://www.uniprot.org/",
    ),
    Source(
        key="uniprot_neg",
        url=_uniprot_url(QUERY_UNIPROT_NEG),
        filename="uniprot_negatives.fasta",
        label=0,
        parser="fasta",
        description=(
            "Peptidos revisados de UniProtKB SIN anotacion antimicrobiana. "
            "ATENCION: son no-anotados, no negativos verificados. Ver docs/01-objetivos.md."
        ),
        license="CC BY 4.0.",
        citation="The UniProt Consortium. UniProt: the Universal Protein Knowledgebase.",
        homepage="https://www.uniprot.org/",
    ),
)

# Fuente planeada que no se pudo usar. Se deja registrada en codigo, y no solo
# en la documentacion, para que reintentarla sea cambiar una linea.
# Estado 2026-10-11: 403 o respuestas de 0 bytes en todas las rutas.
# Ver docs/09-problemas-encontrados.md, P-02.
DEFERRED_SOURCES: tuple[Source, ...] = (
    Source(
        key="dbaasp",
        url="https://dbaasp.org/api/v1/peptides",
        filename="dbaasp_peptides.json",
        label=1,
        parser="json",
        description="DBAASP: unica fuente con valores de MIC. Habilitaria la regresion.",
        license="Uso academico con citacion.",
        citation="Pirtskhalava M. et al. DBAASP v3. Nucleic Acids Research, 2021.",
        homepage="https://dbaasp.org/",
    ),
)


def get_sources(keys: list[str] | None = None) -> list[Source]:
    """Devuelve las fuentes activas, opcionalmente filtradas por clave."""
    if keys is None:
        return list(SOURCES)
    known = {s.key: s for s in SOURCES}
    missing = [k for k in keys if k not in known]
    if missing:
        raise KeyError(f"Fuentes desconocidas: {missing}. Disponibles: {sorted(known)}")
    return [known[k] for k in keys]
