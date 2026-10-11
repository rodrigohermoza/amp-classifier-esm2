"""Rutas, constantes y carga de parametros.

Este modulo no importa nada del proyecto y no toca la red. Todo lo demas puede
importarlo sin riesgo de ciclos.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

# --- Rutas -----------------------------------------------------------------
# Se resuelven desde la ubicacion de este archivo para que la tuberia funcione
# sin importar desde que directorio se invoque el CLI.
PKG_DIR = Path(__file__).resolve().parent
ROOT = PKG_DIR.parent.parent

DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"

RESULTS_DIR = ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
METRICS_DIR = RESULTS_DIR / "metrics"
MODELS_DIR = RESULTS_DIR / "models"
EMBEDDINGS_DIR = RESULTS_DIR / "embeddings"

PARAMS_PATH = ROOT / "params.yaml"

# Artefactos con nombre fijo: son el contrato entre etapas (docs/02).
MANIFEST_PATH = RAW_DIR / "manifest.json"
PEPTIDES_PATH = PROCESSED_DIR / "peptides.csv"
SPLITS_PATH = PROCESSED_DIR / "splits.csv"
CLEAN_REPORT_PATH = PROCESSED_DIR / "clean_report.json"
SPLIT_REPORT_PATH = PROCESSED_DIR / "split_report.json"

ALL_DIRS = (
    RAW_DIR,
    INTERIM_DIR,
    PROCESSED_DIR,
    FIGURES_DIR,
    METRICS_DIR,
    MODELS_DIR,
    EMBEDDINGS_DIR,
)

# --- Constantes biologicas -------------------------------------------------
# Los 20 aminoacidos estandar. El orden es fijo porque define el orden de las
# columnas del vector de composicion; cambiarlo invalidaria modelos guardados.
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# Caracteres que aparecen en FASTA y que NO son aminoacidos estandar:
#   B, J, O, U, Z -> codigos raros o no estandar
#   X             -> residuo desconocido
#   *             -> codon de parada
#   -, .          -> huecos de alineamiento
NON_STANDARD = set("BJOUXZ*-.")

SEED = 42


@dataclass(frozen=True)
class Params:
    """Contenedor tipado de params.yaml.

    Se guarda el diccionario crudo en `raw` y se exponen los bloques como
    atributos. Deliberadamente no se validan todos los campos: los errores de
    configuracion aparecen al usarse, con un KeyError que dice el nombre exacto.
    """

    raw: dict[str, Any] = field(default_factory=dict)

    @property
    def seed(self) -> int:
        return int(self.raw.get("seed", SEED))

    def section(self, name: str) -> dict[str, Any]:
        return dict(self.raw.get(name, {}))

    @property
    def data(self) -> dict[str, Any]:
        return self.section("data")

    @property
    def split(self) -> dict[str, Any]:
        return self.section("split")

    @property
    def features(self) -> dict[str, Any]:
        return self.section("features")

    @property
    def embeddings(self) -> dict[str, Any]:
        return self.section("embeddings")

    @property
    def models(self) -> dict[str, Any]:
        return self.section("models")

    @property
    def evaluate(self) -> dict[str, Any]:
        return self.section("evaluate")

    @property
    def interpret(self) -> dict[str, Any]:
        return self.section("interpret")


def load_params(path: Path | None = None) -> Params:
    """Lee params.yaml. Si no existe, devuelve parametros vacios.

    Que falte el archivo no es un error fatal: cada consumidor define sus
    valores por defecto, asi los tests pueden correr sin el archivo.
    """
    path = path or PARAMS_PATH
    if not path.exists():
        return Params(raw={})
    with open(path, encoding="utf-8") as fh:
        return Params(raw=yaml.safe_load(fh) or {})


def ensure_dirs() -> None:
    """Crea las carpetas de trabajo. Idempotente."""
    for d in ALL_DIRS:
        d.mkdir(parents=True, exist_ok=True)


def hf_cache_dir() -> Path:
    """Donde se guardan los pesos de ESM-2.

    Respeta HF_HOME si el usuario lo definio; si no, usa una carpeta dentro del
    repositorio (ignorada por git) para que borrar el proyecto borre tambien los
    cientos de MB de pesos descargados.
    """
    env = os.environ.get("HF_HOME")
    if env:
        return Path(env)
    return ROOT / ".cache" / "huggingface"
