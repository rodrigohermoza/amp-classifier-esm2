"""Composicion de aminoacidos: 20 features.

El baseline mas simple posible y aun asi competitivo, porque la frecuencia de
K y R codifica la carga, y la de L/I/V/F/W la hidrofobicidad, que son las dos
propiedades que gobiernan el mecanismo de disrupcion de membrana.

No captura nada del orden de los residuos. Esa es exactamente la comparacion
que interesa frente a ESM-2.
"""

from __future__ import annotations

import numpy as np

from ..config import AA_INDEX, AMINO_ACIDS

FEATURE_NAMES: list[str] = [f"aa_{aa}" for aa in AMINO_ACIDS]


def composition_vector(sequence: str) -> np.ndarray:
    """Frecuencia relativa de cada aminoacido estandar.

    Se normaliza por la longitud para que el vector no codifique el tamano del
    peptido: la longitud se anade aparte como descriptor explicito en physchem,
    donde se puede inspeccionar y, si hace falta, quitar.
    """
    vec = np.zeros(len(AMINO_ACIDS), dtype=np.float64)
    if not sequence:
        return vec
    for aa in sequence:
        idx = AA_INDEX.get(aa)
        if idx is not None:
            vec[idx] += 1.0
    return vec / len(sequence)


def composition_matrix(sequences: list[str]) -> np.ndarray:
    """Matriz (n_secuencias, 20)."""
    return np.vstack([composition_vector(s) for s in sequences])
