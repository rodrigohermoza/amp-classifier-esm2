"""Descriptores fisicoquimicos.

Las escalas se escriben explicitas en este archivo, con su fuente, en vez de
importarlas de una libreria. Dos razones: que se puedan auditar sin salir del
repositorio, y que el proyecto no dependa de que una libreria cambie de valores
entre versiones y mueva los numeros publicados.

Justificacion de cada descriptor en docs/05-features-baselines.md.
"""

from __future__ import annotations

import math

import numpy as np

# Kyte & Doolittle (1982), hidropatia. Base del GRAVY.
KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

# Eisenberg et al. (1984), escala de hidrofobicidad consenso normalizada.
# Se usa para el momento hidrofobico, que es el uso para el que se diseno.
EISENBERG = {
    "A": 0.62,
    "R": -2.53,
    "N": -0.78,
    "D": -0.90,
    "C": 0.29,
    "Q": -0.85,
    "E": -0.74,
    "G": 0.48,
    "H": -0.40,
    "I": 1.38,
    "L": 1.06,
    "K": -1.50,
    "M": 0.64,
    "F": 1.19,
    "P": 0.12,
    "S": -0.18,
    "T": -0.05,
    "W": 0.81,
    "Y": 0.26,
    "V": 1.08,
}

# Energias libres de transferencia ciclohexano -> agua de las cadenas laterales
# (kcal/mol), de Radzicka y Wolfenden (1988). Es la escala sobre la que se
# define el indice de Boman. La arginina domina (14.92), seguida de D y E; los
# alifaticos quedan en el extremo negativo.
#
# P-08: aqui habia antes la escala interfacial de Wimley-White con el signo
# invertido (se reconocia por W=1.85, F=1.13, E=-2.02). Con esa tabla la
# arginina quedaba al nivel de N y Q, justo al reves de lo que hace util al
# indice para AMPs cationicos, y LL-37 daba 0.35 en vez de ~3.0.
RADZICKA_WOLFENDEN = {
    "R": 14.92,
    "D": 8.72,
    "E": 6.81,
    "N": 6.64,
    "K": 5.55,
    "Q": 5.54,
    "H": 4.66,
    "S": 3.40,
    "T": 2.57,
    "Y": 0.14,
    "P": 0.00,
    "G": -0.94,
    "C": -1.28,
    "A": -1.81,
    "W": -2.33,
    "M": -2.35,
    "F": -2.98,
    "V": -4.04,
    "I": -4.92,
    "L": -4.92,
}

# Escala interfacial de Wimley-White (kcal/mol), conservada como descriptor
# propio. Mide la preferencia por la interfaz lipidica, no la union a proteinas.
WIMLEY_WHITE_INTERFACIAL = {
    "A": -0.17,
    "R": -0.81,
    "N": -0.42,
    "D": -1.23,
    "C": 0.24,
    "Q": -0.58,
    "E": -2.02,
    "G": -0.01,
    "H": -0.96,
    "I": 0.31,
    "L": 0.56,
    "K": -0.99,
    "M": 0.23,
    "F": 1.13,
    "P": -0.45,
    "S": -0.13,
    "T": -0.14,
    "W": 1.85,
    "Y": 0.94,
    "V": -0.07,
}

# pKa de las cadenas laterales ionizables y de los extremos (valores tipo EMBOSS).
PKA_SIDE_CHAIN = {"D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07, "H": 6.00, "K": 10.53, "R": 12.48}
PKA_N_TERM = 8.60
PKA_C_TERM = 3.60
POSITIVE_RESIDUES = frozenset("KRH")
NEGATIVE_RESIDUES = frozenset("DECY")

HYDROPHOBIC = frozenset("ACFILMVWY")
AROMATIC = frozenset("FWY")

FEATURE_NAMES: list[str] = [
    "length",
    "net_charge",
    "isoelectric_point",
    "gravy",
    "hydrophobic_moment",
    "hydrophobic_fraction",
    "aromaticity",
    "aliphatic_index",
    "boman_index",
    "interfacial_hydrophobicity",
    "charge_density",
]


def net_charge(sequence: str, ph: float = 7.4) -> float:
    """Carga neta por Henderson-Hasselbalch.

    Cada grupo ionizable contribuye con su fraccion protonada o desprotonada al
    pH dado, no con una carga entera. Contar +1 por cada K/R y -1 por cada D/E
    seria mas simple pero equivocado: a pH 7.4, la histidina esta solo
    parcialmente protonada (pKa 6.0), y los extremos de la cadena tambien
    contribuyen.

    La carga positiva neta es la propiedad mas importante de un AMP cationico:
    es lo que lo atrae selectivamente a la membrana bacteriana, mas anionica que
    la de las celulas animales.
    """
    if not sequence:
        return 0.0

    charge = 0.0
    # Extremos de la cadena.
    charge += 1.0 / (1.0 + 10 ** (ph - PKA_N_TERM))
    charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERM - ph))

    for aa in sequence:
        pka = PKA_SIDE_CHAIN.get(aa)
        if pka is None:
            continue
        if aa in POSITIVE_RESIDUES:
            charge += 1.0 / (1.0 + 10 ** (ph - pka))
        elif aa in NEGATIVE_RESIDUES:
            charge -= 1.0 / (1.0 + 10 ** (pka - ph))
    return charge


def isoelectric_point(sequence: str, lo: float = 0.0, hi: float = 20.0) -> float:
    """pH al que la carga neta es cero, por biseccion.

    La curva de carga es monotona decreciente en pH, asi que la biseccion
    converge siempre. 60 iteraciones agotan la precision de doble sobre este
    intervalo; las 100 de antes eran trabajo sin ganancia.

    El intervalo llega a 20 y no a 14 porque un peptido muy basico puede no
    tener su pI por debajo de 14: con 35 argininas en 50 residuos la carga sigue
    siendo positiva a pH 14, y la version anterior devolvia 14.0 en silencio, que
    es el borde del intervalo y no un pI. Si aun asi no hay cambio de signo, se
    devuelve NaN en vez de un numero inventado.
    """
    if not sequence:
        return 0.0
    if net_charge(sequence, lo) <= 0 or net_charge(sequence, hi) >= 0:
        return float("nan")
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if net_charge(sequence, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def gravy(sequence: str) -> float:
    """Hidrofobicidad promedio (grand average of hydropathy)."""
    if not sequence:
        return 0.0
    vals = [KYTE_DOOLITTLE[aa] for aa in sequence if aa in KYTE_DOOLITTLE]
    return float(np.mean(vals)) if vals else 0.0


def hydrophobic_moment(sequence: str, window: int = 11, angle_deg: float = 100.0) -> float:
    """Momento hidrofobico de Eisenberg: la medida de anfipaticidad.

    Cada residuo se trata como un vector: su direccion la da su posicion
    angular en la helice (100 grados de giro por residuo en una helice alfa) y
    su magnitud, su hidrofobicidad. La suma vectorial es grande cuando los
    residuos hidrofobicos quedan agrupados en una cara de la helice y los
    hidrofilicos en la otra.

    Eso es exactamente la arquitectura de la mayoria de los AMPs helicoidales:
    la cara hidrofobica se inserta en el interior lipidico de la membrana y la
    hidrofilica queda hacia el agua. Conceptualmente, es el descriptor mas
    informativo del conjunto.

    Se devuelve el maximo sobre todas las ventanas de tamano `window`, no el
    valor global, porque un peptido puede tener un segmento anfipatico dentro de
    una secuencia que en conjunto no lo es, y es el segmento el que interactua
    con la membrana.
    """
    if not sequence:
        return 0.0

    angle = math.radians(angle_deg)
    values = [EISENBERG.get(aa, 0.0) for aa in sequence]
    win = min(window, len(values))

    best = 0.0
    for start in range(len(values) - win + 1):
        sin_sum = 0.0
        cos_sum = 0.0
        for i in range(win):
            h = values[start + i]
            sin_sum += h * math.sin(i * angle)
            cos_sum += h * math.cos(i * angle)
        moment = math.sqrt(sin_sum**2 + cos_sum**2) / win
        best = max(best, moment)
    return best


def hydrophobic_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    return sum(1 for aa in sequence if aa in HYDROPHOBIC) / len(sequence)


def aromaticity(sequence: str) -> float:
    if not sequence:
        return 0.0
    return sum(1 for aa in sequence if aa in AROMATIC) / len(sequence)


def aliphatic_index(sequence: str) -> float:
    """Indice alifatico de Ikai (1980): volumen relativo de cadenas alifaticas."""
    if not sequence:
        return 0.0
    n = len(sequence)
    pct = {aa: 100.0 * sequence.count(aa) / n for aa in "AVIL"}
    return pct["A"] + 2.9 * pct["V"] + 3.9 * (pct["I"] + pct["L"])


def boman_index(sequence: str) -> float:
    """Indice de Boman: potencial de interaccion con proteinas.

    Suma de las energias libres de transferencia ciclohexano -> agua de las
    cadenas laterales, dividida por la longitud. Valores altos indican tendencia
    a unirse a proteinas; valores bajos, preferencia por la membrana.

    Verificado contra valores publicados: LL-37 ~3.0 kcal/mol (APD3 reporta
    2,9-3,0) y magainina-2 ~0.42.
    """
    if not sequence:
        return 0.0
    total = sum(RADZICKA_WOLFENDEN.get(aa, 0.0) for aa in sequence)
    return total / len(sequence)


def interfacial_hydrophobicity(sequence: str) -> float:
    """Hidrofobicidad interfacial media (escala de Wimley-White).

    Mide la preferencia por la interfaz agua-lipido, que es donde se inserta un
    AMP. Correlaciona con GRAVY, pero pondera distinto los aromaticos: el
    triptofano, que se acumula en las interfaces de membrana, es el residuo con
    el valor mas alto de la escala.
    """
    if not sequence:
        return 0.0
    return sum(WIMLEY_WHITE_INTERFACIAL.get(aa, 0.0) for aa in sequence) / len(sequence)


def physchem_vector(sequence: str, ph: float = 7.4, window: int = 11, angle: float = 100.0):
    """Los 10 descriptores en el orden de FEATURE_NAMES."""
    charge = net_charge(sequence, ph)
    length = len(sequence)
    return np.array(
        [
            float(length),
            charge,
            isoelectric_point(sequence),
            gravy(sequence),
            hydrophobic_moment(sequence, window, angle),
            hydrophobic_fraction(sequence),
            aromaticity(sequence),
            aliphatic_index(sequence),
            boman_index(sequence),
            interfacial_hydrophobicity(sequence),
            charge / length if length else 0.0,
        ],
        dtype=np.float64,
    )


def physchem_matrix(
    sequences: list[str], ph: float = 7.4, window: int = 11, angle: float = 100.0
) -> np.ndarray:
    return np.vstack([physchem_vector(s, ph, window, angle) for s in sequences])
