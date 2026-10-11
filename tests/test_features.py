"""Los descriptores se verifican contra valores publicados, no contra si mismos.

Un test que solo compara la salida de la funcion con la salida de la funcion no
detecta un error en la formula. Los valores de referencia de aqui vienen de la
literatura sobre peptidos conocidos.
"""

import numpy as np

from ampesm.features import composition, physchem

# Peptidos de referencia con propiedades publicadas.
MAGAININ_2 = "GIGKFLHSAKKFGKAFVGEIMNS"  # AMP cationico, carga publicada +3 a +4
LL37 = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES"  # AMP humano, carga publicada +6
POLY_K = "KKKKKKKKKK"
POLY_E = "EEEEEEEEEE"
POLY_L = "LLLLLLLLLLLLLLL"
AMPHIPATHIC = "LKKLLKLLKKLLKLL"  # anfipatico de diseno


class TestComposition:
    def test_frequencies_sum_to_one(self):
        vec = composition.composition_vector(MAGAININ_2)
        assert np.isclose(vec.sum(), 1.0)

    def test_poly_k_is_all_lysine(self):
        vec = composition.composition_vector(POLY_K)
        assert np.isclose(vec[composition.AA_INDEX["K"]], 1.0)
        assert np.isclose(vec.sum(), 1.0)

    def test_empty_sequence_is_zeros(self):
        assert composition.composition_vector("").sum() == 0.0

    def test_length_invariance(self):
        """Duplicar la secuencia no cambia la composicion: esta normalizada."""
        single = composition.composition_vector(MAGAININ_2)
        doubled = composition.composition_vector(MAGAININ_2 * 2)
        assert np.allclose(single, doubled)

    def test_matrix_shape(self):
        matrix = composition.composition_matrix([MAGAININ_2, LL37, POLY_K])
        assert matrix.shape == (3, 20)


class TestNetCharge:
    def test_magainin_2_matches_literature(self):
        """Carga publicada de la magainina 2 a pH neutro: +3 a +4."""
        charge = physchem.net_charge(MAGAININ_2, ph=7.4)
        assert 2.5 < charge < 4.0

    def test_ll37_matches_literature(self):
        """Carga publicada de LL-37: +6."""
        charge = physchem.net_charge(LL37, ph=7.4)
        assert 5.0 < charge < 7.0

    def test_poly_lysine_is_strongly_positive(self):
        assert physchem.net_charge(POLY_K) > 9.0

    def test_poly_glutamate_is_strongly_negative(self):
        assert physchem.net_charge(POLY_E) < -9.0

    def test_charge_decreases_with_ph(self):
        """La curva de carga es monotona decreciente en pH."""
        charges = [physchem.net_charge(MAGAININ_2, ph) for ph in (2.0, 5.0, 7.4, 10.0, 12.0)]
        assert all(a > b for a, b in zip(charges, charges[1:]))

    def test_histidine_is_partially_protonated(self):
        """A pH 7.4 la histidina (pKa 6.0) aporta mucho menos de +1.

        Es el motivo de usar Henderson-Hasselbalch en vez de contar cargas
        enteras.
        """
        with_his = physchem.net_charge("AAHAA", ph=7.4)
        without = physchem.net_charge("AAAAA", ph=7.4)
        contribution = with_his - without
        assert 0.0 < contribution < 0.2


class TestIsoelectricPoint:
    def test_charge_is_zero_at_pi(self):
        pi = physchem.isoelectric_point(MAGAININ_2)
        assert abs(physchem.net_charge(MAGAININ_2, pi)) < 0.01

    def test_cationic_peptide_has_basic_pi(self):
        assert physchem.isoelectric_point(POLY_K) > 10.0

    def test_anionic_peptide_has_acidic_pi(self):
        assert physchem.isoelectric_point(POLY_E) < 5.0


class TestHydrophobicMoment:
    def test_homopolymer_has_near_zero_moment(self):
        """Sin variacion de hidrofobicidad no hay anfipaticidad."""
        assert physchem.hydrophobic_moment(POLY_L) < 0.1

    def test_designed_amphipathic_has_high_moment(self):
        """Un patron L/K que agrupa los hidrofobicos en una cara de la helice."""
        assert physchem.hydrophobic_moment(AMPHIPATHIC) > 0.5

    def test_amphipathic_exceeds_homopolymer(self):
        assert physchem.hydrophobic_moment(AMPHIPATHIC) > physchem.hydrophobic_moment(POLY_L)

    def test_shuffling_changes_moment(self):
        """El momento depende del ORDEN, no solo de la composicion.

        Es la propiedad que lo distingue del vector de composicion: dos
        secuencias con los mismos aminoacidos en otro orden tienen momentos
        distintos.
        """
        ordered = physchem.hydrophobic_moment("LLLLLLLKKKKKKKK")
        alternating = physchem.hydrophobic_moment("LKLKLKLKLKLKLKL")
        assert not np.isclose(ordered, alternating)


class TestOtherDescriptors:
    def test_gravy_sign(self):
        assert physchem.gravy(POLY_L) > 0  # leucina es hidrofobica
        assert physchem.gravy(POLY_E) < 0  # glutamato es hidrofilico

    def test_hydrophobic_fraction_bounds(self):
        assert physchem.hydrophobic_fraction(POLY_L) == 1.0
        assert physchem.hydrophobic_fraction(POLY_K) == 0.0

    def test_aromaticity(self):
        assert physchem.aromaticity("FFWWYY") == 1.0
        assert physchem.aromaticity(POLY_K) == 0.0

    def test_vector_length_matches_names(self):
        vec = physchem.physchem_vector(MAGAININ_2)
        assert len(vec) == len(physchem.FEATURE_NAMES)

    def test_no_nan_or_inf(self):
        """Ningun descriptor debe producir NaN o infinito.

        Importa porque un NaN se propaga en silencio hasta el entrenamiento y
        alli se manifiesta como un error oscuro de scikit-learn.
        """
        for seq in (MAGAININ_2, LL37, POLY_K, POLY_E, POLY_L, "A", "AC", "W" * 50):
            vec = physchem.physchem_vector(seq)
            assert np.all(np.isfinite(vec)), f"valor no finito en {seq!r}: {vec}"

    def test_empty_sequence_does_not_crash(self):
        vec = physchem.physchem_vector("")
        assert np.all(np.isfinite(vec))
