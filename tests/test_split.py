"""Tests de la division de datos.

El mas importante del proyecto es test_no_cluster_spans_splits: si falla,
ninguna metrica publicada es valida, porque habria fuga de informacion por
homologia entre entrenamiento y prueba.
"""

import importlib.util

import numpy as np
import pandas as pd
import pytest

from ampesm.split import (
    assign_splits,
    audit_train_test_identity,
    cluster_builtin,
    kmer_set,
)

FRACTIONS = {"train": 0.70, "val": 0.15, "test": 0.15}


def make_df(cluster_sizes, positive_rate=0.25, seed=0):
    """DataFrame sintetico con clusteres de tamanos dados."""
    rng = np.random.default_rng(seed)
    rows = []
    for cid, size in enumerate(cluster_sizes):
        for j in range(size):
            rows.append(
                {
                    "id": f"seq_{cid}_{j}",
                    "cluster_id": cid,
                    "label": int(rng.random() < positive_rate),
                }
            )
    return pd.DataFrame(rows)


class TestKmers:
    def test_kmer_extraction(self):
        assert kmer_set("ABCDE", 3) == frozenset({"ABC", "BCD", "CDE"})

    def test_short_sequence_returns_itself(self):
        """Una secuencia mas corta que k no tiene k-meros; se usa ella misma.

        Sin este caso, un peptido de 2 residuos tendria un conjunto vacio, su
        Jaccard contra todo seria 0 y nunca se agruparia con nada.
        """
        assert kmer_set("AB", 3) == frozenset({"AB"})

    def test_identical_sequences_have_identical_kmers(self):
        assert kmer_set("GIGKFLHSAKKF", 3) == kmer_set("GIGKFLHSAKKF", 3)


class TestAssignSplits:
    def test_no_cluster_spans_splits(self):
        """EL test del proyecto: ningun cluster puede estar en dos splits.

        Si este test falla, hay fuga por homologia y todas las metricas estan
        infladas. Ver docs/04-division-clusters.md.
        """
        df = make_df([50, 40, 30, 25, 20] + [10] * 20 + [1] * 100)
        assignment = assign_splits(df, FRACTIONS, seed=42)
        df["split"] = df.cluster_id.map(assignment)

        per_cluster = df.groupby("cluster_id").split.nunique()
        assert (per_cluster == 1).all(), f"clusteres repartidos: {per_cluster[per_cluster > 1]}"

    def test_all_clusters_assigned(self):
        df = make_df([10] * 30)
        assignment = assign_splits(df, FRACTIONS, seed=42)
        assert set(assignment) == set(df.cluster_id.unique())
        assert set(assignment.values()) <= set(FRACTIONS)

    def test_fractions_approximately_respected(self):
        """Con muchos clusteres pequenos, el reparto voraz debe acercarse mucho."""
        df = make_df([5] * 200)
        assignment = assign_splits(df, FRACTIONS, seed=42)
        df["split"] = df.cluster_id.map(assignment)

        for split, target in FRACTIONS.items():
            actual = (df.split == split).mean()
            assert abs(actual - target) < 0.05, f"{split}: esperado {target}, obtenido {actual}"

    def test_every_split_receives_data(self):
        df = make_df([20] * 10)
        assignment = assign_splits(df, FRACTIONS, seed=42)
        df["split"] = df.cluster_id.map(assignment)
        for split in FRACTIONS:
            assert (df.split == split).sum() > 0, f"{split} quedo vacio"

    def test_one_huge_cluster_does_not_break_assignment(self):
        """Caso adverso: un cluster con el 40% de los datos.

        Pasa de verdad: una familia de AMPs muy representada puede dominar. El
        reparto tiene que seguir siendo valido aunque las proporciones salgan
        desviadas, porque la alternativa (partir el cluster) seria fuga.
        """
        df = make_df([400] + [5] * 120)
        assignment = assign_splits(df, FRACTIONS, seed=42)
        df["split"] = df.cluster_id.map(assignment)
        assert (df.groupby("cluster_id").split.nunique() == 1).all()
        for split in FRACTIONS:
            assert (df.split == split).sum() > 0

    def test_deterministic_given_seed(self):
        df = make_df([10] * 40)
        a = assign_splits(df, FRACTIONS, seed=7)
        b = assign_splits(df, FRACTIONS, seed=7)
        assert a == b

    def test_class_balance_roughly_preserved(self):
        df = make_df([5] * 200, positive_rate=0.25, seed=3)
        assignment = assign_splits(df, FRACTIONS, seed=42)
        df["split"] = df.cluster_id.map(assignment)
        overall = df.label.mean()
        for split in FRACTIONS:
            rate = df[df.split == split].label.mean()
            assert abs(rate - overall) < 0.12, f"{split}: tasa {rate} frente a {overall}"


@pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython no instalado; el respaldo de clustering lo necesita",
)
class TestBuiltinClustering:
    def test_identical_sequences_cluster_together(self):
        seqs = {"a": "GIGKFLHSAKKFGKAFVGEIMNS", "b": "GIGKFLHSAKKFGKAFVGEIMNS"}
        mapping, _ = cluster_builtin(seqs, identity=0.4, coverage=0.5)
        assert mapping["a"] == mapping["b"]

    def test_near_identical_sequences_cluster_together(self):
        """Dos mutaciones no deben bastar para separar: es el caso que motiva todo.

        Si estas dos acabaran en splits distintos, el modelo veria una en
        entrenamiento y acertaria la otra en prueba por memorizacion.
        """
        seqs = {
            "wt": "GIGKFLHSAKKFGKAFVGEIMNS",
            "mut": "GIGKFLHSAKKFGKAFVGEIMAA",
        }
        mapping, _ = cluster_builtin(seqs, identity=0.4, coverage=0.5)
        assert mapping["wt"] == mapping["mut"]

    def test_unrelated_sequences_stay_apart(self):
        seqs = {"cationic": "KKKKRRRRKKKKRRRR", "acidic": "DDDDEEEEDDDDEEEE"}
        mapping, _ = cluster_builtin(seqs, identity=0.4, coverage=0.5)
        assert mapping["cationic"] != mapping["acidic"]

    def test_all_sequences_receive_a_cluster(self):
        seqs = {
            f"s{i}": seq
            for i, seq in enumerate(
                ["GIGKFLHSAKKF", "GIGKFLHSAKKW", "DDDDEEEE", "WWWWWWWW", "KRKRKRKR"]
            )
        }
        mapping, _ = cluster_builtin(seqs, identity=0.4, coverage=0.5)
        assert set(mapping) == set(seqs)

    def test_deterministic(self):
        seqs = {f"s{i}": s for i, s in enumerate(["GIGKFLHSAKKF", "GIGKFLHSAKKW", "DDDDEEEE"])}
        assert cluster_builtin(seqs, 0.4, 0.5)[0] == cluster_builtin(seqs, 0.4, 0.5)[0]


@pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython no instalado; el respaldo de clustering lo necesita",
)
class TestPrefilterRegression:
    """Regresiones de P-03: el prefiltro anulaba el umbral de identidad.

    El prefiltro viejo (Jaccard de 3-meros >= 0.15) era mas estricto que el
    propio umbral del 40%, asi que los parientes que el umbral pretende agrupar
    nunca llegaban a alinearse. Estos tests fijan el comportamiento corregido.
    """

    def test_short_peptide_contained_in_long_one_clusters(self):
        """El caso que el prefiltro de Jaccard descartaba siempre.

        Un 6-mero identico a un fragmento de un 50-mero tiene identidad 1.00 y
        cobertura 1.00, pero su Jaccard de 3-meros contra el largo es ~0.07,
        asi que con el umbral viejo nunca se alineaban. Con contencion, da 1.00.
        """
        largo = "GIGKFLHSAKKFGKAFVGEIMNSKKWWRRPLFVACDEGHIKLMNPQ"
        corto = largo[10:16]
        mapping, _ = cluster_builtin({"largo": largo, "corto": corto}, 0.4, 0.5)
        assert (
            mapping["largo"] == mapping["corto"]
        ), "un peptido corto contenido entero en uno largo debe agruparse con el"

    def test_moderately_diverged_pair_clusters(self):
        """Dos secuencias con ~40-50% de identidad deben caer en el mismo cluster.

        Es exactamente el rango que el umbral de 0.4 existe para capturar y el
        que el prefiltro viejo descartaba en el 99-100% de los casos.
        """
        base = "GIGKFLHSAKKFGKAFVGEIMNS"
        # Se mutan 11 de 23 posiciones -> ~52% de identidad.
        mutada = list(base)
        for i in range(0, len(base), 2):
            mutada[i] = "A" if base[i] != "A" else "W"
        mapping, _ = cluster_builtin({"base": base, "mut": "".join(mutada)}, 0.4, 0.5)
        assert mapping["base"] == mapping["mut"]

    def test_diagnostics_report_containment_criterion(self):
        _mapping, diag = cluster_builtin({"a": "GIGKFLHSAKKF"}, 0.4, 0.5)
        assert diag["prefilter_criterion"] == "containment"
        assert diag["prefilter"] == 0.05

    def test_truly_unrelated_still_separate(self):
        """Bajar el prefiltro no debe fusionar secuencias sin relacion."""
        seqs = {"cat": "KKKKRRRRKKKKRRRR", "aci": "DDDDEEEEDDDDEEEE", "aro": "WWWWFFFFWWWWFFFF"}
        mapping, _ = cluster_builtin(seqs, 0.4, 0.5)
        assert len(set(mapping.values())) == 3


class TestDegenerateSplits:
    """Regresiones de P-06: splits vacios o monoclase en silencio."""

    def test_correlated_labels_still_give_both_classes(self):
        """El caso real: familias de AMPs donde la etiqueta sigue al cluster.

        Antes salia val con 0% de positivos y test con 100%. La siembra previa
        garantiza ambas clases en cada split siempre que haya material para
        ello, es decir al menos un cluster de cada clase por split.
        """
        rows = []
        for cid in range(12):
            label = 1 if cid < 6 else 0
            for j in range(50):
                rows.append({"id": f"s{cid}_{j}", "cluster_id": cid, "label": label})
        df = pd.DataFrame(rows)

        assignment = assign_splits(df, FRACTIONS, seed=42)
        df["split"] = df.cluster_id.map(assignment)

        for split in FRACTIONS:
            sub = df[df.split == split]
            assert len(sub) > 0, f"{split} quedo vacio"
            assert sub.label.nunique() == 2, (
                f"{split} tiene una sola clase: "
                f"{int((sub.label == 1).sum())}+ / {int((sub.label == 0).sum())}-"
            )

    def test_too_few_positive_clusters_is_impossible_not_hidden(self):
        """Con 2 clusteres positivos y 3 splits, no hay solucion valida.

        Darle positivos a los tres splits exigiria partir un cluster, que es
        precisamente la fuga que se quiere evitar. Lo correcto no es que
        assign_splits haga magia, sino que la division resultante sea
        detectable: la validacion de build_splits aborta antes de escribir
        splits.csv. Este test documenta el limite y comprueba que al menos la
        siembra reparte los positivos disponibles entre splits distintos en vez
        de amontonarlos en uno.
        """
        rows = []
        for cid in range(5):
            label = 1 if cid < 2 else 0
            for j in range(100):
                rows.append({"id": f"s{cid}_{j}", "cluster_id": cid, "label": label})
        df = pd.DataFrame(rows)

        assignment = assign_splits(df, FRACTIONS, seed=42)
        df["split"] = df.cluster_id.map(assignment)

        splits_con_positivos = {
            split for split in FRACTIONS if (df[df.split == split].label == 1).any()
        }
        assert len(splits_con_positivos) == 2, (
            "los 2 clusteres positivos deben quedar en splits distintos; "
            f"quedaron en {splits_con_positivos}"
        )
        # Y el invariante de no fuga se mantiene aunque la division sea invalida.
        assert (df.groupby("cluster_id").split.nunique() == 1).all()

    def test_one_dominant_cluster_still_gives_both_classes(self):
        rows = [{"id": f"big_{j}", "cluster_id": 0, "label": 1} for j in range(600)]
        for cid in range(1, 101):
            rows.append({"id": f"s{cid}", "cluster_id": cid, "label": cid % 3 == 0})
        df = pd.DataFrame(rows)

        assignment = assign_splits(df, FRACTIONS, seed=42)
        df["split"] = df.cluster_id.map(assignment)
        for split in FRACTIONS:
            sub = df[df.split == split]
            assert len(sub) > 0 and sub.label.nunique() == 2, f"{split} degenerado"

    def test_seeding_does_not_break_no_leakage(self):
        """La siembra no puede romper el invariante que todo esto protege."""
        df = make_df([200, 150, 80] + [20] * 30 + [1] * 200)
        assignment = assign_splits(df, FRACTIONS, seed=42)
        df["split"] = df.cluster_id.map(assignment)
        assert (df.groupby("cluster_id").split.nunique() == 1).all()


@pytest.mark.skipif(importlib.util.find_spec("Bio") is None, reason="Biopython no instalado")
class TestIdentityAudit:
    """P-04: la verificacion directa de la premisa del proyecto."""

    def test_detects_planted_leakage(self):
        """Si se mete a mano un pariente identico en train, debe aparecer."""
        seq = "GIGKFLHSAKKFGKAFVGEIMNS"
        df = pd.DataFrame(
            [
                {"id": "tr1", "sequence": seq, "split": "train", "label": 1, "cluster_id": 0},
                {"id": "te1", "sequence": seq, "split": "test", "label": 1, "cluster_id": 1},
            ]
        )
        audit = audit_train_test_identity(df, threshold=0.4)
        assert audit["n_over_threshold"] == 1
        assert audit["fraction_over_threshold"] == 1.0
        assert audit["max_identity_distribution"]["max"] == 1.0

    def test_clean_split_reports_no_leakage(self):
        df = pd.DataFrame(
            [
                {
                    "id": "tr1",
                    "sequence": "KKKKRRRRKKKKRRRR",
                    "split": "train",
                    "label": 1,
                    "cluster_id": 0,
                },
                {
                    "id": "te1",
                    "sequence": "DDDDEEEEDDDDEEEE",
                    "split": "test",
                    "label": 0,
                    "cluster_id": 1,
                },
            ]
        )
        audit = audit_train_test_identity(df, threshold=0.4)
        assert audit["n_over_threshold"] == 0

    def test_empty_split_returns_error_not_crash(self):
        df = pd.DataFrame(
            [{"id": "tr1", "sequence": "GIGKFL", "split": "train", "label": 1, "cluster_id": 0}]
        )
        assert "error" in audit_train_test_identity(df, threshold=0.4)
