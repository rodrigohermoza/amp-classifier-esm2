"""Tests de la division de datos.

El mas importante del proyecto es test_no_cluster_spans_splits: si falla,
ninguna metrica publicada es valida, porque habria fuga de informacion por
homologia entre entrenamiento y prueba.
"""

import importlib.util

import numpy as np
import pandas as pd
import pytest

from ampesm.split import assign_splits, cluster_builtin, kmer_set

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
        mapping = cluster_builtin(seqs, identity=0.4, coverage=0.5)
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
        mapping = cluster_builtin(seqs, identity=0.4, coverage=0.5)
        assert mapping["wt"] == mapping["mut"]

    def test_unrelated_sequences_stay_apart(self):
        seqs = {"cationic": "KKKKRRRRKKKKRRRR", "acidic": "DDDDEEEEDDDDEEEE"}
        mapping = cluster_builtin(seqs, identity=0.4, coverage=0.5)
        assert mapping["cationic"] != mapping["acidic"]

    def test_all_sequences_receive_a_cluster(self):
        seqs = {f"s{i}": seq for i, seq in enumerate(
            ["GIGKFLHSAKKF", "GIGKFLHSAKKW", "DDDDEEEE", "WWWWWWWW", "KRKRKRKR"]
        )}
        mapping = cluster_builtin(seqs, identity=0.4, coverage=0.5)
        assert set(mapping) == set(seqs)

    def test_deterministic(self):
        seqs = {f"s{i}": s for i, s in enumerate(["GIGKFLHSAKKF", "GIGKFLHSAKKW", "DDDDEEEE"])}
        assert cluster_builtin(seqs, 0.4, 0.5) == cluster_builtin(seqs, 0.4, 0.5)
