"""Tests de la division de datos.

El mas importante del proyecto es test_no_cluster_spans_splits: si falla,
ninguna metrica publicada es valida, porque habria fuga de informacion por
homologia entre entrenamiento y prueba.
"""

import importlib.util
import random

import numpy as np
import pandas as pd
import pytest

from ampesm.split import (
    _positive_rate_drift,
    _publication_gates,
    assign_splits,
    audit_train_test_identity,
    cluster_builtin,
    cluster_ids_from_pairs,
    cluster_mmseqs2,
    kmer_set,
    merge_contained_clusters,
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
        audit = audit_train_test_identity(df, threshold=0.4, min_coverage=0.5)
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
        audit = audit_train_test_identity(df, threshold=0.4, min_coverage=0.5)
        assert audit["n_over_threshold"] == 0

    def test_empty_split_returns_error_not_crash(self):
        df = pd.DataFrame(
            [{"id": "tr1", "sequence": "GIGKFL", "split": "train", "label": 1, "cluster_id": 0}]
        )
        assert "error" in audit_train_test_identity(df, threshold=0.4, min_coverage=0.5)


@pytest.mark.skipif(importlib.util.find_spec("Bio") is None, reason="Biopython no instalado")
class TestLeakageDetector:
    """P-14, P-15 y P-18: el detector de fuga debe DETECTAR FUGA.

    Suena obvio y sin embargo las dos primeras versiones no lo hacian, cada una
    por un motivo distinto, y ninguna fallaba de forma visible. Estos tests
    plantan fuga a proposito y exigen que se vea.
    """

    AA = "ACDEFGHIKLMNPQRSTVWY"

    def _frame(self, n_leak, identity, seed=11):
        rng = random.Random(seed)
        rows, train = [], []
        for i in range(150):
            s = "".join(rng.choices(self.AA, k=rng.randint(15, 45)))
            train.append(s)
            rows.append(
                {"id": f"tr{i}", "sequence": s, "split": "train", "label": 1, "cluster_id": i}
            )
        for j in range(60):
            if j < n_leak:
                chars = list(train[j])
                for k in range(len(chars)):
                    if rng.random() > identity:
                        chars[k] = rng.choice(self.AA)
                s = "".join(chars)
            else:
                s = "".join(rng.choices(self.AA, k=rng.randint(15, 45)))
            rows.append(
                {"id": f"te{j}", "sequence": s, "split": "test", "label": 1, "cluster_id": 1000 + j}
            )
        return pd.DataFrame(rows)

    def test_no_leakage_is_not_flagged(self):
        audit = audit_train_test_identity(self._frame(0, 0.0), 0.4, 0.5, null_replicates=20, seed=5)
        assert not audit["leakage"]["significant"], f"falso positivo: {audit['leakage']}"

    def test_planted_leakage_is_detected(self):
        """Regresion de P-18.

        El estadistico anterior (diferencia de proporciones agregadas) no
        detectaba ni una fuga del 17% al 80% de identidad, porque comparaba dos
        proporciones ya saturadas por el ruido de fondo.
        """
        audit = audit_train_test_identity(
            self._frame(12, 0.75), 0.4, 0.5, null_replicates=20, seed=5
        )
        leak = audit["leakage"]
        assert leak["significant"], f"no detecto 12 de 60 parientes al 75%: {leak}"
        assert leak["n_flagged"] > leak["n_expected_by_chance"]

    def test_short_subject_is_not_skipped(self):
        """Regresion de P-14: el corte temprano usaba una cota invertida.

        Una consulta larga con un sujeto corto contenido entero daba identidad
        1.00, pero la cota n_shared/n_kmers_consulta la dejaba por debajo del
        prefiltro y el 'break' descartaba ese sujeto y todos los siguientes.
        """
        largo = "GIGKFLHSAKKFGKAFVGEIMNSKKWWRRPLFVACDEGHIKLMNPQRSTV"
        df = pd.DataFrame(
            [
                {
                    "id": "tr0",
                    "sequence": "QQQQQQQQQQQQQQ",
                    "split": "train",
                    "label": 1,
                    "cluster_id": 0,
                },
                {
                    "id": "tr1",
                    "sequence": largo[10:18],
                    "split": "train",
                    "label": 1,
                    "cluster_id": 1,
                },
                {"id": "te0", "sequence": largo, "split": "test", "label": 1, "cluster_id": 2},
            ]
        )
        audit = audit_train_test_identity(df, 0.4, 0.5, null_model=False)
        assert (
            audit["max_identity_distribution"]["max"] == 1.0
        ), "el fragmento exacto corto debe encontrarse"


class TestValAudit:
    """P-27: la auditoria tiene que poder mirar val, no solo test.

    val fija los hiperparametros y el umbral de decision, y ese umbral se
    aplica a test. Homologia val-train sesga cinco de las siete metricas de
    test sin que test tenga un solo pariente en train.
    """

    def _df(self):
        seq = "GIGKFLHSAKKFGKAFVGEIMNS"
        return pd.DataFrame(
            [
                {"id": "tr1", "sequence": seq, "split": "train", "label": 1, "cluster_id": 0},
                # El pariente plantado esta en val, NO en test.
                {"id": "va1", "sequence": seq, "split": "val", "label": 1, "cluster_id": 1},
                {
                    "id": "te1",
                    "sequence": "DDDDEEEEDDDDEEEEDDDD",
                    "split": "test",
                    "label": 0,
                    "cluster_id": 2,
                },
            ]
        )

    def test_finds_leakage_that_the_test_audit_cannot_see(self):
        df = self._df()

        on_test = audit_train_test_identity(df, threshold=0.4, min_coverage=0.5, null_model=False)
        on_val = audit_train_test_identity(
            df, threshold=0.4, min_coverage=0.5, null_model=False, query_split="val"
        )

        assert on_test["n_over_threshold"] == 0, "test esta limpio, y el audit de test lo dice"
        assert (
            on_val["n_over_threshold"] == 1
        ), "la fuga vive en val; si esto falla, la auditoria no mira donde hace falta"

    def test_report_says_which_splits_it_compared(self):
        audit = audit_train_test_identity(
            self._df(), threshold=0.4, min_coverage=0.5, null_model=False, query_split="val"
        )
        assert audit["query_split"] == "val"
        assert audit["reference_split"] == "train"

    def test_missing_split_returns_error_not_crash(self):
        audit = audit_train_test_identity(
            self._df(), threshold=0.4, min_coverage=0.5, null_model=False, query_split="no_existe"
        )
        assert "error" in audit
        assert audit["query_split"] == "no_existe"


class TestPositiveRateDrift:
    """P-26: hay un caso que pasa TODAS las demas validaciones con train al 6%
    de positivos y test al 29%. PR-AUC, F1 y precision dependen de la
    prevalencia, y el umbral elegido en val queda descalibrado para test."""

    def _splits(self, rates):
        return {k: {"n": 100, "positive_rate": r} for k, r in rates.items()}

    def test_flags_the_known_bad_case(self):
        drift = _positive_rate_drift(
            self._splits({"train": 0.06, "val": 0.10, "test": 0.29}),
            {"positive_rate_tolerance": 0.05},
        )
        assert drift["exceeds_tolerance"]
        assert drift["argmin"] == "train"
        assert drift["argmax"] == "test"
        assert drift["spread"] == pytest.approx(0.23)

    def test_accepts_a_balanced_split(self):
        drift = _positive_rate_drift(
            self._splits({"train": 0.57, "val": 0.56, "test": 0.58}),
            {"positive_rate_tolerance": 0.05},
        )
        assert not drift["exceeds_tolerance"]

    def test_ignores_empty_splits(self):
        splits = self._splits({"train": 0.57, "test": 0.58})
        splits["val"] = {"n": 0, "positive_rate": 0.0}
        drift = _positive_rate_drift(splits, {"positive_rate_tolerance": 0.05})
        assert "val" not in drift["rates"], "un split vacio tiene tasa 0.0 y falsearia la deriva"
        assert not drift["exceeds_tolerance"]


class TestPuertasDePublicacion:
    """P-29 y D-22: que impide publicar y que hay que declarar al publicar.

    Historia de este campo, que explica por que esta separado en dos listas:

    1. `publishable` solo miraba la herramienta de clustering, y en la primera
       corrida real dijo `true` con fuga significativa en test Y en val (P-29);
    2. se arreglo metiendo la fuga en la misma lista, y entonces la puerta no se
       abria nunca: con una fuga real, medida y declarada, el proyecto no
       producia ningun numero. Eso no es honestidad, es paralisis, y empuja a
       retocar parametros hasta que la puerta ceda (D-22).

    Ahora la fuga es **salvedad** (se publica declarandola, junto a M-01) y lo
    que no tiene interpretacion posible es **bloqueo**.
    """

    @staticmethod
    def _audit(query, significant):
        return {
            "query_split": query,
            "reference_split": "train",
            "leakage": {
                "significant": significant,
                "excess": 50.19,
                "n_expected_by_chance": 23.81,
                "global_p_value": 0.0,
            },
        }

    def _report(self, test_sig, val_sig):
        return {
            "train_test_identity_audit": self._audit("test", test_sig),
            "train_val_identity_audit": self._audit("val", val_sig),
        }

    def test_la_fuga_es_salvedad_y_no_bloqueo(self):
        bloqueos, salvedades = _publication_gates(self._report(True, True), "mmseqs2")
        assert bloqueos == [], "una fuga medida y declarada no impide publicar (D-22)"
        assert len(salvedades) == 2

    def test_la_salvedad_lleva_las_cifras_y_la_via_para_cerrarla(self):
        _b, salvedades = _publication_gates(self._report(True, False), "mmseqs2")
        texto = salvedades[0]
        assert "test-vs-train" in texto
        assert "+50.2" in texto, "sin el exceso la salvedad no es cuantitativa"
        assert "M-01" in texto, "declarar un defecto sin decir como se cierra es media salvedad"

    def test_fuga_solo_en_val_tambien_se_declara(self):
        # val fija el umbral que se aplica a test, asi que la fuga val-train
        # sesga las metricas de test aunque test este limpio (P-27).
        _b, salvedades = _publication_gates(self._report(False, True), "mmseqs2")
        assert len(salvedades) == 1
        assert "val-vs-train" in salvedades[0]

    def test_sin_fuga_y_con_mmseqs2_no_hay_nada_que_declarar(self):
        bloqueos, salvedades = _publication_gates(self._report(False, False), "mmseqs2")
        assert bloqueos == []
        assert salvedades == []

    def test_el_respaldo_propio_si_bloquea(self):
        # P-23: su umbral efectivo sigue sin medirse, asi que el numero no tiene
        # interpretacion. No se puede declarar lo que no se conoce.
        bloqueos, _s = _publication_gates(self._report(False, False), "builtin")
        assert len(bloqueos) == 1
        assert "respaldo propio" in bloqueos[0]

    def test_una_auditoria_que_fallo_bloquea(self):
        """Cambio de comportamiento deliberado respecto a la version anterior.

        Antes, una auditoria con error **no aportaba motivo**, asi que
        `publishable` salia `true` sin haberse comprobado la fuga: P-29 otra vez
        por otra puerta. No se puede declarar lo que no se midio.
        """
        report = {
            "train_test_identity_audit": {"error": "falta el split test"},
            "train_val_identity_audit": self._audit("val", False),
        }
        bloqueos, _s = _publication_gates(report, "mmseqs2")
        assert len(bloqueos) == 1
        assert "fallo" in bloqueos[0]

    def test_una_auditoria_que_falta_bloquea(self):
        bloqueos, _s = _publication_gates({}, "mmseqs2")
        assert len(bloqueos) == 2, "faltan las dos auditorias, y ninguna se da por limpia"


class TestClusterModeExplicito:
    """P-29: el modo por defecto de MMseqs2 no da la garantia que el proyecto
    declara.

    El defecto resulto ser el 2 (incremental voraz), que garantiza que cada
    miembro se parece a SU representante y nada entre miembros de clusteres
    distintos. Sobre los datos reales dejo DLPSGWSSYEGH en train y
    DCPSDWSPYEGHCYK en test con 75% de identidad.

    Este test fija la decision: es una comprobacion de la orden construida, no
    del comportamiento de MMseqs2, que no se puede exigir en CI sin instalarlo.
    """

    def test_la_orden_pide_componente_conexa(self, monkeypatch, tmp_path):
        capturado = {}

        def falso_run(cmd, *args, **kwargs):
            capturado["cmd"] = cmd
            raise RuntimeError("corta aqui: solo interesa la orden construida")

        # Sin esto el test pasa en mi Mac y falla en el CI (P-34): la funcion
        # devuelve None antes de llamar a subprocess si no encuentra el binario,
        # asi que donde no hay MMseqs2 instalado no se construye orden alguna.
        # Se finge que esta, porque lo que se comprueba es la orden, no la
        # herramienta.
        monkeypatch.setattr("ampesm.split.shutil.which", lambda _: "/usr/bin/mmseqs")
        monkeypatch.setattr("ampesm.split.subprocess.run", falso_run)
        sequences = {"a": "ACDEFGHIKL", "b": "ACDEFGHIKM"}
        with pytest.raises(RuntimeError):
            cluster_mmseqs2(sequences, identity=0.4, coverage=0.5)

        cmd = capturado["cmd"]
        assert "easy-cluster" in cmd, "se capturo otra orden, no la del agrupamiento"
        assert "--cluster-mode" in cmd, (
            "sin --cluster-mode explicito se usa el defecto, que no garantiza "
            "que nada por encima del umbral cruce la frontera (P-29)"
        )
        assert cmd[cmd.index("--cluster-mode") + 1] == "1"


class TestFusionPorContencion:
    """P-30: el 8.5% del conjunto son subcadenas exactas de otras secuencias.

    Los casos de abajo son reales, sacados de la auditoria de la primera corrida
    con --cluster-mode 1, y MMseqs2 no agrupa ninguno (comprobado preguntandole
    pareja por pareja con todo permisivo: encuentra 3 de 8).
    """

    CASOS_REALES = [
        ("SLDCVLRP", "SKHSSLDCVLRP"),
        ("LFYLVPGPGH", "QLLIRMIYKNILFYLVPGPGHGAEPERRNIKYL"),
        ("ILGKIIKVVK", "NVWKKILGKIIKVVK"),
        ("PFKLSLHL", "TPFKLSLHL"),
        ("KKKKVVAATYVLF", "KKKKVVAATYVLFA"),
    ]

    def test_une_los_casos_reales_que_mmseqs2_no_ve(self):
        sequences, mapping = {}, {}
        for i, (corta, larga) in enumerate(self.CASOS_REALES):
            sequences[f"corta_{i}"] = corta
            sequences[f"larga_{i}"] = larga
            mapping[f"corta_{i}"] = 2 * i
            mapping[f"larga_{i}"] = 2 * i + 1

        merged, report = merge_contained_clusters(sequences, mapping)

        for i in range(len(self.CASOS_REALES)):
            assert (
                merged[f"corta_{i}"] == merged[f"larga_{i}"]
            ), f"el caso {self.CASOS_REALES[i]} sigue partido en dos clusteres"
        assert report["n_contained_sequences"] == len(self.CASOS_REALES)
        assert report["clusters_after"] == len(self.CASOS_REALES)

    def test_no_une_secuencias_sin_relacion(self):
        sequences = {"a": "KWKLFKKIGAVLKVL", "b": "ACDEFGHIKLMNPQRS", "c": "MSTYWVQNDEHGPLK"}
        mapping = {"a": 0, "b": 1, "c": 2}
        merged, report = merge_contained_clusters(sequences, mapping)
        assert merged == mapping
        assert report["n_contained_sequences"] == 0
        assert report["n_cluster_merges"] == 0

    def test_arrastra_a_los_companeros_de_cluster(self):
        # El fallo opuesto al que arregla la funcion: si al fusionar solo se
        # mueven las dos secuencias emparentadas, una tercera que compartia
        # cluster con el fragmento se queda atras y el cluster se PARTE, lo que
        # crearia fuga nueva en vez de quitarla.
        sequences = {
            "frag": "SLDCVLRP",
            "cont": "SKHSSLDCVLRP",
            "vecina": "KWKLFKKIGAVLKVL",
        }
        mapping = {"frag": 7, "vecina": 7, "cont": 3}
        merged, _ = merge_contained_clusters(sequences, mapping)
        assert merged["frag"] == merged["cont"] == merged["vecina"], (
            "la vecina del fragmento quedo en otro cluster: la fusion partio un "
            "cluster que ya estaba bien"
        )

    def test_une_con_todos_los_contenedores_no_solo_el_primero(self):
        # P-31: 803 de las 1 582 secuencias contenidas del conjunto real estan
        # dentro de mas de un contenedor. Cortar en el primero dejaba el resto
        # sin unir y, al depender del orden de iteracion de un set, no era ni
        # reproducible. Aqui el fragmento esta en dos contenedores que no se
        # contienen entre si: los tres clusteres tienen que acabar en uno.
        sequences = {
            "frag": "MPFSTGKRIMLGE",
            "cont_a": "AAVALLPAVLLALLAPMPFSTGKRIMLGE",
            "cont_b": "MPFSTGKRIMLGEWWWWDDDD",
        }
        mapping = {"frag": 10, "cont_a": 20, "cont_b": 30}
        merged, report = merge_contained_clusters(sequences, mapping)
        assert (
            len(set(merged.values())) == 1
        ), "el fragmento solo se unio a uno de sus dos contenedores"
        assert report["n_contained_sequences"] == 1
        assert report["n_containment_pairs"] == 2

    def test_dos_fusiones_que_se_solapan_no_se_sobrescriben(self):
        # P-31b: la primera version hacia union-find sobre SECUENCIAS y luego
        # reetiquetaba grupo a grupo leyendo del mapeo original. Cuando dos
        # grupos compartian un cluster, el segundo deshacia la fusion del
        # primero. Sintoma real: RPPGFSPFR quedo en test y ERPPGFSPFR, que lo
        # contiene, en train.
        #
        # El caso de abajo es el minimo que distingue las dos versiones:
        # cont_b comparte cluster (5) con frag_a, cuyo contenedor esta en el 2.
        # La version mala dejaba frag_a en 5 y cont_a en 2, separando el
        # fragmento de su propio contenedor. Verificado reproduciendo la
        # version antigua: daba 2 clusteres donde debe haber 1.
        sequences = {
            "frag_a": "SLDCVLRP",
            "cont_a": "SKHSSLDCVLRP",
            "frag_b": "RPPGFSPFR",
            "cont_b": "ERPPGFSPFR",
        }
        mapping = {"frag_a": 5, "cont_a": 2, "frag_b": 9, "cont_b": 5}
        merged, _ = merge_contained_clusters(sequences, mapping)
        assert (
            merged["frag_a"] == merged["cont_a"]
        ), "el fragmento quedo separado de su propio contenedor"
        assert (
            len(set(merged.values())) == 1
        ), f"una fusion sobrescribio a la otra: {sorted(set(merged.values()))}"
        assert set(merged.values()) == {2}, "el cluster resultante debe ser el id menor"

    def test_el_resultado_no_depende_del_orden_de_entrada(self):
        # El id del cluster resultante tiene que ser el menor de los fusionados,
        # no el primero que se encuentre: si no, dos corridas con los mismos
        # datos dan informes distintos (P-31).
        casos = {
            "a": "SLDCVLRP",
            "b": "SKHSSLDCVLRP",
            "c": "MMSKHSSLDCVLRPWW",
        }
        m1 = {"a": 9, "b": 4, "c": 7}
        merged1, _ = merge_contained_clusters(casos, m1)
        merged2, _ = merge_contained_clusters(
            dict(reversed(list(casos.items()))), dict(reversed(list(m1.items())))
        )
        assert set(merged1.values()) == {4}
        assert merged1 == merged2

    def test_la_contencion_es_transitiva(self):
        # a dentro de b, b dentro de c: los tres al mismo cluster.
        sequences = {"a": "SLDCVLRP", "b": "SKHSSLDCVLRP", "c": "MMSKHSSLDCVLRPWW"}
        mapping = {"a": 0, "b": 1, "c": 2}
        merged, report = merge_contained_clusters(sequences, mapping)
        assert len(set(merged.values())) == 1
        assert report["clusters_after"] == 1

    def test_secuencias_mas_cortas_que_el_prefiltro_no_rompen_nada(self):
        sequences = {"corta": "AC", "larga": "ACDEFGH"}
        mapping = {"corta": 0, "larga": 1}
        merged, report = merge_contained_clusters(sequences, mapping)
        # No se afirma que las encuentre: se afirma que no revienta. Con
        # min_length 5 en params no llegan aqui.
        assert merged["larga"] == 1
        assert report["n_contained_sequences"] == 0

    def test_el_informe_cuenta_lo_que_dice_contar(self):
        sequences, mapping = {}, {}
        for i, (corta, larga) in enumerate(self.CASOS_REALES[:3]):
            sequences[f"c{i}"] = corta
            sequences[f"l{i}"] = larga
            mapping[f"c{i}"] = 2 * i
            mapping[f"l{i}"] = 2 * i + 1
        _, report = merge_contained_clusters(sequences, mapping)
        assert report["clusters_before"] == 6
        assert report["clusters_after"] == 3
        assert report["n_cluster_merges"] == 3
        assert report["fraction_contained"] == 0.5
        assert report["n_containment_pairs"] == 3
        assert len(report["examples"]) == 3


@pytest.mark.skipif(importlib.util.find_spec("Bio") is None, reason="Biopython no instalado")
class TestAuditoriaExigeCobertura:
    """P-32: la auditoria no puede marcar parejas que el agrupamiento no separa.

    Las dos parejas de abajo son reales, sacadas de las secuencias que marcaba
    la auditoria de val en la cuarta corrida. Un alineador global con huecos
    encuentra cinco residuos dispersos dentro de casi cualquier peptido y
    devuelve identidad 1.000; eso no es homologia, es aritmetica.

    Medido, para que se vea por que el arreglo que P-32 tenia escrito (exigir
    cobertura sobre la secuencia CORTA) no servia:

        30 vs 5  identidad 1.000  /corta 1.000  /larga 0.167
        29 vs 29 identidad 1.000  /corta 1.000  /larga 1.000
    """

    FALSOS_POSITIVOS = [
        ("GIPCGESCVWIPCISGAIGCSCKSKVCYRN", "CCISG"),
        ("ATCDLFSFRSKWVTPNHAGCAAHCIFLGNRG", "VTPFL"),
    ]

    HOMOLOGAS = [
        ("KWKLFKKIGIGKFLHSAKKF", "KWKLFKKIGIGAFLHSAKKF"),
        ("GIGKFLHSAKKFGKAFVGEIMNS", "GIGKFLHSAKKWGKAFVGEIMNS"),
    ]

    @staticmethod
    def _df(largo: str, corto: str) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {"id": "tr0", "sequence": largo, "split": "train", "label": 1, "cluster_id": 0},
                {"id": "te0", "sequence": corto, "split": "test", "label": 1, "cluster_id": 1},
            ]
        )

    @pytest.mark.parametrize(("largo", "corto"), FALSOS_POSITIVOS)
    def test_no_marca_lo_que_el_agrupamiento_no_separa(self, largo, corto):
        audit = audit_train_test_identity(
            self._df(largo, corto), threshold=0.4, min_coverage=0.5, null_model=False
        )
        assert audit["n_over_threshold"] == 0, (
            f"{corto} ({len(corto)} residuos) contra {len(largo)} da cobertura "
            f"{len(corto) / len(largo):.2f} sobre la larga: el agrupamiento con "
            f"-c 0.5 --cov-mode 1 nunca tuvo que separarlos"
        )

    @pytest.mark.parametrize(("a", "b"), HOMOLOGAS)
    def test_sigue_marcando_las_homologas_de_verdad(self, a, b):
        audit = audit_train_test_identity(
            self._df(a, b), threshold=0.4, min_coverage=0.5, null_model=False
        )
        assert audit["n_over_threshold"] == 1, (
            "el filtro de cobertura no puede cegar la auditoria ante parejas de "
            "longitud parecida, que son las que si son homologia"
        )

    def test_marca_la_subcadena_exacta_aunque_no_llegue_a_la_cobertura(self):
        """El segundo termino de la O: la garantia de D-20, no la de MMseqs2.

        Si la auditoria filtrara solo por cobertura, una subcadena exacta corta
        en un split distinto dejaria de marcarse, y eso es exactamente el fallo
        que la fusion por contencion puede tener (P-31, P-31b fueron dos). La
        guarda tiene que poder ver un fallo de D-20, no solo uno de MMseqs2.
        """
        largo = "GIGKFLHSAKKFGKAFVGEIMNSKKWWRRPLFVACDEGHIKLMNPQRSTV"
        audit = audit_train_test_identity(
            self._df(largo, largo[10:18]), threshold=0.4, min_coverage=0.5, null_model=False
        )
        assert audit["n_over_threshold"] == 1, (
            f"cobertura {8 / len(largo):.2f}, por debajo de 0.5, pero es subcadena "
            "exacta: tiene que marcarse por el segundo termino"
        )

    def test_el_informe_deja_escrita_la_definicion_usada(self):
        audit = audit_train_test_identity(
            self._df(*self.HOMOLOGAS[0]), threshold=0.4, min_coverage=0.5, null_model=False
        )
        assert audit["min_coverage"] == 0.5, (
            "el p-valor de fuga solo se puede leer sabiendo que definicion de " "homologo se audito"
        )

    def test_min_coverage_no_tiene_valor_por_defecto(self):
        """P-29 otra vez: un defecto se comporta como una decision del proyecto.

        Quien audite declara que definicion de homologo usa, o no compila.
        """
        import inspect

        firma = inspect.signature(audit_train_test_identity)
        assert firma.parameters["min_coverage"].default is inspect.Parameter.empty


class TestIdsDeClusterDeterministas:
    """P-36: MMseqs2 escribe las filas en un orden que cambia entre corridas.

    Comprobado corriendolo dos veces sobre el mismo FASTA de 18 644 secuencias:
    los dos TSV difieren byte a byte y, ordenados, coinciden en las 18 644 filas
    sin una sola diferencia. Misma particion, orden de filas distinto.

    Asignar los ids por orden de aparicion convertia eso en divisiones
    distintas, porque `assign_splits` baraja con semilla fija partiendo del
    orden de los ids. Medido pasando los dos TSV reales por la cadena entera:
    **6 531 de 18 644 secuencias (35%) cambiaban de split** entre dos corridas
    identicas. Con ids canonicos, cero.
    """

    PARES = [
        ("rep_b", "seq_1"),
        ("rep_b", "seq_2"),
        ("rep_a", "seq_3"),
        ("rep_c", "seq_4"),
    ]

    def test_el_orden_de_las_filas_no_cambia_ni_una_etiqueta(self):
        directo = cluster_ids_from_pairs(self.PARES)
        invertido = cluster_ids_from_pairs(list(reversed(self.PARES)))
        assert directo == invertido, (
            "mismas parejas en otro orden dan otras etiquetas: eso basta para "
            "que la division entera cambie (P-36)"
        )

    def test_el_id_sale_del_representante_ordenado(self):
        ids = cluster_ids_from_pairs(self.PARES)
        # rep_a < rep_b < rep_c alfabeticamente, asi que 0, 1, 2.
        assert ids["seq_3"] == 0
        assert ids["seq_1"] == ids["seq_2"] == 1
        assert ids["seq_4"] == 2

    def test_todas_las_permutaciones_dan_lo_mismo(self):
        import itertools

        esperado = cluster_ids_from_pairs(self.PARES)
        for perm in itertools.permutations(self.PARES):
            assert cluster_ids_from_pairs(list(perm)) == esperado

    def test_la_particion_se_conserva(self):
        """Lo que no debe cambiar el arreglo: quien va con quien."""
        ids = cluster_ids_from_pairs(self.PARES)
        assert ids["seq_1"] == ids["seq_2"], "seq_1 y seq_2 comparten representante"
        assert len(set(ids.values())) == 3
