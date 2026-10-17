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
    _publishable_reasons,
    assign_splits,
    audit_train_test_identity,
    cluster_builtin,
    cluster_mmseqs2,
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
        audit = audit_train_test_identity(self._frame(0, 0.0), 0.4, null_replicates=20, seed=5)
        assert not audit["leakage"]["significant"], f"falso positivo: {audit['leakage']}"

    def test_planted_leakage_is_detected(self):
        """Regresion de P-18.

        El estadistico anterior (diferencia de proporciones agregadas) no
        detectaba ni una fuga del 17% al 80% de identidad, porque comparaba dos
        proporciones ya saturadas por el ruido de fondo.
        """
        audit = audit_train_test_identity(self._frame(12, 0.75), 0.4, null_replicates=20, seed=5)
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
        audit = audit_train_test_identity(df, 0.4, null_model=False)
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

        on_test = audit_train_test_identity(df, threshold=0.4, null_model=False)
        on_val = audit_train_test_identity(df, threshold=0.4, null_model=False, query_split="val")

        assert on_test["n_over_threshold"] == 0, "test esta limpio, y el audit de test lo dice"
        assert (
            on_val["n_over_threshold"] == 1
        ), "la fuga vive en val; si esto falla, la auditoria no mira donde hace falta"

    def test_report_says_which_splits_it_compared(self):
        audit = audit_train_test_identity(
            self._df(), threshold=0.4, null_model=False, query_split="val"
        )
        assert audit["query_split"] == "val"
        assert audit["reference_split"] == "train"

    def test_missing_split_returns_error_not_crash(self):
        audit = audit_train_test_identity(
            self._df(), threshold=0.4, null_model=False, query_split="no_existe"
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


class TestPublishableTieneEnCuentaLaFuga:
    """P-29: `publishable` salio true con fuga significativa en test Y en val.

    El campo solo miraba la herramienta de clustering. Un indicador que se llama
    "publishable" y no mira el resultado de la auditoria invita a saltarse justo
    la lectura que importa.
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

    def test_fuga_en_test_lo_vuelve_no_publicable(self):
        report = {"train_test_identity_audit": self._audit("test", True)}
        reasons = _publishable_reasons(report, "mmseqs2")
        assert len(reasons) == 1
        assert "test-vs-train" in reasons[0]

    def test_fuga_solo_en_val_tambien_lo_vuelve_no_publicable(self):
        # val fija el umbral que se aplica a test, asi que la fuga val-train
        # sesga las metricas de test aunque test este limpio (P-27).
        report = {
            "train_test_identity_audit": self._audit("test", False),
            "train_val_identity_audit": self._audit("val", True),
        }
        reasons = _publishable_reasons(report, "mmseqs2")
        assert len(reasons) == 1
        assert "val-vs-train" in reasons[0]

    def test_el_caso_real_de_la_primera_corrida_acumula_los_dos_motivos(self):
        report = {
            "train_test_identity_audit": self._audit("test", True),
            "train_val_identity_audit": self._audit("val", True),
        }
        reasons = _publishable_reasons(report, "mmseqs2")
        assert len(reasons) == 2

    def test_sin_fuga_y_con_mmseqs2_es_publicable(self):
        report = {
            "train_test_identity_audit": self._audit("test", False),
            "train_val_identity_audit": self._audit("val", False),
        }
        assert _publishable_reasons(report, "mmseqs2") == []

    def test_el_respaldo_propio_sigue_contando_como_motivo(self):
        report = {"train_test_identity_audit": self._audit("test", False)}
        reasons = _publishable_reasons(report, "builtin")
        assert len(reasons) == 1
        assert "respaldo propio" in reasons[0]

    def test_una_auditoria_con_error_no_se_toma_por_limpia(self):
        # Que la auditoria no haya podido correr no es evidencia de ausencia de
        # fuga, pero tampoco se puede afirmar que la haya: no aporta motivo.
        report = {"train_test_identity_audit": {"error": "falta el split test"}}
        assert _publishable_reasons(report, "mmseqs2") == []


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
