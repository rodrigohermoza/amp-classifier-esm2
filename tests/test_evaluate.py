"""Tests de metricas, umbral y bootstrap."""

import json

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

from ampesm.config import Params
from ampesm.evaluate import (
    PREREGISTERED_BASELINE,
    PREREGISTERED_EMBEDDING,
    _select_comparison_pair,
    all_metrics,
    bootstrap_metrics,
    choose_threshold,
    paired_difference,
    threshold_metrics,
)


def make_predictions(n_clusters=40, per_cluster=5, separation=1.5, seed=0):
    """Predicciones sinteticas con estructura de clusteres.

    Las secuencias de un cluster comparten etiqueta y reciben puntuaciones
    correlacionadas, que es como se comportan los datos reales: los parientes
    cercanos se predicen parecido.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for cid in range(n_clusters):
        label = cid % 4 == 0  # ~25% positivos
        center = rng.normal(separation if label else -separation, 0.5)
        for j in range(per_cluster):
            rows.append(
                {
                    "id": f"s_{cid}_{j}",
                    "cluster_id": cid,
                    "split": "test",
                    "label": int(label),
                    "score": 1 / (1 + np.exp(-(center + rng.normal(0, 0.3)))),
                }
            )
    return pd.DataFrame(rows)


class TestMetrics:
    def test_perfect_separation(self):
        y = np.array([0, 0, 1, 1])
        scores = np.array([0.1, 0.2, 0.8, 0.9])
        metrics = all_metrics(y, scores, threshold=0.5)
        assert metrics["roc_auc"] == 1.0
        assert metrics["mcc"] == 1.0

    def test_random_scores_give_chance_auc(self):
        rng = np.random.default_rng(0)
        y = rng.integers(0, 2, 2000)
        scores = rng.random(2000)
        assert abs(all_metrics(y, scores, 0.5)["roc_auc"] - 0.5) < 0.05

    def test_inverted_scores_give_low_auc(self):
        y = np.array([0, 0, 1, 1])
        scores = np.array([0.9, 0.8, 0.2, 0.1])
        assert all_metrics(y, scores, 0.5)["roc_auc"] == 0.0

    def test_single_class_gives_nan_for_every_metric(self):
        """Regresion de P-07.

        Antes solo roc_auc daba NaN, porque es la unica metrica que lanza
        excepcion. PR-AUC devolvia 1.0, MCC devolvia 0.0 y los guardas daban
        0.0, y todos esos valores entraban en los intervalos de confianza de
        las dos metricas principales del proyecto.
        """
        for y in (np.array([1, 1, 1]), np.array([0, 0, 0])):
            metrics = all_metrics(y, np.array([0.2, 0.5, 0.9]), 0.5)
            assert all(
                np.isnan(v) for v in metrics.values()
            ), f"con y={y.tolist()} deberian ser todas NaN, salieron {metrics}"

    def test_empty_input_gives_nan_not_indexerror(self):
        metrics = all_metrics(np.array([]), np.array([]), 0.5)
        assert all(np.isnan(v) for v in metrics.values())

    def test_sensitivity_and_specificity(self):
        y = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        pred = np.array([0, 0, 0, 1, 0, 1, 1, 1])
        metrics = threshold_metrics(y, pred)
        assert metrics["sensitivity"] == 0.75  # 3 de 4 positivos
        assert metrics["specificity"] == 0.75  # 3 de 4 negativos


class TestThreshold:
    def test_finds_separating_threshold(self):
        y = np.array([0, 0, 0, 1, 1, 1])
        scores = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])
        threshold, value = choose_threshold(y, scores, "mcc")
        assert value == 1.0
        assert 0.3 < threshold <= 0.7

    def test_returns_valid_threshold_on_noise(self):
        rng = np.random.default_rng(0)
        y = rng.integers(0, 2, 200)
        threshold, value = choose_threshold(y, rng.random(200), "mcc")
        assert 0.0 <= threshold <= 1.0
        assert -1.0 <= value <= 1.0

    def test_handles_imbalance(self):
        """Con 5% de positivos, el umbral optimo no es 0.5."""
        y = np.array([0] * 95 + [1] * 5)
        scores = np.concatenate([np.linspace(0.0, 0.5, 95), np.linspace(0.55, 0.7, 5)])
        threshold, value = choose_threshold(y, scores, "mcc")
        assert value > 0.9
        assert threshold < 0.6


class TestBootstrap:
    def test_point_estimate_inside_interval(self):
        df = make_predictions()
        result = bootstrap_metrics(df, threshold=0.5, n_boot=200, unit="cluster", seed=42)
        for name, values in result.items():
            if not np.isnan(values["point"]):
                assert values["ci_low"] <= values["point"] <= values["ci_high"], name

    def test_cluster_bootstrap_is_wider_than_sequence(self):
        """La justificacion de D-08, comprobada en vez de afirmada.

        Las secuencias de un cluster no son independientes. Remuestrear
        secuencias individuales finge que si lo son y produce intervalos
        demasiado estrechos. Remuestrear clusteres completos respeta la
        dependencia y da intervalos mas anchos, que son los honestos.
        """
        # separation baja a proposito: con clases perfectamente separables el
        # AUC es 1.0 en todo remuestreo y ambos intervalos tendrian ancho cero,
        # con lo que el test no mediria nada.
        df = make_predictions(n_clusters=30, per_cluster=10, separation=0.35, seed=5)
        assert 0.6 < roc_auc_score(df.label, df.score) < 0.98, "los datos deben solaparse"
        by_cluster = bootstrap_metrics(df, 0.5, n_boot=300, unit="cluster", seed=42)
        by_sequence = bootstrap_metrics(df, 0.5, n_boot=300, unit="sequence", seed=42)

        width_cluster = by_cluster["roc_auc"]["ci_high"] - by_cluster["roc_auc"]["ci_low"]
        width_sequence = by_sequence["roc_auc"]["ci_high"] - by_sequence["roc_auc"]["ci_low"]
        assert width_cluster > width_sequence, (
            f"el IC por cluster ({width_cluster:.4f}) deberia ser mas ancho que "
            f"el de secuencia ({width_sequence:.4f})"
        )

    def test_deterministic_given_seed(self):
        df = make_predictions()
        a = bootstrap_metrics(df, 0.5, n_boot=100, seed=7)
        b = bootstrap_metrics(df, 0.5, n_boot=100, seed=7)
        assert a["roc_auc"] == b["roc_auc"]

    def test_more_data_narrows_interval(self):
        small = make_predictions(n_clusters=20, per_cluster=5, seed=1)
        large = make_predictions(n_clusters=200, per_cluster=5, seed=1)
        w_small = (lambda r: r["roc_auc"]["ci_high"] - r["roc_auc"]["ci_low"])(
            bootstrap_metrics(small, 0.5, n_boot=200, seed=42)
        )
        w_large = (lambda r: r["roc_auc"]["ci_high"] - r["roc_auc"]["ci_low"])(
            bootstrap_metrics(large, 0.5, n_boot=200, seed=42)
        )
        assert w_large < w_small


class TestPairedDifference:
    def test_identical_models_differ_by_zero(self):
        df = make_predictions()
        result = paired_difference(df, df.copy(), 0.5, 0.5, "mcc", n_boot=200, seed=42)
        assert result["difference"] == 0.0
        assert not result["significant"]
        assert result["ci_low"] <= 0 <= result["ci_high"]

    def test_clearly_better_model_is_detected(self):
        df_good = make_predictions(separation=2.5, seed=1)
        df_bad = df_good.copy()
        rng = np.random.default_rng(2)
        df_bad["score"] = rng.random(len(df_bad))  # modelo al azar

        result = paired_difference(df_good, df_bad, 0.5, 0.5, "mcc", n_boot=300, seed=42)
        assert result["difference"] > 0
        assert result["significant"]
        assert result["ci_low"] > 0

    def test_misaligned_predictions_raise(self):
        df = make_predictions()
        try:
            paired_difference(df, df.iloc[:10], 0.5, 0.5, n_boot=10)
        except ValueError as exc:
            assert "alinean" in str(exc)
        else:
            raise AssertionError("se esperaba ValueError al cruzar predicciones desalineadas")


class TestThresholdGuards:
    """Regresiones de P-10d: umbrales degenerados en silencio."""

    def test_empty_validation_raises(self):
        with pytest.raises(ValueError, match="vacio"):
            choose_threshold(np.array([]), np.array([]), "mcc")

    def test_single_class_validation_raises(self):
        """Antes devolvia el score MAS BAJO de val, que en test predice casi
        todo positivo, sin que nada avisara."""
        with pytest.raises(ValueError, match="una sola clase"):
            choose_threshold(np.array([1, 1, 1]), np.array([0.9, 0.8, 0.7]), "mcc")


class TestBootstrapBookkeeping:
    """P-07: descartar remuestreos condiciona el IC, asi que hay que contarlos."""

    def test_reports_valid_bootstrap_count(self):
        df = make_predictions(separation=0.4, seed=3)
        result = bootstrap_metrics(df, 0.5, n_boot=100, unit="cluster", seed=42)
        for name, values in result.items():
            assert values["n_boot_requested"] == 100, name
            assert 0 <= values["n_boot_valid"] <= 100, name


class TestComparisonPairSelection:
    """Regresiones de P-05: el par se elegia por su resultado en test."""

    def _summary(self, tags):
        return pd.DataFrame({"tag": tags, "mcc": np.linspace(0.9, 0.5, len(tags))})

    def test_preregistered_pair_wins_when_present(self):
        tags = [
            PREREGISTERED_EMBEDDING,
            "facebook_esm2_t6_8M_UR50D__mlp",
            PREREGISTERED_BASELINE,
            "composition__logreg",
        ]
        # Val deliberadamente favorable a los NO preregistrados: el preregistro
        # debe ganar de todas formas.
        val = {
            PREREGISTERED_EMBEDDING: 0.70,
            "facebook_esm2_t6_8M_UR50D__mlp": 0.99,
            PREREGISTERED_BASELINE: 0.70,
            "composition__logreg": 0.99,
        }
        emb, base, basis = _select_comparison_pair(self._summary(tags), val)
        assert emb == PREREGISTERED_EMBEDDING
        assert base == PREREGISTERED_BASELINE
        assert "preregistro" in basis

    def test_falls_back_to_validation_not_test(self):
        """Sin el par preregistrado, se elige por val y NUNCA por el orden de
        summary, que esta ordenado por MCC de test."""
        tags = ["facebook_a__logreg", "facebook_b__mlp", "composition__logreg", "physchem__logreg"]
        summary = self._summary(tags)  # mcc de test: a > b > composition > physchem
        val = {
            "facebook_a__logreg": 0.60,
            "facebook_b__mlp": 0.95,  # peor en test, mejor en val
            "composition__logreg": 0.60,
            "physchem__logreg": 0.95,  # peor en test, mejor en val
        }
        emb, base, basis = _select_comparison_pair(summary, val)
        assert emb == "facebook_b__mlp"
        assert base == "physchem__logreg"
        assert "validacion" in basis

    def test_missing_family_returns_none(self):
        emb, base, basis = _select_comparison_pair(
            self._summary(["composition__logreg"]), {"composition__logreg": 0.8}
        )
        assert emb is None and base is None
        assert "no aplicable" in basis


class TestFeaturesHash:
    """P-24: cambiar `features` y correr solo `interpret` no avisaba.

    `error_analysis` recalcula los descriptores con la configuracion actual, asi
    que con otra `moment_window` el analisis describe un modelo que no existe.
    No fallaba nada: salia una figura plausible y equivocada.
    """

    def test_hash_changes_with_the_configuration(self):
        a = Params(raw={"features": {"physchem_ph": 7.4, "moment_window": 11}})
        b = Params(raw={"features": {"physchem_ph": 7.4, "moment_window": 18}})
        assert a.features_hash() != b.features_hash()

    def test_hash_does_not_depend_on_key_order(self):
        a = Params(raw={"features": {"physchem_ph": 7.4, "moment_window": 11}})
        b = Params(raw={"features": {"moment_window": 11, "physchem_ph": 7.4}})
        assert a.features_hash() == b.features_hash()

    def test_interpret_aborts_when_features_changed(self, tmp_path, monkeypatch):
        import ampesm.interpret as interpret_mod

        monkeypatch.setattr(interpret_mod, "MODELS_DIR", tmp_path)
        (tmp_path / "modelo.json").write_text(
            json.dumps({"tag": "modelo", "features_hash": "0000000000000000"})
        )

        with pytest.raises(ValueError, match="features"):
            interpret_mod._check_features_hash(
                "modelo", Params(raw={"features": {"moment_window": 11}})
            )

    def test_interpret_accepts_a_matching_hash(self, tmp_path, monkeypatch):
        import ampesm.interpret as interpret_mod

        params = Params(raw={"features": {"moment_window": 11}})
        monkeypatch.setattr(interpret_mod, "MODELS_DIR", tmp_path)
        (tmp_path / "modelo.json").write_text(
            json.dumps({"tag": "modelo", "features_hash": params.features_hash()})
        )

        interpret_mod._check_features_hash("modelo", params)  # no debe lanzar

    def test_old_model_without_hash_warns_but_does_not_abort(self, tmp_path, monkeypatch):
        """Un modelo entrenado antes de este arreglo no tiene el campo. No se
        puede comprobar, pero romper la etapa por eso seria peor."""
        import ampesm.interpret as interpret_mod

        monkeypatch.setattr(interpret_mod, "MODELS_DIR", tmp_path)
        (tmp_path / "modelo.json").write_text(json.dumps({"tag": "modelo"}))

        interpret_mod._check_features_hash(
            "modelo", Params(raw={"features": {"moment_window": 11}})
        )
