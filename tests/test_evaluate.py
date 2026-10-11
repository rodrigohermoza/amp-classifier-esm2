"""Tests de metricas, umbral y bootstrap."""

import numpy as np
import pandas as pd

from ampesm.evaluate import (
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

    def test_single_class_gives_nan_auc_not_crash(self):
        """Pasa de verdad en un remuestreo bootstrap; no debe reventar."""
        metrics = all_metrics(np.array([1, 1, 1]), np.array([0.2, 0.5, 0.9]), 0.5)
        assert np.isnan(metrics["roc_auc"])

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
        df = make_predictions(n_clusters=30, per_cluster=10)
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
        w_small = (
            lambda r: r["roc_auc"]["ci_high"] - r["roc_auc"]["ci_low"]
        )(bootstrap_metrics(small, 0.5, n_boot=200, seed=42))
        w_large = (
            lambda r: r["roc_auc"]["ci_high"] - r["roc_auc"]["ci_low"]
        )(bootstrap_metrics(large, 0.5, n_boot=200, seed=42))
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
