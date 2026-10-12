"""Metricas, intervalos de confianza por bootstrap y figuras.

Trabaja solo sobre los CSV de predicciones que escribe models.py, asi que
recalcular metricas no obliga a reentrenar.

Dos decisiones metodologicas que estan implementadas aqui y justificadas en
docs/07-evaluacion.md:

1. El umbral se elige en validacion y se aplica a prueba sin tocarlo.
2. El bootstrap remuestrea clusteres completos, no secuencias. Las secuencias
   de un mismo cluster son parientes; tratarlas como independientes produce
   intervalos demasiado estrechos.
"""

from __future__ import annotations

import json
from collections.abc import Callable

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
)

from .config import FIGURES_DIR, METRICS_DIR, ensure_dirs, load_params
from .models import PREDICTIONS_DIR

# Metricas que dependen de un umbral y metricas que no. La distincion importa
# porque las primeras necesitan que el umbral venga de validacion.
RANKING_METRICS: dict[str, Callable] = {
    "roc_auc": roc_auc_score,
    "pr_auc": average_precision_score,
}


def threshold_metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    """Metricas que requieren predicciones binarias."""
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {
        "mcc": float(matthews_corrcoef(y, pred)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "sensitivity": float(tp / (tp + fn)) if (tp + fn) else 0.0,
        "specificity": float(tn / (tn + fp)) if (tn + fp) else 0.0,
        "precision": float(tp / (tp + fp)) if (tp + fp) else 0.0,
    }


METRIC_NAMES = ("roc_auc", "pr_auc", "mcc", "f1", "sensitivity", "specificity", "precision")


def all_metrics(y: np.ndarray, scores: np.ndarray, threshold: float) -> dict[str, float]:
    """Todas las metricas de una vez.

    Si el remuestreo deja una sola clase (o esta vacio), NINGUNA metrica esta
    definida y se devuelven todas como NaN, para que el percentil del bootstrap
    las ignore.

    P-07: antes solo se capturaba la excepcion de ROC-AUC, que es la unica
    metrica que lanza. Las demas devuelven valores finitos pero degenerados y
    entraban en el intervalo de confianza: con `y=[1,1,1]`,
    `average_precision_score` devuelve 1.0, `matthews_corrcoef` devuelve 0.0 por
    denominador nulo, y los guardas de `threshold_metrics` devuelven 0.0. Eso
    contaminaba los intervalos de PR-AUC y MCC, que son las dos metricas
    principales del proyecto, con unos y ceros duros. En `paired_difference` era
    peor: aportaba un cero exacto que contraia el intervalo de la diferencia
    hacia "no distinguible".
    """
    if len(y) == 0 or len(np.unique(y)) < 2:
        return dict.fromkeys(METRIC_NAMES, float("nan"))

    out: dict[str, float] = {}
    for name, fn in RANKING_METRICS.items():
        try:
            out[name] = float(fn(y, scores))
        except ValueError:
            out[name] = float("nan")
    out.update(threshold_metrics(y, (scores >= threshold).astype(int)))
    return out


def choose_threshold(
    y_val: np.ndarray, scores_val: np.ndarray, metric: str = "mcc"
) -> tuple[float, float]:
    """Umbral que maximiza una metrica en validacion.

    Se barren los propios valores de puntuacion como candidatos, en vez de una
    rejilla fija: el optimo siempre esta en uno de ellos, porque la metrica solo
    cambia cuando el umbral cruza una puntuacion observada.

    P-10d: con validacion vacia o monoclase, todos los candidatos empataban y se
    devolvia el score MAS BAJO de val como umbral, con lo que en test se
    predecia casi todo positivo. Y con val vacia se devolvia `-inf`, que
    `json.dumps` escribe como `-Infinity`, que no es JSON valido. Ahora falla
    explicitamente en vez de devolver un umbral sin sentido.
    """
    if len(y_val) == 0:
        raise ValueError("No se puede elegir el umbral: el conjunto de validacion esta vacio")
    if len(np.unique(y_val)) < 2:
        raise ValueError(
            "No se puede elegir el umbral: validacion tiene una sola clase. "
            "Revisa la division (docs/09, P-06)."
        )

    candidates = np.unique(scores_val)
    if len(candidates) > 2000:
        candidates = np.quantile(candidates, np.linspace(0, 1, 2000))

    best_threshold, best_value = 0.5, -np.inf
    for threshold in candidates:
        value = threshold_metrics(y_val, (scores_val >= threshold).astype(int))[metric]
        if value > best_value:
            best_threshold, best_value = float(threshold), value
    return best_threshold, float(best_value)


def bootstrap_metrics(
    df: pd.DataFrame,
    threshold: float,
    n_boot: int = 2000,
    unit: str = "cluster",
    ci: float = 0.95,
    seed: int = 42,
) -> dict[str, dict[str, float]]:
    """Intervalos de confianza por bootstrap.

    Con unit="cluster" se remuestrean clusteres completos con reemplazo, lo que
    respeta la dependencia entre secuencias emparentadas. Con unit="sequence" se
    remuestrean filas, que es lo habitual en la literatura y da intervalos
    optimistamente estrechos; esta disponible solo para poder mostrar la
    diferencia.
    """
    rng = np.random.default_rng(seed)
    y = df.label.to_numpy()
    scores = df.score.to_numpy()

    point = all_metrics(y, scores, threshold)
    samples: dict[str, list[float]] = {k: [] for k in point}

    if unit == "cluster":
        groups = df.groupby("cluster_id").indices
        cluster_ids = np.array(list(groups.keys()))
        index_lists = [groups[c] for c in cluster_ids]

        for _ in range(n_boot):
            picked = rng.integers(0, len(cluster_ids), size=len(cluster_ids))
            idx = np.concatenate([index_lists[i] for i in picked])
            for name, value in all_metrics(y[idx], scores[idx], threshold).items():
                samples[name].append(value)
    else:
        n = len(df)
        for _ in range(n_boot):
            idx = rng.integers(0, n, size=n)
            for name, value in all_metrics(y[idx], scores[idx], threshold).items():
                samples[name].append(value)

    alpha = (1.0 - ci) / 2.0
    out = {}
    for name, values in samples.items():
        arr = np.asarray(values, dtype=float)
        n_total = len(arr)
        arr = arr[~np.isnan(arr)]
        # Descartar remuestreos degenerados condiciona el intervalo, asi que se
        # registra cuantos sobrevivieron: sin este numero no queda rastro de
        # cuantos se fueron y el IC parece mas fiable de lo que es (P-07).
        out[name] = {
            "point": point[name],
            "ci_low": float(np.quantile(arr, alpha)) if len(arr) else float("nan"),
            "ci_high": float(np.quantile(arr, 1 - alpha)) if len(arr) else float("nan"),
            "n_boot_valid": int(len(arr)),
            "n_boot_requested": int(n_total),
        }
        if n_total and len(arr) / n_total < 0.99:
            print(
                f"[eval] AVISO: {n_total - len(arr)} de {n_total} remuestreos degenerados "
                f"descartados para {name}; el IC esta condicionado"
            )
    return out


def paired_difference(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    threshold_a: float,
    threshold_b: float,
    metric: str = "mcc",
    n_boot: int = 2000,
    unit: str = "cluster",
    ci: float = 0.95,
    seed: int = 42,
) -> dict:
    """Bootstrap de la diferencia A - B, con los mismos remuestreos para ambos.

    Esta es la forma correcta de comparar dos modelos. Mirar si los intervalos
    individuales se solapan es estadisticamente incorrecto y es un error
    frecuente: dos intervalos que se solapan pueden corresponder a una
    diferencia significativa, y al contrario.

    Si el intervalo de la diferencia incluye el cero, la conclusion es "no se
    distinguen", no "A es mejor".
    """
    merged = df_a.merge(df_b, on=["id", "cluster_id", "split", "label"], suffixes=("_a", "_b"))
    if len(merged) != len(df_a):
        raise ValueError(
            f"Las predicciones no se alinean: {len(df_a)} frente a {len(merged)} tras el cruce"
        )

    rng = np.random.default_rng(seed)
    y = merged.label.to_numpy()
    sa = merged.score_a.to_numpy()
    sb = merged.score_b.to_numpy()

    def value(idx):
        ma = all_metrics(y[idx], sa[idx], threshold_a)[metric]
        mb = all_metrics(y[idx], sb[idx], threshold_b)[metric]
        return ma - mb

    point = value(np.arange(len(merged)))
    diffs = []

    if unit == "cluster":
        groups = merged.groupby("cluster_id").indices
        index_lists = [groups[c] for c in groups]
        for _ in range(n_boot):
            picked = rng.integers(0, len(index_lists), size=len(index_lists))
            diffs.append(value(np.concatenate([index_lists[i] for i in picked])))
    else:
        for _ in range(n_boot):
            diffs.append(value(rng.integers(0, len(merged), size=len(merged))))

    arr = np.asarray(diffs, dtype=float)
    arr = arr[~np.isnan(arr)]
    alpha = (1.0 - ci) / 2.0
    low = float(np.quantile(arr, alpha))
    high = float(np.quantile(arr, 1 - alpha))

    return {
        "metric": metric,
        "difference": float(point),
        "ci_low": low,
        "ci_high": high,
        "significant": bool(low > 0 or high < 0),
        "interpretation": (
            "La diferencia es distinguible de cero."
            if (low > 0 or high < 0)
            else "El intervalo incluye el cero: los modelos no se distinguen."
        ),
    }


# --- Figuras ---------------------------------------------------------------


def _setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")  # sin ventana: la tuberia corre sin interfaz grafica
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 130,
            "savefig.bbox": "tight",
            "font.size": 10,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    return plt


def plot_curves(results: dict[str, pd.DataFrame]) -> None:
    """Curvas ROC y precision-recall de todos los modelos."""
    from sklearn.metrics import precision_recall_curve, roc_curve

    plt = _setup_matplotlib()

    for kind in ("roc", "pr"):
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        for tag, df in sorted(results.items()):
            y, scores = df.label.to_numpy(), df.score.to_numpy()
            # Si test quedara monoclase, roc_curve lanza ValueError. Aqui ya
            # estan escritas todas las metricas, asi que reventar al dibujar
            # perderia la etapa entera por una figura: se omite la curva y se
            # avisa.
            try:
                if kind == "roc":
                    x_vals, y_vals, _ = roc_curve(y, scores)
                    score = roc_auc_score(y, scores)
                else:
                    y_vals, x_vals, _ = precision_recall_curve(y, scores)
                    score = average_precision_score(y, scores)
            except ValueError as exc:
                print(f"[eval] AVISO: no se pudo dibujar la curva {kind} de {tag}: {exc}")
                continue
            ax.plot(x_vals, y_vals, lw=1.6, label=f"{tag} ({score:.3f})")

        if kind == "roc":
            ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5, label="azar")
            ax.set_xlabel("Tasa de falsos positivos (1 - especificidad)")
            ax.set_ylabel("Tasa de verdaderos positivos (sensibilidad)")
            ax.set_title("Curvas ROC en el conjunto de prueba")
        else:
            ax.set_xlabel("Sensibilidad (recall)")
            ax.set_ylabel("Precision")
            ax.set_title("Curvas precision-recall en el conjunto de prueba")

        ax.legend(fontsize=7, loc="lower right" if kind == "roc" else "lower left")
        fig.savefig(FIGURES_DIR / f"{kind}_curves.png")
        plt.close(fig)


def plot_metrics_comparison(summary: pd.DataFrame) -> None:
    """Barras con barras de error: el grafico que resume el proyecto."""
    plt = _setup_matplotlib()
    metrics = ["roc_auc", "pr_auc", "mcc", "f1"]
    tags = summary.tag.tolist()

    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 5), sharey=False)
    for ax, metric in zip(axes, metrics):
        point = summary[f"{metric}"].to_numpy()
        low = summary[f"{metric}_ci_low"].to_numpy()
        high = summary[f"{metric}_ci_high"].to_numpy()
        # Las barras de error se dan como distancia al punto, no como limites.
        err = np.vstack([point - low, high - point])
        ax.barh(range(len(tags)), point, xerr=err, color="#4C72B0", alpha=0.85, capsize=3)
        ax.set_yticks(range(len(tags)))
        ax.set_yticklabels(tags if metric == metrics[0] else [], fontsize=7)
        ax.set_xlabel(metric)
        ax.set_title(metric)
        ax.set_xlim(0, 1)

    fig.suptitle("Metricas en prueba con IC 95% (bootstrap por cluster)", y=1.02)
    fig.savefig(FIGURES_DIR / "metrics_comparison.png")
    plt.close(fig)


def plot_confusion(results: dict[str, pd.DataFrame], thresholds: dict[str, float]) -> None:
    plt = _setup_matplotlib()
    tags = sorted(results)
    cols = min(len(tags), 4)
    rows = (len(tags) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 3.2 * rows), squeeze=False)

    for ax, tag in zip(axes.ravel(), tags):
        df = results[tag]
        pred = (df.score.to_numpy() >= thresholds[tag]).astype(int)
        cm = confusion_matrix(df.label.to_numpy(), pred, labels=[0, 1])
        ax.imshow(cm, cmap="Blues")
        for i in range(2):
            for j in range(2):
                ax.text(
                    j,
                    i,
                    f"{cm[i, j]:,}",
                    ha="center",
                    va="center",
                    color="white" if cm[i, j] > cm.max() / 2 else "black",
                    fontsize=9,
                )
        ax.set_xticks([0, 1], ["no-AMP", "AMP"])
        ax.set_yticks([0, 1], ["no-AMP", "AMP"])
        ax.set_xlabel("Prediccion")
        ax.set_ylabel("Real")
        ax.set_title(tag, fontsize=7)
        ax.grid(False)

    for ax in axes.ravel()[len(tags) :]:
        ax.axis("off")
    fig.suptitle("Matrices de confusion en prueba", y=1.0)
    fig.savefig(FIGURES_DIR / "confusion_matrices.png")
    plt.close(fig)


def plot_calibration(results: dict[str, pd.DataFrame]) -> None:
    """Curva de calibracion: ¿una probabilidad de 0.8 significa 80%?"""
    from sklearn.calibration import calibration_curve

    plt = _setup_matplotlib()
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5, label="calibracion perfecta")

    for tag, df in sorted(results.items()):
        try:
            frac, mean_pred = calibration_curve(df.label, df.score, n_bins=10, strategy="quantile")
            ax.plot(mean_pred, frac, "o-", lw=1.3, ms=4, label=tag)
        except ValueError:
            continue

    ax.set_xlabel("Probabilidad predicha media")
    ax.set_ylabel("Fraccion real de positivos")
    ax.set_title("Calibracion en el conjunto de prueba")
    ax.legend(fontsize=7)
    fig.savefig(FIGURES_DIR / "calibration.png")
    plt.close(fig)


# --- Orquestacion ---------------------------------------------------------


# Par preregistrado en docs/07, seccion "Compromiso previo". Si ambas etiquetas
# existen entre las predicciones, la comparacion usa exactamente estas, que es
# lo que el preregistro promete. Si no, se cae a la mejor por validacion.
PREREGISTERED_EMBEDDING = "facebook_esm2_t12_35M_UR50D__logreg"
PREREGISTERED_BASELINE = "comp_physchem__xgboost"


def _select_comparison_pair(
    summary: pd.DataFrame, val_scores: dict[str, float]
) -> tuple[str | None, str | None, str]:
    """Elige los dos modelos a comparar SIN mirar test.

    Devuelve (modelo_embeddings, modelo_baseline, criterio_usado).

    Dos criterios, en orden:

    1. Las etiquetas preregistradas en docs/07, si existen. Es lo mas fuerte
       posible: el par se fijo por escrito antes de correr nada.
    2. El mejor de cada familia por ROC-AUC de **validacion**.

    Lo que nunca se usa es el resultado en test, que es lo que hacia la version
    anterior (P-05).
    """
    tags = list(summary.tag)
    emb_tags = [t for t in tags if t.startswith("facebook_")]
    base_tags = [t for t in tags if not t.startswith("facebook_")]

    if not emb_tags or not base_tags:
        return None, None, "no aplicable: falta una de las dos familias"

    if PREREGISTERED_EMBEDDING in emb_tags and PREREGISTERED_BASELINE in base_tags:
        return PREREGISTERED_EMBEDDING, PREREGISTERED_BASELINE, "preregistro (docs/07)"

    best_emb = max(emb_tags, key=lambda t: val_scores.get(t, -1.0))
    best_base = max(base_tags, key=lambda t: val_scores.get(t, -1.0))
    return best_emb, best_base, "mejor ROC-AUC en validacion"


def evaluate_all() -> pd.DataFrame:
    """Evalua todas las predicciones guardadas y escribe metricas y figuras."""
    ensure_dirs()
    params = load_params()
    ecfg = params.evaluate
    n_boot = int(ecfg.get("bootstrap_n", 2000))
    unit = str(ecfg.get("bootstrap_unit", "cluster"))
    ci = float(ecfg.get("ci", 0.95))
    threshold_metric = str(ecfg.get("threshold_metric", "mcc"))

    files = sorted(PREDICTIONS_DIR.glob("*.csv"))
    if not files:
        raise FileNotFoundError(
            f"No hay predicciones en {PREDICTIONS_DIR}. Corre 'make baselines' y 'make train'."
        )

    test_frames: dict[str, pd.DataFrame] = {}
    thresholds: dict[str, float] = {}
    val_scores_by_tag: dict[str, float] = {}
    rows = []

    print(f"[eval] {len(files)} modelos; bootstrap de {n_boot} por {unit}")

    for path in files:
        tag = path.stem
        df = pd.read_csv(path)
        val = df[df.split == "val"]
        test = df[df.split == "test"]

        threshold, val_value = choose_threshold(
            val.label.to_numpy(), val.score.to_numpy(), threshold_metric
        )
        thresholds[tag] = threshold
        test_frames[tag] = test
        # Se guarda el ROC-AUC de VALIDACION, que es el criterio con el que se
        # elige el par a comparar (P-05). Nunca el de test.
        val_scores_by_tag[tag] = float(roc_auc_score(val.label, val.score))

        metrics = bootstrap_metrics(test, threshold, n_boot, unit, ci, params.seed)
        pred = (test.score.to_numpy() >= threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(test.label, pred, labels=[0, 1]).ravel()

        record = {
            "model": tag,
            "split": "test",
            "threshold": threshold,
            "threshold_selected_on": "val",
            "threshold_metric": threshold_metric,
            f"val_{threshold_metric}": val_value,
            "metrics": metrics,
            "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
            "n_test": int(len(test)),
            "n_clusters_test": int(test.cluster_id.nunique()),
            "bootstrap_n": n_boot,
            "bootstrap_unit": unit,
            "ci": ci,
        }
        (METRICS_DIR / f"{tag}.json").write_text(
            json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8"
        )

        row = {"tag": tag, "threshold": threshold, "n_test": len(test)}
        for name, values in metrics.items():
            row[name] = values["point"]
            row[f"{name}_ci_low"] = values["ci_low"]
            row[f"{name}_ci_high"] = values["ci_high"]
        rows.append(row)

        print(
            f"[eval] {tag:48s} ROC-AUC {metrics['roc_auc']['point']:.4f} "
            f"[{metrics['roc_auc']['ci_low']:.4f}, {metrics['roc_auc']['ci_high']:.4f}]  "
            f"MCC {metrics['mcc']['point']:.4f}"
        )

    summary = pd.DataFrame(rows).sort_values("mcc", ascending=False).reset_index(drop=True)
    summary.to_csv(METRICS_DIR / "summary.csv", index=False)

    # Comparacion emparejada: la pregunta central del proyecto (docs/06).
    #
    # P-05: antes se tomaba `summary.iloc[0]`, es decir el mejor de 13 modelos
    # ordenados por MCC **de test**, y se medía la diferencia en ese mismo test.
    # Elegir en test y medir en test es winner's curse: la diferencia puntual
    # sale inflada y el IC bootstrap NO lo corrige, porque el bootstrap no sabe
    # que hubo una seleccion previa. Ademas contradecia dos compromisos escritos
    # en este repositorio (docs/02 y el preregistro de docs/07).
    #
    # Ahora la seleccion es por validacion, y si las etiquetas preregistradas
    # existen se usan esas directamente, que es lo que el preregistro promete.
    best_emb, best_base, selection_basis = _select_comparison_pair(summary, val_scores_by_tag)
    if best_emb and best_base:
        comparison = paired_difference(
            test_frames[best_emb],
            test_frames[best_base],
            thresholds[best_emb],
            thresholds[best_base],
            metric=threshold_metric,
            n_boot=n_boot,
            unit=unit,
            ci=ci,
            seed=params.seed,
        )
        comparison["model_a"] = best_emb
        comparison["model_b"] = best_base
        comparison["selected_on"] = selection_basis
        comparison["selection_note"] = (
            "Los dos modelos comparados se eligieron SIN mirar test. La tabla "
            "summary.csv con los 13 modelos es descriptiva; la unica inferencia "
            "del proyecto es esta comparacion emparejada. Ver docs/09, P-05."
        )
        (METRICS_DIR / "comparison_embeddings_vs_baseline.json").write_text(
            json.dumps(comparison, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        print(f"\n[eval] comparacion preregistrada (seleccion: {selection_basis})")
        print(f"[eval] {best_emb}  frente a  {best_base}")
        print(
            f"[eval] diferencia de {threshold_metric}: {comparison['difference']:+.4f} "
            f"[{comparison['ci_low']:+.4f}, {comparison['ci_high']:+.4f}]"
        )
        print(f"[eval] {comparison['interpretation']}")

    plot_curves(test_frames)
    plot_metrics_comparison(summary)
    plot_confusion(test_frames, thresholds)
    plot_calibration(test_frames)
    print(f"\n[eval] metricas en {METRICS_DIR}, figuras en {FIGURES_DIR}")
    return summary
