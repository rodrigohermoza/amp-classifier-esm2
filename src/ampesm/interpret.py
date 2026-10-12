"""Interpretacion: UMAP, analisis de errores y controles de sesgo.

Esta etapa es la que dice si el modelo aprendio biologia o aprendio procedencia.
Las figuras de aqui pesan mas en las conclusiones que cualquier metrica
agregada, porque una metrica alta es compatible con haber aprendido lo
equivocado.

Detalle de las preguntas que se responden en docs/07-evaluacion.md.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from .config import FIGURES_DIR, METRICS_DIR, MODELS_DIR, ensure_dirs, load_params
from .evaluate import _setup_matplotlib, choose_threshold
from .features import physchem
from .models import PREDICTIONS_DIR, load_splits


def plot_length_distribution(df: pd.DataFrame) -> None:
    """Longitud por clase y por fuente: control de sesgo basico.

    Si las distribuciones de longitud por clase no se solapan, el modelo puede
    separar las clases contando residuos. Esta figura es la comprobacion visual
    de que el emparejado de longitudes de clean.py hizo su trabajo.
    """
    plt = _setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    bins = np.arange(df.length.min(), df.length.max() + 2) - 0.5
    for label, name, color in ((1, "AMP", "#C44E52"), (0, "no-AMP", "#4C72B0")):
        sub = df[df.label == label]
        axes[0].hist(sub.length, bins=bins, alpha=0.6, label=f"{name} (n={len(sub)})", color=color)
    axes[0].set_xlabel("Longitud (residuos)")
    axes[0].set_ylabel("Numero de peptidos")
    axes[0].set_title("Longitud por clase")
    axes[0].legend(fontsize=8)

    for source in sorted(df.source.unique()):
        sub = df[df.source == source]
        axes[1].hist(
            sub.length, bins=bins, histtype="step", lw=1.5, label=f"{source} (n={len(sub)})"
        )
    axes[1].set_xlabel("Longitud (residuos)")
    axes[1].set_title("Longitud por fuente")
    axes[1].legend(fontsize=7)

    fig.suptitle("Control de sesgo: distribucion de longitudes", y=1.02)
    fig.savefig(FIGURES_DIR / "length_distribution.png")
    plt.close(fig)


def plot_umap(df: pd.DataFrame, representation: str) -> dict:
    """UMAP de los embeddings, coloreado por clase, por fuente y por split.

    Colorear por fuente es el control que importa y el que suele omitirse. Si
    los grupos se separan por base de datos en vez de por clase, el modelo tiene
    un atajo disponible y cualquier metrica alta es sospechosa.

    Advertencia de lectura: en un UMAP las distancias no son fiables. Dos puntos
    cercanos en el dibujo no estan necesariamente cercanos en el espacio
    original. Sirve para ver si hay agrupacion, no para medir cuanta.
    """
    import umap

    from .models import get_feature_matrix

    params = load_params()
    x, _ = get_feature_matrix(df, representation, params)
    ucfg = params.interpret.get("umap", {})

    reducer = umap.UMAP(
        n_neighbors=int(ucfg.get("n_neighbors", 15)),
        min_dist=float(ucfg.get("min_dist", 0.1)),
        metric=str(ucfg.get("metric", "euclidean")),
        random_state=params.seed,
    )
    coords = reducer.fit_transform(x)

    plt = _setup_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    for label, name, color in ((0, "no-AMP", "#4C72B0"), (1, "AMP", "#C44E52")):
        mask = (df.label == label).to_numpy()
        axes[0].scatter(
            coords[mask, 0], coords[mask, 1], s=4, alpha=0.5, label=name, color=color, lw=0
        )
    axes[0].set_title("Por clase")
    axes[0].legend(fontsize=8, markerscale=3)

    for source in sorted(df.source.unique()):
        mask = (df.source == source).to_numpy()
        axes[1].scatter(coords[mask, 0], coords[mask, 1], s=4, alpha=0.5, label=source, lw=0)
    axes[1].set_title("Por fuente (el control que importa)")
    axes[1].legend(fontsize=7, markerscale=3)

    for split in ("train", "val", "test"):
        mask = (df.split == split).to_numpy()
        axes[2].scatter(coords[mask, 0], coords[mask, 1], s=4, alpha=0.5, label=split, lw=0)
    axes[2].set_title("Por split")
    axes[2].legend(fontsize=8, markerscale=3)

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)

    fig.suptitle(f"UMAP de los embeddings de {representation}", y=1.02)
    fig.savefig(FIGURES_DIR / "umap_embeddings.png")
    plt.close(fig)

    return {"representation": representation, "n_points": int(len(df))}


def error_analysis(df: pd.DataFrame, tag: str) -> dict:
    """Compara descriptores de aciertos y fallos, y desglosa por fuente.

    Las cuatro preguntas de docs/07 se responden aqui con numeros, no con
    impresiones.
    """
    path = PREDICTIONS_DIR / f"{tag}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Falta {path}")

    preds = pd.read_csv(path)
    val = preds[preds.split == "val"]
    threshold, _ = choose_threshold(val.label.to_numpy(), val.score.to_numpy(), "mcc")

    test = preds[preds.split == "test"].merge(
        df[["id", "sequence", "source", "length"]], on="id", how="left"
    )
    # P-10a: si hay predicciones de una corrida anterior, el cruce deja
    # sequence=NaN para los ids que ya no existen y physchem_matrix revienta con
    # un TypeError opaco ("'float' object is not iterable"). Pasa al re-correr
    # 'make data' o cambiar negative_ratio sin borrar las predicciones.
    if test.sequence.isna().any():
        n = int(test.sequence.isna().sum())
        raise RuntimeError(
            f"{n} predicciones de {tag} corresponden a ids que ya no estan en "
            f"splits.csv: son de una corrida anterior. Corre 'make clean' y "
            f"reentrena ('make baselines' y 'make train')."
        )
    test["pred"] = (test.score >= threshold).astype(int)
    test["correct"] = test.pred == test.label

    # Descriptores fisicoquimicos de cada peptido de prueba.
    # P-10f: con los valores por defecto, si alguien cambia el pH o la ventana
    # en params.yaml el analisis de errores describiria los peptidos con
    # descriptores distintos de los que vio el modelo, sin ningun aviso.
    fcfg = load_params().features
    desc = physchem.physchem_matrix(
        test.sequence.tolist(),
        ph=float(fcfg.get("physchem_ph", 7.4)),
        window=int(fcfg.get("moment_window", 11)),
        angle=float(fcfg.get("moment_angle_deg", 100.0)),
    )
    for i, name in enumerate(physchem.FEATURE_NAMES):
        test[name] = desc[:, i]

    def category(row):
        if row.label == 1:
            return "verdadero_positivo" if row.pred == 1 else "falso_negativo"
        return "falso_positivo" if row.pred == 1 else "verdadero_negativo"

    test["category"] = test.apply(category, axis=1)

    # Pregunta 1 y 2: ¿en que se diferencian los fallos de los aciertos?
    by_category = {
        cat: {
            "n": int(len(sub)),
            **{
                name: round(float(sub[name].mean()), 3)
                for name in ("net_charge", "gravy", "hydrophobic_moment", "length", "aromaticity")
            },
        }
        for cat, sub in test.groupby("category")
    }

    # Pregunta 3: ¿el error se concentra en alguna fuente?
    by_source = {}
    for source, sub in test.groupby("source"):
        entry = {
            "n": int(len(sub)),
            "accuracy": round(float(sub.correct.mean()), 4),
            "label": int(sub.label.iloc[0]),
        }
        if sub.label.nunique() > 1:
            entry["roc_auc"] = round(float(roc_auc_score(sub.label, sub.score)), 4)
        by_source[source] = entry

    # Pregunta 4: ¿el error depende de la longitud?
    # P-10g: derivado de params.yaml. Con el borde cableado en 51, subir
    # max_length dejaria los peptidos mas largos en NaN y groupby(observed=True)
    # los descartaria del reporte en silencio.
    max_len = int(load_params().data.get("max_length", 50))
    bins = [4, 10, 20, 30, 40, max(51, max_len + 1)]
    test["length_bin"] = pd.cut(test.length, bins=bins, right=False)
    by_length = {
        str(interval): {"n": int(len(sub)), "accuracy": round(float(sub.correct.mean()), 4)}
        for interval, sub in test.groupby("length_bin", observed=True)
    }

    # Los falsos positivos de mayor puntuacion son los candidatos a ser
    # negativos mal etiquetados: peptidos cationicos de funcion no anotada que
    # probablemente si son antimicrobianos.
    top_fp = (
        test[test.category == "falso_positivo"]
        .nlargest(15, "score")[["id", "sequence", "score", "net_charge", "length"]]
        .to_dict("records")
    )
    top_fn = (
        test[test.category == "falso_negativo"]
        .nsmallest(15, "score")[["id", "sequence", "score", "net_charge", "length"]]
        .to_dict("records")
    )

    report = {
        "model": tag,
        "threshold": threshold,
        "n_test": int(len(test)),
        "by_category": by_category,
        "by_source": by_source,
        "by_source_note": (
            "Si la exactitud en uniprot_amp es mucho peor que en apd3, el modelo "
            "depende de la procedencia y no de la actividad. Ver docs/01-objetivos.md."
        ),
        "by_length": by_length,
        "top_false_positives": top_fp,
        "top_false_positives_note": (
            "Candidatos a negativos mal etiquetados. Si son cationicos y secretados, "
            "es probable que sean AMPs no anotados, y entonces la precision real es "
            "mejor que la medida."
        ),
        "top_false_negatives": top_fn,
        "top_false_negatives_note": (
            "Si son anionicos o ricos en prolina, son familias de AMPs que no siguen "
            "el mecanismo cationico-anfipatico dominante, y el fallo es explicable."
        ),
    }
    (METRICS_DIR / f"error_analysis_{tag}.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )

    # Figura: descriptores de aciertos frente a fallos.
    plt = _setup_matplotlib()
    features = ["net_charge", "gravy", "hydrophobic_moment", "length"]
    fig, axes = plt.subplots(1, len(features), figsize=(4.0 * len(features), 4.2))
    order = ["verdadero_positivo", "falso_negativo", "falso_positivo", "verdadero_negativo"]
    present = [c for c in order if c in set(test.category)]

    for ax, feature in zip(axes, features):
        ax.boxplot(
            [test[test.category == c][feature].to_numpy() for c in present],
            tick_labels=[c.replace("_", "\n") for c in present],
            showfliers=False,
        )
        ax.set_title(feature)
        ax.tick_params(axis="x", labelsize=7)

    fig.suptitle(f"Analisis de errores: {tag}", y=1.03)
    fig.savefig(FIGURES_DIR / "error_analysis.png")
    plt.close(fig)

    print(f"[interpret] analisis de errores de {tag}:")
    for cat, info in sorted(by_category.items()):
        print(f"  {cat:20s} n={info['n']:5d}  carga media {info['net_charge']:+.2f}")
    print("[interpret] exactitud por fuente:")
    for source, info in sorted(by_source.items()):
        print(f"  {source:14s} n={info['n']:5d}  exactitud {info['accuracy']:.3f}")

    return report


def _best_model_by_validation() -> str:
    """Modelo a interpretar, elegido por VALIDACION.

    P-05: antes se tomaba el primero de summary.csv, que esta ordenado por MCC
    de test. Interpretar el ganador en test no infla una metrica (el analisis de
    errores es descriptivo), pero sesga la narrativa: se acaba contando la
    historia del modelo que mejor le fue al conjunto de prueba. Los ROC-AUC de
    validacion ya estan guardados por save_result en results/models/*.json.
    """
    metas = sorted(MODELS_DIR.glob("*.json"))
    if metas:
        best, best_auc = None, -1.0
        for path in metas:
            try:
                meta = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            auc = float(meta.get("val_roc_auc", -1.0))
            if auc > best_auc and (PREDICTIONS_DIR / f"{path.stem}.csv").exists():
                best, best_auc = path.stem, auc
        if best:
            print(f"[interpret] modelo elegido por val ROC-AUC {best_auc:.4f}")
            return best

    available = sorted(p.stem for p in PREDICTIONS_DIR.glob("*.csv"))
    if not available:
        raise FileNotFoundError("No hay predicciones; corre 'make baselines' y 'make train'")
    return available[0]


def interpret_all(tag: str | None = None) -> None:
    """Corre todas las figuras de interpretacion."""
    ensure_dirs()
    df = load_splits()
    params = load_params()

    plot_length_distribution(df)
    print("[interpret] distribucion de longitudes")

    # Si no se indica modelo, se usa el de mejor ROC-AUC en VALIDACION (P-05).
    if tag is None:
        tag = _best_model_by_validation()

    print(f"[interpret] modelo elegido: {tag}")
    error_analysis(df, tag)

    models = params.embeddings.get("models", [])
    if models:
        representation = models[-1]  # el mas grande de los configurados
        print(f"[interpret] UMAP con {representation}")
        plot_umap(df, representation)

    print(f"[interpret] figuras en {FIGURES_DIR}")
