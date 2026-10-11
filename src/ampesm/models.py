"""Representaciones, modelos y entrenamiento.

Separacion que se mantiene en todo el modulo: aqui se entrena y se guardan
puntuaciones; las metricas se calculan en evaluate.py a partir de esas
puntuaciones. La razon es que recalcular metricas, umbrales o intervalos de
confianza no deberia obligar a reentrenar nada.

Cada combinacion de representacion y modelo escribe un CSV de predicciones con
las columnas id, cluster_id, split, label y score. Ese archivo es el unico
contrato con la etapa de evaluacion.
"""

from __future__ import annotations

import itertools
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import METRICS_DIR, MODELS_DIR, SPLITS_PATH, ensure_dirs, load_params
from .embeddings import embedding_matrix
from .features import composition, physchem

PREDICTIONS_DIR = METRICS_DIR / "predictions"


# --- Representaciones ------------------------------------------------------


def get_feature_matrix(
    df: pd.DataFrame, representation: str, params
) -> tuple[np.ndarray, list[str]]:
    """Construye la matriz de features para una representacion dada.

    `representation` es "composition", "physchem", "comp_physchem" o el nombre
    de un modelo de ESM-2 en HuggingFace.
    """
    fcfg = params.features
    ph = float(fcfg.get("physchem_ph", 7.4))
    window = int(fcfg.get("moment_window", 11))
    angle = float(fcfg.get("moment_angle_deg", 100.0))
    sequences = df.sequence.tolist()

    if representation == "composition":
        return composition.composition_matrix(sequences), composition.FEATURE_NAMES
    if representation == "physchem":
        return physchem.physchem_matrix(sequences, ph, window, angle), physchem.FEATURE_NAMES
    if representation == "comp_physchem":
        x = np.hstack(
            [
                composition.composition_matrix(sequences),
                physchem.physchem_matrix(sequences, ph, window, angle),
            ]
        )
        return x, composition.FEATURE_NAMES + physchem.FEATURE_NAMES

    # Cualquier otra cosa se interpreta como un modelo de ESM-2.
    ecfg = params.embeddings
    x = embedding_matrix(
        df,
        representation,
        layer=ecfg.get("layer", "last"),
        pooling=ecfg.get("pooling", "mean"),
        batch_size=int(ecfg.get("batch_size", 16)),
        device=str(ecfg.get("device", "cpu")),
    )
    return x, [f"emb_{i}" for i in range(x.shape[1])]


# --- Modelos ---------------------------------------------------------------


def expand_grid(grid: dict) -> list[dict]:
    """Producto cartesiano de un diccionario de listas."""
    if not grid:
        return [{}]
    keys = list(grid)
    values = [v if isinstance(v, list) else [v] for v in (grid[k] for k in keys)]
    return [dict(zip(keys, combo)) for combo in itertools.product(*values)]


def make_estimator(kind: str, hp: dict, seed: int):
    """Construye un estimador de scikit-learn ya configurado.

    Los modelos lineales van dentro de un Pipeline con escalado. Los basados en
    arboles no lo necesitan (son invariantes a transformaciones monotonas de
    cada feature) y meterlo solo anadiria una etapa que confunde al leer el
    modelo guardado.
    """
    if kind == "logreg":
        return Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        C=hp.get("C", 1.0),
                        max_iter=int(hp.get("max_iter", 5000)),
                        class_weight="balanced",
                        random_state=seed,
                    ),
                ),
            ]
        )
    if kind == "random_forest":
        return RandomForestClassifier(
            n_estimators=int(hp.get("n_estimators", 300)),
            max_depth=hp.get("max_depth"),
            min_samples_leaf=int(hp.get("min_samples_leaf", 1)),
            class_weight="balanced",
            n_jobs=-1,
            random_state=seed,
        )
    if kind == "xgboost":
        from xgboost import XGBClassifier

        return XGBClassifier(
            n_estimators=int(hp.get("n_estimators", 300)),
            max_depth=int(hp.get("max_depth", 6)),
            learning_rate=float(hp.get("learning_rate", 0.1)),
            subsample=float(hp.get("subsample", 0.8)),
            eval_metric="logloss",
            n_jobs=-1,
            random_state=seed,
        )
    raise ValueError(f"Modelo desconocido: {kind!r}")


# --- MLP -------------------------------------------------------------------


def fit_mlp(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_val: np.ndarray,
    y_val: np.ndarray,
    hp: dict,
    seed: int,
):
    """MLP pequeno en PyTorch, con parada temprana por ROC-AUC en validacion.

    Se mantiene pequeno a proposito: con unos pocos miles de ejemplos y 320-480
    dimensiones, una red grande se sobreajusta sin aportar nada. La parada
    temprana se mide en validacion, nunca en prueba.

    Devuelve (funcion_de_puntuacion, informacion_de_entrenamiento). La funcion
    devuelta encapsula el escalado, para que quien la use no tenga que recordar
    aplicarlo.
    """
    import torch
    from torch import nn

    torch.manual_seed(seed)

    scaler = StandardScaler().fit(x_train)
    xt = torch.tensor(scaler.transform(x_train), dtype=torch.float32)
    yt = torch.tensor(y_train, dtype=torch.float32).unsqueeze(1)
    xv = torch.tensor(scaler.transform(x_val), dtype=torch.float32)

    hidden_sizes = hp.get("hidden_sizes", [256])
    dropout = float(hp.get("dropout", 0.3))

    layers: list[nn.Module] = []
    in_dim = xt.shape[1]
    for size in hidden_sizes:
        layers += [nn.Linear(in_dim, int(size)), nn.ReLU(), nn.Dropout(dropout)]
        in_dim = int(size)
    layers.append(nn.Linear(in_dim, 1))
    model = nn.Sequential(*layers)

    # Peso de la clase positiva para compensar el desbalance, igual que el
    # class_weight="balanced" de los modelos de scikit-learn.
    n_pos = float(y_train.sum())
    n_neg = float(len(y_train) - n_pos)
    pos_weight = torch.tensor([n_neg / max(n_pos, 1.0)])

    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.Adam(model.parameters(), lr=float(hp.get("lr", 1e-3)))

    epochs = int(hp.get("epochs", 60))
    batch_size = int(hp.get("batch_size", 32))
    patience = int(hp.get("patience", 10))

    best_auc = -1.0
    best_state = None
    best_epoch = 0
    generator = torch.Generator().manual_seed(seed)

    for epoch in range(1, epochs + 1):
        model.train()
        perm = torch.randperm(len(xt), generator=generator)
        for start in range(0, len(xt), batch_size):
            idx = perm[start : start + batch_size]
            optimizer.zero_grad()
            loss = criterion(model(xt[idx]), yt[idx])
            loss.backward()
            optimizer.step()

        model.eval()
        with torch.no_grad():
            val_scores = torch.sigmoid(model(xv)).squeeze(1).numpy()
        auc = roc_auc_score(y_val, val_scores)

        if auc > best_auc:
            best_auc = auc
            best_epoch = epoch
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        elif epoch - best_epoch >= patience:
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()

    def score_fn(x: np.ndarray) -> np.ndarray:
        with torch.no_grad():
            xx = torch.tensor(scaler.transform(x), dtype=torch.float32)
            return torch.sigmoid(model(xx)).squeeze(1).numpy()

    return score_fn, {"best_epoch": best_epoch, "val_roc_auc": float(best_auc)}


# --- Entrenamiento ---------------------------------------------------------


@dataclass
class TrainResult:
    tag: str
    representation: str
    model: str
    best_params: dict
    val_roc_auc: float
    val_mcc: float
    predictions: pd.DataFrame


def _split_masks(df: pd.DataFrame):
    return (
        (df.split == "train").to_numpy(),
        (df.split == "val").to_numpy(),
        (df.split == "test").to_numpy(),
    )


def train_one(
    df: pd.DataFrame,
    x: np.ndarray,
    representation: str,
    model_kind: str,
    grid: dict,
    seed: int,
) -> TrainResult:
    """Entrena un modelo, elige hiperparametros en validacion y puntua todo.

    La seleccion es por ROC-AUC en validacion. No se usa validacion cruzada:
    los pliegues tendrian que respetar los clusteres y ya tenemos una particion
    que los respeta, asi que un val fijo es mas simple y mas facil de auditar.
    """
    y = df.label.to_numpy()
    tr, va, _ = _split_masks(df)

    combos = expand_grid(grid)
    best = None

    for hp in combos:
        if model_kind == "mlp":
            score_fn, info = fit_mlp(x[tr], y[tr], x[va], y[va], hp, seed)
            val_scores = score_fn(x[va])
            fitted = score_fn
        else:
            est = make_estimator(model_kind, hp, seed)
            est.fit(x[tr], y[tr])
            val_scores = est.predict_proba(x[va])[:, 1]
            fitted = est
            info = {}

        auc = roc_auc_score(y[va], val_scores)
        if best is None or auc > best["auc"]:
            best = {"auc": auc, "hp": hp, "fitted": fitted, "scores": val_scores, "info": info}

    assert best is not None  # expand_grid siempre devuelve al menos un combo

    fitted = best["fitted"]
    all_scores = fitted(x) if model_kind == "mlp" else fitted.predict_proba(x)[:, 1]

    # El MCC en validacion se calcula con umbral 0.5, solo como referencia
    # rapida. El umbral de verdad se optimiza en evaluate.py, tambien sobre
    # validacion, pero maximizando la metrica configurada.
    val_mcc = matthews_corrcoef(y[va], (best["scores"] >= 0.5).astype(int))

    tag = f"{representation.replace('/', '_')}__{model_kind}"
    predictions = pd.DataFrame(
        {
            "id": df.id.to_numpy(),
            "cluster_id": df.cluster_id.to_numpy(),
            "split": df.split.to_numpy(),
            "label": y,
            "score": all_scores,
        }
    )

    return TrainResult(
        tag=tag,
        representation=representation,
        model=model_kind,
        best_params={**best["hp"], **best["info"]},
        val_roc_auc=float(best["auc"]),
        val_mcc=float(val_mcc),
        predictions=predictions,
    )


def save_result(result: TrainResult) -> Path:
    PREDICTIONS_DIR.mkdir(parents=True, exist_ok=True)
    path = PREDICTIONS_DIR / f"{result.tag}.csv"
    result.predictions.to_csv(path, index=False)

    meta = {
        "tag": result.tag,
        "representation": result.representation,
        "model": result.model,
        "best_params": {k: (None if v is None else v) for k, v in result.best_params.items()},
        "val_roc_auc": result.val_roc_auc,
        "val_mcc_at_0.5": result.val_mcc,
    }
    (MODELS_DIR / f"{result.tag}.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    return path


def load_splits() -> pd.DataFrame:
    if not SPLITS_PATH.exists():
        raise FileNotFoundError(f"Falta {SPLITS_PATH}. Corre 'make split' primero.")
    return pd.read_csv(SPLITS_PATH)


def run_baselines() -> list[TrainResult]:
    """Entrena los baselines tabulares: 3 representaciones x 3 modelos."""
    ensure_dirs()
    params = load_params()
    df = load_splits()
    results = []

    for representation in ("composition", "physchem", "comp_physchem"):
        x, _names = get_feature_matrix(df, representation, params)
        print(f"[baselines] {representation}: matriz {x.shape}")
        for model_kind in ("logreg", "random_forest", "xgboost"):
            grid = params.models.get(model_kind, {})
            result = train_one(df, x, representation, model_kind, grid, params.seed)
            save_result(result)
            print(
                f"[baselines]   {model_kind:14s} val ROC-AUC {result.val_roc_auc:.4f}  "
                f"params {result.best_params}"
            )
            results.append(result)

    return results


def run_embedding_models() -> list[TrainResult]:
    """Entrena logreg y MLP sobre los embeddings de cada modelo de ESM-2."""
    ensure_dirs()
    params = load_params()
    df = load_splits()
    model_names = params.embeddings.get("models", [])
    results = []

    for representation in model_names:
        x, _names = get_feature_matrix(df, representation, params)
        print(f"[train] {representation}: matriz {x.shape}")
        for model_kind in ("logreg", "mlp"):
            grid = params.models.get(model_kind, {})
            result = train_one(df, x, representation, model_kind, grid, params.seed)
            save_result(result)
            print(
                f"[train]   {model_kind:14s} val ROC-AUC {result.val_roc_auc:.4f}  "
                f"params {result.best_params}"
            )
            results.append(result)

    return results
