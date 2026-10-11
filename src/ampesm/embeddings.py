"""Embeddings de ESM-2 con cache en disco.

Es la etapa mas caro de la tuberia, asi que se cachea agresivamente: indexado
por hash de secuencia, de modo que reentrenar un clasificador con otros
hiperparametros tarda segundos en vez de minutos.

Contexto y eleccion de modelo en docs/06-embeddings-esm2.md.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import EMBEDDINGS_DIR, hf_cache_dir

# Dimension del embedding por modelo. Solo informativo, para mensajes y para
# detectar pronto que un modelo devolvio algo con la forma equivocada.
MODEL_DIMS = {
    "facebook/esm2_t6_8M_UR50D": 320,
    "facebook/esm2_t12_35M_UR50D": 480,
    "facebook/esm2_t30_150M_UR50D": 640,
    "facebook/esm2_t33_650M_UR50D": 1280,
}

MODEL_SIZES_MB = {
    "facebook/esm2_t6_8M_UR50D": 31,
    "facebook/esm2_t12_35M_UR50D": 140,
    "facebook/esm2_t30_150M_UR50D": 600,
    "facebook/esm2_t33_650M_UR50D": 2500,
}


def cache_path(model_name: str, layer: str | int, pooling: str) -> Path:
    """Ruta del cache. El nombre codifica la configuracion completa.

    Asi un cambio de capa o de agregacion no reutiliza silenciosamente
    embeddings calculados con otra configuracion, que es un error que no daria
    ningun aviso y falsearia la comparacion.
    """
    slug = model_name.replace("/", "_")
    return EMBEDDINGS_DIR / f"{slug}__layer-{layer}__{pooling}.npz"


def load_cache(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {}
    with np.load(path) as data:
        return {k: data[k] for k in data.files}


def save_cache(path: Path, cache: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Escritura atomica: un .npz truncado por una interrupcion seria peor que
    # ninguno, porque la siguiente corrida lo cargaria como valido.
    tmp = path.with_suffix(".npz.part")
    # Se escribe sobre un descriptor abierto y no sobre la ruta: si se le pasa
    # una ruta, numpy le anade ".npz" cuando no termina en esa extension, y el
    # renombrado posterior apuntaria a un archivo que no existe.
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, **cache)
    tmp.replace(path)


def masked_mean_pool(hidden, attention_mask, input_ids, special_ids):
    """Promedia sobre residuos reales, excluyendo relleno y tokens especiales.

    Este es el punto de la tuberia mas facil de arruinar sin que nada avise. Si
    se promedia sin enmascarar, los vectores de los peptidos cortos quedan
    diluidos por los ceros del relleno, de forma proporcional a su longitud: el
    clasificador acaba leyendo la longitud en el propio embedding y el resultado
    parece bueno por el motivo equivocado.

    Ademas del relleno hay que excluir <cls> y <eos>. ESM-2 se entrena solo con
    enmascaramiento de tokens, sin objetivo a nivel de secuencia, asi que <cls>
    no tiene un significado global garantizado y meterlo en la media solo anade
    una constante que difumina la senal de los peptidos cortos.

    tests/test_embeddings.py verifica que el embedding de una secuencia
    procesada sola coincide con el de la misma dentro de un lote con relleno.
    """
    mask = attention_mask.clone()
    for sid in special_ids:
        if sid is not None:
            mask = mask & (input_ids != sid)
    mask = mask.unsqueeze(-1).to(hidden.dtype)

    summed = (hidden * mask).sum(dim=1)
    counts = mask.sum(dim=1).clamp(min=1.0)  # clamp: evita dividir por cero
    return summed / counts


def masked_max_pool(hidden, attention_mask, input_ids, special_ids):
    """Maximo por dimension sobre residuos reales."""
    import torch

    mask = attention_mask.clone()
    for sid in special_ids:
        if sid is not None:
            mask = mask & (input_ids != sid)
    neg_inf = torch.finfo(hidden.dtype).min
    masked = hidden.masked_fill(~mask.unsqueeze(-1).bool(), neg_inf)
    return masked.max(dim=1).values


def compute_embeddings(
    sequences: dict[str, str],
    model_name: str,
    layer: str | int = "last",
    pooling: str = "mean",
    batch_size: int = 16,
    device: str = "cpu",
    use_cache: bool = True,
) -> dict[str, np.ndarray]:
    """Calcula embeddings para {hash: secuencia}, usando y actualizando el cache.

    Solo se pasa por el modelo lo que falta en el cache. Los lotes se ordenan
    por longitud para que cada lote tenga relleno minimo, lo que en CPU es una
    diferencia de tiempo apreciable.
    """
    import torch
    from transformers import AutoModel, AutoTokenizer

    path = cache_path(model_name, layer, pooling)
    cache = load_cache(path) if use_cache else {}
    todo = {h: s for h, s in sequences.items() if h not in cache}

    if not todo:
        print(f"[embeddings] {len(sequences)} vectores ya en cache: {path.name}")
        return {h: cache[h] for h in sequences}

    size = MODEL_SIZES_MB.get(model_name)
    tamano = f" (~{size} MB la primera vez)" if size else ""
    print(f"[embeddings] {model_name}{tamano}")
    print(f"[embeddings] a calcular: {len(todo)} de {len(sequences)} (resto en cache)")

    tokenizer = AutoTokenizer.from_pretrained(model_name, cache_dir=str(hf_cache_dir()))
    model = AutoModel.from_pretrained(
        model_name,
        cache_dir=str(hf_cache_dir()),
        output_hidden_states=(layer != "last"),
    )
    model.eval()
    model.to(device)

    special_ids = [tokenizer.cls_token_id, tokenizer.eos_token_id, tokenizer.pad_token_id]

    # Ordenar por longitud agrupa secuencias parecidas en el mismo lote y
    # reduce el relleno; en CPU eso se nota.
    ordered = sorted(todo.items(), key=lambda kv: len(kv[1]))

    with torch.no_grad():
        for start in range(0, len(ordered), batch_size):
            chunk = ordered[start : start + batch_size]
            hashes = [h for h, _ in chunk]
            seqs = [s for _, s in chunk]

            batch = tokenizer(seqs, return_tensors="pt", padding=True, truncation=False)
            batch = {k: v.to(device) for k, v in batch.items()}
            out = model(**batch)

            hidden = out.last_hidden_state if layer == "last" else out.hidden_states[int(layer)]

            if pooling == "mean":
                pooled = masked_mean_pool(
                    hidden, batch["attention_mask"], batch["input_ids"], special_ids
                )
            elif pooling == "mean_max":
                mean = masked_mean_pool(
                    hidden, batch["attention_mask"], batch["input_ids"], special_ids
                )
                mx = masked_max_pool(
                    hidden, batch["attention_mask"], batch["input_ids"], special_ids
                )
                pooled = torch.cat([mean, mx], dim=-1)
            else:
                raise ValueError(f"Agregacion desconocida: {pooling!r} (usa mean o mean_max)")

            for h, vec in zip(hashes, pooled.cpu().numpy()):
                cache[h] = vec.astype(np.float32)

            done = min(start + batch_size, len(ordered))
            if done % (batch_size * 20) == 0 or done == len(ordered):
                print(f"[embeddings]   {done}/{len(ordered)}")

    if use_cache:
        save_cache(path, cache)
        print(f"[embeddings] cache en {path} ({path.stat().st_size / 1e6:.1f} MB)")

    return {h: cache[h] for h in sequences}


def embedding_matrix(
    df,
    model_name: str,
    layer: str | int = "last",
    pooling: str = "mean",
    batch_size: int = 16,
    device: str = "cpu",
) -> np.ndarray:
    """Matriz (n_filas_de_df, dim), en el orden de las filas del DataFrame.

    El cache se indexa por hash de secuencia, no por fila, asi que una secuencia
    presente en varias fuentes se calcula una sola vez y aqui se repite su
    vector donde haga falta.
    """
    unique = dict(zip(df.seq_sha1, df.sequence))
    vectors = compute_embeddings(unique, model_name, layer, pooling, batch_size, device)
    return np.vstack([vectors[h] for h in df.seq_sha1])
