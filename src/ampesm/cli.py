"""Interfaz de linea de comandos: un subcomando por etapa de la tuberia.

Las etapas pasan informacion por disco, no en memoria, asi que cada subcomando
se puede correr por separado y la tuberia se puede reanudar en cualquier punto.
Ver docs/02-arquitectura.md.
"""

from __future__ import annotations

import argparse
import sys
import time


def _cmd_download(args: argparse.Namespace) -> None:
    from .data.download import download_all

    download_all(keys=args.sources, force=args.force)


def _cmd_data(args: argparse.Namespace) -> None:
    from .data.clean import build_dataset
    from .data.download import download_all

    if not args.skip_download:
        download_all(keys=args.sources, force=args.force)
    build_dataset(keys=args.sources)


def _cmd_split(args: argparse.Namespace) -> None:
    from .split import build_splits

    build_splits(force_method=args.method)


def _cmd_baselines(_args: argparse.Namespace) -> None:
    from .models import run_baselines

    run_baselines()


def _cmd_embeddings(_args: argparse.Namespace) -> None:
    """Precalcula y cachea embeddings sin entrenar nada.

    Util para separar la parte lenta: se deja calculando y despues el
    entrenamiento es casi instantaneo porque todo sale del cache.
    """
    from .config import load_params
    from .embeddings import compute_embeddings
    from .models import load_splits

    params = load_params()
    df = load_splits()
    ecfg = params.embeddings
    unique = dict(zip(df.seq_sha1, df.sequence))
    print(f"[embeddings] {len(unique)} secuencias unicas de {len(df)} filas")

    for model_name in ecfg.get("models", []):
        compute_embeddings(
            unique,
            model_name,
            layer=ecfg.get("layer", "last"),
            pooling=ecfg.get("pooling", "mean"),
            batch_size=int(ecfg.get("batch_size", 16)),
            device=str(ecfg.get("device", "cpu")),
        )


def _cmd_train(_args: argparse.Namespace) -> None:
    from .models import run_embedding_models

    run_embedding_models()


def _cmd_evaluate(_args: argparse.Namespace) -> None:
    from .evaluate import evaluate_all

    evaluate_all()


def _cmd_interpret(args: argparse.Namespace) -> None:
    from .interpret import interpret_all

    interpret_all(tag=args.model)


def _cmd_all(args: argparse.Namespace) -> None:
    _cmd_data(args)
    _cmd_split(args)
    _cmd_baselines(args)
    _cmd_train(args)
    _cmd_evaluate(args)
    _cmd_interpret(args)


def _cmd_info(_args: argparse.Namespace) -> None:
    """Estado de la tuberia: que artefactos existen y que falta."""
    from .config import (
        CLEAN_REPORT_PATH,
        MANIFEST_PATH,
        PEPTIDES_PATH,
        SPLIT_REPORT_PATH,
        SPLITS_PATH,
    )
    from .data.sources import SOURCES

    print("Fuentes configuradas:")
    for src in SOURCES:
        rol = "positivos" if src.is_positive else "negativos"
        print(f"  {src.key:14s} {rol:10s} {src.license}")

    print("\nArtefactos:")
    for label, path in (
        ("manifiesto", MANIFEST_PATH),
        ("peptidos", PEPTIDES_PATH),
        ("reporte limpieza", CLEAN_REPORT_PATH),
        ("splits", SPLITS_PATH),
        ("reporte split", SPLIT_REPORT_PATH),
    ):
        if path.exists():
            print(f"  [ok]    {label:18s} {path.stat().st_size:>10,} B  {path}")
        else:
            print(f"  [falta] {label:18s} {'':>10}     {path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ampesm",
        description="Clasificador de peptidos antimicrobianos con embeddings de ESM-2.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def add_source_args(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "--sources", nargs="*", default=None,
            help="Claves de fuentes a usar (por defecto, todas).",
        )
        p.add_argument(
            "--force", action="store_true",
            help="Volver a descargar aunque el archivo ya exista.",
        )

    p = sub.add_parser("download", help="Solo descargar las fuentes crudas.")
    add_source_args(p)
    p.set_defaults(func=_cmd_download)

    p = sub.add_parser("data", help="Descargar, limpiar y escribir peptides.csv.")
    add_source_args(p)
    p.add_argument(
        "--skip-download", action="store_true",
        help="Usar los archivos que ya estan en data/raw.",
    )
    p.set_defaults(func=_cmd_data)

    p = sub.add_parser("split", help="Agrupar por identidad y dividir train/val/test.")
    p.add_argument(
        "--method", choices=["mmseqs2", "cdhit", "builtin"], default=None,
        help="Forzar un metodo de clustering en vez de probarlos en cascada.",
    )
    p.set_defaults(func=_cmd_split)

    p = sub.add_parser("baselines", help="Entrenar los baselines tabulares.")
    p.set_defaults(func=_cmd_baselines)

    p = sub.add_parser("embeddings", help="Precalcular y cachear embeddings de ESM-2.")
    p.set_defaults(func=_cmd_embeddings)

    p = sub.add_parser("train", help="Entrenar clasificadores sobre embeddings.")
    p.set_defaults(func=_cmd_train)

    p = sub.add_parser("evaluate", help="Metricas en prueba con IC bootstrap y figuras.")
    p.set_defaults(func=_cmd_evaluate)

    p = sub.add_parser("interpret", help="UMAP y analisis de errores.")
    p.add_argument("--model", default=None, help="Etiqueta del modelo (por defecto, el mejor).")
    p.set_defaults(func=_cmd_interpret)

    p = sub.add_parser("all", help="Correr la tuberia completa.")
    add_source_args(p)
    p.add_argument("--skip-download", action="store_true")
    p.add_argument("--method", choices=["mmseqs2", "cdhit", "builtin"], default=None)
    p.add_argument("--model", default=None)
    p.set_defaults(func=_cmd_all)

    p = sub.add_parser("info", help="Estado de la tuberia y fuentes configuradas.")
    p.set_defaults(func=_cmd_info)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    start = time.perf_counter()
    try:
        args.func(args)
    except (FileNotFoundError, RuntimeError, ValueError, KeyError) as exc:
        # Los errores previstos se muestran como mensaje, no como traza: casi
        # siempre son "falta una etapa anterior" o "una fuente cambio".
        print(f"\nError: {exc}", file=sys.stderr)
        return 1
    print(f"\n[{args.command}] listo en {time.perf_counter() - start:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
