"""Descarga de las fuentes y escritura del manifiesto.

Este es el unico modulo del proyecto que toca la red. Todo lo demas trabaja
sobre archivos en disco, asi que los tests corren sin internet.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import requests

from ..config import MANIFEST_PATH, RAW_DIR, ensure_dirs
from .sources import QUERY_UNIPROT_NEG, QUERY_UNIPROT_POS, Source, get_sources, uniprot_count_url

TIMEOUT = 120
CHUNK = 1 << 16  # 64 KB
RETRIES = 3

# UniProt es lento en generar un stream grande; conviene un User-Agent que
# identifique al cliente, que es lo que piden sus terminos de uso.
HEADERS = {"User-Agent": "ampesm/0.1 (proyecto academico; contacto via GitHub)"}

QUERIES = {"uniprot_amp": QUERY_UNIPROT_POS, "uniprot_neg": QUERY_UNIPROT_NEG}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def _declared_total(query: str) -> int | None:
    """Cuantos registros dice UniProt que hay, leido de x-total-results.

    Es informativo: se guarda en el manifiesto para poder detectar que una
    consulta devuelve algo distinto en una fecha posterior. Un fallo aqui no
    debe interrumpir la descarga, asi que se traga la excepcion.
    """
    try:
        resp = requests.get(uniprot_count_url(query), headers=HEADERS, timeout=30, stream=True)
        resp.close()
        total = resp.headers.get("x-total-results")
        return int(total) if total is not None else None
    except (requests.RequestException, ValueError):
        return None


def download_one(source: Source, dest_dir: Path | None = None, force: bool = False) -> dict:
    """Descarga una fuente y devuelve su entrada de manifiesto.

    Si el archivo ya existe y `force` es falso, no vuelve a descargar pero si
    recalcula el hash: asi una corrida repetida verifica la integridad de lo que
    hay en disco en vez de confiar en que sigue intacto.
    """
    dest_dir = dest_dir or RAW_DIR
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / source.filename

    reused = dest.exists() and not force
    if not reused:
        _fetch(source.url, dest)

    entry = {
        "source": source.key,
        "url": source.url,
        "filename": source.filename,
        "bytes": dest.stat().st_size,
        "sha256": sha256_file(dest),
        "label": source.label,
        "license": source.license,
        "citation": source.citation,
        "homepage": source.homepage,
        "downloaded_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "reused_existing_file": reused,
    }
    if source.key in QUERIES:
        entry["query"] = QUERIES[source.key]
        entry["declared_total"] = _declared_total(QUERIES[source.key])
    return entry


def _fetch(url: str, dest: Path) -> None:
    """GET con reintentos, escribiendo a un temporal y renombrando al final.

    El temporal evita dejar un archivo truncado en `data/raw/` si la descarga se
    corta a medias: un archivo parcial que parece completo es peor que ninguno,
    porque la siguiente corrida lo reutilizaria.
    """
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error: Exception | None = None

    for attempt in range(1, RETRIES + 1):
        try:
            with requests.get(
                url, headers=HEADERS, timeout=TIMEOUT, stream=True, allow_redirects=True
            ) as resp:
                resp.raise_for_status()
                written = 0
                with open(tmp, "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=CHUNK):
                        if chunk:
                            fh.write(chunk)
                            written += len(chunk)
            if written == 0:
                # Lo que hace DBAASP: 200 con cuerpo vacio. Un archivo de 0
                # bytes pasaria los filtros y daria 0 secuencias sin explicar
                # por que, asi que se trata como error explicito.
                raise RuntimeError(f"Respuesta de 0 bytes desde {url}")
            tmp.replace(dest)
            return
        except (requests.RequestException, RuntimeError) as exc:
            last_error = exc
            tmp.unlink(missing_ok=True)
            if attempt < RETRIES:
                print(f"  intento {attempt}/{RETRIES} fallo ({exc}); reintentando")

    raise RuntimeError(f"No se pudo descargar {url} tras {RETRIES} intentos: {last_error}")


def download_all(keys: list[str] | None = None, force: bool = False) -> dict:
    """Descarga todas las fuentes y escribe data/raw/manifest.json."""
    ensure_dirs()
    sources = get_sources(keys)
    entries = []

    for src in sources:
        print(f"[descarga] {src.key}: {src.filename}")
        entry = download_one(src, force=force)
        estado = "reutilizado" if entry["reused_existing_file"] else "descargado"
        print(f"  {estado}: {entry['bytes']:,} bytes  sha256={entry['sha256'][:12]}...")
        entries.append(entry)

    manifest = {
        "generated_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "total_bytes": sum(e["bytes"] for e in entries),
        "entries": entries,
    }
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[descarga] manifiesto en {MANIFEST_PATH}")
    print(f"[descarga] total: {manifest['total_bytes'] / 1e6:.2f} MB")
    return manifest
