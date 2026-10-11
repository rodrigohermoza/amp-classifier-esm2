# Carpeta de datos

Esta carpeta está **vacía en el repositorio**, a propósito.

## Por qué no hay datos aquí

APD3 y DRAMP piden citación y restringen el uso comercial, y ninguno declara una
licencia abierta. Redistribuir sus archivos desde un repositorio público con
licencia MIT sería, como mínimo, discutible. Así que el repositorio incluye los
scripts de descarga y los manifiestos, no los datos.

Efecto secundario conveniente: el repositorio se clona en segundos.

## Cómo reconstruir los datos

```bash
make data
```

Descarga unos 5 MB desde las fuentes originales y escribe
`data/processed/peptides.csv`.

## Qué va en cada subcarpeta

| Carpeta | Contenido | ¿Se versiona? |
|---|---|---|
| `raw/` | Archivos tal como los devuelve cada fuente | No, salvo `manifest.json` |
| `interim/` | Intermedios de depuración | No |
| `processed/` | `peptides.csv`, `splits.csv` y los reportes | No los CSV; sí los reportes |

Los archivos que **sí** se versionan:

- `raw/manifest.json` — URL, fecha, tamaño y `sha256` de cada descarga. Es lo
  que permite verificar que descargaste el mismo conjunto que se usó para los
  números publicados, y detectar que una fuente cambió bajo nuestros pies.
- `processed/clean_report.json` — conteos de cada filtro, duplicados y
  conflictos de etiqueta.
- `processed/split_report.json` — método de clustering, su versión, tamaños de
  clúster y balance de clases por split.

## Si una fuente cambia

El `sha256` del manifiesto no coincidirá con el de la nueva descarga. Eso no es
un error de la tubería: es la señal de que la base de datos se actualizó. Lo que
hay que hacer es anotarlo en
[docs/09-problemas-encontrados.md](../docs/09-problemas-encontrados.md) y decidir
si se vuelven a generar los resultados con la versión nueva.

Detalle de cada fuente, con su licencia y cita, en
[docs/03-datos-fuentes.md](../docs/03-datos-fuentes.md).
