# Arquitectura: como se conecta todo

## Idea general

El proyecto es una tuberia lineal de siete etapas. Cada etapa lee archivos de la
anterior y escribe archivos propios. No hay estado en memoria entre etapas: todo
pasa por disco. Esto hace que se pueda reanudar en cualquier punto y que cada
artefacto intermedio sea inspeccionable.

```
  fuentes web
      |
      v
 [1] descarga  ------> data/raw/*.fasta|.xlsx  + manifest.json (URL, fecha, sha256)
      |
      v
 [2] limpieza  ------> data/processed/peptides.csv        (secuencia, etiqueta, fuente)
      |
      v
 [3] clusteres ------> data/processed/splits.csv          (+ cluster_id, split)
      |
      +----------------------------+
      |                            |
      v                            v
 [4] features              [5] embeddings ESM-2
      |                            |
      v                            v
 [6] modelos  <---------------------+
      |
      v
 [7] evaluacion ------> results/metrics/*.json, results/figures/*.png
```

## Mapa de modulos

| Modulo | Responsabilidad | No hace |
|---|---|---|
| `ampesm/config.py` | Rutas, semilla, constantes (alfabeto, longitudes), carga de `params.yaml` | Logica cientifica |
| `ampesm/data/sources.py` | Registro declarativo de fuentes: URL, licencia, tipo de parseo | Red |
| `ampesm/data/download.py` | Descarga con reintentos, calcula sha256, escribe el manifiesto | Parsear |
| `ampesm/data/parse.py` | FASTA y XLSX -> lista de `PeptideRecord` | Filtrar |
| `ampesm/data/clean.py` | Filtros de alfabeto y longitud, deduplicado, emparejado de longitud | Descargar |
| `ampesm/split.py` | Clustering por identidad, asignacion de clusteres a train/val/test y las auditorias de fuga | Entrenar |
| `ampesm/features/composition.py` | Frecuencia de los 20 aminoacidos | Fisicoquimica |
| `ampesm/features/physchem.py` | Carga neta, GRAVY, momento hidrofobico, pI, aromaticidad | Composicion |
| `ampesm/embeddings.py` | Carga ESM-2, agrega por media de residuos, cachea en disco | Clasificar |
| `ampesm/models.py` | Definicion y ajuste de LR / RF / XGBoost / MLP, y las 4 representaciones tabulares | Metricas |
| `ampesm/evaluate.py` | Metricas, bootstrap, curvas, matriz de confusion | Entrenar |
| `ampesm/interpret.py` | UMAP, analisis de errores | Metricas |
| `ampesm/cli.py` | Un subcomando por etapa | Todo lo demas |

Regla que se sigue en todo el codigo: **los modulos de datos no importan modulos
de modelos, y ningun modulo llama a la red salvo `download.py`.** Asi los tests
corren sin internet.

## Por que pasa todo por disco

Tres razones concretas:

1. Los embeddings de ESM-2 son lo mas caro de calcular (minutos en CPU para miles
   de peptidos). Cachearlos en `.npz` indexado por hash de secuencia significa
   que reentrenar un clasificador toma segundos y no minutos.
2. Si la tuberia falla en la etapa 6, no hay que volver a descargar nada.
3. Un archivo intermedio se puede abrir y mirar. Un objeto en memoria no.

## Contrato de datos entre etapas

`data/processed/peptides.csv`:

| columna | tipo | descripcion |
|---|---|---|
| `id` | str | Identificador estable: `{fuente}:{id_original}` |
| `sequence` | str | Secuencia en mayusculas, solo los 20 aminoacidos estandar |
| `length` | int | Longitud, entre 5 y 50 |
| `label` | int | 1 = AMP, 0 = no-AMP |
| `source` | str | `apd3`, `dramp`, `uniprot_amp`, `uniprot_neg` |
| `seq_sha1` | str | Hash de la secuencia, usado como clave del cache de embeddings |

`data/processed/splits.csv` anade:

| columna | tipo | descripcion |
|---|---|---|
| `cluster_id` | int | Cluster de identidad de secuencia al 40% |
| `split` | str | `train`, `val` o `test` |

**Invariante critico:** ningun `cluster_id` aparece en mas de un `split`. Hay un
test que lo verifica (`tests/test_split.py`). Si ese test falla, todas las
metricas del proyecto son invalidas.

`results/models/*.json`, por modelo entrenado:

| campo | descripcion |
|---|---|
| `tag` | `{representacion}__{modelo}`, la clave que cruza con las predicciones |
| `representation` | `composition`, `physchem`, `comp_physchem`, `comp_physchem_nolen` o un modelo de ESM-2 |
| `best_params` | hiperparametros elegidos **en validacion** |
| `val_roc_auc` | lo que usan `interpret` y el fallback de `evaluate` para elegir modelo |
| `features_hash` | huella de la seccion `features` de `params.yaml` con la que se entreno (P-24) |

## Guardas y diagnosticos por etapa

Esta tabla existe porque la mitad del codigo del proyecto no calcula nada: vigila
que lo calculado signifique algo. Cada guarda nacio de un fallo real que no se
veia.

| Etapa | Guarda | Que pasa si salta | Origen |
|---|---|---|---|
| todas (via `make`) | `check-env`: el paquete se importa | aborta con el comando de arreglo | P-25 |
| `data` | cabeceras vacias de APD3 != 1 | aviso: el formato de origen cambio | P-01 |
| `data` | conjunto vacio o de una sola clase | aborta | — |
| `split` | ningun cluster en dos splits | aborta, no escribe nada | invariante |
| `split` | split degenerado (falta una clase, <20 positivos, <10 clusteres) | aborta | P-06 |
| `split` | `leakage.global_p_value` en test **y en val** | aviso ruidoso | P-04, P-27 |
| `split` | `positive_rate_drift` entre splits | aviso | P-26 |
| `split` | `length_ks_pvalue` por split | diagnostico, se reporta | P-10 |
| `split` | `publishable: false` si se uso el respaldo | aviso ruidoso | D-11, P-03 |
| `split` | borra `predictions/*.csv` rancias | silencioso, a proposito | P-20 |
| `baselines` | XGBoost no disponible (`libomp`) | aviso y sigue sin el | P-09 |
| `evaluate` | ids de predicciones contra `splits.csv` | aborta | P-20 |
| `evaluate` | validacion vacia o de una clase en `choose_threshold` | aborta | P-07 |
| `evaluate` | `n_boot_valid / n_boot_requested` < 99% | aviso | P-16 |
| `evaluate` | no-finitos a `null`, `allow_nan=False` | aborta si quedara un NaN | P-21 |
| `interpret` | `features_hash` distinto del guardado | **aborta** | P-24, D-18 |

El criterio de cuando abortar y cuando avisar esta en D-18: abortar cuando la
salida es un artefacto que viaja solo (una figura, un JSON), avisar cuando el
numero sigue siendo util si se conoce el problema.

## Reproducibilidad

- Dependencias con version exacta en `pyproject.toml`, resueltas por `uv.lock`.
- Python fijado a 3.12 en `.python-version`.
- Una sola semilla, `config.SEED = 42`, propagada explicitamente a cada funcion
  que use aleatoriedad. No se confia en el estado global de numpy.
- El mapeo secuencia -> cluster -> split se guarda como artefacto versionado, de
  modo que la division sea identica aunque la herramienta de clustering cambie de
  version.

## Uso de la prueba final

`test` se toca **una vez**, al final, en la etapa 7. Toda la seleccion de
hiperparametros ocurre contra `val`. Esto esta escrito aqui porque es facil de
violar sin darse cuenta y porque es la diferencia entre un numero reportable y
uno que no lo es.
