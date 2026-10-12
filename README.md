# Clasificador de péptidos antimicrobianos con ESM-2

[English](README.en.md) · [Documentación completa](docs/00-indice.md)

Clasificador de péptidos antimicrobianos (AMP vs no-AMP) que compara embeddings
de un modelo de lenguaje de proteínas (ESM-2) contra descriptores
fisicoquímicos, con división de datos por clústeres de identidad de secuencia
para evitar fuga de información.

> **Estado: código completo, pendiente de ejecución.**
> Las tablas de resultados están vacías a propósito. Se rellenan corriendo la
> tubería, y los números saldrán de `results/metrics/`, no de una estimación.

## Por qué este proyecto

Un péptido antimicrobiano es una cadena corta de aminoácidos que mata bacterias
u hongos, normalmente rompiendo su membrana. Interesan porque ese mecanismo es
más difícil de esquivar que la inhibición de una enzima concreta, y porque
probar candidatos en el laboratorio es lento y caro: un clasificador que ordene
candidatos antes de sintetizarlos tiene valor práctico.

En 2023, en Sup'Biotech (Francia), participé en un proyecto de diseño de AMPs
con un algoritmo genético ([repositorio del
equipo](https://github.com/hgrandjean/Anti-Fungi-Peptide)). Mi aporte fue la
obtención de datos de DRAMP y UniProt. Esta es mi propia versión, donde la parte
de aprendizaje automático es trabajo mío de principio a fin.

## Lo que este proyecto intenta hacer bien

El objetivo no es un ROC-AUC alto. En esta tarea es fácil llegar a 0.95+ y no
haber aprendido nada sobre actividad antimicrobiana. Lo que se intenta hacer bien
es:

**1. Una división de datos sin fuga de información.** Los AMPs vienen en
familias de secuencias casi idénticas: la magainina 1 y la 2 difieren en dos
residuos. Una división aleatoria pone parientes cercanos en entrenamiento y en
prueba, y entonces la métrica mide memorización. Aquí se agrupan todas las
secuencias al 40% de identidad y se asignan **clústeres completos** a
train/val/test. Es la pieza central del proyecto y, en mi opinión, el error
metodológico más común en la literatura de predicción de AMPs.

**2. Baselines honestos.** La actividad antimicrobiana depende sobre todo de dos
propiedades globales: carga positiva neta y una cara hidrofóbica. Eso hace que 30
números fisicoquímicos sean difíciles de superar. Si el baseline gana, el README
lo dirá en el titular.

**3. Incertidumbre reportada.** Métricas con intervalos de confianza por
bootstrap, remuestreando **clústeres** y no secuencias, porque las secuencias
emparentadas no son independientes y tratarlas como si lo fueran da intervalos
demasiado estrechos. Dos modelos cuyos intervalos se solapan no se distinguen.

**4. La debilidad principal dicha en voz alta.** Abajo.

## La debilidad conocida: los negativos

No existe una base de datos de péptidos que se probaron y resultaron **no** ser
antimicrobianos. Lo que se usa aquí, como en toda la literatura, son péptidos de
UniProt sin anotación antimicrobiana. Eso no significa "no antimicrobiano":
significa **"nadie lo ha probado o nadie lo anotó"**.

Dos consecuencias:

- Hay falsos negativos reales entre los negativos, así que el rendimiento
  medible tiene un techo por debajo del real.
- Peor: el modelo puede aprender a distinguir **procedencia** en vez de
  **actividad**. Los AMPs de APD3 son péptidos secretados y catiónicos, muchos de
  piel de anfibio; los fragmentos de UniProt son otra cosa. Un AUC alto puede
  significar "sé de qué base de datos viene" y no "sé si mata bacterias".

Mitigaciones aplicadas: emparejar la distribución de longitud entre clases,
colorear el UMAP por fuente además de por clase, y medir el rendimiento por
fuente por separado en el análisis de errores. Si el acierto en los positivos de
UniProt es mucho peor que en los de APD3, el modelo depende de la procedencia y
el resultado hay que leerlo con eso en mente.

## Datos

El repositorio **no contiene los datos crudos**: APD3 y DRAMP piden citación y
restringen el uso comercial. Contiene los scripts de descarga y un manifiesto con
la URL, la fecha y el `sha256` de cada archivo. `make data` reconstruye el
conjunto.

| Fuente | Rol | Registros | Tamaño | Licencia |
|---|---|---|---|---|
| [APD3](https://aps.unmc.edu/) | Positivos | 3 306 | 147 KB | Académica, con citación |
| [DRAMP 3.0](https://dramp.cpu-bioinfor.org/) | Positivos | — | 2,8 MB | Académica, no comercial |
| [UniProt](https://www.uniprot.org/) KW-0929 | Positivos | 1 585 | ~0,2 MB | CC BY 4.0 |
| UniProt sin KW-0929, revisados | Negativos | 11 982 | ~1,5 MB | CC BY 4.0 |

Tamaños y conteos medidos el 2026-10-11. Total ≈ 5 MB.

Limpieza: solo los 20 aminoácidos estándar, longitud 5–50, deduplicado por
secuencia exacta y emparejado de longitud. Detalle en
[docs/03-datos-fuentes.md](docs/03-datos-fuentes.md).

**DBAASP no está.** Era la cuarta fuente y la única con valores de MIC, lo que
habilitaba una regresión de potencia. Al 2026-10-11 devuelve 403 o respuestas de
0 bytes en todas las rutas probadas, así que la regresión queda diferida.
Registrado en [docs/09-problemas-encontrados.md](docs/09-problemas-encontrados.md), P-02.

## Método

```
descarga → limpieza → clústeres al 40% → ┬→ descriptores ┬→ modelos → evaluación
                                         └→ ESM-2        ┘
```

| Etapa | Qué hace |
|---|---|
| **División** | MMseqs2 al 40% de identidad (obligatorio para resultados publicables); respaldo en CD-HIT y en un clustering propio en Python para que el repo funcione recién clonado. Clústeres completos a train/val/test (70/15/15), con validación que aborta si algún split queda degenerado. |
| **Auditoría de fuga** | Identidad máxima de cada secuencia de test contra **todo** train, comparada contra una distribución nula (las mismas secuencias con los residuos barajados). Lo que se reporta como fuga es el **exceso sobre el nulo**, porque la identidad entre péptidos cortos tiene un ruido de fondo del ~33%. |
| **Baselines** | Composición de aminoácidos (20 features) y 11 descriptores fisicoquímicos (carga neta por Henderson-Hasselbalch, GRAVY, momento hidrofóbico de Eisenberg, pI, índice de Boman…) con regresión logística, random forest y XGBoost. |
| **Embeddings** | ESM-2 8M y 35M, media por residuo con máscara explícita de relleno y tokens especiales, cacheados en disco. Logística y un MLP pequeño encima. |
| **Evaluación** | ROC-AUC, PR-AUC, MCC y F1 con IC del 95% por bootstrap de clústeres. Umbral elegido en validación. Comparación emparejada entre el par **preregistrado**, nunca entre los ganadores en test. |
| **Interpretación** | UMAP por clase, fuente y split; análisis de errores por categoría, fuente y longitud. |

Las escalas fisicoquímicas están escritas explícitas y citadas en
[`physchem.py`](src/ampesm/features/physchem.py), para que se puedan auditar sin
salir del repositorio y para que los números no cambien si una librería cambia
de versión. Verificadas contra valores publicados: magainina-2 da +3,0 de carga
neta a pH 7,4 (publicado +3 a +4), LL-37 da +5,9 (publicado +6) y el índice de
Boman de LL-37 da 3,0 (publicado 2,9–3,0).

## Resultados

Pendientes de ejecución.

| Modelo | ROC-AUC | PR-AUC | MCC | F1 |
|---|---|---|---|---|
| Composición + fisicoquímicos, XGBoost | — | — | — | — |
| ESM-2 35M + logística | — | — | — | — |

### Compromiso previo

Escrito antes de correr nada, para que esto no sea narrar a posteriori lo que
salió bien:

- Modelo principal: **ESM-2 35M + regresión logística**.
- Baseline principal: **composición + fisicoquímicos con XGBoost**.
- Métrica principal: **MCC** en prueba, con IC por bootstrap de clústeres.
- Prueba se evalúa **una sola vez**, al final. La tabla con todos los modelos
  es descriptiva; la única inferencia del proyecto es la comparación emparejada
  entre ese par, elegido por preregistro y no por su resultado en prueba.
- Mi expectativa: los embeddings deberían ganar, pero por un margen modesto
  (quizá 0,02–0,04 de ROC-AUC), y puede que los intervalos se solapen. Si el
  margen sale enorme, lo primero que hay que sospechar es fuga de información,
  no éxito.

## Reproducir

Requiere [uv](https://docs.astral.sh/uv/). Python 3.12 (no 3.14: PyTorch aún no
publica ruedas estables para esa versión).

```bash
git clone https://github.com/rodrigohermoza/amp-classifier-esm2.git
cd amp-classifier-esm2
make setup
make all
```

Etapa por etapa:

```bash
make data        # descarga y limpia  (~5 MB de red)
make split       # clústeres y división
make baselines   # descriptores + LR/RF/XGBoost
make embeddings  # ESM-2: ~31 MB el modelo de 8M, ~140 MB el de 35M
make train       # clasificadores sobre embeddings
make evaluate    # métricas con IC + figuras
make interpret   # UMAP y análisis de errores
```

`ampesm info` muestra qué artefactos existen y qué falta. Todo corre en CPU en un
Mac con Apple Silicon. Hacen falta dos cosas de Homebrew:
`brew install mmseqs2` (~30 MB, obligatorio para resultados publicables) y
`brew install libomp` (~2 MB, sin ella XGBoost no carga en macOS).

## Estructura

```
src/ampesm/
  config.py        rutas, alfabeto, parámetros
  data/            fuentes, descarga, parseo, limpieza
  split.py         clustering por identidad y división sin fuga
  features/        composición y descriptores fisicoquímicos
  embeddings.py    ESM-2 con cache en disco
  models.py        entrenamiento y selección de hiperparámetros
  evaluate.py      métricas, bootstrap, figuras
  interpret.py     UMAP y análisis de errores
docs/              cuaderno de laboratorio: decisiones, errores, mejoras
tests/             117 casos; los descriptores se validan contra literatura
params.yaml        todo lo ajustable en un solo archivo
```

## Documentación

El [cuaderno de laboratorio](docs/00-indice.md) incluye el registro de
decisiones con su justificación ([08](docs/08-decisiones.md)), los problemas
reales encontrados y cómo se resolvieron ([09](docs/09-problemas-encontrados.md)),
las mejoras pendientes ordenadas por esfuerzo/beneficio
([10](docs/10-mejoras-futuras.md)) y el uso de asistencia de IA con tokens y
tiempo ([11](docs/11-agentes.md)).

Si vas a criticar el método, empieza por [04](docs/04-division-clusters.md) y
[09](docs/09-problemas-encontrados.md): ahí está lo interesante.

## Citas

- Wang G., Li X., Wang Z. *APD3: the antimicrobial peptide database.* NAR, 2016.
- Shi G. et al. *DRAMP 3.0.* NAR, 2022.
- The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase.*
- Lin Z. et al. *Evolutionary-scale prediction of atomic-level protein structure
  with a language model.* Science, 2023. (ESM-2)

## Licencia

Código bajo [MIT](LICENSE). **La licencia no cubre los datos**: cada base tiene
sus propios términos, detallados en [docs/03](docs/03-datos-fuentes.md).
