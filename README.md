# Clasificador de péptidos antimicrobianos con ESM-2

[English](README.en.md) · [Documentación completa](docs/00-indice.md)

Clasificador de péptidos antimicrobianos (AMP vs no-AMP) que compara embeddings
de un modelo de lenguaje de proteínas (ESM-2) contra descriptores
fisicoquímicos, con división de datos por clústeres de identidad de secuencia y
**la fuga residual medida y declarada** en lugar de dada por eliminada.

> **Estado: tubería completa corrida (2026-10-17).** 18 644 péptidos, 16 modelos
> evaluados, todas las cifras de este README salen de `results/metrics/`.
> El resultado principal es un empate: la diferencia entre ESM-2 y el baseline
> fisicoquímico no se distingue del cero. La división lleva fuga por homología
> **medida y declarada**, no eliminada.

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

**1. Una división por clústeres con la fuga medida, no supuesta.** Los AMPs
vienen en familias de secuencias casi idénticas: la magainina 1 y la 2 difieren
en dos residuos. Una división aleatoria pone parientes cercanos en entrenamiento
y en prueba, y entonces la métrica mide memorización. Aquí se agrupan todas las
secuencias al 40% de identidad y se asignan **clústeres completos** a
train/val/test.

Y aquí está lo que de verdad distingue a este repositorio: **agrupar no basta, y
el proyecto lo demuestra sobre sí mismo.** Una auditoría independiente mide la
identidad máxima de cada secuencia de prueba contra *todo* entrenamiento, contra
un nulo por barajado, y encuentra fuga residual significativa: exceso de +68
secuencias marcadas sobre las 24 esperadas por azar en test, +62 en validación.
Se concentra en péptidos de menos de 15 residuos, donde el prefiltro de k-meros
de MMseqs2 falla por construcción. Entre el 74% y el 86% de lo marcado tiene
identidad ≥ 0,60 con un pariente en entrenamiento.

Ese número **no está escondido ni redondeado hacia abajo**: va dentro de cada
JSON de métricas y en el pie de cada figura, para que no se pueda separar del
resultado ([D-22](docs/08-decisiones.md)). Lo que se publica no es "métricas
limpias", es métricas con su defecto medido al lado, y la curva de rendimiento
frente al umbral de identidad ([M-01](docs/10-mejoras-futuras.md)) cuantifica
cuánto importa.

La literatura de predicción de AMPs reporta de forma habitual divisiones
"sin fuga" sin ningún instrumento que lo compruebe. La diferencia de este
proyecto no es tener menos fuga: es tener el número.

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
| [DRAMP 3.0](https://dramp.cpu-bioinfor.org/) | Positivos | 12 784 | 2,8 MB | Académica, no comercial |
| [UniProt](https://www.uniprot.org/) KW-0929 | Positivos | 1 585 | ~0,2 MB | CC BY 4.0 |
| UniProt sin KW-0929, revisados | Negativos | 11 982 | ~1,5 MB | CC BY 4.0 |

Tamaños y conteos medidos el 2026-10-11. Total ≈ 5 MB. Tras limpiar quedan
**18 644 péptidos** (10 594 positivos, 8 050 negativos).

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
| **Baselines** | Composición de aminoácidos (20 features) y 11 descriptores fisicoquímicos (carga neta por Henderson-Hasselbalch, GRAVY, momento hidrofóbico de Eisenberg, pI, índice de Boman…) con regresión logística, random forest y XGBoost. Se entrena además una **ablación sin `length` ni `charge_density`**: los negativos son 7,8 residuos más largos de media y el emparejado de longitudes no puede corregirlo, así que la diferencia entre el baseline completo y la ablación mide cuánto de él era longitud. |
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

Corrida del 2026-10-17 sobre la división canónica. **Toda métrica de esta tabla
se lee junto a la fuga declarada arriba**, que va escrita dentro de cada JSON de
métricas y en el pie de cada figura.

### La comparación preregistrada, que es la única inferencia del proyecto

| Modelo | ROC-AUC | PR-AUC | MCC | F1 |
|---|---|---|---|---|
| Composición + fisicoquímicos, XGBoost | 0,909 [0,865–0,943] | 0,931 | 0,637 [0,540–0,716] | 0,827 |
| ESM-2 35M + logística | 0,895 [0,838–0,940] | 0,909 | 0,645 [0,543–0,727] | 0,841 |

**Diferencia de MCC: +0,008, IC 95% [−0,033, +0,051]. El intervalo incluye el
cero: los dos modelos no se distinguen.**

Esto es exactamente lo que el compromiso previo de abajo decía esperar, y algo
más conservador todavía: no es que los embeddings ganen por poco, es que con
estos datos y esta división **no se puede afirmar que ganen**. Un modelo de
lenguaje de proteínas de 35 millones de parámetros empata con 31 números que se
calculan con una tabla de aminoácidos.

No se reescribe el titular para que suene mejor. Si el proyecto solo publicara
la fila que más conviene, no serviría para nada de lo que dice servir.

### Todos los modelos, como tabla descriptiva

Ordenados por MCC. **No son inferencia**: elegir el mejor de 16 en prueba y
presentar su diferencia sería el error que la comparación preregistrada existe
para evitar (P-05).

| Modelo | ROC-AUC | MCC | F1 |
|---|---|---|---|
| ESM-2 35M + MLP | 0,913 [0,862–0,950] | 0,659 | 0,852 |
| ESM-2 35M + logística | 0,895 [0,838–0,940] | 0,645 | 0,841 |
| comp+fisicoquímicos, XGBoost | 0,909 [0,865–0,943] | 0,637 | 0,827 |
| comp+fisicoquímicos, random forest | 0,903 [0,852–0,940] | 0,634 | 0,848 |
| ESM-2 8M + logística | 0,902 [0,860–0,934] | 0,631 | 0,839 |
| comp+fisicoq. sin longitud, XGBoost | 0,907 [0,862–0,941] | 0,626 | 0,818 |
| comp+fisicoq. sin longitud, random forest | 0,901 [0,850–0,938] | 0,620 | 0,846 |
| ESM-2 8M + MLP | 0,891 [0,824–0,938] | 0,614 | 0,837 |
| composición, random forest | 0,891 [0,845–0,926] | 0,598 | 0,835 |
| composición, XGBoost | 0,872 [0,819–0,913] | 0,594 | 0,835 |
| fisicoquímicos, random forest | 0,873 [0,819–0,916] | 0,559 | 0,801 |
| fisicoquímicos, XGBoost | 0,870 [0,817–0,914] | 0,534 | 0,816 |
| comp+fisicoquímicos, logística | 0,821 [0,749–0,877] | 0,481 | 0,786 |
| comp+fisicoq. sin longitud, logística | 0,795 [0,733–0,846] | 0,452 | 0,741 |
| composición, logística | 0,779 [0,720–0,829] | 0,424 | 0,666 |
| fisicoquímicos, logística | 0,755 [0,661–0,838] | 0,377 | 0,567 |

### Tres cosas que esta tabla dice y conviene no pasar por alto

**El baseline no era débil.** 0,909 de ROC-AUC con 31 descriptores
fisicoquímicos. La advertencia escrita en `docs/01` —que esos números son
difíciles de superar porque la actividad antimicrobiana depende sobre todo de
carga y hidrofobicidad— se cumplió.

**La longitud no era el truco.** La ablación `sin longitud` (D-16) baja el MCC
de 0,637 a 0,626, dentro del ruido. El baseline no ganaba midiendo cuán largo
es el péptido, así que la comparación contra ESM-2 es más limpia de lo que la
asimetría de esa ablación hacía temer. Esto sustituye por un número una
advertencia que hasta ahora estaba solo escrita.

**Los intervalos son anchos, y eso es honesto, no un defecto.** El bootstrap
remuestrea **clústeres completos**, no secuencias. Con secuencias daría
intervalos mucho más estrechos y falsos, porque las secuencias de un clúster son
parientes y no son observaciones independientes.

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

## Lo que viene: una fase de generación

Hay una segunda fase diseñada y todavía sin implementar, que cierra el círculo
con el proyecto de 2023: usar el clasificador como juez para que un algoritmo
**diseñe** péptidos nuevos, comparando un algoritmo genético clásico contra una
evolución guiada por ESM-2 —la misma pregunta de la primera mitad, hecha otra
vez.

Con una trampa que es justo lo interesante: optimizar contra tu propio
clasificador no es diseño, es un ataque adversarial. Un algoritmo que maximice
su puntuación converge en poli-arginina, que puntúa altísimo y es tóxica e
inespecífica, porque esa región del espacio no tiene negativos en ninguna base
de datos de AMPs. Así que el resultado declarado de esa fase no será "diseñé
antimicrobianos nuevos" sino **cuánto engaña el generador al clasificador, y qué
guardas hacen falta**. Diseño completo en
[docs/17-diseno-generacion.md](docs/17-diseno-generacion.md).

## Reproducir

Requiere [uv](https://docs.astral.sh/uv/). Python 3.12 (no 3.14: PyTorch aún no
publica ruedas estables para esa versión).

```bash
git clone https://github.com/rodrigohermoza/amp-classifier-esm2.git
cd amp-classifier-esm2
make setup
make all
```

> El entorno virtual se crea en `~/.venvs/ampesm` y `.venv` queda como un enlace
> a esa ruta. No es capricho: en la maquina de desarrollo el proyecto vive en una
> carpeta sincronizada por iCloud, que renombraba archivos dentro del entorno y
> lo dejaba inservible ([P-28](docs/09-problemas-encontrados.md), [D-19](docs/08-decisiones.md)).
> La ruta es un `?=`: exporta `UV_PROJECT_ENVIRONMENT` para cambiarla.

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
tests/             129 casos; los descriptores se validan contra literatura
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
