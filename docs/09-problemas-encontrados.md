# Problemas encontrados

Bitacora de lo que fallo de verdad. Se mantiene porque es la parte de un proyecto
que normalmente se borra y la que mas ensena a quien viene detras.

Formato: sintoma, causa, resolucion, leccion.

---

## P-01 — La primera linea del FASTA de APD3 no es una secuencia
**Fecha:** 2026-10-11 · **Estado:** resuelto en el parser

**Sintoma:** el archivo `naturalAMPs_APD2024a.fasta` empieza asi:

```
>Your search led to 3306 peptides
>AP00001
GLWSKIKEVGKEAAKAAAKAAGKAALGAVSEAV
```

**Causa:** la descarga se genera desde la interfaz de busqueda del sitio y
arrastra el encabezado del resultado como si fuera un registro FASTA.

**Efecto si no se trata:** un parser ingenuo crea un registro con id
`Your search led to 3306 peptides` y secuencia vacia, o peor, le asigna la
secuencia del primer peptido real, desplazando todas las etiquetas una posicion.
Un desplazamiento asi es silencioso y envenena el conjunto entero.

**Resolucion:** `parse.py` descarta todo registro cuya secuencia quede vacia
despues del filtro de alfabeto, y registra cuantos descarto. Si ese numero no es
exactamente 1 para APD3, el formato del archivo cambio y hay que mirarlo.

**Leccion:** nunca confiar en que un archivo de una base de datos biologica
cumple su formato. Validar y contar siempre; un conteo que no cuadra con el
esperado es la senal mas barata de que algo se rompio.

---

## P-02 — DBAASP inaccesible por script
**Fecha:** 2026-10-11 · **Estado:** sin resolver, con rodeo

**Sintoma:** todas las rutas de `dbaasp.org` fallan:

| Ruta | Resultado |
|---|---|
| `https://dbaasp.org/` | HTTP 403 |
| `https://dbaasp.org/home` | HTTP 403 |
| `https://dbaasp.org/api/v1/peptides?limit=1` | HTTP 200, cuerpo de 0 bytes |
| `https://dbaasp.org/api/v1/peptides/all` | HTTP 200, cuerpo de 0 bytes |
| `https://www.dbaasp.org/api/v1/peptides?limit=1` | HTTP 200, cuerpo de 0 bytes |

Se probo con y sin cabecera `User-Agent` de navegador, y con
`Accept: application/json`. Mismo resultado.

**Causa probable:** proteccion anti-bot (tipo Cloudflare) que exige ejecucion de
JavaScript, o filtrado a nivel de red en el entorno de desarrollo. El patron
"200 con cuerpo vacio" apunta mas a lo primero.

**Efecto:** se pierde la unica fuente con valores de MIC, y con ella la regresion
de potencia.

**Rodeo adoptado:** entregar la clasificacion binaria con APD3 + DRAMP + UniProt,
que no necesita DBAASP. Ver [08-decisiones.md](08-decisiones.md), D-07.

**Vias de recuperacion sin resolver:**
1. Descargar a mano desde un navegador y colocar el archivo en `data/raw/`; el
   parser lo aceptaria. Rompe la reproducibilidad de un solo comando, asi que
   seria un camino opcional documentado, no el principal.
2. Buscar un espejo de DBAASP en Zenodo o en el material suplementario de algun
   articulo.
3. Usar los valores de MIC que DRAMP trae en sus metadatos. Son menos y menos
   estructurados, pero existen y vienen en un archivo al que si se accede.

**Leccion:** verificar el acceso real a cada fuente **antes** de escribir el plan
que depende de ella. Comprobar las cuatro fuentes al principio costo unos minutos
y cambio el alcance del proyecto; descubrirlo en la etapa 5 habria costado dias.

---

## P-03 — El prefiltro de k-meros anulaba el umbral de identidad del 40%
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Sintoma:** ninguno visible. Es el peor tipo de error: todo pasa, los tests
pasan, y las metricas saldrian infladas sin ningun aviso.

**Causa:** en `cluster_builtin` (`split.py`), un representante solo se alinea si
comparte al menos el 15% de sus 3-meros (indice de Jaccard). Ese prefiltro se
documento como una optimizacion de coste, pero resulta ser **mucho mas
restrictivo que el propio umbral de identidad**, asi que cambia el resultado y
no solo el tiempo.

Medido de forma independiente por dos auditorias, con pares sinteticos:

| Identidad real | % de pares que llega a alinearse |
|---|---|
| 0.43 – 0.45 | **0 – 1%** |
| 0.52 | 0% |
| 0.63 | 20% |
| 0.72 | 87% |
| 0.81 | 100% |

El umbral **efectivo** del respaldo estaba alrededor de 0.75, no de 0.40. Un par
al 45% de identidad practicamente nunca se agrupaba, con lo que sus miembros
podian caer en train y en test.

Caso aun peor, y es justo el que el docstring de `identity_and_coverage` dice
querer detectar: un peptido corto contenido **entero** dentro de uno largo.

| Peptido corto dentro de un 50-mero | Jaccard de 3-meros | ¿Se agrupaba? | Identidad real |
|---|---|---|---|
| 5 residuos | 0.062 | **No** | 1.00 |
| 8 residuos | 0.125 | **No** | 1.00 |
| 10 residuos | 0.167 | Si | 1.00 |

Con `min_length: 5`, todo peptido de 8 residuos o menos identico a un fragmento
de otro acababa en un cluster distinto.

**Por que los tests no lo detectaron:** `test_no_cluster_spans_splits` verifica
el invariante de la **asignacion** (ningun cluster en dos splits), no el de la
**agrupacion**. Si dos parientes quedan en clusteres distintos, cada uno en su
split, no se viola ninguna regla. El test era correcto y aun asi la premisa del
proyecto estaba rota.

**Correccion pendiente:**
1. Usar **contencion** en vez de Jaccard, `|A∩B| / min(|A|,|B|)`, coherente con
   el resto del modulo, que ya normaliza por la secuencia mas corta. Admite el
   caso corto-dentro-de-largo sin reventar el coste.
2. Bajar el umbral del prefiltro de 0.15 a 0.05. Medido con 8 000 secuencias:
   2,4 s con 0.15 y 6 459 clusteres; 5,2 s con 0.05 y 3 592 clusteres. Casi el
   doble de agrupacion por el doble de unos pocos segundos.
3. Degradar el respaldo a "solo pruebas y CI": exigir MMseqs2 para cualquier
   numero publicable, con un aviso ruidoso y el campo `method` del reporte
   marcado como no publicable.
4. Anadir la verificacion que falta y que vale mas que cualquier ajuste de
   parametros (ver P-04).

**Resuelto asi:** los cuatro puntos. El prefiltro usa contencion con umbral
0.05; `build_splits` avisa de forma ruidosa y marca `publishable: false` en el
reporte cuando se usa el respaldo; se anadieron `-s 7.5` y
`--alignment-mode 3` a MMseqs2 por el mismo motivo (su prefiltro por defecto
tambien pierde pares a 40% en peptidos cortos); y se implemento la auditoria de
P-04. Tests de regresion en `TestPrefilterRegression`, incluido el 6-mero
contenido en un 50-mero.

**Leccion:** un test puede verificar correctamente un invariante y aun asi dejar
pasar el fallo que ese invariante pretendia prevenir, porque el invariante estaba
puesto en el eslabon equivocado de la cadena. Documentar un parametro como
"optimizacion de coste" sin medir su efecto en el resultado fue el error de
fondo: el prefiltro no era una optimizacion, era una segunda regla de agrupacion
escondida y mas estricta que la declarada.

---

## P-04 — Falta la verificacion directa de la premisa del proyecto
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Sintoma:** `docs/04` afirma que "ninguna secuencia de prueba tiene un pariente
cercano en entrenamiento". Nada en el codigo lo comprueba.

**Causa:** el clustering voraz por representante, **con cualquiera de las tres
herramientas** (MMseqs2 y CD-HIT incluidos), no garantiza esa afirmacion. Dos
secuencias A y B al 50% de identidad entre si pueden quedar en clusteres
distintos si A se unio al representante R1 y B a R2, sin que ninguna supere el
umbral contra el representante de la otra. Es el comportamiento conocido de los
algoritmos tipo CD-HIT; no es un bug, es una propiedad.

**Correccion pendiente:** tras asignar los splits, calcular para cada secuencia
de `test` la identidad maxima contra **todo** `train` (el indice invertido de
k-meros ya existe) y guardar en `split_report.json` la distribucion y el numero
de secuencias de test con identidad superior al umbral. Opcionalmente, mover esas
secuencias a train.

**Por que importa mas que el propio clustering:** convierte la premisa central
del proyecto de una afirmacion en un numero auditable y publicable. Es la
diferencia entre "dividimos por clusteres, asi que no hay fuga" y "el 0,3% de
las secuencias de test tiene un pariente por encima del 40% en train, y aqui
esta el histograma".

**Resuelto asi:** `audit_train_test_identity` en `split.py`, que corre dentro de
`build_splits` y escribe su resultado en `split_report.json`: distribucion de
identidades maximas, histograma, lista de casos por encima del umbral y la
distribucion nula de P-12. Al correrla sobre datos sinteticos aparecio
inmediatamente el problema de interpretacion que documenta P-12, asi que la
auditoria ya demostro su valor antes de tocar datos reales.

---

## P-05 — El modelo de la comparacion principal se elegia por su resultado en test
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Sintoma:** ninguno visible.

**Causa:** `evaluate_all` ordena la tabla resumen por MCC **de test** y luego
toma `iloc[0]` de los modelos con embeddings y de los baselines para la
comparacion emparejada. El mejor de 13 modelos, elegido en test y evaluado en
test. El intervalo de confianza del bootstrap **no corrige** ese sesgo de
seleccion: la diferencia puntual queda inflada y el intervalo demasiado estrecho
en la direccion del ganador. El mismo error en `interpret_all`, donde el modelo
que se interpreta tambien se elige por test.

Contradice dos compromisos escritos en este mismo repositorio: el de `docs/02`
("test se toca una vez") y el compromiso previo de `docs/07`, que ya
**preregistra** el par a comparar.

**Correccion pendiente:** elegir por validacion. El dato ya esta guardado:
`results/models/{tag}.json` contiene `val_roc_auc`. Mejor aun, fijar el par a
las dos etiquetas preregistradas y dejar el resto de la tabla como descriptiva,
diciendolo en el README.

**Resuelto asi:** `_select_comparison_pair` en `evaluate.py` elige por este
orden: (1) las dos etiquetas preregistradas en docs/07 si existen entre las
predicciones, (2) el mejor de cada familia por ROC-AUC de **validacion**. El
JSON de la comparacion incluye `selected_on` para que quede por escrito cual se
uso. Lo mismo en `interpret.py`, que elige el modelo a interpretar leyendo los
`val_roc_auc` de `results/models/*.json`. Verificado en la corrida sintetica:
eligio el par preregistrado aunque no era el mejor en test.

**Leccion:** escribir el compromiso previo no basta si el codigo que produce el
numero no lo implementa. El preregistro y la implementacion tienen que apuntar
al mismo sitio, y aqui no lo hacian.

---

## P-06 — `assign_splits` podia dejar un split vacio o de una sola clase, en silencio
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Causa:** la heuristica arranca con deficit 1.0 en los tres splits y pondera el
balance de clases con peso 0.5, asi que los primeros clusteres grandes se
reparten uno por split sin que el balance pueda vencer al deficit.

Ejecutado con distribuciones adversas **realistas** (familias de AMPs, donde la
etiqueta correlaciona con el cluster):

| Caso | Resultado |
|---|---|
| 5 clusteres grandes con etiqueta correlacionada | **val con 0% de positivos, test con 100%** |
| 1 cluster dominante (60% de los datos) | **val y test con 0% de positivos** |
| 2 clusteres | **test vacio** |

**Efecto:** con val monoclase, `roc_auc_score` lanza `ValueError` y revienta a
mitad de la rejilla de hiperparametros, perdiendo todo el entrenamiento. Con test
monoclase es peor: ROC-AUC sale NaN y MCC sale 0.0 **sin aviso**, y se
publicarian numeros basura.

**Resuelto asi:** dos cambios. `assign_splits` siembra primero cada split con un
cluster de cada clase, eligiendo los mas pequenos para no gastar la cuota, y
solo despues reparte por deficit. Y `build_splits` valida antes de escribir:
cada split con ambas clases, al menos 20 positivos y al menos 10 clusteres (por
debajo de eso el bootstrap por cluster no significa nada).

**Limite que no se puede arreglar, y es correcto que no se pueda:** si solo hay
dos clusteres con positivos, no existe reparto valido para tres splits sin
partir un cluster, que seria fuga. En ese caso la validacion aborta con un
mensaje que explica la causa. Un test documenta el caso explicitamente.

---

## P-07 — El manejo de NaN del bootstrap solo cubria ROC-AUC
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Sintoma:** el docstring de `all_metrics` afirma que un remuestreo monoclase
devuelve NaN y el percentil lo ignora. Es falso para cinco de las seis metricas.

**Causa:** solo `roc_auc_score` lanza excepcion. Verificado con `y=[1,1,1]`:

```
{'roc_auc': nan, 'pr_auc': 1.0, 'mcc': 0.0, 'f1': 1.0,
 'sensitivity': 1.0, 'specificity': 0.0, 'precision': 1.0}
```

`average_precision_score` devuelve 1.0 o 0.0 duros, `matthews_corrcoef` devuelve
0.0 por denominador nulo, y los guardas de `threshold_metrics` devuelven 0.0.
Esos valores sobreviven al filtro `arr[~np.isnan(arr)]` y contaminan el intervalo
de PR-AUC y de MCC, que son las dos metricas principales del proyecto. En
`paired_difference` un remuestreo degenerado aporta un cero exacto, que contrae
el intervalo hacia "no distinguible".

**Resuelto asi:** `all_metrics` comprueba `len(y) == 0 or len(np.unique(y)) < 2`
antes de calcular nada y devuelve NaN para las siete metricas. El bootstrap
registra `n_boot_valid` y `n_boot_requested` por metrica y avisa si descarta mas
del 1%. Y `choose_threshold` ahora lanza una excepcion con validacion vacia o
monoclase, en vez de devolver el score mas bajo como umbral.

---

## P-08 — `boman_index` usaba la escala equivocada
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Causa:** la tabla `TRANSFER_FREE_ENERGY` es, residuo a residuo, la escala
**interfacial de Wimley-White** con el signo invertido (se identifica por W=1.85,
F=1.13, E=−2.02). El indice de Boman se define sobre las energias libres de
transferencia ciclohexano→agua de Radzicka y Wolfenden, donde **la arginina
domina** (R=14.92), seguida de D=8.72 y E=6.81. En la tabla del codigo la
arginina queda al nivel de N y Q, exactamente al reves de lo que hace util al
indice para AMPs cationicos ricos en arginina.

**Verificacion:** con la implementacion actual, LL-37 da 0.35 y la magainina-2
da 0.14. Con la escala publicada, LL-37 da ≈3.0 (APD3 reporta 2,9–3,0 kcal/mol)
y la magainina-2 ≈0.42. Un orden de magnitud de diferencia.

**Efecto:** el descriptor que se publica como "potencial de union a proteinas" es
en realidad una hidrofobicidad interfacial media con el signo cambiado, es decir
casi un duplicado de GRAVY.

**Resuelto asi:** ambas cosas. `RADZICKA_WOLFENDEN` sustituye a la tabla
anterior y `boman_index` ya no invierte el signo; verificado, LL-37 da
exactamente **3.00** y la magainina-2 **0.42**, que son los valores publicados.
La escala de Wimley-White se conserva con su nombre correcto como un descriptor
nuevo, `interfacial_hydrophobicity`, asi que los descriptores fisicoquimicos
pasan de 10 a 11. Tests de regresion contra los valores de la literatura y
contra el orden de la escala (la arginina debe dominar).

**Leccion:** el objetivo declarado de escribir las escalas explicitas en el
archivo era que fueran auditables contra la fuente. Funciono: una auditoria las
comparo residuo a residuo y encontro el error. Lo que fallo fue no haber hecho
esa comparacion al escribirlas, y haber verificado solo los descriptores
faciles (carga, momento) y no este.

---

## P-09 — XGBoost no carga: falta libomp
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Sintoma:** al ejecutar la etapa de baselines:

```
XGBoost Library (libxgboost.dylib) could not be loaded.
  Library not loaded: @rpath/libomp.dylib
  tried: '/opt/homebrew/opt/libomp/lib/libomp.dylib' (no such file)
```

**Causa:** la rueda de XGBoost para macOS depende de la biblioteca OpenMP de
Homebrew, que no se instala con `uv` porque no es un paquete de Python.

**Efecto:** `make baselines` muere tras entrenar los dos primeros modelos, y ese
trabajo se pierde. `XGBoostError` hereda de `ValueError`, asi que el CLI lo
captura y termina con codigo 1 y un mensaje poco informativo. Ningun test cubre
XGBoost, por eso los 89 pasaban.

**Resuelto asi:** `brew install libomp` instalado, y el import envuelto en una
excepcion propia, `ModelUnavailable`, que `run_baselines` captura para seguir
con los demas modelos en vez de abortar la etapa. Perder XGBoost debe costar un
baseline, no la corrida entera. Verificado: los 9 baselines entrenan.

---

## P-10 — Guardas que faltan, encontradas ejecutando la tuberia
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

Hallazgos de una corrida real sobre datos sinteticos. Ninguno es grave por si
solo; juntos son la diferencia entre una tuberia que falla con un mensaje util y
una que falla con una traza.

| # | Problema | Efecto |
|---|---|---|
| a | `interpret.py` cruza predicciones con `how="left"`: si hay predicciones rancias, `sequence` queda NaN y `physchem_matrix` lanza `TypeError`, que **no** esta en la tupla capturada por el CLI | Traza cruda al re-correr `make data` sin `make clean` |
| b | Ids duplicados: `dict(zip(df.id, df.sequence))` pisa una de las dos secuencias. La deduplicacion es por secuencia, no por id | Degradacion silenciosa; en el peor caso, error tardio y desconcertante en `evaluate` |
| c | Un `DataFrame` vacio tras los filtros llega a `np.vstack([])` | `ValueError: need at least one array to concatenate`, mensaje opaco |
| d | `choose_threshold` con val vacio devuelve `-inf`, que `json.dumps` escribe como `-Infinity` | JSON invalido para cualquier lector que no sea Python |
| e | `fractions` vacio en `params.yaml` deja la columna `split` en None y la verificacion de fuga pasa | `splits.csv` inservible; falla dos etapas mas tarde |
| f | `interpret.py` calcula los descriptores con los valores por defecto en vez de los de `params.yaml` | Si se cambia el pH o la ventana, el analisis de errores describe los peptidos con descriptores distintos de los que vio el modelo |
| g | Los bins de longitud estan cableados `[4,10,20,30,40,51]` | Si `max_length` sube de 50, los peptidos mas largos desaparecen del reporte en silencio |
| h | `ampesm all` omite la etapa `embeddings` que el `Makefile` si incluye | Cosmetico: una pausa larga e inesperada dentro de `train` |

**Resuelto asi:** las ocho. (a) se valida el cruce y `TypeError` se anade a las
excepciones que captura el CLI; ademas `make split` y `make clean` borran las
predicciones, porque una division nueva las invalida. (b) los ids duplicados se
desambiguan en `clean.py` con los 8 primeros caracteres del hash de la
secuencia, y `split.py` conserva una comprobacion de seguridad. (c) `clean.py`
aborta con un mensaje util si el conjunto queda vacio o de una sola clase.
(d) `choose_threshold` lanza excepcion en vez de devolver `-inf`. (e) se validan
las claves y la suma de `fractions`. (f) `interpret.py` lee `params.yaml`.
(g) los bins de longitud se derivan de `max_length`. (h) `ampesm all` incluye
`embeddings`.

---

## P-11 — Buena noticia: el respaldo de clustering no tarda horas
**Fecha:** 2026-10-11 · **Estado:** resuelto (era una preocupacion infundada)

Se temia que el clustering voraz en Python puro fuera inviable para ~16 000
peptidos. Medido con secuencias sinteticas de 5 a 50 residuos y 15% de homologos
cercanos:

| n | Tiempo | Clusteres | Alineamientos |
|---|---|---|---|
| 4 000 | 1,6 s | 3 279 | 734 |
| 8 000 | 6,4 s | 6 459 | 1 596 |
| 16 000 | **25,9 s** | 12 657 | 3 483 |

Escala cuadratica, pero con una constante pequenisima. Y el perfilado da una
sorpresa util: de 7,5 s a n=8 000, los alineamientos de Biopython son **0,127 s**
y el resto es el bucle de Python recorriendo listas invertidas y calculando
`len(kmers | rep_kmers[...])` para todos los candidatos. El cuello de botella no
es el alineamiento, es el prefiltro. Como los candidatos vienen ordenados por
k-meros compartidos y esa cuenta es una cota superior del Jaccard, se puede
cortar el bucle en cuanto la cota baje del umbral en vez de seguir recorriendo.
Importa poco a 26 s, pero pagara cuando se baje el umbral del prefiltro (P-03).

---

## P-12 — La auditoria de fuga daba un numero ininterpretable sin distribucion nula
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

Este salio de correr la auditoria recien escrita de P-04 sobre datos
sinteticos, y es el hallazgo mas instructivo de la sesion: **una verificacion
puede ser correcta y aun asi llevar a la conclusion contraria si no se compara
contra un control.**

**Sintoma:** la primera corrida de la auditoria dijo que el **29,3%** de las
secuencias de test tenian identidad >= 0.4 contra train. Leido tal cual, eso
significaria que la division por clusteres no sirve de nada y que todo el
proyecto esta roto.

**Causa:** la metrica de identidad tiene un ruido de fondo alto en peptidos
cortos. Con alineamiento global, huecos permitidos y normalizacion por la
secuencia mas corta, dos peptidos **sin ninguna relacion** alinean al 0.33 de
identidad de media, y uno de cada tres supera el 0.4 por puro azar. El umbral
del 40% esta **dentro del ruido** de la propia metrica.

Lo que delataba que algo no cuadraba: la identidad maxima observada era 0.54 y
no habia ni un solo caso por encima de 0.6. Las familias sinteticas se
construyeron con 0 a 2 mutaciones, asi que un pariente real habria dado ~0.95.
Si hubiera habido fuga de verdad, se veria en la cola alta del histograma, y la
cola alta estaba vacia.

**Resuelto asi:** la auditoria calcula ahora una **distribucion nula**. Baraja
los residuos de cada secuencia de test, lo que conserva longitud y composicion
exactas (que es lo que determina la identidad de fondo) y destruye cualquier
parentesco real, y repite la medicion. Se reportan tres numeros:

- `fraction_over_threshold`: el crudo, 29,3%.
- `null_model.fraction_over_threshold`: el ruido de fondo, 30,4%.
- `excess_over_null`: la diferencia, **-1,1%**. Este es el que mide fuga.

Y un `suggested_threshold_above_noise` (el percentil 99 del nulo, 0.62) por
encima del cual un parecido ya no se explica por azar: 0 secuencias de 92.

Conclusion correcta: **no hay fuga detectable**, justo lo contrario de lo que
sugeria el numero crudo.

**Leccion:** un numero sin su control no es evidencia, es una cifra. Es la misma
leccion que el proyecto predica para las metricas de clasificacion (de ahi los
baselines y los intervalos de confianza) y casi se incumple en la propia
herramienta de diagnostico. Si se hubiera publicado el 29% sin el nulo, habria
sido un error vergonzoso y perfectamente evitable.

---

## P-13 — Los 114 tests fallaron de golpe por el entorno, no por el codigo
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Sintoma:** a mitad de la sesion, los tests pasaron de pasar todos a fallar
todos en la coleccion con `ModuleNotFoundError: No module named 'ampesm'`.

**Causa:** la instalacion editable de `uv` deja un archivo
`_editable_impl_ampesm.pth` en `site-packages` que anade `src/` a `sys.path`.
En ese entorno los archivos `.pth` **dejaron de procesarse** por completo: ni
el del proyecto ni los de `virtualenv` o `distutils`. El archivo era correcto y
apuntaba a una carpeta existente; simplemente `site` no lo leia. En un entorno
creado desde cero con `python -m venv` el mismo mecanismo funcionaba, asi que
era ese `.venv` concreto y no la maquina.

**Resuelto asi:** dos cosas, y la segunda importa mas que la primera.

1. Se recreo el entorno, que lo arreglo.
2. Se anadio `pythonpath = ["src"]` a la configuracion de pytest, de modo que
   los tests importen el paquete directamente desde `src/` sin depender del
   instalable editable.

**Leccion:** un fallo del entorno que se presenta como un fallo del codigo cuesta
mucho mas de lo que deberia, porque uno empieza buscando en el sitio
equivocado. Que los tests no dependan del estado de la instalacion es barato y
evita exactamente esa confusion a quien clone el repositorio.

---

## P-14 — El corte temprano de la auditoria usaba una cota invertida
**Fecha:** 2026-10-11 · **Estado:** RESUELTO · **Introducido por el arreglo de P-03**

**Sintoma:** ninguno. La auditoria devolvia numeros plausibles que eran
aproximadamente la mitad de los correctos.

**Causa:** el corte temprano por k-meros usaba `n_shared / n_kmers_consulta` y
el comentario lo llamaba "cota superior" de la contencion. Es al reves. La
contencion es `n_shared / min(n_kmers_consulta, n_kmers_sujeto)`; si el sujeto
tiene **menos** k-meros, el denominador real es menor y la contencion real es
**mayor**, de modo que esa cantidad es una cota **inferior**.

Caso medido: consulta de 50 residuos (48 k-meros) y sujeto de 8 residuos que es
un fragmento exacto suyo (6 k-meros, los 6 compartidos).

| | valor |
|---|---|
| cota usada | 6/48 = **0.125** |
| contencion real | 6/6 = **1.000** |
| identidad real | **1.00** |

Con prefiltro 0.05 la cota no llegaba al umbral, se ejecutaba `break` y se
descartaba ese sujeto **y todos los siguientes**. Es exactamente el caso
corto-dentro-de-largo que P-03 habia arreglado, reintroducido dentro de la
herramienta que debia auditarlo.

En `cluster_builtin` el mismo codigo es inofensivo porque las secuencias se
procesan de mas larga a mas corta y el representante siempre tiene al menos
tantos k-meros que la consulta, con lo que la cota es exacta. En la auditoria
no hay ningun orden: las consultas son de test y los sujetos de train.

**Impacto medido antes del arreglo:** cambiaba la identidad maxima en el 41% de
las consultas, siempre a la baja, y la fraccion por encima del umbral pasaba de
0.80 a 0.37.

**Resuelto asi:** `break` por `continue` en `_max_identity_against`. En
`cluster_builtin` se conserva el corte, pero condicionado explicitamente a
`len(rep_kmers[rep_idx]) >= n_kmers` en vez de confiar en el orden, y con el
comentario corregido. Test de regresion `test_short_subject_is_not_skipped`.

**Leccion:** el comentario equivocado era lo que hacia que el codigo pareciera
correcto. Una cota mal razonada y bien explicada es peor que una sin comentar,
porque desactiva la revision. Y el patron de fondo se repite: **el arreglo de
un problema reintrodujo ese mismo problema en otro sitio**, porque se copio la
idea sin volver a comprobar sus condiciones.

---

## P-15 — El primer control de la auditoria no tenia potencia
**Fecha:** 2026-10-11 · **Estado:** RESUELTO · **Introducido por el arreglo de P-12**

**Causa:** P-12 establecio, con razon, que el numero crudo necesita un control.
El control elegido fue `excess_over_null`: la fraccion de test por encima del
umbral menos la misma fraccion en secuencias barajadas. El problema es que son
**dos proporciones saturadas**: con el ruido de fondo, ambas rondan 0.80, asi
que anadir parientes reales apenas puede moverlas. Hay techo.

Medido plantando 20 parientes en 120 secuencias de test (16,7% de fuga real,
que deberia dar ~+0.167):

| escenario | `excess_over_null` | ¿detecta? |
|---|---|---|
| sin fuga | −0.030 | — |
| 20/120 al 60% de identidad | +0.025 | **no** |
| 20/120 al 80% de identidad | +0.017 | **no** |

Y una sola baraja hacia el estadistico tan ruidoso (sd 0.037 sobre datos sin
fuga) que el aviso configurado en 0.02 se disparaba con ruido puro la mitad de
las veces.

**Resuelto asi:** ver P-18. El agregado se conserva en el reporte como
`excess_over_null_deprecated`, con su advertencia, para poder comparar con
corridas antiguas.

**Leccion:** tener un control no basta; el control tiene que **poder** detectar
el efecto que busca. Es la diferencia entre hacer el gesto de la rigurosidad y
ser riguroso. Un test que no puede fallar no es evidencia de nada.

---

## P-16 — `paired_difference` reventaba con remuestreos degenerados
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Causa:** el arreglo de P-07 se aplico a `bootstrap_metrics` pero no a
`paired_difference`, que es **la unica inferencia del proyecto**. Con
remuestreos degenerados, `arr` podia quedar vacio y `np.quantile` lanzaba
`IndexError: index -1 is out of bounds`. Peor que el fallo: no registraba
`n_boot_valid`, asi que el intervalo podia estar condicionado a una fraccion
arbitraria de remuestreos sin dejar rastro, que es justo lo que P-07 decia
haber resuelto.

**Resuelto asi:** misma guarda y misma contabilidad que en `bootstrap_metrics`,
mas error explicito si no sobrevive ningun remuestreo.

**Leccion:** al arreglar un problema, buscar todos los sitios donde vive el
mismo patron. Un `grep` del sintoma habria encontrado este en diez segundos.

---

## P-17 — Los empates en validacion se rompian por el resultado en test
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

**Causa:** `_select_comparison_pair` construia las listas de candidatos desde
`list(summary.tag)`, y `summary` viene ordenado por MCC **de test**. Como
`max()` conserva el primer maximo, un empate en ROC-AUC de validacion lo
desempataba el orden de test. Con un conjunto de validacion pequeno los empates
son plausibles.

**Resuelto asi:** `sorted(summary.tag)` y desempate por `(val_score, tag)`. El
desempate alfabetico es arbitrario, que es precisamente lo que se quiere:
arbitrario pero independiente de test. Ademas se imprime un aviso cuando no
estan las etiquetas preregistradas y se cae al criterio de validacion, porque
antes esa degradacion era silenciosa.

---

## P-18 — El segundo control tampoco detectaba nada, y por construccion
**Fecha:** 2026-10-11 · **Estado:** RESUELTO

Tercer intento de medir la fuga, y el unico que funciona. Se documenta entero
porque la secuencia de los tres errores es mas instructiva que el resultado.

**Causa:** el reemplazo de P-15 fue un p-valor empirico por secuencia con
correccion de Benjamini-Hochberg. Daba **cero fugas detectadas en todos los
escenarios, incluido uno con parientes al 80% de identidad**.

El motivo no era estadistico sino aritmetico: con R barajas, el p-valor minimo
posible es `1/(R+1)`. Para que sobreviva a BH sobre m pruebas hace falta
`1/(R+1) * m < 0.05`, es decir `R > 20m`. Con 120 secuencias de test harian
falta **mas de 2 400 barajas por secuencia**. Con R=20 el test daba cero
**por construccion**, no por ausencia de fuga.

**Resuelto asi:** en vez de exigir significacion secuencia a secuencia, se
marca una secuencia cuando su identidad observada supera **todas** sus propias
barajas, y se contrasta el **numero** de marcadas contra lo esperado por azar.
Bajo la nula cada secuencia se marca con probabilidad `1/(R+1)`, asi que el
conteo sigue una binomial y un test binomial da un p-valor global. Agregar el
conteo recupera la potencia que se perdia al exigir significacion individual.

Validado contra fuga plantada (20 de 120, R=20):

| escenario | marcadas | esperadas | p global | ¿detecta? |
|---|---|---|---|---|
| sin fuga | 6 | 5,7 | 0,51 | no (correcto) |
| 20/120 al 50% | 15 | 5,7 | 0,0006 | **si** |
| 20/120 al 60% | 27 | 5,7 | <1e-6 | **si** |
| 20/120 al 80% | 25 | 5,7 | <1e-6 | **si** |

Detecta desde el 50% de identidad, donde el estadistico de P-15 no veia nada ni
al 80%. Test de regresion `TestLeakageDetector`.

**Leccion, y es la de toda la sesion:** hicieron falta **tres** intentos para
medir una sola cosa, y los tres primeros fallaban en silencio devolviendo
numeros plausibles. Lo que rompio la racha no fue pensar mejor, fue **plantar
fuga conocida y exigir que el detector la encontrara**. Un diagnostico que
nunca se ha probado contra un caso positivo conocido no es un diagnostico.

---

## P-19 — La etapa `split` paso de segundos a horas
**Fecha:** 2026-10-11 · **Estado:** RESUELTO · **Introducido por el arreglo de P-18**

**Sintoma:** ninguno funcional. Todo seguia correcto, solo que inviable.

**Causa:** el nulo por secuencia de P-18 paso de **1** baraja a **20**, y el
coste de la auditoria se multiplico exactamente por 20. Peor: el coste crece de
forma cuadratica con el tamano del conjunto, porque cada secuencia de test se
compara contra todo train.

Medido por el agente de robustez:

| n peptidos | n_test | etapa `split` |
|---|---|---|
| 600 | 88 | 19,3 s |
| 2 000 | 298 | 152,6 s |
| 4 000 | 598 | **526,2 s** |

Atribucion a n=2 000: auditoria sin nulo 7,9 s, con nulo 160,1 s. El nulo era
practicamente toda la etapa. Extrapolado a los ~16 000 peptidos reales:
**alrededor de 2,3 horas solo en `split`**.

**Resuelto asi:** una seccion `split.audit` en `params.yaml` con
`null_replicates` (20) y `max_test` (500). Topar el conjunto auditado rompe el
crecimiento cuadratico sin perder potencia: el test binomial con 500 secuencias
detecta de sobra el exceso que importa. Medido tras el arreglo con 3 909
peptidos: **113 s frente a 526 s**, 4,6 veces mas rapido, y el coste ya no
crece con el tamano del test.

**Leccion:** un arreglo de correctitud puede romper la viabilidad, y eso
tambien es romperlo. El coste de un diagnostico debe mirarse al escribirlo, no
cuando alguien espera dos horas. Que el proyecto solo se hubiera probado con
600 peptidos sinteticos oculto el problema: **probar a una sola escala no es
probar.**

---

## P-20 — `evaluate` publicaba predicciones rancias en silencio
**Fecha:** 2026-10-11 · **Estado:** RESUELTO · **Arreglo incompleto de P-10a**

**Causa:** P-10a puso `rm -rf results/metrics/predictions` **solo en el target
`split` del Makefile**. `ampesm split` no borraba nada, y por tanto tampoco
`ampesm all`, que P-10h acababa de convertir en la tuberia completa. Ademas
`evaluate_all` no tenia la comprobacion de ids que si tiene `interpret.py`.

**Reproducido por el agente:** entrenar 9 modelos, reducir `peptides.csv` de
586 a 480 filas, correr `ampesm split` y `ampesm evaluate`. Resultado: **9 de 9
ficheros de metricas publicados**, cada uno con 106 ids inexistentes,
**codigo de salida 0 y ningun aviso**. En una variante con XGBoost caido, el
par **preregistrado** salio de un archivo rancio: la unica inferencia del
proyecto calculada sobre un artefacto viejo.

**Resuelto asi:** el borrado se mueve a `build_splits`, que es donde ocurre la
invalidacion, y se avisa de cuantos archivos se borran. Y `evaluate_all`
valida los ids contra `splits.csv` antes de calcular nada, igual que
`interpret.py`.

**Leccion:** poner una salvaguarda en el Makefile y no en el codigo deja fuera
a quien use el CLI directamente. La invariante pertenece al sitio donde se
rompe, no al atajo mas comodo para escribirla.

---

## P-21 — Los JSON de metricas podian no ser JSON valido
**Fecha:** 2026-10-11 · **Estado:** RESUELTO · **Agravado por el arreglo de P-07**

**Causa:** `json.dumps` escribe `NaN` e `Infinity` literales, que Python relee
sin quejarse pero que **no son JSON valido**: cualquier otro lector los
rechaza. P-10d arreglo el `-Infinity` de `choose_threshold` pero no el NaN, y
el arreglo de P-07 **agrando el agujero**: antes solo `roc_auc` podia ser NaN,
ahora pueden serlo las siete metricas.

**Resuelto asi:** `json_safe` convierte los valores no finitos en `null`, y
`write_json` escribe con `allow_nan=False` avisando cuando hubo alguno. Se
mapea a `null` en vez de fallar porque un NaN es informacion legitima (esa
metrica no estaba definida) y perder el reporte entero por eso seria peor.

---

## P-22 — Los bins de longitud solo estaban parametrizados por arriba
**Fecha:** 2026-10-11 · **Estado:** RESUELTO · **Arreglo incompleto de P-10g**

**Causa:** P-10g derivo el borde superior de `max_length` pero dejo el inferior
cableado en 4. Con `min_length: 3`, los peptidos de 3 residuos caian fuera de
todos los bins y `groupby(observed=True)` los descartaba **en silencio**:
exactamente el fallo que P-10g decia arreglar, en la otra direccion.
Verificado: 2 de 4 filas desaparecian. Ademas, con `max_length` grande el
ultimo bin era un unico cubo de 40 residuos y el reporte por longitud dejaba de
informar.

**Resuelto asi:** los bordes se derivan de `min_length` y `max_length`, con
pasos intermedios de 10.

**Leccion:** al parametrizar una constante, mirar si tiene gemela. Dos de los
arreglos de esta tanda (este y P-20) fallaron por quedarse a medio camino.

---

## P-23 — El umbral efectivo del prefiltro sigue por encima del declarado
**Fecha:** 2026-10-11 · **Estado:** ABIERTO, documentado

Medicion independiente sobre el codigo ya corregido:

| Identidad real | Agrupan | Los descarta el prefiltro | Los descarta el test de identidad |
|---|---|---|---|
| 0,43 | 68% | **32%** | 0% |
| 0,45 | 82% | 18% | 0% |
| 0,52 | 96% | 4% | 0% |
| 0,63 | 100% | 0% | 0% |

El prefiltro de k-meros sigue siendo la restriccion que manda: **ni un solo par
lo rechaza el test de identidad**. El umbral efectivo del respaldo esta en
0,55-0,60, no en 0,40. Es una mejora enorme frente al ~0,75 de P-03, que daba
entre 0% y 1% a 0,43, pero **no es el 0,40 declarado**.

Mitigado en la practica por `publishable: false` y por que MMseqs2 es
obligatorio para numeros publicables (D-11). Si alguna vez hiciera falta que el
respaldo fuera fiel, hay que bajar el prefiltro a 0,02-0,03 y volver a medir
el coste, que ya es el factor limitante (P-19).

**Correccion de una afirmacion previa:** el docstring y P-03 decian que bajar
el prefiltro de 0,15 a 0,05 costaba "casi el doble". Medido a 16 000
secuencias: 2,2 s a 22,6 s, es decir **unas 10 veces mas**, no el doble. El
numero absoluto sigue siendo aceptable; la caracterizacion relativa estaba mal.

---

## P-24 — Pendiente: cambiar `features` sin reentrenar no avisa
**Fecha:** 2026-10-11 · **Estado:** ABIERTO

Tras el arreglo de P-10f, `interpret.py` lee `params.yaml`, asi que los
descriptores del analisis de errores coinciden con la configuracion **actual**.
Pero si se cambia `features.moment_window` y se corre solo `interpret` sin
reentrenar, esos descriptores ya **no** son los que vio el modelo, y nada lo
detecta. Un hash de la seccion `features` guardado en `results/models/*.json` y
comparado en `interpret` cerraria el hueco.

## Plantilla para entradas nuevas

```
## P-NN — Titulo corto
**Fecha:** · **Estado:** resuelto | con rodeo | sin resolver

**Sintoma:** que se observo, con el mensaje de error literal si lo hubo.
**Causa:** la causa real, no la primera hipotesis.
**Resolucion:** que se cambio, con el archivo y la funcion.
**Leccion:** que se haria distinto la proxima vez.
```
