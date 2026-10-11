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
**Fecha:** 2026-10-11 · **Estado:** sin resolver, identificado antes de la primera corrida

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

**Leccion:** un test puede verificar correctamente un invariante y aun asi dejar
pasar el fallo que ese invariante pretendia prevenir, porque el invariante estaba
puesto en el eslabon equivocado de la cadena. Documentar un parametro como
"optimizacion de coste" sin medir su efecto en el resultado fue el error de
fondo: el prefiltro no era una optimizacion, era una segunda regla de agrupacion
escondida y mas estricta que la declarada.

---

## P-04 — Falta la verificacion directa de la premisa del proyecto
**Fecha:** 2026-10-11 · **Estado:** sin resolver

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

---

## P-05 — El modelo de la comparacion principal se elegia por su resultado en test
**Fecha:** 2026-10-11 · **Estado:** sin resolver

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

**Leccion:** escribir el compromiso previo no basta si el codigo que produce el
numero no lo implementa. El preregistro y la implementacion tienen que apuntar
al mismo sitio, y aqui no lo hacian.

---

## P-06 — `assign_splits` podia dejar un split vacio o de una sola clase, en silencio
**Fecha:** 2026-10-11 · **Estado:** sin resolver

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

**Correccion pendiente:** validar y abortar antes de escribir `splits.csv` (mismo
patron que ya se usa para la fuga: mejor no escribir el archivo que escribir uno
invalido), y convertir el balance de clases en restriccion y no en penalizacion.

---

## P-07 — El manejo de NaN del bootstrap solo cubria ROC-AUC
**Fecha:** 2026-10-11 · **Estado:** sin resolver

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

**Correccion pendiente:** detectar `len(np.unique(y)) < 2` **antes** de calcular
y devolver NaN para todas las metricas; registrar `n_boot_valid` por metrica y
avisar si se descarta mas del 1%, porque descartar remuestreos condiciona el
intervalo y hoy no queda rastro de cuantos se fueron.

---

## P-08 — `boman_index` usaba la escala equivocada
**Fecha:** 2026-10-11 · **Estado:** sin resolver

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

**Correccion pendiente:** sustituir por la escala de Radzicka y Wolfenden sin
invertir el signo, con un test de regresion de LL-37 ≈ 3.0. Si se quiere
conservar la escala actual como descriptor adicional, renombrarla a
`interfacial_hydrophobicity` y no llamarla Boman.

**Leccion:** el objetivo declarado de escribir las escalas explicitas en el
archivo era que fueran auditables contra la fuente. Funciono: una auditoria las
comparo residuo a residuo y encontro el error. Lo que fallo fue no haber hecho
esa comparacion al escribirlas, y haber verificado solo los descriptores
faciles (carga, momento) y no este.

---

## P-09 — XGBoost no carga: falta libomp
**Fecha:** 2026-10-11 · **Estado:** sin resolver, requiere instalacion

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

**Correccion pendiente:** `brew install libomp` (~2 MB), y envolver el import
para degradar con aviso en vez de abortar toda la etapa.

---

## P-10 — Guardas que faltan, encontradas ejecutando la tuberia
**Fecha:** 2026-10-11 · **Estado:** sin resolver

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

## Plantilla para entradas nuevas

```
## P-NN — Titulo corto
**Fecha:** · **Estado:** resuelto | con rodeo | sin resolver

**Sintoma:** que se observo, con el mensaje de error literal si lo hubo.
**Causa:** la causa real, no la primera hipotesis.
**Resolucion:** que se cambio, con el archivo y la funcion.
**Leccion:** que se haria distinto la proxima vez.
```
