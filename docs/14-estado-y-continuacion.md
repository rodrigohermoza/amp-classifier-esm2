# Estado del proyecto y como continuar

**Ultima actualizacion: 2026-10-17, tarde (tuberia completa corrida, resultados
publicados, reproducibilidad comprobada).** Este archivo es el punto de entrada
para retomar el trabajo. Lo que esta hecho, lo que falta y lo que hay que decidir.

---

## Resumen en una linea

**El proyecto tiene resultado.** La tuberia entera corrio sobre datos reales y la
comparacion preregistrada da un **empate**: ESM-2 35M + logistica frente a
composicion + fisicoquimicos con XGBoost, diferencia de MCC **+0.008 con IC 95%
[-0.033, +0.051]**, que incluye el cero. La division lleva fuga por homologia
**medida y declarada** (no eliminada), y esa declaracion viaja dentro de cada
JSON de metricas y en el pie de cada figura. Los dos README y `docs/13` ya tienen
los numeros reales; no queda ningun `[PENDIENTE]`.

Lo siguiente **no es urgente ni bloqueante**. Se puede retomar cuando se quiera.

---

## Lo que hay que hacer cuando se retome (resumen de 30 segundos)

Nada esta a medias. El orden recomendado, de mas a menos valor:

1. **M-01, la curva de rendimiento frente al umbral de identidad.** Es ahora la
   pieza mas valiosa del proyecto, no un extra: cuantifica cuanto importa la fuga
   declarada, que es la salvedad que acompaña a todos los numeros. Repetir con
   umbrales de 0.30, 0.40, 0.50, 0.70 y 0.90 y graficar MCC y ROC-AUC contra
   umbral. **Coste: 2-3 h**, y los embeddings NO se recalculan (dependen solo de
   las secuencias, no de la division). Ver la seccion propia mas abajo.
2. **La residualizacion por longitud** que cerraria del todo la asimetria de
   D-16 (`docs/10`). Ahora importa menos: la ablacion ya mostro que la longitud
   no sostenia el baseline (MCC 0.637 -> 0.626, dentro del ruido).
3. **M-02, evaluacion cruzada por fuente.** `interpret` ya da una pista de que
   hay señal de procedencia: exactitud 0.872 en APD3 y 0.697 en `uniprot_amp`.
4. **Fase 2, generacion.** Diseñada y sin implementar en
   [17-diseno-generacion.md](17-diseno-generacion.md). Es la continuidad con el
   proyecto de 2023 y la parte mas vistosa para el portafolio.
5. **M-07 (opcional, con etiqueta).** Fusionar las parejas cortas que MMseqs2 no
   ve. Solo como analisis de sensibilidad declarado, nunca como resultado: ver
   D-22 para por que.

Y una cosa que **no** hay que hacer: seguir ajustando el agrupamiento para que
la fuga baje. Esa decision esta tomada y razonada en D-22.

---

## El resultado, para no tener que buscarlo

**Comparacion preregistrada** (la unica inferencia del proyecto):

| Modelo | ROC-AUC | MCC |
|---|---|---|
| comp+fisicoquimicos, XGBoost | 0.909 [0.865-0.943] | 0.637 [0.540-0.716] |
| ESM-2 35M + logistica | 0.895 [0.838-0.940] | 0.645 [0.543-0.727] |

**Diferencia de MCC +0.008, IC 95% [-0.033, +0.051]: incluye el cero.**

Mejor modelo de los 16 por MCC: ESM-2 35M + MLP, ROC-AUC 0.913, MCC 0.659. Eso
es **descriptivo, no inferencia** (P-05).

Tres lecturas que no conviene perder:

- **El baseline no era debil**: 0.909 de ROC-AUC con 31 descriptores. La
  advertencia de `docs/01` se cumplio.
- **La longitud no era el truco**: la ablacion sin longitud da MCC 0.626 frente
  a 0.637, dentro del ruido. Eso convierte en numero una advertencia que estaba
  solo escrita, y hace la comparacion contra ESM-2 mas limpia de lo que se temia.
- **Los intervalos son anchos (±0.09 en MCC) a proposito**: el bootstrap
  remuestrea clusteres completos. Con secuencias saldrian estrechos y falsos.

De `interpret`: los falsos positivos tienen carga media +2.03 frente a +0.36 de
los verdaderos negativos. El modelo se equivoca donde hay carga alta sin
actividad, que es **exactamente** el fallo que `docs/17` predice para la fase 2
(la poliarginina). Vale la pena recordarlo al empezar esa fase.

---

## Estado por etapa

| Etapa | Codigo | Probado con sinteticos | Corrido con datos reales |
|---|---|---|---|
| `data` (descarga + limpieza) | si | si | **si, 2026-10-11** |
| `split` (clusteres + division) | si | si, con MMseqs2 | **si, 4 veces el 2026-10-17** |
| `baselines` (12 modelos) | si | si, XGBoost incluido | no |
| `embeddings` (ESM-2 8M y 35M) | si | si, pesos ya en cache | no |
| `train` (logreg y MLP) | si | si | no |
| `evaluate` (metricas + IC + figuras) | si | si | no |
| `interpret` (UMAP + errores) | si | si, 8 figuras | no |

**`data` ya corrio** y los cuatro controles salieron en orden (detalle en
[09](09-problemas-encontrados.md), seccion "Que dijo la primera corrida real"):
la cabecera espuria de APD3 fue exactamente 1, DRAMP encontro la columna
correcta, 271 conflictos de etiqueta y un `deficit` de emparejado de longitudes
de **23 732**, que es el numero que motivo D-16.

### Que dijo `split` en sus siete corridas (2026-10-17)

Dos de los cuatro controles salieron perfectos y se mantuvieron en todas las
corridas: **deriva de tasa positiva 0.000** (las tres divisiones en 0.568, el
guard de P-26 sin nada que decir) y **MMseqs2 usado** (18-8cc5c). El tercero,
`publishable`, mentia (P-29), luego bloqueaba para siempre, y ahora distingue
bloqueo de salvedad (D-22). El cuarto, la fuga, es la historia de la sesion.

> **Lee esta tabla con cuidado: las dos mitades NO son comparables, y la mitad
> de arriba tampoco lo es consigo misma.** Dos cosas la rompen:
>
> - en la corrida 5 cambio **que se mide** (P-32: la auditoria empezo a exigir la
>   misma cobertura que el agrupamiento), y al limpiar el ruido de fondo la
>   medicion gano potencia, asi que los numeros subieron;
> - hasta la corrida 7, dos corridas identicas daban **divisiones distintas**
>   (P-36: el 35% de las secuencias cambiaba de split), asi que parte de la
>   "tendencia" de las corridas 1 a 6 es ruido de reasignacion, no mejora.
>
> Se conserva entera porque el registro de lo que creimos tambien informa, pero
> **no se puede citar como una progresion**.

| Corrida | Cambio | Clusteres | Exceso test | Exceso val |
|---|---|---|---|---|
| 1 | tal cual estaba | 8 173 | +50.2 (p=0) | +47.2 (p=0) |
| 2 | `--cluster-mode 1` (P-29) | 7 886 | +34.2 (p=0) | +32.2 (p=0) |
| 3 | fusion por contencion (P-30/D-20) | 7 438 | +20.2 (p=8.7e-05) | +30.2 (p=0) |
| 4 | arreglado el `break` de P-31 | 7 178 | +10.2 (p=0.025) | +28.2 (p=0) |
| | *— aqui cambia lo que se mide (P-32) —* | | | |
| 5 | la auditoria exige cobertura (P-32/D-21) | 7 153 | +68.2 (p=0) | +62.2 (p=0) |
| 6 | sin cambios de agrupamiento, solo D-22 | 7 153 | +75.2 (p=0) | +66.2 (p=0) |
| | *— aqui se arregla la reproducibilidad (P-36) —* | | | |
| 7 | ids de cluster canonicos (P-36) | 7 153 | **+54.2 (p=0)** | **+65.2 (p=0)** |

El esperado por azar es 23.81 marcadas de 500 en todas las corridas. **La
corrida 7 es la division canonica en uso**, y es la primera que otra corrida
identica reproduciria: las corridas 5 y 6 son el mismo codigo de agrupamiento y
dan numeros distintos, que es precisamente como se encontro P-36.

### Por que la fuga subio al arreglar la medicion, y por que eso es correcto

Es lo mas importante que entender de esta sesion. Durante cuatro corridas la
auditoria no exigia cobertura, asi que un alineador global con huecos encontraba
cinco residuos dispersos dentro de casi cualquier peptido y devolvia identidad
1.000. Ese ruido **tambien contaminaba el nulo**: su p99 era 1.000, es decir las
barajas tambien saturaban, y la regla de marcado ("la observada supera a sus 20
barajas") no podia dispararse casi nunca. El `+10.2` de la corrida 4 no era poca
fuga: era un test **sin potencia**.

Al exigir cobertura, el ruido de fondo bajo de 0.748 a 0.540 y el p99 del nulo de
1.000 a 0.875, y la fuga real se hizo visible. Se comprobo mirando **que** se
marca, que es la unica forma de saberlo:

| banda de identidad con su pariente en train | test | val |
|---|---|---|
| inequivoca, >= 0.80 | 44% | 56% |
| alta, 0.60-0.80 | 44% | 30% |
| dudosa, 0.40-0.60 | 12% | 14% |

`GILDAIKAIAKAAG` vs `GILDPIKAFAKAAG`, `ILGKIWKGIKSLF` vs `ILGKIWKIKKLF`,
`VQLRIRVRVIRK` vs `VRLRIRVRVIRK`. Eso es homologia, no aritmetica. Entre 19 y 24
de cada 50 marcadas miden **menos de 15 residuos**, que es donde el prefiltro de
k-meros de MMseqs2 falla por construccion (P-30).

Conclusion: la fuga no empeoro, estaba tapada. Y el proyecto no la esconde: se
declara con cada metrica (D-22) y se cuantifica con M-01.

## Numeros actuales

| Metrica | Valor |
|---|---|
| Commits | 36 |
| Codigo en `src/` | 4 892 lineas, 16 modulos |
| Tests | **162 casos, todos pasan** (5,9 s) |
| Documentacion | 5 831 lineas, 18 archivos |
| Datos procesados | 18 644 peptidos (10 594 + / 8 050 -) |
| Clusteres (corrida 7, canonica) | 7 153 (5 615 de 1), 88 con ambas clases |
| Fuga declarada | test +54.2, val +65.2, ambas p = 0, sobre 23.81 por azar |
| Lint | ruff limpio |
| CI | en verde (se cayo dos veces hoy: P-33 y P-34) |
| Repositorio clonable | ~1 MB |
| Entorno virtual | 1,0 GB, **en `~/.venvs/ampesm`, fuera de iCloud** (D-19) |
| Cache de pesos ESM-2 | 160 MB (`.cache/`, ignorada por git) |
| Problemas registrados | 38 (mas P-10a-h y P-25b/P-31b como subapartados); **1 abierto: P-23, a proposito** |
| Decisiones registradas | 22 (D-01 a D-22) |
| Repositorio publico | <https://github.com/rodrigohermoza/amp-classifier-esm2> |

### Coste medido de cada etapa en el M2 (8 nucleos, 24 GB)

Medido, no estimado. Util para decidir que cabe en una sesion:

| Etapa | Tiempo | Nucleos |
|---|---|---|
| `data` | 83 s | 1 |
| `split` | **12,5 min** | 1 (la auditoria no paraleliza) |
| `baselines` | ~3-5 min | 8 |
| `embeddings` | **7 min** (1,2 el de 8M + 5,3 el de 35M) | 8 |
| `train` | ~15-25 min | 8 |
| `evaluate` | ~5-10 min | 1-2 |
| `interpret` | ~5 min | 2-4 |

Pico de memoria muy por debajo de 4 GB: el conjunto son 488 802 residuos y la
matriz de embeddings mas grande son 36 MB. Nada de esto necesita GPU, y el motivo
es que los peptidos son cortos (26 residuos de media). Con proteinas completas
seria otra conversacion.

Palancas si `split` molesta: `split.audit.audit_val: false` quita ~6 min y
`null_replicates` de 20 a 10 otros ~3. Las dos degradan justo las guardas que
distinguen el proyecto, asi que a 12 minutos no vale la pena.

## Dependencias del sistema ya instaladas

```bash
brew install mmseqs2   # 18-8cc5c. OBLIGATORIO para numeros publicables (D-11)
brew install libomp    # sin ella XGBoost no carga en macOS (P-09)
```

---

## Lo que falta, en orden

### 1. `make data` — HECHO (2026-10-11)
18 644 peptidos en 83 s, 5.01 MB de red. Los cuatro controles y su lectura, en
[09](09-problemas-encontrados.md). Lo unico que salio mal fue el emparejado de
longitudes, y la decision al respecto es D-16.

Si hace falta rehacerlo: `make clean-data && make data`. Y si `ampesm` no
importa, es P-25: `uv sync --reinstall-package ampesm && uv sync --extra dev`.

### 1b. Los cuatro arreglos previos a la tuberia — HECHOS (2026-10-11)
Los cuatro estaban pensados para proteger numeros que **todavia no se han
publicado**, asi que se hicieron antes de correr nada:

- **P-27**: la auditoria de fuga mira ahora tambien **val contra train**, no
  solo test. val fija el umbral que se aplica a test, asi que la homologia
  val-train sesga cinco de las siete metricas de test aunque test este limpio.
  El reporte gana `train_val_identity_audit`.
- **P-26**: se vigila la **deriva de la tasa de positivos** entre splits
  (`positive_rate_drift` en el reporte, tolerancia 0.05). Avisa, no aborta.
- **P-24**: hash de la seccion `features` guardado al entrenar y comprobado en
  `interpret`, que **aborta** si cambio (criterio en D-18).
- **P-25**: `make check-env` comprueba el import antes de cada etapa. Se valido
  sola en vivo: el entorno se rompio por cuarta vez mientras se escribia.

### 1c. Lo arreglado el 2026-10-17

- **P-28**: iCloud sincronizaba `~/Desktop` y renombraba archivos dentro del
  entorno virtual (`__init__ 2.py`), dejandolo intacto en disco e invisible para
  Python. 4 785 archivos afectados. **Era la causa real de P-13 y P-25**, que
  reaparecieron cinco veces porque se trataba el sintoma. Segundo efecto, el
  grave: `uv run pytest` caia al pytest de anaconda (Python 3.10) **sin dar
  error**, asi que `make test` medio otro interprete durante un tiempo
  indeterminado. El entorno vive ahora en `~/.venvs/ampesm` con `.venv` como
  enlace (D-19), y `check-env` distingue las dos causas.
- **P-29**: el agrupamiento no pasaba `--cluster-mode` y el defecto resulto ser
  el 2 (incremental voraz), cuya garantia es "cada miembro se parece a **su**
  representante" y no dice nada entre miembros de clusteres distintos. **La
  garantia que el proyecto declara no se deducia del metodo usado.** Puesto
  `--cluster-mode 1`. En el mismo hallazgo: `publishable` decia `true` con fuga
  significativa en los dos lados, porque solo miraba la herramienta; ahora se
  calcula tras las auditorias y enumera los motivos.
- **P-30 / D-20**: el 8.5% del conjunto (1 582 secuencias) son subcadena exacta
  de otra. MMseqs2 no las ve a 8-15 residuos y no hay configuracion que lo
  arregle. Se fusionan por contencion exacta, **no** alineando: igualar la regla
  de la guarda con la del agrupamiento la dejaria ciega a sus propios puntos
  ciegos.
- **P-31 y P-31b**: dos bugs mios en esa fusion, encontrados midiendo. El bucle
  cortaba en el primer contenedor (y 803 de 1 582 estan en varios), y el
  reetiquetado permitia que dos fusiones solapadas se sobrescribieran.

- **P-32 / D-21**: la auditoria descartaba la cobertura que el agrupamiento
  exige, asi que marcaba parejas que no son homologas (`CCISG`, cinco residuos,
  identidad 1.000 por aritmetica). Al arreglarlo, **el arreglo que estaba
  escrito no servia**: la cobertura que ya se calculaba se normaliza por la
  secuencia corta y vale 1.000 justo en los casos malos. La que discrimina es la
  normalizada por la larga, que es la de MMseqs2. La condicion final es una **O**
  (cobertura sobre la larga **o** contencion exacta), con un termino por cada
  garantia que el proyecto afirma.
- **P-33**: el enlace `.venv` se colo en un commit, porque `.gitignore` decia
  `.venv/` con barra y la barra solo casa directorios. **El repo no se podia
  clonar y usar**; el CI se caia en 10 s.
- **P-34**: un test que solo pasaba con MMseqs2 instalado. El CI lo encontro.
- **P-35**: el codigo documentaba **al reves** lo que hace `--cov-mode 1`. La
  frase llevaba semanas ahi y se leyo como premisa al arreglar P-32.
- **P-36, el mas grave**: dos corridas identicas de `split` daban divisiones
  distintas. MMseqs2 escribe las filas del TSV en orden variable (es multihilo)
  y los ids de cluster salian de ese orden, asi que `assign_splits`, que baraja
  con semilla fija partiendo del orden de los ids, asignaba distinto. **El 35% de
  las secuencias cambiaba de split entre dos corridas identicas.** Ningun numero
  anterior era reproducible. Arreglado con ids canonicos y comprobado de punta a
  punta.
- **D-22**: `publishable` separa ahora **bloqueos** (nada publicable, `evaluate`
  aborta) de **salvedades** (se publica declarandolas, y la declaracion va dentro
  de cada JSON y en el pie de cada figura). Con la puerta anterior, una fuga real
  y declarada dejaba el proyecto sin producir **ningun** numero, para siempre.

Queda **un solo problema abierto**:

- **P-23**, deliberadamente abierto: no afecta a ningun numero publicable porque
  MMseqs2 es obligatorio y el respaldo bloquea la publicacion.

### 2. `make split` — HECHO, con la fuga declarada (2026-10-17, corrida 7)

Historia completa en la seccion "Que dijo `split`" de arriba, y el detalle en
P-29, P-30, P-31, P-32, P-35 y P-36.

**Estado: `publishable: true` CON SALVEDADES.** No significa que la division este
limpia: significa que su defecto esta **medido, acotado y declarado**, y que la
declaracion viaja pegada a cada numero (D-22). Fuga en test +54.2 (p = 0) y en
val +65.2 (p = 0), sobre 23.81 marcadas esperadas por azar.

La decision de no seguir ajustando esta tomada y razonada en D-22: la via para
cerrar el residuo es **M-01**, no mas parametros. La alternativa de fusionar las
parejas cortas (M-07) queda como analisis de sensibilidad etiquetado, porque
usaria la misma regla que la auditoria y daria verde por construccion.

Que mirar en `data/processed/split_report.json`:

- **`publication_blockers` debe estar vacio.** Si tiene algo, no hay nada
  publicable y `evaluate` **aborta** (D-18): o se uso el respaldo propio en
  Python, o una auditoria falto o fallo.
- **`publication_caveats`**: lo que hay que declarar al publicar. Hoy son dos,
  las dos de fuga. No es una lista de pendientes: es el texto que `evaluate`
  copia dentro de cada JSON de metricas y en el pie de cada figura.
- `train_test_identity_audit.leakage.global_p_value` y `excess`. **Ese** es el
  numero de fuga, y se lee junto a `min_coverage`, que dice con que definicion
  de homologo se midio (sin eso el p-valor no es interpretable, P-32). NO leer
  `fraction_over_threshold` (dara ~98% y es ruido de fondo de la metrica) ni
  `excess_over_null_deprecated`. Historia en P-12, P-15 y P-18.
- `train_val_identity_audit`: lo mismo para val. Sale peor que test y es
  esperable: a val le toco mas peptido corto.
- `contained_merge`: 1 582 secuencias eran subcadena de otra, 733 fusiones de
  clusteres (P-30, D-20).
- `positive_rate_drift.exceeds_tolerance` debe ser `false`. En las siete
  corridas salio 0.000, perfecto.
- `mixed_label_clusters`: clusteres con ambas clases. Corrida 7: 88.

**Reproducibilidad: resuelta, y conviene no darla por supuesta otra vez.** Hasta
P-36, dos corridas identicas daban divisiones distintas (el 35% de las
secuencias cambiaba de split) porque los ids de cluster salian del orden de las
filas del TSV de MMseqs2, que es multihilo. Si alguna vez un numero de `split`
se mueve sin que nada haya cambiado, **ese** es el sospechoso; y la forma de
comprobarlo sin esperar 12 minutos esta en el script de validacion que se uso:
pasar dos TSV con ordenes distintos por parseo + fusion + asignacion y comparar
las divisiones.

### 3. `make baselines`, `make embeddings`, `make train` — HECHO (2026-10-17)

Corrido entero. Costes reales medidos en esta corrida, muy por debajo de lo
estimado porque los pesos y los embeddings ya estaban en cache:

| Etapa | Tiempo real |
|---|---|
| `baselines` (4 representaciones x 3 modelos) | **39 s** |
| `embeddings` (los dos modelos, ya en cache) | segundos |
| `train` (2 modelos ESM-2 x 2 clasificadores) | **30 s** |
| `evaluate` (16 modelos, bootstrap 2000 por cluster) | **190 s** |
| `interpret` (UMAP + analisis de errores) | **60 s** |

`split` sigue siendo la etapa cara: 12 min en condiciones normales, y se midio
una corrida de 20 min solo por competencia de CPU con el dibujado de pantalla y
la app de Claude. **Mas lento no significa distinto** (P-36).

### 4. `make evaluate` y `make interpret` — HECHO (2026-10-17)

El resultado esta arriba, en "El resultado, para no tener que buscarlo". Lo que
hay en `results/`:

- `metrics/comparison_embeddings_vs_baseline.json` — **la unica inferencia**.
- `metrics/summary.csv` — los 16 modelos, descriptivo.
- `metrics/summary_caveats.json` — las salvedades, porque un CSV no tiene donde
  llevarlas (D-22).
- `metrics/<modelo>.json` — uno por modelo, cada uno con `split_caveats` dentro.
- `figures/` — todas con las salvedades en el pie.

**Reproducibilidad comprobada de punta a punta.** Se guardaron los artefactos, se
relanzo `split` + `baselines` + `train` + `evaluate` y se comparo: `splits.csv` y
la comparacion preregistrada **identicos byte a byte**; en `summary.csv` una sola
diferencia de 1.11e-16 en un umbral, que es el ultimo bit de un flotante. Todas
las metricas coinciden a 12 decimales.

### 5. Rellenar la documentacion — HECHO (2026-10-17)

Los dos README y `13-portafolio.md` tienen los numeros reales. No queda ningun
`[PENDIENTE]`.

`01-objetivos.md` **no se toco**, y es importante que siga como esta: la
expectativa escrita ("los embeddings deberian ganar, pero por margen modesto,
quiza con intervalos solapados") se cumplio, y aun mas conservadora. Si alguna
vez hay que corregirla, se corrige **declarando que se corrige**, nunca
reescribiendola para que encaje.

### 5b. M-01, la curva de umbral de identidad — LO SIGUIENTE

Es lo mas valioso que queda, y ya no es un extra: cada metrica publicada va
acompañada de una salvedad que **remite a esta curva**.

**Que es:** repetir la tuberia con `identity_threshold` en 0.30, 0.40, 0.50,
0.70 y 0.90, y graficar MCC y ROC-AUC contra el umbral. La pendiente **es** la
medida directa de cuanta fuga tenia cada division.

**Coste, con los tiempos ya medidos:** por cada umbral, `split` ~12 min +
`baselines` 39 s + `train` 30 s + `evaluate` 190 s = **~16 min**. Cinco umbrales
= **80 min de maquina**, mas el codigo para orquestarlo y dibujar la curva.
Presupuesto realista con margen: **2-3 h**.

**Lo que NO hay que recalcular:** los embeddings. Dependen solo de las
secuencias, no de la division, y ya estan en cache. Caer en recalcularlos
multiplicaria el coste sin cambiar un numero.

**Que mirar al terminar:** si la curva baja al subir el umbral, la fuga era real
y esta cuantificada. Si es plana, el modelo no dependia de la homologia tanto
como se temia, y eso tambien es un resultado publicable. Las dos salidas sirven;
ninguna obliga a retocar nada.

### 6. Fase 2: generacion de peptidos
Diseñada y aprobada, sin implementar. Todo en
[17-diseno-generacion.md](17-diseno-generacion.md). Lo esencial: el encuadre es
**medir cuanto engaña el generador al clasificador**, no diseñar farmacos
(D-17), hacen falta dos oraculos de representaciones distintas, y no empieza
hasta que la clasificacion este corrida, porque un peptido con puntuacion 0.95
de un clasificador cuyo rendimiento real se desconoce es un numero sin unidades.

### 7. Verificacion: las dos auditorias estan cerradas

**Llegaron los dos informes.** Entre ambos encontraron **11 fallos, todos
introducidos por los arreglos anteriores y no por el codigo original**
(P-14 a P-24). Nueve estan corregidos; dos quedan abiertos a proposito:

- **P-23**: el umbral efectivo del prefiltro del respaldo esta en 0,55-0,60, no
  en el 0,40 declarado. Medido: a 0,43 de identidad real se agrupa el 68% de
  los pares y el 32% restante lo descarta el prefiltro, no el test de
  identidad. Es una mejora enorme frente al ~0,75 anterior, pero no es fiel.
  Mitigado por `publishable: false` y por MMseqs2 obligatorio (D-11). Bajar el
  prefiltro a 0,02-0,03 lo arreglaria, pero hay que medir el coste, que ya es
  el factor limitante.
- ~~**P-24**~~: cerrado el 2026-10-11 con el hash de `features`.

Lo que las auditorias **confirmaron como correcto**, ejecutandolo: la escala de
Boman completa (20 valores, 0 discrepancias), el manejo de NaN en
`all_metrics` y `bootstrap_metrics`, que las etiquetas preregistradas coinciden
con los nombres reales de archivo, la seleccion por validacion en `evaluate` e
`interpret`, que la siembra no rompe el invariante y da 70/15/15 exactos, que
los 5 casos de division degenerada abortan sin escribir archivos, las ocho
guardas de P-10 con mensajes claros y sin trazas, que `ampesm all` y `make all`
corren lo mismo, dos tuberias completas de principio a fin con 16 JSON validos
y 8 figuras, y que los tests pasan aunque se retire el `.pth` del editable.

Los dos puntos de rigor del informe metodologico que quedaban (auditar val, y
validar la deriva de la tasa positiva) **estan hechos**: P-27 y P-26.

---

## Decisiones abiertas para Rodrigo

**Cerradas el 2026-10-11:**
- *Publicar ahora o esperar*: se publico ahora, con las tablas vacias y el
  aviso de "pendiente de ejecucion" en los dos README.
- *Nombre publico*: `amp-classifier-esm2`. La carpeta local sigue llamandose
  `Peptidos-antimicrobianos`, que no afecta a nada.
- *El confusor de longitud* (punto 4 de la lista anterior): ablacion, no
  emparejado por split. Razonamiento completo en D-16.

**Abiertas:**

1. **Lo de GitHub, pendiente de Rodrigo y no de nadie mas.** La descripcion del
   repositorio todavia describe **otro proyecto**: habla de generar peptidos
   para agricultura. Con el README ya con resultados, desentona mas. Se cambia
   desde la propia pagina del repo, en el engranaje de "About".
2. **El nombre del repositorio cuando llegue la fase 2.** `amp-classifier-esm2`
   describe bien la mitad actual y se quedara corto cuando haya generacion.
   Cambiar el nombre en GitHub mantiene las redirecciones, asi que no rompe
   enlaces, pero conviene decidirlo antes de que haya mucho enlazado.
3. **¿Se hace la regresion de MIC con los MIC de DRAMP?** El reexamen de P-02
   cambia el planteamiento: DBAASP no cambio de dominio, su portada responde y
   lo que esta roto es todo lo que devuelve datos (200 con 0 bytes, incluida la
   documentacion de su API). Pero **DRAMP trae MIC en 5 020 de sus 12 784
   filas**, en `Target_Organism` como texto libre, y ese archivo ya se descarga
   y ya se parsea. Lo que cuesta es el parser: regex sobre `MIC=`/`MIC≤`,
   unificar `μg/ml` frente a `µM` (necesita peso molecular) y decidir como
   colapsar varios organismos por peptido. Sigue fuera de alcance (D-07) hasta
   que se decida.

---

## Mejoras priorizadas que siguen sin hacer

Detalle en [10-mejoras-futuras.md](10-mejoras-futuras.md). Las tres de mejor
relacion esfuerzo/beneficio:

1. **Curva de rendimiento frente al umbral de identidad** (30/40/50/70/90%).
   Probablemente la figura mas instructiva del proyecto: la pendiente *es* la
   medida de cuanta fuga tenia la division aleatoria. La tuberia ya esta
   parametrizada.
2. **Negativos barajados** que conserven la composicion de los positivos. Es la
   prueba mas afilada: si ESM-2 no supera a la composicion con esos negativos,
   no esta aportando nada sobre el orden de los residuos. Ahora hay ademas
   codigo que ya baraja secuencias (la nula de la auditoria) y se puede
   reutilizar.
3. **Evaluacion cruzada por fuente** (entrenar en APD3, evaluar en UniProt y al
   reves). Convierte en numero la advertencia sobre procedencia.

Y una de higiene: sacar `jupyter` del extra `dev` a uno propio, para que el CI
no instale ~60 MB de babel y jedi.

---

## Riesgos conocidos que siguen vivos

1. **Los negativos son no-anotados, no verificados.** Es la debilidad central y
   esta escrita en los dos README. No se puede arreglar con codigo.
2. **El respaldo de clustering en Python no sirve para numeros publicables**
   (D-11). Da clusteres distintos de MMseqs2. Mitigado con `publishable: false`
   en el reporte y un aviso ruidoso, pero conviene no olvidarlo.
3. **El respaldo de clustering agrupa la mitad de los pares a 40-45% de
   identidad.** Medido de extremo a extremo. Es la razon de D-11; MMseqs2 es
   obligatorio para publicar. Tabla de recall en docs/04.
4. **La longitud separa las clases y no se puede corregir emparejando.**
   Confirmado con datos reales: positivos 22.9 residuos de media, negativos
   30.6, y un `deficit` de 23 732 negativos cortos que no existen en UniProt.
   Mitigado, no resuelto, con la ablacion `comp_physchem_nolen` (D-16). El KS
   por split sigue reportandose como diagnostico.
5. **La ablacion de D-16 es asimetrica**: quita la longitud al baseline y no a
   ESM-2, asi que esa lectura esta sesgada a favor de los embeddings. Por eso
   la comparacion preregistrada sigue siendo contra `comp_physchem` completo.
6. ~~Las etapas 3 a 7 nunca se han ejecutado con datos reales.~~ **Cerrado el
   2026-10-17**: la tuberia completa corrio, dos veces, con el mismo resultado.
7. **La division tiene fuga por homologia real, medida y declarada.** Test
   +54.2 y val +65.2 sobre 23.81 esperadas por azar, ambas p = 0. No bloquea
   publicar (D-22), **pero ninguna metrica de este proyecto se cita sin ella**.
   La via para cerrarla es M-01, no mas ajustes. Entre el 74% y el 86% de lo
   marcado tiene identidad >= 0.60 con un pariente en train: es homologia, no
   artefacto de la metrica.
8. **La fuga se concentra en peptidos cortos y ahi MMseqs2 no puede verla.** Su
   prefiltro necesita 6-7 k-meros coincidentes y un peptido de 8 residuos tiene
   dos o tres. No es configuracion, es estructural. La fusion por contencion
   (D-20) tapa las subcadenas exactas; las sustituciones de un residuo
   (`ILGKIWKGIKSLF` vs `ILGKIWKIKKLF`) siguen cruzando. Cerrarlas requiere
   M-07, que no puede ser resultado principal (D-22).
8b. **Las cifras de fuga anteriores al 2026-10-17 no son comparables con las de
   despues**, por dos motivos acumulados: cambio lo que se mide (P-32) y hasta
   P-36 dos corridas identicas daban divisiones distintas. La tabla de las siete
   corridas lo marca; **no se puede citar como una progresion**.
9. **El proyecto vive en una carpeta sincronizada por iCloud.** El entorno ya
   esta fuera (D-19) y hoy no hay dano en `src/`, `docs/` ni `.git` —medido, 0
   archivos afectados— pero no hay garantia de que iCloud no toque algo mas
   algun dia. El repositorio remoto es la copia de seguridad real.

---

## Lo que mas conviene recordar

**De la sesion del 11 de octubre.** Tres intentos hicieron falta para construir
un detector de fuga que detectara fuga. Los tres primeros devolvian numeros
plausibles y no encontraban nada, ni con parientes al 80% de identidad plantados
a mano. Ninguno fallaba de forma visible. Y dos de los cinco fallos nuevos los
**introdujeron los arreglos anteriores**. De ahi D-15 y el habito que vale la
pena conservar: **plantar el efecto que buscas y exigir que tu herramienta lo
encuentre**, antes de creerte el numero que te devuelve.

**De la sesion del 17 de octubre**, cuatro cosas que costaron caro:

1. **Un fallo que reaparece cinco veces sin causa identificada no esta
   documentado, esta anotado.** "El `.pth` se desprocesa de vez en cuando" era
   una descripcion con tono de explicacion. La pregunta que lo habria cortado en
   la primera reaparicion es la mas barata que existe: *¿que mas hay en esta
   carpeta que pueda estar escribiendo en ella aparte de mi?*
2. **El valor por defecto de una herramienta se comporta como una decision del
   proyecto aunque no lo sea.** Aqui el defecto de `--cluster-mode` elegia la
   semantica del agrupamiento, que es **la** decision metodologica central del
   repositorio, y no estaba escrita en ninguna parte. Todo parametro del que
   dependa una afirmacion del README se pasa explicito aunque coincida con el
   defecto.
3. **Un caso plantado que no reproduce el mecanismo da un verde que no vale.**
   Paso dos veces el mismo dia: la guarda de P-28 no salto porque `uv run`
   repara el `.pth` por su cuenta, y el primer test de P-31b lo pasaban las dos
   versiones, la buena y la mala. La unica forma de saberlo es **ejecutar la
   version mala y verla fallar**. D-15 sin ese paso es media regla.
4. **Conocer un patron no basta para no repetirlo.** P-31 es P-14 otra vez —un
   `break` que pierde parejas— cometido veinte lineas por debajo del comentario
   de ocho lineas que lo advierte. Hay que buscarlo activamente en cada bucle
   nuevo que recorra candidatos, no confiar en recordarlo.

Y una distincion que conviene no perder, porque es la parte fina de todo esto:
**compartir la definicion no es compartir la regla.** La guarda y lo vigilado
tienen que estar de acuerdo en *que cuenta como homologo* (identidad y
cobertura), o la guarda mide otra cosa que lo garantizado — eso era P-32. Pero no
deben usar *el mismo alineador*, porque entonces la guarda hereda sus puntos
ciegos y no puede contradecirlo — eso es D-20. La primera igualdad es
obligatoria; la segunda, inadmisible.

**De la tarde del 17 de octubre**, cinco mas, todas pagadas:

5. **Lo que no se mide, no se arregla aunque se escriba el arreglo.** El arreglo
   de P-32 llevaba un dia escrito en los docs y **no servia**: al medirlo antes
   de implementarlo, la cobertura que proponia exigir valia 1.000 justo en los
   casos que lo motivaban. Un arreglo razonado sigue siendo una hipotesis.
6. **La documentacion que miente es peor que la que falta**, porque se usa para
   decidir. El docstring de `--cov-mode 1` decia lo contrario de lo que hace
   (P-35), y se leyo como premisa. Un comentario que explica por que un
   parametro vale lo que vale es una afirmacion comprobable, y se comprueba.
7. **Un test que estorba puede estar defendiendo algo.** El filtro de cobertura
   rompio la regresion de P-14. La pregunta correcta no era si borrarlo, sino
   **que garantia defendia**, y la respuesta cambio el diseño: la condicion paso
   a ser una O con un termino por mecanismo (D-21).
8. **El orden es la tercera forma del mismo fallo.** P-14 (una cota mal puesta),
   P-31 (iteracion de un `set`) y P-36 (el orden de las filas de un TSV). La
   pregunta que los detecta: *¿alguien lee esto en orden, y ese orden decide
   algo?* En P-36 el lector era un `rng.shuffle` tres funciones mas abajo.
9. **El CI es la unica maquina que prueba el repo como lo recibe un tercero.**
   Dos fallos del dia (P-33, P-34) pasaban los 145 tests en local y hacian que el
   repositorio **no se pudiera clonar y usar**. Para un repo cuyo proposito es
   que lo lean otros, eso no es un detalle de infraestructura.

Y la que resume la sesion entera: **relanzar algo que no tenia por que cambiar y
comparar el resultado** encontro el fallo mas grave del dia (P-36). No lo
encontro ningun test ni ninguna lectura del codigo. Conviene hacerlo a proposito
de vez en cuando.

## Reglas de trabajo acordadas

- Aprobar cada instalacion y descarga antes de hacerla, diciendo que es y
  cuanto pesa.
- Commits pequenos, mensajes claros, **sin coautoria de Claude**.
- Escribir el codigo completo antes de ejecutarlo; pruebas pequenas que no dejen
  residuos estan bien.
- Todo lo importante va a `docs/`: cambios, errores, logica, decisiones y el
  comportamiento y las estadisticas de los agentes.
- Todo gratis y en CPU de un Mac con Apple Silicon; si hiciera falta GPU, Kaggle
  o Colab.
