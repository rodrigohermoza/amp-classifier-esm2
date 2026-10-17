# Estado del proyecto y como continuar

**Ultima actualizacion: 2026-10-17 (entorno arreglado, `split` corrido cuatro
veces, agrupamiento corregido dos veces).** Este archivo es el punto de entrada
para retomar el trabajo. Lo que esta hecho, lo que falta y lo que hay que decidir.

---

## Resumen en una linea

Codigo completo y auditado, publicado en
<https://github.com/rodrigohermoza/amp-classifier-esm2>. **`data` y `split` ya
corren con fuentes reales.** La auditoria de fuga de `split` **todavia dice
`publishable: false`**, y el ultimo diagnostico (P-32) dice que parte de lo que
marca no es fuga sino un defecto de la propia metrica. **Lo siguiente es aplicar
el arreglo de P-32, correr `split` una vez mas y, si sale limpio, seguir con
`baselines` en adelante.** Hay una fase 2 de generacion diseñada y sin
implementar en [17-diseno-generacion.md](17-diseno-generacion.md).

---

## Lo que hay que hacer primero al retomar (resumen de 30 segundos)

1. Aplicar el arreglo propuesto en **P-32**: que la auditoria exija el mismo
   `min_coverage` que el agrupamiento. Son pocas lineas en
   `_max_identity_against`, en `src/ampesm/split.py`.
2. `make split` (12,5 min) y leer `leakage.global_p_value` en **las dos**
   auditorias y `publishable_reasons`.
3. Si `publishable` es true, seguir con `baselines` -> `embeddings` -> `train`
   -> `evaluate` -> `interpret` (~35 min en total).
4. Si val sigue significativa, **no seguir ajustando**: declararlo en el informe
   y dejar M-01 (curva de umbral de identidad) como la forma publicable de
   cerrarlo. Ajustar parametros hasta que el numero salga bien es exactamente el
   fallo que este repositorio existe para evitar.

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

### Que dijo `split` en sus cuatro corridas (2026-10-17)

Dos de los cuatro controles salieron perfectos y se mantuvieron en las cuatro
corridas: **deriva de tasa positiva 0.000** (las tres divisiones en 0.568, el
guard de P-26 sin nada que decir) y **MMseqs2 usado** (18-8cc5c). El tercero,
`publishable`, mentia (P-29) y ahora funciona. El cuarto, la fuga, es la historia
de la sesion:

| Corrida | Cambio | Clusteres | Exceso test | Exceso val |
|---|---|---|---|---|
| 1 | tal cual estaba | 8 173 | +50.2 (p=0) | +47.2 (p=0) |
| 2 | `--cluster-mode 1` (P-29) | 7 886 | +34.2 (p=0) | +32.2 (p=0) |
| 3 | fusion por contencion (P-30/D-20) | 7 438 | +20.2 (p=8.7e-05) | +30.2 (p=0) |
| 4 | arreglado el `break` de P-31 | 7 178 | **+10.2 (p=0.025)** | **+28.2 (p=0)** |

El esperado por azar es 23.81 marcadas de 500 en todas las corridas. La corrida 4
es **anterior** al arreglo de P-31b (union-find sobre clusteres), que esta en el
codigo y **sin medir**: la proxima corrida dara numeros algo mejores que estos.

Lo que falta por hacer esta diagnosticado y es **P-32**: la auditoria descarta la
cobertura y el agrupamiento la exige, asi que marca parejas como `CCISG` (cinco
residuos) contra cualquier peptido rico en cisteina. No estan midiendo lo mismo.

## Numeros actuales

| Metrica | Valor |
|---|---|
| Commits | 29 |
| Codigo en `src/` | 4 631 lineas, 16 modulos |
| Tests | **145 casos, todos pasan** (7,4 s) |
| Documentacion | 5 308 lineas, 18 archivos |
| Datos procesados | 18 644 peptidos (10 594 + / 8 050 -) |
| Clusteres (ultima corrida) | 7 178 (5 621 de 1), 92 con ambas clases |
| Lint | ruff limpio |
| Repositorio clonable | ~1 MB |
| Entorno virtual | 1,0 GB, **en `~/.venvs/ampesm`, fuera de iCloud** (D-19) |
| Cache de pesos ESM-2 | 160 MB (`.cache/`, ignorada por git) |
| Problemas registrados | 33 (mas P-10a-h y P-25b/P-31b como subapartados); **2 abiertos: P-23 (a proposito) y P-32 (siguiente tarea)** |
| Decisiones registradas | 20 (D-01 a D-20) |
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

Quedan **dos problemas abiertos**:

- **P-23**, deliberadamente abierto: no afecta a ningun numero publicable porque
  MMseqs2 es obligatorio y el respaldo marca `publishable: false`.
- **P-32**, que es la siguiente tarea: la auditoria descarta la cobertura que el
  agrupamiento exige, asi que marca parejas que no son homologas.

### 2. `make split` — CORRIDO 4 VECES, aun no limpio (bloqueante)

Historia completa de las cuatro corridas en la seccion "Que dijo `split`" de
arriba, y el detalle en P-29, P-30, P-31 y P-32.

**Estado:** `publishable: false`, con fuga significativa en test (+10.2,
p = 0.025) y en val (+28.2, p = 0). **Siguiente paso: aplicar P-32** (que la
auditoria exija el mismo `min_coverage` que el agrupamiento) y volver a correr.

Que mirar en `data/processed/split_report.json`:

- `publishable` debe ser `true`. Si es `false`, `publishable_reasons` dice por
  que, en texto. Ese campo ya tiene en cuenta la fuga (antes no, P-29).
- `train_test_identity_audit.leakage.global_p_value` por encima de 0.05 y
  `excess` cerca de 0. **Ese** es el numero de fuga. NO leer
  `fraction_over_threshold` (dara 100% y es ruido de fondo de la metrica) ni
  `excess_over_null_deprecated`. Historia en P-12, P-15 y P-18.
- `train_val_identity_audit.leakage.global_p_value`: lo mismo para val. Si este
  sale significativo y el de test no, el umbral se eligio sobre un conjunto
  optimista y las metricas de test estan sesgadas igualmente (P-27).
- `contained_merge`: cuantas secuencias eran subcadena de otra y cuantos
  clusteres se fusionaron por ello (P-30).
- `positive_rate_drift.exceeds_tolerance` debe ser `false`. En las cuatro
  corridas salio 0.000, perfecto.
- `splits.*.length_ks_pvalue`: si es muy bajo en `test`, la longitud sigue
  separando las clases dentro del split y conviene reportar la ablacion sin
  `length` ni `charge_density`.
- `mixed_label_clusters`: clusteres con ambas clases. Ultima corrida: 92.

**Si val sigue significativa despues de P-32:** parar de ajustar. El residuo ya
es homologia real, la respuesta es declararlo en el informe y usar M-01 (curva de
umbral de identidad) como forma publicable de cerrarlo. Seguir tocando parametros
hasta que el numero salga bien es el fallo que este repositorio existe para
evitar, y seria una ironia terminar cometiendolo.

### 3. `make baselines`, `make embeddings`, `make train`
Los pesos de ESM-2 ya estan en cache. Con 18 644 peptidos, los embeddings
tardaran del orden de 10-20 minutos en CPU para los dos modelos. `baselines`
entrena ahora **4 representaciones x 3 modelos**, porque D-16 añadio
`comp_physchem_nolen`.

Los tiempos de todas las etapas estan medidos en la tabla de "Coste medido" de
arriba. `split` tarda 12,5 min con 18 644 peptidos, y las palancas para bajarlo
(con su coste en calidad) estan ahi mismo.

### 4. `make evaluate` y `make interpret`
Leer `results/metrics/comparison_embeddings_vs_baseline.json`, que es la unica
inferencia del proyecto. El resto de la tabla es descriptiva.

### 5. Rellenar la documentacion con los resultados
- Tablas de resultados de `README.md` y `README.en.md` (hoy con guiones).
- Los `[PENDIENTE]` de `13-portafolio.md`, que son dos: el paso 3 en simple y en
  tecnico.
- Actualizar `01-objetivos.md` si la expectativa escrita (margen modesto, quiza
  intervalos solapados) resulta equivocada. **No reescribir la expectativa para
  que encaje con el resultado**: esta puesta ahi justamente para poder
  contrastarla.

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

**Abierta:**

1. **¿Se hace la regresion de MIC con los MIC de DRAMP?** El reexamen de P-02
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
6. **Las etapas 3 a 7 nunca se han ejecutado con datos reales.** `split` ya
   si, cuatro veces.
7. **La division todavia no esta limpia.** `publishable: false` con fuga
   significativa en val. Diagnosticado (P-32) y con arreglo propuesto, pero
   mientras siga asi **ninguna metrica de test es publicable**, y esa es
   precisamente la funcion del campo.
8. **La metrica de identidad no sirve para peptidos de 5-7 residuos**, y hay 536
   en el conjunto. Un alineador global con huecos da identidad 1.00 entre
   `CCISG` y cualquier peptido rico en cisteina. Afecta a la auditoria, no al
   clasificador, pero contamina el p-valor de fuga (P-32).
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
cobertura), o la guarda mide otra cosa que lo garantizado — eso es P-32 y hay que
arreglarlo. Pero no deben usar *el mismo alineador*, porque entonces la guarda
hereda sus puntos ciegos y no puede contradecirlo — eso es D-20 y hay que
evitarlo.

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
