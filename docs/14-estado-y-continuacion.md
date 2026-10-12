# Estado del proyecto y como continuar

**Ultima actualizacion: 2026-10-11 (repositorio publicado, `data` corrido,
los cuatro arreglos previos a la tuberia hechos, fase 2 diseñada).** Este archivo es el punto de entrada para retomar el trabajo.
Lo que esta hecho, lo que falta y lo que hay que decidir.

---

## Resumen en una linea

Codigo completo y auditado, publicado en
<https://github.com/rodrigohermoza/amp-classifier-esm2>. **`data` ya corrio con
fuentes reales: 18 644 peptidos.** **Lo siguiente es `make split` y el resto de
la tuberia**, que es lo unico que queda para tener resultados. Hay una fase 2 de
generacion diseñada y sin implementar en
[17-diseno-generacion.md](17-diseno-generacion.md).

---

## Estado por etapa

| Etapa | Codigo | Probado con sinteticos | Corrido con datos reales |
|---|---|---|---|
| `data` (descarga + limpieza) | si | si | **si, 2026-10-11** |
| `split` (clusteres + division) | si | si, con MMseqs2 | no |
| `baselines` (9 modelos) | si | si, XGBoost incluido | no |
| `embeddings` (ESM-2 8M y 35M) | si | si, pesos ya en cache | no |
| `train` (logreg y MLP) | si | si | no |
| `evaluate` (metricas + IC + figuras) | si | si | no |
| `interpret` (UMAP + errores) | si | si, 8 figuras | no |

**`data` ya corrio** y los cuatro controles salieron en orden (detalle en
[09](09-problemas-encontrados.md), seccion "Que dijo la primera corrida real"):
la cabecera espuria de APD3 fue exactamente 1, DRAMP encontro la columna
correcta, 271 conflictos de etiqueta y un `deficit` de emparejado de longitudes
de **23 732**, que es el numero que motivo D-16. El riesgo que queda esta ahora
en `split` con 18 644 secuencias: es la primera vez que MMseqs2 y la auditoria
de fuga ven este volumen.

## Numeros actuales

| Metrica | Valor |
|---|---|
| Commits | 22 |
| Codigo en `src/` | 4 445 lineas, 16 modulos |
| Tests | **129 casos, todos pasan** |
| Documentacion | 4 551 lineas, 18 archivos |
| Datos procesados | 18 644 peptidos (10 594 + / 8 050 -) |
| Lint | ruff limpio |
| Repositorio clonable | ~1 MB |
| Entorno `.venv` | 1,0 GB |
| Cache de pesos ESM-2 | 160 MB (`.cache/`, ignorada por git) |
| Problemas registrados | 27 (P-01 a P-27, mas P-25b); **1 abierto: P-23** |
| Decisiones registradas | 18 (D-01 a D-18) |
| Repositorio publico | <https://github.com/rodrigohermoza/amp-classifier-esm2> |

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

Queda **un solo problema abierto, P-23**, y esta deliberadamente abierto: no
afecta a ningun numero publicable porque MMseqs2 es obligatorio y el respaldo
marca `publishable: false`.

### 2. `make split` (siguiente, bloqueante)
Usara MMseqs2. Que mirar en `data/processed/split_report.json`:

- `publishable` debe ser `true`.
- `train_test_identity_audit.leakage.global_p_value` por encima de 0.05 y
  `excess` cerca de 0. **Ese** es el numero de fuga. NO leer
  `fraction_over_threshold` (dara ~30-50% y es ruido de fondo) ni
  `excess_over_null_deprecated`. Historia en P-12, P-15 y P-18.
- `train_val_identity_audit.leakage.global_p_value`: lo mismo para val. Si este
  sale significativo y el de test no, el umbral se eligio sobre un conjunto
  optimista y las metricas de test estan sesgadas igualmente (P-27).
- `positive_rate_drift.exceeds_tolerance` debe ser `false`. Si es `true`,
  PR-AUC, F1 y precision no son comparables entre splits y el umbral de val
  esta descalibrado para test (P-26).
- `splits.*.length_ks_pvalue`: si es muy bajo en `test`, la longitud sigue
  separando las clases dentro del split y conviene reportar la ablacion sin
  `length` ni `charge_density`.
- `mixed_label_clusters`: clusteres con ambas clases.

### 3. `make baselines`, `make embeddings`, `make train`
Los pesos de ESM-2 ya estan en cache. Con 18 644 peptidos, los embeddings
tardaran del orden de 10-20 minutos en CPU para los dos modelos. `baselines`
entrena ahora **4 representaciones x 3 modelos**, porque D-16 añadio
`comp_physchem_nolen`.

**Tiempos medidos de `split`** (incluye la auditoria): 16 s con 600 peptidos,
113 s con 3 900. Con 18 644 habra que ver: es la primera vez a este volumen, y
ademas P-27 **duplica** el coste de la auditoria porque ahora audita val tambien
(se desactiva con `split.audit.audit_val: false` si hiciera falta). El
coste de la auditoria esta topado por `split.audit.max_test` (500) para que no
crezca de forma cuadratica (P-19); si tarda demasiado, bajar `null_replicates`
de 20 a 10 lo reduce a la mitad.

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
6. **Las etapas 2 a 7 nunca se han ejecutado con datos reales.**

---

## Lo que mas conviene recordar de esta sesion

Tres intentos hicieron falta para construir un detector de fuga que detectara
fuga. Los tres primeros devolvian numeros plausibles y no encontraban nada, ni
con parientes al 80% de identidad plantados a mano. Ninguno fallaba de forma
visible.

Y dos de los cinco fallos nuevos los **introdujeron los arreglos anteriores**:
el corte temprano de la auditoria reintrodujo exactamente el bug que P-03
acababa de corregir, porque se copio la idea sin volver a comprobar sus
condiciones.

De ahi D-15 y el habito que vale la pena conservar: **plantar el efecto que
buscas y exigir que tu herramienta lo encuentre**, antes de creerte el numero
que te devuelve.

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
