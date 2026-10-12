# Estado del proyecto y como continuar

**Ultima actualizacion: 2026-10-11 (las dos auditorias cerradas).** Este archivo es el punto de entrada para
retomar el trabajo. Lo que esta hecho, lo que falta y lo que hay que decidir.

---

## Resumen en una linea

Codigo completo y auditado, tuberia verificada de punta a punta con datos
sinteticos. **Nunca se ha corrido contra datos reales.** Falta eso y rellenar
los resultados.

---

## Estado por etapa

| Etapa | Codigo | Probado con sinteticos | Corrido con datos reales |
|---|---|---|---|
| `data` (descarga + limpieza) | si | **no** | **no** |
| `split` (clusteres + division) | si | si, con MMseqs2 | no |
| `baselines` (9 modelos) | si | si, XGBoost incluido | no |
| `embeddings` (ESM-2 8M y 35M) | si | si, pesos ya en cache | no |
| `train` (logreg y MLP) | si | si | no |
| `evaluate` (metricas + IC + figuras) | si | si | no |
| `interpret` (UMAP + errores) | si | si, 8 figuras | no |

**La unica etapa que nunca se ha ejecutado es `data`**, y es justamente donde
esta el riesgo que queda: el parseo del XLSX real de DRAMP (nombres de columna
cambiantes, por eso la deteccion es por patron) y la cabecera espuria del FASTA
de APD3 (P-01, con test de regresion y aviso si el conteo no es exactamente 1).

## Numeros actuales

| Metrica | Valor |
|---|---|
| Commits | 14 |
| Codigo en `src/` | 4 065 lineas, 16 modulos |
| Tests | **117 casos, todos pasan** |
| Documentacion | 2 110 lineas, 14 archivos |
| Lint | ruff limpio |
| Repositorio clonable | ~1 MB |
| Entorno `.venv` | 1,0 GB |
| Cache de pesos ESM-2 | 160 MB (`.cache/`, ignorada por git) |
| Problemas registrados | 24 (P-01 a P-24); 2 abiertos: P-23 y P-24 |
| Decisiones registradas | 15 (D-01 a D-15) |

## Dependencias del sistema ya instaladas

```bash
brew install mmseqs2   # 18-8cc5c. OBLIGATORIO para numeros publicables (D-11)
brew install libomp    # sin ella XGBoost no carga en macOS (P-09)
```

---

## Lo que falta, en orden

### 1. Correr `make data` (bloqueante)
Descarga ~5 MB de APD3, DRAMP y UniProt. Es la primera vez que el parser toca
datos reales. Que mirar en `data/processed/clean_report.json`:

- `parse[apd3].dropped_empty` **debe ser exactamente 1** (P-01). Si no, el
  formato de APD3 cambio.
- `parse[dramp].sequence_column` y `n_columns`: confirmar que encontro la
  columna correcta.
- `deduplicate.label_conflicts`: secuencias presentes como positivas y
  negativas. Es una estimacion por lo bajo del ruido en los negativos y un
  numero interesante por si mismo.
- `length_matching.deficit`: si es alto, el emparejado de longitudes no se
  cumple bien.

### 2. `make split`
Usara MMseqs2. Que mirar en `data/processed/split_report.json`:

- `publishable` debe ser `true`.
- `train_test_identity_audit.leakage.global_p_value` por encima de 0.05 y
  `excess` cerca de 0. **Ese** es el numero de fuga. NO leer
  `fraction_over_threshold` (dara ~30% y es ruido de fondo) ni
  `excess_over_null_deprecated`. Historia en P-12, P-15 y P-18.
- `splits.*.length_ks_pvalue`: si es muy bajo en `test`, la longitud sigue
  separando las clases dentro del split y conviene reportar la ablacion sin
  `length` ni `charge_density`.
- `mixed_label_clusters`: clusteres con ambas clases.

### 3. `make baselines`, `make embeddings`, `make train`
Los pesos de ESM-2 ya estan en cache. Con ~16 000 peptidos, los embeddings
tardaran del orden de 10-15 minutos en CPU para los dos modelos.

**Tiempos medidos de `split`** (incluye la auditoria): 16 s con 600 peptidos,
113 s con 3 900. El coste de la auditoria esta topado por `split.audit.max_test`
(500) para que no crezca de forma cuadratica (P-19); si con el conjunto real
tarda demasiado, bajar `null_replicates` de 20 a 10 lo reduce a la mitad.

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

### 6. Verificacion: las dos auditorias estan cerradas

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
- **P-24**: cambiar `features.*` y correr solo `interpret` sin reentrenar no
  avisa. Un hash de esa seccion en `results/models/*.json` lo cerraria.

Lo que las auditorias **confirmaron como correcto**, ejecutandolo: la escala de
Boman completa (20 valores, 0 discrepancias), el manejo de NaN en
`all_metrics` y `bootstrap_metrics`, que las etiquetas preregistradas coinciden
con los nombres reales de archivo, la seleccion por validacion en `evaluate` e
`interpret`, que la siembra no rompe el invariante y da 70/15/15 exactos, que
los 5 casos de division degenerada abortan sin escribir archivos, las ocho
guardas de P-10 con mensajes claros y sin trazas, que `ampesm all` y `make all`
corren lo mismo, dos tuberias completas de principio a fin con 16 JSON validos
y 8 figuras, y que los tests pasan aunque se retire el `.pth` del editable.

Tambien quedan dos puntos de rigor del informe metodologico, no de correccion:

- **La auditoria solo cubre test, no val.** `val` fija los hiperparametros y el
  umbral, asi que la homologia val-train tambien infla la seleccion.
- **No se valida la deriva de la tasa positiva entre splits.** Hay un caso que
  pasa todas las comprobaciones con train al 6% de positivos y test al 29%.

---

## Decisiones abiertas para Rodrigo

1. **¿Publicar el repositorio ahora o cuando haya resultados?** Recomendacion:
   cuando haya resultados, para que la primera impresion no sea un repositorio
   con tablas vacias.
2. **Nombre publico.** La carpeta es `Peptidos-antimicrobianos`; para un CV en
   ingles, `amp-classifier-esm2` dice mas de un vistazo.
3. **¿Se intenta recuperar DBAASP?** Es lo que habilita la regresion de MIC
   (P-02, D-07). Vias: descarga manual desde un navegador, un espejo en Zenodo,
   o usar los valores de MIC que DRAMP trae en sus metadatos.

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
4. **El emparejado de longitudes es global, no por split.** El reporte incluye
   ahora un KS por split para detectarlo, pero si sale mal hay que decidir entre
   emparejar por split o publicar la ablacion sin `length`.
5. **La etapa `data` nunca se ejecuto.**

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
