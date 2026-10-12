# Evaluacion

## Metricas y por que cada una

| Metrica | Que mide | Por que esta |
|---|---|---|
| **ROC-AUC** | Probabilidad de ordenar un positivo por encima de un negativo | Comparable con la literatura. Insensible al umbral. |
| **PR-AUC** | Area bajo precision-recall | Mas informativa que ROC cuando las clases estan desbalanceadas, y aqui lo estan (~3 negativos por positivo). |
| **MCC** | Correlacion de Matthews | La metrica de un solo numero mas honesta en clasificacion binaria desbalanceada: solo es alta si las cuatro casillas de la matriz de confusion son buenas. |
| **F1** | Media armonica de precision y recall | Comparabilidad; se reporta aunque ignora los verdaderos negativos. |
| **Sensibilidad y especificidad** | Por clase | Para ver *de que lado* falla el modelo. |

El umbral para las metricas que lo necesitan (MCC, F1, matriz de confusion) se
elige **en validacion**, maximizando MCC, y se aplica sin cambios a test. Elegir
el umbral en test es una forma sutil y frecuente de inflar resultados.

## Intervalos de confianza por bootstrap

Un numero sin incertidumbre no permite comparar modelos. Procedimiento:

1. Remuestrear el conjunto de test con reemplazo, 2 000 veces.
2. Recalcular cada metrica en cada remuestreo.
3. Reportar el percentil 2.5 y el 97.5 (IC del 95%).

**Detalle importante:** el remuestreo es **por cluster, no por secuencia**. Las
secuencias dentro de un cluster son parientes; tratarlas como independientes
produce intervalos demasiado estrechos. Remuestrear clusteres completos respeta
la estructura de dependencia de los datos y da intervalos mas anchos y mas
honestos.

Para comparar dos modelos se usa el bootstrap de la **diferencia** emparejada
(mismos remuestreos para ambos modelos). Si el IC del 95% de la diferencia
incluye el cero, la conclusion es "no se distingue", no "el modelo A es mejor".
Comparar los intervalos individuales por separado es estadisticamente incorrecto
y es un error que se ve a menudo.

## Figuras

Todas en `results/figures/`, generadas por codigo, sin retoques manuales:

| Figura | Contenido |
|---|---|
| `roc_curves.png` | ROC de todos los modelos con banda de confianza |
| `pr_curves.png` | Precision-recall equivalente |
| `confusion_matrices.png` | Matriz del mejor baseline y del mejor modelo con ESM-2 |
| `metrics_comparison.png` | Barras con barras de error (IC 95%) por modelo y metrica |
| `calibration.png` | Curva de calibracion: ¿una probabilidad de 0.8 significa 80%? |
| `length_distribution.png` | Longitud por clase y por fuente: control de sesgo |
| `umap_embeddings.png` | UMAP coloreado por clase **y** por fuente |
| `error_analysis.png` | Descriptores de los aciertos frente a los fallos |

## Analisis de errores

Es la seccion que mas dice sobre si el modelo aprendio biologia o procedencia.
Preguntas que se responden con datos:

1. **¿Los falsos negativos son AMPs atipicos?** Se espera que lo sean: peptidos
   anionicos, ricos en prolina, no helicoidales. Son familias reales de AMPs que
   no siguen el mecanismo cationico-anfipatico dominante. Si el modelo falla
   justo ahi, se entiende lo que aprendio.
2. **¿Los falsos positivos son negativos mal etiquetados?** Se inspeccionan a
   mano las entradas de UniProt con mayor puntuacion. Si son peptidos cationicos
   secretados de funcion desconocida, es muy posible que sean antimicrobianos no
   anotados, y entonces la "precision" real es mejor que la medida.
3. **¿El error se concentra por fuente?** Si el rendimiento en los positivos de
   UniProt es mucho peor que en los de APD3, el modelo depende de la procedencia.
   Esta es la prueba mas directa del sesgo descrito en
   [01-objetivos.md](01-objetivos.md).
4. **¿El error depende de la longitud?** Los peptidos muy cortos (5–10 residuos)
   tienen poca senal; se espera peor rendimiento ahi.

## Comparacion con un predictor publicado

Objetivo opcional: pasar el conjunto de test por un servidor web publico
(CAMPR4, AMPfun o similar) y comparar.

Hay que ser claro sobre la limitacion metodologica: esos predictores se
entrenaron con datos que **con casi total seguridad incluyen parte de nuestro
test**, porque todos salen de las mismas bases. Asi que la comparacion es
favorable a ellos por construccion y solo sirve como referencia aproximada, no
como evaluacion justa. Se reporta con esa advertencia escrita al lado, o no se
reporta.

## Formato de los resultados

`results/metrics/` contiene un JSON por modelo con la estructura:

```json
{
  "model": "esm2_t12_35M+logreg",
  "split": "test",
  "threshold": 0.43,
  "threshold_selected_on": "val",
  "metrics": {
    "roc_auc": { "point": 0.000, "ci_low": 0.000, "ci_high": 0.000 },
    "pr_auc":  { "point": 0.000, "ci_low": 0.000, "ci_high": 0.000 },
    "mcc":     { "point": 0.000, "ci_low": 0.000, "ci_high": 0.000 }
  },
  "confusion_matrix": { "tn": 0, "fp": 0, "fn": 0, "tp": 0 },
  "n_test": 0,
  "n_clusters_test": 0,
  "bootstrap_n": 2000,
  "bootstrap_unit": "cluster"
}
```

Y `results/metrics/summary.csv` junta todo en una tabla, que es lo que se pega en
el README.

## Tres cosas que hay que leer antes que la tabla de metricas

Las metricas no se interpretan solas. Antes de mirar `summary.csv` conviene
abrir `split_report.json` y comprobar tres campos, porque los tres pueden
invalidar o reencuadrar la tabla entera:

1. **`train_test_identity_audit.leakage.global_p_value` y el de
   `train_val_identity_audit`.** El segundo importa tanto como el primero y es
   menos obvio: `val` fija el umbral que se aplica a `test`, asi que homologia
   val-train sesga cinco de las siete metricas **aunque test este limpio**
   (P-27). Y en los dos casos el numero que se lee es el p-valor y el exceso,
   nunca `fraction_over_threshold`.
2. **`positive_rate_drift.exceeds_tolerance`.** Si es `true`, PR-AUC, F1 y
   precision **no son comparables entre splits** porque dependen de la
   prevalencia, y el umbral elegido en `val` esta descalibrado para `test`
   (P-26).
3. **`publishable`.** Si es `false` se uso el respaldo de clustering en Python
   y los numeros no son publicables (D-11).

## La ablacion de longitud, y como no leerla al reves

La tabla incluye `comp_physchem_nolen`: el baseline sin `length` ni
`charge_density`. Existe porque los positivos miden 22.9 residuos de media y
los negativos 30.6, y el emparejado de longitudes no puede corregirlo (D-16).

**La diferencia entre `comp_physchem` y `comp_physchem_nolen` es la medida de
cuanto del baseline era longitud.** Esa es su lectura.

Lo que **no** es: un baseline limpio contra el que comparar ESM-2. La ablacion
quita la longitud al baseline y no a los embeddings, asi que esa comparacion
esta sesgada a favor de ESM-2. La comparacion preregistrada sigue siendo contra
`comp_physchem` completo.

## Compromiso previo

Para que esto no sea narrar lo que salio bien a posteriori, queda escrito antes de
correr nada:

- El modelo principal es **ESM-2 35M + regresion logistica**.
- El baseline principal es **composicion + fisicoquimicos con XGBoost**
  (el completo, no la ablacion).
- La metrica principal es **MCC** sobre test, con IC por bootstrap de clusteres.
- Test se evalua **una sola vez**, al final.

Si el resultado es que el baseline gana, el README lo dira en el titular.
