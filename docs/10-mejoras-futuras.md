# Mejoras futuras

Ordenadas por relacion entre esfuerzo y valor. Las primeras son las que de verdad
cambiarian las conclusiones; las ultimas son vistosas pero aportan poco.

> **Estado al 2026-10-11.** Lo que antes estaba en esta lista y ya esta hecho:
> auditar `val` contra `train` (P-27), validar la deriva de la tasa de positivos
> (P-26), el hash de `features` (P-24) y la guarda de entorno (P-25). La
> **fase 2 de generacion** dejo de ser una idea de esta lista y tiene diseño
> propio aprobado en [17-diseno-generacion.md](17-diseno-generacion.md).
>
> De lo que sigue aqui, el orden recomendado despues de tener resultados es
> **M-01, luego los negativos barajados, luego M-02**. Y una que **no** se
> recomienda hacer: bajar el prefiltro del respaldo de clustering (P-23). No
> afecta a ningun numero publicable, porque MMseqs2 es obligatorio y el respaldo
> marca `publishable: false`, y el mismo esfuerzo invertido en M-01 da una figura
> publicable en vez de una nota de higiene. M-04 queda fuera de alcance mientras
> la fase 2 este viva (D-17), y ademas su via buena ya no es recuperar DBAASP
> sino parsear los MIC que DRAMP trae en 5 020 de sus filas (P-02).

### Una mejora que cerraria la asimetria de D-16
La ablacion `comp_physchem_nolen` quita la longitud al baseline y no a ESM-2, lo
que sesga esa lectura a favor de los embeddings. Para cerrarlo haria falta un
control del lado de los embeddings: **residualizar cada dimension contra la
longitud** (ajustar una regresion de la dimension sobre la longitud y quedarse
con el residuo) antes de entrenar. Esfuerzo bajo, y convierte una advertencia
escrita en un numero. Tiene sentido hacerla junto con M-01.

## Lo primero despues de tener resultados

### M-01 es ahora obligatoria, no opcional

La division en uso tiene **fuga por homologia real y declarada** (D-22), y la
salvedad que viaja con cada metrica remite a esta curva. Ya no es "la figura mas
instructiva del proyecto": es la forma en que el proyecto cierra su propio
defecto principal. Coste medido: `split` son ~12 min por umbral y los embeddings
**no** hay que recalcularlos (dependen solo de las secuencias, no de la
division), asi que son cinco corridas de split + entrenamiento + evaluacion.

### M-07 — Fusionar las parejas cortas que MMseqs2 no ve (solo como sensibilidad)

Fuerza bruta sobre las parejas con `min(longitud) < 15` y razon de longitudes
>= 0.5: **25 124 791** parejas, ~13 min medidos a 33 000 alineamientos/s. Fusiona
los clusteres de las que pasen identidad >= 0.4 y cobertura >= 0.5. Cerraria el
punto ciego que deja el prefiltro de k-meros de MMseqs2 en peptidos de 8-14
residuos, que es donde se concentra la fuga que queda.

**No puede ser el resultado principal, y el motivo esta en D-22:** esa fusion
mide con Biopython por fuerza bruta y la auditoria mide con Biopython con
prefiltro, asi que la guarda seria una version mas debil de lo vigilado y daria
verde por construccion. Si se hace, va **etiquetada** como analisis de
sensibilidad: "metricas con la redundancia corta eliminada por nuestra propia
regla, sin comprobacion independiente en esa clase". Informa; no certifica.

## Alto valor, bajo esfuerzo

### M-01 — Curva de rendimiento frente al umbral de identidad
Repetir todo con umbrales de 30%, 40%, 50%, 70% y 90%, y graficar ROC-AUC contra
umbral. La curva descendente **es** la medida directa de cuanta fuga tenia la
division aleatoria. Seria probablemente la figura mas instructiva del proyecto y
cuesta poco: la tuberia ya esta parametrizada.

### M-02 — Evaluacion cruzada por fuente
Entrenar solo con APD3 y evaluar en los positivos de UniProt, y al reves. Si el
rendimiento se cae, el modelo aprendio procedencia. Es la prueba mas directa del
sesgo descrito en [01-objetivos.md](01-objetivos.md) y solo requiere reordenar los
splits.

### M-03 — Capas intermedias de ESM-2
Probar la extraccion desde varias capas, no solo la ultima. Las ultimas capas se
especializan en la tarea de enmascaramiento y a veces pierden informacion util.
Cuesta una pasada por el modelo guardando todos los estados ocultos.

## Alto valor, esfuerzo medio

### M-04 — Recuperar DBAASP y hacer la regresion de MIC
Ver [09-problemas-encontrados.md](09-problemas-encontrados.md), P-02, para las
vias. Habilita la regresion, que es una tarea mas difícil y mas util que la
clasificacion: predecir *cuanto* de potente, no solo *si*.

### M-05 — Negativos mas duros
Los negativos actuales son faciles. Alternativas:
- Peptidos con actividad biologica conocida que **no** sea antimicrobiana
  (hormonas, neuropeptidos, inhibidores de proteasas). Obliga al modelo a
  distinguir "activo antimicrobiano" de "activo", no "activo" de "inerte".
- Secuencias barajadas que conserven la composicion exacta de los positivos.
  Elimina la composicion como atajo y mide si el modelo usa el **orden**. Esta
  es, conceptualmente, la prueba mas afilada: si ESM-2 no supera a la composicion
  con negativos barajados, no esta aportando nada sobre el orden de los residuos.

### M-06 — Calibracion de probabilidades
Platt scaling o regresion isotonica ajustada en validacion. Importa si el modelo
se usara para priorizar candidatos: un ordenamiento sirve, pero una probabilidad
calibrada permite decidir cuantos sintetizar.

## Esfuerzo alto

### M-07 — Fine-tuning de ESM-2
Ajustar los pesos en vez de usarlo congelado. Es la extension con mas potencial
de mejora y la que mas riesgo de sobreajuste tiene con unos pocos miles de
ejemplos. Necesita GPU: Kaggle da ~30 h semanales de T4 gratis. Con LoRA sobre el
modelo de 150M seria viable.

### M-08 — Atencion por residuo para interpretacion
Extraer pesos de atencion o gradientes por posicion para ver que residuos pesan
en la prediccion. Atractivo visualmente, pero hay que tener cuidado: la atencion
no es explicacion, y hay literatura abundante sobre ese malentendido. Como
mucho, generaria hipotesis.

### M-09 — Modelo multi-tarea con hemolisis y toxicidad
Un AMP util tiene que matar al microbio sin matar celulas humanas. Predecir
actividad y toxicidad juntas es mas parecido al problema real. Necesita otras
bases de datos (HemoPI, ToxinPred) y multiplica el trabajo de curado.

## El algoritmo genetico

### M-10 — Reproducir la idea de Sup'Biotech, como experimento negativo
Ver [08-decisiones.md](08-decisiones.md), D-10, para el encuadre completo.

El valor no esta en generar secuencias: esta en **mostrar que optimizar contra un
clasificador produce ejemplos adversarios**. Un algoritmo genetico con el
clasificador como funcion de aptitud converge a secuencias con puntuacion de 0.99
y propiedades absurdas, porque explota los puntos ciegos del modelo en vez de
descubrir biologia.

Lo que se reportaria:
1. La trayectoria de aptitud subiendo mientras el realismo de las secuencias cae.
2. Las mitigaciones probadas: penalizar la distancia a la distribucion de
   entrenamiento, restringir la composicion a rangos observados, usar un comite
   de modelos en vez de uno solo.
3. La conclusion metodologica: un clasificador es un filtro, no una funcion
   objetivo.

Esto no produce candidatos para sintesis y no es el objetivo. Es una
demostracion de los limites del metodo, que es lo que hace falta entender antes
de que un diseno generativo signifique algo.
