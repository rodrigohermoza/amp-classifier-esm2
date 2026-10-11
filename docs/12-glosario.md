# Glosario

Pensado en dos direcciones: biologia para quien viene de ML, y ML para quien
viene de biologia.

## Biologia, para lectores de machine learning

**Aminoacido.** La unidad de una proteina. Hay 20 estandar, cada uno con una
letra: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y. Para un modelo
de ML, un vocabulario de tamano 20.

**Peptido / proteina.** La misma cosa a distinta escala. Por convencion, peptido
es corto (menos de ~50 residuos) y proteina es largo. Aqui se trabaja con 5–50.

**Residuo.** Un aminoacido dentro de la cadena. "Posicion 7" y "residuo 7" son
lo mismo.

**AMP (peptido antimicrobiano).** Peptido corto que mata o inhibe microbios,
normalmente rompiendo su membrana.

**Cationico.** De carga positiva neta. Importa porque la membrana bacteriana es
mas negativa que la de las celulas animales, asi que un peptido positivo es
atraido selectivamente hacia la bacteria.

**Anfipatico.** Con una cara hidrofobica y otra hidrofilica. Es la arquitectura
tipica de los AMPs helicoidales: la cara hidrofobica se mete en el interior
lipidico de la membrana y la hidrofilica se queda en el agua.

**Helice alfa.** Estructura en espiral con un giro de ~100° por residuo. Por eso
el momento hidrofobico se calcula con ese angulo.

**MIC (concentracion minima inhibitoria).** La concentracion mas baja que impide
el crecimiento visible del microbio. Mas baja = mas potente. Es la medida
cuantitativa de actividad, y la variable de la regresion diferida.

**Hemolisis.** Destruccion de globulos rojos. El efecto secundario que hace
inutil a muchos AMPs potentes: matan la bacteria y tambien las celulas del
paciente.

**Homologia.** Parentesco evolutivo entre secuencias. Dos secuencias homologas se
parecen porque descienden de un ancestro comun. Es la razon de que una division
aleatoria filtre informacion.

**Identidad de secuencia.** Porcentaje de posiciones iguales en un alineamiento.
40% de identidad es aproximadamente el limite por debajo del cual el parentesco
deja de ser evidente.

**UniProt / Swiss-Prot / TrEMBL.** UniProt es la base de referencia de secuencias
de proteinas. Swiss-Prot es su parte revisada a mano (`reviewed:true`); TrEMBL es
la parte anotada automaticamente, mucho mas grande y mas ruidosa.

**Palabra clave de UniProt (keyword).** Etiqueta controlada de anotacion
funcional. `KW-0929` es "Antimicrobial". Ausencia de la palabra clave significa
"no anotado", no "no tiene esa propiedad". Toda la seccion de negativos de este
proyecto cuelga de esa distincion.

**FASTA.** Formato de texto para secuencias: una linea de cabecera que empieza
con `>` y luego la secuencia. Simple y universalmente mal implementado.

---

## Machine learning, para lectores de biologia

**Embedding.** Un vector de numeros que representa una secuencia. Lo produce un
modelo entrenado y codifica propiedades aprendidas. La idea clave: secuencias
con funcion parecida acaban cerca en ese espacio, aunque sus letras difieran.

**Modelo de lenguaje de proteinas.** Red neuronal entrenada a predecir
aminoacidos ocultos en millones de secuencias. Para hacerlo bien tiene que
aprender restricciones estructurales y evolutivas, que es lo que se aprovecha.

**Transferencia de aprendizaje.** Usar un modelo entrenado en una tarea grande
(predecir aminoacidos en millones de proteinas) para otra pequena (clasificar
unos miles de peptidos). Es lo que hace viable este proyecto con tan pocos datos.

**Train / validation / test.** Tres conjuntos disjuntos. Entrenamiento ajusta el
modelo; validacion elige hiperparametros y umbral; prueba se mira **una sola vez**
al final. Mirar la prueba varias veces y elegir en consecuencia es, en la
practica, entrenar con ella.

**Fuga de informacion (data leakage).** Cuando informacion del conjunto de prueba
llega al entrenamiento. Aqui la via es la homologia: un pariente cercano en
entrenamiento hace que el acierto en prueba sea memorizacion.

**ROC-AUC.** Probabilidad de que el modelo puntue mas alto un positivo al azar
que un negativo al azar. 0.5 es azar, 1.0 es perfecto.

**PR-AUC.** Area bajo precision-recall. Mas informativa que ROC cuando una clase
es mucho mas abundante que la otra.

**MCC (coeficiente de correlacion de Matthews).** Resumen de la matriz de
confusion entre -1 y 1. Solo sale alto si el modelo acierta en ambas clases, por
lo que es dificil de inflar.

**Matriz de confusion.** Tabla de 2x2: verdaderos positivos, falsos positivos,
falsos negativos, verdaderos negativos.

**Bootstrap.** Remuestrear los datos con reemplazo muchas veces para estimar la
incertidumbre de una metrica. Produce intervalos de confianza sin suponer una
distribucion.

**Sobreajuste (overfitting).** El modelo memoriza el entrenamiento en lugar de
aprender el patron, y falla en datos nuevos.

**Baseline.** Un modelo simple contra el que se compara. Sin baseline, un numero
alto no significa nada.

**Regresion logistica.** Modelo lineal para clasificacion binaria. Simple e
interpretable; aqui es la comparacion justa, porque se usa el mismo modelo sobre
dos representaciones distintas.

**XGBoost.** Conjunto de arboles de decision con refuerzo del gradiente. Suele
ser el mejor en datos tabulares de tamano moderado.

**UMAP.** Metodo para proyectar vectores de muchas dimensiones a dos, para poder
dibujarlos. Util para ver agrupaciones, pero las distancias en el dibujo no son
fiables: dos puntos cercanos en el UMAP no estan necesariamente cerca en el
espacio original.

**Ejemplo adversario.** Una entrada construida para que el modelo se equivoque
con alta confianza. Es lo que produce un algoritmo genetico que optimiza contra un
clasificador, y la razon por la que un clasificador es un filtro y no una funcion
objetivo.
