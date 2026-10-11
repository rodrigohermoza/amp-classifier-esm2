# Embeddings de ESM-2

## Que es ESM-2

ESM-2 (Evolutionary Scale Modeling, version 2) es un modelo de lenguaje de
proteinas de Meta AI. La idea es la misma que en los modelos de lenguaje de
texto: se entrena un transformer con enmascaramiento de tokens sobre millones de
secuencias de UniRef, pidiendole que prediga el aminoacido oculto a partir del
contexto.

Para hacerlo bien, el modelo tiene que aprender implicitamente que posiciones
estan en contacto en la estructura plegada, que sustituciones son conservativas,
y que patrones son caracteristicos de cada familia. Eso queda codificado en los
estados ocultos, y esos estados son los **embeddings** que se usan aqui como
features.

La ventaja para este proyecto es concreta: el conjunto de datos tiene unos pocos
miles de peptidos, insuficiente para entrenar un modelo de secuencias desde cero.
ESM-2 trae conocimiento de millones de proteinas y se usa como extractor de
features congelado. Es transferencia de aprendizaje.

## Modelos considerados

| Modelo | Parametros | Capas | Dim. embedding | Peso en disco |
|---|---|---|---|---|
| `facebook/esm2_t6_8M_UR50D` | 8 M | 6 | 320 | ~31 MB |
| `facebook/esm2_t12_35M_UR50D` | 35 M | 12 | 480 | ~140 MB |
| `facebook/esm2_t30_150M_UR50D` | 150 M | 30 | 640 | ~600 MB |
| `facebook/esm2_t33_650M_UR50D` | 650 M | 33 | 1280 | ~2,5 GB |

**Se empieza por el de 8M y luego el de 35M.** Razones:

- Ambos corren en CPU en un Mac con Apple Silicon en minutos, no horas.
- La comparacion 8M vs 35M responde una pregunta util por si misma: ¿la escala
  del modelo ayuda en esta tarea? Si 35M no mejora a 8M, no hay motivo para
  pensar que 650M lo haria, y eso es un resultado reportable.
- Los peptidos son cortos (5–50 residuos). Los modelos grandes brillan sobre todo
  capturando contactos de largo alcance en proteinas grandes, que aqui no existen.

Los de 150M y 650M quedan como extension opcional, para Kaggle o Colab si hace
falta GPU.

## Agregacion: de residuos a vector

ESM-2 devuelve un vector por residuo: una matriz de `(longitud, dim)`. El
clasificador necesita un vector de tamano fijo por peptido. Opciones:

1. **Media sobre residuos** (la elegida). Se promedia a lo largo del eje de
   longitud, **excluyendo los tokens especiales** `<cls>` y `<eos>` y cualquier
   relleno. Es el estandar en la literatura de ESM y es robusto.
2. Token `<cls>`. Funciona en BERT porque se entrena con un objetivo a nivel de
   secuencia; ESM-2 se entrena solo con enmascaramiento, asi que `<cls>` no tiene
   un significado global garantizado.
3. Concatenar media + maximo. Duplica la dimension; se prueba como variante.

**Detalle que es facil de arruinar:** si se promedia sin enmascarar el relleno,
los vectores de los peptidos cortos quedan diluidos por los ceros del relleno, y
el modelo aprende la longitud en lugar de la secuencia. La implementacion
multiplica por la mascara de atencion antes de promediar y divide por la suma de
la mascara, no por la longitud del tensor. Hay un test que lo verifica comparando
el embedding de una secuencia procesada sola contra la misma dentro de un lote
con relleno: deben coincidir.

## De que capa se extrae

De la **ultima capa oculta** por defecto. Existe evidencia de que las capas
intermedias funcionan mejor para algunas tareas (las ultimas capas se
especializan en la tarea de enmascaramiento), asi que la capa es un parametro en
`params.yaml` y se prueba al menos una intermedia.

## Cache

Calcular los embeddings es lo mas caro de la tuberia. Se cachean en
`results/embeddings/{modelo}_{capa}_{agregacion}.npz`, indexados por el `sha1` de
la secuencia. La consecuencia practica: reentrenar un clasificador con otros
hiperparametros toma segundos, porque los embeddings ya estan en disco.

La clave es el hash de la secuencia y no el `id`, para que una secuencia presente
en dos fuentes se calcule una sola vez.

El cache **no se versiona** (son decenas de MB). Esta en `.gitignore`.

## Clasificadores sobre embeddings

- **Regresion logistica** con escalado. Es la comparacion justa contra el
  baseline lineal: cambia la representacion, no el modelo.
- **MLP pequeno** (una o dos capas ocultas, dropout). Se mantiene pequeno a
  proposito: con unos pocos miles de ejemplos y 320–480 dimensiones, una red
  grande se sobreajusta sin aportar.

No se hace *fine-tuning* de ESM-2. Con este tamano de datos el riesgo de
sobreajuste es alto y el coste en CPU seria prohibitivo. Queda anotado en
[10-mejoras-futuras.md](10-mejoras-futuras.md) como la extension con mas
potencial si hubiera GPU gratis.

## La comparacion que importa

Todo esto existe para responder: **¿un modelo de lenguaje de proteinas supera a
30 numeros fisicoquimicos en esta tarea, con una division honesta?**

Mi expectativa antes de ver los resultados, anotada aqui para no reescribir la
historia despues: los embeddings deberian ganar, pero por un margen modesto
(quizas 0.02–0.04 de ROC-AUC), y puede que los intervalos de confianza se
solapen. La actividad antimicrobiana depende de propiedades globales que los
descriptores ya capturan bien. Si el margen sale enorme, lo primero que hay que
sospechar es fuga de informacion, no exito.
