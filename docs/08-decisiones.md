# Registro de decisiones

Cada entrada: que se decidio, que alternativas habia, por que esta, y que haria
cambiar de opinion. El formato busca que una decision se pueda revisar dentro de
seis meses sin tener que reconstruir el razonamiento.

---

## D-01 — No versionar los datos crudos
**Fecha:** 2026-10-11

APD3 y DRAMP piden citacion y restringen el uso comercial; ninguno declara una
licencia tipo Creative Commons. Subir sus archivos a un repositorio publico MIT
seria, en el mejor caso, discutible.

**Decision:** el repositorio incluye los scripts de descarga y
`data/raw/manifest.json` con URL, fecha, tamano y `sha256` de cada archivo. Los
datos se reconstruyen con `make data`.

**Alternativas:** subir los datos (riesgo de licencia); usar solo UniProt, que es
CC BY 4.0 (perderia los positivos mejor curados).

**Que lo cambiaria:** que APD3 o DRAMP publicaran una licencia explicitamente
permisiva.

---

## D-02 — Clustering al 40% de identidad y division por cluster
**Fecha:** 2026-10-11

Es la decision central del proyecto. Justificacion completa en
[04-division-clusters.md](04-division-clusters.md).

**Decision:** agrupar positivos y negativos juntos al 40% de identidad y asignar
clusteres completos a train/val/test.

**Alternativas:** division aleatoria (infla metricas por homologia); umbral mas
alto, 70–90% (solo quita casi-duplicados); dividir por organismo de origen (mas
realista para algunos usos, pero reduce mucho el conjunto).

**Que lo cambiaria:** nada en el principio. El umbral concreto es discutible y se
podria reportar una curva de rendimiento frente al umbral, que seria un buen
anexo.

---

## D-03 — Agrupar las dos clases juntas, no por separado
**Fecha:** 2026-10-11

Si se agrupan positivos y negativos en ejecuciones separadas, un positivo de test
puede tener un pariente cercano entre los negativos de entrenamiento. Dado que el
modelo ve ambas clases en entrenamiento, eso tambien es fuga.

**Decision:** una sola ejecucion de clustering sobre el conjunto completo.

**Efecto secundario esperado:** apareceran clusteres con ambas clases dentro.
Eso es informacion valiosa, no un problema: son los casos donde secuencias
parecidas tienen etiquetas distintas, y probablemente donde esta el ruido de
etiquetado. El reporte cuenta cuantos clusteres son mixtos.

---

## D-04 — Python 3.12, no 3.14
**Fecha:** 2026-10-11

El sistema tiene Python 3.14.6. PyTorch y las ruedas precompiladas del ecosistema
cientifico suelen tardar meses en cubrir una version nueva de Python.

**Decision:** fijar 3.12 en `.python-version`. `uv` descarga una compilacion
independiente de CPython 3.12, sin tocar el Python del sistema.

**Que lo cambiaria:** que PyTorch publique ruedas estables para 3.14 en arm64.

---

## D-05 — Empezar por ESM-2 8M y 35M, no por 650M
**Fecha:** 2026-10-11

Restriccion dura del proyecto: tiene que correr en CPU en un Mac. El de 650M pesa
~2,5 GB y es lento en CPU.

**Decision:** 8M (~31 MB) y 35M (~140 MB). La comparacion entre los dos responde
por si sola si la escala ayuda en esta tarea.

**Que lo cambiaria:** acceso a GPU gratis en Kaggle o Colab; entonces 150M y 650M
son extensiones directas, porque el codigo solo cambia en el nombre del modelo.

---

## D-06 — Clustering propio en Python como respaldo
**Fecha:** 2026-10-11

MMseqs2 no esta instalado y requiere Homebrew; CD-HIT, lo mismo. Si `make split`
falla en una maquina limpia, el proyecto deja de ser reproducible para quien lo
clone.

**Decision:** implementar un clustering voraz al estilo CD-HIT con prefiltro de
k-meros, usando solo Biopython. Se usa MMseqs2 si esta disponible; si no, el
respaldo. El metodo usado se registra en el reporte.

**Coste:** mas lento y da clusteres algo distintos. Aceptable: el archivo de
splits se versiona, asi que los numeros publicados son reproducibles sin
depender del metodo.

---

## D-07 — La regresion de MIC queda diferida
**Fecha:** 2026-10-11

DBAASP era la unica fuente con valores de MIC y esta inaccesible por script:
403 o respuestas de 0 bytes en todas las rutas probadas. Detalle en
[03-datos-fuentes.md](03-datos-fuentes.md) y [09-problemas-encontrados.md](09-problemas-encontrados.md), P-02.

**Decision:** entregar primero la clasificacion binaria, que no depende de
DBAASP. La regresion se retoma si se consigue acceso.

**Nota adicional:** aunque se consiguiera, la regresion de MIC tiene un problema
propio. Los valores vienen de distintos laboratorios, organismos diana, medios de
cultivo y unidades (µg/mL frente a µM). La variabilidad entre laboratorios para
el mismo peptido puede ser de un factor de 2 a 4. Habria que restringir a un
organismo (p. ej. *E. coli* ATCC 25922), convertir todo a log10 µM, y aun asi el
techo de R² estaria limitado por el ruido experimental, no por el modelo.

---

## D-08 — Bootstrap por cluster, no por secuencia
**Fecha:** 2026-10-11

Las secuencias de un cluster son parientes y no son independientes. Remuestrear
secuencias individuales produce intervalos de confianza demasiado estrechos.

**Decision:** el bootstrap remuestrea clusteres completos.

**Coste:** intervalos mas anchos, o sea numeros menos impresionantes. Es el
precio de que sean correctos.

---

## D-09 — El umbral se elige en validacion
**Fecha:** 2026-10-11

**Decision:** maximizar MCC en `val` y aplicar ese umbral tal cual a `test`.
Test se evalua una sola vez.

Esta escrito como decision explicita porque es facil de violar sin mala
intencion: basta mirar la matriz de confusion de test, pensar "con 0.4 saldria
mejor", y cambiarlo.

---

## D-10 — El algoritmo genetico, si se hace, es un experimento negativo
**Fecha:** 2026-10-11

Usar el clasificador como funcion de aptitud de un algoritmo genetico
reproduciria el proyecto de Sup'Biotech. Pero optimizar contra un clasificador
produce **ejemplos adversarios**: secuencias con puntuacion altisima y propiedades
absurdas (por ejemplo `KKKKKRRRRRKKKKK`), porque el optimizador explota los
artefactos del modelo en vez de descubrir biologia.

**Decision:** si se implementa, el entregable es el analisis de ese fallo y las
mitigaciones (penalizar la distancia a la distribucion de los datos, restringir la
composicion, usar un comite de modelos), no una lista de candidatos.

Es mas honesto y, para una entrevista, mas interesante: demuestra que entiendo el
limite de la herramienta. No se publicaran secuencias propuestas como candidatos
de sintesis.
