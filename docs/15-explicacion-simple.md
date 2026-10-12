# El proyecto explicado en simple

Este documento cuenta todo el proyecto sin dar por sabido nada: que problema
resuelve, que hace cada pieza, por que esta ahi y como se conectan. Si algo de
aqui te suena a jerga, esta en el [glosario](12-glosario.md). La version con
formulas, nombres de funcion y detalles de implementacion es
[16-explicacion-tecnica.md](16-explicacion-tecnica.md).

---

## 1. El problema, en una frase

Hay cadenas cortas de aminoacidos que matan bacterias. Queremos un programa que,
dada una cadena nueva, diga si probablemente mata bacterias o no.

### Por que eso importa

Un **peptido antimicrobiano** (AMP) es una cadena de entre 5 y 50 aminoacidos
que mata bacterias u hongos. La mayoria no lo hace envenenando una enzima
concreta, como hacen muchos antibioticos clasicos, sino **rompiendo la membrana**
de la bacteria. Esa diferencia es la razon de todo el interes: una bacteria
puede mutar una enzima para dejar de ser sensible a un antibiotico, pero
rediseñar su membrana entera es muchisimo mas dificil. Por eso los AMPs se
estudian como respuesta a la resistencia a antibioticos.

El cuello de botella es practico: probar un peptido en el laboratorio es lento y
caro. Sintetizarlo, cultivar bacterias, medir si mueren. Si tienes diez mil
candidatos no puedes probarlos todos. Un programa que los **ordene** de mas
prometedor a menos prometedor, aunque se equivoque a veces, ahorra meses.

### Como mata un AMP (los 30 segundos de biologia que hacen falta)

La membrana de una bacteria tiene carga negativa en su superficie; la de tus
celulas, mucho menos. Un peptido con **carga positiva** es atraido
electrostaticamente hacia la bacteria y no hacia ti. Una vez pegado, si el
peptido tiene una **cara grasa** (hidrofobica) y otra cara que se lleva bien con
el agua, puede clavarse en la membrana como una cuña y abrirla.

Esas dos propiedades — carga positiva y tener una cara grasa — son casi toda la
historia. Guarda eso, porque es la razon por la que este proyecto es
interesante: **se puede predecir bastante bien con dos numeros**, y eso pone la
vara muy alta para cualquier metodo sofisticado.

---

## 2. La pregunta real del proyecto

No es "¿podemos acertar mucho?". Acertar mucho aqui es facil y no demuestra
nada. La pregunta es:

> ¿Un modelo de lenguaje de proteinas, que ha "leido" millones de secuencias y
> aprendido patrones generales, aporta algo **por encima** de medir carga e
> hidrofobicidad con una calculadora?

Dos candidatos compiten:

- **El baseline:** ~31 numeros calculados con formulas de los años 80 y 90
  (carga, hidrofobicidad, longitud, etc.) y un modelo estadistico encima.
- **El contendiente:** **ESM-2**, una red neuronal entrenada con millones de
  proteinas, que convierte cada peptido en una lista de 480 numeros que resumen
  "lo que la red entiende de esa secuencia".

Y la respuesta puede perfectamente ser "no aporta nada". Eso tambien es un
resultado, y el README lo dira en el titular si pasa.

---

## 3. Por que el proyecto se preocupa tanto por hacer trampa sin querer

Esta es la parte que de verdad distingue este proyecto, asi que va con su
propia seccion.

### La trampa

Los AMPs vienen en **familias**. La magainina 1 y la magainina 2 se diferencian
en dos aminoacidos de veintitres. En las bases de datos hay docenas de variantes
del mismo peptido, muchas de ellas sinteticas.

Para medir si un modelo funciona, se parte el conjunto de datos: con una parte
se entrena y con otra, que el modelo **nunca** vio, se mide. Si esa particion se
hace al azar, la magainina 1 cae en entrenamiento y la magainina 2 en medicion.
El modelo acierta la 2 porque practicamente ya habia visto la 2.

Entonces el numero que publicas mide **memoria**, no capacidad de generalizar. Y
mide alto. Un 0.98 de acierto que no significa nada.

Esto, en mi opinion, es el error mas comun en la literatura de prediccion de
AMPs, y buena parte de los numeros espectaculares que se publican vienen de
aqui.

### La solucion

Antes de partir nada, se **agrupan** todas las secuencias por parecido: todo lo
que se parezca mas de un 40% a otra cosa cae en el mismo grupo. Y despues se
reparten **grupos completos**, nunca secuencias sueltas. Las magaininas entran
las dos en el mismo grupo y el grupo entero va a un solo lado.

Resultado: cuando mides, el modelo se enfrenta a familias que no ha visto nunca.
El numero baja, y es el numero honesto.

### La solucion no basta, y eso costo aprenderlo

Las herramientas que agrupan por parecido son rapidas porque hacen trampa a
proposito: no comparan todo contra todo, comparan cada secuencia contra un
"representante" de cada grupo. Eso significa que dos parientes **pueden**
acabar en grupos distintos y colarse a lados opuestos.

Asi que no basta confiar en la agrupacion: hay que **medirlo**. El proyecto
coge cada secuencia del conjunto de medicion y busca su pariente mas parecido
en todo el conjunto de entrenamiento.

Y aqui esta la parte que mas me enseño. La primera vez que corri esa
comprobacion dijo que el **29%** de las secuencias de medicion tenian un
pariente demasiado parecido en entrenamiento. Leido tal cual, eso significaba
que el metodo entero no servia.

Era falso. Resulta que si coges dos peptidos cortos **sin ninguna relacion** y
los alineas, se parecen un 33% de media, solo por azar: son cadenas cortas de un
alfabeto de 20 letras. El 29% no era fuga de informacion, era **el ruido de
fondo de la propia medida**, y de hecho estaba por debajo del azar.

La forma de arreglarlo es tener algo con que comparar. A cada peptido de
medicion se le **barajan los aminoacidos** 20 veces. Barajar conserva la
longitud y la composicion exactas (que es lo que genera el parecido de fondo) y
destruye cualquier parentesco real. Si un peptido se parece a entrenamiento
**mas que sus propias 20 barajas**, eso ya no es azar. Se cuentan cuantos
peptidos superan esa prueba, y se compara con cuantos la superarian por pura
casualidad.

Costo **tres intentos** llegar a una version que funcionara. Las dos primeras
devolvian numeros creibles y **no detectaban nada**, ni cuando les plantaba a
mano veinte parientes al 80% de parecido. Ninguna fallaba de forma visible: las
dos simplemente decian "todo bien".

De ahi la regla que quedo escrita en el proyecto: **un diagnostico que nunca ha
encontrado el problema que busca no es un diagnostico.** Ahora hay una prueba
automatica que planta fuga a proposito y exige que la herramienta la encuentre.

---

## 4. La debilidad del proyecto, dicha en voz alta

Para entrenar hace falta enseñarle ejemplos de las dos clases: peptidos que
matan bacterias y peptidos que no.

Los que matan son faciles: hay bases de datos de AMPs confirmados en
laboratorio.

**Los que no matan no existen como tales.** Nadie publica una base de datos de
"peptidos que probamos y no hacian nada". Lo que se usa aqui, como en toda la
literatura, son peptidos de UniProt que **no estan anotados** como
antimicrobianos. Y eso no significa "no es antimicrobiano": significa **"nadie
lo ha probado, o nadie lo anoto"**.

Dos consecuencias, y la segunda es la grave:

1. Hay AMPs reales escondidos entre los negativos, asi que el modelo es
   castigado por acertar. El rendimiento medible tiene un techo por debajo del
   real. En la corrida real aparecieron **271 secuencias que figuran a la vez
   como positivas en una base y como negativas en otra**, lo cual pone un suelo
   de ~3% a ese ruido. Es una estimacion por lo bajo.
2. El modelo puede aprender a distinguir **de donde viene el dato** en vez de
   **que hace el peptido**. Los AMPs de las bases especializadas son peptidos
   secretados, cationicos, muchos de piel de rana; los fragmentos de UniProt son
   otra cosa. Un acierto altisimo puede significar "se de que archivo salio
   esto", no "se si mata bacterias".

Lo que se hace al respecto: colorear los graficos por fuente de datos ademas de
por clase, medir el acierto en cada fuente por separado, y emparejar las
longitudes entre clases. Si el modelo acierta mucho peor en los positivos de
UniProt que en los de las bases especializadas, depende de la procedencia y hay
que leer el resultado con eso en mente.

Esto no se arregla con codigo. Se declara.

### Y un confusor concreto que aparecio al correr los datos

El emparejado de longitudes del punto anterior **no se pudo cumplir**. Los
positivos miden 22.9 aminoacidos de media y los negativos 30.6. El programa
intento compensarlo y le faltaron 23 732 negativos cortos que simplemente **no
existen** en UniProt.

Asi que la longitud es un atajo disponible: un modelo podria acertar bastante
diciendo "si es corto, es AMP", y eso no es biologia, es un artefacto de como se
construyo el conjunto.

La decision tomada: entrenar el baseline **dos veces**, una con la longitud y
otra sin ella (ni la longitud ni la carga-por-residuo, que tambien la lleva
dentro). La diferencia entre los dos numeros **es** la medida de cuanto del
baseline era longitud disfrazada. Se reportan los dos. Detalle del razonamiento
en [08-decisiones.md](08-decisiones.md), D-16.

---

## 5. El recorrido completo, etapa por etapa

```
descarga → limpieza → agrupar por parecido → partir → ┬→ calculadora → modelos → medir
                                                      └→ ESM-2       ┘
```

Siete etapas. Cada una deja sus archivos en disco, asi que se pueden correr
sueltas y repetir solo la que cambio.

### Etapa 1: `data` — conseguir los peptidos

Baja cuatro archivos de tres sitios publicos (APD3, DRAMP y UniProt), unos 5 MB
en total, y de cada uno guarda la direccion, la fecha y una **huella digital**
(un `sha256`, un codigo que cambia si el archivo cambia aunque sea un caracter).
Eso es lo que permite que otra persona compruebe que bajo exactamente lo mismo.

El repositorio **no incluye los datos**, solo los scripts para bajarlos: dos de
las bases piden citacion y restringen el uso comercial, y redistribuirlas no
seria limpio.

Luego limpia:

- Fuera todo lo que tenga letras raras (aminoacidos no estandar, que las formulas
  no saben medir).
- Fuera lo mas corto de 5 y lo mas largo de 50 aminoacidos.
- Fuera los duplicados exactos. Y cuando una misma secuencia aparece como
  positiva en un sitio y negativa en otro, **gana la positiva**: alguien la
  probo en un laboratorio, lo cual pesa mas que una ausencia de anotacion.

De 25 557 registros quedaron **18 644 peptidos** (10 594 positivos, 8 050
negativos).

### Etapa 2: `split` — agrupar y partir

Lo del apartado 3: agrupar al 40% de parecido con una herramienta llamada
MMseqs2, repartir grupos completos en 70% entrenamiento / 15% validacion / 15%
medicion, y despues **auditar** que no se colo nadie.

Hay tres cosas que vigilar aqui y el programa las vigila solo: que ningun grupo
aparezca en dos lados a la vez (si eso falla, **ningun numero del proyecto
vale**), que los tres trozos tengan ejemplos de las dos clases, y la auditoria
de parentesco con su baraja de control.

**¿Por que tres trozos y no dos?** Porque hay decisiones que hay que tomar
mirando datos: cuanta regularizacion, donde poner el corte entre "si" y "no".
Si esas decisiones se toman mirando el trozo de medicion, el trozo de medicion
ya dejo de ser limpio. Entonces: se decide mirando **validacion**, y
**medicion** se mira una sola vez, al final del todo.

### Etapa 3: `baselines` — la calculadora

Dos familias de numeros:

- **Composicion** (20 numeros): simplemente que porcentaje de la secuencia es
  cada aminoacido. Tonto de simple, y sorprendentemente fuerte, porque la carga
  y la grasitud estan ahi dentro de forma indirecta.
- **Fisicoquimicos** (11 numeros): carga neta, punto isoelectrico, grasitud
  media, **momento hidrofobico** (lo explico abajo porque es el bonito), fraccion
  apolar, aromaticidad, indice alifatico, indice de Boman, hidrofobicidad
  interfacial, longitud y carga por residuo.

El **momento hidrofobico** es el descriptor mas elegante del conjunto. Muchos
AMPs se enrollan en helice, y en una helice cada aminoacido esta girado 100°
respecto al anterior. Imagina clavar en cada aminoacido una flecha que apunta
hacia afuera en su angulo, y cuya longitud es lo graso que es. Si sumas todas
las flechas y el resultado es grande, significa que los grasos estan todos en
**una cara** de la helice y los demas en la otra. Eso es exactamente la forma de
cuña que necesita para clavarse en la membrana. Un numero que mide una forma.

Estos numeros alimentan tres modelos: regresion logistica (lineal e
interpretable), random forest y XGBoost. Nueve combinaciones, mas tres de la
ablacion sin longitud.

Un detalle de higiene que importa: las formulas y sus tablas de valores estan
**escritas a mano dentro del repositorio** y con su cita, en vez de llamar a una
libreria. Dos razones: se pueden auditar sin salir del proyecto, y los numeros
no cambian solos cuando una libreria se actualiza.

Y una razon mas, aprendida a golpes: una de esas tablas estaba **invertida**. El
indice de Boman de un AMP conocido daba 0.35 en vez de 3.0, un orden de
magnitud, y el aminoacido que deberia dominar la escala quedaba en medio de la
tabla. No fallaba nada: solo empeoraba el modelo en silencio. Desde entonces
**cada descriptor entra con una prueba que compara su valor contra el publicado
para un peptido conocido**.

### Etapa 4: `embeddings` — ESM-2

ESM-2 es un modelo de lenguaje, pariente lejano de los que escriben texto, pero
entrenado con secuencias de proteinas en vez de frases. Se entreno tapando
aminoacidos al azar y pidiendole adivinar cual falta, millones de veces. Para
hacerlo bien tuvo que aprender, sin que nadie se lo dijera, cosas sobre
estructura y funcion.

Aqui **no se reentrena nada**: se usa como una funcion fija que convierte un
peptido en una lista de numeros (320 o 480, segun el tamaño del modelo). Se usan
los dos modelos mas pequeños a proposito: todo tiene que correr en un portatil
sin tarjeta grafica, y el modelo de 35 millones de parametros pesa 140 MB.
Los resultados se guardan en disco, porque calcularlos otra vez es lo mas lento
de todo el proyecto.

### Etapa 5: `train` — modelos sobre ESM-2

Regresion logistica y una red neuronal pequeña encima de esos numeros. Los
ajustes se eligen mirando validacion.

### Etapa 6: `evaluate` — medir, con honestidad

Se miden varias cosas (ROC-AUC, PR-AUC, MCC, F1...) porque cada una esconde algo
distinto, pero lo importante es **como** se miden.

**Con intervalo de confianza.** Un numero suelto como "0.84" no dice si es
0.84 solido o 0.84 de suerte. La tecnica, llamada bootstrap, es: vuelve a medir
2 000 veces sobre muestras ligeramente distintas de los mismos datos y mira
cuanto baila el resultado. Si baila entre 0.79 y 0.89, el numero es 0.84 con esa
horquilla.

**Y el bootstrap remuestrea grupos, no secuencias.** Es el mismo razonamiento
del apartado 3: las secuencias emparentadas no son observaciones independientes,
y tratarlas como si lo fueran da horquillas artificialmente estrechas, que te
hacen creer que distingues dos modelos que en realidad no distingues.

**Dos modelos cuyas horquillas se solapan no se distinguen.** Esa frase es la
que decide la conclusion del proyecto, no quien saca el numero mas alto.

**Y la pareja a comparar estaba elegida de antemano.** Escrito en el README
antes de correr nada: ESM-2 de 35M con logistica contra composicion mas
fisicoquimicos con XGBoost. ¿Por que fijarlo antes? Porque si entrenas trece
modelos y luego comparas el mejor de un grupo contra el mejor del otro, estas
comparando dos ganadores de concurso, y los ganadores de concurso ganan en parte
por suerte. La diferencia sale inflada **siempre**, incluso si los dos metodos
son identicos. La tabla con los trece modelos se publica como material
descriptivo; **la unica conclusion del proyecto es la comparacion de esa pareja
fijada**.

Tambien esta escrita la expectativa: que ESM-2 gane, pero poco (0.02-0.04), y
puede que con horquillas solapadas. Si gana por mucho, lo primero que hay que
sospechar no es exito, es fuga de informacion.

### Etapa 7: `interpret` — mirar los errores

Un mapa en dos dimensiones (UMAP) de todos los peptidos, coloreado por clase,
por fuente y por trozo, que deja ver si las clases se separan o si lo que se
separan son las bases de datos. Y una lista de en que se equivoca el modelo:
¿falla mas en los cortos? ¿en los de una fuente concreta? Un error que se
concentra en un sitio suele ser una pista; un error repartido, el techo real
del metodo.

---

## 6. Como esta construido por dentro

- **Todo lo ajustable vive en un solo archivo** (`params.yaml`): los umbrales,
  los tamaños, los modelos. Asi un cambio de configuracion aparece en el
  historial de git como lo que es, y no escondido en medio del codigo.
- **118 pruebas automaticas.** No son adorno: una de ellas comprueba que ningun
  grupo este en dos trozos a la vez, que es la pieza de la que depende la
  validez de todo lo demas. Otras comparan cada formula contra un valor
  publicado. Otra planta fuga de informacion a proposito y exige que la
  auditoria la encuentre.
- **Un cuaderno de laboratorio** (esta carpeta `docs/`), con las decisiones y su
  justificacion, y con los 24 problemas reales que aparecieron y como se
  resolvieron. Incluidos los que me costo tres intentos y los que me cause yo
  mismo arreglando otra cosa.

## 7. En que punto esta

Bajados y limpiados los datos reales: **18 644 peptidos**. El codigo esta
completo y auditado, pero las etapas de la 2 a la 7 todavia no han corrido con
datos reales, asi que las tablas de resultados estan **vacias a proposito**. Se
rellenaran con lo que salga, no con lo que convenga.

Lo que falta, lo que esta decidido y lo que queda por decidir esta en
[14-estado-y-continuacion.md](14-estado-y-continuacion.md).
