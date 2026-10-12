# Fase 2: diseño de la parte de generacion

**Estado: diseñado, aprobado el encuadre, sin escribir una linea de codigo.**
Esta fase no empieza hasta que la clasificacion este corrida y publicada, y la
razon es tecnica, no de gestion: esta en la seccion "Por que no antes".

Este documento existe para que el diseño no se pierda entre sesiones. Lo que
hay aqui es el **plan**, con sus trampas identificadas y las decisiones ya
tomadas. Lo que no hay es implementacion.

---

## Por que esta fase existe

En 2023, en Sup'Biotech, el proyecto del equipo diseñaba AMPs con un algoritmo
genetico ([repositorio del equipo](https://github.com/hgrandjean/Anti-Fungi-Peptide)).
El aporte de Rodrigo fue la obtencion de datos de DRAMP y UniProt. La primera
mitad de este proyecto hace que el **clasificador** sea suyo de principio a
fin; esta fase cierra el arco, porque un clasificador es exactamente la funcion
de aptitud que un algoritmo genetico necesita.

Dicho de otra forma: en 2023 dio los datos para que el generador de otros
puntuara candidatos. Al terminar esta fase, el oraculo, el generador y **la
critica del generador** son suyos.

### La simetria que justifica el diseño

El proyecto actual hace una pregunta: *¿aporta ESM-2 algo sobre medir carga e
hidrofobicidad?* Esta fase hace **la misma pregunta otra vez** en otro
contexto: *¿aporta ESM-2 algo sobre un algoritmo genetico clasico?*

| | Metodo simple | Metodo con modelo de lenguaje |
|---|---|---|
| **Clasificacion** | 31 descriptores fisicoquimicos | embeddings de ESM-2 |
| **Generacion** | algoritmo genetico con mutacion aleatoria | evolucion guiada por ESM-2 |

Un proyecto con una pregunta hecha dos veces se lee como una tesis. Dos mitades
sin relacion se leen como una coleccion de modulos.

---

## La trampa central, que es el resultado

**Un algoritmo genetico que maximice la puntuacion del clasificador va a
converger en algo como `RRRRRRRRRRRRRR`.**

No es una exageracion retorica. La poli-arginina saca puntuacion altisima en
casi cualquier clasificador de AMPs, porque el clasificador aprendio que la
carga positiva predice actividad, **y es verdad que la predice**. Pero la
poli-arginina es conocida por ser inespecifica y toxica para celulas animales:
inutil como farmaco. El clasificador no lo sabe porque **nadie publica
poli-arginina en una base de datos de AMPs**, asi que esa region del espacio de
secuencias nunca tuvo ejemplos negativos. El modelo no es que se equivoque: es
que se le esta preguntando fuera de donde tiene datos.

El generador, por tanto, no encuentra peptidos buenos. Encuentra **los puntos
ciegos del clasificador**. En terminos de ML es un ataque adversarial contra el
propio modelo; fuera de ML es la ley de Goodhart: una metrica usada como
objetivo deja de ser una metrica.

### La consecuencia de encuadre, que es la decision mas importante del diseño

**El resultado honesto de esta fase no es "diseñe antimicrobianos nuevos".** Es:

> Construi un generador contra mi propio clasificador, medi cuanto lo engaña, y
> estas son las guardas que hacen falta para que no lo haga.

Eso es publicable, honesto y no necesita laboratorio. La version ingenua —"mira,
genere 50 peptidos con puntuacion 0.99"— no es defendible ante nadie que
conozca el campo, y es lo que hace la mayoria de los repositorios publicos sobre
el tema.

Esta inversion del encuadre es el nucleo del diseño. Todo lo demas se deriva de
ella: si el objetivo es **medir el engaño**, entonces hacen falta dos oraculos,
hace falta un control positivo, y hace falta que las metricas de exito no sean
la puntuacion del oraculo que se esta optimizando.

---

## Bloque 1 — Endurecer el oraculo antes de atacarlo

No se puede optimizar contra un clasificador y evaluar con el mismo
clasificador: no habria forma de distinguir "mejoro" de "aprendio a engañarlo".

### 1.1 Dos oraculos independientes

- **Oraculo de optimizacion:** ESM-2 35M + logistica. Es el que ve el
  generador.
- **Oraculo de evaluacion:** composicion + fisicoquimicos + XGBoost. Entrenado
  en el mismo `train`, **nunca** visto por el generador.

La **brecha entre las dos puntuaciones es la medida del engaño**, y es el
numero principal de esta fase. Si lo generado saca 0.99 en el suyo y 0.40 en el
independiente, ya esta el resultado.

Que sean de representaciones distintas no es un detalle: dos modelos sobre la
misma representacion comparten los mismos puntos ciegos, asi que la brecha
saldria artificialmente pequeña. Los dos oraculos son justamente el par
preregistrado de la primera mitad del proyecto, lo que ahorra entrenarlos.

### 1.2 Un penalizador de toxicidad

Hallazgo util, ya verificado sobre los datos descargados: el XLSX de DRAMP
tiene una columna `Hemolytic_activity` con **informacion real en 5 818 de sus
12 784 filas**. Ejemplos literales:

```
[Ref.22467870] HD50 = 11.6 μM against human red blood cells (Type A).
Tn-AFP1 has no adverse effect on red blood cells at concentration of 32 μg mL−1
[Ref.12054867] No hemolytic activity was detected when purified aureocin A53...
```

Las otras 6 966 filas dicen literalmente "No hemolysis information or data
found in the reference(s) presented in this entry", asi que son ausencia de
dato y no ausencia de hemolisis: **el mismo problema que los negativos de la
primera mitad**, y hay que tratarlo igual (declararlo, no ignorarlo).

Da para un clasificador binario grueso de hemolisis que sirva de penalizador.
Es la guarda que mata a la poli-arginina, porque la poli-arginina **si** esta
caracterizada como hemolitica en la literatura. Seria aproximado y hay que
decirlo en el README.

**No hace falta descargar nada:** ese archivo ya se baja, ya esta en el
manifiesto con su `sha256` y ya lo parsea `parse_xlsx_source`.

### 1.3 Restricciones duras

Gratis, y evitan la mitad de las degeneraciones antes de que aparezcan:

- longitud entre 5 y 50, como en el resto del proyecto;
- entropia de Shannon minima de la secuencia;
- ningun aminoacido por encima de un porcentaje del total.

El tercero por si solo elimina la poli-arginina. No sustituye al penalizador
aprendido: una restriccion dura dice "esto no", y el penalizador dice "esto
cuanto".

---

## Bloque 2 — Dos generadores

### 2.1 El baseline: algoritmo genetico

Deliberadamente parecido al de 2023: poblacion, seleccion por torneo, mutacion
puntual, cruce, elitismo.

```
aptitud = puntuacion(oraculo_opt) - λ_hemo · penalizacion_hemolisis
                                  - λ_nov · penalizacion_novedad
```

La penalizacion de novedad es importante y se explica en 3.2: sin ella, el
optimo trivial es copiar un AMP del conjunto de entrenamiento.

Esta aqui por dos razones: es la continuidad con el proyecto de 2023, y es la
vara contra la que se mide el metodo con modelo de lenguaje.

### 2.2 El contendiente: evolucion guiada por ESM-2

En vez de mutar al azar, se **enmascara una posicion y se le pregunta a ESM-2
que aminoacido iria ahi**, que es literalmente la tarea con la que se entreno
(*masked language modeling*). Las mutaciones propuestas son entonces plausibles
como proteina, y el oraculo solo selecciona entre propuestas razonables.

Esto ataca la degeneracion por construccion: ESM-2 no propone la decimocuarta
arginina seguida, porque eso no se parece a nada de lo que leyo en UniRef. El
prior del modelo de lenguaje hace de regularizador sobre el espacio de
busqueda.

**No necesita nada nuevo:** el modelo ya esta en cache (160 MB). Lo unico que
cambia es usar su cabeza de lenguaje enmascarado (`EsmForMaskedLM`) en vez de
las capas intermedias que usa `embeddings.py`.

### 2.3 Hipotesis preregistrada

A escribir **antes** de correr, igual que la de la primera mitad:

> El GA ganara en la puntuacion del oraculo de optimizacion y perdera en el
> oraculo independiente. La version guiada por ESM-2 puntuara mas bajo en el
> suyo y aguantara mejor en el otro. Es decir: la brecha entre oraculos sera
> **mayor para el GA**.

Si sale al reves, es un hallazgo y se reporta tal cual. Si el GA aguanta igual
de bien en el oraculo independiente, entonces el prior de ESM-2 no aporta nada
aqui, que es una respuesta legitima a la pregunta del proyecto.

---

## Bloque 3 — Evaluar sin laboratorio

Cuatro medidas. Tres reutilizan codigo que ya existe y esta probado.

### 3.1 Brecha entre oraculos
La medida del engaño. Es el numero principal.

### 3.2 Novedad
Identidad maxima de cada peptido generado contra **todo** el conjunto de
entrenamiento. Un generador que escupe magainina-2 con dos letras cambiadas no
ha generado nada: ha recuperado el optimo trivial de su funcion de aptitud.

**Reutiliza `_max_identity_against` de `split.py` sin tocarla.** Es la misma
funcion de la auditoria de fuga aplicada a otra pregunta, lo cual es una señal
de que el modulo estaba bien partido.

### 3.3 Realismo distribucional
¿Los 11 descriptores fisicoquimicos caen dentro del rango de los AMPs reales o
se van a las colas? ¿La composicion es degenerada? Los descriptores ya estan
escritos y validados contra literatura (D-14).

### 3.4 Control positivo: recuperacion de una familia escondida

Esta es la prueba que de verdad convence, y es **D-15 aplicado a la
generacion**.

1. Se retira **entera** una familia de AMPs reales del entrenamiento. DRAMP
   tiene candidatas naturales: 393 peptidos de la subfamilia de las brevininas
   (piel de rana), 358 de la familia DEFL.
2. Se entrenan los oraculos sin ella.
3. Se corre el generador.
4. Se mide si lo generado se parece a lo escondido, con la misma metrica de
   identidad.

Si el generador **no** redescubre nada parecido a una familia de AMPs reales
que existe y que no vio, no esta diseñando AMPs: esta puntuando alto. Igual que
el detector de fuga, un generador que nunca ha encontrado lo que busca no esta
demostrado.

---

## Por que no antes de cerrar la clasificacion

No es prioridad: es dependencia. **El generador no puede interpretarse sin un
oraculo medido.**

Ahora mismo no se sabe si el clasificador saca 0.90 o 0.70 en test, si los
embeddings ganaron al baseline, ni si la longitud era el atajo que los numeros
de `data` sugieren. Un peptido generado con puntuacion 0.95 de un clasificador
cuyo rendimiento real se desconoce es **un numero sin unidades**. Y la brecha
entre oraculos no se puede leer si no se conoce la brecha que ya existe entre
los dos oraculos sobre datos reales: parte de ella es desacuerdo normal entre
modelos, no engaño del generador. **Ese desacuerdo de base es la linea cero de
esta fase**, y sale de la primera mitad.

---

## Coste estimado

**Descargas: ninguna.** Pesos de ESM-2 ya en cache, DRAMP ya bajado, codigo de
alineamiento y descriptores ya escrito.

**CPU:** cada candidato nuevo necesita su embedding, que es lo mas lento del
proyecto. Una corrida de 200 individuos x 100 generaciones son ~20 000
evaluaciones, del mismo orden que embeber el conjunto entero (18 644).

Mitigaciones previstas:
- modelo de 8M en el bucle interno, 35M solo para la puntuacion final;
- el cache por `seq_sha1` que ya existe absorbe los candidatos repetidos, que
  en un GA con elitismo son muchos.

Estimacion honesta: **10-25 minutos por corrida** con el modelo pequeño,
multiplicado por el numero de semillas. Minimo tres semillas, porque un GA con
una sola semilla no dice nada.

**Codigo nuevo:** un modulo de generacion, uno de penalizadores, uno de
evaluacion de lo generado, y sus tests. Del orden de `split.py`, que es el
modulo mas grande del proyecto. No es poco.

---

## Decisiones ya tomadas

| # | Decision | Estado |
|---|---|---|
| 1 | El encuadre es "mido cuanto engaña mi generador a mi clasificador", no "diseñe antimicrobianos" | **aprobada** |
| 2 | Dos generadores (GA y guiado por ESM-2), para mantener la pregunta del proyecto | **aprobada** |
| 3 | Fase 2: despues de cerrar y publicar la clasificacion | **aprobada** |
| 4 | La regresion de MIC sale del alcance mientras esta fase este viva | **aprobada** |

Sobre la 4: un MIC predicho de mala calidad metido como funcion de aptitud
haria el generador peor y mucho mas dificil de interpretar. Si mas adelante se
quiere potencia, entra mejor como filtro posterior que como objetivo.

## Decision pendiente

**El nombre del repositorio.** `amp-classifier-esm2` deja de describir el
proyecto cuando esta fase entre. Opciones: renombrar (GitHub redirige los
enlaces viejos) o dejarlo y que el README lo explique. No es urgente, pero
conviene decidirlo antes de que el enlace circule en solicitudes.

## Lo primero que habria que hacer al retomar esta fase

1. Releer la seccion de la trampa central. Si el encuadre se pierde, la fase se
   convierte en el repositorio ingenuo del que intenta distinguirse.
2. Medir el **desacuerdo de base entre los dos oraculos** sobre el conjunto de
   test real. Es la linea cero contra la que se lee todo lo demas, y sale casi
   gratis de los resultados de la primera mitad.
3. Escribir el control positivo (3.4) **antes** del generador, no despues. Si
   se escribe despues, se acabara ajustando hasta que el generador lo pase.
