# Division sin fuga de informacion

Esta es la etapa que decide si los numeros del proyecto significan algo. Si se
hace mal, todo lo demas es decoracion.

## El problema

Las secuencias biologicas no son observaciones independientes. Los AMPs llegan en
familias: la magainina 1 y la magainina 2 difieren en dos residuos; hay decenas de
variantes de cecropina publicadas; APD3 contiene muchos analogos sinteticos de un
mismo peptido natural.

Una division aleatoria 80/20 pone casi con certeza miembros de la misma familia en
entrenamiento y en prueba. El modelo puede entonces acertar en prueba por
reconocer una secuencia que ya vio con dos letras cambiadas. La metrica mide
memorizacion, no generalizacion, y queda inflada.

Este es, en mi opinion, el error metodologico mas comun en los articulos de
prediccion de AMPs. Buena parte de los AUC de 0.98 que se reportan en la
literatura vienen de aqui.

## La solucion

1. Agrupar **todas** las secuencias (positivas y negativas juntas) en clusteres
   con un umbral de identidad del 40%.
2. Asignar **clusteres completos**, no secuencias, a `train` / `val` / `test`.

Asi, ninguna secuencia de prueba tiene un pariente cercano en entrenamiento. El
modelo tiene que generalizar a familias que no vio.

## Por que 40%

Es el umbral convencional en este campo y marca aproximadamente el limite por
debajo del cual dos proteinas ya no se consideran homologas de forma obvia. Mas
alto (70–90%) solo elimina duplicados casi exactos y deja fuga por homologia.
Mas bajo (25–30%) es territorio de "zona de penumbra", donde el alineamiento
local es poco fiable para secuencias tan cortas.

Hay que decir que 40% de identidad sobre un peptido de 15 residuos es una medida
ruidosa: son 6 residuos coincidentes. Para secuencias cortas, cualquier umbral de
identidad es aproximado. Es mejor que dividir al azar, no es perfecto.

## Herramientas, en orden de preferencia

El codigo intenta tres caminos y usa el primero disponible. El metodo realmente
usado se registra en `data/processed/split_report.json`, porque afecta a la
reproducibilidad.

### 1. MMseqs2 (preferido)

```
mmseqs easy-cluster entrada.fasta prefijo tmp --min-seq-id 0.4 -c 0.5 --cov-mode 1 \
    -s 7.5 --alignment-mode 3 --cluster-mode 1
```

Rapido y es el estandar actual. Se instala con `brew install mmseqs2`
(~30 MB). `-c 0.5 --cov-mode 1` exige que el alineamiento cubra al menos la
mitad de la secuencia mas corta, lo que importa cuando se mezclan peptidos de 5 y
de 50 residuos.

**`--cluster-mode 1` es obligatorio y no es el valor por defecto** (P-29). Los
tres modos dan garantias distintas, y solo una sirve aqui:

| Modo | Garantiza | Sirve para dividir? |
|---|---|---|
| 2, incremental voraz (**el defecto**) | cada miembro se parece a **su** representante | **No** |
| 0, set-cover | idem, con otra heuristica de eleccion | No |
| 1, componente conexa | si dos secuencias se parecen, estan en el mismo grupo, por transitividad | **Si** |

La diferencia no es academica. Con el modo por defecto, dos miembros de
clusteres distintos pueden parecerse entre si todo lo que quieran mientras
ninguno sea representante del otro, y eso es exactamente fuga. Medido sobre los
datos reales: `DLPSGWSSYEGH` quedo en train y `DCPSDWSPYEGHCYK` en test con
**75% de identidad** y cobertura 1.00. MMseqs2 veia el parecido perfectamente
(valor E 3.5e-08 al preguntarselo directamente); el modo de agrupamiento no lo
usaba.

Coste de cambiar al modo 1 sobre los 18 644 peptidos: 8 173 -> 7 886 clusteres,
el mayor pasa de 379 a 454 secuencias (2.4% del total). No aparecen cadenas
degeneradas, asi que la division 70/15/15 sigue siendo posible.

### 2. CD-HIT

```
cd-hit -i entrada.fasta -o salida -c 0.4 -n 2 -l 4
```

Respaldo. Ojo: CD-HIT con `-c 0.4` necesita `-n 2` (palabras de 2 caracteres) y
avisa de que por debajo de 0.4 su heuristica pierde fiabilidad. Para peptidos
cortos hay que bajar `-l` (longitud minima) del valor por defecto de 10.

### 3. Fallback propio en Python (sin dependencias externas)

Clustering incremental voraz al estilo CD-HIT, implementado en
`ampesm/split.py`. Existe para que `make split` funcione en una maquina limpia
sin instalar nada.

Algoritmo:

1. Ordenar las secuencias de mas larga a mas corta.
2. Para cada secuencia, compararla con los representantes de los clusteres ya
   creados. Si la identidad con alguno supera el umbral, se une a ese cluster.
   Si no, se vuelve representante de un cluster nuevo.
3. **Prefiltro por k-meros** para evitar el coste cuadratico: solo se alinea
   contra representantes que comparten suficientes 3-meros. Un alineamiento
   completo contra los miles de representantes seria inviable; el prefiltro
   descarta la mayoria con operaciones de conjuntos.

   El criterio es la **contencion**, `|A inter B| / min(|A|,|B|)`, con umbral
   0.05. Esto no es un detalle de implementacion: la primera version usaba el
   indice de Jaccard con umbral 0.15 y eso **anulaba el umbral de identidad**.
   A 45% de identidad real, entre el 0% y el 1% de los pares llegaba a
   alinearse, asi que el umbral efectivo era ~0.75 en vez de 0.40. Jaccard
   castiga la diferencia de tamano, de modo que un peptido de 6 residuos
   contenido **entero** dentro de uno de 50 daba 0.06 y se descartaba, pese a
   tener identidad y cobertura de 1.00. Historia completa en
   [09-problemas-encontrados.md](09-problemas-encontrados.md), P-03.
4. La identidad se calcula con alineamiento global de Biopython
   (`PairwiseAligner`, matriz BLOSUM62, modo global) y se define como
   `residuos identicos / longitud de la secuencia mas corta`.

Es mas lento que MMseqs2 y da clusteres ligeramente distintos. Para ~16 000
peptidos cortos corre en pocos minutos en un Mac con Apple Silicon.

## Proporciones y como se reparten

Objetivo: **70% train / 15% val / 15% test**, contado en secuencias, no en
clusteres (los clusteres tienen tamanos muy desiguales).

El reparto no es trivial: hay que llenar cuotas de secuencias asignando bloques
indivisibles de tamano variable, manteniendo a la vez el balance de clases. El
procedimiento usado:

1. Ordenar los clusteres por tamano, de mayor a menor.
2. Recorrerlos y asignar cada uno al split que esta mas lejos (en terminos
   relativos) de cumplir su cuota.
3. Dentro de ese criterio, preferir el split cuyo balance de clases se degrada
   menos.

Es una heuristica voraz. Se usa porque el reparto exacto es un problema de
particion (NP-duro) y porque la aproximacion es suficiente: lo que importa es que
no haya fuga, no que los porcentajes sean exactos. El reporte guarda los
porcentajes reales obtenidos.

Los clusteres grandes van primero a proposito: si se dejan para el final, uno
solo puede desbordar la cuota de `test`.

## Invariante que se verifica

`tests/test_split.py` comprueba que la interseccion de `cluster_id` entre
cualquier par de splits sea vacia. **Si ese test falla, ninguna metrica del
proyecto es valida.** Esta puesto como test y no como comentario a proposito.

Pero ese invariante **no basta**, y conviene entender por que. Verifica la
*asignacion*, no la *agrupacion*: si dos parientes acaban en clusteres
distintos, cada uno puede irse a un split distinto sin violar ninguna regla. Es
exactamente lo que pasaba con el prefiltro viejo. Un test puede ser correcto y
aun asi vigilar el eslabon equivocado de la cadena.

## La verificacion que de verdad respalda la premisa

Ninguna herramienta de clustering voraz garantiza que no haya parientes entre
train y test. Todas, incluidas MMseqs2 y CD-HIT, comparan cada secuencia solo
contra **representantes**, asi que A y B al 50% de identidad pueden acabar
separadas si cada una se unio a un representante distinto. No es un fallo, es
una propiedad del algoritmo.

Asi que en vez de confiar en la agrupacion, se mide el resultado: para cada
secuencia de test, su **identidad maxima contra todo el conjunto de
entrenamiento**. El resultado va a `split_report.json` con su histograma.

### Y ese numero se lee contra una distribucion nula

La primera corrida de esta auditoria dijo que el 29% de las secuencias de test
tenian identidad >= 0.4 contra train, lo que leido tal cual significaria que
todo el metodo no sirve. Es falso: **la identidad de secuencia tiene un ruido
de fondo alto en peptidos cortos**. Con alineamiento global y huecos
permitidos, dos peptidos sin ninguna relacion alinean al 0.33 de media.

El control es **por secuencia**: cada peptido de test se baraja 20 veces, lo
que conserva su longitud y su composicion exactas (que es lo que determina su
identidad de fondo) y destruye cualquier parentesco, y se mide la identidad
maxima de cada baraja contra train. Una secuencia se **marca** cuando su
identidad observada supera todas sus propias barajas.

Bajo la hipotesis de que no hay fuga, cada secuencia se marca con probabilidad
1/(R+1), asi que el numero de marcadas sigue una binomial y un test binomial da
un p-valor global. Lo que se reporta como fuga es ese p-valor y el exceso de
marcadas sobre lo esperado por azar:

```json
"leakage": {
  "n_flagged": 6,
  "n_expected_by_chance": 5.7,
  "excess": 0.3,
  "global_p_value": 0.51,
  "significant": false
}
```

**La auditoria corre dos veces: test contra train y `val` contra train.** Lo
natural es auditar solo test, que es el que da el numero publicado, pero `val`
fija los hiperparametros y el umbral de decision, y ese umbral se aplica despues
a test: entra en cinco de las siete metricas. Si `val` tiene parientes en train,
el umbral se eligio sobre un conjunto optimista y las metricas de test salen
sesgadas **aunque test este limpio**. La fuga no necesita tocar test para
contaminar test (P-27).

Validado plantando fuga conocida: detecta 20 parientes de 120 desde el 50% de
identidad, y hay un caso que la planta **solo en val** para comprobar que la
auditoria de val la ve y la de test, correctamente, no. Costo tres intentos llegar a un estadistico que funcionara, y los
dos primeros devolvian numeros plausibles sin detectar nada; la historia
completa esta en [09-problemas-encontrados.md](09-problemas-encontrados.md),
P-12, P-15 y P-18. La leccion que dejo: **un diagnostico que nunca se ha
probado contra un caso positivo conocido no es un diagnostico.**

### Honestidad sobre el respaldo en Python

Medido de extremo a extremo (pares que de verdad cumplen identidad >= 0.4 y
cobertura >= 0.5 y acaban en el mismo cluster):

| Identidad real | Recall del respaldo |
|---|---|
| 0.40 | 47% |
| 0.45 | 53% |
| 0.60 | 80% |
| 0.80 | 85% |

Es decir: **a 40-45% los pares se agrupan la mitad de las veces**. Unos 15
puntos los pierde el prefiltro y el resto es la voracidad por representante,
que es intrinseca al algoritmo. Por eso el respaldo marca `publishable: false`
y MMseqs2 es obligatorio para numeros publicables (D-11).
