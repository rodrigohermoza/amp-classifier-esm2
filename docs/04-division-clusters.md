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
mmseqs easy-cluster entrada.fasta prefijo tmp --min-seq-id 0.4 -c 0.5 --cov-mode 1
```

Rapido y es el estandar actual. Se instala con `brew install mmseqs2`
(~30 MB). `-c 0.5 --cov-mode 1` exige que el alineamiento cubra al menos la
mitad de la secuencia mas corta, lo que importa cuando se mezclan peptidos de 5 y
de 50 residuos.

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
   contra representantes que comparten al menos una fraccion minima de 3-meros
   (indice de Jaccard sobre el conjunto de 3-meros). Un alineamiento completo
   contra los ~miles de representantes seria inviable; el prefiltro descarta la
   gran mayoria con operaciones de conjuntos.
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

## Lo que se guarda

`data/processed/splits.csv` con `cluster_id` y `split` por secuencia, y
`data/processed/split_report.json` con: metodo usado, version de la herramienta,
umbral, numero de clusteres, tamanos, porcentajes reales y balance de clases por
split.

El archivo de splits se versiona como artefacto. Asi, aunque MMseqs2 cambie de
version y produzca clusteres distintos, los numeros publicados se pueden
reproducir exactamente.
