# Baselines: composicion y fisicoquimica

Los baselines no son un tramite. Son la vara con la que se mide si ESM-2 aporta
algo, y en esta tarea tienen fama de ser dificiles de superar, porque la actividad
antimicrobiana depende mucho de dos propiedades globales y poco de la posicion
exacta de cada residuo.

## Por que estos descriptores

El mecanismo dominante de los AMPs es la disrupcion de membrana. Para eso un
peptido necesita basicamente dos cosas:

1. **Carga positiva neta**, para ser atraido por la membrana bacteriana, que es
   mas aniónica que la de las celulas animales (fosfatidilglicerol,
   cardiolipina, acidos teicoicos).
2. **Una cara hidrofobica**, para insertarse en el interior lipidico.

Casi todos los descriptores de abajo son formas de medir una de esas dos cosas.
Esto es lo que hace que un modelo lineal sobre 25 numeros sea competitivo: los
numeros codifican el mecanismo.

## Conjunto A: composicion de aminoacidos (20 features)

Frecuencia relativa de cada uno de los 20 aminoacidos estandar:
`count(aa) / length`.

Simple hasta lo tonto, y aun asi fuerte: captura indirectamente la carga (K, R),
la hidrofobicidad (L, I, V, F, W) y rasgos estructurales (alto contenido de P o
de G es caracteristico de ciertas familias). No captura nada del orden.

## Conjunto B: descriptores fisicoquimicos (11 features)

| Descriptor | Que mide | Como se calcula |
|---|---|---|
| `length` | Longitud | Numero de residuos |
| `net_charge` | Carga neta a pH 7.4 | Ecuacion de Henderson-Hasselbalch sobre los pKa de las cadenas laterales ionizables y los extremos |
| `isoelectric_point` | pH de carga nula | Busqueda binaria sobre la curva de carga |
| `gravy` | Hidrofobicidad promedio | Media de la escala de Kyte-Doolittle |
| `hydrophobic_moment` | Anfipaticidad | Momento de Eisenberg con angulo de 100° (helice alfa) |
| `hydrophobic_fraction` | Proporcion apolar | Fraccion de A, C, F, I, L, M, V, W, Y |
| `aromaticity` | Contenido aromatico | Fraccion de F, W, Y |
| `aliphatic_index` | Volumen alifatico | Indice de Ikai |
| `boman_index` | Potencial de union a proteinas | Media de las energias de transferencia ciclohexano-agua (Radzicka y Wolfenden) |
| `interfacial_hydrophobicity` | Preferencia por la interfaz lipidica | Media de la escala interfacial de Wimley-White |
| `charge_density` | Carga por residuo | `net_charge / length` |

El **indice de Boman** merece una advertencia, porque este proyecto ya se
equivoco con el una vez. Se define sobre las energias libres de transferencia
**ciclohexano -> agua** de Radzicka y Wolfenden, donde la arginina domina la
escala con 14.92 kcal/mol. La primera version uso por error la escala
interfacial de Wimley-White con el signo invertido: LL-37 daba 0.35 en vez de
3.0, un orden de magnitud, y la arginina quedaba al nivel de la asparagina,
justo al reves de lo que hace util al indice para AMPs cationicos. Detalle en
[09-problemas-encontrados.md](09-problemas-encontrados.md), P-08. La escala de
Wimley-White se conserva, pero como descriptor propio y con su nombre correcto.

Desde entonces, **todo descriptor entra con un test que compara su valor contra
el publicado** para al menos un peptido conocido (D-14).

El **momento hidrofobico** merece una nota, porque es el descriptor conceptualmente
mas interesante. Trata cada residuo como un vector: direccion dada por su angulo
en la helice (100° de giro por residuo) y magnitud dada por su hidrofobicidad.
La suma vectorial es grande cuando los residuos hidrofobicos se agrupan en una
cara de la helice y los hidrofilicos en la otra, es decir cuando el peptido es
anfipatico. Esa es precisamente la arquitectura de la mayoria de los AMPs
helicoidales.

Se calcula sobre la ventana de 11 residuos que maximiza el momento, no sobre toda
la secuencia, porque un peptido puede tener un segmento anfipatico dentro de una
secuencia que globalmente no lo es.

## Conjunto C: A + B (31 features)

La combinacion. Es el baseline principal contra el que se compara ESM-2.

## Modelos

| Modelo | Por que esta |
|---|---|
| Regresion logistica | Lineal, interpretable, con escalado de features. El coeficiente de `net_charge` deberia salir positivo y grande; si no, hay un error en la tuberia. Funciona como verificacion de cordura. |
| Random Forest | Captura interacciones y no lineales sin ajuste fino. |
| XGBoost | Suele ser el mejor en datos tabulares de este tamano. |

Los hiperparametros se eligen por busqueda en rejilla contra `val`, nunca contra
`test`. El espacio de busqueda esta en `params.yaml` para que quede a la vista.

## Interpretacion de coeficientes: una advertencia

La tentacion es leer los coeficientes de la regresion logistica como
"importancia biologica". No se puede, por dos razones:

1. Los descriptores estan correlacionados entre si (GRAVY y la fraccion
   hidrofobica miden casi lo mismo). La regularizacion reparte el peso entre
   features correlacionadas de forma arbitraria.
2. Un coeficiente alto puede reflejar el sesgo de procedencia de los datos y no
   biologia. Si los negativos de UniProt son sistematicamente mas largos, el
   coeficiente de `length` sera grande y no significara nada sobre actividad
   antimicrobiana.

Por eso el analisis de errores de [07-evaluacion.md](07-evaluacion.md) pesa mas
en las conclusiones que la importancia de features.
