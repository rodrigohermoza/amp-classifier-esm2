# Objetivos

## El problema

Un peptido antimicrobiano (AMP) es una cadena corta de aminoacidos, tipicamente
de 10 a 50 residuos, que mata o inhibe bacterias, hongos o virus. Son de interes
porque muchos actuan rompiendo la membrana celular en vez de inhibir una enzima
concreta, lo que hace mas difícil que el microbio desarrolle resistencia.

El problema practico: probar un peptido en el laboratorio cuesta tiempo y dinero,
y el espacio de secuencias posibles es absurdamente grande (20^30 para un peptido
de 30 residuos). Un clasificador que ordene candidatos antes de sintetizarlos
tiene valor real.

## Pregunta concreta de este proyecto

> Dada una secuencia de aminoacidos, ¿es un peptido antimicrobiano?

Es clasificacion binaria sobre secuencias. La variante de regresion (predecir la
concentracion minima inhibitoria, MIC) queda diferida: ver
[08-decisiones.md](08-decisiones.md), decision D-07.

## Que cuenta como exito

Deliberadamente, **el exito no es un numero alto de ROC-AUC**. En esta tarea es
facil llegar a 0.95+ y que el modelo no haya aprendido nada sobre actividad
antimicrobiana. El exito es:

1. **Una division de datos que no filtra informacion.** Los AMPs vienen en
   familias de secuencias casi identicas. Si una division al azar pone
   `GLWSKIKEVGKEAAKAAAKAAGKAALGAVSEAV` en entrenamiento y una variante con dos
   mutaciones en prueba, la metrica mide memorizacion, no generalizacion. La
   division por clusteres de identidad es la pieza central del proyecto.
2. **Baselines honestos.** Si una regresion logistica sobre 20 numeros
   (composicion de aminoacidos) iguala a ESM-2, eso es el resultado, y hay que
   publicarlo asi. Un proyecto que solo reporta el modelo grande no demuestra
   criterio.
3. **Incertidumbre reportada.** Metricas con intervalos de confianza por
   bootstrap. Una diferencia de 0.01 en AUC con intervalos que se solapan no es
   una diferencia.
4. **Las debilidades dichas en voz alta.** En particular el problema de los
   negativos (abajo). Un repositorio que nombra su propio punto flaco es mas
   creible que uno que presume el numero.

## La debilidad conocida: los negativos

No existe una base de datos de "peptidos que se probaron y no son
antimicrobianos". Lo que se usa en la literatura, y lo que se usa aqui, son
peptidos de UniProt sin anotacion antimicrobiana. Eso **no** significa
"no antimicrobiano": significa "nadie lo ha probado o nadie lo anoto".

Consecuencias que hay que vigilar:

- Hay falsos negativos reales en el conjunto de negativos. El techo de
  rendimiento medible esta por debajo del rendimiento verdadero.
- Peor: el modelo puede aprender a distinguir **procedencia** en vez de
  **actividad**. Los AMPs curados en APD3 son peptidos secretados, cationicos,
  muchos de piel de anfibio; los fragmentos de UniProt son otra cosa. Un AUC alto
  puede significar "se cual base de datos es" y no "se si mata bacterias".

Mitigaciones aplicadas: emparejar la distribucion de longitud entre clases,
colorear el UMAP por fuente y no solo por clase (si los grupos se separan por
base de datos, hay un problema), y dejar el punto escrito en el README.

## Fuera de alcance

- Predecir estructura tridimensional.
- Predecir la MIC (regresion de potencia). Diferida, ver D-07. La via viable
  ya no es recuperar DBAASP sino parsear los MIC que DRAMP trae en 5 020 de
  sus filas (P-02).
- **Producir una lista de candidatos para sintetizar.** Esto merece precision,
  porque el proyecto va a tener una fase de generacion y la frontera es justo
  aqui: se generan peptidos, pero el resultado declarado es **cuanto engaña el
  generador al clasificador**, no una lista de candidatos. Un clasificador es
  un mal criterio de optimizacion por si solo, y demostrarlo con numeros es el
  objetivo. Diseño completo en
  [17-diseno-generacion.md](17-diseno-generacion.md) y encuadre en D-17.
- Predecir hemolisis **como objetivo propio**. Entra en la fase 2 solo como
  penalizador aproximado del generador, construido con la columna
  `Hemolytic_activity` de DRAMP (informacion real en 5 818 filas), y con la
  misma advertencia que los negativos: las filas sin dato son ausencia de
  informacion, no ausencia de hemolisis.

## Contexto personal

En 2023, en Sup'Biotech (Francia), participe en un proyecto de diseno de AMPs con
un algoritmo genetico ([repositorio del
equipo](https://github.com/hgrandjean/Anti-Fungi-Peptide)). Mi aporte fue la
obtencion de datos de DRAMP y UniProt. Este proyecto es mi propia version, donde
la parte de aprendizaje automatico es trabajo mio de principio a fin.

Con la fase 2 el arco se cierra: en 2023 di los datos para que el generador de
otros puntuara candidatos; al terminar, el oraculo, el generador y **la critica
del generador** son mios.
