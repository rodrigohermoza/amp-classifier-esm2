# Texto para el portafolio

Borrador para `rodrigohermoza.pages.dev`, siguiendo
`pagina-web/docs/11-guia-historias.md`.

Por ser proyecto personal, va en `organizacion: personal` (sección aparte al
final de `/proyectos`) y usa `historiaCorta` de 3 pasos: Antes → Lo que hice →
Ahora, con el título de cada paso como gancho.

**Estado (2026-10-17):** la tubería corrió y los huecos estan rellenos con los
numeros reales. Listo para publicar.

---

## resumen (simple)

Hay bacterias que ya no responden a los antibióticos de siempre. Una vía de
respuesta son unos fragmentos muy cortos de proteína que matan microbios
rompiéndoles la pared, en lugar de bloquearles una pieza interna; al microbio le
cuesta más acostumbrarse a eso. El problema es que probar cada fragmento en el
laboratorio toma semanas, y las combinaciones posibles son más que todos los
granos de arena de la Tierra. Hice un programa que lee la receta de un fragmento
y estima si va a funcionar, para que alguien pruebe primero los que más
prometen.

## resumenTecnico

Clasificador binario AMP vs no-AMP sobre embeddings de ESM-2 (8M y 35M, media
por residuo), comparado contra baselines de composición de aminoácidos y
descriptores fisicoquímicos (carga neta por Henderson-Hasselbalch, momento
hidrofóbico de Eisenberg, GRAVY, pI, índice de Boman) con regresión logística,
random forest y XGBoost. Positivos de APD3, DRAMP 3.0 y UniProt (KW-0929);
negativos de UniProt revisado sin anotación antimicrobiana, emparejados por
longitud. División train/val/test por clústeres de identidad al 40% (MMseqs2, con
respaldo propio en Python) para evitar fuga por homología. Evaluación con
ROC-AUC, PR-AUC, MCC y F1 con IC del 95% por bootstrap de clústeres, y
comparación emparejada entre modelos. Python, PyTorch, scikit-learn, uv.

---

## Paso 1 — «Probar uno por uno, y el reloj corriendo»

**simple.** En un proyecto en Francia, en 2023, el equipo quería encontrar
fragmentos de proteína capaces de frenar un hongo. Cada candidato había que
encargarlo, esperarlo y probarlo en el laboratorio: semanas por intento. Y los
candidatos posibles no son miles ni millones, son una cifra que no tiene nombre
corto. Mi parte entonces fue reunir los datos de los fragmentos que ya se
conocían. Me quedó la pregunta de si se podía adivinar antes, en la computadora,
cuáles valía la pena probar.

**tecnico.** En 2023, en Sup'Biotech, participé en un proyecto que diseñaba
péptidos antifúngicos con un algoritmo genético. Mi aporte fue la obtención y
curado de datos desde DRAMP y UniProt; el algoritmo era del equipo. El cuello de
botella era la función de aptitud: sin un predictor de actividad, la búsqueda en
un espacio de 20^30 secuencias se guía por heurísticas fisicoquímicas sueltas.
Quedó pendiente construir el predictor.

## Paso 2 — «Lo fácil era hacer trampa sin darme cuenta»

**simple.** Rehice el proyecto por mi cuenta, con la parte de aprendizaje
automático hecha por mí. Lo difícil no fue el programa: fue darme cuenta de que
era muy fácil engañarme solo. Estos fragmentos vienen en familias que se parecen
muchísimo entre sí, como hermanos. Si pongo a un hermano en los ejemplos de
estudio y al otro en el examen, el programa aprueba por haberlo visto antes, no
por haber entendido. Tuve que agrupar primero todas las familias y repartirlas
enteras, cada familia completa a estudio o a examen, nunca partida. Con eso las
notas bajan, pero son de verdad. También usé un modelo de lenguaje, parecido a
los que completan frases, pero entrenado con millones de proteínas en vez de
texto, y lo comparé contra un método simple que solo mide dos cosas: si el
fragmento tiene carga positiva y si tiene un lado aceitoso. Ese método simple es
más difícil de superar de lo que parece.

**tecnico.** Reimplementé la tarea como clasificación binaria con el ML propio.
La decisión central fue la división: clustering de todas las secuencias (ambas
clases juntas) al 40% de identidad y asignación de clústeres completos a
train/val/test, para eliminar la fuga por homología. Agrupar las clases por
separado seguiría filtrando, porque el modelo ve ambas en entrenamiento.
Representaciones comparadas: embeddings de ESM-2 frente a 30 features tabulares.
El bootstrap de los intervalos de confianza remuestrea clústeres, no secuencias,
porque las secuencias emparentadas no son independientes. La limitación
principal, escrita en el README, es que los negativos son péptidos **no
anotados** y no negativos verificados: el riesgo real es que el modelo aprenda
procedencia de base de datos en vez de actividad, y por eso el UMAP y el análisis
de errores se desglosan por fuente.

## Paso 3 — «Un filtro, no un oráculo»

**simple.** El resultado honesto es que **el modelo grande no gana**: acierta
igual que medir carga e hidrofobicidad con una tabla de toda la vida, y la
diferencia entre los dos es tan pequena que no se distingue del azar. Eso no es
un fracaso del proyecto, es su resultado: ahora hay un numero donde antes habia
una suposicion.
El programa no decide: ordena una lista para que quien tenga el laboratorio
empiece por arriba. Y sirve igual cuando se equivoca: los fragmentos que marca
como buenos sin estar en las listas conocidas son justo los que nadie ha probado
todavía, y esos son candidatos interesantes por sí mismos. Todo está publicado
con las cuentas a la vista, incluidos los errores que cometí por el camino, para
que cualquiera pueda repetirlo o contradecirme.

**tecnico.** En prueba, ESM-2 35M + logistica da MCC 0,645 [0,543-0,727] y
ROC-AUC 0,895 [0,838-0,940]; el baseline preregistrado (composicion +
fisicoquimicos con XGBoost) da MCC 0,637 [0,540-0,716] y ROC-AUC 0,909
[0,865-0,943]. La **comparacion emparejada** del par preregistrado da una
diferencia de MCC de **+0,008 con IC 95% [-0,033, +0,051]**: incluye el cero, los
modelos no se distinguen. Los intervalos son por bootstrap de **clusteres
completos**, no de secuencias. La division lleva fuga por homologia **medida y
declarada** (exceso +54,2 en test y +65,2 en val sobre 23,81 esperadas por azar),
escrita dentro de cada JSON de metricas y en el pie de cada figura. El repositorio es reproducible con un comando, con
dependencias de versión fija, datos reconstruibles desde la fuente (no
redistribuidos, por licencia) y el mapeo secuencia→clúster→split versionado como
artefacto. La documentación incluye registro de decisiones, bitácora de errores
reales y mejoras pendientes ordenadas por esfuerzo/beneficio. Extensión natural:
usar el clasificador como función de aptitud del algoritmo genético original,
pero el entregable ahí sería mostrar que optimizar contra un clasificador produce
ejemplos adversarios, es decir que un clasificador es un filtro y no una función
objetivo.

---

## cierre (opcional)

**simple.** Lo siguiente sería predecir no solo si funciona, sino cuánto; para
eso hace falta una base de datos que ahora mismo no me deja entrar.

**tecnico.** La regresión de MIC queda diferida por el bloqueo de acceso a
DBAASP (P-02). Vías alternativas: los valores de MIC de los metadatos de DRAMP,
o un espejo en Zenodo.

---

## Notas para la animación

La guía pide animación nueva, de 3 escenas sincronizadas con los pasos, pixel art
128×72, tramado Bayer, paleta del sitio. Nunca reciclada. Idea, a revisar:

1. **Escena 1:** una rejilla enorme de candidatos; el cursor prueba uno, espera
   (barra lenta), falla. Se repite dos veces. Transmite el coste por intento.
2. **Escena 2:** la idea clave, que es visual y por eso vale la pena dibujarla.
   Los candidatos se agrupan en racimos de parecidos; primero se ve una línea de
   reparto que corta un racimo por la mitad (el error) y luego los racimos
   enteros caen a un lado o al otro, sin partirse. Es literalmente la figura del
   método.
3. **Escena 3:** la rejilla se reordena de mejor a peor y se ilumina la franja de
   arriba: no es un veredicto, es una cola de prioridad.

Si se quiere una cuarta escena, la honesta sería el caso adversario: una
secuencia absurda tipo `KKKKKKKK` con puntuación máxima, para ilustrar que el
clasificador es un filtro y no un oráculo. Encaja con el tono de "sin vender
humo" de la guía.

## Checklist de la guía

- [x] `resumen` y `resumenTecnico`
- [x] simple + técnico en cada paso
- [x] Títulos con gancho que se entienden solos
- [x] Empieza en un problema humano (el equipo esperando semanas por intento)
- [x] Tono sin soberbia; se dice que el algoritmo genético era del equipo y que
      mi aporte de 2023 fueron los datos
- [x] Nada confidencial; el repositorio del equipo es académico y público
- [x] Ningún dato inventado → los `[PENDIENTE]` se rellenaron con la corrida del 2026-10-17
- [ ] Animación nueva de 3 escenas y miniatura nueva
- [ ] `npm run build` sin errores
- [ ] Actualizar `05-animaciones.md`, `06-decisiones.md` y `07-pendientes.md`
