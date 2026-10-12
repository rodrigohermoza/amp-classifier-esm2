# Uso de agentes de IA

Este archivo registra como se uso asistencia de IA en el proyecto: que hizo cada
agente, que decidio una persona, y el coste en tokens y tiempo.

Existe por una razon practica: si este repositorio va en un CV y en
postulaciones a maestrias, la pregunta "¿que parte es tuya?" es legitima y
merece una respuesta verificable en vez de una vaga.

## Convencion

- **Agente**: una ejecucion delegada y aislada, con su propio contexto, lanzada
  para una subtarea concreta.
- **En linea**: trabajo hecho en la conversacion principal, sin delegar.

Los tokens son del contexto acumulado de la sesion segun lo reporta la
herramienta de uso; no son una factura. El tiempo es de reloj, incluyendo las
esperas de red.

---

## Sesion 1 — 2026-10-11 · Planificacion y construccion inicial

**Agentes lanzados: ninguno.** Todo el trabajo fue en linea. La razon: cada
agente arranca sin contexto y tiene que volver a derivar lo que la conversacion
principal ya sabe. Para un proyecto donde las decisiones se encadenan (la
eleccion de fuentes condiciona el alcance, que condiciona la arquitectura),
delegar habria costado mas tokens y habria producido decisiones peor informadas.

### Trabajo hecho en linea

| Tarea | Resultado |
|---|---|
| Verificacion de las 4 fuentes de datos | Tamanos reales medidos por HTTP; DBAASP detectado como inaccesible (P-02) |
| Decision de alcance | Regresion de MIC diferida; clasificacion binaria como entrega principal |
| Estructura del repositorio | Arbol de carpetas, `pyproject.toml`, `Makefile`, `.gitignore`, licencia |
| Documentacion | 13 archivos en `docs/` |
| Codigo | Modulos de `src/ampesm/` |

### Coste

| Metrica | Valor |
|---|---|
| Tokens de contexto acumulados | ~107 000 |
| Ventana de contexto | 1 000 000 (11% usado) |
| Modelo | Claude Opus 5 |
| Tiempo de reloj aproximado | ~35 min |
| Llamadas de red a fuentes de datos | 18 (solo cabeceras y sondas, sin descargas) |

### Decisiones que tomo la persona, no la IA

- Alcance del proyecto y encuadre como pieza de portafolio.
- El limite de 150 MB para los datos, y la instruccion de revisar los tamanos
  antes de comprometerse.
- La exigencia de `docs/` como cuaderno de laboratorio y no como README ampliado.
- La regla de escribir todo antes de ejecutar la tuberia completa.
- La aprobacion explicita de cada instalacion y descarga.

### Lo que la IA aporto como juicio propio, no como ejecucion

Vale la pena separarlo porque es donde un lector esceptico mirara:

- Senalar que el problema de los negativos (no-anotado ≠ no-antimicrobiano) es la
  debilidad central del campo, y que debe ir escrita en el README y no escondida.
- Proponer el bootstrap por cluster en vez de por secuencia.
- Proponer la agrupacion conjunta de ambas clases (D-03), que es un detalle facil
  de pasar por alto y que abre una via de fuga.
- Reencuadrar el algoritmo genetico como experimento negativo (D-10).

---

## Sesion 2 — 2026-10-11 · Auditoria por agentes

Primera vez que se delegan subtareas. El motivo: el codigo estaba escrito
completo pero sin ejecutar, y revisar tu propio codigo recien escrito tiene un
punto ciego conocido — se verifica lo que se penso al escribirlo, no lo que se
dio por supuesto.

### Agentes lanzados

| Agente | Funcion | Que hizo | Tokens | Tiempo | Herramientas |
|---|---|---|---|---|---|
| Auditoria metodologica | Validez estadistica: fuga, seleccion de modelo, bootstrap, correccion numerica de los descriptores | Leyo docs/02, 04 y 07 y comparo la intencion declarada contra el codigo. Midio el prefiltro de k-meros con pares sinteticos, ejecuto `assign_splits` con distribuciones adversas, comparo las escalas fisicoquimicas residuo a residuo contra las publicadas | 152 454 | 9 min 51 s | 24 llamadas |
| Auditoria de robustez | Fallos que apareceran en la primera corrida real | Copio `src/` al scratchpad y **ejecuto de verdad** `split → baselines → evaluate → interpret` con 368 peptidos sinteticos. Construyo un `EsmModel` local desde config para probar los embeddings sin red. Midio el clustering con hasta 16 000 secuencias y lo perfilo con cProfile | 161 569 | 12 min 46 s | 38 llamadas |
| **Total** | | | **314 023** | **22 min 37 s** (en paralelo) | 62 |

Ambos corrieron en paralelo y en modo solo lectura: ninguno modifico un archivo
del proyecto. Trabajaron sobre copias en un directorio temporal.

### Que encontraron

12 hallazgos, registrados como P-03 a P-11 en
[09-problemas-encontrados.md](09-problemas-encontrados.md). Los dos criticos:

1. **El prefiltro de k-meros anulaba el umbral de identidad del 40%** (P-03). El
   umbral efectivo era ~0.75. Es la fuga de datos que todo el diseno del proyecto
   existe para evitar. **Los dos agentes lo encontraron de forma independiente**,
   con metodologias distintas y numeros que coinciden, lo que sube mucho la
   confianza en que es real y no un artefacto de un analisis.
2. **El modelo de la comparacion principal se elegia por su resultado en test**
   (P-05), contradiciendo dos compromisos escritos en este mismo repositorio.

Tambien confirmaron explicitamente que varias cosas estan bien, lo cual tiene
valor propio: no hay fuga por escalado, `paired_difference` si usa los mismos
remuestreos para ambos modelos, la agregacion enmascarada de ESM-2 es correcta
(diferencia individual-vs-lote de 4,8e-07), ninguna columna de procedencia entra
como feature, y el contrato de columnas entre etapas se cumple en una corrida
real de punta a punta.

Y una preocupacion resulto infundada (P-11): se temia que el clustering en
Python puro fuera inviable; son **26 segundos** para 16 000 peptidos.

### Coste total de la sesion

| Metrica | Valor |
|---|---|
| Tokens de los agentes | 314 023 |
| Tokens de la conversacion principal | 269 161 (27% de una ventana de 1 M) |
| Modelo | Claude Opus 5 |
| Tiempo de reloj de la sesion | ~2 h |

### Valoracion del uso de agentes

Merecio la pena, y se puede decir por que con precision: los dos hallazgos
criticos eran invisibles desde dentro. P-03 no lo detecta ningun test porque el
test verifica el invariante correcto en el eslabon equivocado; P-05 requiere
cruzar lo que dice la documentacion con lo que hace el codigo, que es
exactamente lo que no haces cuando escribiste los dos.

El patron que funciono: **dar a cada agente una hipotesis que contrastar y
permiso para ejecutar**, no "revisa este codigo". El segundo agente encontro el
fallo de XGBoost (P-09) porque intento correr la etapa de verdad; una lectura del
codigo no lo habria visto nunca, porque el codigo es correcto y lo que falta es
una biblioteca del sistema.

Lo que **no** se delego, y a proposito: las decisiones. Los agentes informan,
proponen correcciones y las ordenan por prioridad; que se arregla, en que orden y
que se documenta como limitacion conocida se decide en la conversacion principal,
con el autor del proyecto.

---

## Sesion 3 — 2026-10-11 · Arreglos y verificacion

### Trabajo en linea

Se aplicaron los 12 hallazgos de la sesion 2 (P-03 a P-13 en
[09-problemas-encontrados.md](09-problemas-encontrados.md)), se instalaron
MMseqs2 y libomp, y se verifico la tuberia completa con 623 peptidos
sinteticos: split con MMseqs2, 9 baselines, embeddings de los dos modelos de
ESM-2, entrenamiento, evaluacion con la comparacion preregistrada e
interpretacion con 8 figuras. 114 tests, ruff limpio.

Dos problemas nuevos aparecieron **durante** los arreglos, no en la auditoria:

- **P-12**, el mas instructivo de todo el proyecto: la auditoria de fuga recien
  escrita dio un 29% que parecia catastrofico y resulto ser ruido de fondo de
  la metrica. Lo delato que la cola alta del histograma estaba vacia. Se
  resolvio anadiendo una distribucion nula.
- **P-13**: los 114 tests fallaron de golpe por un `.pth` del instalable
  editable que dejo de procesarse. Un fallo del entorno disfrazado de fallo del
  codigo.

### Agentes de verificacion

Se lanzaron dos agentes en paralelo para comprobar que los arreglos eran
correctos y suficientes, con el encargo explicito de buscar **problemas nuevos
introducidos por los arreglos**, que es el riesgo real de una tanda grande de
cambios:

| Agente | Funcion |
|---|---|
| Verificacion metodologica | Medir el prefiltro corregido a varios niveles de identidad; comprobar si el corte temprano por cota superior puede descartar pares validos; validar la distribucion nula con fuga plantada; comprobar que test ya no influye en la seleccion; revisar la tabla de Radzicka y Wolfenden residuo a residuo |
| Verificacion de robustez | Provocar las ocho guardas una por una; comprobar que XGBoost degrada sin tumbar la etapa; forzar divisiones degeneradas; corrida completa en una copia; medir el rendimiento del clustering con el prefiltro mas laxo |

**Estado al compactar la sesion: ambos agentes seguian corriendo y sus
informes no llegaron a registrarse aqui.** Si no se revisaron despues, conviene
relanzarlos antes de publicar resultados. Ver
[14-estado-y-continuacion.md](14-estado-y-continuacion.md), punto 6.

### Que se aprendio sobre el uso de agentes

Lo de la sesion 2 se confirmo: el patron util es **dar una hipotesis concreta
que contrastar y permiso para ejecutar**, no "revisa este codigo". Los dos
hallazgos criticos de la sesion 2 (P-03 y P-05) eran invisibles desde dentro,
y los encontraron porque midieron y ejecutaron en vez de leer.

Un matiz nuevo: al encargar la verificacion conviene pedir explicitamente que
digan **que confirman como correcto**, no solo lo que esta mal. Un informe que
solo lista problemas no permite distinguir "lo revise y esta bien" de "no lo
revise".

## Plantilla para sesiones nuevas

```
## Sesion N — FECHA · Titulo

### Agentes lanzados
| Agente | Funcion | Que hizo | Tokens | Tiempo |
|---|---|---|---|---|

### Trabajo en linea
### Coste
### Decisiones de la persona
```

## Politica para este proyecto

1. Se declara el uso de IA. No se oculta y no se exagera.
2. Ningun numero del README sale de una estimacion de un modelo: todos se generan
   corriendo el codigo y quedan en `results/metrics/`.
3. El codigo se lee y se entiende antes de aceptarlo. Un repositorio de
   portafolio que el autor no puede explicar en una entrevista es un pasivo, no
   un activo.
4. Los errores se registran en
   [09-problemas-encontrados.md](09-problemas-encontrados.md) aunque los haya
   causado la IA. Ahi esta la parte educativa.
