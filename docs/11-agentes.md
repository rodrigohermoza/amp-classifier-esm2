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
