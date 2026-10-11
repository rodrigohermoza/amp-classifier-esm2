# Indice de la documentacion

Esta carpeta es el cuaderno de laboratorio del proyecto. La idea es que alguien
que llegue sin contexto pueda leer en orden y entender **que** se hizo, **por que**
se hizo asi, y **que salio mal** en el camino.

| Archivo | Que contiene |
|---|---|
| [01-objetivos.md](01-objetivos.md) | Que problema se resuelve, que cuenta como exito, que queda fuera de alcance |
| [02-arquitectura.md](02-arquitectura.md) | Como se conecta todo: modulos, flujo de datos, artefactos de cada etapa |
| [03-datos-fuentes.md](03-datos-fuentes.md) | Cada fuente: URL, fecha de descarga, licencia, tamano, conteos |
| [04-division-clusters.md](04-division-clusters.md) | Por que no se puede dividir al azar y como se agrupa por similitud |
| [05-features-baselines.md](05-features-baselines.md) | Composicion de aminoacidos y descriptores fisicoquimicos |
| [06-embeddings-esm2.md](06-embeddings-esm2.md) | Que es ESM-2, que modelo se usa y como se agregan los embeddings |
| [07-evaluacion.md](07-evaluacion.md) | Metricas, intervalos de confianza por bootstrap, analisis de errores |
| [08-decisiones.md](08-decisiones.md) | Registro de decisiones tecnicas con su justificacion y fecha |
| [09-problemas-encontrados.md](09-problemas-encontrados.md) | Errores reales, que los causo y como se resolvieron |
| [10-mejoras-futuras.md](10-mejoras-futuras.md) | Lo que falta, ordenado por relacion esfuerzo/beneficio |
| [11-agentes.md](11-agentes.md) | Uso de agentes de IA: funcion, trabajo hecho, tokens y tiempo |
| [12-glosario.md](12-glosario.md) | Terminos de biologia para lectores de ML y viceversa |

## Como leer esto segun quien eres

- **Reclutador o comite de admision**: el [README](../README.md) y luego
  [07-evaluacion.md](07-evaluacion.md).
- **Alguien que quiere reproducir el trabajo**: [02-arquitectura.md](02-arquitectura.md)
  y [03-datos-fuentes.md](03-datos-fuentes.md).
- **Alguien que quiere criticar el metodo**: [04-division-clusters.md](04-division-clusters.md)
  y [09-problemas-encontrados.md](09-problemas-encontrados.md), que es donde esta
  lo interesante.
