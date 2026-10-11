# Problemas encontrados

Bitacora de lo que fallo de verdad. Se mantiene porque es la parte de un proyecto
que normalmente se borra y la que mas ensena a quien viene detras.

Formato: sintoma, causa, resolucion, leccion.

---

## P-01 — La primera linea del FASTA de APD3 no es una secuencia
**Fecha:** 2026-10-11 · **Estado:** resuelto en el parser

**Sintoma:** el archivo `naturalAMPs_APD2024a.fasta` empieza asi:

```
>Your search led to 3306 peptides
>AP00001
GLWSKIKEVGKEAAKAAAKAAGKAALGAVSEAV
```

**Causa:** la descarga se genera desde la interfaz de busqueda del sitio y
arrastra el encabezado del resultado como si fuera un registro FASTA.

**Efecto si no se trata:** un parser ingenuo crea un registro con id
`Your search led to 3306 peptides` y secuencia vacia, o peor, le asigna la
secuencia del primer peptido real, desplazando todas las etiquetas una posicion.
Un desplazamiento asi es silencioso y envenena el conjunto entero.

**Resolucion:** `parse.py` descarta todo registro cuya secuencia quede vacia
despues del filtro de alfabeto, y registra cuantos descarto. Si ese numero no es
exactamente 1 para APD3, el formato del archivo cambio y hay que mirarlo.

**Leccion:** nunca confiar en que un archivo de una base de datos biologica
cumple su formato. Validar y contar siempre; un conteo que no cuadra con el
esperado es la senal mas barata de que algo se rompio.

---

## P-02 — DBAASP inaccesible por script
**Fecha:** 2026-10-11 · **Estado:** sin resolver, con rodeo

**Sintoma:** todas las rutas de `dbaasp.org` fallan:

| Ruta | Resultado |
|---|---|
| `https://dbaasp.org/` | HTTP 403 |
| `https://dbaasp.org/home` | HTTP 403 |
| `https://dbaasp.org/api/v1/peptides?limit=1` | HTTP 200, cuerpo de 0 bytes |
| `https://dbaasp.org/api/v1/peptides/all` | HTTP 200, cuerpo de 0 bytes |
| `https://www.dbaasp.org/api/v1/peptides?limit=1` | HTTP 200, cuerpo de 0 bytes |

Se probo con y sin cabecera `User-Agent` de navegador, y con
`Accept: application/json`. Mismo resultado.

**Causa probable:** proteccion anti-bot (tipo Cloudflare) que exige ejecucion de
JavaScript, o filtrado a nivel de red en el entorno de desarrollo. El patron
"200 con cuerpo vacio" apunta mas a lo primero.

**Efecto:** se pierde la unica fuente con valores de MIC, y con ella la regresion
de potencia.

**Rodeo adoptado:** entregar la clasificacion binaria con APD3 + DRAMP + UniProt,
que no necesita DBAASP. Ver [08-decisiones.md](08-decisiones.md), D-07.

**Vias de recuperacion sin resolver:**
1. Descargar a mano desde un navegador y colocar el archivo en `data/raw/`; el
   parser lo aceptaria. Rompe la reproducibilidad de un solo comando, asi que
   seria un camino opcional documentado, no el principal.
2. Buscar un espejo de DBAASP en Zenodo o en el material suplementario de algun
   articulo.
3. Usar los valores de MIC que DRAMP trae en sus metadatos. Son menos y menos
   estructurados, pero existen y vienen en un archivo al que si se accede.

**Leccion:** verificar el acceso real a cada fuente **antes** de escribir el plan
que depende de ella. Comprobar las cuatro fuentes al principio costo unos minutos
y cambio el alcance del proyecto; descubrirlo en la etapa 5 habria costado dias.

---

## Plantilla para entradas nuevas

```
## P-NN — Titulo corto
**Fecha:** · **Estado:** resuelto | con rodeo | sin resolver

**Sintoma:** que se observo, con el mensaje de error literal si lo hubo.
**Causa:** la causa real, no la primera hipotesis.
**Resolucion:** que se cambio, con el archivo y la funcion.
**Leccion:** que se haria distinto la proxima vez.
```
