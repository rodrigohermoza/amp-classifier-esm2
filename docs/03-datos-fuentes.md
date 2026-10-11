# Fuentes de datos

**Fecha de verificacion de URLs y tamanos: 2026-10-11.**

## Politica de redistribucion

Este repositorio **no contiene los datos crudos**. Contiene los scripts que los
descargan, la URL exacta, la fecha y el `sha256` de cada archivo en
`data/raw/manifest.json`. Quien clone el repositorio corre `make data` y
reconstruye el conjunto.

Dos razones: APD3 y DRAMP piden citacion y restringen el uso comercial, asi que
re-subir sus archivos a GitHub seria discutible; y un repositorio sin blobs de
datos se clona en segundos.

## Tabla resumen

| Fuente | Rol | Archivo | Tamano medido | Registros |
|---|---|---|---|---|
| APD3 | Positivos | `naturalAMPs_APD2024a.fasta` | 150 144 B (147 KB) | 3 306 |
| DRAMP 3.0 | Positivos | `general_amps.xlsx` | 2 960 502 B (2,8 MB) | por confirmar al parsear |
| UniProt (KW-0929) | Positivos | FASTA via API REST | ~0,2 MB | 1 585 |
| UniProt (sin KW-0929) | Negativos | FASTA via API REST | ~1,5 MB | 11 982 |
| **Total** | | | **~5 MB** | |

Los conteos de UniProt vienen de la cabecera `x-total-results` de la API,
consultada el 2026-10-11. Pueden variar en releases posteriores; el manifiesto
guarda el conteo real de cada corrida.

---

## 1. APD3 — Antimicrobial Peptide Database

- **URL base**: https://aps.unmc.edu/
- **Descarga**: https://aps.unmc.edu/assets/sequences/naturalAMPs_APD2024a.fasta
- **Formato**: FASTA. Cabeceras del tipo `>AP00001`. La primera linea del archivo
  es `>Your search led to 3306 peptides`, un artefacto de la interfaz web que hay
  que descartar al parsear (ver [09-problemas-encontrados.md](09-problemas-encontrados.md), P-01).
- **Rol**: positivos. Es la base mas curada de las tres: entradas con actividad
  experimental verificada.
- **Licencia**: uso academico libre con citacion. No se declara una licencia
  tipo Creative Commons en el sitio. **Por eso no se redistribuye.**
- **Cita**: Wang G., Li X., Wang Z. *APD3: the antimicrobial peptide database as
  a tool for research and education.* Nucleic Acids Research, 2016.
- **Subconjuntos adicionales disponibles** (no usados por defecto, utiles para
  analisis por origen): `animalAMPs_APD2024a.fasta` (113 KB),
  `plantAMPs_APD2024.fasta` (13 KB), `bacteriaAMPs_APD2024.fasta` (18 KB),
  `humanAMPs_APD2024.fasta`, `amphibianAMPs_APD2024.fasta`,
  `insectAMPs_APD2024.fasta`.

## 2. DRAMP 3.0 — Data Repository of Antimicrobial Peptides

- **URL base**: https://dramp.cpu-bioinfor.org/
- **Descarga**: `https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.xlsx`
  (la URL con `http://` redirige a `https://` con un 301; seguir redirecciones).
- **Formato**: XLSX. Mas rico en metadatos que APD3 (organismo diana, actividad,
  referencias), lo que lo hace util para el analisis de errores.
- **Rol**: positivos.
- **Licencia**: uso academico y no comercial. No se redistribuye.
- **Cita**: Shi G. et al. *DRAMP 3.0: an enhanced comprehensive data repository
  of antimicrobial peptides.* Nucleic Acids Research, 2022.
- **Nota**: hay solapamiento sustancial con APD3. El deduplicado por secuencia
  exacta en la etapa de limpieza lo resuelve; el conteo de cuantos se pierden se
  registra en el reporte de limpieza.

## 3. UniProtKB — positivos anotados

- **Endpoint**: `https://rest.uniprot.org/uniprotkb/search`
- **Consulta**: `(keyword:KW-0929) AND (length:[5 TO 50])`
- **Formato**: FASTA via `format=fasta`, paginado con `cursor`.
- **Rol**: positivos adicionales, y mas importante, **control de sesgo de fuente**:
  si el modelo rinde bien en los AMPs de APD3 pero mal en los de UniProt, esta
  aprendiendo procedencia.
- `KW-0929` es la palabra clave "Antimicrobial" de UniProt.
- **Licencia**: CC BY 4.0. Es la unica de las tres que se podria redistribuir sin
  problema; se trata igual que las demas por consistencia.

## 4. UniProtKB — negativos

- **Consulta**: `(length:[5 TO 50]) AND (reviewed:true) NOT (keyword:KW-0929) NOT (cc_function:antimicrobial)`
- **Rol**: negativos.
- **Decisiones de la consulta, y por que**:
  - `reviewed:true` (solo Swiss-Prot): la alternativa, incluir TrEMBL, da cientos
    de miles de secuencias en su mayoria predichas automaticamente y sin curar.
    Mas datos pero mucho mas ruido.
  - Se excluye `KW-0929` **y** la mencion de "antimicrobial" en el campo de
    funcion. Lo segundo captura entradas donde la actividad esta descrita en
    texto libre pero nadie asigno la palabra clave. Es un filtro imperfecto.
  - La longitud 5–50 empareja el rango de los positivos. Sin esto, el modelo
    aprende a usar la longitud como atajo.
- **Advertencia principal del proyecto**: estos no son negativos verificados, son
  no-anotados. Ver [01-objetivos.md](01-objetivos.md), seccion "La debilidad
  conocida".

---

## DBAASP: por que no esta

DBAASP (https://dbaasp.org) era la cuarta fuente planeada y la unica con valores
de MIC, es decir la unica que habilitaba la regresion de potencia.

**Estado al 2026-10-11: inaccesible por script.** Todas las rutas probadas
(`/api/v1/peptides`, `/api/v1/peptides/all`, `/home`, la raiz) devuelven o un
403 o un 200 con cuerpo de 0 bytes. El patron sugiere proteccion anti-bot tipo
Cloudflare que exige ejecucion de JavaScript, o un filtrado a nivel de red.

Consecuencia: **la regresion de MIC queda diferida**, no cancelada. La
clasificacion binaria no se afecta, porque los positivos de APD3 + DRAMP +
UniProt son suficientes. Ver [08-decisiones.md](08-decisiones.md), D-07, y
[10-mejoras-futuras.md](10-mejoras-futuras.md) para las vias de recuperacion.

## El manifiesto

Cada corrida de `make data` escribe `data/raw/manifest.json`:

```json
{
  "generated_at": "2026-10-11T12:00:00Z",
  "entries": [
    {
      "source": "apd3",
      "url": "https://aps.unmc.edu/assets/sequences/naturalAMPs_APD2024a.fasta",
      "filename": "naturalAMPs_APD2024a.fasta",
      "bytes": 150144,
      "sha256": "...",
      "downloaded_at": "2026-10-11T12:00:03Z"
    }
  ]
}
```

Este archivo **si** se versiona. Es lo que permite a alguien verificar que
descargo el mismo conjunto que se uso para los numeros publicados, y detectar
que una fuente cambio bajo nuestros pies.
