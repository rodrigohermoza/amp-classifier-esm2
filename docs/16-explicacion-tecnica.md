# El proyecto explicado en tecnico

Recorrido completo con nombres de modulo, formulas, artefactos y decisiones de
implementacion. La version sin jerga es
[15-explicacion-simple.md](15-explicacion-simple.md). Este documento se solapa a
proposito con [02-arquitectura.md](02-arquitectura.md) (que es la referencia
seca de modulos y artefactos) y con [04](04-division-clusters.md),
[05](05-features-baselines.md), [06](06-embeddings-esm2.md) y
[07](07-evaluacion.md), que entran mas a fondo en su etapa. Aqui esta el hilo
completo de principio a fin, con el **por que** de cada eleccion pegado al
**que**.

---

## 0. La tarea, planteada formalmente

Clasificacion binaria sobre secuencias. Entrada: una cadena
`s ∈ Σ*` con `Σ` = los 20 aminoacidos estandar, `5 ≤ |s| ≤ 50`. Salida: una
puntuacion en `[0,1]` interpretable como `P(AMP | s)`.

La pregunta cientifica no es el valor absoluto de la metrica, sino una
comparacion de **representaciones** con todo lo demas fijo:

```
          ┌─ φ_pc(s)  ∈ R^31   descriptores tabulares
   s ──── ┤
          └─ φ_esm(s) ∈ R^480  ESM-2 congelado, pooling por media
```

y la hipotesis es `MCC(φ_esm) > MCC(φ_pc)` en test, con intervalos por
bootstrap de clusteres. Esta registrada antes de correr nada, con su par fijado
(D-12) y su tamaño de efecto esperado (0.02-0.04 de ROC-AUC).

Dos propiedades del dominio condicionan todo el diseño:

1. **Las observaciones no son i.i.d.** Las secuencias estan relacionadas por
   homologia, en familias de identidad >90%. Esto rompe tanto la particion
   (seccion 2) como los intervalos de confianza (seccion 6).
2. **Los negativos no son negativos verificados**, son no-anotados. El problema
   es formalmente *PU learning* (positive-unlabeled) tratado como clasificacion
   binaria, que es lo que hace la literatura. No se corrige; se declara y se
   acota empiricamente (271 conflictos de etiqueta observados, ~3% suelo).

---

## 1. Etapa `data` — adquisicion y limpieza

**Modulos:** `data/sources.py`, `data/download.py`, `data/parse.py`,
`data/clean.py`. **Artefactos:** `data/raw/*`, `data/raw/manifest.json`,
`data/processed/peptides.csv`, `data/processed/clean_report.json`.

### Fuentes

| Clave | Fuente | Formato | Rol | Parseados |
|---|---|---|---|---|
| `apd3` | APD3 | FASTA | positivos | 3 306 |
| `dramp` | DRAMP 3.0 general | XLSX (28 columnas) | positivos | 12 784 |
| `uniprot_amp` | UniProt KW-0929, `length:[5 TO 50]` | FASTA | positivos | 1 585 |
| `uniprot_neg` | UniProt revisados sin KW-0929 | FASTA | negativos | 11 982 |

Total de red: **5.01 MB**. Las consultas a UniProt estan construidas en
`sources.py` como cadenas explicitas (`QUERY_UNIPROT_POS`, `QUERY_UNIPROT_NEG`)
contra `rest.uniprot.org/uniprotkb/stream`, no copiadas de la interfaz web, para
que la consulta sea parte del codigo versionado.

`download.py` escribe un **manifiesto** con URL, fecha, tamaño y `sha256` por
archivo. El repositorio no versiona los datos (APD3 y DRAMP piden citacion y
restringen uso comercial); el manifiesto es lo que hace verificable que otra
persona bajo el mismo bytestream. `_declared_total()` compara el conteo que
declara la API de UniProt con el numero de registros recibidos, para detectar
una descarga truncada, que de otro modo pasaria por un conjunto mas pequeño.

### Parseo

`iter_fasta` es deliberadamente estricta y devuelve `(header, sequence)` sin
interpretar. Dos trampas reales:

- **APD3 incluye una cabecera espuria** sin secuencia: literalmente
  `"Your search led to 3306 peptides"`. Si se parsea con un lector permisivo
  entra como un registro con secuencia vacia. El codigo la descarta y **avisa si
  el numero de descartadas no es exactamente 1**, porque un cambio de formato en
  origen debe aparecer como aviso y no como silencio (P-01). Comprobado en la
  corrida real: `dropped_empty = 1`, cabecera exactamente esa.
- **Los nombres de columna del XLSX de DRAMP cambian entre versiones.**
  `parse_xlsx_source` los detecta por patron (`SEQ_COL_PATTERNS`,
  `ID_COL_PATTERNS`) y registra en el reporte que columna eligio, en vez de
  fijar un nombre literal que se rompa en la siguiente version. En la corrida
  real: `sequence_column = "Sequence"`, `id_column = "DRAMP_ID"`, 28 columnas.

### Limpieza

`clean_records` aplica, en orden:

1. **Alfabeto.** Fuera cualquier secuencia con un caracter de
   `NON_STANDARD = BJOUXZ*-.`. No es purismo: las escalas fisicoquimicas no
   tienen valor para esos simbolos, y una `X` propagaria un NaN hasta la matriz
   de features.
2. **Longitud.** `5 ≤ |s| ≤ 50` (`params.data.min_length/max_length`).
3. **Deduplicado por secuencia exacta** (`deduplicate`). Dos reglas que son
   decisiones, no detalles:
   - **Conflicto de etiqueta → gana el positivo.** Una secuencia anotada como
     AMP en una base y no anotada en UniProt es, con alta probabilidad, un AMP
     que UniProt no anoto. La asimetria de evidencia es real: un positivo
     implica un experimento, un negativo implica una ausencia.
   - El numero de conflictos se reporta (`label_conflicts`) porque **es una
     estimacion por lo bajo de la tasa de falsos negativos** entre los
     negativos. Observado: **271**, sobre 8 050 negativos finales ≈ 3.4%.
   - Entre duplicados de la misma clase, `SOURCE_PRIORITY` decide que registro
     se conserva, de modo que el `source` resultante sea determinista.
4. **Desambiguacion de ids repetidos** añadiendo `seq_sha1[:8]`, porque un id
   duplicado haria que `splits.csv` y las predicciones no se pudieran cruzar.
5. **Emparejado de longitudes** (`match_lengths`): submuestrea negativos
   **por valor de longitud** hasta `negative_ratio = 3.0` positivos por
   negativo en cada longitud.
6. **Abortos**: conjunto vacio o de una sola clase corta la etapa.

### Lo que dijo la primera corrida real

```
tras filtros de alfabeto y longitud: 25 557
tras deduplicado:                    18 683  (6 874 duplicados, 271 conflictos)
emparejado de longitudes:             8 089 -> 8 050 negativos
final:                               18 644  (10 594 +, 8 050 -)
longitud media:                      positivos 22.9, negativos 30.6
```

El numero que importa aqui es `length_matching.deficit = 23 732`: negativos que
el emparejado habria querido en alguna longitud concreta y **no existen en la
fuente**. Es decir, el emparejado de longitudes es estructuralmente inviable con
estas fuentes, y la longitud queda como confusor disponible. Consecuencia
tomada: **D-16**, ablacion `comp_physchem_nolen`, seccion 4.

Los 6 874 duplicados son esperables: DRAMP y APD3 se solapan mucho, y dedupicar
*antes* de clusterizar evita que un peptido repetido infle su propio cluster.

---

## 2. Etapa `split` — clustering por identidad y particion

**Modulo:** `split.py` (1 023 lineas, el modulo critico del proyecto).
**Artefactos:** `data/processed/splits.csv`, `data/processed/split_report.json`.

### El problema formal

Con secuencias homologas, una particion aleatoria produce un estimador sesgado
de la generalizacion: `E[métrica_test]` sobreestima el rendimiento en familias
no vistas, porque parte del test es casi-duplicado del train. El sesgo no es
pequeño; es la diferencia entre 0.98 y lo que de verdad hay.

### La solucion y su orden de preferencia

Clusterizar **todas** las secuencias (ambas clases juntas) al 40% de identidad y
asignar **clusteres completos**. Tres implementaciones, se usa la primera
disponible y el metodo real se registra en el reporte:

1. **MMseqs2** (`cluster_mmseqs2`), preferido:
   ```
   mmseqs easy-cluster in.fasta pref tmp --min-seq-id 0.4 -c 0.5 --cov-mode 1
                                         -s 7.5 --alignment-mode 3
   ```
   `--cov-mode 1` exige cobertura sobre la secuencia **mas corta**, que es lo
   correcto al mezclar peptidos de 5 y 50 residuos. `-s 7.5` sube la
   sensibilidad del prefiltro y `--alignment-mode 3` fuerza alineamiento
   completo: los dos son necesarios porque los valores por defecto estan
   ajustados a proteinas, no a peptidos cortos.
2. **CD-HIT** (`cluster_cdhit`): `-c 0.4 -n 2 -l 4`. `-n 2` es obligatorio por
   debajo de `-c 0.5`, y `-l` hay que bajarlo del defecto 10.
3. **Respaldo propio en Python** (`cluster_builtin`): clustering voraz
   incremental, para que `make split` funcione en una maquina limpia. Ordena por
   longitud descendente, compara cada secuencia contra representantes y usa un
   **prefiltro por k-meros** para evitar el coste cuadratico de alinear todo
   contra todo.

   El criterio del prefiltro es la **contencion**,
   `|A ∩ B| / min(|A|,|B|)`, con umbral 0.05, y esto no es un detalle de
   implementacion: la primera version uso **Jaccard con umbral 0.15** y eso
   **anulaba el umbral de identidad**. Jaccard castiga la diferencia de tamaño,
   asi que un peptido de 6 residuos contenido **entero** en uno de 50 daba 0.06 y
   se descartaba, pese a tener identidad y cobertura 1.00. El umbral efectivo
   era ~0.75 en vez de 0.40 (P-03).

   Sigue abierto que el umbral efectivo del respaldo esta en 0.55-0.60 y no en
   0.40 (P-23). Por eso el respaldo marca **`publishable: false`** en el reporte
   y MMseqs2 es obligatorio para numeros publicables (D-11). La tabla de recall
   medido esta en [04](04-division-clusters.md): a 0.40-0.45 de identidad real
   agrupa el ~50% de los pares.

La identidad se calcula con `Bio.Align.PairwiseAligner` en modo global con
BLOSUM62, y se define como `residuos identicos / longitud de la mas corta`.

### Asignacion de clusteres a splits

`assign_splits` resuelve un problema de particion con bloques indivisibles de
tamaño muy desigual y dos objetivos a la vez (cuotas 70/15/15 en **secuencias**
y balance de clases). Es NP-duro, asi que se usa una heuristica voraz:

1. **Siembra**: a cada split se le da primero un cluster pequeño de cada clase.
   Sin esto se producian splits degenerados con 0% de positivos en `val`
   (P-06), que es un fallo que no rompe nada visiblemente y hace que
   `choose_threshold` devuelva basura.
2. Ordenar el resto de clusteres de mayor a menor y asignar cada uno al split
   **relativamente** mas lejos de su cuota. Los grandes primero a proposito: al
   final, uno solo puede desbordar la cuota de `test`.
3. A igualdad de deficit, preferir el split cuyo balance de clases se degrada
   menos.

### Validaciones que abortan la etapa

`build_splits` no escribe nada si:

- las fracciones no suman 1 o les faltan claves;
- hay ids duplicados;
- **algun `cluster_id` aparece en dos splits** (el invariante del que depende la
  validez de todo lo demas, tambien como test en `tests/test_split.py`);
- algun split queda degenerado: le falta una clase, tiene <20 positivos o <10
  clusteres.

Y borra `results/metrics/predictions/*.csv`, porque una division nueva invalida
las predicciones viejas. Estaba **solo en el Makefile**, asi que `ampesm split`
y `ampesm all` publicaban metricas calculadas sobre ids inexistentes sin un solo
aviso, y en una variante la comparacion preregistrada salio de un archivo
obsoleto (P-20). La salvaguarda vive ahora donde ocurre la invalidacion, no en
el orquestador.

### La auditoria de fuga, que es el aporte metodologico del proyecto

**Ningun** algoritmo de clustering voraz garantiza ausencia de parientes entre
splits: todos, MMseqs2 incluido, comparan contra **representantes**, asi que A y
B al 50% de identidad pueden separarse si cada una se unio a un representante
distinto. No es un bug, es una propiedad del algoritmo. Por tanto la premisa del
proyecto hay que **medirla**, no asumirla.

`audit_train_test_identity` calcula, para cada secuencia de test, su identidad
maxima contra **todo** train.

**Y ese numero no se puede leer crudo.** La primera corrida dio 29.3% de
secuencias por encima de 0.4, lo que leido tal cual invalidaria el metodo
entero. Es espurio: con alineamiento global y huecos permitidos, dos peptidos
cortos **sin relacion** alinean al ~0.33 de media. El 29% estaba por debajo del
azar (P-12).

El control es **por secuencia**, y la eleccion importa: cada peptido de test se
baraja `R = 20` veces. Barajar conserva **longitud y composicion exactas** —que
es justo lo que determina la identidad de fondo— y destruye el parentesco. Una
secuencia se **marca** si su identidad observada supera **todas** sus propias
barajas.

Bajo H0 (sin fuga), cada secuencia se marca con probabilidad `1/(R+1)`, luego
el numero de marcadas `~ Binomial(m, 1/(R+1))` y
`scipy.stats.binomtest(n_flagged, m, 1/(R+1), alternative="greater")` da un
p-valor global:

```json
"leakage": {
  "n_flagged": 6, "n_expected_by_chance": 5.7, "excess": 0.3,
  "global_p_value": 0.51, "significant": false
}
```

**Lo que hay que leer es `global_p_value` y `excess`.** No
`fraction_over_threshold`, que dara ~30-50% y es ruido de fondo, ni
`excess_over_null_deprecated`.

**Y la auditoria corre dos veces: test contra train y `val` contra train**
(P-27). `val` fija los hiperparametros **y el umbral de decision**, y ese
umbral se aplica despues a test, asi que entra en MCC, F1, sensibilidad,
especificidad y precision: cinco de las siete metricas. Homologia val-train
sesga las metricas de test **aunque test este limpio**; la fuga no necesita
tocar test para contaminar test. `query_split` y `reference_split` son
parametros de la misma funcion, y `max_val: null` significa "el mismo tope que
test", no "audita val entero", que reintroduciria el coste cuadratico de P-19.

Costo tres intentos y los dos primeros **no tenian potencia**: devolvian numeros
plausibles y no detectaban nada con 20 parientes al 80% plantados a mano. El
segundo fue peor que inutil: p-valor por secuencia con correccion
Benjamini-Hochberg, que con `R = 20` tiene p minimo `1/21` y **no puede
sobrevivir** a BH sobre `m` tests salvo que `R > 20m`; daba cero detecciones en
todos los escenarios por aritmetica, no por ausencia de fuga (P-15, P-18). La
version actual esta validada plantando fuga: `p = 0.51` sin fuga, `p = 0.0006`
al 50% de identidad, `p ≈ 0` al 60% y 80%. De ahi **D-15**.

### Coste

La auditoria es `O(m · R · n_train)`. Pasar `R` de 1 a 20 multiplico por 20 el
coste de la etapa y la hizo cuadratica en el tamaño del conjunto: 526 s con
4 000 peptidos, ~2.3 h extrapoladas a 16 000 (P-19). Arreglar la correctitud y
romper la viabilidad es romperlo igual. Se topa con
`split.audit.max_test = 500`, suficiente de sobra para el test binomial, que
rompe el crecimiento cuadratico: medido 113 s con 3 909 peptidos. Si con el
conjunto real sigue siendo lento, `null_replicates: 20 → 10` lo reduce a la
mitad a cambio de resolucion en el p-valor.

El reporte incluye tambien `length_ks_pvalue` por split (Kolmogorov-Smirnov de
la longitud entre clases dentro de cada split), porque el emparejado de
longitudes es **global y no por split**, y un split puede quedar desbalanceado
aunque el conjunto global no lo este.

Y `positive_rate_drift` (P-26): la tasa de positivos no puede derivar entre
splits, porque PR-AUC, F1 y precision dependen de la prevalencia y el umbral
elegido en `val` se descalibra si `test` tiene otra. Habia un caso que pasaba
**todas** las demas validaciones con `train` al 6% y `test` al 29%. Avisa y no
aborta, por el criterio de D-18. Los splits vacios se excluyen: tienen
`positive_rate = 0.0` por construccion y falsearian la deriva.

---

## 3. Features tabulares

**Modulos:** `features/composition.py` (20 features),
`features/physchem.py` (11 features).

### Composicion

`count(aa) / |s|` para los 20 aminoacidos. Captura carga (K, R) e
hidrofobicidad (L, I, V, F, W) indirectamente y **nada del orden**. Es la
referencia correcta para preguntar si ESM-2 aporta algo sobre el orden de los
residuos.

### Fisicoquimicos

| Feature | Definicion |
|---|---|
| `length` | `|s|` |
| `net_charge` | Henderson-Hasselbalch sobre cadenas laterales ionizables + extremos, a pH 7.4 |
| `isoelectric_point` | biseccion sobre `net_charge(s, pH) = 0`, 60 iteraciones en `[0, 20]` |
| `gravy` | media de Kyte-Doolittle |
| `hydrophobic_moment` | momento de Eisenberg, 100°/residuo, **maximo sobre ventanas de 11** |
| `hydrophobic_fraction` | fraccion de `ACFILMVWY` |
| `aromaticity` | fraccion de `FWY` |
| `aliphatic_index` | indice de Ikai |
| `boman_index` | media de `ΔG` ciclohexano→agua (Radzicka y Wolfenden) |
| `interfacial_hydrophobicity` | media de la escala interfacial de Wimley-White |
| `charge_density` | `net_charge / length` |

El **momento hidrofobico** es el unico descriptor que mide geometria:
`μ_H = |Σ_i H(a_i) · e^{i·θ·i}|`, con `θ = 100°` (giro por residuo de la helice
α). Es grande cuando los residuos hidrofobicos se concentran en una cara, que es
la arquitectura anfipatica de la mayoria de AMPs helicoidales. Se toma el
**maximo sobre ventanas de 11 residuos**, no el valor global, porque un peptido
puede tener un segmento anfipatico dentro de una secuencia que globalmente no lo
sea; el global lo promediaria hasta desaparecer.

Las escalas estan **escritas literalmente en el modulo con su cita**, no
importadas de una libreria: se auditan sin salir del repositorio y no cambian
cuando una dependencia sube de version.

Esto importa porque una de ellas estuvo **mal**: `boman_index` usaba la escala
interfacial de Wimley-White con el signo invertido. LL-37 daba 0.35 en vez de
3.0 (un orden de magnitud) y la arginina —que debe dominar la escala con
14.92 kcal/mol y es justo lo que hace util al indice para AMPs cationicos—
quedaba al nivel de la asparagina. No fallaba nada: degradaba el modelo en
silencio y envenenaba la interpretacion (P-08). De ahi **D-14**: todo descriptor
entra con un test contra un valor publicado. Hoy verificados: carga de
magainina-2 (+2.98 vs +3 a +4 publicado) y LL-37 (+5.94 vs +6), Boman de LL-37
(3.00) y magainina-2 (0.42), GRAVY contra ExPASy, momento ≈ 0 en homopolimeros.

---

## 4. Etapa `baselines` — representaciones tabulares

**Modulo:** `models.py`. **Artefactos:** `results/models/*.json`,
`results/metrics/predictions/*.csv`.

Cuatro representaciones x tres modelos:

| Representacion | Dim | Para que esta |
|---|---|---|
| `composition` | 20 | ¿cuanto se consigue sin orden y sin fisica? |
| `physchem` | 11 | ¿cuanto consigue el mecanismo solo? |
| `comp_physchem` | 31 | **baseline principal preregistrado** |
| `comp_physchem_nolen` | 29 | ablacion sin `length` ni `charge_density` (D-16) |

Modelos: regresion logistica (con escalado; el coeficiente de `net_charge` debe
salir positivo y grande, y si no, hay un error en la tuberia —funciona como
verificacion de cordura), random forest y XGBoost. Rejilla en `params.yaml`,
seleccion **contra `val`**, nunca contra `test`.

`ModelUnavailable` envuelve el import de xgboost: en macOS sin `libomp` no carga
(P-09), y perder un baseline no debe tumbar la etapa entera, asi que
`run_baselines` avisa una vez y sigue.

### Sobre la ablacion y su asimetria

`comp_physchem_nolen` existe porque el `deficit` de 23 732 de la seccion 1
convierte la longitud en un atajo real (22.9 vs 30.6 residuos de media). La
alternativa era emparejar longitudes **dentro de cada split**, y se descarto:
los clusteres son atomicos, asi que hacerlo obliga a descartar secuencias
despues de la division, y hacerlo **solo si el KS sale bajo** convierte una
decision de diseño en un filtro guiado por un estadistico, que es la misma clase
de error que P-05.

La ablacion entra en `run_baselines` **siempre**, no condicionada al KS,
precisamente para que no sea una decision post-hoc.

**Asimetria que hay que declarar:** la ablacion quita la longitud al baseline,
no a ESM-2. El pooling por media atenua la longitud en los embeddings pero no la
elimina, asi que `comp_physchem_nolen` vs ESM-2 esta sesgado **a favor de
ESM-2**. Por eso la comparacion preregistrada sigue siendo contra
`comp_physchem` completo y la ablacion se reporta al lado como lectura del
confusor. Cerrar la asimetria pediria residualizar cada dimension del embedding
contra la longitud, y eso esta en [10](10-mejoras-futuras.md).

---

## 5. Etapa `embeddings` — ESM-2 congelado

**Modulo:** `embeddings.py`. **Artefactos:** `results/embeddings/*.npz`.

| Modelo | Capas | Dim | Peso |
|---|---|---|---|
| `facebook/esm2_t6_8M_UR50D` | 6 | 320 | ~31 MB |
| `facebook/esm2_t12_35M_UR50D` | 12 | 480 | ~140 MB |

ESM-2 es un transformer entrenado con *masked language modeling* sobre UniRef.
Aqui se usa **congelado**, como extractor de features: sin fine-tuning. Dos
razones: el presupuesto es CPU de un portatil, y fijar la representacion es lo
que hace que la comparacion con los descriptores sea una comparacion de
**representaciones** y no de presupuestos de computo.

Se usan los dos modelos mas pequeños a proposito. Los grandes (650M, 3B) no son
una opcion en CPU, y comparar la variante de 8M con la de 35M da gratis una
señal de si la escala del modelo mueve algo en esta tarea.

**El pooling es donde se cometen los errores silenciosos.** `masked_mean_pool`
promedia solo sobre posiciones reales: excluye padding **y** tokens especiales
(`cls`, `eos`). Promediar el padding diluye los peptidos cortos en proporcion a
cuanto padding tienen, lo que introduce una dependencia de la longitud
exactamente donde menos conviene (seccion 4).

El **cache** esta indexado por `seq_sha1` y el nombre de archivo codifica
modelo, capa y pooling, de modo que cambiar cualquiera de los tres no reutiliza
un cache incompatible.

Y el cache **nunca se escribia**. `np.savez_compressed` añade `.npz` al path si
no lo tiene, asi que el temporal `cache.npz.part` se guardaba como
`cache.npz.part.npz` y el `tmp.replace(path)` posterior lanzaba
`FileNotFoundError`. Habria perdido todos los embeddings en la primera corrida
real. Se arregla escribiendo a traves de un **descriptor de archivo abierto**,
que desactiva ese añadido de sufijo.

---

## 6. Etapas `train` y `evaluate`

**Modulos:** `models.py` (`run_embedding_models`), `evaluate.py`.
**Artefactos:** `results/metrics/summary.csv`,
`results/metrics/comparison_embeddings_vs_baseline.json`,
`results/figures/*.png`.

Sobre los embeddings: regresion logistica y un MLP pequeño (`fit_mlp`, PyTorch,
early stopping contra `val`).

### Metricas

`METRIC_NAMES = (roc_auc, pr_auc, mcc, f1, sensitivity, specificity, precision)`.
Varias a proposito: ROC-AUC es optimista con clases desbalanceadas, PR-AUC es
sensible a la prevalencia, MCC es la que usa toda la matriz de confusion y es la
**metrica principal preregistrada**.

El umbral se fija maximizando MCC **en `val`** (`choose_threshold`), que es la
razon de que existan tres splits y no dos. `choose_threshold` lanza `ValueError`
con validacion vacia o de una sola clase, en vez de devolver un umbral
arbitrario.

`all_metrics` devuelve `NaN` en las **siete** metricas cuando
`len(y) == 0` o hay una sola clase. Antes solo `roc_auc` lo hacia y las demas
devolvian numeros calculados sobre un remuestreo degenerado (P-07).

### Bootstrap de clusteres

`bootstrap_metrics` con `bootstrap_n = 2000` y **`bootstrap_unit = cluster`**.
Remuestrear secuencias trata como independientes a secuencias homologas y
produce intervalos demasiado estrechos; remuestrear **clusteres** respeta la
unidad real de independencia. El intervalo sale mas ancho, y es el ancho
correcto: la conclusion del proyecto depende de si dos intervalos se solapan.

Se reportan `n_boot_valid` / `n_boot_requested` y se avisa por debajo del 99%,
porque un remuestreo puede caer degenerado y un intervalo calculado sobre 1 400
de 2 000 replicas no es el intervalo que dice ser. `paired_difference` tiene la
misma contabilidad (P-16).

`json_safe` / `write_json` mapean no-finitos a `null` y escriben con
`allow_nan=False`: `json.dumps` por defecto emite `NaN` literal, que **no es
JSON valido** y revienta en cualquier consumidor estricto. El arreglo de P-07
agrando el agujero de una metrica a siete (P-21).

### La comparacion, y por que esta preregistrada

```python
PREREGISTERED_EMBEDDING = "facebook_esm2_t12_35M_UR50D__logreg"
PREREGISTERED_BASELINE  = "comp_physchem__xgboost"
```

La primera version elegia el par **por su MCC en test**. Eso es la maldicion del
ganador: comparar el mejor de 13 modelos contra el mejor del otro grupo infla la
diferencia **siempre**, incluso bajo H0 cierta, porque los maximos de muestra
estan sesgados hacia arriba. El par esta fijado por escrito en el README antes de
correr nada; `_select_comparison_pair` solo recurre a un fallback **por
validacion** (con desempate por nombre ordenado, P-17) si una etiqueta
preregistrada falta.

`evaluate_all` valida los ids de las predicciones contra `splits.csv` **antes**
de calcular nada, que es la otra mitad del arreglo de P-20.

La tabla completa de modelos es **descriptiva**. La unica inferencia del
proyecto es la diferencia emparejada de esa pareja, con su intervalo.

---

## 7. Etapa `interpret`

**Modulo:** `interpret.py`.

UMAP de las representaciones coloreado por **clase, fuente y split**: el coloreo
por fuente es el que delata si lo que separa el espacio es la biologia o la
procedencia del archivo.

`error_analysis` desagrega errores por categoria, **por fuente** y por longitud
(bins derivados de `min_length` y `max_length` en pasos de 10; estaban
parametrizados solo por el borde superior, P-22). Si el acierto en los positivos
de UniProt es mucho peor que en los de APD3, el modelo depende de la procedencia
y el resultado se lee con eso delante.

`_best_model_by_validation` elige el modelo a interpretar por `val_roc_auc`
leido de `results/models/*.json`, nunca por test.

`_check_features_hash` **aborta** la etapa si la seccion `features` de
`params.yaml` cambio desde que se entreno el modelo (P-24). `error_analysis`
recalcula los descriptores con la configuracion actual, asi que con otra
`moment_window` el analisis describiria un modelo que no existe, y no fallaba
nada: salia una figura plausible y equivocada. Aborta en vez de avisar porque
la salida es un PNG que se consume despues y viaja sin la consola (D-18).

---

## 8. Ingenieria

- **Entorno:** `uv`, Python `>=3.12,<3.13` (PyTorch no publica ruedas estables
  para 3.14), backend `hatchling`, `ruff`, `pytest`, `pre-commit`.
- **`params.yaml`** concentra todo lo ajustable, para que un cambio de
  configuracion aparezca en el diff y no enterrado en el codigo.
- **`pythonpath = ["src"]`** en `pyproject.toml`: el `.pth` del instalable
  editable dejo de procesarse en el entorno y los 117 tests pasaron a fallar en
  la coleccion con `No module named ampesm`, que parece un fallo del proyecto y
  es del entorno (P-13; reaparecio dos veces, y se arregla con
  `uv sync --reinstall-package ampesm`).
- **129 tests.** Los que de verdad sostienen el proyecto: el invariante de
  clusteres entre splits, los descriptores contra valores publicados, el
  detector de fuga contra fuga plantada, la ablacion contra el renombrado de
  sus propias columnas, la auditoria de val contra fuga plantada **solo en
  val**, y la deriva de prevalencia contra el caso conocido 6%/29%.
- **`make check-env`** comprueba que el paquete se importa antes de cada etapa
  (P-25, cuarta reaparicion). Vive en el Makefile y no dentro de `ampesm info`
  porque cuando el import falla, `info` tampoco arranca: una guarda que no
  puede ejecutarse cuando ocurre el fallo que vigila no es una guarda.
- **Dependencias del sistema:** `mmseqs2` (~30 MB, obligatorio para numeros
  publicables) y `libomp` (~2 MB, sin ella XGBoost no carga en macOS).

### Los tres patrones de fallo que dejo este proyecto

De los 24 problemas registrados, **11 los introdujeron los arreglos de los
anteriores**, y se agrupan en tres patrones que vale mas recordar que los bugs
concretos:

1. **Reintroducir el bug que acabas de arreglar, en otro sitio.** El corte
   temprano de `_max_identity_against` reintrodujo exactamente el caso
   "corto contenido en largo" que P-03 acababa de corregir, porque se copio la
   idea del prefiltro sin volver a comprobar sus condiciones: alli las
   secuencias estaban ordenadas por longitud y aqui no, asi que la cota
   superior era en realidad una cota inferior (P-14). Cambio el 41% de las
   identidades maximas medidas.
2. **Arreglar a medias.** P-20 puso la salvaguarda en el Makefile y no en el
   codigo; P-22 parametrizo un solo borde de los bins. Las dos se veian
   completas desde el sitio donde se hizo el cambio.
3. **Romper la viabilidad al arreglar la correctitud.** P-19: la auditoria
   correcta tardaba 2.3 horas. Un diagnostico que nadie va a correr no es un
   diagnostico.

---

## 9. Estado

`data` corrido con fuentes reales: **18 644 peptidos**, los cuatro controles del
reporte en orden. Las etapas 2-7 todavia no han visto datos reales y las tablas
de resultados estan vacias a proposito.

Hay una **fase 2 de generacion de peptidos** diseñada y sin implementar, cuyo
encuadre es medir cuanto engaña un generador al clasificador propio y no diseñar
farmacos: [17-diseno-generacion.md](17-diseno-generacion.md) y D-17.

Lo pendiente y lo que queda por decidir, en
[14-estado-y-continuacion.md](14-estado-y-continuacion.md).
