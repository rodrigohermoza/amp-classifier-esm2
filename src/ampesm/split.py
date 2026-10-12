"""Clustering por identidad de secuencia y division train/val/test.

Es la etapa que decide si los numeros del proyecto significan algo. La
justificacion completa esta en docs/04-division-clusters.md; aqui solo el como.

Invariante que debe cumplirse siempre: ningun cluster aparece en mas de un
split. Se verifica en tests/test_split.py. Si ese test falla, todas las
metricas del proyecto son invalidas.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from .config import PEPTIDES_PATH, SPLIT_REPORT_PATH, SPLITS_PATH, ensure_dirs, load_params

# --- Utilidades ------------------------------------------------------------


def write_fasta(sequences: dict[str, str], path: Path) -> None:
    """Escribe un FASTA con los ids como cabecera."""
    with open(path, "w", encoding="utf-8") as fh:
        for seq_id, seq in sequences.items():
            fh.write(f">{seq_id}\n{seq}\n")


def tool_version(cmd: list[str]) -> str | None:
    """Version de una herramienta externa, o None si no esta disponible."""
    if shutil.which(cmd[0]) is None:
        return None
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        text = (out.stdout or out.stderr).strip()
        return text.splitlines()[0] if text else "desconocida"
    except (subprocess.SubprocessError, OSError):
        return None


# --- Metodo 1: MMseqs2 -----------------------------------------------------


def cluster_mmseqs2(
    sequences: dict[str, str], identity: float, coverage: float
) -> dict[str, int] | None:
    """Agrupa con MMseqs2. Devuelve id -> cluster, o None si no esta instalado.

    --cov-mode 1 mide la cobertura respecto a la secuencia mas corta, que es lo
    que interesa al mezclar peptidos de 5 y de 50 residuos: con el modo por
    defecto, un peptido corto nunca alcanzaria la cobertura exigida frente a uno
    largo y quedaria en su propio cluster, perdiendo la agrupacion.
    """
    if shutil.which("mmseqs") is None:
        return None

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        fasta = tmpdir / "input.fasta"
        prefix = tmpdir / "clust"
        write_fasta(sequences, fasta)

        cmd = [
            "mmseqs",
            "easy-cluster",
            str(fasta),
            str(prefix),
            str(tmpdir / "tmp"),
            "--min-seq-id",
            str(identity),
            "-c",
            str(coverage),
            "--cov-mode",
            "1",
            # -s 7.5 sube la sensibilidad del prefiltro: con el valor por
            # defecto, a 40% de identidad sobre peptidos cortos MMseqs2 pierde
            # pares reales. --alignment-mode 3 calcula identidad y cobertura en
            # vez de estimarlas. Ambos cuestan tiempo, y aqui el conjunto es
            # pequeno, asi que la sensibilidad vale mas que la velocidad.
            "-s",
            "7.5",
            "--alignment-mode",
            "3",
            "-v",
            "1",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            stderr = getattr(exc, "stderr", "") or ""
            print(f"[split] MMseqs2 fallo: {stderr[:500]}")
            return None

        tsv = prefix.parent / f"{prefix.name}_cluster.tsv"
        if not tsv.exists():
            print("[split] MMseqs2 no produjo el archivo de clusteres esperado")
            return None

        # Formato: representante<TAB>miembro, una linea por miembro.
        rep_to_id: dict[str, int] = {}
        mapping: dict[str, int] = {}
        with open(tsv, encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 2:
                    continue
                rep, member = parts
                if rep not in rep_to_id:
                    rep_to_id[rep] = len(rep_to_id)
                mapping[member] = rep_to_id[rep]

    missing = set(sequences) - set(mapping)
    if missing:
        # No deberia pasar, pero si pasa es mejor darles cluster propio que
        # perderlos en silencio.
        next_id = max(mapping.values(), default=-1) + 1
        for i, seq_id in enumerate(sorted(missing)):
            mapping[seq_id] = next_id + i
        print(f"[split] AVISO: {len(missing)} secuencias sin cluster; se les asigno uno propio")

    return mapping


# --- Metodo 2: CD-HIT ------------------------------------------------------


def cluster_cdhit(
    sequences: dict[str, str], identity: float, coverage: float
) -> dict[str, int] | None:
    """Agrupa con CD-HIT. Devuelve id -> cluster, o None si no esta instalado.

    Dos ajustes obligados para peptidos cortos: el tamano de palabra tiene que
    bajar con el umbral de identidad (-n 2 para 0.4, que es el minimo que acepta
    CD-HIT), y -l 4 permite secuencias de 5 residuos, que el valor por defecto
    de 10 descartaria.
    """
    if shutil.which("cd-hit") is None:
        return None

    word_size = 5 if identity >= 0.7 else 4 if identity >= 0.6 else 3 if identity >= 0.5 else 2

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        fasta = tmpdir / "input.fasta"
        out = tmpdir / "out"
        write_fasta(sequences, fasta)

        cmd = [
            "cd-hit",
            "-i",
            str(fasta),
            "-o",
            str(out),
            "-c",
            str(identity),
            "-n",
            str(word_size),
            "-l",
            "4",
            "-aS",
            str(coverage),
            "-d",
            "0",
            "-M",
            "2000",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            stderr = getattr(exc, "stderr", "") or ""
            print(f"[split] CD-HIT fallo: {stderr[:500]}")
            return None

        clstr = Path(str(out) + ".clstr")
        if not clstr.exists():
            return None

        mapping: dict[str, int] = {}
        current = -1
        with open(clstr, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">Cluster"):
                    current = int(line.split()[1])
                elif ">" in line:
                    # '0  33aa, >AP00001... *'
                    name = line.split(">", 1)[1].split("...")[0].strip()
                    mapping[name] = current

    return mapping if mapping else None


# --- Metodo 3: respaldo propio --------------------------------------------


def kmer_set(sequence: str, k: int) -> frozenset[str]:
    """Conjunto de k-meros. Si la secuencia es mas corta que k, ella misma."""
    if len(sequence) < k:
        return frozenset({sequence})
    return frozenset(sequence[i : i + k] for i in range(len(sequence) - k + 1))


def _make_aligner():
    """Alineador global con parametros tipo BLASTP.

    Se importa Biopython aqui y no arriba para que el modulo se pueda importar
    sin Bio instalado cuando se usa MMseqs2.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def identity_and_coverage(seq_a: str, seq_b: str, aligner) -> tuple[float, float]:
    """Identidad y cobertura de dos secuencias, por alineamiento global.

    Identidad = residuos identicos / longitud de la mas corta.
    Cobertura  = columnas alineadas / longitud de la mas corta.

    Se normaliza por la mas corta y no por la longitud del alineamiento porque
    con peptidos de longitudes muy distintas lo segundo castigaria a un peptido
    corto que esta contenido entero en uno largo, que es precisamente el caso
    que queremos detectar como relacionado.
    """
    try:
        alignment = aligner.align(seq_a, seq_b)[0]
    except (ValueError, IndexError):
        return 0.0, 0.0

    shorter = min(len(seq_a), len(seq_b))
    matches = 0
    aligned_cols = 0
    for (a_start, a_end), (b_start, _b_end) in zip(*alignment.aligned):
        block = a_end - a_start
        aligned_cols += block
        matches += sum(1 for i in range(block) if seq_a[a_start + i] == seq_b[b_start + i])

    return matches / shorter, aligned_cols / shorter


def cluster_builtin(
    sequences: dict[str, str],
    identity: float,
    coverage: float,
    kmer_size: int = 3,
    prefilter: float = 0.05,
) -> tuple[dict[str, int], dict]:
    """Clustering voraz al estilo CD-HIT, sin herramientas externas.

    1. Ordenar de mas larga a mas corta (las largas son mejores representantes).
    2. Para cada secuencia, alinear solo contra los representantes que comparten
       suficientes k-meros, usando un indice invertido para no recorrer todos.
    3. Si alguno pasa el umbral de identidad y cobertura, unirse a ese cluster;
       si no, volverse representante de un cluster nuevo.

    El prefiltro hace esto viable: alinear cada secuencia contra miles de
    representantes en Python puro seria inviable, y las operaciones de conjuntos
    descartan la mayoria de pares antes de alinear.

    IMPORTANTE (P-03). El prefiltro se definia antes como un indice de Jaccard
    de 3-meros con umbral 0.15, y eso **anulaba el umbral de identidad**: a 45%
    de identidad real, entre el 0% y el 1% de los pares llegaba a alinearse, asi
    que el umbral efectivo era ~0.75 en vez de 0.40 y los parientes cercanos
    acababan en clusteres distintos, es decir en splits distintos. Dos cambios
    lo corrigen:

    - Se usa **contencion** `|A∩B| / min(|A|,|B|)` en vez de Jaccard. Jaccard
      castiga la diferencia de tamano, asi que un peptido de 5 residuos
      contenido ENTERO en uno de 50 daba 0.062 y se descartaba, pese a tener
      identidad y cobertura de 1.00. La contencion es ademas coherente con
      `identity_and_coverage`, que ya normaliza por la secuencia mas corta.
    - El umbral baja de 0.15 a 0.05. Medido con 8 000 secuencias: 2,4 s y 6 459
      clusteres con el valor viejo, 5,2 s y 3 592 clusteres con el nuevo. Casi
      el doble de agrupacion por unos segundos mas.

    Aun asi, este respaldo **no sirve para numeros publicables**: ver la nota de
    `build_splits` y docs/04.

    Devuelve (mapeo, diagnostico).
    """
    aligner = _make_aligner()
    order = sorted(sequences, key=lambda i: (-len(sequences[i]), i))

    kmer_index: dict[str, set[int]] = defaultdict(set)  # k-mero -> indices de repr.
    rep_ids: list[str] = []
    rep_kmers: list[frozenset[str]] = []
    mapping: dict[str, int] = {}
    alignments_done = 0
    candidates_skipped = 0

    for seq_id in order:
        seq = sequences[seq_id]
        kmers = kmer_set(seq, kmer_size)
        n_kmers = len(kmers)

        # Candidatos: representantes que comparten al menos un k-mero.
        shared: dict[int, int] = defaultdict(int)
        for km in kmers:
            for rep_idx in kmer_index.get(km, ()):
                shared[rep_idx] += 1

        # Se prueban primero los de mayor solapamiento: si uno cercano pasa el
        # umbral, se corta y se ahorran los alineamientos restantes.
        candidates = sorted(shared.items(), key=lambda kv: -kv[1])
        assigned: int | None = None

        for rep_idx, n_shared in candidates:
            # Corte temprano. OJO con la direccion de la cota (P-14):
            # contencion = n_shared / min(n_kmers, n_kmers_rep). Si el
            # representante tiene MENOS k-meros que la consulta, el denominador
            # real es menor y la contencion real es MAYOR que n_shared/n_kmers,
            # asi que esa cantidad es una cota INFERIOR, no superior, y cortar
            # con ella descartaria pares validos.
            #
            # Aqui si es seguro cortar porque las secuencias se procesan de mas
            # larga a mas corta, de modo que todo representante tiene al menos
            # tantos k-meros como la consulta y la cota es exacta. La condicion
            # se comprueba explicitamente en vez de confiar en el orden.
            if len(rep_kmers[rep_idx]) >= n_kmers and (n_shared / max(n_kmers, 1)) < prefilter:
                candidates_skipped += len(candidates) - candidates.index((rep_idx, n_shared))
                break

            denom = min(n_kmers, len(rep_kmers[rep_idx]))
            if denom == 0 or (n_shared / denom) < prefilter:
                continue

            ident, cov = identity_and_coverage(sequences[rep_ids[rep_idx]], seq, aligner)
            alignments_done += 1
            if ident >= identity and cov >= coverage:
                assigned = rep_idx
                break

        if assigned is None:
            assigned = len(rep_ids)
            rep_ids.append(seq_id)
            rep_kmers.append(kmers)
            for km in kmers:
                kmer_index[km].add(assigned)

        mapping[seq_id] = assigned

    diagnostics = {
        "n_clusters": len(rep_ids),
        "alignments": alignments_done,
        "candidates_skipped_by_bound": candidates_skipped,
        "kmer_size": kmer_size,
        "prefilter": prefilter,
        "prefilter_criterion": "containment",
    }
    print(
        f"[split] respaldo propio: {len(rep_ids)} clusteres, "
        f"{alignments_done:,} alineamientos realizados"
    )
    return mapping, diagnostics


# --- Asignacion de clusteres a splits -------------------------------------


def assign_splits(df: pd.DataFrame, fractions: dict[str, float], seed: int) -> dict[int, str]:
    """Reparte clusteres completos entre train/val/test.

    Heuristica voraz: se recorren los clusteres de mayor a menor tamano y cada
    uno va al split que esta relativamente mas lejos de cumplir su cuota; a
    igualdad, al que degrada menos el balance de clases.

    Los grandes van primero a proposito. Si se dejan para el final, un unico
    cluster grande puede desbordar la cuota de test y descuadrar las
    proporciones por completo.

    El reparto exacto es un problema de particion (NP-duro) y no hace falta
    resolverlo: lo que importa es que no haya fuga, no que los porcentajes
    salgan exactos. Los porcentajes reales se registran en el reporte.

    SIEMBRA PREVIA (P-06). La version anterior solo penalizaba el desbalance de
    clases con peso 0.5, y eso no bastaba: arrancando con deficit 1.0 en los
    tres splits, los primeros clusteres grandes se repartian uno por split sin
    que el balance pudiera vencer al deficit. Con familias de AMPs, donde la
    etiqueta correlaciona con el cluster (el caso real, no uno adverso
    inventado), salia val con 0% de positivos y test con 100%. Ahora el balance
    es **restriccion y no penalizacion**: primero se siembra cada split con al
    menos un cluster de cada clase, y solo despues se aplica el criterio de
    deficit.
    """
    rng = np.random.default_rng(seed)
    total = len(df)
    targets = {k: v * total for k, v in fractions.items()}

    sizes = df.groupby("cluster_id").size()
    pos_counts = df[df.label == 1].groupby("cluster_id").size()

    clusters = [
        {
            "id": int(cid),
            "n": int(sizes[cid]),
            "n_pos": int(pos_counts.get(cid, 0)),
        }
        for cid in sizes.index
    ]
    # El desempate aleatorio evita que el orden dependa del id del cluster, que
    # a su vez depende del orden del archivo de entrada.
    rng.shuffle(clusters)
    clusters.sort(key=lambda c: -c["n"])

    assigned: dict[int, str] = {}
    current = {k: 0 for k in fractions}
    current_pos = {k: 0 for k in fractions}
    global_pos_rate = float(df.label.mean())

    # --- Siembra: garantizar ambas clases en cada split -------------------
    # Se reparten los clusteres que contienen positivos, y luego los que
    # contienen negativos, uno por split, de menor a mayor tamano. Se eligen los
    # mas pequenos para no gastar la cuota de val y test en un cluster enorme:
    # la siembra solo tiene que asegurar presencia, no volumen.
    with_pos = sorted((c for c in clusters if c["n_pos"] > 0), key=lambda c: c["n"])
    with_neg = sorted((c for c in clusters if c["n"] - c["n_pos"] > 0), key=lambda c: c["n"])

    for pool in (with_pos, with_neg):
        for split in fractions:
            # Si este split ya tiene la clase cubierta por la siembra anterior,
            # no hace falta sembrarla otra vez.
            already = (
                current_pos[split] > 0
                if pool is with_pos
                else (current[split] - current_pos[split]) > 0
            )
            if already:
                continue
            seed_cluster = next((c for c in pool if c["id"] not in assigned), None)
            if seed_cluster is None:
                continue  # no hay material; la validacion posterior lo detectara
            assigned[seed_cluster["id"]] = split
            current[split] += seed_cluster["n"]
            current_pos[split] += seed_cluster["n_pos"]

    # --- Reparto por deficit del resto ------------------------------------
    for cluster in clusters:
        if cluster["id"] in assigned:
            continue
        best_split, best_score = None, None
        for split in fractions:
            # Deficit relativo: que fraccion de su cuota le falta a este split.
            deficit = (targets[split] - current[split]) / max(targets[split], 1.0)

            # Desviacion del balance de clases si se le anade este cluster.
            n_after = current[split] + cluster["n"]
            pos_after = current_pos[split] + cluster["n_pos"]
            balance_penalty = abs((pos_after / n_after) - global_pos_rate)

            score = deficit - 0.5 * balance_penalty
            if best_score is None or score > best_score:
                best_split, best_score = split, score

        assigned[cluster["id"]] = best_split
        current[best_split] += cluster["n"]
        current_pos[best_split] += cluster["n_pos"]

    return assigned


# --- Verificacion directa de la premisa del proyecto ----------------------


def _max_identity_against(
    queries: list[str],
    subject_seqs: list[str],
    subject_kmers: list[frozenset[str]],
    kmer_index: dict[str, set[int]],
    aligner,
    kmer_size: int,
    prefilter: float,
) -> list[tuple[float, int | None]]:
    """Para cada consulta, su identidad maxima contra el conjunto sujeto."""
    out: list[tuple[float, int | None]] = []
    for query in queries:
        kmers = kmer_set(query, kmer_size)
        n_kmers = len(kmers)

        shared: dict[int, int] = defaultdict(int)
        for km in kmers:
            for idx in kmer_index.get(km, ()):
                shared[idx] += 1

        best_identity = 0.0
        best_partner: int | None = None
        for idx, n_shared in sorted(shared.items(), key=lambda kv: -kv[1]):
            # Sin 'break' aqui, a diferencia de cluster_builtin (P-14). Las
            # consultas son secuencias de test y los sujetos de train, en
            # cualquier orden de longitud, asi que n_shared/n_kmers NO acota la
            # contencion por arriba: un sujeto corto puede tener contencion 1.00
            # con una cota de 0.12. Cortar el bucle perdia ademas todos los
            # sujetos cortos restantes y falseaba la auditoria a la baja.
            denom = min(n_kmers, len(subject_kmers[idx]))
            if denom == 0 or (n_shared / denom) < prefilter:
                continue
            ident, _cov = identity_and_coverage(subject_seqs[idx], query, aligner)
            if ident > best_identity:
                best_identity = ident
                best_partner = idx
        out.append((best_identity, best_partner))
    return out


def audit_train_test_identity(
    df: pd.DataFrame,
    threshold: float,
    kmer_size: int = 3,
    prefilter: float = 0.05,
    max_test: int | None = None,
    seed: int = 42,
    null_model: bool = True,
    null_replicates: int = 20,
) -> dict:
    """Identidad maxima de cada secuencia de test contra TODO train.

    Esta es la verificacion que faltaba (P-04), y vale mas que cualquier ajuste
    de parametros del clustering.

    docs/04 afirma que "ninguna secuencia de prueba tiene un pariente cercano en
    entrenamiento". Ningun algoritmo de clustering voraz garantiza eso, ni el
    propio ni MMseqs2 ni CD-HIT: todos comparan cada secuencia solo contra
    **representantes**, asi que A y B al 50% de identidad pueden acabar en
    clusteres distintos si A se unio a R1 y B a R2 sin superar el umbral contra
    el representante de la otra. No es un bug, es una propiedad conocida de los
    algoritmos tipo CD-HIT.

    Asi que en vez de confiar en la agrupacion, se mide el resultado
    directamente: para cada secuencia de test, su identidad maxima contra el
    conjunto de entrenamiento completo. Eso convierte la premisa central del
    proyecto de una afirmacion en un numero publicable: no "dividimos por
    clusteres, asi que no hay fuga", sino "el X% de las secuencias de test tiene
    un pariente por encima del umbral en train, y aqui esta la distribucion".

    Se reutiliza el indice invertido de k-meros con un prefiltro laxo. **El
    prefiltro no es inocuo**: medido contra fuerza bruta, cambia la identidad
    maxima de una fraccion apreciable de las consultas y siempre a la baja, asi
    que el resultado es una cota inferior, no una medida exacta. Se acepta
    porque la fuerza bruta sobre miles de secuencias es inviable en Python, y
    porque el nulo se calcula con el mismo prefiltro: el sesgo afecta por igual
    a la medicion y a su control, que es lo que hace valida la comparacion.
    """
    aligner = _make_aligner()
    train = df[df.split == "train"]
    test = df[df.split == "test"]

    if train.empty or test.empty:
        return {"error": "train o test vacio; no se puede auditar"}

    subject_seqs = train.sequence.tolist()
    subject_kmers = [kmer_set(s, kmer_size) for s in subject_seqs]

    kmer_index: dict[str, set[int]] = defaultdict(set)
    for idx, kms in enumerate(subject_kmers):
        for km in kms:
            kmer_index[km].add(idx)

    test_rows = test
    sampled = False
    if max_test is not None and len(test) > max_test:
        test_rows = test.sample(n=max_test, random_state=seed)
        sampled = True

    test_seqs = test_rows.sequence.tolist()
    results = _max_identity_against(
        test_seqs, subject_seqs, subject_kmers, kmer_index, aligner, kmer_size, prefilter
    )
    arr = np.asarray([r[0] for r in results])

    offenders = []
    for row, (ident, partner) in zip(test_rows.itertuples(), results):
        if ident >= threshold and len(offenders) < 25:
            offenders.append(
                {
                    "test_id": row.id,
                    "test_sequence": row.sequence,
                    "max_identity": round(ident, 4),
                    "closest_train_id": train.iloc[partner].id if partner is not None else None,
                }
            )

    report = {
        "threshold": threshold,
        "n_test_audited": int(len(arr)),
        "sampled": sampled,
        "n_over_threshold": int((arr >= threshold).sum()),
        "fraction_over_threshold": round(float((arr >= threshold).mean()), 4) if len(arr) else 0.0,
        "max_identity_distribution": {
            "mean": round(float(arr.mean()), 4),
            "median": round(float(np.median(arr)), 4),
            "p90": round(float(np.quantile(arr, 0.90)), 4),
            "p99": round(float(np.quantile(arr, 0.99)), 4),
            "max": round(float(arr.max()), 4),
        },
        "histogram": {
            f"{lo:.1f}-{lo + 0.1:.1f}": int(((arr >= lo) & (arr < lo + 0.1)).sum())
            for lo in np.arange(0.0, 1.0, 0.1)
        },
        "offenders": offenders,
    }

    # --- Distribucion nula, por secuencia -------------------------------
    # P-12 establecio que el numero crudo es ininterpretable sin control, pero
    # el primer control elegido (la diferencia de proporciones agregadas) no
    # tenia potencia (P-15): con el ruido de fondo, ambas proporciones estan
    # saturadas cerca de 0.8, asi que plantar un 17% de parientes al 80% de
    # identidad movia el exceso solo +0.017, indistinguible de cero.
    #
    # El control correcto es **por secuencia**: para cada peptido de test se
    # baraja su propia secuencia R veces (conservando longitud y composicion
    # exactas, que es lo que determina su identidad de fondo) y se mide la
    # identidad maxima de cada baraja contra train. Eso da la distribucion nula
    # DE ESA SECUENCIA, y su p-valor empirico es la fraccion de barajas que
    # igualan o superan lo observado. Asi cada peptido se compara consigo mismo
    # y no contra un promedio del conjunto.
    if null_model:
        rng = np.random.default_rng(seed)
        per_seq_null_max: list[list[float]] = [[] for _ in test_seqs]

        for _ in range(null_replicates):
            shuffled = []
            for s in test_seqs:
                chars = list(s)
                rng.shuffle(chars)
                shuffled.append("".join(chars))
            rep = _max_identity_against(
                shuffled, subject_seqs, subject_kmers, kmer_index, aligner, kmer_size, prefilter
            )
            for i, (ident, _p) in enumerate(rep):
                per_seq_null_max[i].append(ident)

        null_all = np.asarray([v for lst in per_seq_null_max for v in lst])
        m = len(arr)

        # Como se mide la fuga (P-18). El primer intento fue un p-valor por
        # secuencia con correccion de Benjamini-Hochberg, y era INCAPAZ de
        # detectar nada: con R barajas el p-valor minimo posible es 1/(R+1), y
        # tras corregir por m pruebas hace falta 1/(R+1) * m < 0.05, es decir
        # R > 20*m. Con m=120 secuencias harian falta mas de 2 400 barajas por
        # secuencia. El test daba cero por construccion, no por ausencia de
        # fuga, que es la peor forma posible de equivocarse.
        #
        # Lo que si funciona: marcar una secuencia cuando su identidad observada
        # supera TODAS sus propias barajas, y contrastar el NUMERO de marcadas
        # contra lo que se espera por azar. Bajo la hipotesis nula, cada
        # secuencia se marca con probabilidad 1/(R+1), asi que el conteo sigue
        # una binomial y un test binomial da un p-valor global con potencia
        # real. Se gana potencia agregando el conteo en vez de exigir
        # significacion secuencia a secuencia.
        flagged = np.array(
            [arr[i] > max(per_seq_null_max[i], default=0.0) for i in range(len(arr))]
        )
        n_flagged = int(flagged.sum())
        p_flag = 1.0 / (null_replicates + 1)
        expected = m * p_flag

        from scipy.stats import binomtest

        global_p = float(binomtest(n_flagged, m, p_flag, alternative="greater").pvalue)
        noise_cut = float(np.quantile(null_all, 0.99))

        report["null_model"] = {
            "description": (
                "Por secuencia: cada peptido de test se baraja null_replicates veces "
                "(conserva longitud y composicion, destruye el parentesco) y se mide "
                "su identidad maxima contra train."
            ),
            "replicates": null_replicates,
            "mean": round(float(null_all.mean()), 4),
            "p90": round(float(np.quantile(null_all, 0.90)), 4),
            "p99": round(noise_cut, 4),
            "max": round(float(null_all.max()), 4),
            "fraction_over_threshold": round(float((null_all >= threshold).mean()), 4),
        }

        # --- El numero que hay que leer -----------------------------------
        report["leakage"] = {
            "n_flagged": n_flagged,
            "n_expected_by_chance": round(expected, 2),
            "excess": round(n_flagged - expected, 2),
            "fraction_flagged": round(float(n_flagged / m), 4) if m else 0.0,
            "global_p_value": round(global_p, 6),
            "significant": bool(global_p < 0.05),
            "method": (
                "Una secuencia se marca si su identidad maxima supera las "
                f"{null_replicates} barajas de si misma. Bajo la nula eso ocurre con "
                "probabilidad 1/(R+1); el conteo se contrasta con un test binomial."
            ),
            "note": (
                "ESTE es el numero de fuga. n_flagged por si solo no significa nada: "
                "se espera n_expected_by_chance por azar. Lo que importa es el exceso "
                "y su p-valor global. Verificado contra fuga plantada; ver docs/09, P-18."
            ),
        }
        report["n_above_noise_threshold"] = int((arr > noise_cut).sum())
        report["suggested_threshold_above_noise"] = round(noise_cut, 4)
        report["flagged_ids"] = [test_rows.iloc[i].id for i in range(len(arr)) if flagged[i]][:50]

        # Se conserva el agregado por transparencia, marcado como lo que es.
        report["excess_over_null_deprecated"] = {
            "value": round(float((arr >= threshold).mean() - (null_all >= threshold).mean()), 4),
            "warning": (
                "NO usar como medida de fuga: es la diferencia de dos proporciones "
                "saturadas y no detecta fuga plantada ni al 80% de identidad. Se "
                "conserva solo para poder comparar con corridas antiguas. Ver P-15."
            ),
        }

    report["note"] = (
        "Identidad maxima de cada secuencia de test contra todo train. El numero "
        "de fuga es leakage.global_p_value y su exceso, NO fraction_over_threshold: "
        "la metrica de identidad tiene un ruido de fondo alto en peptidos cortos. "
        "Ver docs/09, P-04, P-12, P-15 y P-18."
    )
    return report


# --- Orquestacion ---------------------------------------------------------


def build_splits(force_method: str | None = None) -> pd.DataFrame:
    """Agrupa y divide. Escribe splits.csv y split_report.json."""
    ensure_dirs()
    params = load_params()
    scfg = params.split

    identity = float(scfg.get("identity_threshold", 0.4))
    coverage = float(scfg.get("min_coverage", 0.5))
    fractions = dict(scfg.get("fractions", {"train": 0.70, "val": 0.15, "test": 0.15}))
    methods = list(scfg.get("methods", ["mmseqs2", "cdhit", "builtin"]))
    bcfg = dict(scfg.get("builtin", {}))

    # P-10e: unas fracciones mal puestas dejarian la columna split en None, la
    # verificacion de fuga pasaria (0 valores unicos por cluster) y el fallo
    # aparaceria dos etapas mas tarde como un error opaco de sklearn.
    if set(fractions) != {"train", "val", "test"}:
        raise ValueError(
            f"split.fractions debe tener exactamente train, val y test; hay: {sorted(fractions)}"
        )
    if abs(sum(fractions.values()) - 1.0) > 0.01:
        raise ValueError(f"split.fractions debe sumar 1.0; suma {sum(fractions.values())}")

    if not PEPTIDES_PATH.exists():
        raise FileNotFoundError(f"Falta {PEPTIDES_PATH}. Corre 'make data' primero.")
    df = pd.read_csv(PEPTIDES_PATH)

    # P-10b: la deduplicacion de clean.py es por SECUENCIA, no por id, asi que
    # dos secuencias distintas bajo el mismo accession producen ids repetidos.
    # El dict de abajo pisaria una de las dos: no fallaria nada, pero una
    # secuencia no se agruparia, ambas filas heredarian el mismo cluster_id y el
    # cruce de interpret.py duplicaria filas inflando los conteos.
    if df.id.duplicated().any():
        dups = df.id[df.id.duplicated()].unique()[:5].tolist()
        raise RuntimeError(
            f"{int(df.id.duplicated().sum())} ids duplicados en {PEPTIDES_PATH.name}, "
            f"por ejemplo {dups}. Desambigua el id en clean.py (p. ej. con seq_sha1)."
        )

    sequences = dict(zip(df.id, df.sequence))
    print(f"[split] {len(sequences)} secuencias, umbral de identidad {identity}")

    if force_method:
        methods = [force_method]

    mapping: dict[str, int] | None = None
    used_method = None
    used_version = None
    cluster_diagnostics: dict = {}

    for method in methods:
        if method == "mmseqs2":
            version = tool_version(["mmseqs", "version"])
            if version is None:
                print("[split] mmseqs no disponible")
                continue
            print(f"[split] usando MMseqs2 ({version})")
            mapping = cluster_mmseqs2(sequences, identity, coverage)
        elif method == "cdhit":
            version = tool_version(["cd-hit", "-h"])
            if version is None:
                print("[split] cd-hit no disponible")
                continue
            print(f"[split] usando CD-HIT ({version})")
            mapping = cluster_cdhit(sequences, identity, coverage)
        elif method == "builtin":
            version = "interno"
            print(
                "[split] AVISO: usando el respaldo propio en Python.\n"
                "[split] Este metodo NO sirve para numeros publicables: agrupa menos\n"
                "[split] que MMseqs2 y deja fuga residual por homologia. Instala\n"
                "[split] MMseqs2 ('brew install mmseqs2') antes de publicar resultados.\n"
                "[split] Ver docs/09-problemas-encontrados.md, P-03."
            )
            mapping, cluster_diagnostics = cluster_builtin(
                sequences,
                identity,
                coverage,
                kmer_size=int(bcfg.get("kmer_size", 3)),
                prefilter=float(bcfg.get("kmer_prefilter", 0.05)),
            )
        else:
            raise ValueError(f"Metodo de clustering desconocido: {method!r}")

        if mapping:
            used_method, used_version = method, version
            break

    if not mapping:
        raise RuntimeError(f"Ningun metodo de clustering funciono. Probados: {methods}")

    df["cluster_id"] = df.id.map(mapping)
    if df.cluster_id.isna().any():
        missing = int(df.cluster_id.isna().sum())
        raise RuntimeError(f"{missing} secuencias quedaron sin cluster: hay un bug en el mapeo")
    df["cluster_id"] = df.cluster_id.astype(int)

    assignment = assign_splits(df, fractions, params.seed)
    df["split"] = df.cluster_id.map(assignment)

    if df.split.isna().any():
        raise RuntimeError(f"{int(df.split.isna().sum())} filas sin split asignado")

    # Verificacion del invariante, aqui y no solo en los tests: si esto falla,
    # es mejor no escribir el archivo que escribir uno invalido.
    per_cluster = df.groupby("cluster_id").split.nunique()
    if (per_cluster > 1).any():
        culprits = per_cluster[per_cluster > 1].index.tolist()[:10]
        raise RuntimeError(f"FUGA: clusteres repartidos entre splits: {culprits}")

    # P-06: splits degenerados. Mismo principio que arriba, mejor abortar que
    # escribir un archivo que produciria metricas sin sentido. Un val monoclase
    # revienta roc_auc_score a mitad de la rejilla de hiperparametros; un test
    # monoclase es peor, porque da NaN y MCC 0.0 SIN aviso y se publicarian
    # numeros basura.
    MIN_POSITIVES = 20
    MIN_CLUSTERS = 10
    problems: list[str] = []
    for split in fractions:
        sub = df[df.split == split]
        n_pos = int((sub.label == 1).sum())
        n_neg = int((sub.label == 0).sum())
        if len(sub) == 0:
            problems.append(f"{split} esta vacio")
            continue
        if n_pos == 0 or n_neg == 0:
            problems.append(f"{split} tiene una sola clase ({n_pos} positivos, {n_neg} negativos)")
        elif n_pos < MIN_POSITIVES:
            problems.append(f"{split} tiene solo {n_pos} positivos (minimo {MIN_POSITIVES})")
        n_clusters = int(sub.cluster_id.nunique())
        if n_clusters < MIN_CLUSTERS:
            problems.append(
                f"{split} tiene solo {n_clusters} clusteres (minimo {MIN_CLUSTERS}); "
                "el bootstrap por cluster seria inutil"
            )
    if problems:
        raise RuntimeError(
            "Division degenerada, no se escribe splits.csv:\n  - "
            + "\n  - ".join(problems)
            + "\n\nCausa habitual: muy pocos clusteres, o unos pocos clusteres enormes que "
            "concentran una clase. Revisa data/processed/clean_report.json y el umbral de "
            "identidad. Ver docs/09-problemas-encontrados.md, P-06."
        )

    df.to_csv(SPLITS_PATH, index=False)

    cluster_sizes = df.groupby("cluster_id").size()
    mixed = df.groupby("cluster_id").label.nunique()
    report = {
        "method": used_method,
        "tool_version": used_version,
        "identity_threshold": identity,
        "min_coverage": coverage,
        "seed": params.seed,
        "n_sequences": int(len(df)),
        "n_clusters": int(df.cluster_id.nunique()),
        "cluster_size": {
            "max": int(cluster_sizes.max()),
            "mean": float(cluster_sizes.mean()),
            "median": float(cluster_sizes.median()),
            "singletons": int((cluster_sizes == 1).sum()),
        },
        "mixed_label_clusters": int((mixed > 1).sum()),
        "mixed_label_note": (
            "Clusteres con secuencias de ambas clases: secuencias parecidas con "
            "etiquetas distintas. Probablemente es donde se concentra el ruido de "
            "etiquetado. Ver docs/08-decisiones.md, D-03."
        ),
        "publishable": used_method != "builtin",
        "publishable_note": (
            "false significa que se uso el respaldo propio en Python, que agrupa "
            "menos que MMseqs2 y deja fuga residual. No publiques metricas con "
            "publishable=false. Ver docs/09, P-03."
        ),
        "cluster_diagnostics": cluster_diagnostics,
        "splits": {},
    }
    for split in fractions:
        sub = df[df.split == split]
        pos = sub[sub.label == 1]
        neg = sub[sub.label == 0]
        report["splits"][split] = {
            "n": int(len(sub)),
            "fraction_actual": round(len(sub) / len(df), 4),
            "fraction_target": fractions[split],
            "n_positive": int(len(pos)),
            "n_negative": int(len(neg)),
            "positive_rate": round(float(sub.label.mean()), 4) if len(sub) else 0.0,
            "n_clusters": int(sub.cluster_id.nunique()),
            # P-05 del informe de metodologia: el emparejado de longitudes es
            # global, asi que puede no cumplirse DENTRO de cada split. Si las
            # medias se separan en test, la longitud vuelve a ser un atajo y
            # `length` es literalmente la primera feature fisicoquimica.
            "length_mean_positive": round(float(pos.length.mean()), 2) if len(pos) else None,
            "length_mean_negative": round(float(neg.length.mean()), 2) if len(neg) else None,
            "length_ks_pvalue": _length_ks_pvalue(pos, neg),
        }

    report["length_matching_note"] = (
        "length_ks_pvalue es un test de Kolmogorov-Smirnov de dos muestras entre "
        "las longitudes de positivos y negativos DENTRO de cada split. Un p-valor "
        "bajo en test significa que la longitud sigue separando las clases y que "
        "parte del rendimiento puede venir de ahi. Ver docs/09, P-10."
    )

    # P-04: la verificacion directa de la premisa del proyecto.
    print("[split] auditando identidad maxima test-vs-train...")
    report["train_test_identity_audit"] = audit_train_test_identity(
        df,
        threshold=identity,
        kmer_size=int(bcfg.get("kmer_size", 3)),
        prefilter=float(bcfg.get("kmer_prefilter", 0.05)),
        seed=params.seed,
    )

    SPLIT_REPORT_PATH.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"[split] {report['n_clusters']} clusteres ({report['cluster_size']['singletons']} de 1)")
    print(f"[split] {report['mixed_label_clusters']} clusteres con ambas clases")
    for split, info in report["splits"].items():
        print(
            f"[split] {split:5s}: {info['n']:6d} ({info['fraction_actual']:.1%}) "
            f"{info['n_positive']}+ / {info['n_negative']}-  "
            f"tasa positiva {info['positive_rate']:.3f}"
        )

    audit = report["train_test_identity_audit"]
    if "error" not in audit:
        print(
            f"[split] identidad test-vs-train: {audit['fraction_over_threshold']:.1%} "
            f"por encima de {identity} (crudo)"
        )
        if "null_model" in audit:
            null = audit["null_model"]
            leak = audit["leakage"]
            print(
                f"[split]   ruido de fondo ({null['replicates']} barajas por secuencia): "
                f"media {null['mean']:.3f}, p99 {null['p99']:.3f}"
            )
            print(
                f"[split]   marcadas: {leak['n_flagged']} de {audit['n_test_audited']} "
                f"(esperadas por azar {leak['n_expected_by_chance']}), "
                f"exceso {leak['excess']:+.1f}"
            )
            print(
                f"[split]   FUGA: p global = {leak['global_p_value']:.4g} "
                f"({'SIGNIFICATIVA' if leak['significant'] else 'no significativa'})"
            )
            if leak["significant"]:
                print(
                    "[split] AVISO: hay fuga por homologia estadisticamente "
                    "significativa. Revisa el metodo de clustering antes de publicar."
                )

    if not report["publishable"]:
        print("[split] AVISO: publishable=false (se uso el respaldo propio). Ver docs/09, P-03.")

    print(f"[split] {SPLITS_PATH}")
    return df


def _length_ks_pvalue(pos: pd.DataFrame, neg: pd.DataFrame) -> float | None:
    """p-valor de Kolmogorov-Smirnov entre las longitudes de ambas clases."""
    if len(pos) < 2 or len(neg) < 2:
        return None
    from scipy.stats import ks_2samp

    return round(float(ks_2samp(pos.length.to_numpy(), neg.length.to_numpy()).pvalue), 6)
