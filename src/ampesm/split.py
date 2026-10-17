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

from .config import (
    METRICS_DIR,
    PEPTIDES_PATH,
    SPLIT_REPORT_PATH,
    SPLITS_PATH,
    ensure_dirs,
    load_params,
)

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

    `--cov-mode 1` mide la cobertura respecto al **objetivo**, que en el
    agrupamiento es el representante, es decir la secuencia larga. Esta frase
    decia antes lo contrario ("respecto a la mas corta") y afirmaba que por eso
    un peptido corto se agrupa con uno largo. Es falso, y medido (P-35):

        larga            GIPCGESCVWIPCISGAIGCSCKSKVCYRN   30 residuos
        corta_subcadena  IPCISGAIGC                       subcadena exacta
        -> con -c 0.5 --cov-mode 1 quedan en clusteres DISTINTOS

    Es decir, la definicion efectiva de homologo de este agrupamiento es
    identidad >= `identity` **y** alineamiento que cubra >= `coverage` de la
    secuencia larga. Las subcadenas cortas no las agrupa, y de ahi vino P-30:
    las tapa aparte la fusion por contencion de D-20, no este modo de cobertura.

    La auditoria de fuga exige la misma definicion, por eso le pasa el mismo
    `min_coverage` (P-32).
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
            # --cluster-mode 1 (componente conexa) NO es el valor por defecto, y
            # el defecto es incorrecto para lo que este proyecto afirma (P-29).
            #
            # El modo por defecto resulto ser el 2, incremental voraz: cada
            # secuencia se asigna al primer representante con el que alinea. Eso
            # garantiza que todo miembro se parece a SU representante, y no
            # garantiza nada entre miembros de clusteres distintos. Medido sobre
            # los datos reales: DLPSGWSSYEGH quedo en train y DCPSDWSPYEGHCYK en
            # test con 75% de identidad y cobertura 1.00, con ambos tools de
            # acuerdo en esa cifra. MMseqs2 veia el parecido (valor E 3.5e-08);
            # el modo de agrupamiento no lo usaba.
            #
            # La componente conexa une por transitividad, que es la unica
            # garantia equivalente a "nada por encima del umbral cruza la
            # frontera". Coste medido: 8173 -> 7886 clusteres, el mayor pasa de
            # 379 a 454 secuencias (2.4% del total), sin cadenas degeneradas.
            "--cluster-mode",
            "1",
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


def identity_and_coverages(seq_a: str, seq_b: str, aligner) -> tuple[float, float, float]:
    """Identidad, cobertura sobre la corta y cobertura sobre la larga.

    Las dos coberturas existen porque el proyecto necesita las dos definiciones
    y confundirlas fue P-32:

    - sobre la **corta**, la que usa `cluster_builtin`: no castiga al peptido
      corto contenido en uno largo, que es el caso que al agrupar interesa unir;
    - sobre la **larga**, la que usa MMseqs2 con `--cov-mode 1` y por tanto la
      que define que parejas el agrupamiento se comprometio a separar. Es la que
      debe exigir la auditoria de fuga.

    Que la diferencia no es cosmetica, medido sobre las parejas reales que
    marcaba la auditoria:

        30 vs 5  (CCISG)   identidad 1.000  /corta 1.000  /larga 0.167
        29 vs 29 homologas identidad 1.000  /corta 1.000  /larga 1.000

    Exigir 0.5 a la cobertura sobre la corta no habria quitado ni uno de los
    falsos positivos, que era el arreglo que P-32 tenia escrito.
    """
    ident, cov_shorter = identity_and_coverage(seq_a, seq_b, aligner)
    shorter = min(len(seq_a), len(seq_b))
    longer = max(len(seq_a), len(seq_b))
    if longer == 0:
        return ident, cov_shorter, 0.0
    return ident, cov_shorter, cov_shorter * shorter / longer


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
    *,
    min_coverage: float,
) -> list[tuple[float, int | None]]:
    """Para cada consulta, su identidad maxima contra el conjunto sujeto.

    Una pareja solo cuenta si pasa identidad **y** cobertura, con la misma
    definicion de cobertura que exige el agrupamiento (P-32): sobre la secuencia
    larga, que es lo que mide MMseqs2 con `--cov-mode 1`. Sin esto la auditoria
    marcaba parejas que el agrupamiento nunca tuvo la obligacion de separar,
    como `CCISG` (cinco residuos) contra cualquier peptido rico en cisteina, con
    identidad 1.000 por aritmetica y cobertura 0.167 sobre la larga.

    `min_coverage` no tiene valor por defecto a proposito. Es el parametro del
    que depende el significado del p-valor de fuga, y P-29 ya enseno que un
    defecto se comporta como una decision del proyecto aunque nadie la haya
    tomado: quien llame a esta funcion declara que definicion de homologo audita.
    """
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
            ident, _cov_corta, cov_larga = identity_and_coverages(subject_seqs[idx], query, aligner)
            # La condicion es una O, y los dos terminos vigilan garantias
            # distintas (P-32):
            #
            # - cobertura sobre la larga: la garantia de MMseqs2 con
            #   `-c --cov-mode 1`, es decir que parejas el agrupamiento se
            #   comprometio a separar;
            # - contencion exacta: la garantia de D-20, la fusion por
            #   contencion. Si se auditara solo por cobertura, una subcadena
            #   corta en un split distinto (8 residuos dentro de 50 da cobertura
            #   0.16) dejaria de marcarse, y la auditoria quedaria ciega
            #   justamente ante el mecanismo mas joven y el que ya fallo dos
            #   veces (P-31, P-31b). La guarda tiene que ver las dos cosas.
            otra = subject_seqs[idx]
            contenida = otra in query or query in otra
            if cov_larga < min_coverage and not contenida:
                continue
            if ident > best_identity:
                best_identity = ident
                best_partner = idx
        out.append((best_identity, best_partner))
    return out


def audit_train_test_identity(
    df: pd.DataFrame,
    threshold: float,
    min_coverage: float,
    kmer_size: int = 3,
    prefilter: float = 0.05,
    max_test: int | None = None,
    seed: int = 42,
    null_model: bool = True,
    null_replicates: int = 20,
    query_split: str = "test",
    reference_split: str = "train",
) -> dict:
    """Identidad maxima de cada secuencia de `query_split` contra TODO `reference_split`.

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

    `query_split` esta parametrizado porque **test no es el unico split que
    importa**. `val` fija los hiperparametros y el umbral de decision, y ese
    umbral se aplica despues a test, asi que entra en cinco de las siete
    metricas (MCC incluida, que es la principal preregistrada). Homologia
    val-train infla la seleccion y contamina las metricas de test **aunque test
    este limpio**: la fuga no necesita tocar test para sesgar test.
    """
    aligner = _make_aligner()
    reference = df[df.split == reference_split]
    query = df[df.split == query_split]

    if reference.empty or query.empty:
        return {
            "error": (f"{reference_split} o {query_split} vacio; no se puede auditar"),
            "query_split": query_split,
            "reference_split": reference_split,
        }

    subject_seqs = reference.sequence.tolist()
    subject_kmers = [kmer_set(s, kmer_size) for s in subject_seqs]

    kmer_index: dict[str, set[int]] = defaultdict(set)
    for idx, kms in enumerate(subject_kmers):
        for km in kms:
            kmer_index[km].add(idx)

    query_rows = query
    sampled = False
    if max_test is not None and len(query) > max_test:
        query_rows = query.sample(n=max_test, random_state=seed)
        sampled = True

    test_seqs = query_rows.sequence.tolist()
    results = _max_identity_against(
        test_seqs,
        subject_seqs,
        subject_kmers,
        kmer_index,
        aligner,
        kmer_size,
        prefilter,
        min_coverage=min_coverage,
    )
    arr = np.asarray([r[0] for r in results])

    offenders = []
    for row, (ident, partner) in zip(query_rows.itertuples(), results):
        if ident >= threshold and len(offenders) < 25:
            offenders.append(
                {
                    "test_id": row.id,
                    "test_sequence": row.sequence,
                    "max_identity": round(ident, 4),
                    "closest_train_id": (
                        reference.iloc[partner].id if partner is not None else None
                    ),
                }
            )

    report = {
        "threshold": threshold,
        # Queda escrito en el informe porque el p-valor de fuga solo se puede
        # leer sabiendo que definicion de homologo se audito (P-32).
        "min_coverage": min_coverage,
        "query_split": query_split,
        "reference_split": reference_split,
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
                shuffled,
                subject_seqs,
                subject_kmers,
                kmer_index,
                aligner,
                kmer_size,
                prefilter,
                min_coverage=min_coverage,
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
        report["flagged_ids"] = [query_rows.iloc[i].id for i in range(len(arr)) if flagged[i]][:50]

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
        f"Identidad maxima de cada secuencia de {query_split} contra todo "
        f"{reference_split}. El numero de fuga es leakage.global_p_value y su "
        "exceso, NO fraction_over_threshold: la metrica de identidad tiene un "
        "ruido de fondo alto en peptidos cortos. Ver docs/09, P-04, P-12, P-15 "
        "y P-18."
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
    acfg = dict(scfg.get("audit", {}))

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

    # P-30: la fusion por contencion va DESPUES del agrupamiento y antes de
    # repartir, y se aplica cualquiera sea el metodo usado: el punto ciego que
    # tapa (fragmentos de 8-15 residuos) lo tienen los tres.
    mapping, contained_report = merge_contained_clusters(sequences, mapping)
    print(
        f"[split] fusion por contencion: {contained_report['n_contained_sequences']} "
        f"secuencias son subcadena de otra, "
        f"{contained_report['clusters_before']} -> {contained_report['clusters_after']} "
        f"clusteres ({contained_report['n_cluster_merges']} fusiones)"
    )

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

    # P-20: una division nueva invalida cualquier prediccion anterior. El
    # borrado vivia solo en el target del Makefile, asi que 'ampesm split' (y
    # por tanto 'ampesm all') dejaba predicciones rancias que evaluate publicaba
    # sin avisar. Se borra aqui, que es donde ocurre la invalidacion.
    predictions_dir = METRICS_DIR / "predictions"
    if predictions_dir.exists():
        stale = list(predictions_dir.glob("*.csv"))
        if stale:
            for f in stale:
                f.unlink()
            print(
                f"[split] {len(stale)} archivos de predicciones borrados: la division "
                f"nueva los invalida. Reentrena con 'make baselines' y 'make train'."
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
        # Valor provisional: solo mira la herramienta. Se recalcula al final de
        # build_splits, cuando las auditorias ya han corrido, porque un campo
        # llamado "publishable" que ignora el resultado de la auditoria de fuga
        # es peor que no tenerlo: en la primera corrida real dijo true con fuga
        # significativa en test Y en val (P-29).
        "publishable": used_method != "builtin",
        "publishable_note": (
            "false significa que no se cumple alguna condicion para publicar "
            "metricas con esta division: o se uso el respaldo propio en Python, "
            "que agrupa menos que MMseqs2 (docs/09, P-03), o alguna auditoria "
            "detecto fuga por homologia significativa (docs/09, P-29). El motivo "
            "exacto esta en publishable_reasons."
        ),
        "publishable_reasons": [],
        "contained_merge": contained_report,
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

    # P-26: la tasa de positivos no puede derivar entre splits.
    report["positive_rate_drift"] = _positive_rate_drift(report["splits"], acfg)

    # P-04: la verificacion directa de la premisa del proyecto.
    print("[split] auditando identidad maxima test-vs-train...")
    max_test = acfg.get("max_test", 500)
    report["train_test_identity_audit"] = audit_train_test_identity(
        df,
        threshold=identity,
        min_coverage=coverage,
        kmer_size=int(bcfg.get("kmer_size", 3)),
        prefilter=float(bcfg.get("kmer_prefilter", 0.05)),
        max_test=None if max_test is None else int(max_test),
        seed=params.seed,
        null_replicates=int(acfg.get("null_replicates", 20)),
    )

    # P-27: val fija los hiperparametros y el umbral, y ese umbral se aplica a
    # test. Homologia val-train sesga las metricas de test sin tocar test.
    if acfg.get("audit_val", True):
        print("[split] auditando identidad maxima val-vs-train...")
        # `max_val: null` en params significa "usa el mismo tope que test", no
        # "audita val entero": lo segundo reintroduciria el coste cuadratico de
        # P-19 por la puerta de atras.
        max_val = acfg.get("max_val") if acfg.get("max_val") is not None else max_test
        report["train_val_identity_audit"] = audit_train_test_identity(
            df,
            threshold=identity,
            min_coverage=coverage,
            kmer_size=int(bcfg.get("kmer_size", 3)),
            prefilter=float(bcfg.get("kmer_prefilter", 0.05)),
            max_test=None if max_val is None else int(max_val),
            seed=params.seed,
            null_replicates=int(acfg.get("null_replicates", 20)),
            query_split="val",
            reference_split="train",
        )

    # P-29: `publishable` se decide aqui y no arriba, porque hasta este punto no
    # se sabe el resultado de las auditorias de fuga.
    report["publishable_reasons"] = _publishable_reasons(report, used_method)
    report["publishable"] = not report["publishable_reasons"]

    SPLIT_REPORT_PATH.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"[split] {report['n_clusters']} clusteres ({report['cluster_size']['singletons']} de 1)")
    print(f"[split] {report['mixed_label_clusters']} clusteres con ambas clases")
    for split, info in report["splits"].items():
        print(
            f"[split] {split:5s}: {info['n']:6d} ({info['fraction_actual']:.1%}) "
            f"{info['n_positive']}+ / {info['n_negative']}-  "
            f"tasa positiva {info['positive_rate']:.3f}"
        )

    for key in ("train_test_identity_audit", "train_val_identity_audit"):
        audit = report.get(key)
        if audit is None or "error" in audit:
            continue
        label = f"{audit['query_split']}-vs-{audit['reference_split']}"
        print(
            f"[split] identidad {label}: {audit['fraction_over_threshold']:.1%} "
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
                extra = (
                    " val fija el umbral que se aplica a test, asi que esto sesga "
                    "las metricas de test aunque test este limpio."
                    if audit["query_split"] == "val"
                    else " Revisa el metodo de clustering antes de publicar."
                )
                print(f"[split] AVISO: fuga por homologia significativa en {label}.{extra}")

    drift = report["positive_rate_drift"]
    print(
        f"[split] deriva de tasa positiva: {drift['spread']:.3f} "
        f"(min {drift['min']:.3f} en {drift['argmin']}, "
        f"max {drift['max']:.3f} en {drift['argmax']})"
    )
    if drift["exceeds_tolerance"]:
        print(f"[split] AVISO: {drift['note']}")

    if not report["publishable"]:
        print("[split] AVISO: publishable=false. No publiques metricas con esta division:")
        for reason in report["publishable_reasons"]:
            print(f"[split]   - {reason}")

    print(f"[split] {SPLITS_PATH}")
    return df


def _length_ks_pvalue(pos: pd.DataFrame, neg: pd.DataFrame) -> float | None:
    """p-valor de Kolmogorov-Smirnov entre las longitudes de ambas clases."""
    if len(pos) < 2 or len(neg) < 2:
        return None
    from scipy.stats import ks_2samp

    return round(float(ks_2samp(pos.length.to_numpy(), neg.length.to_numpy()).pvalue), 6)


def _positive_rate_drift(splits_report: dict, audit_cfg) -> dict:
    """Comprueba que la tasa de positivos no derive entre splits.

    Por que importa (P-26): PR-AUC, F1 y precision dependen de la prevalencia,
    asi que con tasas distintas las metricas de los tres splits no son
    comparables entre si. Peor: el umbral se elige maximizando MCC en `val` y
    se aplica a `test`, y un umbral calibrado al 6% de positivos aplicado a un
    conjunto al 29% esta simplemente mal calibrado. Hay un caso conocido que
    pasa TODAS las demas validaciones con exactamente esas dos tasas.

    No aborta, avisa. A diferencia de un split degenerado (que no tiene
    arreglo), una deriva moderada sigue siendo utilizable si se reporta: lo que
    no se puede es no saberlo. La tolerancia sale de `split.audit.
    positive_rate_tolerance`.
    """
    rates = {k: v["positive_rate"] for k, v in splits_report.items() if v["n"] > 0}
    tolerance = float(audit_cfg.get("positive_rate_tolerance", 0.05))

    if len(rates) < 2:
        return {
            "rates": rates,
            "tolerance": tolerance,
            "error": "hacen falta al menos dos splits no vacios",
        }

    argmin = min(rates, key=lambda k: rates[k])
    argmax = max(rates, key=lambda k: rates[k])
    spread = rates[argmax] - rates[argmin]
    exceeds = spread > tolerance

    return {
        "rates": rates,
        "min": rates[argmin],
        "argmin": argmin,
        "max": rates[argmax],
        "argmax": argmax,
        "spread": round(spread, 4),
        "tolerance": tolerance,
        "exceeds_tolerance": bool(exceeds),
        "note": (
            f"la tasa de positivos varia {spread:.3f} entre splits "
            f"({argmin} {rates[argmin]:.3f} vs {argmax} {rates[argmax]:.3f}), "
            f"por encima de la tolerancia {tolerance}. PR-AUC, F1 y precision no "
            "son comparables entre splits, y el umbral elegido en val esta "
            "descalibrado para test. Ver docs/09, P-26."
            if exceeds
            else "dentro de tolerancia"
        ),
    }


def _publishable_reasons(report: dict, used_method: str) -> list[str]:
    """Motivos por los que NO se pueden publicar metricas con esta division.

    Lista vacia = publicable. Se calcula despues de las auditorias a proposito:
    la version anterior solo miraba la herramienta de clustering y en la primera
    corrida real declaro `publishable: true` con fuga significativa en test y en
    val a la vez (P-29). Un campo que se llama "publishable" y no mira el
    resultado de la auditoria de fuga es peor que no tenerlo, porque invita a
    saltarse justo la lectura que importa.
    """
    reasons: list[str] = []
    if used_method == "builtin":
        reasons.append(
            "se uso el respaldo propio en Python, que agrupa menos que MMseqs2 "
            "(docs/09, P-03, P-23)"
        )
    for key in ("train_test_identity_audit", "train_val_identity_audit"):
        audit = report.get(key)
        if not audit or "error" in audit or "leakage" not in audit:
            continue
        leak = audit["leakage"]
        if leak.get("significant"):
            reasons.append(
                f"fuga por homologia significativa en "
                f"{audit['query_split']}-vs-{audit['reference_split']} "
                f"(exceso {leak['excess']:+.1f} sobre {leak['n_expected_by_chance']} "
                f"esperadas por azar, p = {leak['global_p_value']:.3g})"
            )
    return reasons


def merge_contained_clusters(
    sequences: dict[str, str],
    mapping: dict[str, int],
    min_kmer: int = 4,
    max_examples: int = 10,
) -> tuple[dict[str, int], dict]:
    """Fusiona los clusteres de toda secuencia que sea subcadena exacta de otra.

    Motivo (P-30): el 8.5% del conjunto son fragmentos de otras secuencias -el
    peptido maduro dentro de su precursor, o el mismo registro con los bordes
    cortados de otra forma en DRAMP y APD3- y MMseqs2 no los agrupa. No es
    configuracion: su prefiltro necesita k-meros de 6-7 residuos y un peptido de
    8 tiene dos o tres, asi que desplazar una letra los destruye. Preguntandole
    pareja por pareja con todo permisivo solo encuentra 3 de 8 casos reales.

    Se eligio este criterio y no una pasada de alineamiento por pares (P-30,
    opcion 2) porque la contencion de subcadenas es **exacta e independiente de
    cualquier alineador**. La auditoria de fuga mide con Biopython; si el
    agrupamiento midiera con la misma regla, la auditoria dejaria de poder
    detectar los puntos ciegos de esa regla. Una guarda que comparte regla con
    lo que vigila no es una guarda. Ver D-20.

    El prefiltro de k-meros es una interseccion, no una union: una secuencia solo
    puede estar contenida en otra que contenga **todos** sus k-meros. Eso es una
    condicion necesaria exacta, asi que el prefiltro no descarta ningun caso
    verdadero; solo evita las 174 millones de comparaciones ciegas.
    """
    ids = list(sequences)
    seq_of = sequences
    # Indice de k-meros -> ids que lo contienen, solo de las secuencias que
    # pueden ser el contenedor (cualquiera: la longitud se comprueba despues).
    index: dict[str, set[str]] = defaultdict(set)
    for sid in ids:
        s = seq_of[sid]
        for i in range(len(s) - min_kmer + 1):
            index[s[i : i + min_kmer]].add(sid)

    # Union-find sobre CLUSTERES, no sobre secuencias (P-31b). La primera
    # version unia secuencias y luego reetiquetaba grupo a grupo leyendo del
    # mapeo original: si dos grupos de contencion tocaban el mismo cluster, el
    # segundo sobrescribia la fusion del primero. Resultado medido: RPPGFSPFR
    # quedo en test y ERPPGFSPFR en train, siendo subcadena limpia.
    parent: dict[int, int] = {cid: cid for cid in set(mapping.values())}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            # Al menor, para que el id del cluster resultante no dependa del
            # orden en que se encuentren las parejas.
            lo, hi = min(ra, rb), max(ra, rb)
            parent[hi] = lo

    contained: list[tuple[str, str]] = []
    for sid in ids:
        s = seq_of[sid]
        if len(s) < min_kmer:
            # Demasiado corta para el prefiltro. `min_length: 5` en params la
            # excluye del conjunto, pero la funcion no puede dar eso por hecho.
            continue
        kmers = [s[i : i + min_kmer] for i in range(len(s) - min_kmer + 1)]
        candidates = set(index[kmers[0]])
        for km in kmers[1:]:
            candidates &= index[km]
            if not candidates:
                break
        # Sin 'break' en este bucle: 803 de las 1 582 secuencias contenidas lo
        # estan en MAS DE UN contenedor (P-31). Cortar en el primero dejaba el
        # resto sin unir y hacia el resultado dependiente del orden de iteracion
        # de un set, que no es estable entre procesos. Es la trampa de P-14.
        containers = sorted(
            cid
            for cid in candidates
            if cid != sid and len(seq_of[cid]) > len(s) and s in seq_of[cid]
        )
        for cid in containers:
            contained.append((sid, cid))
            union(mapping[cid], mapping[sid])

    merged = {sid: find(cid) for sid, cid in mapping.items()}
    n_merges = len(set(mapping.values())) - len(set(merged.values()))

    report = {
        # Por secuencia, no por par: una secuencia dentro de tres contenedores
        # aporta tres pares y una sola secuencia redundante.
        "n_contained_sequences": len({a for a, _ in contained}),
        "n_containment_pairs": len(contained),
        "fraction_contained": round(len({a for a, _ in contained}) / max(len(ids), 1), 4),
        "n_cluster_merges": n_merges,
        "clusters_before": len(set(mapping.values())),
        "clusters_after": len(set(merged.values())),
        "examples": [
            {"contained": seq_of[a], "container": seq_of[b], "ids": [a, b]}
            for a, b in contained[:max_examples]
        ],
        "note": (
            "Secuencias que son subcadena exacta de otra mas larga. MMseqs2 no "
            "las agrupa a 8-15 residuos (P-30), asi que sin esta fusion acaban "
            "en splits distintos con identidad ~1.00. El criterio es exacto y no "
            "usa alineador, a proposito: ver D-20."
        ),
    }
    return merged, report
