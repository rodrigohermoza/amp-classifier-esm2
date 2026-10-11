# Antimicrobial peptide classifier with ESM-2

[Español](README.md) · [Full documentation (Spanish)](docs/00-indice.md)

An antimicrobial peptide classifier (AMP vs non-AMP) that compares protein
language model embeddings (ESM-2) against physicochemical descriptors, using
sequence-identity cluster splits to prevent data leakage.

> **Status: code complete, not yet run.**
> The results tables are deliberately empty. They get filled by running the
> pipeline; the numbers will come from `results/metrics/`, not from an estimate.

## Why this project

An antimicrobial peptide is a short amino acid chain that kills bacteria or
fungi, usually by disrupting their membrane. They matter because that mechanism
is harder to evade than inhibiting a specific enzyme, and because lab testing is
slow and expensive: a classifier that ranks candidates before synthesis has
practical value.

In 2023, at Sup'Biotech (France), I took part in a project that designed AMPs
with a genetic algorithm ([team
repository](https://github.com/hgrandjean/Anti-Fungi-Peptide)). My contribution
was data acquisition from DRAMP and UniProt. This is my own version, where the
machine learning is my work end to end.

## What this project tries to get right

The goal is not a high ROC-AUC. In this task it is easy to reach 0.95+ while
having learned nothing about antimicrobial activity. What it tries to get right:

**1. A split with no data leakage.** AMPs come in families of near-identical
sequences: magainin 1 and 2 differ by two residues. A random split puts close
relatives in both train and test, and then the metric measures memorisation.
Here all sequences are clustered at 40% identity and **whole clusters** are
assigned to train/val/test. This is the core of the project and, in my view, the
most common methodological error in the AMP prediction literature.

**2. Honest baselines.** Antimicrobial activity depends mostly on two global
properties: net positive charge and a hydrophobic face. That makes 30
physicochemical numbers hard to beat. If the baseline wins, the README will say
so in the headline.

**3. Reported uncertainty.** Metrics with bootstrap confidence intervals,
resampling **clusters** rather than sequences, because related sequences are not
independent and treating them as if they were produces intervals that are too
narrow. Two models with overlapping intervals are not distinguishable.

**4. The main weakness stated out loud.** Below.

## The known weakness: the negatives

There is no database of peptides that were tested and found **not** to be
antimicrobial. What is used here, as throughout the literature, are UniProt
peptides without antimicrobial annotation. That does not mean "not
antimicrobial": it means **"nobody tested it, or nobody annotated it"**.

Two consequences:

- There are genuine false negatives among the negatives, so measurable
  performance is capped below true performance.
- Worse: the model can learn to distinguish **provenance** rather than
  **activity**. APD3 AMPs are secreted cationic peptides, many from amphibian
  skin; UniProt fragments are something else. A high AUC can mean "I know which
  database this came from" rather than "I know whether it kills bacteria".

Mitigations applied: matching the length distribution across classes, colouring
the UMAP by source as well as by class, and measuring per-source performance
separately in the error analysis. If accuracy on UniProt positives is much worse
than on APD3 positives, the model depends on provenance and the result has to be
read with that in mind.

## Data

The repository **does not contain raw data**: APD3 and DRAMP require citation and
restrict commercial use. It contains download scripts and a manifest with each
file's URL, date and `sha256`. `make data` rebuilds the dataset.

| Source | Role | Records | Size | License |
|---|---|---|---|---|
| [APD3](https://aps.unmc.edu/) | Positives | 3,306 | 147 KB | Academic, with citation |
| [DRAMP 3.0](https://dramp.cpu-bioinfor.org/) | Positives | — | 2.8 MB | Academic, non-commercial |
| [UniProt](https://www.uniprot.org/) KW-0929 | Positives | 1,585 | ~0.2 MB | CC BY 4.0 |
| UniProt without KW-0929, reviewed | Negatives | 11,982 | ~1.5 MB | CC BY 4.0 |

Sizes and counts measured on 2026-10-11. Total ≈ 5 MB.

Cleaning: standard 20 amino acids only, length 5–50, exact-sequence deduplication
and length matching. Details in [docs/03-datos-fuentes.md](docs/03-datos-fuentes.md).

**DBAASP is absent.** It was the fourth source and the only one with MIC values,
which would have enabled a potency regression. As of 2026-10-11 it returns 403 or
zero-byte responses on every route tried, so the regression is deferred. Logged
in [docs/09-problemas-encontrados.md](docs/09-problemas-encontrados.md), P-02.

## Method

```
download → cleaning → 40% clusters → ┬→ descriptors ┬→ models → evaluation
                                     └→ ESM-2       ┘
```

| Stage | What it does |
|---|---|
| **Split** | MMseqs2 at 40% identity; falls back to CD-HIT, then to a built-in greedy Python clustering with a k-mer prefilter. Whole clusters to train/val/test (70/15/15). |
| **Baselines** | Amino acid composition (20 features) and 10 physicochemical descriptors (net charge via Henderson-Hasselbalch, GRAVY, Eisenberg hydrophobic moment, pI, Boman index…) with logistic regression, random forest and XGBoost. |
| **Embeddings** | ESM-2 8M and 35M, per-residue mean pooling with explicit masking of padding and special tokens, cached on disk. Logistic regression and a small MLP on top. |
| **Evaluation** | ROC-AUC, PR-AUC, MCC and F1 with 95% CIs from cluster bootstrap. Threshold selected on validation. Paired comparison between the best embedding model and the best baseline. |
| **Interpretation** | UMAP by class, source and split; error analysis by category, source and length. |

The physicochemical scales are written out explicitly and cited in
[`physchem.py`](src/ampesm/features/physchem.py), so they can be audited without
leaving the repository and so the numbers do not shift when a library changes
version. Verified against published values: magainin-2 gives +3.0 net charge at
pH 7.4 (published +3 to +4) and LL-37 gives +5.9 (published +6).

## Results

Pending execution.

| Model | ROC-AUC | PR-AUC | MCC | F1 |
|---|---|---|---|---|
| Composition + physicochemical, XGBoost | — | — | — | — |
| ESM-2 35M + logistic regression | — | — | — | — |

### Pre-registration

Written before running anything, so this is not a post-hoc account of whatever
turned out well:

- Primary model: **ESM-2 35M + logistic regression**.
- Primary baseline: **composition + physicochemical with XGBoost**.
- Primary metric: **MCC** on test, with cluster bootstrap CIs.
- Test is evaluated **once**, at the end.
- My expectation: embeddings should win, but by a modest margin (perhaps
  0.02–0.04 ROC-AUC), and the intervals may overlap. If the margin comes out
  large, the first thing to suspect is leakage, not success.

## Reproducing

Requires [uv](https://docs.astral.sh/uv/). Python 3.12 (not 3.14: PyTorch does
not yet ship stable wheels for it).

```bash
git clone <url> && cd Peptidos-antimicrobianos
make setup
make all
```

Stage by stage:

```bash
make data        # download and clean  (~5 MB over the network)
make split       # clustering and split
make baselines   # descriptors + LR/RF/XGBoost
make embeddings  # ESM-2: ~31 MB for the 8M model, ~140 MB for 35M
make train       # classifiers on embeddings
make evaluate    # metrics with CIs + figures
make interpret   # UMAP and error analysis
```

`ampesm info` shows which artifacts exist and which are missing. Everything runs
on CPU on an Apple Silicon Mac. Clustering is optionally faster with
`brew install mmseqs2` (~30 MB).

## Documentation

The lab notebook in [`docs/`](docs/00-indice.md) is in Spanish. It includes the
decision log with rationale ([08](docs/08-decisiones.md)), the real problems hit
and how they were resolved ([09](docs/09-problemas-encontrados.md)), pending
improvements ranked by effort/benefit ([10](docs/10-mejoras-futuras.md)) and AI
assistance usage with tokens and time ([11](docs/11-agentes.md)).

If you want to attack the method, start with [04](docs/04-division-clusters.md)
and [09](docs/09-problemas-encontrados.md): that is where the interesting part is.

## Citations

- Wang G., Li X., Wang Z. *APD3: the antimicrobial peptide database.* NAR, 2016.
- Shi G. et al. *DRAMP 3.0.* NAR, 2022.
- The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase.*
- Lin Z. et al. *Evolutionary-scale prediction of atomic-level protein structure
  with a language model.* Science, 2023. (ESM-2)

## License

Code under [MIT](LICENSE). **The license does not cover the data**: each database
has its own terms, detailed in [docs/03](docs/03-datos-fuentes.md).
