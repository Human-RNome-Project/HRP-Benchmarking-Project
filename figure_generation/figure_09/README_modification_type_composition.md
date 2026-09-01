# Modification type composition (original cell-based script)

This is the README for the **original** interactive analysis script,
`modification_type_composition.py` — the cell-based (`#%%`) notebook for the
**RNome** project: a multi-lab comparison of RNA modification calls from Oxford
Nanopore (ONT) direct-RNA sequencing of GM12878 poly(A) RNA, rRNA and tRNA. It
loads per-site modification BED tables, annotates each modified site with a
genomic feature, filters to a high-confidence call set (via an IVT control test
and a composite confidence score), and produces the publication figures
(modification-by-feature heatmap, poly(A) tail-length violins, and per-dataset
PCA plots).

The script is written as a cell-based notebook (`#%%` cells) meant to be stepped
through interactively (VS Code / Jupyter), **not** run end-to-end as a single
`python` invocation. (A refactored command-line version with an argument parser,
manifest-driven inputs and no `dmode` dependency also exists — see
`README.md` — but this README documents the original notebook as written.)

## What it does

The pipeline runs in numbered steps, each an `#%%` cell:

| Step | Purpose |
|------|---------|
| 0 | Imports (pandas + polars used side by side; polars for the large per-site joins). |
| 1 | Central **RNome/Nature plotting theme** — fixed modification/platform colours, Nature figure sizing, `set_rnome_theme()`, panel labels, vector-only saving. |
| 2 | Load the GENCODE **GTF** and build collapsed gene-body intervals (`dmode`). |
| 3 | Load the combined ONT per-site modification BED; drop chrY; split into **genome** vs spiked-in **rRNA** contigs. |
| 4 | Overlap every site with gene bodies (**PyRanges** left join). |
| 5 | **Feature-classification decision tree** → label each site as 5'UTR / CDS / 3'UTR / intronic / rRNA species (18S/28S/5S/5.8S) / lncRNA-etc / intergenic. |
| 6 | Load the **15 individual lab datasets** (re-basecalled per-site BEDs). |
| 7 | Compute a composite `weighted_sum` **confidence score** per site (coverage + frequency + threshold score, each min-max scaled). |
| 8 | **Modification-name harmonisation** — collapse synonyms / ChEBI IDs / substitutions onto canonical short names (`_MOD_MAP`). |
| 9 | Load the **native-vs-IVT** Fisher table (IVT = in-vitro-transcribed, unmodified false-positive baseline). |
| 10 | **Native-vs-IVT enrichment test** — vectorised one-sided hypergeometric survival function + Benjamini-Hochberg FDR correction. |
| 11 | **High-confidence site filter** (significant vs IVT *or* very high confidence, with coverage/frequency floors). |
| 12–13 | Summarise modification × feature distribution and draw **Figure 1: heatmap**. |
| 14–15 | Load per-lab poly(A) **tail-length** parquet tables and draw the **tail-length violin plot**. |
| 16 | Re-imports + `pca_plot()` helper (samples × positions matrix → PCA). |
| 17–19 | **Figure 1S PCA** plots for the polyA-RNA, rRNA, and tRNA datasets. |

### High-confidence site definition (Step 11)

A site is kept if **either**:

- it is significantly enriched vs IVT: `padj ≤ 0.05` **and** `ivt_frequency / frequency ≤ 0.5`, **or**
- it has very high composite confidence: `weighted_sum ≥ 0.9`

**and**, in both cases, has solid support: `coverage ≥ 30` **and** `frequency ≥ 3%`.

## Enviroment setup

Python 3.10+ is required (the script uses `dict[str, str]` / `list[...] | None`
syntax). The following packages are needed:

```
pandas
polars
numpy
matplotlib
seaborn
pyranges
scipy
statsmodels
scikit-learn
upsetplot
tqdm
pyarrow
joblib
dmode          # in-house project helper package (metagene / gene-body utilities)
```

The **in-house `dmode` package** is required for
`dmode.metagene_plot.prepare_gene_body_coverage()`; it is not available on
PyPI/conda and must be installed from the project's own source. (The refactored
CLI version removes this dependency.)

Fonts: the theme targets **Helvetica/Arial**; without them matplotlib falls back
to DejaVu Sans.

### Using Conda

There is no packaged `environment.yml` for the original notebook (chiefly because
`dmode` is not on any public channel). A minimal environment can be created
manually, after which `dmode` is installed separately from its own source:

```bash
conda create -n rnome python=3.10 pandas polars numpy matplotlib seaborn \
    pyranges scipy statsmodels scikit-learn upsetplot tqdm pyarrow joblib
conda activate rnome
# then install the in-house dmode package from the project source, e.g.:
# pip install /path/to/dmode
```

## Input data

All paths are **hard-coded** in the script and point at the project's storage (a
mix of `/global/cfs/cdirs/m5243/...` on NERSC and `/home/stefan/Synology/...`).
You will need to edit these for your environment. Expected inputs:

- **GTF annotation** — `gencode.v49.primary_assembly.annotation.gtf`
- **Combined ONT BED** — `ONT_polyARNA_rRNA_combined.filtered.bed`
  (columns: chrom, chromStart, chromEnd, name, score, strand, thickStart,
  thickEnd, itemRgb, coverage, frequency, single_letter_code, mod_id)
- **15 per-lab genome BEDs** — one `*.merged.aligned.sorted.genome.bed` per lab
  (schema has split per-call counts: n_mod, count_canonical, count_other_mod,
  count_delete, count_fail, count_diff, count_nocall)
- **Native-vs-IVT Fisher BED** — `polyA.native_vs_ivt_fisher.bed`
- **Poly(A) tail-length parquet tables** — one per lab plus a merged file
- **rRNA / tRNA transcriptome-aligned BEDs** — for the PCA steps
  (note: the rRNA/tRNA schema uses `count_modified` rather than `n_mod`)

Datasets are identified by codes like `HRP_A_007_1` (lab / donor / replicate).

## Usage

Run interactively cell by cell (do **not** run the whole file at once).
Recommended order:

1. Run Steps 0–1 to set up imports and the plotting theme.
2. Run Steps 2–10 to build and test the annotated table
   (slow: has per-row annotation loops and a large IVT join).
3. **Checkpoint:** after Step 10 the table is written to
   `combined_labs_with_ivt_fisher.csv`. On later runs you can skip Steps 6–10 and
   read that CSV back in (Step 11) instead of recomputing.
4. Run Steps 11–13 for the heatmap, 14–15 for the tail-length violins, and
   16–19 for the PCA figures.

## Output

- `combined_labs_with_ivt_fisher.csv` — checkpoint of the fully annotated + tested
  table (write once, then re-read to resume from Step 11 without recomputing).
- `single_datasets_and_merged_taillengths.parquet` — stacked tail-length checkpoint.
- **Figures** (all vector, for publication):
  - `Figure1_modification_types_on_features_heatmap.pdf`
  - `taillengths_violin.svg`
  - `Figure1S_polyARNA_PCA_mod_frequency.pdf`
  - `Figure1S_rRNA_PCA_mod_frequency.pdf`
  - `Figure1S_tRNA_PCA_mod_frequency.pdf`

## Importent Notes

These are carried over from comments in the script and are worth checking before
relying on results:

- **Font sizes** — the docstring cites the Nature 5–7 pt guideline, but the values
  actually applied in `set_rnome_theme()` are ~20–28 pt. Lower them if you need
  true Nature sizing.
- **Lab-code mismatches** — a couple of per-lab comments don't match the assigned
  `lab` string (e.g. Lab12 comment says `A016` but `lab="HRP_A_011_1"`; Lab15
  comment says `A_015` but `lab="HRP_A_003_1"`). Double-check the mapping.
- **Feature tie-break** — in the classification loop, `candidate_gene[...][-1]`
  takes the **last** matching GTF record when several annotations overlap. This is
  deliberate but arbitrary.
- **tRNA PCA nesting (Step 19)** — unlike Steps 17/18, `bed_files_list` and
  `names_list` are wrapped in an extra `[ ]`, so `pca_plot` receives a
  single-element list-of-lists. If the tRNA PCA errors or shows a single point,
  that nesting is the likely cause. Left as-is.
- **chrY** is dropped in several places (GM12878 is female; chrY calls would be
  spurious).
- **IVT** = in-vitro-transcribed control (no real modifications → false-positive
  baseline). `score == 0` marks native (non-IVT) rows in the annotated tables.
