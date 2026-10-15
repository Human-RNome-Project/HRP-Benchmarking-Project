# wf-ont-bedrmod-native-vs-ivt

Snakemake workflow for identifying RNA modifications from Oxford Nanopore
direct-RNA sequencing data by comparing native and in-vitro transcribed (IVT)
samples.

For every site called in a native sample, the workflow tests whether the
modification rate is significantly higher than in a matched IVT control
(unmodified by design), applies a set of effect-size and frequency filters,
removes Nanopore cross-reactivity artefacts, and produces diagnostic plots and
supplementary tables at every stage.

The workflow is driven by a single config file and runs over an arbitrary
number of native/IVT pairs ("biotypes"); the shipped config defines three:
`polyA`, `rRNA`, and `tRNA`.

---

## Repository layout

```text
config/config.yml        Configuration: input pairs, thresholds, colors
workflow/Snakefile       All rules
workflow/rules/          Helper rules (bgzip/tabix)
workflow/envs/           Conda environment definitions
workflow/profile/        Snakemake profile for the NIH Biowulf SLURM cluster
scripts/                 Python scripts invoked by the rules
snakemake.sh             Cluster launcher
raw/                     Input BEDRMod files (not tracked)
analysis/                All generated output (not tracked)
```

`analysis/`, `raw/`, `logs/`, and `.snakemake/` are listed in `.gitignore` —
the repository tracks code and configuration only.

---

## Environment setup

The only hard requirement is Snakemake; all rule-level dependencies are
declared per rule and installed by Snakemake itself via `--use-conda`.

- `workflow/envs/base.yml` — Python ≥ 3.13 with `pandas`, `numpy`,
  `matplotlib`, `seaborn`, `statsmodels`, `openpyxl`. Used by every Python
  rule.
- `workflow/envs/samtools.yml` — `samtools` ≥ 1.21, providing `bgzip` and
  `tabix`.

---

## Quick start

Place the input BEDRMod files under `raw/`, point `config/config.yml` at them,
then run:

```bash
./snakemake.sh
```

`snakemake.sh` loads the `snakemake/7.32.3` module and runs the workflow
through the Biowulf profile in `workflow/profile/`, submitting each rule as a
SLURM job. To run locally instead:

```bash
snakemake -s workflow/Snakefile --use-conda --cores 8
```

An alternative config can be supplied with `--configfile`; the `Snakefile`
reads that flag directly and falls back to `config/config.yml`.

---

## Configuration

All settings live in `config/config.yml`.

| Key | Meaning |
|---|---|
| `NATIVE_VS_IVT` | Table of comparisons. Columns: `name`, `native_path`, `ivt_path`, `min_native_coverage`, `min_ivt_coverage`. One row per biotype; `name` becomes the wildcard used throughout the output file names. |
| `ILLUMINA_VALIDATION.bed` | Orthogonal Illumina call set used to benchmark the Nanopore predictions. |
| `SANITIZE.drop_chroms` | Chromosomes excluded during sanitization (e.g. `chrY` for female samples). |
| `ANALYSIS_DIR` / `CONFIG_DIR` / `SCRIPTS_DIR` | Directory locations. |
| `PLOT_COLUMNS` | BEDRMod columns to produce distribution and scatter plots for (`coverage`, `frequency`, `n_mod`, `score`). |
| `VOLCANO.y_cap` | Clip −log10(padj) at this value in the volcano plots. |
| `SCORE_DENSITY.sig_padj` | padj threshold defining the "significant" group in the score-density plots. |
| `FINAL_FILTER` | Per-biotype site filters: `min_log2fc`, `max_padj`, `min_ivt_absent_score`, `min_native_freq`. |
| `PROXIMITY_FILTER` | Per-biotype `max_proximity` (nt) for cross-reactivity deduplication. |
| `MODIFICATION_COLORS` | Fixed color per modification name, applied consistently across every figure. |

Adding a biotype means adding one row to `NATIVE_VS_IVT.data` plus matching
entries under `FINAL_FILTER` and `PROXIMITY_FILTER`; every rule expands over
those names automatically.

---

## Input files

Inputs are **BEDRMod** files produced upstream by the modification caller, one
native and one IVT file per biotype. Column names are read from the header
line, so extra caller-specific columns are tolerated; the workflow requires at
least:

```text
chrom  chromStart  chromEnd  name  score  strand  thickStart  thickEnd  color
coverage  frequency  n_mod  count_canonical  ...
```

- `name` — modification name (`m6A`, `Y`, `Am`, …)
- `score` — caller confidence, 0–1
- `coverage` / `n_mod` / `count_canonical` — read counts at the position
- `frequency` — modified fraction, in percent

Two header conventions are accepted: the bedRModv2 tab-separated column line
(`#chrom<TAB>chromStart…`) and the Python-list form
(`#['chrom', 'chromStart', …]`) used internally. The sanitization step converts
the former to the latter, so every downstream script can auto-detect columns
from a single header format.

---

## Workflow

### Step 1 — Sanitization

#### `rule sanitize_bedrmod` → `scripts/sanitize_bedrmod.py`

Normalizes the header of each raw BEDRMod file to the `#[col_list]` format,
strips a spurious `score_nmod_thresholds` field that the raw files inject into
the metadata lines, removes repeated headers, and drops unwanted chromosomes.
Accepts plain or gzipped input; always writes plain text. All other comment
lines and every data row are passed through unchanged.

**Inputs**
- `NATIVE_VS_IVT.native_path` / `ivt_path` for each biotype

**Outputs**
- `analysis/sanitized/{name}_native.bed`
- `analysis/sanitized/{name}_ivt.bed`

Every subsequent rule reads the sanitized files, never the raw ones.

---

### Step 2 — Input QC

#### `rule plot_bedrmod_col_distribution`
→ `scripts/plot_bedrmod_col_distribution.py`

Distribution of a single BEDRMod column, run for every sanitized file ×
every column in `PLOT_COLUMNS`.

The PDF is a 2×2 figure — cumulative fraction (ECDF), non-cumulative fraction,
cumulative count, non-cumulative count — with one curve per modification plus
an overall curve. For position-level columns such as `coverage` the overall
curve is deduplicated by position. For the coverage distribution the rule also
draws a vertical line at the configured `min_{side}_coverage` and clamps the
x-axis at 1.

**Inputs**
- `analysis/sanitized/{stem}.bed`

**Outputs**
- `analysis/{stem}.{column}_dist.tsv` — percentile table, overall and per
  modification
- `analysis/{stem}.{column}_dist.pdf`

#### `rule native_vs_ivt_scatter` → `scripts/plot_native_vs_ivt_scatter.py`

Hexbin scatter of one column, IVT (x) vs native (y), over sites present in both
files (inner join on `chrom`/`chromStart`/`chromEnd`/`strand`). Page 1 pools
all modifications; one further page per modification type.

**Inputs**
- `analysis/sanitized/{name}_native.bed`
- `analysis/sanitized/{name}_ivt.bed`

**Outputs**
- `analysis/{name}.native_vs_ivt.{column}_scatter.pdf`
- `analysis/{name}.native_vs_ivt.{column}_scatter.tsv` — Pearson and Spearman
  correlations, overall and per modification

---

### Step 3 — Native vs IVT statistical test

#### `rule native_vs_ivt_fisher` → `scripts/compare_native_vs_ivt.py`

The core step. Each native site is annotated with the matching IVT counts and a
one-sided Fisher's exact p-value (native > IVT), computed from `n_mod` vs
`count_canonical` in both conditions. The test is vectorized through the
hypergeometric distribution, and Benjamini–Hochberg FDR correction is applied
jointly over all testable sites.

Sites are retained regardless of testability:

- native coverage ≤ `min_native_coverage` → kept, all IVT and statistics fields
  set to `NA`
- IVT coverage ≤ `min_ivt_coverage`, or the site absent from the IVT file →
  treated as IVT-absent, `NA` statistics

Only the columns needed for the test are loaded into memory; the output is then
streamed line-by-line from the native file.

**Inputs**
- `analysis/sanitized/{name}_native.bed`
- `analysis/sanitized/{name}_ivt.bed`

**Output**
- `analysis/{name}.native_vs_ivt_fisher.bed.gz` (+ `.tbi`)

**Output format**
All native BEDRMod columns, unchanged, plus five appended columns:

| Column | Meaning |
|---|---|
| `ivt_frequency` | Modified fraction at the same position in IVT, or `NA` |
| `ivt_coverage` | IVT read depth, or `NA` |
| `ivt_n_mod` | IVT modified-read count, or `NA` |
| `pvalue` | One-sided Fisher's exact p-value (native > IVT), or `NA` |
| `padj` | Benjamini–Hochberg adjusted p-value, or `NA` |

---

### Step 4 — Statistical QC

All four rules read `analysis/{name}.native_vs_ivt_fisher.bed.gz`.

#### `rule pvalue_distribution` → `scripts/plot_pvalue_distribution.py`

Multi-page PDF, one page per modification plus an overall page. Each page has
three panels: raw p-value histogram with the uniform expectation line, BH
adjusted p-value histogram, and a native frequency vs −log10(padj) hexbin.

**Outputs**
- `analysis/{name}.native_vs_ivt_pvalue_dist.pdf`
- `analysis/{name}.native_vs_ivt_pvalue_dist.tsv` — total, testable, and
  IVT-absent counts, plus sites passing each significance threshold, overall
  and per modification

#### `rule volcano_plot` → `scripts/plot_volcano.py`

log2(native frequency / IVT frequency) vs −log10(padj), testable sites only.
One overall page plus one hexbin page per modification.

**Outputs**
- `analysis/{name}.native_vs_ivt_volcano.pdf`
- `analysis/{name}.native_vs_ivt_volcano.tsv` — site counts at common padj
  thresholds

#### `rule score_density` → `scripts/plot_score_density.py`

KDE of the BEDRMod `score` across four groups: native significant
(padj < `SCORE_DENSITY.sig_padj`), IVT significant (same sites, IVT score),
IVT all, and native IVT-absent. Page 1 overlays all modifications on one axis
(native-significant solid, IVT-all dashed, x restricted to score > 0.5);
subsequent pages give the overall and per-modification breakdown.

**Inputs**
- `analysis/{name}.native_vs_ivt_fisher.bed.gz`
- `analysis/sanitized/{name}_ivt.bed`

**Output**
- `analysis/{name}.native_vs_ivt_score_density.pdf`

#### `rule threshold_sweep` → `scripts/plot_threshold_sweep.py`

Surviving site counts as a function of the padj threshold (testable sites) and
the score threshold (IVT-absent sites), for all biotypes at once. One row per
dataset × two columns; one curve per modification plus a dashed total on a log
y-axis, with vertical lines at the chosen operating thresholds. This is the
plot used to justify the values in `FINAL_FILTER`.

**Output**
- `analysis/threshold_sweep.pdf`

---

### Step 5 — Site filtering

#### `rule filter_sites` → `scripts/filter_sites.py`

Applies the per-biotype `FINAL_FILTER` thresholds. A site is kept when

```text
testable (padj != NA):
    ( padj <= max_padj AND log2(native_freq / ivt_freq) > min_log2fc
      OR score > min_ivt_absent_score )
    AND native_freq > min_native_freq

IVT-absent (padj == NA):
    score > min_ivt_absent_score AND native_freq > min_native_freq
```

Header lines are preserved and the input is streamed in chunks, so memory use
is constant regardless of file size.

**Input**
- `analysis/{name}.native_vs_ivt_fisher.bed.gz`

**Outputs**
- `analysis/{name}.native_vs_ivt_fisher.filtered.bed.gz` (+ `.tbi`)
- `analysis/{name}.native_vs_ivt_fisher.filter_counts.tsv`

**Output columns — `*.filter_counts.tsv`**

The first line is a comment recording the thresholds used
(`#max_padj=…  min_log2fc=…  min_ivt_absent_score=…  min_native_freq=…`).
Column labels embed the actual threshold values, so the table is
self-documenting; one row per modification.

| Column | Meaning |
|---|---|
| `modification` | Modification name |
| `total` | Sites before filtering |
| `ivt_testable` / `ivt_absent` | Split by testability |
| `padj_le_<P>` | Sites passing the padj cutoff |
| `padj_le_<P>_and_log2fc_gt_<F>` | …and the fold-change cutoff |
| `padj_le_<P>_and_log2fc_gt_<F>_freq_gt_<R>` | …and the frequency cutoff |
| `padj_le_<P>_log2fc_gt_<F>_or_score_gt_<S>` | Testable-site OR condition |
| `absent_and_score_gt_<S>` | IVT-absent sites passing the score cutoff |
| `*_freq_gt_<R>` variants | The two conditions above, with the frequency cutoff added |
| `after_filter` | Sites retained |

#### `rule proximity_filter` → `scripts/proximity_filter.py`

[Currently not used as part of the RNOME project]
Removes Nanopore cross-reactivity artefacts: a site is discarded when another
site on the same (`chrom`, `strand`) at a **different** `chromStart` within
`max_proximity` nt carries a strictly higher score. Different modification
types called at the same position are never compared against each other. Input
order is preserved; the input must be coordinate-sorted, which the `bgzip_bed`
rule guarantees.

**Input**
- `analysis/{name}.native_vs_ivt_fisher.filtered.bed.gz`

**Outputs**
- `analysis/{name}.native_vs_ivt_fisher.filtered.deduped.bed.gz` (+ `.tbi`)
- `analysis/{name}.native_vs_ivt_fisher.dedup_counts.tsv` — `modification`,
  `before_dedup`, `after_dedup`, with a `#max_proximity=…` comment line

#### `rule swap_scores` → `scripts/swap_scores.py`

Browser-visualisation copy of the filtered set. Every row and coordinate is
preserved; two columns are rewritten:

- `score` ← −log10(padj). `NA` padj becomes 0; padj = 0 becomes a fixed
  constant above the ~323 ceiling of −log10(padj), so the strongest sites get
  the top code consistently across biotypes.
- `color` ← the original score (0–1) encoded as red intensity, `R,0,0` with
  R = score × 255.

The original score is retained verbatim in an appended `composite_score`
column, so no information is lost. In a genome browser the track is then shaded
by statistical significance while the caller score remains readable as redness.

**Input**
- `analysis/{name}.native_vs_ivt_fisher.filtered.bed.gz`

**Output**
- `analysis/{name}.native_vs_ivt_fisher.filtered.scores_swapped.bed.gz`
  (+ `.tbi`)

---

### Step 6 — Post-filter summaries

#### `rule plot_filtered_summary` → `scripts/plot_filtered_summary.py`

Three-page cross-biotype summary, generated twice — once for the filtered set
and once for the deduplicated set:

1. Stacked bar — site counts per modification and biotype, split into
   testable and IVT-absent panels
2. Violin — log2(native_freq / IVT_freq) per modification, faceted by biotype
   (testable sites only)
3. Heatmap — modifications × biotypes, colour = log10(count), annotated with
   raw counts

**Inputs**
- `analysis/{name}.native_vs_ivt_fisher.{variant}.bed.gz` for all biotypes,
  `variant` ∈ {`filtered`, `filtered.deduped`}

**Outputs**
- `analysis/filtered_summary.pdf`
- `analysis/filtered.deduped_summary.pdf`

#### `rule plot_filtered_scatter` → `scripts/plot_filtered_scatter.py`

Four-page per-biotype view of the filtered sites; pages 1–3 are grids with one
panel per modification, each in its configured colour:

1. Native frequency vs IVT frequency (testable sites, log-log)
2. Score vs native frequency (all sites)
3. Score vs log2(native / IVT frequency) (testable sites)
4. Native frequency distribution per modification (violin / strip)

Threshold lines are drawn at the biotype's `min_ivt_absent_score` and
`min_log2fc`.

**Input**
- `analysis/{name}.native_vs_ivt_fisher.filtered.bed.gz`

**Output**
- `analysis/{name}.native_vs_ivt_fisher.filtered_scatter.pdf`

#### `rule ivt_site_counts`

Total IVT sites per modification, counted directly from the sanitized IVT BED
with `awk`. Provides the denominator for the supplementary tables.

**Input**
- `analysis/sanitized/{name}_ivt.bed`

**Output**
- `analysis/{name}.ivt_counts.tsv` — `modification`, `ivt_total`

#### `rule supplementary_filter_table` / `..._deduped`
→ `scripts/build_filter_table.py`

Assembles the per-biotype count TSVs into a supplementary Excel workbook, one
sheet per biotype, with columns derived from the TSV metadata header and a
`Total` row appended at the bottom of each sheet. The `_deduped` variant also
folds in the proximity-filter counts.

**Inputs**
- `analysis/{name}.native_vs_ivt_fisher.filter_counts.tsv`
- `analysis/{name}.ivt_counts.tsv`
- `analysis/{name}.native_vs_ivt_fisher.dedup_counts.tsv` (deduped variant only)

**Outputs**
- `analysis/supplementary_filter_counts.xlsx`
- `analysis/supplementary_filter_counts_deduped.xlsx`

---

### Step 7 — Orthogonal validation

#### `rule illumina_validation` → `scripts/plot_illumina_validation.py`

Benchmarks the Nanopore calls against an independent Illumina call set. One
subplot per modification type; curves show Illumina overlap (%) as a function
of site rank, for the unfiltered set (sites above the native coverage
threshold) and the filtered set (IVT-testable sites only). A site counts as
overlapping when `chrom`, `chromStart`, `strand`, and modification name all
match.

**Inputs**
- `analysis/{name}.native_vs_ivt_fisher.bed.gz`
- `analysis/{name}.native_vs_ivt_fisher.filtered.bed.gz`
- `ILLUMINA_VALIDATION.bed`

**Output**
- `analysis/{name}.illumina_validation.pdf`

---

### Helper rules — `workflow/rules/samtools.smk`

#### `rule bgzip_bed`

Coordinate-sorts a `.bed` file (`sort -k1,1 -k2,2n`, `#` header lines kept on
top) and compresses it with BGZF, producing a tabix-indexable `.bed.gz`. Every
BED-producing rule writes plain text and lets this rule take over, which is
also what guarantees the sort order `proximity_filter` depends on.

#### `rule tabix_bed`

Indexes a BGZF-compressed BED, producing `.bed.gz.tbi`.

---

## Cluster execution

`workflow/profile/` is a Snakemake profile for the NIH Biowulf SLURM cluster.
It derives all `sbatch` options from each rule's `threads` and `resources`
(`mem_mb`, `runtime`, optionally `disk_mb`, `gpu`, `gpu_model`), infers the
partition from the requested resources, and polls job status from dashboard
data rather than querying SLURM directly. See
[workflow/profile/Readme.md](workflow/profile/Readme.md) for details.

Every rule in this workflow declares `mem_mb` and `runtime`; the heaviest is
`native_vs_ivt_fisher` at 64 GB / 8 h.

---

## Methods

### Identification of RNA modification sites

To distinguish genuine RNA modifications from sequencing noise, modification
calls derived from native direct-RNA sequencing were compared against calls from
an in-vitro transcribed (IVT) control, which is unmodified by design. For each
genomic position present in the native BEDRMod file, a one-sided Fisher's exact
test was performed to assess whether the modification rate in the native sample
was significantly greater than in the IVT control. Specifically, for a given
site, a 2×2 contingency table was constructed from the counts of modified and
canonical reads in each condition, and the p-value was computed as P(X ≥
n<sub>mod,nat</sub>), where X follows a Hypergeometric(N, K, n) distribution
with N = total reads across both conditions, K = total modified reads, and n =
total native reads. Only sites with a minimum coverage of 30 reads in both the
native and IVT samples were considered testable; sites below this threshold were
retained in the output but assigned NA for all statistical fields. Sites absent
from the IVT sample were similarly assigned NA p-values. Multiple testing
correction was applied to all testable sites jointly using the
Benjamini–Hochberg false discovery rate (FDR) procedure.

### Site filtering

Testable sites were retained when they satisfied an adjusted p-value cutoff
together with a minimum log2 fold change of the native over the IVT
modification frequency, or, alternatively, a high caller score; sites absent
from the IVT control were retained on the caller score alone. Both classes were
additionally required to exceed a minimum native modification frequency.
Thresholds were set per biotype (polyA, rRNA, tRNA) from the threshold sweep,
with the more stringent fold-change and frequency cutoffs applied to tRNA.

