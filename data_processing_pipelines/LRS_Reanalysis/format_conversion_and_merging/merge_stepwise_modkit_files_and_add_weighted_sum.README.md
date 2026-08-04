# Stepwise modkit Merge & Weighted-Sum Scoring Pipeline

This directory contains: `merge_stepwise_modkit_files_and_add_weighted_sum.py`, a script that merges a series of [modkit](https://github.com/nanoporetech/modkit)-derived bedRMod files produced at increasing model-confidence thresholds (0.85 → 0.99), records the highest confidence threshold at which each modification site is still called, and computes a composite **weighted-sum score** per site together with empirical, permutation-based p-values (Benjamini–Hochberg corrected). It processes a native dataset and, optionally, an in-vitro-transcribed (IVT) control dataset used as an alternative null.

> **Note:** This is a step-wise / notebook-style script (organised into `#%%` cells) rather than an argparse CLI. Input file names and thresholds are hardcoded (see *Inputs* below); edit the paths in the script before running.

## Scripts/Workflows

| File | Purpose |
|------|---------|
| `merge_stepwise_modkit_files_and_add_weighted_sum.py` | Loads bedRMod files across confidence thresholds 0.85–0.99, stamps each site with the maximal threshold at which it is called (`update_score`), computes a standardised weighted-sum `score`, and derives multiple empirical p-values via subsampled permutation null distributions with BH correction (`calculate_bedrmod_score_weighted_sum`). Writes scored native and IVT tables. |

### Key functions

| Function | Purpose |
|----------|---------|
| `load_bedrmod_polars` | Reads a bedRMod file into Polars, parses `#key=value` comments, auto-assigns column names by count, rescales frequency to 0–100, and aligns modification/chromosome names. |
| `shift_coordinates` | Adds a coordinate offset to start/end columns (used for adapter correction). |
| `align_modification_names_polars` / `align_chromosome_names_polars` | Standardise the `name` and `chrom` columns via `_MOD_MAP` / `_CHROM_MAP`. |
| `update_score` | Left-joins a source and target table on site keys and takes the maximal `score_nmod_thresholds`, propagating the highest confidence boundary. |
| `calculate_bedrmod_score_weighted_sum` | Filters by coverage, min-max standardises features (log₂ coverage, frequency/100, `score_nmod_thresholds`), computes the weighted sum, and builds empirical p-values from permutation nulls (row-shuffle, rank-based, and IVT-based) with BH correction. |
| `write_bedrmod_polars` | Writes a Polars DataFrame back to bedRMod format with metadata header. |

## Requirements & Setup

### Environment setup

To run the script you will need a Python 3.9+ environment with the required dependencies. An `environment.yml` file is provided for Conda users and a `requirements.txt` for pip users.

Core dependencies: `polars`, `pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-learn`, `matplotlib`, `seaborn`, `pyranges`, `upsetplot`, `joblib`, `pyarrow`, `tqdm`.

**Option A: Using Conda**

```bash
conda env create -f environment.yml
conda activate modkit-score
```

**Option B: Using Pip**

```bash
python -m venv .venv
source .venv/bin/activate
pip install polars pandas numpy scipy statsmodels scikit-learn \
            matplotlib seaborn pyranges upsetplot joblib pyarrow tqdm
```

## Input data

### Input files

| Input File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `0.85.bed` … `0.99.bed` | One bedRMod file per modkit confidence threshold, in 0.01 steps (15 files: `0.85`–`0.99`). Each is read for both the native and IVT tracks; rows on `rRNA`/`oligo` references and rows with `frequency == 0` are filtered out. | Tab-separated bedRMod (`.bed`) | `merge_stepwise_modkit_files_and_add_weighted_sum.py` |

### Inputs / configuration

Because this is a notebook-style script, key settings are hardcoded and must be edited in place:

- **Threshold range:** `range(85, 100, 1)` → files `0.85.bed` … `0.99.bed`.
- **Filtering:** references containing `rRNA` or `oligo` are dropped; `frequency > 0` is required.
- **Weighting:** the weighted sum uses three equally weighted (0.3333) standardised features — `log₂(coverage)` (min-max), `frequency / 100`, and `score_nmod_thresholds` (min-max). `n_mod` is available but excluded from the active weighting.
- **Permutation settings:** `n_permutations` (default 1000), `null_sample_size` (default 500), `min_coverage` (default 1); RNG seeded at 42 for reproducibility.
- **tRNA adapter correction:** an optional commented-out block subtracts the 24 nt 5′ adapter length and restricts sites to `0 ≤ chromStart ≤ 96`. Uncomment only for tRNA data.

### References/Accessory files

No external genomic references are required — modification and chromosome/reference mapping tables (`_MOD_MAP`, `_CHROM_MAP`, including human rRNA and tRNA reference names) are hardcoded in the script.

## Usage/Step-wise Execution

### Step 1: Prepare inputs

Place the 15 threshold bedRMod files (`0.85.bed` … `0.99.bed`) in the script's working directory, or edit the file-path template in the loading loop.

### Step 2: Run the script

```bash
python merge_stepwise_modkit_files_and_add_weighted_sum.py
```

The script (a) loads and filters all threshold files into native and IVT lists, (b) iteratively runs `update_score` to record the maximal confidence boundary per site, (c) calls `calculate_bedrmod_score_weighted_sum` to compute the weighted-sum score and permutation p-values, and (d) writes the scored output. It can also be run cell-by-cell (`#%%`) in an interactive environment (VS Code / Jupyter).

## Outputs

| Output File | Format | Script | Description |
|-------------|--------|--------|-------------|
| `combined_and_scored.bed` | Tab-separated bedRMod with `#key=value` metadata header | `merge_stepwise_modkit_files_and_add_weighted_sum.py` | Merged, threshold-annotated sites with standardised feature columns (`std_coverage`, `std_n_mod`, `std_score_nmod_thresholds`, `std_frequency`), the composite `score`, and empirical p-value columns (`pvalue_shuffle_axis_1`, `pvalue_adj_shuffle_axis_1`, `pvalue_shuffle_axis_0`, `pvalue_ranks`, and, when an IVT dataset is supplied, `pvalue_shuffle_axis_1_IVT` / `pvalue_shuffle_axis_0_IVT`). Native and IVT tables are written. |

> **Caveat:** As written, the script writes both the native and IVT tables to the same filename (`combined_and_scored.bed`), so the second write overwrites the first. Give the native and IVT outputs distinct filenames before running if you need to keep both.
