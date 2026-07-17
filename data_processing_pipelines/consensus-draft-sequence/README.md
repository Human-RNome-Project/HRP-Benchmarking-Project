# Draft Reference Pipeline — Reproducible Subset

This folder contains a clean, reproducible subset of the Human RNome Project
draft-reference construction pipeline and plotting process (steps 1–6). It is organized into:

- `inputs/` – small input BEDs/TSVs, plus reference indices if provided.
- `scripts/` – all Python scripts, patched to use `paths.py`.
- `outputs/` – regenerated results from steps 1–5, plus the manually
  curated `consensus_draft_sequence.bed`.
- `scripts/panels/` – optional figure-panel scripts.
- `figures/` – generated panel figures.
- `run_pipeline.sh` – orchestration for steps 1–5.

## Environment setup

The pipeline requires `pandas`, `numpy`, and `matplotlib` (for the optional panel-figure scripts).

### Option 1: Conda

```bash
cd draft_reference_github
conda env create -f environment.yml
conda activate draft-ref-pipeline
```

### Option 2: pip

```bash
cd draft_reference_github
pip install -r requirements.txt
```

## Quick start

Once the environment is active:

```bash
cd draft_reference_github
./run_pipeline.sh
```

The pipeline executes steps 1–5:

1. MS harmonization (`scripts/harmonize_massspec.py`)
2. Parameter-range generation (`scripts/generate_parameter_ranges.py`,
   `scripts/generate_tRNA_parameter_ranges.py`)
3. Grid search (`scripts/run_rrna_grid_from_ranges.py`,
   `scripts/run_3way_nm_from_ranges.py`, `scripts/run_3way_m5C_from_ranges.py`,
   `scripts/run_3way_m6A_from_ranges.py`, `scripts/run_polya_grid_fast.py`,
   `scripts/run_trna_grid_from_ranges.py`)
4. Filtering (`scripts/run_filtering.py`, `scripts/run_tRNA_filtering.py`)
5. Tiered list generation (`scripts/create_tiered_mod_lists.py`,
   `scripts/create_tRNA_tiered_table.py`)

Individual panel plotting scripts (Step 6) are detailed further below.

---

## Input files

All files are in `inputs/`.

| File | Source / purpose |
|---|---|
| `MS_rRNA_tRNA.bed` | rmchrY MassSpec BED. Uses `mxA/mxC/mxG/mxU` for ambiguous calls. |
| `H.sapiens_rRNA_ref_mods.bed` | Curated human rRNA reference modifications. |
| `Illumina_combined_polyARNA_tRNA_rRNA.bed` | Combined Illumina short-read modification calls (rRNA, tRNA, polyA). |
| `ONT_polyARNA_rRNA_combined.filtered.bed` | Filtered ONT direct-RNA calls used for rRNA and polyA analyses. |
| `ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed` | Filtered ONT calls including tRNA, used for MS harmonization and tRNA analyses. |

### External References and Large Files

Certain panel plotting scripts (Step 6) require large genomic references that are 
not included in this repository to save space. These include:
- `gencode.v49.primary_assembly.annotation.gtf.gz`
- `GRCh38.primary_assembly.genome.fa.fai`
- `hs_rRNAs_NR_046235.fa.fai`

By default, the panel scripts will first look for these files inside the `inputs/` folder. 
If they are not found, they gracefully fall back to alternative local paths (e.g., `~/ref/`).

**Adjusting these paths**: If you need to change where the scripts look for these files, 
you can open the respective script in `scripts/panels/` and modify the paths defined 
in the `# ── Paths ───` section at the top of the file.

### Input BED format

Input BEDs follow the **bedRModv2** format with 13 tab-separated columns:

```text
chrom  chromStart  chromEnd  name  score  strand  thickStart  thickEnd  itemRgb  coverage  frequency  single_letter_code  mod_id
```

Header comment lines (`#fileformat=...`, `#modification_names=...`, etc.) are
preserved where relevant.

---

## Scripts

### Shared infrastructure

#### `scripts/paths.py`
Central path definitions. All other scripts import from this file so that
only one place needs editing if the folder layout changes.

| Constant | Points to |
|---|---|
| `ROOT` | `draft_reference_github/` |
| `INPUTS` | `inputs/` |
| `OUTPUTS` | `outputs/` |
| `THRESHOLDS` | `outputs/thresholds/` |
| `FILTERED` | `outputs/filtered_platform_beds/` |
| `TIERED` | `outputs/tiered_lists/` |
| `TIERED_TRNA` | `outputs/tiered_tRNA/` |

#### `scripts/utils.py`
Shared helper functions imported by most scripts.

- `parse_bed(path, platform_label)` – parses a bedRModv2 BED into a DataFrame
  with columns `chrom`, `start`, `end`, `name`, `score`, `strand`,
  `coverage`, `level`, `platform`. The 10th BED column (`frequency`) is
  renamed to `level`.
- `make_grid(series, n_points=7)` – generates up to 7 quantile-based
  thresholds, snapped to actual values in the input series.

---

### Step 1 — MS harmonization

#### `scripts/harmonize_massspec.py`
Resolves ambiguous MassSpec modification calls (`mxA/mxC/mxG/mxU` and legacy
`mA?/mC?/mG?/mU?`) to specific modifications by cross-referencing:

1. MS-self (other non-ambiguous mods at the same position)
2. Reference BED (`H.sapiens_rRNA_ref_mods.bed`) — rRNA only
3. Illumina combined BED — rRNA and tRNA
4. ONT tRNA-inclusive BED — rRNA and tRNA

Priority is MS-self > Reference (rRNA) > Illumina > ONT. Pseudouridine
resolution (`mxU`/`mU?` → `Y`/`pseU`/`psU`) is blocked; `Um` is allowed.

**Inputs**
- `inputs/MS_rRNA_tRNA.bed`
- `inputs/H.sapiens_rRNA_ref_mods.bed`
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed`

**Output**
- `outputs/MS_rRNA_tRNA_harmonized.bed`

**Output format**
bedRModv2 with 13 columns. Column 4 (`name`) is updated when a call is
resolved to a specific modification (e.g. `mxU` → `Um`). Unresolved ambiguous
calls are written as `mA?/mC?/mG?/mU?` in column 4 while the original raw
columns (12–13) preserve the input symbols.

---

### Step 2 — Parameter ranges

#### `scripts/generate_parameter_ranges.py`
Builds quantile-based parameter grids for rRNA and polyA. For each platform
and modification type it computes grids over `score`, `coverage`, and `level`
(`score` and `level` only for MassSpec, because MassSpec coverage is uniformly
0).

**Inputs**
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_combined.filtered.bed`
- `outputs/MS_rRNA_tRNA_harmonized.bed`

**Outputs**
- `outputs/thresholds/rRNA_grid_search_parameter_ranges.tsv`
- `outputs/thresholds/polyA_grid_search_parameter_ranges.tsv`
- `outputs/thresholds/grid_search_space_sizes.tsv`

**Output columns — `*_grid_search_parameter_ranges.tsv`**

| Column | Meaning |
|---|---|
| `biotype` | `rRNA` or `polyA` |
| `platform` | `Illumina`, `ONT`, or `MassSpec` |
| `mod_type` | Modification name (e.g. `Am`, `Y`, `m6A`) |
| `parameter` | `score`, `coverage`, or `level` |
| `n_values` | Number of non-missing values for this parameter |
| `min` / `max` / `mean` / `median` | Descriptive statistics |
| `n_grid_points` | Number of thresholds in the grid |
| `grid_thresholds` | Comma-separated list of grid thresholds |

**Output columns — `grid_search_space_sizes.tsv`**

| Column | Meaning |
|---|---|
| `biotype` | `rRNA` or `polyA` |
| `mod_type` | Modification name |
| `comparison` | Platform comparison (e.g. `Illumina × ONT`) |
| `p1_combos` / `p2_combos` / `p3_combos` | Number of threshold combinations per platform |
| `total_combos` | Total number of parameter combinations searched |
| `search_type` | `2-way` or `3-way` |

#### `scripts/generate_tRNA_parameter_ranges.py`
Same as above but restricted to tRNA rows (`hs_tRNA*` / `hs_mttRNA*`) and
using the tRNA-inclusive ONT BED.

**Inputs**
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed`
- `outputs/MS_rRNA_tRNA_harmonized.bed`

**Outputs**
- `outputs/thresholds/tRNA_grid_search_parameter_ranges.tsv`
- `outputs/thresholds/tRNA_grid_search_space_sizes.tsv`

---

### Step 3 — Grid search

Grid search reads the parameter-range TSVs, enumerates all threshold
combinations, and selects the combination that maximizes the Jaccard index
between the classified position sets of the compared platforms.

Classification rule for a single platform/mod/threshold combo:

```
score >= score_threshold AND coverage >= coverage_threshold AND level >= level_threshold
```

For MassSpec, coverage is ignored (set to 0).

#### `scripts/run_rrna_grid_from_ranges.py`
Pairwise (2-way) rRNA grid search across all platform pairs present for each
mod.

**Inputs**
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_combined.filtered.bed`
- `outputs/MS_rRNA_tRNA_harmonized.bed`
- `outputs/thresholds/rRNA_grid_search_parameter_ranges.tsv`

**Outputs**
- `outputs/thresholds/rRNA_grid_search_all_results.tsv` – every combo evaluated
- `outputs/thresholds/rRNA_grid_search_best.tsv` – best combo per mod/comparison

**Output columns — `rRNA_grid_search_best.tsv`**

| Column | Meaning |
|---|---|
| `mod_type` | Modification name |
| `comparison` | Pairwise comparison (e.g. `Illumina-ONT`) |
| `jaccard` | Best Jaccard index |
| `n_intersection` / `n_union` | Size of intersection / union at best thresholds |
| `n_plat1` / `n_plat2` | Number of positions passing thresholds in each platform |
| `<platform>_score` / `<platform>_coverage` / `<platform>_level` | Best thresholds per platform |

#### `scripts/run_3way_nm_from_ranges.py`
3-way grid search for rRNA Nm modifications: `Am`, `Cm`, `Gm`, `Um`.

**Output**
- `outputs/thresholds/rRNA_3way_Nm_best.tsv`

#### `scripts/run_3way_m5C_from_ranges.py`
3-way grid search for rRNA `m5C`.

**Output**
- `outputs/thresholds/rRNA_3way_m5C_best.tsv`

#### `scripts/run_3way_m6A_from_ranges.py`
3-way grid search for rRNA `m6A`.

**Output**
- `outputs/thresholds/rRNA_3way_m6A_best.tsv`

**Output columns — `rRNA_3way_*_best.tsv`**

| Column | Meaning |
|---|---|
| `mod_type` | Modification name |
| `best_jaccard` | Best 3-way Jaccard index |
| `n_intersection` / `n_union` | 3-way intersection / union sizes |
| `n_Illumina` / `n_ONT` / `n_MassSpec` | Positions passing thresholds per platform |
| `<platform>_score` / `<platform>_coverage` / `<platform>_level` | Best thresholds per platform |

#### `scripts/run_polya_grid_fast.py`
Pairwise polyA grid search for mods present in both Illumina and ONT. Uses
numpy boolean arrays for speed.

**Inputs**
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_combined.filtered.bed`
- `outputs/thresholds/polyA_grid_search_parameter_ranges.tsv`

**Outputs**
- `outputs/thresholds/polyA_grid_search_all_results.tsv`
- `outputs/thresholds/polyA_grid_search_best.tsv`

**Output columns** are the same as for `rRNA_grid_search_best.tsv`.

#### `scripts/run_trna_grid_from_ranges.py`
tRNA grid search. Runs 3-way search for mods present in all three platforms,
and pairwise search for mods present in exactly two platforms.

**Inputs**
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed`
- `outputs/MS_rRNA_tRNA_harmonized.bed`
- `outputs/thresholds/tRNA_grid_search_parameter_ranges.tsv`

**Outputs**
- `outputs/thresholds/tRNA_grid_search_all_results.tsv`
- `outputs/thresholds/tRNA_grid_search_best.tsv` (2-way results)
- `outputs/thresholds/tRNA_3way_best.tsv`

**Output columns** follow the same conventions as the rRNA grid-search files.

---

### Step 4 — Filtering

#### `scripts/run_filtering.py`
Applies the optimal grid-search thresholds to the combined BEDs to produce
filtered subset BEDs for rRNA and polyA. Nm mods (`Am`, `Cm`, `Gm`, `Um`) use
the 3-way thresholds; `m5C` and `m6A` use their 3-way thresholds; all other
rRNA mods use the best pairwise thresholds. PolyA mods use the best
Illumina-ONT thresholds.

**Inputs**
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_combined.filtered.bed`
- `outputs/MS_rRNA_tRNA_harmonized.bed`
- `outputs/thresholds/rRNA_grid_search_best.tsv`
- `outputs/thresholds/rRNA_3way_Nm_best.tsv`
- `outputs/thresholds/rRNA_3way_m5C_best.tsv`
- `outputs/thresholds/rRNA_3way_m6A_best.tsv`
- `outputs/thresholds/polyA_grid_search_best.tsv`

**Outputs**
- `outputs/filtered_platform_beds/rRNA_illumina_filtered.bed`
- `outputs/filtered_platform_beds/rRNA_ont_filtered.bed`
- `outputs/filtered_platform_beds/rRNA_massspec_filtered.bed`
- `outputs/filtered_platform_beds/polyA_illumina_filtered.bed`
- `outputs/filtered_platform_beds/polyA_ont_filtered.bed`
- `outputs/filtered_platform_beds/filtering_summary.tsv`

**Output format**
bedRModv2 with all 13 original columns preserved. Rows are a strict subset of
the input BEDs; only the rRNA MassSpec raw-name columns are normalized from
`mx*` to `m?` to match the legacy filtered BED.

**Output columns — `filtering_summary.tsv`**

| Column | Meaning |
|---|---|
| `biotype` | `rRNA` or `polyA` |
| `platform` | Platform label |
| `n_in` | Rows before filtering |
| `n_out` | Rows after filtering |
| `pct` | Percentage retained |

#### `scripts/run_tRNA_filtering.py`
Same idea for tRNA. Uses 3-way thresholds when available and 2-way thresholds
otherwise.

**Inputs**
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed`
- `outputs/MS_rRNA_tRNA_harmonized.bed`
- `outputs/thresholds/tRNA_grid_search_best.tsv`
- `outputs/thresholds/tRNA_3way_best.tsv`

**Outputs**
- `outputs/filtered_platform_beds/tRNA_illumina_filtered.bed`
- `outputs/filtered_platform_beds/tRNA_ont_filtered.bed`
- `outputs/filtered_platform_beds/tRNA_massspec_filtered.bed`
- `outputs/filtered_platform_beds/tRNA_filtering_summary.tsv`

---

### Step 5 — Tiered lists

#### `scripts/create_tiered_mod_lists.py`
Builds tiered rRNA and polyA site lists.

- **Tier1** = union of grid-search-optimized comparison intersections.
  - rRNA Nm mods use 3-way intersections.
  - rRNA `m5C` and `m6A` use 3-way intersections.
  - rRNA other mods (e.g. `Y`) use all 2-way intersections.
  - polyA uses Illumina-ONT 2-way intersections.
- **Tier2** = sites present in ≥2 raw combined BEDs but not in tier1.

**Inputs**
- `inputs/H.sapiens_rRNA_ref_mods.bed`
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_combined.filtered.bed`
- `outputs/MS_rRNA_tRNA_harmonized.bed`
- `outputs/thresholds/rRNA_grid_search_best.tsv`
- `outputs/thresholds/rRNA_3way_Nm_best.tsv`
- `outputs/thresholds/rRNA_3way_m5C_best.tsv`
- `outputs/thresholds/rRNA_3way_m6A_best.tsv`
- `outputs/thresholds/polyA_grid_search_best.tsv`

**Outputs**
- `outputs/tiered_lists/tiered_rRNA_only.tsv`
- `outputs/tiered_lists/tiered_rRNA_only.bed` (BED6 version)
- `outputs/tiered_lists/tiered_rRNA_shared_mod_types.tsv`
- `outputs/tiered_lists/tiered_polyA.tsv`
- `outputs/tiered_lists/tiered_polyA_<mod>.tsv` per polyA mod

**Output columns — `tiered_rRNA_only.tsv` / `tiered_polyA.tsv`**

| Column | Meaning |
|---|---|
| `chr` | Chromosome |
| `start` / `end` | 0-based BED coordinates |
| `name` | Modification name |
| `tier` | `tier1` or `tier2` |
| `strand` | Strand |
| `in_ref` | `TRUE` if position is in `H.sapiens_rRNA_ref_mods.bed` |
| `ref_mod_type` | Reference modification at that position, or `NA` |

**Output columns — `tiered_rRNA_shared_mod_types.tsv`**
Same columns as above, restricted to rRNA sites whose mod type is also found
in tRNA MassSpec data.

#### `scripts/create_tRNA_tiered_table.py`
Builds tiered tRNA site lists using the same tier1/tier2 logic as the rRNA
script.

**Inputs**
- `inputs/Illumina_combined_polyARNA_tRNA_rRNA.bed`
- `inputs/ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed`
- `outputs/MS_rRNA_tRNA_harmonized.bed`
- `outputs/thresholds/tRNA_grid_search_best.tsv`
- `outputs/thresholds/tRNA_3way_best.tsv`

**Outputs**
- `outputs/tiered_tRNA/tiered_tRNA.tsv`
- `outputs/tiered_tRNA/tiered_tRNA_counts.tsv`

**Output columns — `tiered_tRNA.tsv`**

| Column | Meaning |
|---|---|
| `chr` / `start` / `end` | Coordinates |
| `name` | Modification name |
| `tier` | `tier1` or `tier2` |
| `strand` | Strand |
| `platforms` | Comma-separated supporting platforms |
| `comparison` | Grid-search comparison or `raw_overlap` |
| `jaccard` | Jaccard index of the comparison, or `NA` |

**Output columns — `tiered_tRNA_counts.tsv`**

| Column | Meaning |
|---|---|
| `mod` | Modification name |
| `tier1` / `tier2` / `total` | Site counts per tier |

---

## Step 6 — Panel figures

Optional plotting scripts live in `scripts/panels/` and reproduce individual
panels from the draft-reference figures using the pipeline outputs. All panel
scripts resolve paths relative to the repository root and output their
PDF and PNG results to dedicated subdirectories inside `figures/`.

### `scripts/panels/plot_panel_a_polyA_manhattan_density.py`

Reproduces Panel A: a Manhattan-style plot of polyA modification density in
1 Mb bins across the primary chromosomes.

The script derives the bin counts directly from the tiered polyA site list, so
it does not depend on any pre-computed count file.

**Inputs**
- `outputs/tiered_lists/tiered_polyA.tsv`
- `inputs/GRCh38.primary_assembly.genome.fa.fai` (preferred) or
  `~/ref/GRCh38.primary_assembly.genome.fa.fai` (fallback)

**Outputs**
- `figures/panel_a_polyA_manhattan_density/panel_a_polyA_manhattan_density.pdf`
- `figures/panel_a_polyA_manhattan_density/panel_a_polyA_manhattan_density.png`
- `figures/panel_a_polyA_manhattan_density/panel_a_polyA_1Mb_bin_counts.tsv`

### `scripts/panels/plot_panel_b_chr1_160M_region_zoom.py`

Reproduces Panel B: a zoomed locus view of the chr1:160.45-160.85 Mb region
(SLAM/CD2 gene family) showing protein-coding gene models and polyA
modification sites colored by modification type.

**Inputs**
- `outputs/tiered_lists/tiered_polyA.tsv`
- `inputs/gencode.v49.primary_assembly.annotation.gtf.gz` (preferred) or
  `~/ref/gencode.v49.primary_assembly.annotation.gtf.gz` (fallback)

**Outputs**
- `figures/panel_b_chr1_160M_region_zoom/panel_b_chr1_160M_region_zoom.pdf`
- `figures/panel_b_chr1_160M_region_zoom/panel_b_chr1_160M_region_zoom.png`
- `figures/panel_b_chr1_160M_region_zoom/panel_b_region_modifications.tsv`
- `figures/panel_b_chr1_160M_region_zoom/panel_b_region_genes.tsv`

### `scripts/panels/plot_panel_c_polyA_per_mb_vs_genes.py`

Reproduces Panel C: a scatter plot of polyA modification density vs protein-coding 
gene density per chromosome. Includes a linear fit line and chromosome labels.

**Inputs**
- `outputs/tiered_lists/tiered_polyA.tsv`
- `inputs/gencode.v49.primary_assembly.annotation.gtf.gz` (preferred fallback)
- `inputs/GRCh38.primary_assembly.genome.fa.fai` (preferred fallback)

**Outputs**
- `figures/panel_c_polyA_per_mb_vs_genes/panel_c_polyA_perMb_vs_genes_perMb.pdf`
- `figures/panel_c_polyA_per_mb_vs_genes/panel_c_polyA_perMb_vs_genes_perMb.png`

### `scripts/panels/plot_panel_d_polyA_sites_per_gene_histogram.py`

Reproduces Panel D: a histogram of polyA modification sites per protein-coding gene.

**Inputs**
- `outputs/tiered_lists/tiered_polyA.tsv`
- `inputs/gencode.v49.primary_assembly.annotation.gtf.gz` (preferred fallback)

**Outputs**
- `figures/panel_d_polyA_sites_per_gene_histogram/panel_d_polyA_sites_per_genes_histogram.pdf`
- `figures/panel_d_polyA_sites_per_gene_histogram/panel_d_polyA_sites_per_genes_histogram.png`
- `figures/panel_d_polyA_sites_per_gene_histogram/panel_d_top100_protein_coding_genes_by_polyA_sites.tsv`

### `scripts/panels/plot_panel_e_observed_vs_expected_mod_load.py`

Reproduces Panel E: a plot comparing observed vs expected polyA modification load 
per gene based on a negative-binomial GLM regressing observed counts on gene length 
and transcript abundance (TPM).

**Inputs**
- `outputs/tiered_lists/tiered_polyA.tsv`
- `inputs/gencode.v49.primary_assembly.annotation.gtf.gz` (preferred fallback)
- `inputs/OUT.gene_tpm.tsv` (preferred fallback)

**Outputs**
- `figures/panel_e_observed_vs_expected_mod_load/panel_e_observed_vs_expected_mod_load.pdf`
- `figures/panel_e_observed_vs_expected_mod_load/panel_e_observed_vs_expected_mod_load.png`

### `scripts/panels/plot_panel_f_composite_metagene.py`

Reproduces Panel F: a composite metagene plot showing modification density across 
transcript regions (5' UTR, CDS, 3' UTR) for multiple modification types.

**Inputs**
- `outputs/tiered_lists/tiered_polyA.tsv`
- `inputs/gencode.v49.primary_assembly.annotation.gtf.gz` (preferred fallback)

**Outputs**
- `figures/panel_f_composite_metagene/panel_f_composite_metagene_overlay.pdf` (and `.png`)
- `figures/panel_f_composite_metagene/panel_f_composite_metagene_grid.pdf` (and `.png`)
- `figures/panel_f_composite_metagene/panel_f_composite_metagene_horizontal_3col.pdf` (and `.png`)

### `scripts/panels/plot_panel_g_rRNA_region_PTC.py`

Reproduces Panel G: a regional plot of rRNA modifications across human rRNA sequences.

**Inputs**
- `outputs/tiered_lists/tiered_rRNA_only.tsv`
- `inputs/hs_rRNAs_NR_046235.fa.fai` (preferred fallback)

**Outputs**
- `figures/panel_g_rRNA_region_PTC/panel_g_rRNA_region_PTC.pdf`
- `figures/panel_g_rRNA_region_PTC/panel_g_rRNA_region_PTC.png`

### `scripts/panels/plot_panel_h_tiered_rRNA_only_9o3v.py`

Reproduces Panel H: a 3D structural visualization of human rRNA modifications mapped 
onto the 9o3v ribosome structure using PyMOL.

*Note: This script uses PyMOL. Ensure PyMOL is installed in your environment or run it 
using `pymol -cq scripts/panels/plot_panel_h_tiered_rRNA_only_9o3v.py` to generate the image.*

**Inputs**
- `outputs/tiered_lists/tiered_rRNA_only.tsv`
- `inputs/9o3v.cif` (preferred fallback)
- `inputs/rcsb_pdb_9O3V.fasta` (preferred fallback)

**Outputs**
- `figures/panel_h_tiered_rRNA_only_9o3v/panel_h_tiered_rRNA_only_9o3v.png`

---

## Outputs and verification

All generated outputs were compared against the original files in
`draft_reference/` / `draft_reference_github_backup/`. See
`outputs/VERIFICATION_REPORT.md` for the detailed comparison.

Summary of matches:

- **Exact match:** harmonized MS BED, all grid-search best-threshold TSVs,
  all filtered BEDs (data rows), and tRNA tiered tables.
- **Data identical / order differs:** rRNA and polyA tiered lists (sorting
  order changed).
- **Manually placed:** `outputs/consensus_draft_sequence.bed` is not generated
  by this subset; it is kept here for convenience and matches the original.
