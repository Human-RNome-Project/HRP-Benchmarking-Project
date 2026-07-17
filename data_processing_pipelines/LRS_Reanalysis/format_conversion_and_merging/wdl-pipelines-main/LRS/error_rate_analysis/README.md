# Error Rate Analysis WDL Pipelines

This directory contains WDL workflows and helper scripts for calculating direct-RNA sequencing (DRS) error rates at two resolutions:

- **Native mismatch rate** — per-position error rates computed from pysamstats pileups, with masking and optional Modkit / IVT stratification.
- **Read error rate** — per-read error rates sampled from an aligned BAM, aggregated genome-wide.

A standalone helper script for coverage-binned error-rate analysis is also included.

## Directory contents

| Path | Type | Description |
|---|---|---|
| `pipeline_native_mismatch_rate.wdl` | WDL | Per-sample, per-chromosome native mismatch rate workflow. |
| `pipeline_read_error_rates.wdl` | WDL | Per-read error rate workflow with chromosome scatter/merge. |
| `calculate_native_mismatch_rates_fast.py` | Script | Chromosome-level and merge-level native mismatch calculation. Invoked inside `NativeMismatchRate` tasks. |
| `calculate_read_error_rates_parallel.py` | Script | Chromosome-level and merge-level per-read error calculation. Invoked inside `ReadErrorRates` tasks. |
| `calc_error_rate_inputs/` | Inputs | One JSON per sample for `NativeMismatchRate`. |
| `read_error_calc_inputs/` | Inputs | One JSON per sample for `ReadErrorRates`. |
| `shared_masks/` | Accessory | Mask BEDs used by `NativeMismatchRate` (analysis regions, variants, junctions, transcript strand, GIAB-derived beds). |
| `HG001_GRCh38_1_22_v4.2.1_benchmark.bed` | Accessory | GIAB HG001 high-confidence benchmark regions for GRCh38 chromosomes 1–22. |
| `ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed` | Accessory | Modkit modification pileup BED used for modification stratification (optional). |

## Workflows

| WDL file | Workflow namespace | Purpose |
|---|---|---|
| `pipeline_native_mismatch_rate.wdl` | `NativeMismatchRate` | Per-sample, per-chromosome native mismatch rate calculation from pysamstats TSVs. Scatters over chromosomes, then merges into sample-level summaries. |
| `pipeline_read_error_rates.wdl` | `ReadErrorRates` | Per-read error rate calculation from an aligned BAM. Scatters over chromosomes using reservoir sampling, then merges into genome-wide summaries. |

## Helper scripts

| Script | Purpose |
|---|---|
| `calculate_native_mismatch_rates_fast.py` | Chromosome-level and merge-level native mismatch rate calculation. Invoked inside `NativeMismatchRate` tasks. |
| `calculate_read_error_rates_parallel.py` | Chromosome-level and merge-level per-read error rate calculation. Invoked inside `ReadErrorRates` tasks. |

## Accessory files (`shared_masks/`)

| File | Description |
|---|---|
| `analysis_regions.bed` | Final high-confidence analysis regions after subtracting variants and splice junctions from GIAB regions. |
| `giab_filtered.bed` | GIAB benchmark regions filtered to autosomes and chrX. |
| `giab_minus_variants.bed` | GIAB regions with known variant positions removed. |
| `junctions.bed` | Splice-junction flanking regions excluded from error-rate calculations. |
| `transcript_strand.bed` | Transcript strand annotations used for strand-aware counting. |
| `variants.bed` | Known variant positions (e.g., GIAB HG001) excluded from error-rate calculations. |

## Requirements

### Workflow engine

The production runs used **JAWS** on Perlmutter, but the workflows can also be run with **Cromwell** or **miniwdl** locally.

Examples:
- JAWS: `jaws submit --no-cache <wdl> <inputs.json> perlmutter --tag <tag>`
- Cromwell: `java -jar cromwell.jar run <wdl> -i <inputs.json>`
- miniwdl: `miniwdl run <wdl> --input <inputs.json>`

### Docker

All WDL tasks run inside Docker containers. The executor must be able to pull the images listed below. If you are running offline, pull them beforehand.

## Docker images

```bash
# Native mismatch rate calculation
# Used by ChromTask and MergeTask in pipeline_native_mismatch_rate.wdl
docker pull kandarprj/drs-native-mismatch-1.5

# Per-read error rate calculation
# Used by ChromTask and MergeTask in pipeline_read_error_rates.wdl
docker pull kandarprj/drs-read-error-rates-1.0
```

> **Note:** These images are referenced by tag. If reproducibility is critical, consider rebuilding the images from the Dockerfiles and pinning by digest.

## Input JSON preparation

### Native mismatch rate (`pipeline_native_mismatch_rate.wdl`)

Ready-to-use input JSONs are provided under `calc_error_rate_inputs/`, one per sample (e.g. `HRP_A_001_native_polyA_RNA_001_inputs.json`).

Each JSON contains:
- `pysamstats_tsvs_per_sample`: per-chromosome pysamstats TSVs for the sample
- `sample_labels`: the sample name
- `chroms` / `chroms_str`: autosomes chr1–chr22
- `analysis_regions_bed`, `variants_bed`, `junctions_bed`, `transcript_strand_bed`: point to `shared_masks/`
- `giab_bed`: path to the GIAB benchmark BED (`HG001_GRCh38_1_22_v4.2.1_benchmark.bed`)
- `modkit_bed`: path to a Modkit pileup BED (`ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed`, optional; used for modification stratification)
- Parameter defaults: `min_coverage=10`, `junction_bp=50`, `min_mod_fraction=0.5`, `min_mod_coverage=20`, `cpus=4`, `memory_gb=40`

> **Note:** The provided JSONs currently reference these BEDs via absolute NERSC paths. If you want to use the copies now included in this directory, update the paths in the JSONs, e.g.:
> ```bash
> sed -i 's|/pscratch/sd/k/kandarpj/HRP_benchmarking_project/error_rate_calc/HG001_GRCh38_1_22_v4.2.1_benchmark.bed|HG001_GRCh38_1_22_v4.2.1_benchmark.bed|g' calc_error_rate_inputs/*.json
> sed -i 's|/global/cfs/cdirs/m5243/final_bedRmods/ont/ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed|ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed|g' calc_error_rate_inputs/*.json
> ```

Example for a single sample:

```bash
jaws submit \
  --no-cache \
  pipeline_native_mismatch_rate.wdl \
  calc_error_rate_inputs/HRP_A_001_native_polyA_RNA_001_inputs.json \
  perlmutter \
  --tag "HRP_A_001_native_mismatch"
```

To run the aggregate IVT/no-Supplement control:

```bash
jaws submit \
  --no-cache \
  pipeline_native_mismatch_rate.wdl \
  calc_error_rate_inputs/IVT_noSupp_inputs.json \
  perlmutter \
  --tag "IVT_noSupp_mismatch"
```

### Read error rates (`pipeline_read_error_rates.wdl`)

Ready-to-use input JSONs are provided under `read_error_calc_inputs/`, one per sample (e.g. `HRP_A_001_native_polyA_RNA_001_read_error_rates_inputs.json`).

Each JSON contains:
- `bam` / `bai`: aligned/sorted BAM and index
- `sample_label`: sample name
- `chroms`: chr1–chr22, chrX, chrY, chrM
- `sample_n_per_chrom`: 40000 reads per chromosome
- `min_mapq=20`, `min_length=200`, `seed=1000`, `cpus=4`

Example:

```bash
jaws submit \
  --no-cache \
  pipeline_read_error_rates.wdl \
  read_error_calc_inputs/HRP_A_001_native_polyA_RNA_001_read_error_rates_inputs.json \
  perlmutter \
  --tag "HRP_A_001_read_errors"
```

### Coverage-binned error rate (local helper)

```bash
python calculate_coverage_error_rate.py \
  --pysamstats   /path/to/native_variation_per_position.tsv \
  --outdir       results/coverage_error \
  --sample-label GM12878_native \
  --giab-bed     /path/to/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
  --chroms       chr19
```

## Outputs

### `pipeline_native_mismatch_rate.wdl`

Per sample:

| Output | File | Description |
|---|---|---|
| `mismatch_summaries` | `<sample>_mismatch_rate_summary.tsv` | Aggregated mismatch/error rates across regions and masks. |
| `masking_reports` | `<sample>_masking_report.tsv` | Counts and fractions of positions removed by each mask. |
| `per_position_tsvs` | `<sample>_per_position_mismatch.tsv` | Concatenated per-position mismatch table. |
| `substitution_spectra` | `<sample>_substitution_spectrum.tsv` | 12-class substitution spectrum (optional). |
| `modkit_stratified` | `<sample>_modkit_stratified_rates.tsv` | Error rates stratified by Modkit modification calls (optional). |
| `native_ivt_deltas` | `<sample>_native_ivt_delta.tsv` | Position-wise native-vs-IVT error-rate deltas (optional). |
| `chrom_position_tsvs` | `<sample>_<chrom>_per_position_mismatch.tsv` | Per-chromosome per-position tables. |

### `pipeline_read_error_rates.wdl`

| Output | File | Description |
|---|---|---|
| `chrom_per_read_tsvs` | `<sample>_<chrom>_per_read.tsv` | Per-read metrics per chromosome. |
| `chrom_counts_npzs` | `<sample>_<chrom>_counts.npz` | Raw count arrays per chromosome for merging. |
| `per_read_tsv` | `<sample>_per_read_error_rates.tsv` | Concatenated per-read error-rate table. |
| `summary_tsv` | `<sample>_read_error_rate_summary.tsv` | Genome-wide summary statistics. |
| `spectrum_tsv` | `<sample>_read_substitution_spectrum.tsv` | 12-class substitution spectrum. |
| `positional_bias_tsv` | `<sample>_read_positional_bias.tsv` | Error rates binned by read position. |
