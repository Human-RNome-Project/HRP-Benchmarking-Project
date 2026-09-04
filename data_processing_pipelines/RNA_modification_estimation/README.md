# Error Rate and RNA Modification Estimation Panel Plotting

This directory contains the plotting scripts and documentation for reproducing Panels B, C, D, E, F, and G of the RNA modification estimation section (`Figure 13`). 

## Directory Structure
- `scripts/panels/`: Contains the plotting scripts that generate the panels.
- `scripts/generation/`: Upstream scripts that generate the data required by the panels.
- `inputs/`: Intended directory for the required input data files (not tracked in Git due to size).
- `figures/`: Output directory where the generated plots will be saved.

## Environment Setup

To run the panel scripts, you will need a Python environment with the required dependencies (`pandas`, `numpy`, `matplotlib`, `scipy`). We have provided both an `environment.yml` file for Conda users and a `requirements.txt` file for pip users.

### Option A: Using Conda (Recommended)
1. Ensure you have Conda (or Miniconda/Mamba) installed.
2. Create the environment by running:
   ```bash
   conda env create -f environment.yml
   ```
3. Activate the environment:
   ```bash
   conda activate hrp_panels
   ```

### Option B: Using Pip
1. Create and activate a standard virtual environment:
   ```bash
   python3 -m venv hrp_env
   source hrp_env/bin/activate
   ```
2. Install the required packages:
   ```bash
   pip install -r requirements.txt
   ```

## Data Download Instructions

The preprocessed data files required to run the panel plotting scripts are hosted on the DOE Data Explorer.
- **Data Repository:** [DOE Data Explorer (doi:10.25585/DOE-HRP/3377574)](https://doi.org/10.25585/DOE-HRP/3377574)
- **Directory Path:** `processed_data/RNA_modificaiton_estimation`

To download the required input files and set up your environment:

1. Create the `inputs/` folder in this directory:
   ```bash
   mkdir -p inputs
   ```

2. Download all data files from the `processed_data/RNA_modificaiton_estimation` folder on the DOE Data Explorer and place them inside the `inputs/` folder. The required file structure is described in detail below.

> [!NOTE]
> **Upstream Input Generation (`pysamstats`):**
> Prior to the estimation pipeline, initial position-level count tables (`*_variation.tsv` / `*_variation.tsv.gz`) were generated upstream from aligned BAM files and a reference FASTA:
> 1. **Alignment Filtering (`samtools`):** Primary alignments were extracted per chromosome or target sequence from the alignment BAM using `samtools view -F 2304` to remove non-primary (secondary and supplementary) alignments.
> 2. **Count Table Generation (`pysamstats`):** Count tables were created using `pysamstats` with the following parameters:
>    - `--type variation`: outputs position-by-position base count and variation statistics.
>    - `--fasta <ref_genome>`: specifies the reference genome/transcriptome FASTA file.
>    - `--min-mapq 20`: filters out reads with mapping quality below 20.
>    - `--pad`: pads unread or missing positions across the target sequence.
>    - `-D 500000`: sets the max coverage depth limit to 500,000.
>    - `--chromosome <chrom>`: specifies the chromosome or target sequence being processed.
> 
> The resulting per-chromosome `*_variation.tsv` files serve as the input count tables processed by the upstream generation scripts (`calculate_native_mismatch_rates_memeff_fast.py`, `calculate_coverage_error_rate.py`, and `rdd_genomewide.py`) located in `scripts/generation/`.


## Panels and Required Inputs

The scripts in `scripts/panels/` require the following pre-calculated data files to be present in the `inputs/` folder.

### Panel B: Per-Read Total Error Rate (Separated)
**Script:** `scripts/panels/panel_b_per_read_error.py`
**Output:** `figures/panel_b_per_read_error.png` (and `.pdf`)
**Required Inputs:**
- `*_per_read_error_rates.tsv` (for all native and IVT samples)

### Panel C: Coverage Mismatch Error Rate Combined
**Script:** `scripts/panels/panel_c_coverage_mismatch.py`
**Output:** `figures/panel_c_coverage_mismatch.png` (and `.pdf`)
**Required Inputs:**
- `*_chr*_coverage_binned_error_rates.tsv` (for all native and IVT samples)

### Panel D: Substitution Error Rate
**Script:** `scripts/panels/panel_d_substitution_error.py`
**Output:** `figures/panel_d_substitution_error.png` (and `.pdf`)
**Required Inputs:**
- `no_chrY_sub_summary.tsv`

### Panel E: Per-Site Error Rate by Modification
**Script:** `scripts/panels/panel_e_per_site_error.py`
**Output:** `figures/panel_e_per_site_error.png` (and `.pdf`)
**Required Inputs:**
- `*_modkit_stratified_rates.tsv` (for all native samples)

### Panel F: Substitution Spectrum
**Script:** `scripts/panels/panel_f_substitution_spectrum.py`
**Output:** `figures/panel_f_substitution_spectrum.png` (and `.pdf`)
**Required Inputs:**
- `*_substitution_spectrum.tsv` (for native samples)

### Panel G: RDD Sites vs Threshold
**Script:** `scripts/panels/panel_g_rdd_sites.py`
**Output:** `figures/panel_g_rdd_sites.png` (and `.pdf`)
**Required Inputs:**
- `genomewide_rdd.json` (Note: Other JSON outputs like `genomewide_matrix.json` and `genomewide_permod.json` are generated by the upstream script `rdd_genomewide.py` but are not directly required by this panel script)

---

*Note: The scripts were adapted to use the local `inputs` and `figures` directory relative to the repository root. All scripts will automatically generate outputs into the `figures/` directory.*

## Upstream Data Generation

The input files needed by the panel scripts are generated upstream from aligned BAM files or intermediate pileups. The scripts responsible for this generation have been included in `scripts/generation/`.

### 1. Per-Read Error Rates
**Script:** `scripts/generation/calculate_read_error_rates_parallel.py`

Parses BAM alignments and calculates the total error rate per read (required for Panel B).

**Parameters:**
- `--mode`: Select `chrom` (process a single BAM/chromosome) or `merge` (aggregate chromosome results).
- *In `chrom` mode:*
  - `--bam`: Path to the indexed BAM file.
  - `--chrom`: Chromosome to process (e.g., `chr1`).
- *In `merge` mode:*
  - `--per-read-tsvs`: List of TSV files generated from `chrom` mode.
  - `--npz-files`: List of `.npz` files generated from `chrom` mode.

### 2. Per-Position Mismatch and Modkit Rates
**Script:** `scripts/generation/calculate_native_mismatch_rates_memeff_fast.py`

Takes `pysamstats` output and generates the `*_per_position_mismatch.tsv`, `*_mismatch_rate_summary.tsv`, and `*_modkit_stratified_rates.tsv` files required for Panels D and F. 

**Parameters:**
- `--mode`: Select `chrom` (process a single pysamstats file) or `merge` (aggregate).
- *In `chrom` mode:*
  - `--pysamstats`: Path to the per-chromosome pysamstats TSV output (can be gzipped).
  - `--masks-dir`: Directory containing `analysis_regions.bed`, `variants.bed`, and `junctions.bed` for site masking.
  - `--chrom`: Chromosome being processed (e.g., `chr1`).
- *In `merge` mode:*
  - `--per-pos-files`, `--rates-files`, `--summary-files`: Lists of intermediate files from `chrom` mode to be aggregated.
  - `--sample-name`: Base name for the final merged output files.

### 3. Calculate Coverage Error Rate
**Script:** `scripts/generation/calculate_coverage_error_rate.py`

Calculates per-position error rates binned by coverage depth directly from `pysamstats` TSVs (required for Panel C).

**Parameters:**
- `--pysamstats` (Required): Path to the `pysamstats --type variation` TSV output.
- `--giab-bed`: Path to GIAB high-confidence BED for position filtering (optional).
- `--chroms`: Comma-separated list of chromosomes to analyse (default: all autosomes + chrX).
- `--outdir`: Directory to save the output files (default: `results/coverage_error`).
- `--sample-label`: Label for the sample in the output TSV (default: `sample`).
- `--n-deciles`: Number of depth deciles to bin the output into (default: 20).

### 4. Genome-wide RDD Analysis (Alternative)
**Script:** `scripts/generation/rdd_genomewide.py`

Generates the JSON summary files (`genomewide_matrix.json`, `genomewide_permod.json`, `genomewide_rdd.json`) used to generate alternative RDD panels, including Panel G, from ONT per-position pileups.

**Parameters (Positional):**
- `MODE`: Must be one of `matrix`, `rdd`, or `permod`.
  - `matrix` or `rdd`: Provide the mode followed by a list of `*variation_tsv.gz` pileups.
  - `permod`: Provide `permod`, followed by a `BEDRMOD.bed` file, followed by the pileups.
- *Internal Constants (edit within script):*
  - `MINCOV`: Minimum coverage (default: 20).
  - `SNV_VAF`: Variant Allele Frequency threshold to exclude SNVs (default: 0.35).
