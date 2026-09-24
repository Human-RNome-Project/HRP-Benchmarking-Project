# RNome m6A Pipeline

This directory contains the R pipeline used to detect and quantify N6-methyladenosine (m6A)
sites from eTAM-seq metadata tables, and to compare methylation levels between mRNA samples
(`mRNA`) and their *in vitro* transcribed control (`IVT`), which is unmethylated by design.

The pipeline loads per-sample site tables, annotates each A site with its sequence motif
(DRACH / DRAC / RAC / RRAT), computes a per-site methylation ratio and a binomial p-value
against the sequencing error rate, runs a differential methylation analysis with `edgeR`
(quasi-likelihood framework applied to methylation counts), and exports the significant sites
as tables and bedRMod files.

## Directory Structure

- `pipeline.R` # Main script: per-site binomial test + differential analysis + bedRMod export
- `README.md` # Contains project documentation

Output directories are created by the script itself inside `PATH` (see *Parameters*):

- `plots/` # QC figures (number of A sites, MDS/PCA)
- `input_DE/` # Per-sample Bismark-like count files consumed by edgeR
- `tables/` # Result tables

## Scripts/Workflows

| File | Purpose |
|------|---------|
| `pipeline.R` | Loads eTAM metadata tables, annotates motifs, computes coverage / conversion rate / binomial p-value / FDR per site, exports per-sample bedRMod files, merges samples, runs the mRNA vs IVT differential test and exports significant sites. |

## Requirements & Setup

### Environment setup

To run the scripts you need an **R (≥ 4.3)** environment with Bioconductor. Required packages:

- Bioconductor: `edgeR`, `limma`, `GenomicRanges`, `GenomicFeatures`, `Biostrings`,
  `rtracklayer`, `Rsamtools`, `txdbmaker`, `org.Dm.eg.db`, `AnnotationDbi`
- CRAN: `tidyverse` (incl. `dplyr`, `stringr`, `purrr`, `ggplot2`), `data.table`, `writexl`

Installation:

```r
install.packages(c("tidyverse", "data.table", "writexl", "BiocManager"))
BiocManager::install(c("edgeR", "limma", "GenomicRanges", "GenomicFeatures",
                       "Biostrings", "rtracklayer", "Rsamtools", "txdbmaker",
                       "org.Dm.eg.db", "AnnotationDbi"))
```

No Conda environment or Docker image is used; the pipeline runs directly in R/RStudio.

### Resources

The `estimateDisp` and `glmQLFit` steps run on several million sites and take **tens of minutes
to a few hours** with a few GB of RAM. Run them on a machine with enough memory, or on a subset
of chromosomes for testing.

## Input data

### Data Download Instructions

Input tables are produced upstream by the eTAM-seq pre-processing pipeline (base-calling +
per-position A/G counting) in /data_processing_pipelines/SRS_workflows/eTAM-seq/PreProcess/. Point the `PATH`
parameter to the folder containing the two sub-folders `ETAM_RNA_mrna/metadata/` and
`ETAM_RNA_IVT/metadata/`.

| Input File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `ETAM_RNA_mrna/metadata/metadata_form_mRNA_<N>.tsv` | Per-site counts for the mRNA (methylated) replicates, pre-filtered to coverage ≥ 10 reads. | Tab-separated, with header | `pipeline.R` |
| `ETAM_RNA_IVT/metadata/metadata_form_IVT_<N>.tsv` | Per-site counts for the IVT (unmethylated control) replicates. | Tab-separated, with header | `pipeline.R` |

Expected columns in the metadata tables:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | 1-based genomic position of the A site |
| `nt` | Reference nucleotide (`A` → + strand, `T` → − strand) |
| `context` | Local sequence context; the motif is extracted from positions 2..n-1 |
| `nbA` | Number of reads still reading A (protected = methylated) |
| `nbG` | Number of reads converted to G (unmethylated) |
| `nbA.G` | Total coverage A+G at the site |
| `kind`, `gene` | Genomic feature and gene annotation, carried through to the result tables |

Files whose names contain hyphens are automatically renamed with underscores
(`rename_files_with_hyphens`) so they can be used as R object names.

### Parameters to adapt

All parameters are at the top of `pipeline.R`:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `PATH` | local path | Root folder containing the data; also the working directory and where outputs are written |
| `CTRL` / `nCTRL` | `mRNA` / 3 | Name and number of replicates of the biological condition |
| `IVT` / `nIVT` | `IVT` / 2 | Name and number of replicates of the unmethylated control |
| `EXT` | `.tsv` | Extension of the metadata files |
| `pattern_to_remove` | `metadata_form_` | Prefix stripped to build clean sample names |
| `MIN_COVERAGE` | 10 | Minimum read coverage per site |
| `MIN_SAMPLES` | `min(nIVT, nCTRL)` | Minimum number of samples reaching `MIN_COVERAGE` |
| `FDR_THRESHOLD` | 0.05 | Significance threshold |
| `DELTA_THRESHOLD` | 0.1 | Minimum methylation difference (mRNA − IVT) |
| `p_err` | 0.001 | Assumed sequencing/conversion error rate for the per-site binomial test |

Motif definitions (`DRACH`, `DRAC`, `RAC`, `RRAT`) and their regex patterns are also defined
there, as well as the plotting style (`theme_custom`).

### References/Accessory files

- Genome/annotation: the pipeline reuses the annotation columns already present in the
  metadata tables (`chr`, `pos`, `gene`, `kind`); `org.Dm.eg.db` is loaded for optional
  gene-symbol conversion (*GRCh38*).
- The bedRMod header written by `write_bedRMod` is currently hardcoded to
  `#organism=9606` / `#assembly=GRCh38` / Ensembl 93.

## Usage/Step-wise Execution

The script is meant to be run interactively (RStudio) or in one go:

```bash
Rscript pipeline.R
```
—
Step-by-step:

1. **Setup** — set `PATH`, `CTRL`/`IVT` and the replicate numbers; the script sets the working
   directory and creates `plots/`, `input_DE/` and `tables/`.
2. **Load data** — `rename_files_with_hyphens()` then `import_annot_files()` load every
   `.tsv` of both metadata folders into the global environment, one data frame per sample.
3. **Motif annotation** — `add_motif_columns()` adds `protected` (= `nbA / nbA.G`, the
   methylation ratio), `motif`, the four motif flags, `coord` (`chr-pos`) and `strand`.
4. **Per-site statistics** — coverage, conversion rate, one-sided binomial p-value against
   `p_err` and BH-adjusted FDR are computed for every site of every sample.
5. **QC** — barplot of the number of A sites per sample (`plots/nb_A_sites_per_sample.pdf`)
   and MDS/PCA plot of the `protected` matrix (`plots/PCA_plot.png`).
6. **Per-sample export** — bedRMod files for all sites and for sites with p ≤ 0.05.
7. **Merge and export for edgeR** — `merge_protected_columns()` keeps sites present in all
   samples; `export_for_edgeR()` writes one Bismark-like file per sample in `input_DE/`.
8. **Differential analysis** — `pdata`/`pdata_sub` define the design; edgeR reads the
   `input_DE/` files, filters on coverage, normalises library sizes, fits a quasi-likelihood
   model and tests the `IVT - mRNA` contrast (this block lives in `pipeline_all.R`, see
   *Known issues*).
9. **Annotation and filtering** — results are joined back to the motif/gene annotation,
   `mean_ctrl`, `mean_IVT` and `delta_ivt = mean_ctrl - mean_IVT` are computed, and sites
   are kept when `delta_ivt ≥ DELTA_THRESHOLD` and `PValue ≤ FDR_THRESHOLD`.
10. **Export** — final bedRMod files per condition.

## Outputs

| Output File | Format | Script | Description |
|-------------|--------|--------|-------------|
| `plots/nb_A_sites_per_sample.pdf` | PDF | `pipeline.R` | Number of A sites with coverage ≥ 10 reads per sample |
| `plots/PCA_plot.png` | PNG | `pipeline.R` | MDS/PCA of samples based on the methylation ratios |
| `pdata.csv` | CSV | `pipeline.R` | Sample sheet: sample name, number of sites, group |
| `input_DE/<sample>.txt` | TSV, no header | `pipeline.R` | edgeR input: chr, start, end, methylation ratio, nbA, nbG |
| `<sample>_all_bedRMod.bed` | bedRMod v2 | `pipeline.R` | All sites of a sample, score = binomial p-value, frequency = methylation ratio |
| `<sample>_pval0.05_bedRMod.bed` | bedRMod v2 | `pipeline.R` | Same, restricted to sites with p ≤ 0.05 |
| `mRNA_all_bedRMod.bed` / `mRNA_pV0.05_bedRMod.bed` | bedRMod v2 | `pipeline.R` | Differential results for the mRNA condition (all sites / significant sites) |
| `IVT_all_bedRMod.bed` / `IVT_pV0.05_bedRMod.bed` | bedRMod v2 | `pipeline.R` | Differential results for the IVT condition (all sites / significant sites) |
| `fit_<CTRL>_<IVT>.RData` | RData | `pipeline_all.R` | Saved `glmQLFit` object, to avoid re-running the long fitting step |
