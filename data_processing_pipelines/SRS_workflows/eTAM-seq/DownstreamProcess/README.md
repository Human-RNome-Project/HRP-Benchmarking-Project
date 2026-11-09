# Rnome m6A Pipeline

This directory contains the R pipeline used to detect and quantify N6-methyladenosine (m6A)
sites from eTAM-seq metadata tables (human, GRCh38), and to compare methylation levels between
mRNA samples (`mRNA`) and their *in vitro* transcribed control (`IVT`), which is unmethylated
by design.

The pipeline loads per-sample site tables, annotates each A site with its sequence motif
(DRACH / DRAC / RAC / RRAT), computes a per-site methylation ratio and a binomial p-value
against the sequencing error rate, runs a differential methylation analysis with `edgeR`
(quasi-likelihood framework applied to methylation counts), and exports the results as TSV
tables and bedRMod v2 files.

## Directory Structure

- `pipeline.R` # Self-contained script: per-site statistics, edgeR differential analysis, table and bedRMod export
- `README.md` # Contains project documentation

Output directories are created by the script inside `out_dir` (see *Usage*):

- `plots/` # QC figures (number of A sites, MDS/PCA)
- `input_DE/` # Per-sample Bismark-like count files consumed by edgeR
- `tables/` # Result tables (TSV)
- `bedRMod/` # bedRMod v2 files

## Scripts/Workflows

| File | Purpose |
|------|---------|
| `pipeline.R` | Loads the eTAM metadata tables, annotates motifs, computes coverage / conversion rate / binomial p-value / FDR per site, exports per-sample bedRMod files, runs the edgeR test `IVT - mRNA` (`readBismark2DGE` → `estimateDisp` → `glmQLFit` → `glmQLFTest`), and exports the result tables and per-condition bedRMod files. |

## Requirements & Setup

### Environment setup

To run the script you need an **R (≥ 4.3)** environment with Bioconductor. Required packages:

- Bioconductor: `edgeR`, `limma`
- CRAN: `dplyr`, `stringr`, `tibble`, `ggplot2`, `data.table`

Installation:

```r
install.packages(c("dplyr", "stringr", "tibble", "ggplot2", "data.table", "BiocManager"))
BiocManager::install(c("edgeR", "limma"))
```

No Conda environment or Docker image is used.

### Resources

The `estimateDisp` and `glmQLFit` steps run on several million sites and take **tens of minutes
to a few hours** with a few GB of RAM. Run them on a machine with enough memory, or on a subset
of chromosomes for testing. The fitted model is saved as `fit_<CTRL>_<IVT>.RData`.

## Input data

### Data Download Instructions

Input tables are produced upstream by the eTAM-seq pre-processing pipeline (see `/HRP-Benchmarking-Project/data_processing_pipelines/SRS_workflows/eTAM-seq/PreProcess/`:
trimming, mapping, per-position A/G counting). They are not included in this repository.
They are the `<sample>.annot.txt` files written by `/HRP-Benchmarking-Project/data_processing_pipelines/SRS_workflows/eTAM-seq/PreProcess/run`, unchanged, copied under the
names below by `/HRP-Benchmarking-Project/data_processing_pipelines/SRS_workflows/eTAM-seq/PreProcess/3_make_metadata.sh <annot_dir> <data_dir>`
(`Rnome-mRNA-<N>` → `metadata_form_mRNA_<N>.tsv`, `Rnome-IVT-<N>` → `metadata_form_IVT_<N>.tsv`).
`data_dir` must contain the two sub-folders below; every `*.tsv` file in them is loaded.

| Input File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `<data_dir>/ETAM_RNA_mrna/metadata/metadata_form_mRNA_<N>.tsv` | Per-site counts for the mRNA (methylated) replicates, pre-filtered to coverage ≥ 10 reads. | Tab-separated, with header | `pipeline.R` |
| `<data_dir>/ETAM_RNA_IVT/metadata/metadata_form_IVT_<N>.tsv` | Per-site counts for the IVT (unmethylated control) replicates. | Tab-separated, with header | `pipeline.R` |

Sample names are taken from the file names: the `.tsv` extension and the `metadata_form_` prefix
are removed and hyphens are replaced by underscores (the files themselves are not renamed).
The group is the part before the first `_` and must be `mRNA` or `IVT`; at least two replicates
per group are required. The number of replicates is read from the data, not set by hand.

Expected columns in the metadata tables:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | 1-based genomic position of the A site |
| `nt` | Reference nucleotide (`A` → + strand, `T` → − strand) |
| `context` | Local sequence context; the motif is extracted from positions 2..n-1 |
| `nbA` | Number of reads still reading A (protected = methylated) |
| `nbG` | Number of reads converted to G (unmethylated) |
| `nbA.G` | Total coverage A+G at the site (header `nbA+G` in the file, renamed by R at reading) |
| `kind`, `gene` | Genomic feature and gene annotation, carried through to the result tables |

The files also hold `nbStart`, `nbStop`, `aa`, `frame`, `donor`, `acceptor` (not used here).
`kind`, `context` and `gene` are `NA` for sites outside annotated transcripts; their motif flags are then `NA`.

### Parameters to adapt

Parameters are at the top of `pipeline.R`:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `CTRL` | `mRNA` | Name of the biological condition |
| `IVT` | `IVT` | Name of the unmethylated control |
| `EXT` | `.tsv` | Extension of the metadata files |
| `pattern_to_remove` | `metadata_form_` | Prefix stripped to build clean sample names |
| `MIN_COVERAGE` | 10 | Minimum read coverage per site |
| `MIN_SAMPLES` | `min(nIVT, nCTRL)` | Minimum number of samples reaching `MIN_COVERAGE` (computed from the data) |
| `FDR_THRESHOLD` | 0.05 | Significance threshold (applied to the edgeR `PValue`) |
| `DELTA_THRESHOLD` | 0.1 | Minimum methylation difference (mRNA − IVT) |
| `P_NONCONV` | `NA` (estimated) | Background non-conversion rate for the per-site binomial test; `NA` = pooled `sum(nbA) / sum(nbA.G)` over the IVT samples |
| `BEDRMOD_*` | see below | bedRMod header fields |

Motif definitions (`DRACH`, `DRAC`, `RAC`, `RRAT`) and their regex patterns are also defined
there, as well as the plotting style (`theme_custom`).

### References/Accessory files

- Genome/annotation: the pipeline reuses the annotation columns already present in the
  metadata tables (`chr`, `pos`, `gene`, `kind`), which come from the RefSeq GenBank records of
  **GRCh38.p12** (`GCF_000001405.38`) used in `/HRP-Benchmarking-Project/data_processing_pipelines/SRS_workflows/eTAM-seq/PreProcess/`. No annotation package is loaded. These can be found in the DOE Data Explorer https://doi.org/10.25585/DOE-HRP/3377574 Raw data References folder.
- bedRMod header: `#organism=9606`, `#assembly=GRCh38.p12`, `#annotation_source=RefSeq`,
  `#annotation_version=GCF_000001405.38_GRCh38.p12`, `#sequencing_platform=Illumina NextSeq 2000`
  (parameters `BEDRMOD_*`).

## Usage/Step-wise Execution

```bash
Rscript pipeline.R <data_dir> [out_dir]
```

- `data_dir` (required): folder containing `ETAM_RNA_mrna/metadata/` and `ETAM_RNA_IVT/metadata/`.
- `out_dir` (optional): where the outputs are written; defaults to `data_dir`.

In RStudio, set `DATA_DIR` (and optionally `OUT_DIR`) in the console, then `source("pipeline.R")`.

Step-by-step:

1. **Setup**: read the arguments, create `plots/`, `input_DE/`, `tables/` and `bedRMod/` in `out_dir`.
2. **Load data**: every `.tsv` of both metadata folders is read into a list, one data frame per sample.
3. **Motif annotation**: `add_motif_columns()` adds `protected` (= `nbA / nbA.G`, the
   methylation ratio), `motif`, the four motif flags, `coord` (`chr-pos`) and `strand`.
4. **Per-site statistics**: coverage, conversion rate, one-sided binomial p-value of `nbA`
   (A protected from deamination, i.e. m6A candidate) against the background non-conversion
   rate `P_NONCONV`, and BH-adjusted FDR, for every site of every sample. `P_NONCONV` is
   estimated from the IVT samples (unmethylated by design) unless set by hand.
5. **Sample sheet and QC**: `pdata.csv`, per-sample summary in the console, barplot of the
   number of A sites (`plots/nb_A_sites_per_sample.pdf`), MDS plot of the `protected` matrix on
   the sites shared by all samples (`plots/PCA_plot.png`).
6. **Per-sample export**: bedRMod files for all sites and for sites with p ≤ 0.05.
7. **Export for edgeR**: one Bismark-like file per sample in `input_DE/`
   (chr, start, end, ratio, methylated = `nbA`, unmethylated = `nbG`).
8. **Differential analysis**: edgeR reads the `input_DE/` files, keeps sites with coverage
   ≥ `MIN_COVERAGE` in ≥ `MIN_SAMPLES` samples, sets one library size per sample, fits a
   quasi-likelihood model (`modelMatrixMeth`) and tests the `IVT - mRNA` contrast.
9. **Annotation and filtering**: results are joined to the motif/gene annotation; `mean_ctrl`,
   `mean_IVT` (mean methylation ratio), `cov_ctrl`, `cov_IVT` (mean coverage) and
   `delta_ivt = mean_ctrl - mean_IVT` are computed. Sites are significant when
   `delta_ivt ≥ DELTA_THRESHOLD` and `PValue ≤ FDR_THRESHOLD`.
10. **Export**: result tables and per-condition bedRMod files.

## Outputs

All paths are relative to `out_dir`.

| Output File | Format | Description |
|-------------|--------|-------------|
| `plots/nb_A_sites_per_sample.pdf` | PDF | Number of A sites with coverage ≥ 10 reads per sample |
| `plots/PCA_plot.png` | PNG | MDS of samples based on the methylation ratios |
| `pdata.csv` | CSV | Sample sheet: sample name, number of sites, group, edgeR input file |
| `input_DE/<sample>.txt` | TSV, no header | edgeR input: chr, start, end, methylation ratio, nbA, nbG |
| `fit_<CTRL>_<IVT>.RData` | RData | Saved `glmQLFit` object |
| `tables/<CTRL>_vs_<IVT>_all.tsv` | TSV | All tested sites with edgeR statistics, annotation, group means and `delta_ivt` |
| `tables/<CTRL>_vs_<IVT>_delta0.1_p0.05.tsv` | TSV | Significant sites (`delta_ivt ≥ 0.1`, `PValue ≤ 0.05`) |
| `bedRMod/<sample>_all_bedRMod.bed` | bedRMod v2 | All sites of a sample; score from the binomial p-value (`nbA` vs `P_NONCONV`) |
| `bedRMod/<sample>_pval0.05_bedRMod.bed` | bedRMod v2 | Same, restricted to sites with binomial p ≤ 0.05 |
| `bedRMod/mRNA_all_bedRMod.bed` / `mRNA_pV0.05_bedRMod.bed` | bedRMod v2 | mRNA condition: all tested sites / sites with edgeR `PValue ≤ 0.05`; score from the edgeR p-value |
| `bedRMod/IVT_all_bedRMod.bed` / `IVT_pV0.05_bedRMod.bed` | bedRMod v2 | Same for the IVT condition |

### bedRMod columns

The header lines are followed by a tab-separated column header line:
`#chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb coverage frequency`.

| Column | Content |
|--------|---------|
| `name` | `m6A` |
| `score` | Integer 0 to 1000: `round(-10 × log10(p))`, capped at 1000 (p = 0.05 → 13, p = 1e-10 → 100) |
| `strand` | `+`, `-`, or `.` when unknown |
| `coverage` | Read coverage (`nbA.G`); for the per-condition files, mean over the replicates of the condition |
| `frequency` | Methylation level as an integer percentage 0 to 100 (`round(100 × protected)`); for the per-condition files, mean over the replicates |

Sites with no frequency or coverage for a condition (site not covered in every replicate of that
condition) are left out of that condition's bedRMod file.
