# Figure 11. Short-read sequencing (SRS) of polyA-enriched RNAs

This directory contains the plotting scripts used to generate Figure11 in the manuscript.

## Scripts

- 0.0_anno_to_gene.R
  - Annotate modifications sites to canonical ensembl transcript of genes
  - Calculate relative location in 5' UTR, CDS and 3' UTR for metagene plot

- 1.0_compare_replicates.R
  - UpSet plot to compare modification sites in replicates

- 2.0_compare_methods.R
  - UpSet plot to compare modifications sites from different methods

- 3.0_integrated_metagene.R
  - Plot number of modifications sites for polyA-enriched RNAs
  - Plot distribution of modification levels for each method
  - Plot region distrition and metagene plot for relative location in transcriptome
  - Plot mean modification levels along the transcriptome

- 3.1_integrated_top_Gene.R
  - Get genes with highest average number of modification sites per transcript

- 3.2_integrated_motif.R
  - Motif plot


## Dependencies

These scripts require the following R packages:

- tidyverse
- GenomicFeatures
- clusterProfiler
- org.Hs.eg.db
- UpSetR
- patchwork
- ggseqlogo


## Inputs

- Gencode GTF annotation: gencode.v49.primary_assembly.annotation.gtf
  - Canonical ensembl transcript for each gene are labelled by tag "Ensembl_canonical" in the GTP file
- Combined bedRmod file for SRS modification sites


## Typical Workflow

1. Prepare the combined bedRmod file for SRS modification sites and download Gencode GTF file
2. Install the depencencies
3. Modify the file path and run each script

