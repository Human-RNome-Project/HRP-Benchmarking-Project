# Figure 11. Short-read sequencing (SRS) of polyA-enriched RNAs

This directory contains the plotting scripts used to generate Figure11 in the manuscript.

## Scripts

| File | Purpose |
|------------|-------------|
| `0.0_anno_to_gene.R`  | Annotate modifications sites to 5' UTR, CDS and 3' UTR of canonical ensembl transcript of genes |
| `1.0_integrated_metagene.R` | Plot number of modifications sites for polyA-enriched RNAs; Plot region distrition and metagene plot for relative location in transcriptome;Plot mean modification levels along the transcriptome |
| `2.0_integrated_motif.R` | Seq log plot for polyA modifications |
| `3.0_integrated_top_Gene.R` | Barplot for top20 genes with highest average number of modifications sites |


## Dependencies

These scripts require the following R packages:

- R>=4.3.3
- tidyverse>=2.0.0
- GenomicFeatures>=1.58.0
- clusterProfiler>=4.14.0
- org.Hs.eg.db>=3.20.0
- UpSetR>=1.4.0
- patchwork>=1.3.2
- ggseqlogo>=0.2.2
- data.table>=1.15.4
- rtracklayer>=1.66.0
- qs>=0.26.3
- ggsci>=3.2.0

bedtools >=2.31.1 is also needed to get motif from genome FASTAs. 


# Usage/Step-wise Execution

## Step 1: prepare GTF file, hg38 FASTA and merged bedRmod for Illumina

```sh
indir="test_input"
outdir="test_output"

## 1. download GTF file
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
mv gencode.v49.primary_assembly.annotation.gtf.gz ${indir}

## 2. download genome reference used to get motif
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa ${indir}

## 3. get `Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed` from upstream pipelines, could also be downloaded from our UCSC genome browser
```

## Step 2: Annotate modifications sites to 5' UTR, CDS and 3' UTR of canonical ensembl transcript of genes

### Input files

| Input File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `${indir}/Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed`  | Merged bedRmod for Illumina platform derived from upstream pipelines, could also be downloaded from our UCSC genome browser | Tab-separated bedRMod (`.bed`) | `0.0_anno_to_gene.R` |
| `${indir}/gencode.v49.primary_assembly.annotation.gtf.gz` | GTF file from Gencode | GTF format |  `0.0_anno_to_gene.R` |


### Usage

```sh
indir="test_input"
outdir="test_output"
Rscript 0.0_anno_to_gene.R -i ${indir}/Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed -o ${outdir} \
  --gtf ${indir}/gencode.v49.primary_assembly.annotation.gtf.gz

## Note: 
# 1. This step will consume several minutes as we need to read GTF files
# 2. If chrY sites are incluced in the input file, it will also be removed!
```

### Output files

| Output File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `${outdir}/Anno_database_canonical_transcripts.qs`  | Annotation files created from GTF, used to annotate modification sites to genes | Binary file created by R qs package | `0.0_anno_to_gene.R` |
| `${outdir}/Illumina_polyA_mod_annotated.tsv`  | Modification sites were annotated to ensemble canonical transcripts of genes. Sites in exon regions are also annotated by its relative location in the transcript: utr5 (0-1), cds (1-2), utr3 (2-3) | Tab-separated files | `0.0_anno_to_gene.R` |


## Step 3: Plot number of modification sites and metagene plot of polyA modifications

### Input files

| Input File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `${indir}/Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed`  | Merged bedRmod for Illumina platform derived from upstream pipelines, could also be downloaded from our UCSC genome browser | Tab-separated bedRMod (`.bed`) | `1.0_integrated_metagene.R` |
| `${outdir}/Illumina_polyA_mod_annotated.tsv` | Modification sites with annotated from last step | Tab-separated files |  `1.0_integrated_metagene.R` |

### Usage

```sh
indir="test_input"
outdir="test_output"
Rscript 1.0_integrated_metagene.R -i ${indir}/Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed \
  -a ${outdir}/Illumina_polyA_mod_annotated.tsv \
  -o ${outdir}
```

### Output files

| Output File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `${outdir}/Integrated_barplot_num_of_sites.pdf`  | Barplot showing number of sites for each modification type and each RNA type | PDF | `1.0_integrated_metagene.R` |
| `${outdir}/Integrated_region_proportion.pdf`  | Proportion of poluA modification sites in exonic, intronic, and intergenic regions for each modification. | PDF | `1.0_integrated_metagene.R` |
| `${outdir}/Integrated_metagene.pdf`  | Metagene distribution of polyA modifications across 5′UTR, CDS, and 3′UTR | PDF | `1.0_integrated_metagene.R` |
| `${outdir}/Mean_ratio_along_transripts.pdf`  | Mean modification level across meta-transcript regions | PDF | `1.0_integrated_metagene.R` |


## Step 4: Plot motif logos

### Input files

| Input File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `${indir}/Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed`  | Merged bedRmod for Illumina platform derived from upstream pipelines, could also be downloaded from our UCSC genome browser | Tab-separated bedRMod (`.bed`) | `2.0_integrated_motif.R` |
| `${indir}/GRCh38.primary_assembly.genome.fa` | Genome reference file used to get motif | FASTA | `2.0_integrated_motif.R` |


### Usage

```sh
indir="test_input"
outdir="test_output"
Rscript 2.0_integrated_motif.R -i ${indir}/Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed \
  -f ${indir}/GRCh38.primary_assembly.genome.fa \
  -o ${outdir}

# Note: 
# 1. bedtools gesfasta is called within the R script
# 2. Several minutes needed to get motif of Inosine sites
```


### Output files

| Output File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `${outdir}/polyA_motif_bits.pdf`  | Seq log plot for polyA modifications | PDF | `2.0_integrated_motif.R` |


## Step 5: Plot top genes

### Input files

| Input File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `${indir}/Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed`  | Merged bedRmod for Illumina platform derived from upstream pipelines, could also be downloaded from our UCSC genome browser | Tab-separated bedRMod (`.bed`) | `3.0_integrated_top_Gene.R` |
| `${outdir}/Illumina_polyA_mod_annotated.tsv` | Modification sites with annotated from step 2 | Tab-separated files |  `3.0_integrated_top_Gene.R` |


### Usage

```sh
indir="test_input"
outdir="test_output"
Rscript 3.0_integrated_top_Gene.R -i ${indir}/Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed \
  -a ${outdir}/Illumina_polyA_mod_annotated.tsv \
  -o ${outdir}

# Note: 
# 1. Several minutes needed to get motif of Inosine sites
```


### Output files

| Output File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `${outdir}/Integrated_top20_Coding_gene_by_mod_levels.pdf`  | Barplot for top20 genes with highest average number of modifications sites | PDF | `3.0_integrated_top_Gene.R` |