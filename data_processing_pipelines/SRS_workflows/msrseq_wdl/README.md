# MSR-seq tRNA Modification WDL Pipeline

This directory contains the WDL workflow and input templates for processing MSR-seq (Multiplex Small RNA Sequencing) tRNA modification samples.

## Workflows

| WDL file | Workflow namespace | Purpose |
|---|---|---|
| `msrseq_pipeline.wdl` | `msrseq_pipeline` | Full pipeline for tRNA modification detection: Reference preparation & Bowtie2 indexing → Read 2 sense local alignment → SAM extension binning → BAM conversion & IGV base pileup → Multi-sample aggregation (`tRNAmod_cleaner.R`) → Stoichiometry binomial scoring & bedRMod conversion (`bedrmod_converter.R`) |

## Requirements

### Workflow engine
You need a WDL-compatible execution backend. The production runs used **JAWS** on Perlmutter, but the workflows can also be run with **Cromwell** or **miniwdl** locally.

Examples:
- JAWS: `jaws submit --no-cache <wdl> <inputs.json> perlmutter --tag <tag>`
- Cromwell: `java -jar cromwell.jar run <wdl> -i <inputs.json>`
- miniwdl: `miniwdl run <wdl> --input <inputs.json>`

### Docker
All tasks run inside Docker containers. The executor must be able to pull or load the container image listed below. If you are running offline, pull or build it beforehand.

### Input data
- For `msrseq_pipeline.wdl`:
  - Mature tRNA reference FASTA (`mature_trna_fasta`, e.g. `Step1_260_seq_hg38-mature-tRNAs.fa`)
  - Chromosomal tRNA reference FASTA (`chromosomal_trna_fasta`, e.g. `hg38_chromosomal_tRNA_genes_high_confidence_intro_remove+CCA_upper_case2.fa`)
  - Isodecoder annotation table (`isodecoder_table`, e.g. `Table S3_cyto all isodecoders3.xlsx`)
  - Array of `MSRSample` structs specifying `sample_id`, `treatment`, `replicate`, `barcode`, and demultiplexed Read 2 FASTQ (`fastq_read2`)

## Docker images

The workflow uses a single unified container image containing all required bioinformatics tools (Bowtie2, Samtools, IGVTools/Java 17), Python with pysam/numpy/pandas, R with tidyverse/openxlsx/Biostrings, and all packaged pipeline scripts.

```bash
# Unified MSR-seq pipeline image
# Used by all tasks: PrepareReference, ProcessSample, AggregateTsvs, BedrmodConversion
docker pull <your_dockerhub_user>/msrseq-pipeline:v1

# Or build locally from the provided Dockerfile:
docker build -t msrseq-pipeline:v1 docker/
```

> **Note:** On NERSC Perlmutter, Shifter automatically pulls and caches the image from Docker Hub or Quay.io when submitted via JAWS or Cromwell.

## Input JSON preparation

### NERSC Perlmutter (CFS paths)

Use the provided template `inputs_msrseq_perlmutter.json`:

```json
{
  "msrseq_pipeline.mature_trna_fasta": "/global/cfs/cdirs/m5243/raw_files/MSRseq/Step1_260_seq_hg38-mature-tRNAs.fa",
  "msrseq_pipeline.chromosomal_trna_fasta": "/global/cfs/cdirs/m5243/raw_files/MSRseq/hg38_chromosomal_tRNA_genes_high_confidence_intro_remove+CCA_upper_case2.fa",
  "msrseq_pipeline.isodecoder_table": "/global/cfs/cdirs/m5243/raw_files/MSRseq/Table_S3_cyto_all_isodecoders3.xlsx",
  "msrseq_pipeline.docker_image": "<your_dockerhub_user>/msrseq-pipeline:v1",
  "msrseq_pipeline.ref_cpus": 8,
  "msrseq_pipeline.sample_cpus": 16,
  "msrseq_pipeline.agg_cpus": 4,
  "msrseq_pipeline.convert_cpus": 4,
  "msrseq_pipeline.samples": [
    {
      "sample_id": "TP-AB-19s-R-Con_S12_L008_L1_bc8_GGTA",
      "treatment": "HRPC_ctrl",
      "replicate": 1,
      "barcode": "bc8",
      "fastq_read2": "/global/cfs/cdirs/m5243/raw_files/MSRseq/NovaSeqX_GM12878-1_MSRseq_260211_R1/L1_bc8_GGTA_2.txt.gz"
    }
  ]
}
```

Edit the CFS file paths and Docker image tag before submitting.

### Local testing

Use the provided local test input JSON `inputs_msrseq_local.json`, which points to local datasets in the repository and uses the local image tag `msrseq-pipeline:v1`.

## Usage

### Run on NERSC Perlmutter with JAWS

```bash
jaws submit \
  --no-cache \
  msrseq_pipeline.wdl \
  inputs_msrseq_perlmutter.json \
  perlmutter \
  --tag "msrseq_tRNA_run"
```

### Run on NERSC with Cromwell (Shifter + Slurm)

```bash
# Pull image into Shifter on Perlmutter login node
shifterimg pull <your_dockerhub_user>/msrseq-pipeline:v1

# Submit workflow
java -Dconfig.file=cromwell.nersc.conf -jar /usr/local/bin/cromwell.jar run \
  msrseq_pipeline.wdl \
  -i inputs_msrseq_perlmutter.json
```

### Run locally with miniwdl

```bash
miniwdl run \
  msrseq_pipeline.wdl \
  -i inputs_msrseq_local.json \
  --dir run_output
```

## Outputs

### `msrseq_pipeline.wdl`
- `cleaned_ref_fasta` — header-cleaned mature tRNA reference FASTA
- `bt2_index_files` — array of generated Bowtie2 index files
- `all_tsvs` — array of all binned base-level coverage TSV files across all samples
- `data_cleaned_zip` — aggregated intermediate tRNA modification dataset (`data_cleaned_5_HRPC.csv.zip`)
- `bedrmod_files` — array of all 12 final `bedRModv2` modification files
- `bs_rep1`, `bs_rep2`, `bs_rep3` — Bisulfite treatment deletion modification files
- `cbh_rep1`, `cbh_rep2`, `cbh_rep3` — Cyanoborohydride treatment mutation modification files
- `ctrl_mut_rep1`, `ctrl_mut_rep2`, `ctrl_mut_rep3` — Untreated control mutation modification files
- `ctrl_del_rep1`, `ctrl_del_rep2`, `ctrl_del_rep3` — Untreated control deletion modification files
