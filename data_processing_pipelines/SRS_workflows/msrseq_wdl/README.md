# MSR-seq tRNA Modification WDL Pipeline

This directory contains the WDL workflow, Docker container configuration, and input templates for processing MSR-seq (Multiplex Small RNA Sequencing) tRNA modification samples.

---

## Workflows

| WDL file | Workflow namespace | Purpose |
|---|---|---|
| `msrseq_pipeline.wdl` | `msrseq_pipeline` | Full pipeline for tRNA modification detection: Reference preparation & Bowtie2 indexing → Read 2 sense local alignment → SAM extension binning → BAM conversion & IGV base pileup → Multi-sample aggregation (`tRNAmod_cleaner.R`) → Stoichiometry binomial scoring & bedRMod conversion (`bedrmod_converter.R`) |

---

## Requirements & Software Setup

### Running with miniwdl (Recommended for Local Workstations & Cloud VMs)
For users without a NERSC account, the workflow can be executed locally or on cloud instances (e.g., Google Cloud Platform / GCP) using **miniwdl**:

1. **[Docker](https://docs.docker.com/engine/install/)**: Required to execute containerized pipeline tasks.
2. **[miniwdl](https://miniwdl.readthedocs.io/en/latest/getting_started.html#install-miniwdl)**: Python-based WDL runner.
   ```bash
   pip install miniwdl
   ```

### Running on HPC Clusters
* **JAWS / Perlmutter**: `jaws submit --no-cache msrseq_pipeline.wdl inputs_msrseq_perlmutter.json perlmutter --tag <tag>`
* **Cromwell**: `java -jar cromwell.jar run msrseq_pipeline.wdl -i inputs_msrseq_local.json`

---

## Docker Image Setup

Build the container image locally from the provided `docker/` directory:

```bash
docker build -t msrseq-pipeline:v1 docker/
```

This packages all required tools (Bowtie2 v2.4.4, Samtools v1.13, IGVTools v2.8.0, Java 17, Python with pysam/numpy/pandas, R v4.2.2 with tidyverse/openxlsx/Biostrings) and custom pipeline scripts.

> **Note:** Building the image locally is the simplest method. `miniwdl` will automatically use this local image. Optionally, you can tag and push this image to your container registry (Docker Hub / Quay.io) if deploying to remote clusters.

---

## Directory Layout & Input Data

The pipeline configuration uses local relative paths so that no absolute paths need to be configured:

```
msrseq_wdl/
├── msrseq_pipeline.wdl             # Master WDL workflow
├── docker/                         # Dockerfile and pipeline scripts
├── ref/                            # Bundled tRNA reference files
│   ├── Step1_260_seq_hg38-mature-tRNAs.fa
│   ├── hg38_chromosomal_tRNA_genes_high_confidence_intro_remove+CCA_upper_case2.fa
│   └── Table S3_cyto all isodecoders3.xlsx
├── raw_msrseq/                     # User-created directory for raw FASTQ files
│   ├── HRP_B_012_1_R2.fastq.gz
│   ├── ...
│   └── HRP_B_012_9_R2.fastq.gz
├── inputs_msrseq_local.json        # Template with relative paths
├── inputs_msrseq_gcp.json          # Template for GCP / cloud execution
└── inputs_msrseq_perlmutter.json   # Template with NERSC CFS paths
```

### 1. Reference Files (`ref/`)
The three required reference files are already bundled in the `ref/` subdirectory:
- `ref/Step1_260_seq_hg38-mature-tRNAs.fa`
- `ref/hg38_chromosomal_tRNA_genes_high_confidence_intro_remove+CCA_upper_case2.fa`
- `ref/Table S3_cyto all isodecoders3.xlsx`

### 2. FASTQ Files (`raw_msrseq/`)
Create a `raw_msrseq/` directory inside `msrseq_wdl/` and save the 9 demultiplexed Read 2 FASTQ files from the project data portal:
```bash
mkdir -p raw_msrseq
```
Save files named `HRP_B_012_1_R2.fastq.gz` through `HRP_B_012_9_R2.fastq.gz`.

---

## Input JSON Preparation

### Local / Cloud Execution (`inputs_msrseq_local.json` or `inputs_msrseq_gcp.json`)
Because reference files are in `ref/` and FASTQs are in `raw_msrseq/`, the provided JSON works immediately out-of-the-box without modifying any paths:

```json
{
  "msrseq_pipeline.mature_trna_fasta": "ref/Step1_260_seq_hg38-mature-tRNAs.fa",
  "msrseq_pipeline.chromosomal_trna_fasta": "ref/hg38_chromosomal_tRNA_genes_high_confidence_intro_remove+CCA_upper_case2.fa",
  "msrseq_pipeline.isodecoder_table": "ref/Table S3_cyto all isodecoders3.xlsx",
  "msrseq_pipeline.docker_image": "msrseq-pipeline:v1",
  "msrseq_pipeline.ref_cpus": 8,
  "msrseq_pipeline.sample_cpus": 8,
  "msrseq_pipeline.agg_cpus": 4,
  "msrseq_pipeline.convert_cpus": 4,
  "msrseq_pipeline.samples": [
    {
      "sample_id": "HRP_B_012_tRNA_005",
      "treatment": "HRPC_ctrl",
      "replicate": 1,
      "barcode": "bc8",
      "fastq_read2": "raw_msrseq/HRP_B_012_5_R2.fastq.gz"
    },
    {
      "sample_id": "HRP_B_012_tRNA_007",
      "treatment": "HRPC_ctrl",
      "replicate": 2,
      "barcode": "bc9",
      "fastq_read2": "raw_msrseq/HRP_B_012_7_R2.fastq.gz"
    },
    {
      "sample_id": "HRP_B_012_tRNA_009",
      "treatment": "HRPC_ctrl",
      "replicate": 3,
      "barcode": "bc10",
      "fastq_read2": "raw_msrseq/HRP_B_012_9_R2.fastq.gz"
    },
    {
      "sample_id": "HRP_B_012_tRNA_004",
      "treatment": "HRPC_BS",
      "replicate": 1,
      "barcode": "bc8",
      "fastq_read2": "raw_msrseq/HRP_B_012_4_R2.fastq.gz"
    },
    {
      "sample_id": "HRP_B_012_tRNA_006",
      "treatment": "HRPC_BS",
      "replicate": 2,
      "barcode": "bc9",
      "fastq_read2": "raw_msrseq/HRP_B_012_6_R2.fastq.gz"
    },
    {
      "sample_id": "HRP_B_012_tRNA_008",
      "treatment": "HRPC_BS",
      "replicate": 3,
      "barcode": "bc10",
      "fastq_read2": "raw_msrseq/HRP_B_012_8_R2.fastq.gz"
    },
    {
      "sample_id": "HRP_B_012_tRNA_001",
      "treatment": "HRPC_CBH",
      "replicate": 1,
      "barcode": "bc8",
      "fastq_read2": "raw_msrseq/HRP_B_012_1_R2.fastq.gz"
    },
    {
      "sample_id": "HRP_B_012_tRNA_003",
      "treatment": "HRPC_CBH",
      "replicate": 2,
      "barcode": "bc9",
      "fastq_read2": "raw_msrseq/HRP_B_012_3_R2.fastq.gz"
    },
    {
      "sample_id": "HRP_B_012_tRNA_002",
      "treatment": "HRPC_CBH",
      "replicate": 3,
      "barcode": "bc10",
      "fastq_read2": "raw_msrseq/HRP_B_012_2_R2.fastq.gz"
    }
  ]
}
```

### NERSC Perlmutter (CFS paths)
Use `inputs_msrseq_perlmutter.json` if running via JAWS or Cromwell on NERSC Perlmutter.

---

## Usage

### Run Locally with miniwdl
```bash
miniwdl run \
  msrseq_pipeline.wdl \
  -i inputs_msrseq_local.json \
  --dir run_output \
  --verbose
```
*(Or use `-i inputs_msrseq_gcp.json`)*

### Run on NERSC Perlmutter with JAWS
```bash
jaws submit \
  --no-cache \
  msrseq_pipeline.wdl \
  inputs_msrseq_perlmutter.json \
  perlmutter \
  --tag "msrseq_tRNA_run"
```

---

## Outputs

- `cleaned_ref_fasta` — header-cleaned mature tRNA reference FASTA
- `bt2_index_files` — Bowtie2 index files
- `all_tsvs` — binned base-level coverage TSV files across all samples
- `data_cleaned_zip` — aggregated intermediate tRNA modification dataset (`data_cleaned_5_HRPC.csv.zip`)
- `bedrmod_files` — array of all 12 final `bedRModv2` modification files:
  - `bs_rep1`, `bs_rep2`, `bs_rep3` — Bisulfite deletion modification files ($\Psi, m^7G, Gm$)
  - `cbh_rep1`, `cbh_rep2`, `cbh_rep3` — Cyanoborohydride mutation modification files ($ac^4C, f^5C$)
  - `ctrl_mut_rep1`, `ctrl_mut_rep2`, `ctrl_mut_rep3` — Untreated control mutation modification files ($m^1A, m^1G, m^2_2G, \dots$)
  - `ctrl_del_rep1`, `ctrl_del_rep2`, `ctrl_del_rep3` — Untreated control deletion modification files
