# tRNA WDL Pipeline

This directory contains the WDL workflows and helper scripts used for the tRNA analysis. The workflow chain is:

1. **Basecall + merge POD5 files** → `pipeline_SCATTER_jaws_tRNA_KJ.wdl`
2. **Realign merged unaligned BAMs to the transcriptome** → `tRNA_realign_JM.wdl`
3. **Merge aligned BAMs and run Modkit across thresholds** → `pipeline_tRNA_bam_merge_modkit_SCATTER.wdl`

## Workflows

| WDL file | Workflow namespace | Purpose |
|---|---|---|
| `pipeline_SCATTER_jaws_tRNA_KJ.wdl` | `ont_tRNA` | Basecall raw POD5 files (scattered), merge per-POD5 BAMs into one sample-level unaligned BAM, and produce a barcode table. |
| `tRNA_realign_JM.wdl` | `ont_tRNA` | Take the merged unaligned BAM from step 1 and align it to the tRNA transcriptome reference, then run Modkit pileup. |
| `pipeline_tRNA_bam_merge_modkit_SCATTER.wdl` | `ont_tRNA_merge_modkit_grid` | Merge all sample-level aligned BAMs from step 2 and run Modkit pileup across a grid of modification probability thresholds. |

## Helper scripts

| Script | Purpose |
|---|---|
| `generate_tRNA_full_pipeline_inputs.sh` | Generates a single input JSON for `pipeline_SCATTER_jaws_tRNA_KJ.wdl` |
| `generate_tRNA_full_pipeline_inputs_batch.sh` | Bulk generates input JSONs for HRP_A_020 / HRP_A_021 barcode samples |
| `submit_tRNA_realign_jobs.sh` | Submits all `tRNA_realign_JM.wdl` jobs for the realignment samples |

## Requirements

### Workflow engine
You need a WDL-compatible execution backend. The production runs used **JAWS** on Perlmutter, but the workflows can also be run with **Cromwell** or **miniwdl** locally.

Examples:
- JAWS: `jaws submit --no-cache <wdl> <inputs.json> perlmutter --tag <tag>`
- Cromwell: `java -jar cromwell.jar run <wdl> -i <inputs.json>`
- miniwdl: `miniwdl run <wdl> --input <inputs.json>`

### Docker
All tasks run inside Docker containers. The executor must be able to pull the images listed below. If you are running offline, pull them beforehand.

## Docker images

```bash
# Dorado basecaller
# Used by task DoradoBasecall in pipeline_SCATTER_jaws_tRNA_KJ.wdl
docker pull ontresearch/dorado@sha256:c8f356489fa8b44b31beba841b84d2879de2088e

# Seqtagger demultiplexing
# Used by task Seqtagger in pipeline_SCATTER_jaws_tRNA_KJ.wdl
docker pull lpryszcz/seqtagger:latest

# Minimap2 + samtools alignment utilities
# Used by MergeBams and MinimapTranscriptome tasks
docker pull nanozoo/minimap2:2.28--9e3bd01

# Modkit modification calling and pileup
# Used by ModkitPileup tasks
docker pull ontresearch/modkit@sha256:489d708a48c66368e5d1e118538e5dca68203a64
```

> **Note:** The WDLs reference images by digest (`sha256:...`) or by tag (`latest`). For reproducibility, prefer the digest-pinned versions shown above.

## Data flow

```
Raw POD5 files
    │
    ▼
pipeline_SCATTER_jaws_tRNA_KJ.wdl
    │
    ├──► barcode_table, barcode_pdfs (if barcoded)
    └──► merged_bam (unaligned, per-sample)
             │
             ▼
    tRNA_realign_JM.wdl
             │
             └──► transcriptome_bam + transcriptome_bed (per-sample)
                          │
                          ▼
        pipeline_tRNA_bam_merge_modkit_SCATTER.wdl
                          │
                          └──► merged_bam + modkit_beds (one per threshold)
```

## Input JSON preparation

### Step 1: Basecall + merge POD5 files

Use the provided generator script:

```bash
./generate_tRNA_full_pipeline_inputs.sh <POD5_DIR> <SAMPLE_ID>
```

Example:

```bash
./generate_tRNA_full_pipeline_inputs.sh \
  /path/to/pod5_dir \
  HRP_A_020_1_native_tRNA_001_FAST
```

This writes `tRNA/inputs_<SAMPLE_ID>.json` (relative to the parent `wdl-pipelines/` directory). The script hardcodes tRNA reference paths and the Dorado model directory; edit it if your environment differs.

For bulk HRP_A_020 / HRP_A_021 samples:

```bash
./generate_tRNA_full_pipeline_inputs_batch.sh
```

### Step 2: Transcriptome realignment

There is no generator for this workflow. The input JSON takes the merged unaligned BAM produced in step 1:

```json
{
  "ont_tRNA.unaligned_bam": "/path/to/HRP_A_020_1_native_tRNA_001_FAST.merged.bam",
  "ont_tRNA.sample_id": "HRP_A_020_1_native_tRNA_001_FAST",
  "ont_tRNA.ref_transcriptome": "/path/to/hg38.rtRNA.with_oligos.fa",
  "ont_tRNA.basecalling_type": "fast",
  "ont_tRNA.cpus": 32
}
```

### Step 3: Merge aligned BAMs + Modkit threshold grid

There is no generator for this workflow. Hand-write an input JSON listing all aligned BAMs from step 2:

```json
{
  "ont_tRNA_merge_modkit_grid.bamfiles": [
    "/path/to/sample1.transcriptome.aligned.sorted.bam",
    "/path/to/sample2.transcriptome.aligned.sorted.bam"
  ],
  "ont_tRNA_merge_modkit_grid.bamindices": [
    "/path/to/sample1.transcriptome.aligned.sorted.bam.bai",
    "/path/to/sample2.transcriptome.aligned.sorted.bam.bai"
  ],
  "ont_tRNA_merge_modkit_grid.sample_id": "merged_tRNA_sample",
  "ont_tRNA_merge_modkit_grid.reference": "/path/to/hg38.rtRNA.with_oligos.fa",
  "ont_tRNA_merge_modkit_grid.cpus": 12,
  "ont_tRNA_merge_modkit_grid.mod_thresholds": [
    0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99
  ]
}
```

> **Note:** In this WDL, `bamfiles` and `bamindices` are declared as `File` (not `Array[File]`), and the `MergeBams` task is commented out. Provide a single merged BAM or use an external merge step.

## Usage

### Step 1: Basecall + merge

```bash
jaws submit \
  --no-cache \
  pipeline_SCATTER_jaws_tRNA_KJ.wdl \
  tRNA/inputs_HRP_A_020_1_native_tRNA_001_FAST.json \
  perlmutter \
  --tag "HRP_A_020_1_native_tRNA_001_FAST_basecall"
```

### Step 2: Transcriptome realignment

```bash
jaws submit \
  --no-cache \
  tRNA_realign_JM.wdl \
  tRNA_REALGINMENT/realignment_all/inputs_HRP_A_020_1_native_tRNA_001_FAST.json \
  perlmutter \
  --tag "HRP_A_020_1_native_tRNA_001_FAST_realign"
```

Use `submit_tRNA_realign_jobs.sh` to submit all realignment jobs at once.

### Step 3: Merge + Modkit threshold grid

```bash
jaws submit \
  --no-cache \
  pipeline_tRNA_bam_merge_modkit_SCATTER.wdl \
  inputs_merge_modkit_grid_tRNA.json \
  perlmutter \
  --tag "tRNA_merge_grid"
```

## Outputs

### `pipeline_SCATTER_jaws_tRNA_KJ.wdl`
- `merged_bam`, `merged_bai` — merged unaligned basecalled BAM
- `barcode_table` — demux table (if barcoded)
- `barcode_pdfs` — demux PDFs (if barcoded)

### `tRNA_realign_JM.wdl`
- `transcriptome_bam`, `transcriptome_bai` — transcriptome-aligned BAM
- `transcriptome_bed`, `transcriptome_log` — Modkit transcriptome pileup

### `pipeline_tRNA_bam_merge_modkit_SCATTER.wdl`
- `modkit_beds` — array of Modkit pileup BEDs, one per threshold
- `modkit_logs` — array of Modkit logs, one per threshold
