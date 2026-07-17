# polyA/mRNA WDL Pipeline

This directory contains the WDL workflows and input generator for processing ONT direct-RNA polyA/mRNA samples.

## Workflows

| WDL file | Workflow namespace | Purpose |
|---|---|---|
| `pipeline_SCATTER_jaws.wdl` | `ont_mRNA_pilot` | Full pipeline from POD5 to aligned BAMs and Modkit pileups: Seqtagger demux → Dorado basecalling → BAM merge → Minimap2 genome + transcriptome alignment → NanoComp QC → Modkit pileup |
| `pipeline_polyA_bam_merge_modkit_SCATTER.wdl` | `ont_mRNA_pilot_merge_modkit_grid` | Merge multiple aligned BAMs and run Modkit pileup across a grid of modification probability thresholds |

## Requirements

### Workflow engine
You need a WDL-compatible execution backend. The production runs used **JAWS** on Perlmutter, but the workflows can also be run with **Cromwell** or **miniwdl** locally.

Examples:
- JAWS: `jaws submit --no-cache <wdl> <inputs.json> perlmutter --tag <tag>`
- Cromwell: `java -jar cromwell.jar run <wdl> -i <inputs.json>`
- miniwdl: `miniwdl run <wdl> --input <inputs.json>`

### Docker
All tasks run inside Docker containers. The executor must be able to pull the images listed below. If you are running offline, pull them beforehand.

### Input data
- For `pipeline_SCATTER_jaws.wdl`: POD5 files, reference genome FASTA, transcriptome FASTA, Dorado model directory
- For `pipeline_polyA_bam_merge_modkit_SCATTER.wdl`: aligned/sorted BAM files + indices, reference genome FASTA, array of Modkit probability thresholds

## Docker images

The workflows use the following container images. Pull commands are provided for convenience.

```bash
# Dorado basecaller
# Image used by task DoradoBasecall in pipeline_SCATTER_jaws.wdl
docker pull ontresearch/dorado@sha256:c8f356489fa8b44b31beba841b84d2879de2088e

# Seqtagger demultiplexing
# Image used by task Seqtagger in pipeline_SCATTER_jaws.wdl
docker pull lpryszcz/seqtagger:latest

# Minimap2 + samtools alignment utilities
# Used by MergeBams, MinimapGenome, MinimapTranscriptome tasks
docker pull nanozoo/minimap2:2.28--9e3bd01

# NanoComp QC
# Used by NanoCompQC tasks
docker pull luxendr13/nanocomp:0.6.0

# Modkit modification calling and pileup
# Used by ModkitPileup tasks
docker pull ontresearch/modkit@sha256:489d708a48c66368e5d1e118538e5dca68203a64
```

> **Note:** The WDLs reference images by digest (`sha256:...`) or by tag (`latest`). For reproducibility, prefer the digest-pinned versions shown above.

## Input JSON preparation

### Full polyA/mRNA pipeline

Use the provided generator script:

```bash
./generate_inputs.sh
```

This writes a single-sample input JSON under `polyA_RNA/inputs_<sample_id>.json`. Edit the variables at the top of `generate_inputs.sh` (POD5 directory, sample ID, reference paths, GPU string, etc.) before running.

### Merge + threshold-grid Modkit

There is no generator script for this workflow. Hand-write an input JSON using this template:

```json
{
  "ont_mRNA_pilot_merge_modkit_grid.bamfiles": [
    "/path/to/sample1.aligned.sorted.bam",
    "/path/to/sample2.aligned.sorted.bam"
  ],
  "ont_mRNA_pilot_merge_modkit_grid.bamindices": [
    "/path/to/sample1.aligned.sorted.bam.bai",
    "/path/to/sample2.aligned.sorted.bam.bai"
  ],
  "ont_mRNA_pilot_merge_modkit_grid.sample_id": "my_merged_sample",
  "ont_mRNA_pilot_merge_modkit_grid.reference": "/path/to/GRCh38.primary_assembly.genome.fa",
  "ont_mRNA_pilot_merge_modkit_grid.cpus": 12,
  "ont_mRNA_pilot_merge_modkit_grid.mod_thresholds": [
    0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99
  ]
}
```

The example above uses 15 modification probability thresholds.

## Usage

### Run the full pipeline from POD5

```bash
jaws submit \
  --no-cache \
  pipeline_SCATTER_jaws.wdl \
  polyA_RNA/inputs_polyA_SUP_<sample>.json \
  perlmutter \
  --tag "<sample>_polyA_run"
```

### Run the merge + threshold-grid workflow

```bash
jaws submit \
  --no-cache \
  pipeline_polyA_bam_merge_modkit_SCATTER.wdl \
  inputs_merge_modkit_grid.json \
  perlmutter \
  --tag "my_merge_grid"
```

## Outputs

### `pipeline_SCATTER_jaws.wdl`
- `merged_bam`, `merged_bai` — merged basecalled BAM
- `genome_bam`, `genome_bai` — genome-aligned BAM
- `transcriptome_bam`, `transcriptome_bai` — transcriptome-aligned BAM
- `genome_bed`, `genome_log` — Modkit genome pileup
- `transcriptome_bed`, `transcriptome_log` — Modkit transcriptome pileup
- `nanocomp_report` — NanoComp QC report tarball

### `pipeline_polyA_bam_merge_modkit_SCATTER.wdl`
- `merged_bam`, `merged_bai` — merged input BAMs
- `modkit_beds` — array of Modkit pileup BEDs, one per threshold
- `modkit_logs` — array of Modkit logs, one per threshold
- `nanocomp_report` — NanoComp QC report tarball
