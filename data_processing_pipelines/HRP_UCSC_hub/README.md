# Human RNome Project - UCSC Genome Browser Hub Code (Fig. 14)

This repository contains the code used to generate the UCSC Genome Browser hub
for the Human RNome Project. The hub visualizes RNA modifications across poly(A)
RNA (hg38), ribosomal RNA (rRNA), and transfer RNA (tRNA), detected by short-read
sequencing (Illumina), long-read sequencing (Oxford Nanopore), and mass spectrometry.

The input data and the ready-to-load hub files are deposited separately at the Human RNome Project Data Portal. 
This repository contains only the code.

---

## What is where

| Location | Contents |
|---|---|
| This repository (GitHub) | Pipeline code to generate the hub from input data |
| Human RNome Project Data Portal — `input_data/` | Modification BED files, sequence alignments, references, variant calls |
| Human RNome Project Data Portal — `ucsc_hub/` | The complete, ready-to-load UCSC hub (track files + configuration) |

To reproduce the hub, download `input_data/` from the Human RNome Project Data Portal and run
`scripts/pipeline.py` (see below). To browse the tracks without reproducing them,
load `ucsc_hub/` from the Human RNome Project Data Portal directly.

---

## Scripts

| Script | Purpose |
|---|---|
| `scripts/pipeline.py` | End-to-end pipeline: splits modification BED files by RNA type, converts to bigBed, generates hub configuration and track description pages |
| `scripts/make_decorator.py` | Builds the modification decorator track overlaid on the RNA sequence track |
| `scripts/BedPyLift.py` | Transcriptome-to-genome coordinate lifting utility |

---

## Running the pipeline

### Requirements

- Python 3.8+
- UCSC command-line tools (`bedToBigBed`, `faToTwoBit`) placed in `utils/`
- `openpyxl` (optional, for Excel manifest): `pip install openpyxl`

### Setup

Download `input_data/` from the Human RNome Project Data Portal and place it at the root of this
repository. Then update the input file paths at the top of `scripts/pipeline.py`:


### Run

```bash
python scripts/pipeline.py
```

This produces a complete `ucsc_hub/` directory. To rebuild only the modification
decorator track:

```bash
python scripts/make_decorator.py
```
