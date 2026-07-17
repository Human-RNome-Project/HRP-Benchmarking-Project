# Human RNome Project UCSC Genome Browser Hub (Fig. 14)

UCSC Genome Browser assembly hub visualizing RNA modifications across the human
RNome: poly(A) RNA (hg38), ribosomal RNA (rRNA), and transfer RNA (tRNA).
Modifications are detected by three modalities — short-read sequencing (SRS),
long-read sequencing (LRS), and mass spectrometry (MS). Each modality is
rendered as a single flat track per RNA type, colored by modification type with
saturation proportional to frequency. Each RNA type also has a tiered
cross-platform consensus track (tier1/tier2 confidence). The poly(A) RNA
assembly additionally includes a native RNA sequence track (bigPSL) with
consensus modification sites overlaid as decorators.

The processed data files backing these tracks are deposited at NERSC (see manuscript).

---

## Contents

```
├── README.md
├── scripts/
│   ├── pipeline.py          # End-to-end hub generation pipeline
│   ├── make_decorator.py    # Builds modification decorator bigBed for the sequence track
│   └── BedPyLift.py         # Transcriptome → genome coordinate lifting utility
└── ucsc_hub/                # UCSC assembly hub configuration (no data files)
    ├── hub.txt
    ├── genomes.txt
    ├── hubDescription.html
    ├── hg38_polyA-RNA/      # poly(A) RNA — trackDb + track description pages
    ├── rrna/                # rRNA custom assembly — trackDb + description pages
    └── trna/                # tRNA custom assembly — trackDb + description pages
```

---

## Pipeline

### Overview (`scripts/pipeline.py`)

The pipeline takes pre-merged, genome-mapped BED files (one per modality) and
produces a complete UCSC assembly hub:

1. **Split** — split combined BED files by RNA type
   (chrom containing "rRNA" → rRNA, "tRNA" → tRNA, else → poly(A) RNA hg38)
2. **Convert** — process split BED files → bigBed (per platform), plus tiered
   consensus beds → bigBed (per RNA type)
3. **HTML** — generate track description pages from `RNA_modifications_manifest.tsv`
4. **Manifest** — generate Excel manifest of all tracks
5. **Hub config** — write hub.txt, genomes.txt, trackDb.txt for all assemblies
6. **Decorator** — build modification decorator bigBed for the RNA sequence track

### Input files

Update `COMBINED_FILES` in `scripts/pipeline.py` to point to your input data:

| Key | Expected file | Description |
|---|---|---|
| `SRS` | `final_bedRmods/SRS_modifications.bed` | Short-read sequencing consensus calls |
| `LRS` | `final_bedRmods/LRS_modifications.bed` | Long-read sequencing consensus calls |
| `MS`  | `final_bedRmods/MS_modifications.bed`  | Mass spectrometry consensus calls |

These are pre-merged bedRmod-format files (13 columns) with all RNA types
combined. The pipeline splits them by chromosome name.

The tiered consensus input (`TIERED_COMBINED_FILE`) should be a 6-column BED
with header `chr/start/end/name/tier/strand`.

### Requirements

- Python 3.8+
- UCSC command-line tools (in `utils/`): `bedToBigBed`, `faToTwoBit`
- `openpyxl` (optional, for Excel manifest): `pip install openpyxl`

### Running

```bash
# From the repository root:
python scripts/pipeline.py
```

To rebuild only the modification decorator:
```bash
python scripts/make_decorator.py
```

### Configuration

Before running, update these constants at the top of `scripts/pipeline.py`:

| Constant | Purpose |
|---|---|
| `COMBINED_FILES` | Paths to input BED files |
| `TIERED_COMBINED_FILE` | Path to combined tiered consensus BED |
| `GITHUB_USER`, `GITHUB_REPO` | GitHub repository for hub hosting |
| `BIGPSL_GCS_URL` | Cloud storage URL for the bigPSL sequence file |
| `PAPER_URL` | Published DOI |

---

## UCSC Hub

The `ucsc_hub/` directory contains the complete hub configuration:
- `hub.txt` / `genomes.txt` — hub metadata and genome assembly registrations
- Per-assembly `trackDb.txt` — track stanzas for all modification and sequence tracks
- HTML pages — description pages for every track (required by UCSC)

The hub covers three assemblies:
- **hg38** — poly(A) RNA modifications + native RNA sequence track + NA12878 WGS variants
- **hs_rRNA** — custom rRNA assembly; rRNA modifications by SRS, LRS, and MS
- **hs_tRNA** — custom tRNA assembly; tRNA modifications by SRS, LRS, and MS

Large binary data files (bigBed, bigPSL, VCF, 2bit) are not included here;
they are deposited at NERSC and hosted for UCSC browsing on cloud storage.

---

## Data availability

Processed data files are deposited at NERSC under the Human RNome Project data directory.
See the manuscript for accession details.
