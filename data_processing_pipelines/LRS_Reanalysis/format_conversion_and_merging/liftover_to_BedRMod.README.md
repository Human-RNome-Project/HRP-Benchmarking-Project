# bedRMod Liftover & Merge Pipeline

This directory contains: `liftover_to_BedRMod.py`, a command-line tool that harmonises and merges one or more [bedRMod](https://github.com/anmol25081/bedrmod) files coming from different RNA types (mRNA, tRNA, rRNA) and sequencing technologies (Illumina, ONT, MS) into a single, standardised bedRMod file. It normalises modification names, remaps chromosome/reference names to UCSC-style identifiers, optionally corrects tRNA coordinates for ONT adapter offsets, appends single-letter modification codes and numeric modification IDs, and merges the per-file metadata headers.

## Scripts/Workflows

| File | Purpose |
|------|---------|
| `liftover_to_BedRMod.py` | Loads each input bedRMod file, standardises modification names (`_MOD_MAP`) and chromosome names (`_CHROM_MAP`), applies an optional coordinate offset to ONT tRNA files, vertically concatenates all files, appends `single_letter_code` (`_MOD_TO_SYMBOL`) and `mod_id` (`_MOD_TO_ID`) columns, merges the comment headers, and writes one combined bedRMod file. |

## Requirements & Setup

### Environment setup

To run the script you will need a Python 3.9+ environment with the required dependencies (`polars`, `numpy`, `tqdm`; `argparse` is part of the standard library). An `environment.yml` file is provided for Conda users and a `requirements.txt` for pip users.

**Option A: Using Conda**

```bash
conda env create -f environment.yml
conda activate bedrmod
```

**Option B: Using Pip**

```bash
python -m venv .venv
source .venv/bin/activate
pip install polars numpy tqdm
```

## Input data

### Input files

| Input File | Description | File Type/Format | Script Usage |
|------------|-------------|------------------|--------------|
| `*.bed` (bedRMod files) | One bedRMod file per RNA type / condition. Comment lines are `#key=value` metadata; data rows follow the bedRMod column layout. The loader auto-detects the column count (supports 11–28 columns) and reads the first 11 standard columns. Frequency values ≤ 1.0 are automatically rescaled to a 0–100 percentage. | Tab-separated bedRMod (`.bed`) | `liftover_to_BedRMod.py` |

### Command-line arguments

| Flag | Argument | Description |
|------|----------|-------------|
| `-b`, `--bed_files` | one or more paths | Input bedRMod files to merge (space-separated). |
| `-r`, `--rna_types` | one or more labels | RNA type label for each input file, in the same order (e.g. `mRNA tRNA rRNA`). |
| `-t`, `--technology` | string | Sequencing technology in use: `Illumina`, `ONT`, or `MS`. Written to the `sequencing_plattform` header field. |
| `-s`, `--tRNA_adapter_offset` | int (default `0`) | Coordinate offset added to tRNA files to correct for adapter-inclusive alignment when using ONT (e.g. `-3`). Only applied when `--rna_types` is `tRNA` **and** `--technology` is `ONT`. |
| `-o`, `--output_file` | path | Destination path for the merged bedRMod output. |

### References/Accessory files

No external genomic references are required — the chromosome and modification mapping tables (`_CHROM_MAP`, `_MOD_MAP`, `_MOD_TO_SYMBOL`, `_MOD_TO_ID`) are hardcoded in the script and cover GRCh37/GRCh38 primary chromosomes, scaffolds/ALT contigs, and human rRNA/tRNA reference names.

## Usage/Step-wise Execution

### Step 1: Run the liftover/merge

```bash
python liftover_to_BedRMod.py \
  -b mRNA.bed tRNA.bed rRNA.bed \
  -r mRNA tRNA rRNA \
  -t ONT \
  -s -3 \
  -o combined.bedrmod
```

The number of `--bed_files` must match the number of `--rna_types`, given in the same order. For ONT tRNA inputs, pass `--tRNA_adapter_offset` to shift coordinates back into the reference frame.

## Outputs

| Output File | Format | Script | Description |
|-------------|--------|--------|-------------|
| `<output_file>` (e.g. `combined.bedrmod`) | Tab-separated bedRMod (`bedRModv2`) with `#key=value` metadata header | `liftover_to_BedRMod.py` | Single merged file containing all input rows with standardised modification and chromosome names, plus two appended columns: `single_letter_code` (Modomics-style symbol) and `mod_id` (numeric modification ID). The header aggregates each input file's `assembly`, `annotation_source`, `basecalling`, `bioinformatics_workflow`, `experiment`, and `external_source` fields, tagged by RNA type. |
