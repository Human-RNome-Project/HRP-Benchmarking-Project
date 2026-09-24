# bedRMod Liftover & Merge Pipeline

This directory contains `liftover_to_BedRMod.py`, a command-line tool that harmonises and merges one or more [bedRMod](https://github.com/anmol25081/bedrmod) files coming from different RNA types (mRNA, tRNA, rRNA) and sequencing technologies (Illumina, ONT, MS) into a single, standardised bedRMod file. It loads each input file (auto-detecting the column layout), normalises modification names to a canonical vocabulary, remaps chromosome/reference names to UCSC-style identifiers, optionally corrects tRNA coordinates for ONT adapter offsets, concatenates all files, appends single-letter modification codes and numeric modification IDs, merges the per-file metadata headers, and writes one combined bedRMod file.

## Enviroment setup

The tool requires a Python 3.9+ environment with `polars`, `numpy` and `tqdm` (`argparse` is part of the standard library). A Conda environment file, `liftover_environment.yml`, is provided.

### Using Conda
1. Ensure you have Conda (or Miniconda/Mamba) installed.
2. Create the environment by running:
   ```bash
   conda env create -f liftover_environment.yml
   ```
3. Activate the environment:
   ```bash
   conda activate bedrmod
   ```

Alternatively, with pip: `pip install polars numpy tqdm`.

## Input data

The tool takes one or more **bedRMod** files (`.bed`/`.bedrmod`), one per RNA type / condition, all referring to the same assembly. Each file is tab-separated: comment lines are `#key=value` metadata (parsed and carried into the merged header), and data rows follow the bedRMod column layout. The loader auto-detects the column count (layouts with 11–28 columns are recognised) and keeps the first 11 standard columns: `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`, `thickStart`, `thickEnd`, `itemRgb`, `coverage`, `frequency`. Frequency values on a 0–1 scale (max ≤ 1.0) are automatically rescaled to a 0–100 percentage.

No external genomic references are required — the chromosome and modification mapping tables (`_CHROM_MAP`, `_MOD_MAP`, `_MOD_TO_SYMBOL`, `_MOD_TO_ID`) are hard-coded in the script and cover GRCh37/GRCh38 primary chromosomes, scaffolds/ALT contigs, and human rRNA/tRNA reference names.

## Usage

```bash
python liftover_to_BedRMod.py \
    -b <file1.bed> [<file2.bed> ...] \
    -r <rna_type1> [<rna_type2> ...] \
    -t <Illumina|ONT|MS> \
    -o <output.bedrmod> \
    [-s <tRNA_adapter_offset>]
```

The number of `--bed_files` must match the number of `--rna_types`, given in the same order. Run `python liftover_to_BedRMod.py --help` at any time to see this list with current defaults.

### Arguments

| Argument | Default | Description |
|---|---|---|
| `-b`, `--bed_files` | *(required)* | One or more input bedRMod files to merge (space-separated). |
| `-r`, `--rna_types` | *(required)* | RNA type label for each input file, in the same order (e.g. `mRNA tRNA rRNA`). Written into the merged header and used to decide whether the tRNA offset applies. |
| `-t`, `--technology` | *(none)* | Sequencing technology in use: `Illumina`, `ONT`, or `MS`. Written to the `sequencing_plattform` header field. |
| `-s`, `--tRNA_adapter_offset` | `0` | Integer coordinate offset added to `chromStart`/`chromEnd` of tRNA files to correct for adapter-inclusive ONT alignment (a typical value is `-3`). Applied **only** when the file's `--rna_types` label is `tRNA` **and** `--technology` is `ONT`. |
| `-o`, `--output_file` | *(none)* | Destination path for the merged bedRMod output. |

### Examples

**ONT run** (mRNA + tRNA + rRNA, tRNA coordinates shifted back by 3 to undo adapter-inclusive alignment):
```bash
python liftover_to_BedRMod.py \
    -b mRNA.bed tRNA.bed rRNA.bed \
    -r mRNA tRNA rRNA \
    -t ONT \
    -s -3 \
    -o combined.bedrmod
```

**Illumina run** (no coordinate offset needed):
```bash
python liftover_to_BedRMod.py \
    -b illumina_mRNA.bed \
    -r mRNA \
    -t Illumina \
    -o illumina_combined.bedrmod
```

## Output

One combined **bedRMod** file is written to `--output_file`:

- **`<output_file>`** — a tab-separated `bedRModv2` file with a `#key=value` metadata header, then a `#`-prefixed column-header line, then the merged data rows. It contains all input rows with standardised modification and chromosome names, plus two appended columns: `single_letter_code` (a Modomics-style symbol from `_MOD_TO_SYMBOL`) and `mod_id` (a numeric modification ID from `_MOD_TO_ID`). The header aggregates each input file's `assembly`, `annotation_source`, `basecalling`, `bioinformatics_workflow`, `experiment` and `external_source` fields, each tagged by RNA type (`<rna_type>:<value>;`), and records the technology in `sequencing_plattform`.

## Important Notes

- **`--bed_files` and `--rna_types` are positional partners.** They must have the same length and the same order; each RNA-type label describes the file at the same position.
- **The tRNA offset is conditional.** `--tRNA_adapter_offset` is applied to a file only when its RNA-type label is exactly `tRNA` and `--technology` is `ONT`; otherwise it is ignored. It defaults to `0` (no shift) and should be given as a (typically negative) integer such as `-3`.
- **Frequency is normalised to a percentage.** If a file's `frequency` column is on a 0–1 scale (maximum ≤ 1.0), all values are multiplied by 100; rows with `frequency` ≤ 0 or an empty `name` are dropped.
- **Modification names must be resolvable.** Names are first standardised via `_MOD_MAP` (unmapped names are passed through), then converted to a symbol and a numeric ID via `_MOD_TO_SYMBOL` / `_MOD_TO_ID`. Any post-standardisation name absent from those tables is printed and will stop the run — extend the mapping dictionaries in the script to add new modifications.
- **Chromosome names must be in `_CHROM_MAP`.** Chromosome remapping is strict: names not present in the table are printed and will stop the run. The table already covers GRCh37/GRCh38 primary chromosomes, scaffolds/ALT contigs and human rRNA/tRNA reference names.
- **Only the first 11 columns are kept.** Extra columns from wider bedRMod layouts (12–28 columns) are recognised for parsing but not carried into the output.
