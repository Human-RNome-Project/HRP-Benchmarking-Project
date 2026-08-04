# MS-seq Consensus creation
This folder contains the command-line pipeline for building a consensus bedrmod file from multiple mass-spectrometry (OpenMS) sample bedrmod files. It loads each sample, applies score/frequency/mapping filters, merges fragment positions within each sample, and then combines all samples into a single consensus set of modification sites.

## Enviroment setup
The pipeline requires pandas.

### Using Conda
1. Ensure you have Conda (or Miniconda/Mamba) installed.
2. Create the environment by running:
   ```bash
   conda env create -f environment.yml
   ```
3. Activate the environment:
   ```bash
   conda activate MS-seq_consensus_creation
   ```

## Input data

Put all sample files for a single consensus run into one folder. The script automatically discovers every file ending in `.bed` or `.bedrmod` inside that folder and uses them all as samples.

## Usage

```bash
python consensus_creation_main.py \
    --input-folder <path> \
    --out-file <path> \
    [--min-samples INT] \
    [--q-score FLOAT] \
    [--freq FLOAT] \
    [--min-overlap INT] \
    [--unique-mapping | --no-unique-mapping] \
    [--tRNA]
```

Run `python consensus_creation_main.py --help` at any time to see this list with current defaults.

### Arguments

| Argument | Default | Description |
|---|---|---|
| `--input-folder` | *(required)* | Folder containing sample files; all `.bed`/`.bedrmod` files in it are used as samples. |
| `--out-file` | *(required)* | Path to the output consensus bedrmod file. |
| `--min-samples` | `1` | Minimum number of samples covering a modification site required for accepting it into the consensus. |
| `--q-score` | `0.05` | Score/q-value threshold used for filtering. |
| `--freq` | `10` | Frequency threshold used for filtering. |
| `--min-overlap` | `100` | Minimum percentage of samples that must have a modification at a specific site for it to be included in the final consensus. |
| `--unique-mapping` / `--no-unique-mapping` | `--unique-mapping` | Whether to require unique mapping during filtering. |
| `--tRNA` | off | Treat samples as tRNA data (applies `remove_base_ends` before merging). |

### Examples

**rRNA run** (unique mapping required, default settings):
```bash
python consensus_creation_main.py \
    --input-folder /path/to/rRNA_data/ \
    --out-file test_consensus_rRNA.bed \
    --min-samples 1 --q-score 0.05 --freq 10 --min-overlap 100 --unique-mapping
```

**tRNA run** (unique mapping not required, base-end trimming applied):
```bash
python consensus_creation_main.py \
    --input-folder /path/to/tRNA_data/ \
    --out-file test_consensus_tRNA.bed \
    --min-samples 1 --q-score 0.05 --freq 10 --min-overlap 100 --no-unique-mapping --tRNA
```

## Output

Two files are written per run, based on `--out-file` (e.g. `test_consensus.bed`):

- **`<name>.bed`** — the final consensus bedrmod file, reshaped into standard BED format.
- **`<name>_statistics.bed`** — an information-rich version of the consensus before final reshaping, including per-site sample-overlap statistics.

Both files carry forward the original comment header (`#key=value`) lines from the input samples, followed by a column-header line and tab-separated BED data.

## Importent Notes

- Coverage values are set to `0` in the output, since mass-spec fragment intensity is not equivalent to sequencing read coverage and is not comparable across fragments.
- Inosine (`I`) sites are dropped from each sample before merging. Since Inosine is only 1Da heavier than adenosine, its mass overlaps with the M+1 isotope peak of adenosine, making its detection with the used technology unreliable. 

