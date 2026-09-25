# MS-seq Consensus creation
This folder contains the command-line pipeline for building a consensus bedrmod file from multiple mass-spectrometry (OpenMS) sample bedrmod files. It loads each sample, applies score/frequency/mapping filters, merges fragment positions within each sample, and then combines all samples into a single consensus set of modification sites.
The thresholds used in our paper are the default values of the pipeline.

## Enviroment setup

### Using Docker
Start this from MS-seq main directory.
```bash
docker build -t hrp-ms-seq:latest .
```

## Input data
This analysis step depends on the outputs of the [`./NASE_rRNA/`](NASE_rRNA/) and [`./NASE_tRNA/`](NASE_tRNA/) workflows. 
Place all sample files you want to merge into one folder. The script automatically discovers every file ending in `.bed` or `.bedrmod` inside that folder and uses them all as samples.

## Usage

```bash
docker run --rm -v "$PWD":/data hrp-ms-seq:latest consensus \
    --input-folder /data/samples \
    --out-file /data/consensus.bed \
    [--min-samples INT] \
    [--q-score FLOAT] \
    [--freq FLOAT] \
    [--min-overlap INT] \
    [--unique-mapping | --no-unique-mapping] \
    [--tRNA]
```


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

**rRNA run** (unique mapping was applied in the paper):
```bash
docker run --rm -v "$PWD":/data hrp-ms-seq:latest consensus \
    --input-folder /data/rRNA_data \
    --out-file /data/consensus_rRNA.bed \
    --min-samples 1 \
    --q-score 0.05 \
    --freq 10 \
    --min-overlap 100 \
    --unique-mapping
```

**tRNA run** (unique mapping was not applied for the paper, base-end trimming applied):
```bash
docker run --rm -v "$PWD":/data hrp-ms-seq:latest consensus \
    --input-folder /data/tRNA_data \
    --out-file /data/consensus_tRNA.bed \
    --min-samples 1 \
    --q-score 0.05 \
    --freq 10 \
    --min-overlap 100 \
    --no-unique-mapping \
    --tRNA
```

## Output

Two files are written per run, based on `--out-file` (e.g. `test_consensus.bed`):

- **`<name>.bed`** — the final consensus bed file.
- **`<name>_statistics.bed`** — an information-rich version of the consensus before final reshaping, including per-site sample-overlap statistics.


## Importent Notes

- Coverage values are set to `0` in the output, since mass-spec fragment intensity is not equivalent to sequencing read coverage and is not comparable across fragments.
- Inosine (`I`) sites are dropped from each sample before merging. Since Inosine is only 1Da heavier than adenosine, its mass overlaps with the M+1 isotope peak of adenosine, making its detection with the used technology unreliable. 

