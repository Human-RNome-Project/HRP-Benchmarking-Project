# RNA Biotype Composition Analysis

This folder contains the command-line pipeline for computing and plotting the RNA **biotype composition** (Ensembl/GENCODE `gene_type`) of long-read (e.g. Oxford Nanopore) sequencing samples. It loads one `featureCounts` read-count table per sample, merges them into a gene × sample matrix, normalizes to CPM (counts per million), attaches each gene's biotype from a GTF annotation, optionally relabels biotypes, writes a per-biotype composition table, and draws stacked horizontal bar charts of the relative composition per sample.

## Enviroment setup

The pipeline requires `pandas`, `matplotlib` and `pyyaml`. `polars` is optional and only used for the internal `dmode` GTF-parsing path; without it a built-in GTF parser is used automatically.

### Using Conda
1. Ensure you have Conda (or Miniconda/Mamba) installed.
2. Create the environment by running:
   ```bash
   conda env create -f environment.yml
   ```
3. Activate the environment:
   ```bash
   conda activate biotype-analysis
   ```

The internal `dmode` package is not on Conda/PyPI and is therefore not listed in `environment.yml`. It is optional: install it separately if you want the `dmode` GTF-parsing path, or run with `--no-dmode` to force the built-in parser.

## Input data

Two kinds of input are needed:

- **Count tables** — one `featureCounts` output table (`.tsv`) per sample. Each table is tab-separated with `#` comment header lines, a single column-header row, a gene-id column (`Geneid` by default), and the **read counts in the last column**. An optional extra "combined"/merged table can be supplied and is used as the default sorting column.
- **GTF annotation** — the annotation used for quantification, from which the gene-level `gene_type`, `gene_name` and `gene_id` are read (`feature == "gene"`). GENCODE `gene_type` and Ensembl `gene_biotype` attributes are both supported.


## Read Alignment and Biotype Quantification

This section starts from unaligned BAM (uBAM) files produced by basecalling. Biotype analysis needs **two separate alignment strategies**:

- **rRNA and poly(A) RNA**: standard spliced long-read alignment to the genome.
- **tRNA**: short, highly structured reads need much more sensitive minimap2 settings. They are aligned first to a combined rRNA/tRNA reference. Reads that don't align there are then aligned to the genome.

### Requirements

- [samtools](http://www.htslib.org/)
- [minimap2](https://github.com/lh3/minimap2)
- [featureCounts](https://subread.sourceforge.net/) (Subread package)
- References:
  - `hg38.fa`: human genome
  - `hg38.rtRNA.with_oligos.fa`: rRNA + tRNA sequences, including adapter oligos
  - `combined_v49_rRNAs.gtf`: annotation including rRNA features

---

### 1. Alignment for rRNA and poly(A) libraries

```bash
samtools fastq -T "*" in.ubam \
  | minimap2 -y --MD -t 5 -ax splice-ont -uf -k14 $NERSC/references/hg38.fa - \
  | samtools sort -@ 2 \
  | samtools view --threads 2 -hb -F 260 \
  > out.aligned.sorted.bam
```

| Option | Purpose |
|---|---|
| `samtools fastq -T "*"` | Converts uBAM to FASTQ and keeps all read tags |
| `minimap2 -y` | Copies the tags into the alignment output |
| `-ax splice-ont -uf -k14` | Spliced alignment for ONT direct RNA, forward-strand only |
| `-F 260` | Removes unmapped (4) and secondary (256) alignments |

---

### 2. Alignment for tRNA libraries

#### 2.1 Align to the rRNA/tRNA reference

```bash
samtools fastq -T "*" in.ubam \
  | minimap2 -y -ax splice -t 5 -k 7 -w 3 -A 2 -B 1 -O1,32 -E2,0 -n 1 -m 13 -s 30 \
      --secondary=no --MD $NERSC/references/hg38.rtRNA.with_oligos.fa - \
  | samtools sort -@ 5 \
  | samtools view --threads 5 -hb \
  > out.rtRNA.all.sorted.bam
```

The small k-mer/window sizes (`-k 7 -w 3`), relaxed scoring and low chaining/score thresholds (`-m 13 -s 30`) let minimap2 seed and align short tRNA reads.

#### 2.2 Split aligned and unaligned reads

```bash
# Reads aligned to rRNA/tRNA
samtools view -hb -F 4 out.rtRNA.all.sorted.bam > out.rtRNA.aligned.sorted.bam

# Unaligned reads, passed on to genome alignment
samtools view -hb -f 4 out.rtRNA.all.sorted.bam > out.rtRNA.unaligned.sorted.bam
```

#### 2.3 Align the remaining reads to the genome

```bash
samtools fastq -T "*" out.rtRNA.unaligned.sorted.bam \
  | minimap2 -y --MD -t 5 -ax splice-ont -uf -k14 $NERSC/references/hg38.fa - \
  | samtools sort -@ 2 \
  | samtools view --threads 2 -hb -F 260 \
  > out.genome.aligned.sorted.bam
```

#### 2.4 Merge both alignments

```bash
samtools merge -@ 5 -f out.aligned.sorted.bam \
  out.rtRNA.aligned.sorted.bam \
  out.genome.aligned.sorted.bam
```

---

### 3. Quantification with featureCounts

Use the same command for every library type:

```bash
featureCounts -L -T 50 \
  -a combined_v49_rRNAs.gtf \
  -o out.aligned.sorted.tsv \
  out.aligned.sorted.bam
```

`-L` turns on long-read mode, which is required for nanopore data.

---

### Output overview

| File | Description |
|---|---|
| `out.aligned.sorted.bam` | Final alignments used for quantification |
| `out.rtRNA.aligned.sorted.bam` | tRNA workflow only: reads aligned to the rRNA/tRNA reference |
| `out.genome.aligned.sorted.bam` | tRNA workflow only: remaining reads aligned to the genome |
| `out.aligned.sorted.tsv` | featureCounts count table (plus a `.summary` file) |

## Usage

The tool can be run in two modes: a single run driven by command-line flags, or a batch run driven by a YAML config file.

**Single run:**
```bash
python biotype_analysis_executable.py \
    --name <label> \
    --counts <sample1.tsv> <sample2.tsv> ... \
    --gtf <annotation.gtf> \
    --output-dir <path> \
    [--names S1 S2 ...] \
    [--combined <merged.tsv>] [--combined-name combined] \
    [--relabel <relabel.yaml>] [--palette gene_type|transcript|<file>] \
    [--select fraction|top_n] [--threshold FLOAT] [--top-n INT] \
    [--sort-by NAME] [--combined-only] [--reverse-samples] \
    [--hide-yticks] [--figsize W H] [--font-preset small|large] \
    [--no-dmode] [--show] [-v]
```

**Batch run (recommended for reproducibility):**
```bash
python biotype_analysis_executable.py --config config.example.yaml
```

Run `python biotype_analysis_executable.py --help` at any time to see this list with current defaults.

### Arguments

| Argument | Default | Description |
|---|---|---|
| `--config` | *(none)* | YAML config describing one or more analyses (batch mode). If given, the single-run flags below are ignored. |
| `--counts` | *(required in single-run)* | One or more `featureCounts` `.tsv` tables, one per sample. |
| `--names` | *(derived)* | Sample names, in the same order as `--counts`. If omitted, names are derived from the file names. |
| `--combined` | *(none)* | Optional merged/combined table, added as an extra column and used as the default `--sort-by` column. |
| `--combined-name` | `combined` | Column name for the `--combined` table. |
| `--gtf` | *(required in single-run)* | GTF annotation file. |
| `--id-column` | `Geneid` | Gene-id column in the count tables. |
| `--output-dir` / `-o` | `biotype_results` | Directory for the composition table and figures (created if needed). |
| `--name` | `analysis` | Short label used in output file names and as the default y-axis label. |
| `--relabel` | *(none)* | YAML/JSON file with an `{old_biotype: new_biotype}` map (e.g. collapse `rRNA_5S` → `rRNA`). |
| `--palette` | `gene_type` | Built-in palette name (`gene_type` or `transcript`), or a YAML/JSON file mapping biotypes to hex colours. |
| `--select` | `fraction` | Biotype selection: keep those above `--threshold` in any sample (`fraction`), or the `--top-n` most abundant (`top_n`). |
| `--threshold` | `0.01` | Relative-share threshold for `--select fraction`. |
| `--top-n` | `8` | Number of biotypes to keep for `--select top_n`. |
| `--sort-by` | `combined` | Sample column to sort the composition table by (falls back to total signal if absent). |
| `--y-label` | *(= `--name`)* | Y-axis label. |
| `--x-label` | `Relative composition (CPM / 1M)` | X-axis label. |
| `--hide-yticks` | off | Hide per-sample y tick labels. |
| `--figsize` | `16 8` | Figure size in inches (width height). |
| `--reverse-samples` | off | Reverse the order samples are stacked in the plot. |
| `--combined-only` | off | Plot only the `--combined` column instead of all samples. |
| `--font-preset` | `large` | Matplotlib font-size preset (`large` or `small`). |
| `--no-dmode` | off | Force the built-in GTF parser even if `dmode` is installed. |
| `--show` | off | Display figures interactively (default: save only, headless-safe). |
| `-v` / `--verbose` | off | Verbose logging. |

For batch mode, each analysis in the config accepts the same options as keys, plus a `plots` list (one entry per figure) with keys `sample_cols` (list, or `null` = all non-combined samples), `y_label`, `hide_yticks`, `figsize`, `method`, `threshold`, `top_n`, `reverse` and `outfile`. Values under `defaults` are inherited by every analysis and can be overridden per analysis. See `config.example.yaml`.

### Examples

**Single polyA run** (explicit sample names, combined column, default palette):
```bash
python biotype_analysis_executable.py \
    --name polyA \
    --counts S001.tsv S002.tsv S003.tsv \
    --names HRP_A_001_1 HRP_A_002_1 HRP_A_003_1 \
    --combined polyA_native_merged.tsv \
    --gtf gencode.v49.basic.annotation.gtf \
    --output-dir results/
```

**tRNA run** (built-in parser, relabel map, top-8 biotypes, small fonts):
```bash
python biotype_analysis_executable.py \
    --name tRNA \
    --counts tRNA_*.tsv \
    --combined combined_dataset_final_merge.tsv \
    --gtf combined_v49_rRNAs_tRNAs.gtf \
    --relabel trna_relabel.yaml \
    --palette transcript \
    --select top_n --top-n 8 \
    --font-preset small \
    --no-dmode \
    --output-dir results/tRNA/
```

**Batch run** reproducing the polyA, rRNA and tRNA analyses at once:
```bash
python biotype_analysis_executable.py --config config.example.yaml
```

## Output

For each analysis, two kinds of files are written into `--output-dir`:

- **`gene_type_composition_<name>.tsv`** — the biotype composition table: rows are biotypes, columns are samples, values are the percentage of total CPM signal contributed by that biotype (each column sums to 100). The `table_name` config key overrides this name.
- **One or more figures** — stacked horizontal bar charts of the relative biotype composition (one bar per sample). File names and formats follow `--name` in single-run mode, or the `outfile` of each plot spec in batch mode; the format (`.pdf`, `.png`, `.svg`, …) is taken from the file extension.

The tool is also importable as a module: `load_count_matrix`, `compute_cpm`, `load_gene_annotation`, `annotate_counts`, `composition_table`, `plot_composition` and `run_analysis` can be reused directly.

## Important Notes

- **CPM normalization.** Each sample is scaled by its own total counts / 1e6, so columns are comparable across sequencing depths. The composition table then expresses each biotype as a percentage of the per-sample total.
- **GTF parsing is portable.** `dmode.utility.gtf_to_df` is used when `dmode` is importable; otherwise a built-in parser reads gene-level `gene_type`/`gene_biotype`, `gene_name` and `gene_id`. Use `--no-dmode` to always use the built-in parser.
- **Missing genes are treated as zero.** Count tables are merged with an outer join on the gene id; a gene absent from a sample contributes 0 counts.
- **Headless by default.** Unless `--show` is passed, a non-interactive matplotlib backend is used so the pipeline runs over SSH / on a cluster with no display.
- **Behavioural difference from the original script.** The composition table normalizes each column by its own total (identical to the original `sum / 1e6 * 100` for CPM input, but correct when some genes are missing from a sample), and the leftover `key_0` merge artifact is no longer written.
