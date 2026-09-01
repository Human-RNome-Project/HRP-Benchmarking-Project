# Biotype Composition Analysis (original script)

This folder contains the analysis script (`biotype_analysis.py`) used to compute and plot the RNA **biotype composition** (Ensembl/GENCODE `gene_type`) of long-read (Oxford Nanopore) sequencing samples for three RNA classes: **polyA RNA**, **rRNA** and **tRNA**. For each class it loads per-sample `featureCounts` read-count tables, merges them into a gene × sample matrix, normalizes to CPM (counts per million), attaches each gene's biotype from a GTF annotation, aggregates the composition per biotype, and draws stacked horizontal bar charts.

The script is written as a **notebook-style cell script** (`#%%` cell separators, intended for an interactive IDE such as VS Code or Jupyter). Input paths, sample lists and output paths are **hard-coded** in the source; it is not a command-line tool and takes no arguments. A parameterized, command-line version of the same analysis is documented in `README.md`.

## Enviroment setup

The script requires `pandas`, `polars`, `matplotlib` and the internal `dmode` package (imported as `import dmode` and used via `dmode.metagene_plot.prepare_gene_body_coverage`). `dmode` in turn pulls in additional dependencies (`numpy`, `scipy`, `scikit-learn`, `seaborn`, `umap-learn`, `pyranges`, `tqdm`).

### Using Conda
1. Ensure you have Conda (or Miniconda/Mamba) installed.
2. Create an environment with the required packages, for example:
   ```bash
   conda create -n biotype-analysis -c conda-forge \
       python pandas polars matplotlib numpy scipy scikit-learn seaborn pyranges tqdm
   conda activate biotype-analysis
   pip install umap-learn
   ```
3. Install the internal `dmode` package into the environment (it is not on Conda/PyPI):
   ```bash
   pip install -e /path/to/dmode
   ```

## Input data

The script reads two kinds of input, both referenced by absolute paths written directly in the source:

- **Count tables** — `featureCounts` output tables (`.tsv`), tab-separated with `#` comment header lines, a `Geneid` column and the read counts in the last column. One table per sample, plus a merged "combined" table per RNA class. The preprocessing that produces these (alignment with `minimap2`, sorting/filtering with `samtools`, quantification with `featureCounts`) is documented in the docstring at the top of the script.
- **GTF annotation** — a GENCODE/rRNA/tRNA GTF, parsed via `dmode` to obtain the gene-level `gene_type`, `gene_name` and `gene_id`.

The three RNA classes read from different base directories and GTFs, all hard-coded in the source (e.g. `/global/cfs/cdirs/m5243/...`, `/home/stefan/Synology/...`, `/mnt/data01/stpastore/...`). **These paths must be edited to match your system before running.**

## Usage

There are no command-line arguments. The script is meant to be edited and run cell by cell in an interactive session, or executed top to bottom:

```bash
python biotype_analysis.py
```

Before running, edit the hard-coded values described below so they point at your own data.

### Configuration (edit in the source)

| Variable / location | Section | Description |
|---|---|---|
| `base` | each of polyA / rRNA / tRNA blocks | Base directory the count tables and outputs are read from / written to. Reassigned at the start of each RNA-class block. |
| `HRP_A_*_df = pd.read_csv(...)` | each block | Individual per-sample table paths. One `read_csv` line per sample. |
| `combined_df = pd.read_csv(...)` | each block | Path to the merged/combined table for that RNA class. |
| `dfs = { ... }` | each block | Mapping of sample label → dataframe; controls which samples are included and their display names. |
| `prepare_gene_body_coverage("...")` | each block | Absolute path to the GTF annotation for that RNA class. |
| `sample_cols = [ ... ]` | each block | Which columns are plotted (e.g. `["combined"]` for the single-dataset figure, or the full sample list for the per-sample figure). |
| `.replace({ ... })` | rRNA / tRNA blocks | Biotype relabelling maps (e.g. `rRNA_5S` → `rRNA`, `tRNAAla_CGC` → `Ala_CGC`). |
| `gene_type_colors`, `transcript_colors` | top of file / tRNA block | Biotype → hex colour palettes. |
| `keep = ...` / `top_n` | each block | Threshold (`rel > 0.01`) or top-N selection controlling which biotypes are shown before collapsing the rest into `"other"`. |
| `plt.savefig("...")` / `fig.savefig("...")` | each block | Output figure paths and formats. |

### Examples

The script runs the same three blocks in sequence; there is no per-run parameterization. A typical workflow is:

1. Open `biotype_analysis.py` in VS Code / Jupyter.
2. In the **polyA** block, set `base`, the per-sample `read_csv` paths, `combined_df`, the `dfs` mapping and the GTF path; run the cells to produce the polyA table and figures.
3. Repeat for the **rRNA** block (note it additionally collapses `rRNA_5S/5.8S/18S/28S` into `rRNA`).
4. Repeat for the **tRNA** block (note it relabels tRNA/rRNA gene types to short amino-acid_anticodon names, uses a smaller font preset, and adds an anticodon/amino-acid distribution figure with `transcript_colors`).

## Output

Running all three blocks writes composition tables and figures (paths depend on the hard-coded `base` and `savefig` locations):

- **polyA:** `gene_type_composition_polyA_combined_only.tsv`; figures `Figure1_gene_type_composition_combined_only_CPMS.pdf`, `Figure1_gene_type_composition.png` / `.svg` / `_CPMS.pdf`.
- **rRNA:** `gene_type_composition_rRNA_combined_only.tsv`; figures `Figure1_rRNA_biotype_CPMS_single_dataset.pdf`, `Figure1_rRNA_biotype_CPMS.pdf`.
- **tRNA:** `gene_type_composition_tRNA_combined_only.tsv`; figures `Figure1E_tRNA_biotypes.pdf`, `Figure1_tRNA_biotype_CPMS.pdf`, `Figure1_tRNA_anticodon_aminoacid_distribution_CPMS.pdf`.

Composition tables give each biotype as a percentage of total signal per sample; figures are stacked horizontal bar charts of the relative composition, with low-abundance biotypes collapsed into `"other"`.

## Important Notes

- **Hard-coded, system-specific paths.** Input and output paths point at specific NERSC / Synology / local mounts and mix several machines; they must be edited before the script will run elsewhere.
- **Requires the internal `dmode` package.** GTF parsing goes through `dmode.metagene_plot.prepare_gene_body_coverage`; the script cannot run without `dmode` installed.
- **CPM normalization.** Each sample column is divided by its total counts / 1e6, and the composition table is then scaled to a percentage (`/ 1e6 * 100`).
- **Interactive by design.** The script calls `plt.show()` and is structured around `#%%` cells, so it expects an interactive/display environment rather than a headless run.
- **Biotype relabelling differs per class.** The rRNA block collapses rRNA subunits into a single `rRNA` biotype, while the tRNA block expands raw tRNA gene types into amino-acid_anticodon labels — check you are applying the intended map for your data.
