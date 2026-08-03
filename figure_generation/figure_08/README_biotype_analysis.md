# biotype_analysis.py

Quantify and visualize the RNA biotype composition of Oxford Nanopore (ONT) direct-RNA
sequencing samples. The script reads per-gene count tables produced by `featureCounts`,
normalizes them to CPM, annotates each gene with its biotype from a GENCODE GTF, and
produces stacked horizontal bar charts of relative composition for three RNA classes:
**polyA RNA**, **rRNA**, and **tRNA**.

## Overview

The script is organized as a set of `#%%` cells (intended for interactive execution in
VS Code / Spyder / Jupyter) and is split into three largely parallel blocks:

1. **polyA RNA** — genome-aligned samples annotated against `gencode.v49.basic.annotation.gtf`.
2. **rRNA** — reads aligned to an rRNA reference; the individual 5S / 5.8S / 18S / 28S
   biotypes are collapsed into a single `rRNA` category.
3. **tRNA** — reads aligned to an rRNA+tRNA reference, with a per-isodecoder breakdown
   (amino-acid / anticodon) for the cytoplasmic and mitochondrial tRNA families.

Each block follows the same pattern: load count tables → merge into a single matrix →
CPM-normalize → join biotype annotation → group by biotype → keep biotypes above a 1%
relative share (collapsing the rest into `other`) → plot and export figures.

## Upstream preprocessing

The docstring at the top of the file documents the alignment and quantification steps
that generate the input `.tsv` tables (not run by this script):

- **polyA / rRNA:** `minimap2` splice-aware alignment (`-ax splice-ont`) to `hg38.fa`,
  sorted and filtered with `samtools`.
- **tRNA:** a two-pass strategy — first align to an rRNA+tRNA reference with tuned
  short-read `minimap2` parameters, then re-align the unaligned remainder to the genome
  and merge.
- **Quantification:** `featureCounts` in long-read mode (`-L`) against a combined GTF.

## Inputs

- Per-sample count tables (`*.merged.aligned.sorted.tsv` / `*_final_merge.tsv`) from
  `featureCounts`, one column of counts per sample plus a `Geneid` column.
- A combined/merged count table used as the `combined` column.
- A GENCODE GTF (e.g. `gencode.v49.basic.annotation.gtf`) for biotype annotation,
  loaded via `dmode.metagene_plot.prepare_gene_body_coverage`.

**Note:** all input/output paths are hard-coded (e.g. `/global/cfs/cdirs/m5243/...`) and
must be edited to match your environment. Some paths within the file are inconsistent
(mixing cluster, Synology, and local mounts) and will need to be reconciled before running.

## Outputs

Stacked bar charts and source tables written next to the inputs, including:

- `gene_type_composition_*_combined_only.tsv` — biotype composition tables (%).
- `Figure1_gene_type_composition*.{pdf,png,svg}` — polyA biotype composition.
- `Figure1_rRNA_biotype_CPMS*.pdf` — rRNA composition.
- `Figure1_tRNA_biotype_CPMS.pdf` and
  `Figure1_tRNA_anticodon_aminoacid_distribution_CPMS.pdf` — tRNA composition and
  per-isodecoder breakdown.

## Key details

- **Normalization:** counts are scaled to CPM (`counts / (column_sum / 1e6)`). The
  reads-per-kilobase (length) normalization is present but commented out.
- **Biotype filtering:** biotypes exceeding 1% relative share in at least one sample are
  retained; everything else is summed into an `other` category and bars are ordered by
  total signal.
- **Color maps:** two curated palettes are defined —
  `gene_type_colors` (an earth-tone palette keyed by GENCODE biotype) and
  `transcript_colors` (blues/greens/rusts keyed by individual rRNA and tRNA isodecoders).
- **tRNA relabeling:** raw reference names (e.g. `hs_tRNAAla_CGC`, `hs_mttRNAAla_TGC`) are
  mapped to compact labels (`Ala_CGC`, `mt-Ala_TGC`) for plotting.

## Dependencies

- `pandas`, `polars`, `matplotlib`
- `dmode` — an internal/project module providing
  `dmode.metagene_plot.prepare_gene_body_coverage` for GTF parsing and gene-body coverage.

## Usage

The script is meant to be run cell-by-cell in an interactive session rather than as a
single batch job. Before running:

1. Install the dependencies and ensure the `dmode` package is importable. 
A gtf parser is imported from dmode, which transforms the gtf into a table format and 
extracts 3'UTR, CDS and 5'UTR for the MANE transcript of a gene in a second table. 
2. Update the `base` path and every input filename to point at your data.
3. Confirm the GTF path(s) and reference names match your annotation.
4. Execute the cells for the RNA class(es) of interest (polyA, rRNA, tRNA).

## Caveats

- Paths are hard-coded and environment-specific.
- The three blocks reuse variable names (`base`, `merged`, `cpm`, `sample_cols`,
  `collapse`, etc.), so cells must be run in order within each block.
- The `rel_plot` / `abs_plot` variables in the final tRNA isodecoder plot are carried over
  from the preceding cell; double-check they reflect the intended grouping before use.
