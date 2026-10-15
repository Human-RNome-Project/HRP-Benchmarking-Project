#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
biotype_analysis.py
===================

Compute and visualise the RNA **biotype composition** (Ensembl/GENCODE
``gene_type``) of long-read (e.g. Oxford Nanopore) sequencing samples,
starting from per-sample ``featureCounts`` read-count tables.

The original version of this script was a single, hard-coded notebook that
analysed three specific datasets (polyA RNA, rRNA and tRNA) with absolute
file paths baked in.  This rewrite keeps *exactly the same analysis logic*
but turns it into:

* a small library of reusable, documented functions, and
* a command-line tool (``argparse`` + ``if __name__ == "__main__"``) that
  can be pointed at **any** set of featureCounts tables and a GTF
  annotation -- either through command-line flags for a single run, or
  through a YAML config file describing several analyses at once.

--------------------------------------------------------------------------
Analysis pipeline (unchanged from the original)
--------------------------------------------------------------------------
1. Load one featureCounts table per sample. featureCounts writes a
   commented header (``#`` lines) followed by a table whose **last column**
   holds the read counts; the gene identifier is in the ``Geneid`` column.
2. Merge all samples into one gene x sample count matrix (outer join on the
   gene id; genes absent from a sample become 0).
3. Normalise to **CPM** (counts per million): divide each sample column by
   its own total counts / 1e6.
4. Attach the biotype of every gene from the GTF annotation
   (``gene_type``, ``gene_name``, ``gene_id`` for ``feature == "gene"``).
5. Optionally relabel biotype names (e.g. collapse ``rRNA_5S`` -> ``rRNA``).
6. Aggregate CPM per biotype, write a composition table (percent of total
   signal per sample) to TSV.
7. Draw stacked horizontal bar charts of the relative biotype composition,
   keeping the most abundant biotypes and collapsing the rest into
   ``"other"``.

--------------------------------------------------------------------------
GTF parsing
--------------------------------------------------------------------------
The original relied on the internal ``dmode`` package
(``dmode.metagene_plot.prepare_gene_body_coverage`` -> ``dmode.utility.gtf_to_df``).
Only the gene-level ``gene_type`` / ``gene_name`` / ``gene_id`` columns are
needed here.  To keep the script runnable on machines that do not have
``dmode`` installed, ``load_gene_annotation`` uses ``dmode`` when it is
importable and otherwise falls back to a self-contained GTF parser that
extracts the same columns.  Use ``--no-dmode`` to force the built-in parser.

--------------------------------------------------------------------------
Examples
--------------------------------------------------------------------------
Single run from the command line::

    python biotype_analysis.py \\
        --name polyA \\
        --counts sample1.tsv sample2.tsv sample3.tsv \\
        --names  HRP_A_001_1 HRP_A_002_1 HRP_A_003_1 \\
        --combined polyA_native_merged.tsv \\
        --gtf gencode.v49.basic.annotation.gtf \\
        --output-dir results/ \\
        --palette gene_type

Reproducible batch run from a config file (see ``config.example.yaml``)::

    python biotype_analysis.py --config config.example.yaml

The module is also importable, e.g.::

    from biotype_analysis import (
        load_count_matrix, compute_cpm, load_gene_annotation,
        annotate_counts, composition_table, plot_composition, run_analysis,
    )
"""

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

logger = logging.getLogger("biotype_analysis")

# A path may be given as a string or a pathlib.Path throughout.
PathLike = Union[str, Path]

# =========================================================================
# Default colour palettes
# -------------------------------------------------------------------------
# These are the exact palettes from the original script.  They map a biotype
# label to a hex colour.  Any label not present in the palette is drawn in
# ``DEFAULT_OTHER_COLOR``.  Palettes can be overridden from a config file or
# ``--palette`` (see ``resolve_palette``).
# =========================================================================

DEFAULT_OTHER_COLOR = "#CCCCCC"

# Gene-level biotypes (protein_coding, rRNA, tRNA families, IG/TR genes,
# pseudogenes, ...).  Used for the polyA and rRNA figures in the original.
DEFAULT_GENE_TYPE_COLORS: Dict[str, str] = {
    # Protein-coding -- deep forest greens
    "protein_coding": "#2E5E3A",
    "TEC":            "#6B8E5A",
    # tRNA family -- warm earth/bark tones
    "tRNA":           "#A0522D",
    "Mt_tRNA":        "#8B4513",
    # rRNA family -- clay / terracotta
    "rRNA":           "#C1693C",
    "Mt_rRNA":        "#9C4A2A",
    "rRNA_pseudogene": "#D4926A",
    "18S":            "#E08D5B",
    "28S":            "#C76B3A",
    "5.8S":           "#E8A87C",
    "5S":             "#B5651D",
    # Small/regulatory RNAs -- wildflower & moss tones
    "miRNA":          "#7B68A6",   # lavender
    "snoRNA":         "#9370B5",   # mauve-violet
    "snRNA":          "#5D8AA8",   # cornflower
    "scaRNA":         "#6CA0B5",   # dusty teal
    "sRNA":           "#88B04B",   # moss green
    "misc_RNA":       "#A3A847",   # olive
    "vault_RNA":      "#4F7942",   # fern
    "ribozyme":       "#3D6B4A",   # pine
    "lncRNA":         "#4A7C8C",   # slate teal
    # IG genes -- autumn reds/oranges
    "IG_C_gene":      "#B22222",
    "IG_V_gene":      "#C1432E",
    "IG_D_gene":      "#D35438",
    "IG_J_gene":      "#A83232",
    # IG pseudogenes -- faded autumn
    "IG_C_pseudogene": "#D98C7A",
    "IG_V_pseudogene": "#E0A08F",
    "IG_J_pseudogene": "#CD9384",
    "IG_pseudogene":  "#C99182",
    # TR genes -- golden / amber
    "TR_C_gene":      "#DAA520",
    "TR_V_gene":      "#E0B23A",
    "TR_D_gene":      "#C99700",
    "TR_J_gene":      "#EAC54F",
    # TR pseudogenes -- muted gold
    "TR_V_pseudogene": "#D4BC72",
    "TR_J_pseudogene": "#CBB76A",
    # General pseudogenes -- stone / lichen greys
    "processed_pseudogene":               "#9A9B7A",
    "unprocessed_pseudogene":             "#8A8B6C",
    "transcribed_processed_pseudogene":   "#A8A98A",
    "transcribed_unprocessed_pseudogene": "#7C7D60",
    "transcribed_unitary_pseudogene":     "#B0B190",
    "unitary_pseudogene":                 "#90916F",
    "translated_processed_pseudogene":    "#A0A17E",
    # Other
    "artifact":       "#6E6E6E",   # neutral grey
}

# Per-transcript rRNA / tRNA isoacceptor palette.  Used for the tRNA
# anticodon/amino-acid figure in the original.
DEFAULT_TRANSCRIPT_COLORS: Dict[str, str] = {
    # rRNA -- ocean blues
    "5S": "#1F4E79", "5.8S": "#2E6CA4", "18S": "#3D8BC0", "28S": "#5BA3D0",
    "mt-12S": "#4A6670", "mt-16S": "#637E88",
    # cytoplasmic tRNA -- green -> earth amino-acid spectrum
    "Ala_CGC": "#2E5E3A", "Ala_TGC": "#3A6E45", "Ala_AGC": "#477E50",
    "Arg_ACG": "#548D5B", "Arg_TCG": "#629C66", "Arg_CCG2": "#5E9B5E",
    "Arg_TCT2": "#6FA86B", "Arg_CCT": "#7DB073", "Arg_CCG1": "#8AB87B",
    "Arg_TCT1": "#97C084", "Asn_GTT": "#6B9362", "Asp_GTC": "#759C5A",
    "Cys_GCA": "#83A551", "Gln_CTG_TTG": "#90AD49", "Glu_CTC": "#9DB541",
    "Glu_TTC": "#A8BC45", "Gly_GCC": "#B0C04A", "Gly_TCC": "#A3A847",
    "Gly_CCC": "#959E42", "His_GTG": "#AEB158", "Ile_AAT": "#B8B05A",
    "Ile_TAT": "#C2AF55", "Ile_GAT": "#CBAE4E", "Leu_AAG": "#D4AD47",
    "Leu_CAG": "#DAA520", "Leu_CAA": "#D6A032", "Leu_TAA": "#CC9A3A",
    "Lys_CTT": "#C29440", "Lys_TTT": "#B98C3E", "Met_CAT": "#B5651D",
    "Phe_GAA": "#AD7A35", "Pro_AGG_CGG_TGG": "#A86E2E", "SeC_TCA": "#9C6328",
    "Ser_AGA": "#A0522D", "Ser_CGA": "#985534", "Ser_GCT": "#8F5A3C",
    "Thr_AGT": "#8B4F2E", "Thr_CGT1": "#7E5536", "Thr_TGT": "#76502F",
    "Thr_CGT2": "#6E4A2B", "Trp_CCA": "#7A5A3A", "Tyr_GTA2": "#6B4F38",
    "Tyr_ATA": "#5F4A33", "Tyr_GTA1": "#54442F", "Val_AAC_CAC": "#4E5340",
    "Val_TAC": "#5A6049",
    # mitochondrial tRNA -- terracotta / rust
    "mt-Ala_TGC": "#C1693C", "mt-Arg_TCG": "#CD7548", "mt-Asn_GTT": "#D98153",
    "mt-Asp_GTC": "#B85A30", "mt-Cys_GCA": "#C26538", "mt-Gln_TTG": "#CE7040",
    "mt-Glu_TTC": "#A84E2A", "mt-Gly_TCC": "#B45932", "mt-His_GTG": "#C0643A",
    "mt-Ile_GAT": "#9C4A2A", "mt-Leu_TAA": "#A85432", "mt-Leu_TAG": "#B45F3A",
    "mt-Lys_TTT": "#92442A", "mt-Met_CAT": "#9E4E32", "mt-Phe_GAA": "#AA583A",
    "mt-Pro_TGG": "#883E28", "mt-Ser_GCT": "#944830", "mt-Ser_TGA": "#A05238",
    "mt-Thr_TGT": "#7E3826", "mt-Trp_TCA": "#8A422E", "mt-Tyr_GTA": "#964C36",
    "mt-Val_TAC": "#723424",
}

# Named palettes usable from the config / CLI (``palette: gene_type``).
BUILTIN_PALETTES: Dict[str, Dict[str, str]] = {
    "gene_type": DEFAULT_GENE_TYPE_COLORS,
    "transcript": DEFAULT_TRANSCRIPT_COLORS,
}

# Matplotlib rcParams presets matching the two font sizes used in the
# original script ("large" for the polyA/rRNA figures, "small" for tRNA).
FONT_PRESETS: Dict[str, Dict[str, float]] = {
    "large": {
        "font.size": 20, "axes.titlesize": 24, "axes.labelsize": 28,
        "xtick.labelsize": 20, "ytick.labelsize": 20,
        "legend.fontsize": 20, "legend.title_fontsize": 22,
    },
    "small": {
        "font.size": 11, "axes.titlesize": 14, "axes.labelsize": 13,
        "xtick.labelsize": 11, "ytick.labelsize": 11,
        "legend.fontsize": 11, "legend.title_fontsize": 12,
    },
}


# =========================================================================
# Loading and normalising count tables
# =========================================================================

def load_count_table(
    path: PathLike,
    sample_name: Optional[str] = None,
    id_column: str = "Geneid",
) -> pd.Series:
    """Load a single ``featureCounts`` table and return its count column.

    Parameters
    ----------
    path :
        Path to a featureCounts ``.tsv`` (tab separated, ``#`` comment lines,
        one header row).
    sample_name :
        Name to give the returned count column.  If ``None`` the original
        last-column header (usually the input BAM path) is kept.
    id_column :
        Column holding the gene identifier used as the index. Default
        ``"Geneid"`` (the featureCounts default).

    Returns
    -------
    pandas.Series
        Read counts indexed by gene id, named ``sample_name``.
    """
    df = pd.read_csv(path, sep="\t", header=0, comment="#")
    if id_column not in df.columns:
        raise KeyError(
            f"{path}: expected an id column {id_column!r}; found {list(df.columns)}"
        )
    # featureCounts always writes the read counts in the LAST column.
    counts_col = df.columns[-1]
    series = df.set_index(id_column)[counts_col]
    series.name = sample_name if sample_name is not None else counts_col
    return series


def load_count_matrix(
    counts: Union[Mapping[str, PathLike], Sequence[PathLike]],
    id_column: str = "Geneid",
) -> pd.DataFrame:
    """Load several featureCounts tables into one gene x sample matrix.

    Parameters
    ----------
    counts :
        Either a mapping ``{sample_name: path}`` (preserves the given order
        and names) or a plain sequence of paths (sample names are derived
        from the file names).
    id_column :
        Gene-id column, passed through to :func:`load_count_table`.

    Returns
    -------
    pandas.DataFrame
        Counts indexed by gene id, one column per sample. Genes missing in a
        sample (outer join) are filled with 0.
    """
    named = _as_named_mapping(counts)
    series = [
        load_count_table(path, sample_name=name, id_column=id_column)
        for name, path in named.items()
    ]
    if not series:
        raise ValueError("No count tables were provided.")
    # Outer join keeps every gene seen in any sample; absent -> 0.
    matrix = pd.concat(series, axis=1).fillna(0)
    logger.info(
        "Loaded %d samples x %d genes.", matrix.shape[1], matrix.shape[0]
    )
    return matrix


def compute_cpm(matrix: pd.DataFrame) -> pd.DataFrame:
    """Convert a raw count matrix to CPM (counts per million).

    Each sample (column) is divided by its own total counts divided by 1e6,
    so every column sums to 1,000,000. This is library-size normalisation and
    makes samples of different sequencing depth directly comparable.
    """
    scaling = matrix.sum(axis=0) / 1e6
    if (scaling == 0).any():
        empty = list(scaling.index[scaling == 0])
        raise ValueError(f"These samples have zero total counts: {empty}")
    return matrix.div(scaling, axis=1)


# =========================================================================
# Gene annotation (biotype) from a GTF
# =========================================================================

# Matches GTF attribute entries of the form:  key "value";
_GTF_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _extract_gtf_attr(attr_series: pd.Series, key: str) -> pd.Series:
    """Vectorised extraction of a single GTF attribute value."""
    return attr_series.str.extract(rf'{key}\s+"([^"]+)"', expand=False)


def _load_gene_annotation_builtin(gtf_path: PathLike) -> pd.DataFrame:
    """Standalone GTF parser (no external deps) for gene-level biotypes.

    Reads the 9-column GTF, keeps ``feature == "gene"`` rows and pulls
    ``gene_id``, ``gene_type`` (falling back to ``gene_biotype`` as used by
    Ensembl) and ``gene_name`` out of the attribute field.
    """
    cols = [
        "seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attribute",
    ]
    df = pd.read_csv(
        gtf_path, sep="\t", comment="#", header=None, names=cols,
        dtype=str, quoting=3,  # quoting=3 == csv.QUOTE_NONE: keep the quotes
    )
    genes = df[df["feature"] == "gene"].copy()
    if genes.empty:
        raise ValueError(
            f"{gtf_path}: no rows with feature == 'gene' were found."
        )
    attr = genes["attribute"]
    gene_type = _extract_gtf_attr(attr, "gene_type")
    # Ensembl GTFs use gene_biotype instead of GENCODE's gene_type.
    gene_type = gene_type.fillna(_extract_gtf_attr(attr, "gene_biotype"))
    out = pd.DataFrame({
        "gene_type": gene_type,
        "gene_name": _extract_gtf_attr(attr, "gene_name"),
        "gene_id": _extract_gtf_attr(attr, "gene_id"),
    })
    out = out.dropna(subset=["gene_id"]).reset_index(drop=True)
    logger.info("Parsed %d genes from GTF (built-in parser).", len(out))
    return out


def _load_gene_annotation_dmode(gtf_path: PathLike) -> pd.DataFrame:
    """Gene-level biotypes via the internal ``dmode`` package.

    Mirrors what the original script did through
    ``dmode.metagene_plot.prepare_gene_body_coverage`` but only touches the
    part that this analysis needs: ``dmode.utility.gtf_to_df`` followed by a
    ``feature == "gene"`` selection.
    """
    import polars as pl  # imported lazily; only needed on the dmode path
    from dmode.utility import gtf_to_df

    gtf_df = gtf_to_df(str(gtf_path))
    if not isinstance(gtf_df, pl.DataFrame):
        gtf_df = pl.from_pandas(gtf_df)
    genes = (
        gtf_df.filter(pl.col("feature") == "gene")
              .select(["gene_type", "gene_name", "gene_id"])
              .to_pandas()
    )
    logger.info("Parsed %d genes from GTF (dmode).", len(genes))
    return genes


def load_gene_annotation(
    gtf_path: PathLike,
    prefer_dmode: bool = True,
) -> pd.DataFrame:
    """Return a gene-level annotation table from a GTF.

    Tries the internal ``dmode`` package first (when ``prefer_dmode`` is
    True and it is importable), and otherwise -- or on any failure -- falls
    back to the built-in parser, so the script runs anywhere.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_type``, ``gene_name``, ``gene_id`` (one row per gene).
    """
    if prefer_dmode:
        try:
            return _load_gene_annotation_dmode(gtf_path)
        except Exception as exc:  # ImportError, or any dmode/polars issue
            logger.warning(
                "dmode GTF parsing unavailable (%s); using built-in parser.",
                exc,
            )
    return _load_gene_annotation_builtin(gtf_path)


def annotate_counts(
    cpm: pd.DataFrame,
    annotation: pd.DataFrame,
    relabel: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """Attach the biotype of each gene to the CPM matrix.

    Parameters
    ----------
    cpm :
        CPM matrix indexed by gene id (output of :func:`compute_cpm`).
    annotation :
        Gene annotation table with a ``gene_id`` column
        (output of :func:`load_gene_annotation`).
    relabel :
        Optional ``{old_gene_type: new_gene_type}`` map, applied after the
        merge (e.g. collapse ``rRNA_5S`` -> ``rRNA``).

    Returns
    -------
    pandas.DataFrame
        The CPM sample columns plus ``gene_type``, ``gene_name``, ``gene_id``.
        Inner join: genes absent from the annotation are dropped.
    """
    merged = cpm.merge(
        annotation, left_index=True, right_on="gene_id", how="inner"
    ).reset_index(drop=True)
    if relabel:
        merged["gene_type"] = merged["gene_type"].replace(dict(relabel))
    logger.info(
        "Annotated %d genes across %d biotypes.",
        len(merged), merged["gene_type"].nunique(),
    )
    return merged


# =========================================================================
# Composition table and plotting
# =========================================================================

def composition_table(
    annotated: pd.DataFrame,
    sample_cols: Sequence[str],
    sort_by: Optional[str] = None,
) -> pd.DataFrame:
    """Aggregate CPM per biotype and express it as a percentage.

    Each biotype's summed CPM is divided by the sample total and multiplied
    by 100, so every column sums to 100 %.  (For CPM input this is identical
    to the original ``sum / 1e6 * 100``.)

    Parameters
    ----------
    annotated :
        Output of :func:`annotate_counts`.
    sample_cols :
        Columns (samples) to include.
    sort_by :
        Sample column to sort biotypes by (descending). If ``None`` or not
        present, biotypes are ordered by total signal across all columns.
    """
    comp = annotated.groupby("gene_type")[list(sample_cols)].sum()
    pct = comp.div(comp.sum(axis=0), axis=1) * 100.0
    if sort_by and sort_by in pct.columns:
        pct = pct.sort_values(by=sort_by, ascending=False)
    else:
        order = pct.sum(axis=1).sort_values(ascending=False).index
        pct = pct.loc[order]
    return pct


def select_biotypes(
    rel: pd.DataFrame,
    method: str = "fraction",
    threshold: float = 0.01,
    top_n: int = 8,
) -> pd.Index:
    """Choose which biotypes to show explicitly (rest -> "other").

    ``method="fraction"`` keeps biotypes whose relative share exceeds
    ``threshold`` in at least one sample; ``method="top_n"`` keeps the
    ``top_n`` biotypes with the largest total signal.  Either way, the kept
    biotypes are returned ordered by total signal (descending) so the
    "other" segment stacks consistently at the end.
    """
    if method == "top_n":
        return rel.sum(axis=1).sort_values(ascending=False).index[:top_n]
    if method == "fraction":
        keep = rel.index[(rel > threshold).any(axis=1)]
        return rel.loc[keep].sum(axis=1).sort_values(ascending=False).index
    raise ValueError(f"Unknown selection method: {method!r}")


def collapse_other(
    df: pd.DataFrame,
    keep: pd.Index,
    other_label: str = "other",
) -> pd.DataFrame:
    """Keep ``keep`` rows and sum everything else into a single row."""
    out = df.loc[keep].copy()
    dropped = df.drop(index=keep)
    if not dropped.empty:
        out.loc[other_label] = dropped.sum(axis=0)
    return out


def plot_composition(
    annotated: pd.DataFrame,
    sample_cols: Sequence[str],
    colors: Optional[Mapping[str, str]] = None,
    *,
    method: str = "fraction",
    threshold: float = 0.01,
    top_n: int = 8,
    y_label: str = "Sample",
    x_label: str = "Relative composition (CPM / 1M)",
    hide_yticks: bool = False,
    figsize: Tuple[float, float] = (16, 8),
    other_color: str = DEFAULT_OTHER_COLOR,
    other_label: str = "other",
    ax: "Optional[matplotlib.axes.Axes]" = None,
) -> Tuple["matplotlib.figure.Figure", "matplotlib.axes.Axes", pd.DataFrame]:
    """Draw a stacked horizontal bar chart of relative biotype composition.

    One horizontal bar per sample; segment widths are the biotype's relative
    share (0-1). Returns ``(figure, axes, plotted_dataframe)``.
    """
    comp = annotated.groupby("gene_type")[list(sample_cols)].sum()
    rel = comp.div(comp.sum(axis=0), axis=1)
    keep = select_biotypes(rel, method=method, threshold=threshold, top_n=top_n)
    rel_plot = collapse_other(rel, keep, other_label=other_label)

    palette = dict(colors or {})
    bar_colors = [palette.get(b, other_color) for b in rel_plot.index]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    rel_plot.T.plot(kind="barh", stacked=True, color=bar_colors, width=0.8, ax=ax)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(0, 1)
    if hide_yticks:
        ax.set_yticklabels("")
    ax.legend(title="Biotype", bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    return fig, ax, rel_plot


# =========================================================================
# High-level driver: one analysis
# =========================================================================

def run_analysis(
    name: str,
    counts: Union[Mapping[str, PathLike], Sequence[PathLike]],
    gtf: PathLike,
    output_dir: PathLike,
    *,
    combined: Union[Mapping[str, PathLike], PathLike, None] = None,
    combined_name: str = "combined",
    relabel: Optional[Mapping[str, str]] = None,
    colors: Optional[Mapping[str, str]] = None,
    plots: Optional[List[Dict[str, Any]]] = None,
    id_column: str = "Geneid",
    sort_by: Optional[str] = "combined",
    prefer_dmode: bool = True,
    annotation: Optional[pd.DataFrame] = None,
    show: bool = False,
    table_name: Optional[str] = None,
) -> Dict[str, Any]:
    """Run the full pipeline for one dataset and write its outputs.

    Parameters
    ----------
    name :
        Short label for this analysis (used in default output filenames).
    counts :
        featureCounts tables: ``{sample_name: path}`` or a list of paths.
    gtf :
        GTF annotation file.
    output_dir :
        Directory for the composition TSV and figures (created if needed).
    combined :
        Optional extra "merged" table, added as one more sample column and
        used as the default ``sort_by`` column.
    combined_name :
        Column name for ``combined`` (default ``"combined"``).
    relabel :
        Optional biotype relabelling map (see :func:`annotate_counts`).
    colors :
        Biotype -> hex colour palette for the figures.
    plots :
        List of plot-spec dicts. Each may set: ``sample_cols`` (list, or
        ``None``/empty = all non-combined samples), ``y_label``,
        ``hide_yticks``, ``figsize``, ``method``, ``threshold``, ``top_n``,
        ``reverse`` (flip bar order) and ``outfile``. If ``None``, one plot
        with all samples is drawn.
    annotation :
        Pre-loaded gene annotation (reused across analyses sharing a GTF).
    show :
        Call ``plt.show()`` for each figure (interactive use).

    Returns
    -------
    dict
        ``{"table": Path, "figures": [Path, ...], "annotated": DataFrame}``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Assemble all count sources (samples + optional combined column).
    sources: Dict[str, PathLike] = dict(_as_named_mapping(counts))
    if combined is not None:
        extra = _as_named_mapping(combined, default_name=combined_name)
        sources.update(extra)

    matrix = load_count_matrix(sources, id_column=id_column)
    cpm = compute_cpm(matrix)

    if annotation is None:
        annotation = load_gene_annotation(gtf, prefer_dmode=prefer_dmode)
    annotated = annotate_counts(cpm, annotation, relabel=relabel)

    all_cols = list(matrix.columns)
    sort_col = sort_by if (sort_by and sort_by in all_cols) else None

    # --- composition table (over all samples) ---
    table = composition_table(annotated, all_cols, sort_by=sort_col)
    table_path = output_dir / (table_name or f"gene_type_composition_{name}.tsv")
    with pd.option_context("display.float_format", "{:.6f}".format):
        table.to_csv(table_path, sep="\t", index=True, header=True)
    logger.info("Wrote composition table: %s", table_path)

    # --- figures ---
    if plots is None:
        plots = [{"sample_cols": None, "y_label": name}]

    figures: List[Path] = []
    for i, spec in enumerate(plots):
        sample_cols = spec.get("sample_cols")
        if not sample_cols:  # None or [] -> all non-combined samples
            sample_cols = [c for c in all_cols if c != combined_name] or all_cols
        if spec.get("reverse"):
            sample_cols = list(sample_cols)[::-1]

        fig, _ax, _plotted = plot_composition(
            annotated,
            sample_cols,
            colors=colors,
            method=spec.get("method", "fraction"),
            threshold=spec.get("threshold", 0.01),
            top_n=spec.get("top_n", 8),
            y_label=spec.get("y_label", name),
            x_label=spec.get("x_label", "Relative composition (CPM / 1M)"),
            hide_yticks=spec.get("hide_yticks", False),
            figsize=tuple(spec.get("figsize", (16, 8))),
        )
        outfile = spec.get("outfile", f"{name}_composition_{i + 1}.pdf")
        outpath = output_dir / outfile
        fig.savefig(outpath, bbox_inches="tight")
        logger.info("Wrote figure: %s", outpath)
        if show:
            plt.show()
        else:
            plt.close(fig)
        figures.append(outpath)

    return {"table": table_path, "figures": figures, "annotated": annotated}


# =========================================================================
# Config-file driven runs
# =========================================================================

def _as_named_mapping(
    obj: Union[Mapping[str, PathLike], Sequence[PathLike], PathLike],
    default_name: Optional[str] = None,
) -> "Dict[str, PathLike]":
    """Normalise a count/combined spec into an ordered ``{name: path}`` dict.

    Accepts a mapping (returned as-is), a single path, or a sequence of
    paths (names derived from file stems). If ``default_name`` is given and
    exactly one path results, that name is used.
    """
    if isinstance(obj, Mapping):
        return dict(obj)
    if isinstance(obj, (str, Path)):
        items: List[PathLike] = [obj]
    else:
        items = list(obj)
    if default_name and len(items) == 1:
        return {default_name: items[0]}
    out: "Dict[str, PathLike]" = {}
    for path in items:
        out[Path(path).stem] = path
    return out


def resolve_palette(
    value: Union[str, Mapping[str, str], None],
) -> Dict[str, str]:
    """Resolve a palette spec into a ``{biotype: colour}`` dict.

    ``value`` may be: ``None`` (empty palette), a mapping (used directly),
    the name of a built-in palette (``"gene_type"`` / ``"transcript"``), or
    a path to a YAML/JSON file mapping biotypes to colours.
    """
    if value is None:
        return {}
    if isinstance(value, Mapping):
        return dict(value)
    if isinstance(value, str):
        if value in BUILTIN_PALETTES:
            return dict(BUILTIN_PALETTES[value])
        loaded = _load_mapping_file(value)
        if loaded is not None:
            return loaded
        logger.warning("Palette %r is not a known name or file; using empty.", value)
    return {}


def resolve_mapping(
    value: Union[str, Mapping[str, str], None],
) -> Dict[str, str]:
    """Resolve a relabel spec (inline mapping or YAML/JSON file path)."""
    if value is None:
        return {}
    if isinstance(value, Mapping):
        return dict(value)
    loaded = _load_mapping_file(str(value))
    return loaded or {}


def _load_mapping_file(path: str) -> Optional[Dict[str, str]]:
    """Load a ``{str: str}`` mapping from a YAML or JSON file, if it exists."""
    p = Path(path)
    if not p.exists():
        return None
    text = p.read_text()
    if p.suffix.lower() in (".yaml", ".yml"):
        import yaml
        return dict(yaml.safe_load(text) or {})
    if p.suffix.lower() == ".json":
        return dict(json.loads(text))
    # Unknown extension: try YAML (a JSON file is valid YAML too).
    import yaml
    return dict(yaml.safe_load(text) or {})


def load_config(path: PathLike) -> Dict[str, Any]:
    """Load a YAML analysis config file into a dict."""
    import yaml
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict) or "analyses" not in cfg:
        raise ValueError(
            "Config must be a mapping with an 'analyses' list "
            "(see config.example.yaml)."
        )
    return cfg


def run_from_config(config: Mapping[str, Any], show: bool = False) -> List[Dict[str, Any]]:
    """Run every analysis described in a config mapping.

    Top-level keys: ``defaults`` (values inherited by each analysis) and
    ``analyses`` (a list of analysis specs). Each analysis spec accepts the
    same keys as :func:`run_analysis` plus ``palette``/``relabel`` (name,
    inline mapping, or file path). GTF annotations are cached and reused
    across analyses that share the same GTF path.
    """
    defaults: Dict[str, Any] = dict(config.get("defaults", {}))
    analyses = config.get("analyses", [])
    if not analyses:
        raise ValueError("Config contains no analyses.")

    annotation_cache: Dict[str, pd.DataFrame] = {}
    results: List[Dict[str, Any]] = []

    for spec in analyses:
        merged = {**defaults, **spec}  # per-analysis keys win over defaults
        gtf = merged["gtf"]
        prefer_dmode = merged.get("prefer_dmode", True)

        cache_key = f"{gtf}|{prefer_dmode}"
        if cache_key not in annotation_cache:
            annotation_cache[cache_key] = load_gene_annotation(
                gtf, prefer_dmode=prefer_dmode
            )
        annotation = annotation_cache[cache_key]

        result = run_analysis(
            name=merged.get("name", "analysis"),
            counts=merged["counts"],
            gtf=gtf,
            output_dir=merged.get("output_dir", "."),
            combined=merged.get("combined"),
            combined_name=merged.get("combined_name", "combined"),
            relabel=resolve_mapping(merged.get("relabel")),
            colors=resolve_palette(merged.get("palette", "gene_type")),
            plots=merged.get("plots"),
            id_column=merged.get("id_column", "Geneid"),
            sort_by=merged.get("sort_by", "combined"),
            prefer_dmode=prefer_dmode,
            annotation=annotation,
            show=show,
            table_name=merged.get("table_name"),
        )
        results.append(result)

    return results


# =========================================================================
# Command-line interface
# =========================================================================

def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the command-line interface."""
    parser = argparse.ArgumentParser(
        prog="biotype_analysis.py",
        description=(
            "Compute and plot the RNA biotype (gene_type) composition of "
            "featureCounts tables. Provide either --config for a batch run, "
            "or --counts + --gtf for a single run."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="See config.example.yaml for the batch-config format.",
    )

    # --- batch mode ---
    parser.add_argument(
        "--config",
        help="YAML config describing one or more analyses (batch mode).",
    )

    # --- single-run inputs ---
    grp_in = parser.add_argument_group("single-run inputs")
    grp_in.add_argument(
        "--counts", nargs="+",
        help="featureCounts .tsv tables, one per sample.",
    )
    grp_in.add_argument(
        "--names", nargs="+",
        help="Sample names, in the same order as --counts (optional).",
    )
    grp_in.add_argument(
        "--combined",
        help="Optional merged/combined table (added as an extra column and "
             "used as the default --sort-by column).",
    )
    grp_in.add_argument(
        "--combined-name", default="combined",
        help="Column name for the --combined table.",
    )
    grp_in.add_argument("--gtf", help="GTF annotation file.")
    grp_in.add_argument(
        "--id-column", default="Geneid",
        help="Gene-id column in the count tables.",
    )

    # --- outputs ---
    grp_out = parser.add_argument_group("outputs")
    grp_out.add_argument(
        "--output-dir", "-o", default="biotype_results",
        help="Directory for the composition table and figures.",
    )
    grp_out.add_argument(
        "--name", default="analysis",
        help="Short label used in output file names.",
    )

    # --- annotation / relabelling / colours ---
    grp_ann = parser.add_argument_group("annotation and appearance")
    grp_ann.add_argument(
        "--relabel",
        help="YAML/JSON file with a {old_biotype: new_biotype} map.",
    )
    grp_ann.add_argument(
        "--palette", default="gene_type",
        help="Built-in palette name ('gene_type' or 'transcript'), or a "
             "YAML/JSON file mapping biotypes to hex colours.",
    )
    grp_ann.add_argument(
        "--no-dmode", action="store_true",
        help="Force the built-in GTF parser even if dmode is installed.",
    )

    # --- biotype selection ---
    grp_sel = parser.add_argument_group("biotype selection")
    grp_sel.add_argument(
        "--select", choices=["fraction", "top_n"], default="fraction",
        help="Keep biotypes above --threshold in any sample, or the "
             "--top-n most abundant.",
    )
    grp_sel.add_argument(
        "--threshold", type=float, default=0.01,
        help="Relative-share threshold for --select fraction.",
    )
    grp_sel.add_argument(
        "--top-n", type=int, default=8,
        help="Number of biotypes to keep for --select top_n.",
    )
    grp_sel.add_argument(
        "--sort-by", default="combined",
        help="Sample column to sort the composition table by.",
    )

    # --- plot appearance ---
    grp_plot = parser.add_argument_group("plot appearance")
    grp_plot.add_argument("--y-label", help="Y-axis label (default: --name).")
    grp_plot.add_argument(
        "--x-label", default="Relative composition (CPM / 1M)",
        help="X-axis label.",
    )
    grp_plot.add_argument(
        "--hide-yticks", action="store_true",
        help="Hide per-sample y tick labels.",
    )
    grp_plot.add_argument(
        "--figsize", nargs=2, type=float, default=[16, 8],
        metavar=("W", "H"), help="Figure size in inches.",
    )
    grp_plot.add_argument(
        "--reverse-samples", action="store_true",
        help="Reverse the order samples are stacked in the plot.",
    )
    grp_plot.add_argument(
        "--combined-only", action="store_true",
        help="Plot only the --combined column instead of all samples.",
    )
    grp_plot.add_argument(
        "--font-preset", choices=["small", "large"], default="large",
        help="Matplotlib font-size preset.",
    )

    # --- misc ---
    parser.add_argument(
        "--show", action="store_true",
        help="Display figures interactively (default: save only).",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging.",
    )
    return parser


def _to_named_counts(
    counts: Sequence[str], names: Optional[Sequence[str]]
) -> Union[Dict[str, str], List[str]]:
    """Turn --counts/--names into a mapping (or a bare list if no names)."""
    if names:
        if len(names) != len(counts):
            raise SystemExit(
                f"--names ({len(names)}) must match --counts ({len(counts)})."
            )
        return dict(zip(names, counts))
    return list(counts)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    # Use a non-interactive backend unless the user asked to view figures,
    # so the script works over SSH / on a cluster with no display.
    if not args.show:
        matplotlib.use("Agg", force=True)

    # --- batch mode -------------------------------------------------------
    if args.config:
        plt.rcParams.update(FONT_PRESETS[args.font_preset])
        config = load_config(args.config)
        run_from_config(config, show=args.show)
        return 0

    # --- single-run mode --------------------------------------------------
    if not args.counts or not args.gtf:
        parser.error("--counts and --gtf are required unless --config is given.")

    plt.rcParams.update(FONT_PRESETS[args.font_preset])

    counts = _to_named_counts(args.counts, args.names)
    relabel = resolve_mapping(args.relabel)
    colors = resolve_palette(args.palette)

    plot_spec: Dict[str, Any] = {
        "y_label": args.y_label if args.y_label is not None else args.name,
        "x_label": args.x_label,
        "hide_yticks": args.hide_yticks,
        "figsize": tuple(args.figsize),
        "method": args.select,
        "threshold": args.threshold,
        "top_n": args.top_n,
        "reverse": args.reverse_samples,
        "outfile": f"{args.name}_biotype_composition.pdf",
    }
    if args.combined_only:
        plot_spec["sample_cols"] = [args.combined_name]

    result = run_analysis(
        name=args.name,
        counts=counts,
        gtf=args.gtf,
        output_dir=args.output_dir,
        combined=args.combined,
        combined_name=args.combined_name,
        relabel=relabel,
        colors=colors,
        plots=[plot_spec],
        id_column=args.id_column,
        sort_by=args.sort_by,
        prefer_dmode=not args.no_dmode,
        show=args.show,
    )
    logger.info(
        "Done. Table: %s | Figures: %s",
        result["table"],
        ", ".join(str(p) for p in result["figures"]),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
