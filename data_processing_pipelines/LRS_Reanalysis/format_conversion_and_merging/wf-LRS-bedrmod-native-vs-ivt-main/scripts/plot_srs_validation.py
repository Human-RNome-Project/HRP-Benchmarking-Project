#!/usr/bin/env python3
"""
Validate Nanopore modification predictions against an Illumina reference set.

Single-page PDF with one subplot per modification type.
Unfiltered set: sites with native coverage above the minimum threshold.
Filtered set: IVT-testable sites only (padj present).
Curves show Illumina overlap (%) vs site rank.

Overlap criterion: same chrom, chromStart, strand, and modification name.
"""

import ast
import gzip
import sys
import argparse
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams.update({
    "font.size":        16,
    "axes.titlesize":   18,
    "axes.labelsize":   16,
    "xtick.labelsize":  14,
    "ytick.labelsize":  14,
    "legend.fontsize":  14,
    "figure.dpi":       150,
})

_LINE_UNF = dict(linewidth=2.2, linestyle="--", alpha=0.90)  # unfiltered: dashed
_LINE_FLT = dict(linewidth=2.2, linestyle="-",  alpha=0.90)  # filtered:   solid

_MAX_COLS = 4


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--unfiltered", required=True,
                   help="Unfiltered native_vs_ivt_fisher BED (plain or .gz)")
    p.add_argument("--filtered",   required=True,
                   help="Filtered BED (plain or .gz)")
    p.add_argument("--illumina",   required=True,
                   help="Illumina validation BED")
    p.add_argument("--name",       required=True,
                   help="Biotype label for plot titles")
    p.add_argument(
        "--color", nargs=2, metavar=("MOD", "HEX"),
        action="append", default=[],
        help="Color for a modification (repeatable)",
    )
    p.add_argument("--min-native-coverage", type=int, default=30, metavar="N",
                   help="Minimum native coverage for unfiltered sites (default: 30)")
    p.add_argument("--output-pdf", required=True, help="Output PDF path")
    return p.parse_args()


# ── I/O ──────────────────────────────────────────────────────────────────────

def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _parse_header(path):
    """Return (col_names_or_None, n_header_lines).

    Understands two BEDRMod header styles:
      #['chrom', ...]          — Python-list format (sanitized)
      #chrom\tchromStart\t...  — BEDRMod v2 tab format (raw)
    """
    col_names = None
    n_skip = 0
    with _open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            n_skip += 1
            stripped = line.lstrip("#").strip()
            if stripped.startswith("["):
                try:
                    col_names = ast.literal_eval(stripped)
                except (ValueError, SyntaxError):
                    pass
            elif "\t" in stripped and not col_names:
                cols = stripped.split("\t")
                if cols[0] in ("chrom", "chromStart", "chr"):
                    col_names = cols
    return col_names, n_skip


def load_nanopore(path):
    needed = ["chrom", "chromStart", "strand", "name", "score", "frequency", "coverage", "padj"]
    col_names, n_skip = _parse_header(path)
    if col_names is None:
        raise ValueError(f"No column header found in {path}")
    available = [c for c in needed if c in col_names]
    df = pd.read_csv(
        path, sep="\t", skiprows=n_skip, header=None,
        names=col_names, usecols=available, dtype=str,
    )
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    if "frequency" in df.columns:
        df["frequency"] = pd.to_numeric(df["frequency"], errors="coerce")
    if "coverage" in df.columns:
        df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    if "padj" in df.columns:
        df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    return df


def load_illumina(path):
    """Return a set of 'chrom|chromStart|strand|name' match keys."""
    needed = ["chrom", "chromStart", "strand", "name"]
    col_names, n_skip = _parse_header(path)
    if col_names is None:
        col_names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
        n_skip = 0
    missing = [c for c in needed if c not in col_names]
    if missing:
        raise ValueError(f"Illumina BED missing columns: {missing}")
    df = pd.read_csv(
        path, sep="\t", skiprows=n_skip, header=None,
        names=col_names, usecols=needed, dtype=str,
        comment="#",
    )
    return set(df["chrom"] + "|" + df["chromStart"] + "|" + df["strand"] + "|" + df["name"])


def label_sites(df, illumina_keys):
    """Vectorised membership test; returns bool numpy array."""
    keys = df["chrom"] + "|" + df["chromStart"] + "|" + df["strand"] + "|" + df["name"]
    return keys.isin(illumina_keys).values


# ── Curve computation ─────────────────────────────────────────────────────────

def compute_cumprec(labels, scores):
    """Cumulative precision at each rank (score descending)."""
    if len(labels) == 0:
        return np.array([0]), np.array([0.0])
    order = np.argsort(scores)[::-1]
    lab   = labels[order]
    tp    = np.cumsum(lab)
    prec  = tp / np.arange(1, len(lab) + 1)
    return np.arange(1, len(lab) + 1), prec


# ── Plotting ──────────────────────────────────────────────────────────────────

def _clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _plot_mod_ax(ax, color, labels_u, freqs_u, labels_f, neglog_padj_f, n_illumina):
    """Fill one axes with the two curves for a single modification."""
    if n_illumina == 0:
        ax.text(0.5, 0.5, "no Illumina reference\nsites for this biotype",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=12, color="#999999", style="italic")
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)
        return

    n_flt = max(len(labels_f), 1)

    for labels, ranks_col, line_kw in [
        (labels_u, freqs_u,       _LINE_UNF),
        (labels_f, neglog_padj_f, _LINE_FLT),
    ]:
        if len(labels) == 0:
            continue
        ranks, cp = compute_cumprec(labels, ranks_col)
        ax.plot(ranks, cp * 100, color=color, **line_kw)

    ax.set_xlim(1, n_flt)
    ax.set_ylim(bottom=0)

    ax.set_xlabel("site rank", labelpad=6)
    ax.set_ylabel("Illumina overlap (%)", labelpad=6)
    _clean(ax)

    ax.text(
        0.02, 0.04,
        f"Illumina: {n_illumina:,}  |  unf: {len(labels_u):,}  |  flt: {len(labels_f):,}",
        transform=ax.transAxes, fontsize=9,
        ha="left", va="bottom", color="#555555",
    )


def plot_all_mods(fig, biotype, mods_data):
    """
    Draw all modifications onto a grid of subplots in *fig*.

    mods_data: list of (mod, color, labels_u, freqs_u, labels_f, neglog_padj_f, n_illumina)
    """
    n = len(mods_data)
    n_cols = min(n, _MAX_COLS)
    n_rows = (n + n_cols - 1) // n_cols

    axes = fig.subplots(n_rows, n_cols, squeeze=False)

    for idx, (mod, color, labels_u, freqs_u, labels_f, neglog_padj_f, n_ill) in \
            enumerate(mods_data):
        row, col = divmod(idx, n_cols)
        ax = axes[row][col]
        ax.set_title(mod, pad=8)
        _plot_mod_ax(ax, color, labels_u, freqs_u, labels_f, neglog_padj_f, n_ill)

    for idx in range(n, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row][col].set_visible(False)

    # Figure-level title + subtitle
    fig.suptitle(
        f"{biotype}  —  Illumina overlap vs site rank",
        fontsize=18, y=0.99,
    )

    # Shared legend at the bottom using neutral-colour proxy lines
    leg_handles = [
        plt.Line2D([0], [0], color="#555555", linewidth=2.2, linestyle="--", alpha=0.90,
                   label="unfiltered  (ranked by frequency)"),
        plt.Line2D([0], [0], color="#555555", linewidth=2.2, linestyle="-",  alpha=0.90,
                   label="filtered    (ranked by −log₁₀(padj))"),
    ]
    fig.legend(handles=leg_handles, loc="lower center", ncol=2,
               frameon=False, fontsize=13, bbox_to_anchor=(0.5, 0.0))


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args      = parse_args()
    color_map = {mod: hex_ for mod, hex_ in args.color}

    print(f"[info] Loading Illumina: {args.illumina}", file=sys.stderr)
    illumina_keys = load_illumina(args.illumina)
    print(f"  {len(illumina_keys):,} sites", file=sys.stderr)

    print(f"[info] Loading unfiltered: {args.unfiltered}", file=sys.stderr)
    df_unf = load_nanopore(args.unfiltered)
    n_total_unf = len(df_unf)
    # Chroms interrogated by this biotype's assay.  Biotypes occupy disjoint
    # chrom namespaces (genomic for polyA, hs_rRNA_* for rRNA, hs_tRNA_* for
    # tRNA), so this restricts the Illumina reference to the matching biotype.
    biotype_chroms = set(df_unf["chrom"].unique())
    df_unf = df_unf[df_unf["coverage"] > args.min_native_coverage].copy()
    print(f"  {n_total_unf:,} sites total → {len(df_unf):,} "
          f"with native coverage > {args.min_native_coverage}", file=sys.stderr)

    print(f"[info] Loading filtered: {args.filtered}", file=sys.stderr)
    df_flt_all = load_nanopore(args.filtered)
    df_flt = df_flt_all[df_flt_all["padj"].notna()].copy()
    # The most significant sites can have padj underflow to exactly 0 (the
    # smallest positive padj is already denormal, so flooring padj itself
    # underflows).  Work in -log10 space instead: compute it for positive padj
    # and rank the zeros just above the largest finite value so they sit at the
    # top rather than being dropped.
    padj    = df_flt["padj"].to_numpy(dtype=float)
    neglog  = np.empty(len(padj))
    posmask = padj > 0
    neglog[posmask] = -np.log10(padj[posmask])
    max_finite = neglog[posmask].max() if posmask.any() else 0.0
    neglog[~posmask] = max_finite + 1.0
    df_flt["neglog_padj"] = neglog
    print(f"  {len(df_flt_all):,} sites total → {len(df_flt):,} IVT-testable", file=sys.stderr)

    print("[info] Labelling ...", file=sys.stderr)
    df_unf["match"] = label_sites(df_unf, illumina_keys)
    df_flt["match"] = label_sites(df_flt, illumina_keys)

    # Illumina reference counts per mod, restricted to this biotype's chroms.
    # A mod can be present genome-wide but have zero sites in this biotype;
    # such mods are still shown (with a placeholder) rather than dropped.
    illumina_counts = Counter(
        k.split("|")[3] for k in illumina_keys if k.split("|")[0] in biotype_chroms
    )
    illumina_mods = set(k.split("|")[3] for k in illumina_keys)
    nanopore_mods = set(df_unf["name"].unique())
    mods          = sorted(illumina_mods & nanopore_mods)

    if not mods:
        print("[warn] No shared modification types — writing empty PDF.", file=sys.stderr)
        with PdfPages(args.output_pdf) as pdf:
            fig, ax = plt.subplots(figsize=(6, 3))
            ax.text(0.5, 0.5, "No shared modification types",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=14, color="grey")
            ax.axis("off")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
        return

    print(f"[info] Mods to validate: {', '.join(mods)}", file=sys.stderr)

    mods_data = []
    for mod in mods:
        sub_u     = df_unf[df_unf["name"] == mod]
        sub_f     = df_flt[df_flt["name"] == mod]
        sub_f_all = df_flt_all[df_flt_all["name"] == mod]
        n_ill = illumina_counts[mod]   # biotype-restricted Illumina reference sites
        print(f"  {mod}: illumina={n_ill:,}  unf={len(sub_u):,}"
              f"  flt={len(sub_f_all):,} (testable={len(sub_f):,})"
              f"  TP_flt={int(sub_f['match'].sum())}",
              file=sys.stderr)
        mods_data.append((
            mod, color_map.get(mod, "#444444"),
            sub_u["match"].values, sub_u["frequency"].values,
            sub_f["match"].values, sub_f["neglog_padj"].values,
            n_ill,
        ))

    n      = len(mods_data)
    n_cols = min(n, _MAX_COLS)
    n_rows = (n + n_cols - 1) // n_cols
    fig    = plt.figure(figsize=(4.5 * n_cols, 5.5 * n_rows))
    plot_all_mods(fig, args.name, mods_data)
    fig.tight_layout(rect=[0, 0.04, 1, 0.99])

    with PdfPages(args.output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"[info] Written → {args.output_pdf}", file=sys.stderr)


if __name__ == "__main__":
    main()
