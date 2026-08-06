#!/usr/bin/env python3
"""
Plot the value distribution of a single column from a BEDRMod file.

Outputs:
  - TSV: percentile table overall and per modification.
  - PDF: 2×2 figure — cumulative fraction (ECDF), non-cumulative fraction,
         cumulative count, and non-cumulative count. One curve per
         modification plus an overall curve. For position-level columns
         (e.g. coverage) the overall curve is deduplicated by position.

Column names are auto-detected from the #[...] header line.
"""

import ast
import gzip
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages


# Columns stored once per genomic position, not once per modification.
# The overall curves for these use deduplicated positions.
POSITION_LEVEL_COLS = {"coverage", "std_coverage"}

KEY_COLS    = ["chrom", "chromStart", "chromEnd", "strand"]
PERCENTILES = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 99.9, 100]
ECDF_POINTS = 2000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",        help="Input BEDRMod file")
    p.add_argument("--column",     required=True, help="Column to plot")
    p.add_argument("--output-tsv", required=True, help="Output TSV file")
    p.add_argument("--output-pdf", required=True, help="Output PDF file")
    p.add_argument("--bins",       type=int, default=100,
                   help="Histogram bins (default: 100)")
    p.add_argument("--color", nargs=2, metavar=("MOD", "HEX"),
                   action="append", default=[],
                   help="Colour for a modification (repeat for each mod)")
    p.add_argument("--vline", type=float, default=None,
                   help="Draw a dashed vertical reference line at this x value")
    p.add_argument("--xmin", type=float, default=None,
                   help="Lower limit for the x-axis")
    return p.parse_args()


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def detect_header(path):
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
    if col_names is None:
        raise ValueError(f"Could not detect column names from header of {path}")
    return col_names, n_skip


def load_data(path, column):
    col_names, n_skip = detect_header(path)
    if column not in col_names:
        raise ValueError(
            f"Column '{column}' not found in {path}.\nAvailable: {col_names}"
        )
    load_cols = [c for c in KEY_COLS + ["name", column] if c in col_names]
    df = pd.read_csv(
        path, sep="\t", skiprows=n_skip, header=None,
        names=col_names, usecols=load_cols,
        dtype={column: "float32", "name": "category",
               "chrom": "category", "strand": "category"},
    )
    return df.dropna(subset=[column])


def ecdf_xy(arr):
    q = np.linspace(0, 1, min(ECDF_POINTS, len(arr)))
    return np.quantile(arr, q), q


def make_bin_edges(arr, n_bins, use_log):
    lo, hi = arr.min(), arr.max()
    if use_log:
        lo = max(lo, 1e-9)
        return np.logspace(np.log10(lo), np.log10(hi), n_bins + 1)
    return np.linspace(lo, hi, n_bins + 1)


def write_tsv(df, column, is_pos_level, path):
    def block(arr, label):
        rows = [("section", label, ""), ("stat", "value", "")]
        for p in PERCENTILES:
            name = f"p{p:g}" if p not in (0, 100) else ("min" if p == 0 else "max")
            rows.append((name, f"{np.percentile(arr, p):.4g}", ""))
        rows += [
            ("mean",  f"{arr.mean():.4g}", ""),
            ("std",   f"{arr.std():.4g}",  ""),
            ("count", str(len(arr)),        ""),
            ("", "", ""),
        ]
        return rows

    all_arr = (
        df.drop_duplicates(subset=KEY_COLS)[column].values.astype(float)
        if is_pos_level else df[column].values.astype(float)
    )
    rows = block(all_arr, "all_positions" if is_pos_level else "all")
    for mod, grp in df.groupby("name", observed=True):
        rows += block(grp[column].values.astype(float), str(mod))

    with open(path, "w") as fh:
        for r in rows:
            if r[2]:
                fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")
            elif r[0]:
                fh.write(f"{r[0]}\t{r[1]}\n")
            else:
                fh.write("\n")


def _fmt(ax, column, use_log):
    if use_log:
        ax.set_xscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel(column, fontsize=11)
    ax.tick_params(labelsize=10)


def draw_cumulative(axes, groups, bin_edges, use_log, column):
    """Page 1 — left: ECDF, right: cumulative count."""
    for label, color, arr, is_all in groups:
        lw = 2.5 if is_all else 1.5
        ls = "--" if is_all else "-"
        n  = len(arr)
        x, y_frac = ecdf_xy(arr)
        axes[0].plot(x, y_frac,     color=color, lw=lw, ls=ls,
                     label=f"{label}  (n={n:,})")
        axes[1].plot(x, y_frac * n, color=color, lw=lw, ls=ls)

    for ax in axes:
        _fmt(ax, column, use_log)

    axes[0].set_ylabel("cumulative fraction", fontsize=11)
    axes[0].set_ylim(0, 1)
    axes[0].set_title("ECDF", fontsize=11, fontweight="bold")

    axes[1].set_ylabel("cumulative count", fontsize=11)
    axes[1].set_title("cumulative count", fontsize=11, fontweight="bold")
    axes[1].yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{int(x):,}")
    )


def draw_noncumulative(axes, groups, bin_edges, use_log, column):
    """Page 2 — left: fraction per bin, right: count per bin."""
    for label, color, arr, is_all in groups:
        lw = 2.5 if is_all else 1.5
        ls = "--" if is_all else "-"
        n  = len(arr)
        counts, _ = np.histogram(arr, bins=bin_edges)
        axes[0].stairs(counts / n, bin_edges, color=color, lw=lw, ls=ls,
                       label=f"{label}  (n={n:,})")
        axes[1].stairs(counts,     bin_edges, color=color, lw=lw, ls=ls)

    for ax in axes:
        _fmt(ax, column, use_log)

    axes[0].set_ylabel("fraction per bin", fontsize=11)
    axes[0].set_title("distribution (fraction)", fontsize=11, fontweight="bold")

    axes[1].set_ylabel("count per bin", fontsize=11)
    axes[1].set_title("count (histogram)", fontsize=11, fontweight="bold")
    axes[1].yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{int(x):,}")
    )


def write_pdf(df, column, is_pos_level, color_map, n_bins, input_path, path,
              vline=None, xmin=None):
    stem = input_path.rsplit("/", 1)[-1]
    mods = sorted(df["name"].astype(str).unique())

    # Fallback palette for mods not in the user-supplied color map
    cmap = plt.get_cmap("tab10")
    fallback_idx = 0
    resolved = {}
    for mod in mods:
        if mod in color_map:
            resolved[mod] = color_map[mod]
        else:
            resolved[mod] = cmap(fallback_idx % 10)
            fallback_idx += 1

    use_log = is_pos_level or (
        df[column].min() > 0 and df[column].max() > df[column].min() * 100
    )

    all_arr = (
        df.drop_duplicates(subset=KEY_COLS)[column].values.astype(float)
        if is_pos_level else df[column].values.astype(float)
    )
    all_label = "all  [pos-deduped]" if is_pos_level else "all"

    # Build group list: per-mod first, then "all" on top
    groups = [
        (mod, resolved[mod],
         df[df["name"].astype(str) == mod][column].values.astype(float),
         False)
        for mod in mods
    ]
    groups.append((all_label, "black", all_arr, True))

    bin_edges = make_bin_edges(all_arr, n_bins, use_log)

    def _save_page(pdf, draw_fn, suptitle):
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle(suptitle, fontsize=12)
        draw_fn(axes, groups, bin_edges, use_log, column)
        for ax in axes:
            if xmin is not None:
                ax.set_xlim(left=xmin)
            if vline is not None:
                ax.axvline(vline, color="#888888", lw=1.3, ls="--", zorder=1)
                ax.text(vline, 1.0, f" ≥ {vline:g}",
                        transform=ax.get_xaxis_transform(),
                        va="top", ha="left", fontsize=9, color="#666666")
        handles, labels = axes[0].get_legend_handles_labels()
        ncol = min(len(handles), 5)
        fig.legend(handles, labels, loc="lower center", ncol=ncol,
                   fontsize=10, frameon=False, bbox_to_anchor=(0.5, 0.0))
        fig.tight_layout(rect=[0, 0.10, 1, 1])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    with PdfPages(path) as pdf:
        _save_page(pdf, draw_cumulative,    f"{stem}  ·  {column}  —  cumulative")
        _save_page(pdf, draw_noncumulative, f"{stem}  ·  {column}  —  distribution")


def main():
    args = parse_args()
    is_pos_level = args.column in POSITION_LEVEL_COLS
    color_map = {mod: hex_ for mod, hex_ in args.color}

    print(f"[info] Reading '{args.column}' from {args.input} ...", file=sys.stderr)
    df = load_data(args.input, args.column)
    mods = sorted(df["name"].astype(str).unique())
    print(
        f"[info] {len(df):,} rows | modifications: {', '.join(mods)}",
        file=sys.stderr,
    )
    if is_pos_level:
        n_pos = len(df.drop_duplicates(subset=KEY_COLS))
        print(
            f"[info] Position-level column — overall curves use "
            f"{n_pos:,} unique positions",
            file=sys.stderr,
        )

    print(f"[info] Writing TSV → {args.output_tsv}", file=sys.stderr)
    write_tsv(df, args.column, is_pos_level, args.output_tsv)

    print(f"[info] Writing PDF → {args.output_pdf}", file=sys.stderr)
    write_pdf(df, args.column, is_pos_level, color_map,
              args.bins, args.input, args.output_pdf,
              vline=args.vline, xmin=args.xmin)

    print("[info] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
