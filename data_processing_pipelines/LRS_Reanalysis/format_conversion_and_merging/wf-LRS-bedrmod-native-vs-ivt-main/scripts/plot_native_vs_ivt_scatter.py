#!/usr/bin/env python3
"""
Scatter (hexbin) plots comparing a single column between native and IVT BEDRMod files.

For each site present in both files (inner join on chrom/start/end/strand), the
column value in IVT (x) is plotted against the value in native (y).

Outputs one PDF per call (parameterised by --column), with:
  - Page 1: all modifications combined
  - One additional page per modification type

TSV contains Pearson and Spearman correlation statistics, overall and per modification.

Column names are auto-detected from the #[...] header line.
"""

import ast
import gzip
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm


KEY_COLS = ["chrom", "chromStart", "chromEnd", "strand", "name"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("native",             help="Native BEDRMod file")
    p.add_argument("ivt",               help="IVT BEDRMod file")
    p.add_argument("--column",          required=True, help="Column to compare")
    p.add_argument("--output-pdf",      required=True)
    p.add_argument("--output-tsv",      required=True)
    p.add_argument("--min-native-coverage", type=int, default=0, metavar="N")
    p.add_argument("--min-ivt-coverage",    type=int, default=0, metavar="N")
    p.add_argument("--bins",            type=int, default=80,
                   help="Hexbin grid size (default: 80)")
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


def load_file(path, column, min_coverage):
    col_names, n_skip = detect_header(path)
    for c in (column, "name", "coverage"):
        if c not in col_names:
            raise ValueError(f"Column '{c}' not found in {path}. Available: {col_names}")

    dtypes = {
        "chrom":    "category",
        "chromStart": "int32",
        "chromEnd":   "int32",
        "strand":   "category",
        "name":     "category",
        "coverage": "int32",
        column:     "float32",
    }
    load_cols = KEY_COLS + ["coverage", column]
    df = pd.read_csv(
        path, sep="\t", skiprows=n_skip, header=None,
        names=col_names, usecols=load_cols, dtype=dtypes,
    )
    if min_coverage > 0:
        df = df[df["coverage"] > min_coverage]
    if column != "coverage":
        df = df.drop(columns=["coverage"])
    return df


def corr_stats(x, y):
    """Return (pearson_r, pearson_p, spearman_r, spearman_p) or NaNs if too few points."""
    if len(x) < 3:
        return np.nan, np.nan, np.nan, np.nan
    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)
    return pr, pp, sr, sp


def draw_scatter(ax, x, y, column, title, bins, use_log):
    n = len(x)
    if n == 0:
        ax.set_title(f"{title}  (no data)")
        return

    if use_log:
        x_plot = np.log10(np.clip(x, 1e-9, None))
        y_plot = np.log10(np.clip(y, 1e-9, None))
        label  = f"{column}  (log₁₀)"
    else:
        x_plot, y_plot = x, y
        label = column

    hb = ax.hexbin(x_plot, y_plot, gridsize=bins, cmap="YlOrRd",
                   mincnt=1, norm=LogNorm())
    plt.colorbar(hb, ax=ax, label="count")

    # Diagonal reference line
    lo = min(x_plot.min(), y_plot.min())
    hi = max(x_plot.max(), y_plot.max())
    ax.plot([lo, hi], [lo, hi], color="steelblue", linestyle="--",
            linewidth=0.9, alpha=0.7, label="y = x")

    pr, _, sr, _ = corr_stats(x_plot, y_plot)
    ax.set_xlabel(f"IVT  {label}")
    ax.set_ylabel(f"native  {label}")
    ax.set_title(f"{title}\nn={n:,}  |  Pearson r={pr:.3f}  |  Spearman r={sr:.3f}",
                 fontsize=9)
    ax.legend(fontsize=7)


def write_pdf(merged, column, native_stem, ivt_stem, bins, path):
    use_log = (
        merged[f"{column}_nat"].gt(0).all()
        and merged[f"{column}_ivt"].gt(0).all()
        and merged[f"{column}_nat"].max() > merged[f"{column}_nat"].min() * 100
    )

    mods = sorted(merged["name"].unique())
    n_pages = 1 + len(mods)

    with PdfPages(path) as pdf:
        # Overall page
        fig, ax = plt.subplots(figsize=(7, 6))
        fig.suptitle(f"{native_stem}  vs  {ivt_stem}  ·  {column}", fontsize=10)
        draw_scatter(
            ax,
            merged[f"{column}_ivt"].values,
            merged[f"{column}_nat"].values,
            column,
            f"all modifications  (n={len(merged):,})",
            bins, use_log,
        )
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Per-modification pages
        for mod in mods:
            grp = merged[merged["name"] == mod]
            fig, ax = plt.subplots(figsize=(7, 6))
            fig.suptitle(f"{native_stem}  vs  {ivt_stem}  ·  {column}  ·  {mod}",
                         fontsize=10)
            draw_scatter(
                ax,
                grp[f"{column}_ivt"].values,
                grp[f"{column}_nat"].values,
                column, mod, bins, use_log,
            )
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"[info] PDF written ({n_pages} pages): {path}", file=sys.stderr)


def write_tsv(merged, column, path):
    rows = [("modification", "n_sites",
             "pearson_r", "pearson_pval",
             "spearman_r", "spearman_pval")]

    def add_row(label, sub):
        x = sub[f"{column}_ivt"].values.astype(float)
        y = sub[f"{column}_nat"].values.astype(float)
        pr, pp, sr, sp = corr_stats(x, y)
        rows.append((label, len(sub),
                     f"{pr:.4f}", f"{pp:.3e}",
                     f"{sr:.4f}", f"{sp:.3e}"))

    add_row("all", merged)
    for mod, grp in merged.groupby("name", observed=True):
        add_row(mod, grp)

    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")

    print(f"[info] TSV written: {path}", file=sys.stderr)


def main():
    args = parse_args()

    nat_stem = args.native.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    ivt_stem = args.ivt.rsplit("/", 1)[-1].rsplit(".", 1)[0]

    print(f"[info] Loading native ({args.column}) ...", file=sys.stderr)
    nat = load_file(args.native, args.column, args.min_native_coverage)

    print(f"[info] Loading IVT ({args.column}) ...", file=sys.stderr)
    ivt = load_file(args.ivt, args.column, args.min_ivt_coverage)

    print("[info] Merging (inner join) ...", file=sys.stderr)
    merged = nat.merge(
        ivt[KEY_COLS + [args.column]],
        on=KEY_COLS,
        suffixes=("_nat", "_ivt"),
    )
    print(f"[info] Shared sites: {len(merged):,}", file=sys.stderr)

    write_tsv(merged, args.column, args.output_tsv)
    write_pdf(merged, args.column, nat_stem, ivt_stem, args.bins, args.output_pdf)

    print("[info] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
