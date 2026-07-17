#!/usr/bin/env python3
"""
Plot p-value and adjusted p-value distributions from a native-vs-IVT BEDRMod
output file (produced by filter_native_vs_ivt.py).

Outputs:
  - TSV: summary counts (total, testable, IVT-absent; sites passing each
         significance threshold) — overall and per modification.
  - PDF (multi-page): one page per modification + an overall summary page.
         Each page has three panels:
           (1) raw p-value histogram with uniform expectation line,
           (2) adjusted p-value (BH) histogram,
           (3) native frequency vs -log10(padj) hexbin.

Column names are auto-detected from the #[...] header line or, as a fallback,
from the first plain tab-separated header line.
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
from matplotlib.colors import LogNorm


THRESHOLDS = [0.001, 0.01, 0.05, 0.1, 0.2]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",        help="native_vs_ivt_fisher BEDRMod file")
    p.add_argument("--output-tsv", required=True, help="Output TSV file")
    p.add_argument("--output-pdf", required=True, help="Output PDF file")
    p.add_argument("--bins", type=int, default=50,
                   help="Number of histogram bins for p-values (default: 50)")
    return p.parse_args()


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def detect_header(path):
    """Return (col_names, n_skip).

    Accepts:
      1. BEDRMod #[...] comment line (preferred).
      2. Plain tab-separated first non-comment line (fallback).
    """
    col_names = None
    n_skip = 0
    with _open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                if col_names is None:
                    col_names = line.rstrip("\n").split("\t")
                    n_skip += 1
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


def threshold_rows(pval_clean, padj_clean, n_testable):
    rows = []
    for t in THRESHOLDS:
        n_pv = int((pval_clean < t).sum())
        n_pj = int((padj_clean < t).sum())
        rows.append((
            f"<{t}",
            f"{n_pv} ({n_pv/n_testable:.4f})",
            f"{n_pj} ({n_pj/n_testable:.4f})",
        ))
    return rows


def write_tsv(df, total_sites, path):
    pval = pd.to_numeric(df["pvalue"],    errors="coerce")
    padj = pd.to_numeric(df["padj"],      errors="coerce")

    n_testable   = int(pval.notna().sum())
    n_ivt_absent = total_sites - n_testable

    lines = [
        ("section", "overall", ""),
        ("stat",    "count",   "fraction_of_total"),
        ("total_sites",  str(total_sites),  "1.0000"),
        ("testable",     str(n_testable),   f"{n_testable/total_sites:.4f}"),
        ("ivt_absent",   str(n_ivt_absent), f"{n_ivt_absent/total_sites:.4f}"),
        ("", "", ""),
        ("threshold", "pvalue_lt", "padj_lt"),
        *threshold_rows(pval.dropna(), padj.dropna(), n_testable),
    ]

    for mod, grp in df.groupby("name"):
        pv = pd.to_numeric(grp["pvalue"], errors="coerce")
        pj = pd.to_numeric(grp["padj"],   errors="coerce")
        n_tot  = len(grp)
        n_test = int(pv.notna().sum())
        lines += [
            ("", "", ""),
            ("section", str(mod), ""),
            ("stat", "count", "fraction_of_mod_total"),
            ("total_sites",  str(n_tot),  "1.0000"),
            ("testable",     str(n_test), f"{n_test/n_tot:.4f}"),
            ("ivt_absent",   str(n_tot - n_test), f"{(n_tot-n_test)/n_tot:.4f}"),
            ("", "", ""),
            ("threshold", "pvalue_lt", "padj_lt"),
            *threshold_rows(pv.dropna(), pj.dropna(), max(n_test, 1)),
        ]

    with open(path, "w") as fh:
        for row in lines:
            if row[2]:
                fh.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")
            elif row[0]:
                fh.write(f"{row[0]}\t{row[1]}\n")
            else:
                fh.write("\n")


XMAX = 0.9  # x-axis upper limit for p-value / padj plots


def _cumcount(values, x_max):
    """Return (x, y) where y is the cumulative count of values <= x, capped at x_max."""
    x = np.sort(values[values <= x_max])
    y = np.arange(1, len(x) + 1, dtype=np.int64)
    return x, y


def draw_page(pval, padj, freq, score, title, n_bins, fig):
    """Fill a 2×3 figure for one group of sites.

    Layout:
      [0,0] p-value histogram (x ≤ 0.9)
      [0,1] padj histogram    (x ≤ 0.9)
      [0,2] frequency vs −log10(padj)
      [1,0] p-value ECDF      (x ≤ 0.9)
      [1,1] padj ECDF         (x ≤ 0.9)
      [1,2] score — IVT-absent sites
    """
    axes = fig.subplots(2, 3)
    fig.suptitle(title, fontsize=11)

    pval_clean = pval.dropna().values
    padj_clean = padj.dropna().values
    n_testable = len(pval_clean)

    # bins limited to [0, XMAX]
    bin_edges = np.linspace(0, XMAX, n_bins + 1)

    # --- [0,0]: p-value histogram (x ≤ XMAX) ---
    ax = axes[0, 0]
    ax.hist(pval_clean[pval_clean <= XMAX], bins=bin_edges,
            color="#4C72B0", edgecolor="none", alpha=0.85)
    n_shown = int((pval_clean <= XMAX).sum())
    expected = n_shown / n_bins
    ax.axhline(expected, color="#d62728", linestyle="--",
               linewidth=1.2, label=f"uniform ({expected:,.0f})")
    for t in [0.05, 0.1]:
        ax.axvline(t, color="grey", linestyle=":", linewidth=0.8, alpha=0.7)
    ax.set_xlim(0, XMAX)
    ax.set_xlabel("p-value")
    ax.set_ylabel("count")
    ax.set_title(f"p-value histogram  (n={n_testable:,})")
    ax.legend(fontsize=7)

    # --- [0,1]: padj histogram (x ≤ XMAX) ---
    ax = axes[0, 1]
    ax.hist(padj_clean[padj_clean <= XMAX], bins=bin_edges,
            color="#55A868", edgecolor="none", alpha=0.85)
    for t in [0.05, 0.1]:
        ax.axvline(t, color="grey", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.text(t + 0.01, ax.get_ylim()[1] * 0.97, f"{t}",
                fontsize=7, va="top", color="grey")
    n_sig = int((padj_clean < 0.05).sum())
    ax.set_xlim(0, XMAX)
    ax.set_xlabel("padj  (BH)")
    ax.set_ylabel("count")
    ax.set_title(f"padj histogram  (n={n_testable:,},  <0.05: {n_sig:,})")

    # --- [0,2]: frequency vs -log10(padj) ---
    ax = axes[0, 2]
    mask = padj.notna() & freq.notna()
    x = freq[mask].values
    y = -np.log10(padj[mask].values.clip(1e-300))
    if len(x) > 1:
        hb = ax.hexbin(x, y, gridsize=80, cmap="YlOrRd",
                       mincnt=1, norm=LogNorm())
        fig.colorbar(hb, ax=ax, label="count")
    ax.axhline(-np.log10(0.05), color="#d62728", linestyle="--",
               linewidth=0.9, alpha=0.8, label="padj=0.05")
    ax.set_xlabel("native frequency (%)")
    ax.set_ylabel("−log₁₀(padj)")
    ax.set_title("frequency vs significance")
    ax.legend(fontsize=7)

    # --- [1,0]: p-value cumulative count (x ≤ XMAX) ---
    ax = axes[1, 0]
    if n_testable > 0:
        x_pv, y_pv = _cumcount(pval_clean, XMAX)
        ax.plot(x_pv, y_pv, color="#4C72B0", linewidth=1.2)
    for t in [0.05, 0.1]:
        ax.axvline(t, color="grey", linestyle=":", linewidth=0.8, alpha=0.7)
        n_lt = int((pval_clean < t).sum())
        ax.annotate(f"{n_lt:,}", xy=(t, n_lt), xytext=(t + 0.01, n_lt),
                    fontsize=7, color="grey", va="bottom")
    ax.set_xlim(0, XMAX)
    ax.set_xlabel("p-value")
    ax.set_ylabel("cumulative count")
    ax.set_title("p-value cumulative count")

    # --- [1,1]: padj cumulative count (x ≤ XMAX) ---
    ax = axes[1, 1]
    if n_testable > 0:
        x_pj, y_pj = _cumcount(padj_clean, XMAX)
        ax.plot(x_pj, y_pj, color="#55A868", linewidth=1.2)
    for t in [0.05, 0.1]:
        ax.axvline(t, color="grey", linestyle=":", linewidth=0.8, alpha=0.7)
        n_lt = int((padj_clean < t).sum())
        ax.annotate(f"{n_lt:,}", xy=(t, n_lt), xytext=(t + 0.01, n_lt),
                    fontsize=7, color="grey", va="bottom")
    ax.set_xlim(0, XMAX)
    ax.set_xlabel("padj  (BH)")
    ax.set_ylabel("cumulative count")
    ax.set_title("padj cumulative count")

    # --- [1,2]: score distribution for IVT-absent sites ---
    ax = axes[1, 2]
    ivt_absent_score = score[pval.isna()].dropna().values
    n_absent = len(ivt_absent_score)
    if n_absent > 0:
        ax.hist(ivt_absent_score, bins=n_bins, range=(0, 1),
                color="#8172B2", edgecolor="none", alpha=0.85)
    ax.set_xlabel("score")
    ax.set_ylabel("count")
    ax.set_title(f"score  —  IVT-absent sites  (n={n_absent:,})")

    fig.tight_layout()


def write_pdf(df, stem, n_bins, path):
    mods = sorted(df["name"].unique())

    with PdfPages(path) as pdf:
        # Overall summary page
        fig = plt.figure(figsize=(15, 9))
        draw_page(
            pd.to_numeric(df["pvalue"],    errors="coerce"),
            pd.to_numeric(df["padj"],      errors="coerce"),
            pd.to_numeric(df["frequency"], errors="coerce"),
            pd.to_numeric(df["score"],     errors="coerce"),
            f"{stem}  ·  all modifications  (n={len(df):,})",
            n_bins, fig,
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # One page per modification
        for mod in mods:
            grp = df[df["name"] == mod]
            fig = plt.figure(figsize=(15, 9))
            draw_page(
                pd.to_numeric(grp["pvalue"],    errors="coerce"),
                pd.to_numeric(grp["padj"],      errors="coerce"),
                pd.to_numeric(grp["frequency"], errors="coerce"),
                pd.to_numeric(grp["score"],     errors="coerce"),
                f"{stem}  ·  {mod}  (n={len(grp):,})",
                n_bins, fig,
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def main():
    args = parse_args()

    col_names, n_skip = detect_header(args.input)
    for col in ("pvalue", "padj", "frequency", "score", "name"):
        if col not in col_names:
            raise ValueError(
                f"Column '{col}' not found in {args.input}.\nAvailable: {col_names}"
            )

    stem = args.input.rsplit("/", 1)[-1]

    print(f"[info] Reading {args.input} ...", file=sys.stderr)
    df = pd.read_csv(
        args.input,
        sep="\t",
        skiprows=n_skip,
        header=None,
        names=col_names,
        usecols=["name", "score", "pvalue", "padj", "frequency"],
        dtype=str,
    )
    # Drop any stray header row that survived as data (plain-header fallback edge case)
    df = df[df["name"] != "name"].reset_index(drop=True)

    total_sites = len(df)
    mods = sorted(df["name"].unique())
    print(
        f"[info] {total_sites:,} sites | modifications: {', '.join(mods)}",
        file=sys.stderr,
    )

    print(f"[info] Writing TSV → {args.output_tsv}", file=sys.stderr)
    write_tsv(df, total_sites, args.output_tsv)

    print(f"[info] Writing PDF → {args.output_pdf}", file=sys.stderr)
    write_pdf(df, stem, args.bins, args.output_pdf)

    print("[info] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
