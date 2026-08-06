#!/usr/bin/env python3
"""
Volcano plot: log2(native frequency / IVT frequency) vs −log10(padj).

Input: native_vs_ivt_fisher BEDRMod file (output of filter_native_vs_ivt.py).
Only testable sites (padj ≠ NA) are plotted.

Outputs:
  - PDF: one overall page + one page per modification type (hexbin volcano).
  - TSV: site counts at common padj thresholds, overall and per modification.

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
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm


PADJ_THRESHOLDS = [0.001, 0.01, 0.05, 0.1, 0.2]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",        help="native_vs_ivt_fisher BEDRMod file")
    p.add_argument("--output-pdf", required=True)
    p.add_argument("--output-tsv", required=True)
    p.add_argument(
        "--y-cap", type=float, default=30.0, metavar="N",
        help="Cap -log10(padj) at this value for display (default: 30); "
             "y-axis floor is fixed at padj=0.1 (-log10=1)",
    )
    p.add_argument("--bins", type=int, default=100,
                   help="Hexbin grid size (default: 100)")
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


def load_data(path):
    col_names, n_skip = detect_header(path)
    for col in ("name", "frequency", "ivt_frequency", "padj"):
        if col not in col_names:
            raise ValueError(
                f"Column '{col}' not found in {path}.\nAvailable: {col_names}"
            )

    df = pd.read_csv(
        path, sep="\t", skiprows=n_skip, header=None,
        names=col_names,
        usecols=["name", "frequency", "ivt_frequency", "padj"],
        dtype=str,
    )
    df = df[df["name"] != "name"].reset_index(drop=True)

    df["frequency"]     = pd.to_numeric(df["frequency"],     errors="coerce")
    df["ivt_frequency"] = pd.to_numeric(df["ivt_frequency"], errors="coerce")
    df["padj"]          = pd.to_numeric(df["padj"],          errors="coerce")

    n_total      = len(df)
    n_ivt_absent = int(df["padj"].isna().sum())

    testable = df.dropna(subset=["padj", "frequency", "ivt_frequency"]).copy()

    # Remove pvalue=1 (uninformative, dominate colour scale) and pvalue=0
    # (floating-point underflow — handled separately via y-cap).
    n_pv1 = int((testable["padj"] == 1).sum())
    testable = testable[testable["padj"] < 1].copy()

    # Drop sites where log2FC is undefined (zero in either frequency).
    n_before = len(testable)
    testable = testable[
        (testable["frequency"] > 0) & (testable["ivt_frequency"] > 0)
    ].copy()
    n_undef = n_before - len(testable)

    testable["log2fc"]        = np.log2(testable["frequency"] / testable["ivt_frequency"])
    testable["neg_log_padj"]  = -np.log10(testable["padj"].clip(lower=1e-300))
    testable["name"]          = testable["name"].astype("category")

    print(
        f"[info] Total sites: {n_total:,} | "
        f"Testable: {len(testable):,} | "
        f"IVT-absent (excluded): {n_ivt_absent:,} | "
        f"padj=1 (excluded): {n_pv1:,} | "
        f"log2FC undefined (excluded): {n_undef:,}",
        file=sys.stderr,
    )
    return testable


def draw_volcano(ax, x, y_raw, y_cap, title, bins):
    y = np.minimum(y_raw, y_cap)
    n = len(x)

    if n == 0:
        ax.set_title(f"{title}  (no data)")
        return

    hb = ax.hexbin(x, y, gridsize=bins, cmap="YlOrRd", mincnt=1, norm=LogNorm())
    cb = plt.colorbar(hb, ax=ax)
    cb.set_label("count", fontsize=12)
    cb.ax.tick_params(labelsize=11)

    yt_05 = -np.log10(0.05)
    if yt_05 <= y_cap:
        ax.axhline(yt_05, color="#ff7f0e", linestyle="--", linewidth=2.0,
                   alpha=0.9, label="padj = 0.05")

    ax.axvline(0, color="grey", linestyle="-", linewidth=0.8, alpha=0.4)

    n_clipped = int((y_raw > y_cap).sum())
    if n_clipped:
        ax.text(0.98, 0.98, f"{n_clipped:,} sites above cap",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=12, color="grey")

    x_lo, x_hi = np.percentile(x, [0.5, 99.5])
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(-np.log10(0.1), y_cap)
    ax.set_xlabel("log₂(native frequency / IVT frequency)", fontsize=14)
    ax.set_ylabel(f"−log₁₀(padj)  [padj < 0.1 only, capped at {y_cap}]", fontsize=14)
    ax.set_title(f"{title}\nn={n:,}", fontsize=14)
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=12, loc="upper left")


def write_pdf(df, stem, y_cap, bins, path):
    mods = sorted(df["name"].unique())
    n_pages = 1 + len(mods)

    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(9, 7))
        fig.suptitle(f"{stem}  ·  volcano  (all modifications)", fontsize=14)
        draw_volcano(ax, df["log2fc"].values, df["neg_log_padj"].values,
                     y_cap, "all modifications", bins)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for mod in mods:
            grp = df[df["name"] == mod]
            fig, ax = plt.subplots(figsize=(9, 7))
            fig.suptitle(f"{stem}  ·  volcano  ·  {mod}", fontsize=14)
            draw_volcano(ax, grp["log2fc"].values, grp["neg_log_padj"].values,
                         y_cap, mod, bins)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"[info] PDF written ({n_pages} pages): {path}", file=sys.stderr)


def write_tsv(df, path):
    header = ("modification", "n_testable",
              *[f"padj<{t}" for t in PADJ_THRESHOLDS])

    def make_row(label, sub):
        counts = tuple(int((sub["padj"] < t).sum()) for t in PADJ_THRESHOLDS)
        return (label, len(sub), *counts)

    rows = [header, make_row("all", df)]
    for mod, grp in df.groupby("name", observed=True):
        rows.append(make_row(mod, grp))

    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")

    print(f"[info] TSV written: {path}", file=sys.stderr)


def main():
    args = parse_args()

    stem = args.input.rsplit("/", 1)[-1]
    print(f"[info] Reading {args.input} ...", file=sys.stderr)
    df = load_data(args.input)
    mods = sorted(df["name"].unique())
    print(f"[info] Modifications: {', '.join(mods)}", file=sys.stderr)

    write_tsv(df, args.output_tsv)
    write_pdf(df, stem, args.y_cap, args.bins, args.output_pdf)

    print("[info] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
