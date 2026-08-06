#!/usr/bin/env python3
"""
Score density plot: KDE of BEDRMod score across four site groups.

Page 1 (summary) — all modifications overlaid on one axis:
  solid lines  : native, significant  (padj < --sig-padj)
  dashed lines : IVT, all
  x-axis restricted to scores > 0.5; each mod coloured by --color.

Remaining pages (overall + per modification):
  1. native, significant  — native score for sites with padj < --sig-padj
  2. IVT, significant     — IVT score for those same sites (from IVT BEDRMod)
  3. IVT, all             — IVT score for all sites in the IVT BEDRMod file
  4. native, IVT-absent   — native score for IVT-absent sites (padj = NA)

Column names are auto-detected from the #[...] header line.
"""

import ast
import gzip
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


KEY_COLS = ["chrom", "chromStart", "chromEnd", "strand", "name"]

GROUPS_STYLE = {
    "native, significant": "#d62728",
    "IVT, significant":    "#ff7f0e",
    "IVT, all":            "#1f77b4",
    "native, IVT-absent":  "#7f7f7f",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fisher",       help="native_vs_ivt_fisher BEDRMod file")
    p.add_argument("ivt",          help="IVT BEDRMod file")
    p.add_argument("--output-pdf", required=True)
    p.add_argument(
        "--sig-padj", type=float, default=0.01,
        help="padj threshold for 'significant' group (default: 0.01)",
    )
    p.add_argument(
        "--color", nargs=2, metavar=("MOD", "HEX"),
        action="append", dest="colors", default=[],
        help="Modification colour for the summary page (repeatable)",
    )
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


def load_fisher(path, sig_padj):
    """Return (sig_native, na_native) — KEY_COLS + score subsets."""
    col_names, n_skip = detect_header(path)
    for c in KEY_COLS + ["score", "padj"]:
        if c not in col_names:
            raise ValueError(f"Column '{c}' not found in {path}.")

    df = pd.read_csv(
        path, sep="\t", skiprows=n_skip, header=None,
        names=col_names,
        usecols=KEY_COLS + ["score", "padj"],
        dtype={
            "chrom": "category", "strand": "category", "name": "category",
            "chromStart": "int32", "chromEnd": "int32",
            "score": "float32", "padj": str,
        },
    )

    padj_num   = pd.to_numeric(df["padj"], errors="coerce")
    sig_native = df[padj_num < sig_padj][KEY_COLS + ["score"]].copy()
    na_native  = df[padj_num.isna()][KEY_COLS + ["score"]].copy()

    print(
        f"[info] Fisher: {len(df):,} total | "
        f"significant (padj<{sig_padj}): {len(sig_native):,} | "
        f"IVT-absent: {len(na_native):,}",
        file=sys.stderr,
    )
    return sig_native, na_native


def load_ivt(path):
    """Return DataFrame with KEY_COLS + score for all IVT sites."""
    col_names, n_skip = detect_header(path)
    for c in KEY_COLS + ["score"]:
        if c not in col_names:
            raise ValueError(f"Column '{c}' not found in {path}.")

    df = pd.read_csv(
        path, sep="\t", skiprows=n_skip, header=None,
        names=col_names,
        usecols=KEY_COLS + ["score"],
        dtype={
            "chrom": "category", "strand": "category", "name": "category",
            "chromStart": "int32", "chromEnd": "int32",
            "score": "float32",
        },
    )
    print(f"[info] IVT: {len(df):,} sites loaded.", file=sys.stderr)
    return df


def get_ivt_sig_scores(sig_native, ivt_all):
    """Join significant native sites with IVT BEDRMod to retrieve IVT scores."""
    merged = sig_native[KEY_COLS].merge(
        ivt_all[KEY_COLS + ["score"]].rename(columns={"score": "ivt_score"}),
        on=KEY_COLS,
        how="inner",
    )
    print(
        f"[info] Significant sites matched in IVT: {len(merged):,} "
        f"(of {len(sig_native):,})",
        file=sys.stderr,
    )
    return merged


def draw_density(ax, groups, title, sig_padj):
    any_drawn = False
    for label, color, arr in groups:
        arr = arr[np.isfinite(arr)]
        if len(arr) < 5:
            continue
        kde = gaussian_kde(arr, bw_method="scott")
        xs  = np.linspace(0, 1, 500)
        ys  = kde(xs)
        ax.plot(xs, ys, color=color, linewidth=2.0,
                label=f"{label}  (n={len(arr):,})")
        ax.fill_between(xs, ys, alpha=0.1, color=color)
        any_drawn = True

    if not any_drawn:
        ax.set_title(f"{title}  (no data)")
        return

    ax.set_xlim(0, 1)
    ax.set_xlabel("score", fontsize=14)
    ax.set_ylabel("density", fontsize=14)
    ax.set_title(f"{title}\n(significant: padj < {sig_padj})", fontsize=14)
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=11, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_histogram(ax, groups, title, sig_padj):
    any_drawn = False
    for label, color, arr in groups:
        arr = arr[np.isfinite(arr)]
        arr = arr[(arr >= 0.7) & (arr <= 1.0)]
        if len(arr) == 0:
            continue
        ax.hist(arr, bins=30, range=(0.7, 1.0), color=color,
                histtype="stepfilled", alpha=0.5,
                label=f"{label}  (n={len(arr):,})")
        any_drawn = True

    if not any_drawn:
        ax.set_title(f"{title}  (no data)")
        return

    ax.set_yscale("log")
    ax.set_xlim(0.7, 1.0)
    ax.set_xlabel("score", fontsize=14)
    ax.set_ylabel("count (log scale)", fontsize=14)
    ax.set_title(f"{title}\n(score ≥ 0.7)", fontsize=14)
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=11, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def build_groups(sig_nat, sig_ivt, ivt_all, na_nat, mod=None):
    def scores(df, col="score"):
        if mod is not None:
            df = df[df["name"].astype(str) == mod]
        return df[col].values.astype(float)

    return [
        ("native, significant", GROUPS_STYLE["native, significant"], scores(sig_nat)),
        ("IVT, significant",    GROUPS_STYLE["IVT, significant"],    scores(sig_ivt, "ivt_score")),
        ("IVT, all",            GROUPS_STYLE["IVT, all"],            scores(ivt_all)),
        ("native, IVT-absent",  GROUPS_STYLE["native, IVT-absent"],  scores(na_nat)),
    ]


_SUMMARY_BINS = np.linspace(0.7, 1.0, 31)   # 30 bins, score > 0.7


def _hist_steps(arr, bins):
    """Return (x, y) suitable for step/fill_between; closes at bins[-1] with 0."""
    counts, edges = np.histogram(arr, bins=bins)
    x = np.append(edges[:-1], edges[-1])
    y = np.append(counts, 0)
    return x, y


def draw_summary_page(sig_nat, ivt_all, color_map, sig_padj, stem):
    """
    Returns a Figure with one row per modification (score > 0.7).
    Native (padj < sig_padj): filled step histogram.
    IVT all: dashed step outline.
    Works even when a mod has very few sites.
    """
    mods = sorted(
        set(sig_nat["name"].astype(str).unique())
        | set(ivt_all["name"].astype(str).unique())
    )
    n = len(mods)

    fig, axes = plt.subplots(n, 1, figsize=(7, 1.2 * n), sharex=True)
    if n == 1:
        axes = [axes]
    fig.suptitle(f"{stem}  ·  score distribution", fontsize=16)

    for ax, mod in zip(axes, mods):
        color = color_map.get(mod, "#aaaaaa")

        arr_nat = sig_nat[sig_nat["name"].astype(str) == mod]["score"].values.astype(float)
        arr_nat = arr_nat[np.isfinite(arr_nat) & (arr_nat > 0.7)]

        arr_ivt = ivt_all[ivt_all["name"].astype(str) == mod]["score"].values.astype(float)
        arr_ivt = arr_ivt[np.isfinite(arr_ivt) & (arr_ivt > 0.7)]

        vmax = 1.0
        if len(arr_nat) > 0:
            x, y = _hist_steps(arr_nat, _SUMMARY_BINS)
            vmax = max(vmax, float(y.max()))
            ax.fill_between(x, y, step="post", color=color, alpha=0.35,
                            linewidth=0)
            ax.step(x, y, where="post", color=color, linewidth=1.5)

        if len(arr_ivt) > 0:
            x, y = _hist_steps(arr_ivt, _SUMMARY_BINS)
            vmax = max(vmax, float(y.max()))
            ax.step(x, y, where="post", color=color, linewidth=1.5,
                    linestyle="--")

        ax.axvline(0.9, color="#888888", linewidth=1.0, linestyle=":",
                   zorder=0)
        ax.set_title(mod, fontsize=13, color=color, loc="left", pad=2)
        ax.set_xlim(0.7, 1.0)
        ax.set_yscale("log")
        # Powers of 10 already give a readable scale once the counts span a
        # decade; only sub-10 panels need finer 1-2-5 ticks to avoid showing a
        # lone "1".  Headroom guarantees at least two ticks are visible.
        ax.set_ylim(0.7, max(vmax * 1.8, 3.0))
        subs = (1.0,) if vmax >= 10 else (1.0, 2.0, 5.0)
        ax.yaxis.set_major_locator(plt.LogLocator(base=10, subs=subs))
        ax.yaxis.set_minor_locator(plt.NullLocator())
        ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda v, _: f"{int(round(v)):,}" if v >= 1 else "")
        )
        ax.tick_params(labelsize=11)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel("score", fontsize=14)
    fig.supylabel("site count", fontsize=14)

    style_handles = [
        Patch(facecolor="black", alpha=0.35, edgecolor="black", linewidth=1.5,
              label=f"native  (padj < {sig_padj})"),
        Line2D([0], [0], color="black", linewidth=1.5, linestyle="--",
               label="IVT, all"),
    ]
    axes[0].legend(handles=style_handles,
                   bbox_to_anchor=(1.02, 1.0), loc="upper left",
                   frameon=False, fontsize=12)

    fig.tight_layout()
    return fig


def write_pdf(sig_nat, sig_ivt, ivt_all, na_nat, stem, sig_padj, color_map, path):
    mods = sorted(
        set(sig_nat["name"].astype(str).unique())
        | set(na_nat["name"].astype(str).unique())
        | set(ivt_all["name"].astype(str).unique())
    )
    n_pages = 2 + len(mods)

    with PdfPages(path) as pdf:
        # ── Page 1: summary — vertically faceted by modification, score > 0.7
        fig = draw_summary_page(sig_nat, ivt_all, color_map, sig_padj, stem)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig, (ax_kde, ax_hist) = plt.subplots(1, 2, figsize=(18, 6))
        fig.suptitle(f"{stem}  ·  score distributions  (all modifications)",
                     fontsize=14)
        grps = build_groups(sig_nat, sig_ivt, ivt_all, na_nat)
        draw_density(ax_kde, grps, "all modifications", sig_padj)
        draw_histogram(ax_hist, grps, "all modifications", sig_padj)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for mod in mods:
            fig, (ax_kde, ax_hist) = plt.subplots(1, 2, figsize=(18, 6))
            fig.suptitle(f"{stem}  ·  score distributions  ·  {mod}", fontsize=14)
            grps = build_groups(sig_nat, sig_ivt, ivt_all, na_nat, mod=mod)
            draw_density(ax_kde, grps, mod, sig_padj)
            draw_histogram(ax_hist, grps, mod, sig_padj)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"[info] PDF written ({n_pages} pages): {path}", file=sys.stderr)


def main():
    args      = parse_args()
    stem      = args.fisher.rsplit("/", 1)[-1]
    color_map = dict(args.colors)

    sig_nat, na_nat = load_fisher(args.fisher, args.sig_padj)
    ivt_all         = load_ivt(args.ivt)
    sig_ivt         = get_ivt_sig_scores(sig_nat, ivt_all)

    write_pdf(sig_nat, sig_ivt, ivt_all, na_nat, stem, args.sig_padj,
              color_map, args.output_pdf)
    print("[info] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
