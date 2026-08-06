#!/usr/bin/env python3
"""
Per-biotype scatter plots for filtered modification sites.

Four-page PDF; pages 1-3 are grids with one panel per modification, each in
its own colour:
  Page 1 — Native frequency vs IVT frequency (testable sites; log-log)
  Page 2 — Score vs native frequency (all sites)
  Page 3 — Score vs log₂(native / IVT freq) (testable sites)
  Page 4 — Native frequency distribution per modification (violin / strip)
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
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.backends.backend_pdf import PdfPages


plt.rcParams.update({
    "font.size":        17,
    "axes.titlesize":   20,
    "axes.labelsize":   18,
    "xtick.labelsize":  16,
    "ytick.labelsize":  16,
    "legend.fontsize":  16,
    "figure.dpi":       150,
})

MAX_SCATTER = 20_000
MAX_VIOLIN  = 50_000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",  required=True,
                   help="Filtered BEDRMod file (plain or .gz)")
    p.add_argument("--name",   required=True,
                   help="Biotype label for plot titles")
    p.add_argument(
        "--color", nargs=2, metavar=("MOD", "HEX"),
        action="append", default=[],
        help="Color for a modification (repeatable)",
    )
    p.add_argument("--score-threshold", type=float, default=None, metavar="S",
                   help="Show a vertical reference line at this score value")
    p.add_argument("--log2fc-threshold", type=float, default=None, metavar="L",
                   help="Show a horizontal reference line at this log2FC value")
    p.add_argument("--output-pdf", required=True, help="Output PDF path")
    return p.parse_args()


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


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
        raise ValueError(f"No column header found in {path}")
    return col_names, n_skip


def load_data(path):
    col_names, n_skip = detect_header(path)
    needed = ["name", "frequency", "ivt_frequency", "padj", "score"]
    missing = [c for c in needed if c not in col_names]
    if missing:
        raise ValueError(f"Columns missing from {path}: {missing}")
    df = pd.read_csv(
        path, sep="\t", skiprows=n_skip, header=None,
        names=col_names, usecols=needed, dtype=str,
    )
    for col in ("frequency", "ivt_frequency", "padj", "score"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["is_testable"] = df["padj"].notna()
    with np.errstate(divide="ignore", invalid="ignore"):
        df["log2fc"] = np.where(
            df["is_testable"] & (df["frequency"] > 0) & (df["ivt_frequency"] > 0),
            np.log2(df["frequency"].values / df["ivt_frequency"].values),
            np.nan,
        )
    return df


def canonical_order(color_map, all_mods):
    ordered = [m for m in color_map if m in all_mods]
    ordered += sorted(all_mods - set(ordered))
    return ordered


def _clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _subsample(df, n, rng):
    return df.sample(n, random_state=rng) if len(df) > n else df


# ── Per-modification grid helpers ─────────────────────────────────────────────

_MAX_COLS = 4


def _make_grid(n, sharex=False, sharey=False):
    n_cols = min(n, _MAX_COLS)
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(4.3 * n_cols, 4.8 * n_rows),
        squeeze=False, sharex=sharex, sharey=sharey,
    )
    flat = [axes[r][c] for r in range(n_rows) for c in range(n_cols)]
    return fig, flat


def _finish_grid(axes, n_used, xlabel, ylabel):
    n_cols = min(n_used, _MAX_COLS)
    for i, ax in enumerate(axes[:n_used]):
        ax.set_box_aspect(1)                       # square plotting area
        _clean(ax)
        if i % n_cols == 0:                        # left column
            ax.set_ylabel(ylabel, fontsize=13, labelpad=5)
        if i + n_cols >= n_used:                   # bottom edge (incl. ragged row)
            ax.set_xlabel(xlabel, fontsize=13, labelpad=5)
    for ax in axes[n_used:]:
        ax.set_visible(False)


# ── Page 1: native freq vs IVT freq (log-log), one panel per mod ──────────────

def page_freq_scatter(df, color_map, present, name):
    sel = df[
        df["is_testable"]
        & df["frequency"].notna() & df["ivt_frequency"].notna()
        & (df["frequency"] > 0) & (df["ivt_frequency"] > 0)
    ]
    if not sel.empty:
        # Independent per-axis ranges so each axis starts near its own data
        # minimum.  Native frequency (y) is floored by the filter, so extra
        # bottom margin is added to the y-axis to lengthen the y = x diagonal.
        xlo = max(sel["ivt_frequency"].min() * 0.7, 1e-3)
        xhi = sel["ivt_frequency"].max() * 1.4
        ylo = max(sel["frequency"].min() * 0.3, 1e-3)
        yhi = sel["frequency"].max() * 1.4
    else:
        xlo, xhi, ylo, yhi = 1e-3, 1.0, 1e-3, 1.0
    d_lo, d_hi = max(xlo, ylo), min(xhi, yhi)   # y = x over the overlapping range

    here = [m for m in present if (sel["name"] == m).any()]
    fig, axes = _make_grid(len(here))
    rng = np.random.default_rng(42)
    for ax, mod in zip(axes, here):
        sub = _subsample(sel[sel["name"] == mod], MAX_SCATTER, rng)
        ax.scatter(sub["ivt_frequency"], sub["frequency"],
                   color=color_map.get(mod, "#888888"),
                   s=16, alpha=0.35, linewidths=0)
        if d_lo < d_hi:
            ax.plot([d_lo, d_hi], [d_lo, d_hi], color="#aaaaaa",
                    linewidth=1.0, linestyle="--", zorder=0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
        ax.set_title(mod, color=color_map.get(mod, "#333333"),
                     fontsize=15, pad=8)

    _finish_grid(axes, len(here),
                 "IVT frequency (%)", "native frequency (%)")
    n_absent = int((~df["is_testable"]).sum())
    fig.suptitle(f"{name}  —  native vs IVT frequency  (testable sites)",
                 fontsize=19, y=0.995)
    fig.text(0.5, 0.955, f"IVT-absent: {n_absent:,} sites not shown",
             ha="center", va="top", fontsize=12, color="#777777")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


# ── Page 2: score vs native frequency, one panel per mod ──────────────────────

def page_score_vs_freq(df, color_map, present, name, score_thr):
    here = [m for m in present if (df["name"] == m).any()]
    fig, axes = _make_grid(len(here), sharex=True, sharey=True)
    rng = np.random.default_rng(42)
    for ax, mod in zip(axes, here):
        sub = _subsample(df[df["name"] == mod], MAX_SCATTER, rng)
        col = color_map.get(mod, "#888888")
        t = sub[sub["is_testable"]]
        a = sub[~sub["is_testable"]]
        if not t.empty:
            ax.scatter(t["score"], t["frequency"], color=col, s=16,
                       alpha=0.35, linewidths=0, marker="o")
        if not a.empty:
            ax.scatter(a["score"], a["frequency"], color=col, s=22,
                       alpha=0.5, linewidths=0.6, edgecolors="white", marker="^")
        if score_thr is not None:
            ax.axvline(score_thr, color="#555555", linewidth=1.1,
                       linestyle=":", zorder=0)
        ax.set_xlim(0, 1)
        ax.set_title(mod, color=col, fontsize=15, pad=8)

    _finish_grid(axes, len(here), "score", "native frequency (%)")
    fig.suptitle(f"{name}  —  score vs native frequency", fontsize=19, y=0.995)
    ht = mlines.Line2D([], [], color="#555555", marker="o", linestyle="None",
                       markersize=9, label="IVT-testable")
    ha = mlines.Line2D([], [], color="#555555", marker="^", linestyle="None",
                       markersize=9, label="IVT-absent")
    fig.legend(handles=[ht, ha], loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    return fig


# ── Page 3: score vs log2FC, one panel per mod ────────────────────────────────

def page_score_vs_log2fc(df, color_map, present, name, score_thr, log2fc_thr):
    sel = df[df["is_testable"] & df["log2fc"].notna()]
    here = [m for m in present if (sel["name"] == m).any()]
    fig, axes = _make_grid(len(here), sharex=True, sharey=True)
    rng = np.random.default_rng(42)
    for ax, mod in zip(axes, here):
        sub = _subsample(sel[sel["name"] == mod], MAX_SCATTER, rng)
        col = color_map.get(mod, "#888888")
        ax.scatter(sub["score"], sub["log2fc"], color=col, s=16,
                   alpha=0.35, linewidths=0)
        if score_thr is not None:
            ax.axvline(score_thr, color="#555555", linewidth=1.1,
                       linestyle=":", zorder=0)
        if log2fc_thr is not None:
            ax.axhline(log2fc_thr, color="#555555", linewidth=1.1,
                       linestyle="--", zorder=0)
        ax.set_xlim(0, 1)
        ax.set_title(mod, color=col, fontsize=15, pad=8)

    _finish_grid(axes, len(here), "score", "log₂(native / IVT freq)")
    fig.suptitle(f"{name}  —  score vs log₂FC", fontsize=19, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


# ── Page 4: native + IVT frequency distributions ─────────────────────────────

_W       = 0.20   # violin half-width (must be < _O)
_O       = 0.26   # offset from group centre for each pair
_SPACING = 1.3    # distance between consecutive modification groups


def _draw_violin_or_strip(ax, arr, pos, color, alpha, rng):
    if len(arr) >= 5 and np.std(arr) > 1e-9:
        parts = ax.violinplot(
            [arr], positions=[pos],
            showmedians=True, showextrema=False, widths=_W * 2,
        )
        parts["bodies"][0].set_facecolor(color)
        parts["bodies"][0].set_edgecolor("none")
        parts["bodies"][0].set_alpha(alpha)
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(1.8)
    elif len(arr) > 0:
        jitter = rng.uniform(-_O * 0.6, _O * 0.6, len(arr))
        ax.scatter(pos + jitter, arr, color=color,
                   s=16, alpha=alpha, linewidths=0)


def plot_freq_distributions(ax, df, color_map, mods_ordered, title):
    present = [m for m in mods_ordered if m in df["name"].astype(str).values]
    rng = np.random.default_rng(42)

    for idx, mod in enumerate(present):
        mod_df = df[df["name"].astype(str) == mod]
        col    = color_map.get(mod, "#888888")
        xc     = idx * _SPACING

        nat_arr = mod_df["frequency"].dropna().values
        if len(nat_arr) > MAX_VIOLIN:
            nat_arr = rng.choice(nat_arr, MAX_VIOLIN, replace=False)

        ivt_arr = mod_df.loc[
            mod_df["is_testable"] & mod_df["ivt_frequency"].notna(),
            "ivt_frequency",
        ].values
        if len(ivt_arr) > MAX_VIOLIN:
            ivt_arr = rng.choice(ivt_arr, MAX_VIOLIN, replace=False)

        _draw_violin_or_strip(ax, nat_arr, xc - _O, col, 0.82, rng)
        _draw_violin_or_strip(ax, ivt_arr, xc + _O, col, 0.38, rng)

    centres = [i * _SPACING for i in range(len(present))]
    ax.set_xticks(centres)
    ax.set_xticklabels(present, rotation=45, ha="right")
    if centres:
        ax.set_xlim(centres[0] - _SPACING * 0.6, centres[-1] + _SPACING * 0.6)
    ax.set_ylabel("frequency (%)", labelpad=8)
    ax.set_title(title, pad=12)
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f"{v:.0f}")
    )
    _clean(ax)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    color_map = {mod: hex_ for mod, hex_ in args.color}

    print(f"[info] Loading {args.input} ...", file=sys.stderr)
    df = load_data(args.input)
    n_test   = int(df["is_testable"].sum())
    n_absent = int((~df["is_testable"]).sum())
    all_mods     = set(df["name"].astype(str).unique())
    mods_ordered = canonical_order(color_map, all_mods)
    present_mods = [m for m in mods_ordered if m in all_mods]
    print(
        f"  {len(df):,} sites | testable={n_test:,} IVT-absent={n_absent:,} "
        f"| mods: {', '.join(present_mods)}",
        file=sys.stderr,
    )

    name = args.name

    with PdfPages(args.output_pdf) as pdf:

        # ── Pages 1-3: per-modification scatter grids ─────────────────────
        for fig in (
            page_freq_scatter(df, color_map, present_mods, name),
            page_score_vs_freq(
                df, color_map, present_mods, name, args.score_threshold),
            page_score_vs_log2fc(
                df, color_map, present_mods, name,
                args.score_threshold, args.log2fc_threshold),
        ):
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # ── Page 4: native + IVT freq distributions ───────────────────────
        fig, ax = plt.subplots(
            figsize=(max(4.5, (_SPACING * len(present_mods) * 0.72 + 1.5) * 0.9), 5)
        )
        plot_freq_distributions(
            ax, df, color_map, mods_ordered,
            f"{name}  —  native vs IVT frequency per modification",
        )
        h_nat = mpatches.Patch(color="#555555", alpha=0.82, label="native")
        h_ivt = mpatches.Patch(color="#555555", alpha=0.38, label="IVT")
        fig.legend(
            handles=[h_nat, h_ivt], loc="lower center",
            ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0),
        )
        fig.tight_layout(rect=[0, 0.06, 1, 1])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"[info] Written → {args.output_pdf}", file=sys.stderr)


if __name__ == "__main__":
    main()
