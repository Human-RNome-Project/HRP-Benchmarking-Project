#!/usr/bin/env python3
"""
Publication-quality summary of filtered native vs IVT sites.

Three-page PDF:
  Page 1 — Stacked bar: site counts per modification and biotype, split
            into testable (padj-passing) and IVT-absent panels.
  Page 2 — Violin: log2(native_freq / IVT_freq) per modification, faceted
            by biotype (testable sites only).
  Page 3 — Heatmap: modifications x biotypes, cell colour = log10(count),
            annotated with raw counts.
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
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpecFromSubplotSpec


MAX_VIOLIN_POINTS = 50_000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--dataset", nargs=2, metavar=("NAME", "PATH"),
        action="append", default=[],
        help="Dataset name and filtered BED/BED.gz path (repeatable)",
    )
    p.add_argument(
        "--color", nargs=2, metavar=("MOD", "HEX"),
        action="append", default=[],
        help="Colour for a modification (repeatable)",
    )
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


def load_filtered(path):
    col_names, n_skip = detect_header(path)
    needed = ["name", "frequency", "ivt_frequency", "padj", "score"]
    missing = [c for c in needed if c not in col_names]
    if missing:
        raise ValueError(f"Columns missing from {path}: {missing}")
    df = pd.read_csv(
        path, sep="\t", skiprows=n_skip, header=None,
        names=col_names, usecols=needed, dtype=str,
    )
    df["frequency"]     = pd.to_numeric(df["frequency"],     errors="coerce")
    df["ivt_frequency"] = pd.to_numeric(df["ivt_frequency"], errors="coerce")
    df["padj"]          = pd.to_numeric(df["padj"],          errors="coerce")
    df["score"]         = pd.to_numeric(df["score"],         errors="coerce")
    return df


def canonical_order(color_map, all_mods):
    """Mods in color_map insertion order, then any extras alphabetically."""
    ordered = [m for m in color_map if m in all_mods]
    ordered += sorted(all_mods - set(ordered))
    return ordered


# ── Page 1 top: per-biotype mod bars ────────────────────────────────────────

def plot_mod_bars(axes, datasets, color_map, mods_ordered):
    """One bar per modification per biotype; count annotated above each bar."""
    for ax, (bio, df) in zip(axes, datasets.items()):
        present = [m for m in mods_ordered
                   if m in df["name"].astype(str).values]
        if not present:
            ax.set_title(bio, fontsize=14, fontweight="bold")
            continue

        x      = np.arange(len(present))
        counts = [int((df["name"].astype(str) == m).sum()) for m in present]
        colors = [color_map.get(m, "#888888") for m in present]

        bars = ax.bar(x, counts, color=colors, width=0.7)

        top = max(counts) if counts else 1
        for bar, count in zip(bars, counts):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + top * 0.02,
                f"{count:,}", ha="center", va="bottom", fontsize=11,
            )

        ax.set_xticks(x)
        ax.set_xticklabels(present, rotation=45, ha="right", fontsize=12)
        ax.set_ylabel("site count", fontsize=13)
        ax.set_title(bio, fontsize=14, fontweight="bold")
        ax.margins(y=0.15)
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _: f"{int(v):,}")
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)


# ── Page 1 bottom: stacked bar counts ───────────────────────────────────────

def plot_counts(fig, axes, datasets, color_map, mods_ordered):
    biotypes = list(datasets.keys())
    x = np.arange(len(biotypes))

    for ax, is_testable, title in [
        (axes[0], True,  "Testable sites  (padj ≤ threshold)"),
        (axes[1], False, "IVT-absent sites  (score > threshold)"),
    ]:
        bottoms = np.zeros(len(biotypes))
        for mod in mods_ordered:
            heights = []
            for bio in biotypes:
                df   = datasets[bio]
                mask = df["name"].astype(str) == mod
                mask &= df["padj"].notna() if is_testable else df["padj"].isna()
                heights.append(int(mask.sum()))
            heights = np.array(heights)
            if heights.sum() == 0:
                continue
            ax.bar(
                x, heights, bottom=bottoms,
                color=color_map.get(mod, "#888888"),
                label=mod, width=0.65,
            )
            bottoms += heights

        top = bottoms.max() if bottoms.max() > 0 else 1
        for xi, total in zip(x, bottoms):
            if total > 0:
                ax.text(xi, total + top * 0.02, f"{int(total):,}",
                        ha="center", va="bottom", fontsize=11)

        ax.set_xticks(x)
        ax.set_xticklabels(biotypes, fontsize=13)
        ax.set_ylabel("site count", fontsize=13)
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.margins(y=0.15)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _: f"{int(v):,}")
        )


# ── Page 2: log2FC violin ────────────────────────────────────────────────────

def _color_violins(parts, colors):
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("none")
        body.set_alpha(0.82)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.8)


def plot_log2fc(fig, axes, datasets, color_map, mods_ordered):
    biotypes = list(datasets.keys())
    rng = np.random.default_rng(seed=42)

    for ax, bio in zip(axes, biotypes):
        df = datasets[bio]
        testable = df[
            df["padj"].notna()
            & df["ivt_frequency"].notna()
            & (df["frequency"]     > 0)
            & (df["ivt_frequency"] > 0)
        ].copy()

        if testable.empty:
            ax.set_title(bio, fontsize=14, fontweight="bold")
            ax.text(0.5, 0.5, "no testable sites", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="grey")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            continue

        testable["log2fc"] = np.log2(
            testable["frequency"] / testable["ivt_frequency"]
        )

        present = [m for m in mods_ordered
                   if m in testable["name"].astype(str).values]
        if not present:
            continue

        data_per_mod, colors_list = [], []
        for mod in present:
            arr = testable.loc[
                testable["name"].astype(str) == mod, "log2fc"
            ].values
            if len(arr) > MAX_VIOLIN_POINTS:
                arr = rng.choice(arr, MAX_VIOLIN_POINTS, replace=False)
            data_per_mod.append(arr)
            colors_list.append(color_map.get(mod, "#888888"))

        positions = np.arange(len(present))

        violin_idx = [i for i, a in enumerate(data_per_mod)
                      if len(a) >= 5 and np.std(a) > 1e-9]
        strip_idx  = [i for i in range(len(present))
                      if i not in violin_idx]

        if violin_idx:
            parts = ax.violinplot(
                [data_per_mod[i] for i in violin_idx],
                positions=[positions[i] for i in violin_idx],
                showmedians=True, showextrema=False,
            )
            _color_violins(parts, [colors_list[i] for i in violin_idx])

        for i in strip_idx:
            ax.scatter(
                np.full(len(data_per_mod[i]), positions[i]),
                data_per_mod[i],
                color=colors_list[i], s=20, alpha=0.7,
            )

        ax.axhline(1.0, color="#aaaaaa", linestyle="--", linewidth=1)
        ax.set_xticks(positions)
        ax.set_xticklabels(present, rotation=45, ha="right", fontsize=12)
        ax.set_ylabel("log₂(native / IVT freq)", fontsize=13)
        ax.set_title(bio, fontsize=14, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)


# ── Page 3: heatmap ──────────────────────────────────────────────────────────

def plot_heatmap(ax, datasets, color_map, mods_ordered):
    biotypes = list(datasets.keys())
    present  = [m for m in mods_ordered
                if any(m in df["name"].astype(str).values
                       for df in datasets.values())]

    counts = np.zeros((len(present), len(biotypes)), dtype=int)
    for j, bio in enumerate(biotypes):
        for i, mod in enumerate(present):
            counts[i, j] = int(
                (datasets[bio]["name"].astype(str) == mod).sum()
            )

    log_counts = np.log10(counts.astype(float) + 1)
    vmax = max(log_counts.max(), 1.0)

    im = ax.imshow(log_counts, cmap="YlOrRd", aspect="auto",
                   vmin=0, vmax=vmax)

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.03)
    cbar.set_label("log₁₀(count + 1)", fontsize=12)

    # White cell separators
    ax.set_xticks(np.arange(-0.5, len(biotypes)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(present)),  minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)

    for i in range(len(present)):
        for j in range(len(biotypes)):
            val = counts[i, j]
            txt_color = "white" if log_counts[i, j] / vmax > 0.6 else "black"
            ax.text(j, i, f"{val:,}", ha="center", va="center",
                    fontsize=12, color=txt_color, fontweight="bold")

    ax.set_xticks(range(len(biotypes)))
    ax.set_yticks(range(len(present)))
    ax.set_xticklabels(biotypes, fontsize=13)
    ax.set_yticklabels(present, fontsize=12)
    ax.set_title("filtered sites per modification and biotype",
                 fontsize=14, fontweight="bold", pad=12)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    if not args.dataset:
        sys.exit("error: provide at least one --dataset NAME PATH")

    color_map = {mod: hex_ for mod, hex_ in args.color}

    print("[info] Loading datasets ...", file=sys.stderr)
    datasets = {}
    for name, path in args.dataset:
        print(f"  {name}: {path}", file=sys.stderr)
        df = load_filtered(path)
        n_test   = int(df["padj"].notna().sum())
        n_absent = int(df["padj"].isna().sum())
        mods     = sorted(df["name"].astype(str).unique())
        print(
            f"    {len(df):,} sites | testable={n_test:,} "
            f"IVT-absent={n_absent:,} | mods: {', '.join(mods)}",
            file=sys.stderr,
        )
        datasets[name] = df

    all_mods     = set().union(*(df["name"].astype(str).unique()
                                 for df in datasets.values()))
    mods_ordered = canonical_order(color_map, all_mods)
    present_mods = [m for m in mods_ordered if m in all_mods]
    biotypes     = list(datasets.keys())
    n_bio        = len(biotypes)

    mod_handles = [
        mpatches.Patch(color=color_map.get(m, "#888888"), label=m)
        for m in present_mods
    ]

    with PdfPages(args.output_pdf) as pdf:

        # ── Page 1: per-biotype mod bars (top) + stacked counts (bottom) ──
        fig = plt.figure(figsize=(max(5 * n_bio, 14), 12))
        fig.suptitle("Filtered site counts", fontsize=16, fontweight="bold")
        outer = fig.add_gridspec(2, 1, hspace=0.55, height_ratios=[0.6, 1])

        top_gs = GridSpecFromSubplotSpec(
            1, n_bio, subplot_spec=outer[0], wspace=0.4,
        )
        axes_top = [fig.add_subplot(top_gs[i]) for i in range(n_bio)]
        plot_mod_bars(axes_top, datasets, color_map, mods_ordered)

        bot_gs = GridSpecFromSubplotSpec(
            1, 2, subplot_spec=outer[1], wspace=0.35,
        )
        axes_bot = [fig.add_subplot(bot_gs[i]) for i in range(2)]
        plot_counts(fig, axes_bot, datasets, color_map, mods_ordered)

        fig.legend(
            handles=mod_handles, loc="lower center",
            ncol=min(len(mod_handles), 7), fontsize=12,
            frameon=False, bbox_to_anchor=(0.5, 0.0),
        )
        fig.tight_layout(rect=[0, 0.06, 1, 1])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ── Page 2: log2FC violin ─────────────────────────────────────────
        fig, axes = plt.subplots(
            1, n_bio,
            figsize=(max(5 * n_bio, 10), 6),
            sharey=True,
        )
        if n_bio == 1:
            axes = [axes]
        fig.suptitle(
            "log₂(native / IVT frequency)  —  testable sites",
            fontsize=16, fontweight="bold",
        )
        plot_log2fc(fig, axes, datasets, color_map, mods_ordered)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ── Page 3: heatmap ───────────────────────────────────────────────
        n_mods = len(present_mods)
        fig, ax = plt.subplots(
            figsize=(max(5, 2.8 * n_bio), max(4, 0.65 * n_mods + 2))
        )
        plot_heatmap(ax, datasets, color_map, mods_ordered)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"[info] Written → {args.output_pdf}", file=sys.stderr)


if __name__ == "__main__":
    main()
