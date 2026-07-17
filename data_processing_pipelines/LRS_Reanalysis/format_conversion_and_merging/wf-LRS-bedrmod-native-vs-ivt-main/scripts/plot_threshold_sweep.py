#!/usr/bin/env python3
"""
Threshold sweep: for each padj (testable sites) and score (IVT-absent sites)
threshold, plot the number of surviving sites per modification.

One figure with one row per dataset × 2 columns (padj sweep | score sweep).
Each subplot has one curve per modification plus a dashed total; y-axis is log.
A vertical reference line marks the chosen operating thresholds.
"""

import ast
import gzip
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PADJ_THRESHOLDS  = np.linspace(0.001, 0.50, 400)
SCORE_THRESHOLDS = np.linspace(0.50,  0.99, 400)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset", nargs=2, metavar=("NAME", "PATH"),
                   action="append", default=[],
                   help="Dataset name and path to native_vs_ivt_fisher.bed "
                        "(repeat for multiple datasets)")
    p.add_argument("--output-pdf",  required=True)
    return p.parse_args()


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def detect_header(path):
    col_names = None
    with _open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            stripped = line.lstrip("#").strip()
            if stripped.startswith("["):
                try:
                    col_names = ast.literal_eval(stripped)
                except (ValueError, SyntaxError):
                    pass
    if col_names is None:
        raise ValueError(f"Could not detect header from {path}")
    return col_names


def load_data(path):
    """Return dict: mod -> {'padj': float32 array, 'score': float32 array}."""
    col_names = detect_header(path)
    name_idx  = col_names.index("name")
    score_idx = col_names.index("score")
    padj_idx  = col_names.index("padj")

    padj_lists  = {}
    score_lists = {}

    print(f"[info] Reading {path} ...", file=sys.stderr)
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= padj_idx or fields[padj_idx] == "padj":
                continue

            mod      = fields[name_idx]
            padj_val = fields[padj_idx]

            if padj_val == "NA":
                try:
                    score_lists.setdefault(mod, []).append(float(fields[score_idx]))
                except ValueError:
                    pass
            else:
                try:
                    padj_lists.setdefault(mod, []).append(float(padj_val))
                except ValueError:
                    pass

    mods = sorted(set(padj_lists) | set(score_lists))
    result = {}
    for mod in mods:
        result[mod] = {
            "padj":  np.array(padj_lists.get(mod,  []), dtype=np.float32),
            "score": np.array(score_lists.get(mod, []), dtype=np.float32),
        }
        print(
            f"  {mod}: {len(result[mod]['padj']):,} testable  "
            f"{len(result[mod]['score']):,} IVT-absent",
            file=sys.stderr,
        )
    return result


def sweep(arr, thresholds, less_than):
    arr2 = arr[:, None]
    t    = thresholds[None, :]
    return (arr2 < t).sum(axis=0) if less_than else (arr2 > t).sum(axis=0)


def draw_panel(ax, data, thresholds, key, less_than, title, colors):
    total = np.zeros(len(thresholds), dtype=np.int64)
    for mod, color in zip(sorted(data), colors):
        arr = data[mod][key]
        if len(arr) == 0:
            continue
        counts = sweep(arr, thresholds, less_than)
        total += counts
        ax.plot(thresholds, np.maximum(counts, 1), label=mod,
                color=color, linewidth=1.5)

    ax.plot(thresholds, np.maximum(total, 1), color="black", linewidth=2,
            linestyle="--", label="Total")

    ax.set_yscale("log")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Threshold", fontsize=8)
    ax.set_ylabel("Sites surviving", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}")
    )


def main():
    args = parse_args()
    if not args.dataset:
        sys.exit("Error: provide at least one --dataset NAME PATH pair")

    all_mods = set()
    datasets = {}
    for name, path in args.dataset:
        d = load_data(path)
        datasets[name] = d
        all_mods |= set(d)

    mods_sorted = sorted(all_mods)
    cmap   = plt.get_cmap("tab10")
    colors = {mod: cmap(i % 10) for i, mod in enumerate(mods_sorted)}

    n_rows = len(datasets)
    fig, axes = plt.subplots(n_rows, 2, figsize=(13, 4.5 * n_rows))
    if n_rows == 1:
        axes = axes[np.newaxis, :]
    fig.suptitle("Sites surviving per threshold — by modification  (log scale)",
                 fontsize=11)

    for row, (name, data) in enumerate(datasets.items()):
        mod_colors = [colors[m] for m in sorted(data)]

        draw_panel(
            axes[row, 0], data, PADJ_THRESHOLDS, "padj", True,
            f"{name}  ·  testable sites  (padj < threshold)",
            mod_colors,
        )
        draw_panel(
            axes[row, 1], data, SCORE_THRESHOLDS, "score", False,
            f"{name}  ·  IVT-absent sites  (score > threshold)",
            mod_colors,
        )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6,
               fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.04, 1, 1])

    with PdfPages(args.output_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    print(f"[info] Written: {args.output_pdf}", file=sys.stderr)


if __name__ == "__main__":
    main()
