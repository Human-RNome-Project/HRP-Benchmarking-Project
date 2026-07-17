#!/usr/bin/env python3
"""
Panel A: Manhattan-like plot of polyA modification density across 1 Mb bins.

This is a GitHub-ready adaptation of:
  draft_reference/new_run/plot_polyA_manhattan_density.py

Inputs:
  - outputs/tiered_lists/tiered_polyA.tsv
  - GRCh38.primary_assembly.genome.fa.fai (1 Mb bin definitions)

Outputs:
  - figures/panel_a_polyA_manhattan_density/panel_a_polyA_manhattan_density.pdf
  - figures/panel_a_polyA_manhattan_density/panel_a_polyA_manhattan_density.png
  - figures/panel_a_polyA_manhattan_density/panel_a_polyA_1Mb_bin_counts.tsv
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
GITHUB_ROOT = Path(__file__).resolve().parents[2]
TIERED_POLYA = GITHUB_ROOT / "outputs" / "tiered_lists" / "tiered_polyA.tsv"
OUTDIR = GITHUB_ROOT / "figures" / "panel_a_polyA_manhattan_density"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Locate the GRCh38 .fai file: prefer a local copy, fall back to the reference area
LOCAL_FAI = GITHUB_ROOT / "inputs" / "GRCh38.primary_assembly.genome.fa.fai"
FALLBACK_FAI = Path.home() / "ref" / "GRCh38.primary_assembly.genome.fa.fai"
FAI = LOCAL_FAI if LOCAL_FAI.exists() else FALLBACK_FAI
if not FAI.exists():
    raise FileNotFoundError(
        f"GRCh38 .fai not found at {LOCAL_FAI} or {FALLBACK_FAI}. "
        "Place the index in inputs/ or update the path in this script."
    )

# ── Parameters ───────────────────────────────────────────────────────────────
BIN_SIZE = 1_000_000  # 1 Mb
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
BAR_COLOR = "#0D3B6E"

# ── Load chromosome lengths ──────────────────────────────────────────────────
CHR_LENGTHS = {}
with open(FAI) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            CHR_LENGTHS[p[0]] = int(p[1])

# Keep only requested chromosomes that are present in the .fai
CHR_ORDER = [c for c in CHR_ORDER if c in CHR_LENGTHS]

# ── Build 1 Mb bin counts from tiered polyA sites ────────────────────────────
tiered = pd.read_csv(TIERED_POLYA, sep="\t")
# Restrict to primary chromosomes and count sites per 1 Mb bin
tiered = tiered[tiered["chr"].isin(CHR_ORDER)].copy()
tiered["bin_index"] = tiered["start"] // BIN_SIZE

observed_counts = (
    tiered.groupby(["chr", "bin_index"])
    .size()
    .reset_index(name="n_mods")
)

# Build a complete set of 1 Mb bins from the .fai (including zero-count bins)
full_bins = []
for chrom in CHR_ORDER:
    chrom_len = CHR_LENGTHS[chrom]
    n_bins = int(np.ceil(chrom_len / BIN_SIZE))
    for b in range(n_bins):
        bin_start = b * BIN_SIZE
        bin_end = min((b + 1) * BIN_SIZE, chrom_len)
        full_bins.append({
            "chr": chrom,
            "bin_index": b,
            "bin_start": bin_start,
            "bin_end": bin_end,
        })

bin_counts = pd.DataFrame(full_bins)
bin_counts = bin_counts.merge(observed_counts, on=["chr", "bin_index"], how="left")
bin_counts["n_mods"] = bin_counts["n_mods"].fillna(0).astype(int)

# Save the generated bin counts alongside the figure
bin_counts.to_csv(OUTDIR / "panel_a_polyA_1Mb_bin_counts.tsv", sep="\t", index=False)

# ── Compute cumulative genomic positions for Manhattan layout ──────────────────
cumulative_offset = {}
offset = 0
for chrom in CHR_ORDER:
    cumulative_offset[chrom] = offset
    offset += CHR_LENGTHS[chrom]

total_length = offset

bin_counts["x_mid"] = (
    bin_counts["bin_index"] * BIN_SIZE
    + (bin_counts["bin_end"] - bin_counts["bin_start"]) / 2
    + bin_counts["chr"].map(cumulative_offset)
)

# ── Figure setup ─────────────────────────────────────────────────────────────
fig_width_mm = 380
fig_height_mm = 100
fig_width_in = fig_width_mm / 25.4
fig_height_in = fig_height_mm / 25.4

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 8
plt.rcParams["axes.labelsize"] = 16
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 8
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in))

# ── Alternating chromosome background shades ─────────────────────────────────
for i, chrom in enumerate(CHR_ORDER):
    start = cumulative_offset[chrom]
    end = start + CHR_LENGTHS[chrom]
    shade = "#F0F0F0" if i % 2 == 0 else "#FFFFFF"
    ax.axvspan(start, end, facecolor=shade, edgecolor="none", alpha=1.0, zorder=0)

# ── Plot bars ────────────────────────────────────────────────────────────────
for chrom in CHR_ORDER:
    sub = bin_counts[bin_counts["chr"] == chrom]
    ax.bar(
        sub["x_mid"],
        sub["n_mods"],
        width=BIN_SIZE,
        color=BAR_COLOR,
        alpha=1.0,
        edgecolor="none",
        zorder=2,
    )

# ── Chromosome boundary lines and labels ─────────────────────────────────────
prev_end = 0
xticks = []
xtick_labels = []
for chrom in CHR_ORDER:
    start = cumulative_offset[chrom]
    end = start + CHR_LENGTHS[chrom]
    mid = (start + end) / 2
    xticks.append(mid)
    xtick_labels.append(chrom.replace("chr", ""))
    if prev_end > 0:
        ax.axvline(x=prev_end, color="#CCCCCC", linewidth=0.4, zorder=0)
    prev_end = end

ax.set_xticks(xticks)
ax.set_xticklabels(xtick_labels, fontweight="bold")
plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
for label in ax.get_yticklabels():
    label.set_fontweight("bold")

# ── Axis labels and limits ───────────────────────────────────────────────────
ax.set_xlim(-BIN_SIZE / 2, total_length + BIN_SIZE / 2)
ax.set_ylim(0, 100)
ax.set_xlabel("Chromosome", fontweight="bold")
ax.set_ylabel("Modifications per 1 Mb bin", fontweight="bold", fontsize=12)

# ── Annotate top bins on specific chromosomes ──────────────────────────────────
target_chroms = ["chr1", "chr4", "chr19"]
for t_chrom in target_chroms:
    sub = bin_counts[bin_counts["chr"] == t_chrom]
    if sub.empty:
        continue
    top_bin = sub.loc[sub["n_mods"].idxmax()]
    bin_start = int(top_bin["bin_index"] * BIN_SIZE / 1e6)
    bin_end = int((top_bin["bin_index"] + 1) * BIN_SIZE / 1e6)
    label_text = f"{t_chrom.replace('chr', '')}:{bin_start}-{bin_end} Mb"
    y_pos = min(top_bin["n_mods"], 95)

    ax.annotate(
        label_text,
        xy=(top_bin["x_mid"], top_bin["n_mods"]),
        xytext=(0, 15),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color="black",
        arrowprops=dict(arrowstyle="-|>", color="black", lw=1.0),
    )

# ── Clean spines ─────────────────────────────────────────────────────────────
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.8)
ax.spines["bottom"].set_linewidth(0.8)
ax.tick_params(axis="both", which="major", width=0.8)

plt.tight_layout()
plt.savefig(OUTDIR / "panel_a_polyA_manhattan_density.pdf", format="pdf", bbox_inches="tight")
plt.savefig(OUTDIR / "panel_a_polyA_manhattan_density.png", format="png", bbox_inches="tight", dpi=300)
plt.close()

print(f"Saved panel A plots to {OUTDIR}")
