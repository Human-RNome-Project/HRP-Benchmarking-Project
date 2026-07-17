#!/usr/bin/env python3
"""
Regional plot of rRNA modifications across human rRNA sequences.

For each rRNA sequence that has at least one modification, draws a panel showing
modification positions as colored vertical bars along the rRNA sequence bar.
Only sequences with modifications are plotted.

Inputs:
  - draft_reference/tiered_lists/tiered_rRNA_only.tsv
  - ~/ref/hs_rRNAs_NR_046235.fa.fai

# Outputs:
#   - figures/panel_g_rRNA_region_PTC/panel_g_rRNA_region_PTC.pdf
#   - figures/panel_g_rRNA_region_PTC/panel_g_rRNA_region_PTC.png
"""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
GITHUB_ROOT = Path(__file__).resolve().parents[2]
OUTDIR = GITHUB_ROOT / "figures" / "panel_g_rRNA_region_PTC"
OUTDIR.mkdir(parents=True, exist_ok=True)
RRNA_MODS = GITHUB_ROOT / "outputs" / "tiered_lists" / "tiered_rRNA_only.tsv"

LOCAL_FAI = GITHUB_ROOT / "inputs" / "hs_rRNAs_NR_046235.fa.fai"
FALLBACK_FAI = Path.home() / "ref" / "hs_rRNAs_NR_046235.fa.fai"
FAI = LOCAL_FAI if LOCAL_FAI.exists() else FALLBACK_FAI

if not FAI.exists():
    raise FileNotFoundError(f"Reference FAI not found at {LOCAL_FAI} or {FALLBACK_FAI}.")

# ── Plotting setup ───────────────────────────────────────────────────────────
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 14
plt.rcParams["axes.labelsize"] = 16
plt.rcParams["axes.titlesize"] = 18
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["text.color"] = "black"
plt.rcParams["axes.labelcolor"] = "black"
plt.rcParams["xtick.color"] = "black"
plt.rcParams["ytick.color"] = "black"

# ── Official project RNA-modification colors ─────────────────────────────────
MOD_COLORS = {
    "Am": "#D44F3E",
    "Cm": "#0D3B6E",
    "Gm": "#74B354",
    "m5C": "#001427",
    "m6A": "#721817",
    "Um": "#C47A02",
    "Y": "#F0A202",
}

# ── rRNA structural domain annotations ───────────────────────────────────────
# Human rRNA domain boundaries (approximate, based on secondary structure):
#   18S: standard four-domain model (5′, central, 3′ major, 3′ minor)
#   28S: six canonical domains mapped to human 28S structural regions
#        (Gilmore et al., 2017, PMC5737216; RiboVision human 28S model)
#   5.8S: approximate three-part division for visualization
RRNA_DOMAINS = {
    "hs_rRNA_18S": [
        ("5′ domain", 1, 564, "#C8E6C9"),
        ("Central domain", 565, 926, "#BBDEFB"),
        ("3′ major domain", 927, 1393, "#FFE0B2"),
        ("3′ minor domain", 1394, 1869, "#F8BBD0"),
    ],
    "hs_rRNA_28S": [
        ("Domain I", 1, 1107, "#C8E6C9"),
        ("Domain II", 1108, 1686, "#BBDEFB"),
        ("Domain III", 1687, 2805, "#FFE0B2"),
        ("Domain IV", 2806, 3899, "#F8BBD0"),
        ("Domain V", 3900, 4597, "#E1BEE7"),
        ("Domain VI", 4598, 5070, "#B2DFDB"),
    ],
    "hs_rRNA_5.8S": [
        ("5′ region", 1, 60, "#C8E6C9"),
        ("Central region", 61, 110, "#BBDEFB"),
        ("3′ region", 111, 157, "#FFE0B2"),
    ],
}

# ── Load rRNA sequence lengths ───────────────────────────────────────────────
seq_lengths = {}
with open(FAI) as fh:
    for line in fh:
        p = line.strip().split("\t")
        if len(p) >= 2:
            seq_lengths[p[0]] = int(p[1])
print("rRNA sequence lengths:")
for name, length in seq_lengths.items():
    print(f"  {name}: {length:,} bp")

# ── Load rRNA modifications ──────────────────────────────────────────────────
mods = pd.read_csv(RRNA_MODS, sep="\t")
mods["pos"] = mods["start"] + 0.5  # midpoint of BED interval
print(f"\nTotal rRNA modifications: {len(mods)}")
print(mods["name"].value_counts().to_string())

# Keep only sequences that have modifications, ordered by sequence length
seqs_with_mods = sorted(mods["chr"].unique(), key=lambda s: seq_lengths.get(s, 0))
print(f"\nSequences with modifications: {seqs_with_mods}")

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(
    nrows=len(seqs_with_mods), ncols=1,
    figsize=(14, 1.2 * len(seqs_with_mods)),
    squeeze=False,
)
axes = axes.flatten()

for idx, seq_name in enumerate(seqs_with_mods):
    ax = axes[idx]
    seq_len = seq_lengths.get(seq_name, 0)
    seq_mods = mods[mods["chr"] == seq_name].copy()
    display_name = seq_name.replace("hs_rRNA_", "")

    # Plot subtle annotated domain backgrounds behind everything
    domains = RRNA_DOMAINS.get(seq_name, [])
    for dom_name, dom_start, dom_end, dom_color in domains:
        rect = mpatches.Rectangle(
            (dom_start, -0.08),
            dom_end - dom_start,
            0.16,
            facecolor=dom_color,
            edgecolor="none",
            linewidth=0,
            alpha=0.45,
            zorder=0,
        )
        ax.add_patch(rect)
        # Draw thin vertical lines at domain boundaries above the position line
        ax.vlines(
            [dom_start, dom_end],
            0.0,
            0.28,
            colors="#666666",
            linewidth=1.0,
            zorder=1,
        )
        # Domain label at center of domain, near the top
        label_x = (dom_start + dom_end) / 2
        label_y = 0.09
        ax.text(
            label_x, label_y, dom_name,
            ha="center", va="bottom",
            fontsize=16, fontweight="bold",
            color="black", zorder=2,
        )

        # Highlight PTC region in Domain V of human 28S rRNA
        if seq_name == "hs_rRNA_28S" and dom_name == "Domain V":
            ptc_rect = mpatches.Rectangle(
                (dom_start, -0.09),
                dom_end - dom_start,
                0.18,
                facecolor="none",
                edgecolor="red",
                linewidth=2.0,
                linestyle="--",
                zorder=4,
            )
            ax.add_patch(ptc_rect)

    # Plot rRNA sequence bar as a thin reference baseline
    ax.plot([0, seq_len], [0, 0], color="#263238", linewidth=1.5, zorder=1)

    # Plot modifications as dots colored by type with alpha
    for mod_type in sorted(seq_mods["name"].unique()):
        sub = seq_mods[seq_mods["name"] == mod_type]
        color = MOD_COLORS.get(mod_type, "#333333")
        ax.plot(
            sub["pos"],
            [0] * len(sub),
            marker="o",
            linestyle="None",
            color=color,
            markersize=8,
            markeredgecolor="none",
            alpha=0.7,
            label=mod_type,
            zorder=3,
        )

    ax.set_xlim(0, seq_len)
    ax.set_ylim(-0.15, 0.40)
    ax.set_yticks([])

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="y", which="both", left=False)

    # x-axis ticks
    ax.set_xlabel(f"{display_name} position (nt)", fontweight="bold", fontsize=16)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))

# Single figure-level legend at the bottom
legend_elements = [
    plt.Line2D([0], [0], color="none", label="Modification type:")
] + [
    plt.Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=color, markeredgecolor="none",
               markersize=10, label=mod_type)
    for mod_type, color in sorted(MOD_COLORS.items())
]
leg = fig.legend(
    handles=legend_elements,
    loc="lower center",
    ncol=len(MOD_COLORS) + 1,
    fontsize=13,
    frameon=False,
    bbox_to_anchor=(0.5, -0.02),
    handlelength=0.6,
    handletextpad=0.4,
    columnspacing=1.2,
)
if leg:
    texts = leg.get_texts()
    if texts:
        texts[0].set_fontweight("bold")
        texts[0].set_fontsize(15)

fig.text(0.01, 0.5, "Modifications", ha="center", va="center",
         rotation="vertical", fontsize=18, fontweight="bold")

plt.tight_layout(rect=[0.04, 0.08, 0.98, 0.94], h_pad=0.8)
plt.savefig(OUTDIR / "panel_g_rRNA_region_PTC.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.savefig(OUTDIR / "panel_g_rRNA_region_PTC.png", format="png", bbox_inches="tight", dpi=300)
plt.close()

print(f"\nSaved: {OUTDIR / 'panel_g_rRNA_region_PTC.pdf'}")
print(f"Saved: {OUTDIR / 'panel_g_rRNA_region_PTC.png'}")
