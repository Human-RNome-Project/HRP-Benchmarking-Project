#!/usr/bin/env python3
"""
Panel B: Zoomed locus plot for chr1:160.45-160.85 Mb showing gene models and
polyA modification sites (m6A and Y).

This is a GitHub-ready adaptation of:
  draft_reference/new_run/plot_chr1_160M_region_zoom.py

Inputs:
  - outputs/tiered_lists/tiered_polyA.tsv
  - gencode.v49.primary_assembly.annotation.gtf.gz (gene models)

Outputs:
  - figures/panel_b_chr1_160M_region_zoom/panel_b_chr1_160M_region_zoom.pdf
  - figures/panel_b_chr1_160M_region_zoom/panel_b_chr1_160M_region_zoom.png
  - figures/panel_b_chr1_160M_region_zoom/panel_b_region_modifications.tsv
  - figures/panel_b_chr1_160M_region_zoom/panel_b_region_genes.tsv
"""

import gzip
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
GITHUB_ROOT = Path(__file__).resolve().parents[2]
OUTDIR = GITHUB_ROOT / "figures" / "panel_b_chr1_160M_region_zoom"
OUTDIR.mkdir(parents=True, exist_ok=True)

POLYA = GITHUB_ROOT / "outputs" / "tiered_lists" / "tiered_polyA.tsv"

LOCAL_GTF = GITHUB_ROOT / "inputs" / "gencode.v49.primary_assembly.annotation.gtf.gz"
FALLBACK_GTF = Path.home() / "ref" / "gencode.v49.primary_assembly.annotation.gtf.gz"
GTF = LOCAL_GTF if LOCAL_GTF.exists() else FALLBACK_GTF
if not GTF.exists():
    raise FileNotFoundError(
        f"GENCODE GTF not found at {LOCAL_GTF} or {FALLBACK_GTF}. "
        "Place the GTF in inputs/ or update the path in this script."
    )

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

# ── Region ───────────────────────────────────────────────────────────────────
CHROM = "chr1"
START = 160_450_000
END = 160_850_000

# ── Official project RNA-modification colors ─────────────────────────────────
MOD_COLORS = {
    "m6A": "#721817",  # dark red
    "m5C": "#001427",  # very dark blue
    "I": "#2D6E1E",    # green
    "Y": "#F0A202",    # orange/yellow
}

# ── Load polyA sites in region ───────────────────────────────────────────────
polya = pd.read_csv(POLYA, sep="\t")
region_mods = polya[
    (polya["chr"] == CHROM) &
    (polya["start"] >= START) &
    (polya["start"] <= END)
].copy()
region_mods = region_mods.sort_values("start")
print(f"Modifications in {CHROM}:{START:,}-{END:,}: {len(region_mods)}")
print(region_mods["name"].value_counts().to_string())

# Save the modifications used in the panel
region_mods.to_csv(OUTDIR / "panel_b_region_modifications.tsv", sep="\t", index=False)

# ── Parse GTF for transcripts and exons in region ────────────────────────────
transcripts = {}
exons = []
genes = []

with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = p
        if chrom != CHROM:
            continue
        start = int(start)
        end = int(end)
        if end < START or start > END:
            continue
        m_gene_id = re.search(r'gene_id "([^"]+)"', attrs)
        m_gene_name = re.search(r'gene_name "([^"]+)"', attrs)
        m_trans_id = re.search(r'transcript_id "([^"]+)"', attrs)
        m_gene_type = re.search(r'gene_type "([^"]+)"', attrs)
        if not m_gene_id:
            continue
        gene_id = m_gene_id.group(1)
        gene_name = m_gene_name.group(1) if m_gene_name else gene_id
        gene_type = m_gene_type.group(1) if m_gene_type else ""
        trans_id = m_trans_id.group(1) if m_trans_id else None

        if feature == "gene":
            genes.append({
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_type": gene_type,
                "start": start,
                "end": end,
                "strand": strand,
            })
        elif feature == "transcript" and trans_id:
            transcripts[trans_id] = {
                "transcript_id": trans_id,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "start": start,
                "end": end,
                "strand": strand,
            }
        elif feature == "exon" and trans_id:
            exons.append({
                "transcript_id": trans_id,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "start": start,
                "end": end,
                "strand": strand,
            })

genes_df = pd.DataFrame(genes)
exons_df = pd.DataFrame(exons)
print(f"\nGenes in region: {len(genes_df)}")
print(f"Transcripts in region: {len(transcripts)}")
print(f"Exons in region: {len(exons_df)}")
print(genes_df[["gene_name", "gene_type", "start", "end", "strand"]].to_string(index=False))

# Save the genes displayed in the panel
genes_df.to_csv(OUTDIR / "panel_b_region_genes.tsv", sep="\t", index=False)

# Determine which genes overlap at least one modification (strand-matched)
def gene_overlaps_mods(gene):
    sub = region_mods[
        (region_mods["chr"] == CHROM) &
        (region_mods["start"] >= gene["start"]) &
        (region_mods["start"] <= gene["end"]) &
        (region_mods["strand"] == gene["strand"])
    ]
    return len(sub) > 0

genes_df["has_mods"] = genes_df.apply(gene_overlaps_mods, axis=1)

# Show all protein-coding genes; show lncRNAs only if they have modifications
genes_to_show = genes_df[
    (genes_df["gene_type"] == "protein_coding") |
    ((genes_df["gene_type"] == "lncRNA") & (genes_df["has_mods"]))
].copy().reset_index(drop=True)
exons_df = exons_df[exons_df["gene_id"].isin(genes_to_show["gene_id"])]

# Pick one representative transcript per gene (longest cumulative exon length)
transcript_lengths = exons_df.groupby("transcript_id").apply(
    lambda x: (x["end"] - x["start"]).sum(), include_groups=False
).reset_index(name="length")
transcript_lengths["gene_id"] = transcript_lengths["transcript_id"].map(
    {t["transcript_id"]: t["gene_id"] for t in transcripts.values()}
)
transcript_lengths = transcript_lengths.sort_values("length", ascending=False)
rep_transcripts = transcript_lengths.drop_duplicates("gene_id", keep="first")["transcript_id"].tolist()

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 2.5))

# Plot genomic sequence bar as a thin reference baseline
ax.plot([START, END], [0, 0], color="#263238", linewidth=1.5, zorder=1)

# ── Gene models ──────────────────────────────────────────────────────────────
for _, gene in genes_to_show.iterrows():
    gene_name = gene["gene_name"]
    gene_start = gene["start"]
    gene_end = gene["end"]
    gene_type = gene["gene_type"]

    visible_start = max(START, gene_start)
    visible_end = min(END, gene_end)

    box_color = "#BBDEFB" if gene_type == "protein_coding" else "#E1BEE7"
    rect = mpatches.Rectangle(
        (visible_start, -0.08),
        visible_end - visible_start,
        0.16,
        facecolor=box_color,
        edgecolor="none",
        alpha=0.6,
        zorder=0,
    )
    ax.add_patch(rect)

    label_x = (visible_start + visible_end) / 2
    label_text = f"{gene_name} →" if gene["strand"] == "+" else f"← {gene_name}"
    ax.text(
        label_x,
        0.12,
        label_text,
        ha="center",
        va="bottom",
        fontsize=15,
        fontweight="bold",
        color="black",
        zorder=4,
    )

# ── Modifications ────────────────────────────────────────────────────────────
mod_order = ["m6A", "m5C", "I", "Y"]
present_mod_types = []

for mod_type in mod_order:
    sub = region_mods[region_mods["name"] == mod_type]
    if sub.empty:
        continue
    present_mod_types.append(mod_type)
    color = MOD_COLORS.get(mod_type, "#333333")
    ax.plot(
        sub["start"],
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

# ── Axis Formatting ──────────────────────────────────────────────────────────
ax.set_xlim(START, END)
ax.set_ylim(-0.25, 0.35)
ax.set_yticks([])
ax.set_xlabel("Genomic position (Mb)", fontweight="bold", fontsize=16)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(True)
ax.spines["bottom"].set_color("#263238")
ax.spines["bottom"].set_linewidth(1.5)
ax.tick_params(axis="x", which="both", bottom=True)
ax.tick_params(axis="y", which="both", left=False)

ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x/1e6:.2f}M"))

# Single figure-level legend at the bottom for modifications
mod_legend_elements = [
    plt.Line2D([0], [0], color="none", label="Modification type:")
] + [
    plt.Line2D([0], [0], marker="o", linestyle="None",
               markerfacecolor=MOD_COLORS[mt], markeredgecolor="none",
               markersize=10, label=mt)
    for mt in present_mod_types
]
leg = fig.legend(
    handles=mod_legend_elements,
    loc="lower center",
    ncol=len(present_mod_types) + 1,
    fontsize=13,
    frameon=False,
    bbox_to_anchor=(0.5, -0.05),
    handlelength=0.6,
    handletextpad=0.4,
    columnspacing=1.2,
)
if leg:
    texts = leg.get_texts()
    if texts:
        texts[0].set_fontweight("bold")
        texts[0].set_fontsize(15)

plt.tight_layout(rect=[0.02, 0.08, 0.98, 0.98])
out_pdf = OUTDIR / "panel_b_chr1_160M_region_zoom.pdf"
out_png = OUTDIR / "panel_b_chr1_160M_region_zoom.png"
plt.savefig(out_pdf, format="pdf", bbox_inches="tight", dpi=300)
plt.savefig(out_png, format="png", bbox_inches="tight", dpi=300)
plt.close()

print(f"\nSaved: {out_pdf}")
print(f"Saved: {out_png}")
