#!/usr/bin/env python3
"""
Histogram of polyA modification sites per protein-coding gene.

Inputs:
  - draft_reference/tiered_lists/tiered_polyA.tsv
  - ~/ref/gencode.v49.primary_assembly.annotation.gtf.gz

Outputs:
  - figures/panel_d_polyA_sites_per_gene_histogram/panel_d_polyA_sites_per_genes_histogram.pdf
  - figures/panel_d_polyA_sites_per_gene_histogram/panel_d_polyA_sites_per_genes_histogram.png
  - figures/panel_d_polyA_sites_per_gene_histogram/panel_d_top100_protein_coding_genes_by_polyA_sites.tsv
"""

import gzip
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
GITHUB_ROOT = Path(__file__).resolve().parents[2]
OUTDIR = GITHUB_ROOT / "figures" / "panel_d_polyA_sites_per_gene_histogram"
OUTDIR.mkdir(parents=True, exist_ok=True)
POLYA = GITHUB_ROOT / "outputs" / "tiered_lists" / "tiered_polyA.tsv"

# Locate the references: prefer a local copy, fall back to the reference area
LOCAL_GTF = GITHUB_ROOT / "inputs" / "gencode.v49.primary_assembly.annotation.gtf.gz"
FALLBACK_GTF = Path.home() / "ref" / "gencode.v49.primary_assembly.annotation.gtf.gz"
GTF = LOCAL_GTF if LOCAL_GTF.exists() else FALLBACK_GTF

if not GTF.exists():
    raise FileNotFoundError(f"GENCODE GTF not found at {LOCAL_GTF} or {FALLBACK_GTF}.")

# ── Plotting setup ───────────────────────────────────────────────────────────
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 4
plt.rcParams["axes.labelsize"] = 3.5
plt.rcParams["axes.titlesize"] = 4
plt.rcParams["xtick.labelsize"] = 3.5
plt.rcParams["ytick.labelsize"] = 3.5
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["font.weight"] = "bold"

# ── Load data ────────────────────────────────────────────────────────────────
polya = pd.read_csv(POLYA, sep="\t")
print(f"Total polyA sites: {len(polya):,}")

# Parse protein-coding genes
genes = []
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = p
        if feature != "gene":
            continue
        if not chrom.startswith("chr"):
            continue
        if not re.search(r'gene_type "protein_coding"', attrs):
            continue
        m_id = re.search(r'gene_id "([^"]+)"', attrs)
        m_name = re.search(r'gene_name "([^"]+)"', attrs)
        genes.append({
            "gene_id": m_id.group(1) if m_id else "",
            "gene_name": m_name.group(1) if m_name else "",
            "chrom": chrom,
            "start": int(start),
            "end": int(end),
            "strand": strand,
        })

genes_df = pd.DataFrame(genes)
print(f"Total protein-coding genes: {len(genes_df):,}")

# ── Count sites per gene (strand-matched) ────────────────────────────────────
sites_per_gene = []
for _, gene in genes_df.iterrows():
    sub = polya[polya["chr"] == gene["chrom"]]
    n = (
        (sub["start"] >= gene["start"]) &
        (sub["start"] <= gene["end"]) &
        (sub["strand"] == gene["strand"])
    ).sum()
    sites_per_gene.append(int(n))

genes_df["n_sites"] = sites_per_gene

n_with_sites = int((genes_df["n_sites"] > 0).sum())
n_without_sites = int((genes_df["n_sites"] == 0).sum())

print(f"Genes with 0 sites: {n_without_sites:,}")
print(f"Genes with ≥1 site: {n_with_sites:,}")
print(f"Mean sites per gene (all genes): {genes_df['n_sites'].mean():.2f}")
print(f"Median sites per gene (all genes): {genes_df['n_sites'].median():.1f}")

# ── Histogram ────────────────────────────────────────────────────────────────
genes_with_sites = genes_df[genes_df["n_sites"] > 0].copy()
values = genes_with_sites["n_sites"].values
max_display = 40
# Restrict to genes with >= 1 site and cap display bin
values_pos = values[values >= 1]
binned = np.minimum(values_pos, max_display)

fig, ax = plt.subplots(figsize=(1.2, 1.2))
bins = np.arange(1, max_display + 2) - 0.5
ax.hist(binned, bins=bins, color="#0D3B6E", edgecolor="white", linewidth=0.2, alpha=0.8, zorder=2)

# Show every 10th tick to avoid crowding, plus the 40+ label
xticks = list(range(10, max_display, 10)) + [max_display]
ax.set_xticks(xticks)
ax.set_xticklabels([str(i) if i < max_display else f"{max_display}+" for i in xticks])
ax.set_xlabel("Modification sites per gene", fontweight="bold", labelpad=1)
ax.set_ylabel("Number of genes", fontweight="bold", labelpad=1)



ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.3)
ax.spines["bottom"].set_linewidth(0.3)
ax.tick_params(axis="both", which="major", width=0.3, pad=1)

# ── Inset: histogram for 10–40+ mods ─────────────────────────────────────────
inset_ax = ax.inset_axes([0.40, 0.42, 0.58, 0.56])
inset_values = values[values >= 10]
inset_binned = np.minimum(inset_values, max_display)
inset_bins = np.arange(10, max_display + 2) - 0.5
inset_ax.hist(inset_binned, bins=inset_bins, color="#0D3B6E", edgecolor="white", linewidth=0.2, alpha=0.8, zorder=2)

# Show every 10th tick in the inset, with the final bin labeled 40+
inset_xticks = list(range(10, max_display, 10)) + [max_display]
inset_ax.set_xticks(inset_xticks)
inset_ax.set_xticklabels([str(i) if i < max_display else f"{max_display}+" for i in inset_xticks], fontsize=3.5)
inset_ax.set_ylabel("Genes", fontweight="bold", fontsize=3.5, labelpad=1)

inset_ax.tick_params(axis="both", labelsize=3.5, pad=1)
inset_ax.spines["top"].set_visible(True)
inset_ax.spines["right"].set_visible(True)
inset_ax.spines["left"].set_linewidth(0.3)
inset_ax.spines["bottom"].set_linewidth(0.3)
inset_ax.spines["top"].set_linewidth(0.3)
inset_ax.spines["right"].set_linewidth(0.3)

plt.subplots_adjust(left=0.25, right=0.95, bottom=0.25, top=0.95)
plt.savefig(OUTDIR / "panel_d_polyA_sites_per_genes_histogram.pdf", format="pdf", dpi=300)
plt.savefig(OUTDIR / "panel_d_polyA_sites_per_genes_histogram.png", format="png", dpi=300)
plt.close()

print(f"Saved: {OUTDIR / 'panel_d_polyA_sites_per_genes_histogram.pdf'}")
print(f"Saved: {OUTDIR / 'panel_d_polyA_sites_per_genes_histogram.png'}")

# ── Top 100 genes table ──────────────────────────────────────────────────────
top100 = genes_with_sites.sort_values("n_sites", ascending=False).head(100)

# Unique sites in top 100 genes
polya_by_chrom = {}
for idx, row in polya.iterrows():
    polya_by_chrom.setdefault(row["chr"], []).append((idx, int(row["start"])))

unique_sites_in_top100 = set()
for _, gene in top100.iterrows():
    for site_idx, pos in polya_by_chrom.get(gene["chrom"], []):
        if gene["start"] <= pos <= gene["end"]:
            unique_sites_in_top100.add(site_idx)

pct = 100 * len(unique_sites_in_top100) / len(polya)
print(f"\nUnique sites in top 100 protein-coding genes: {len(unique_sites_in_top100):,}")
print(f"Percentage of all polyA sites: {pct:.1f}%")

top100_out = OUTDIR / "panel_d_top100_protein_coding_genes_by_polyA_sites.tsv"
top100[["gene_name", "gene_id", "chrom", "start", "end", "n_sites"]].to_csv(top100_out, sep="\t", index=False)
print(f"Saved: {top100_out}")
