#!/usr/bin/env python3
"""
Scatter plot of polyA modification density vs protein-coding gene density
per chromosome. Includes linear fit line and chromosome labels.

Inputs:
  - draft_reference/tiered_lists/tiered_polyA.tsv
  - ~/ref/gencode.v49.primary_assembly.annotation.gtf.gz
  - ~/ref/GRCh38.primary_assembly.genome.fa.fai

Outputs:
  - figures/panel_c_polyA_per_mb_vs_genes/panel_c_polyA_perMb_vs_genes_perMb.pdf
  - figures/panel_c_polyA_per_mb_vs_genes/panel_c_polyA_perMb_vs_genes_perMb.png
"""

import gzip
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
GITHUB_ROOT = Path(__file__).resolve().parents[2]
OUTDIR = GITHUB_ROOT / "figures" / "panel_c_polyA_per_mb_vs_genes"
OUTDIR.mkdir(parents=True, exist_ok=True)
POLYA = GITHUB_ROOT / "outputs" / "tiered_lists" / "tiered_polyA.tsv"

# Locate the references: prefer a local copy, fall back to the reference area
LOCAL_GTF = GITHUB_ROOT / "inputs" / "gencode.v49.primary_assembly.annotation.gtf.gz"
FALLBACK_GTF = Path.home() / "ref" / "gencode.v49.primary_assembly.annotation.gtf.gz"
GTF = LOCAL_GTF if LOCAL_GTF.exists() else FALLBACK_GTF

LOCAL_FAI = GITHUB_ROOT / "inputs" / "GRCh38.primary_assembly.genome.fa.fai"
FALLBACK_FAI = Path.home() / "ref" / "GRCh38.primary_assembly.genome.fa.fai"
FAI = LOCAL_FAI if LOCAL_FAI.exists() else FALLBACK_FAI

if not FAI.exists():
    raise FileNotFoundError(f"GRCh38 .fai not found at {LOCAL_FAI} or {FALLBACK_FAI}.")
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

# ── Load chromosome lengths ──────────────────────────────────────────────────
chrom_lengths = {}
with open(FAI) as fh:
    for line in fh:
        p = line.strip().split("\t")
        if len(p) >= 2:
            chrom_lengths[p[0]] = int(p[1])

# Keep main nuclear chromosomes
main_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
main_chroms = [c for c in main_chroms if c in chrom_lengths]

# ── Load and count polyA modifications per chromosome ────────────────────────
polya = pd.read_csv(POLYA, sep="\t")
mod_counts = polya.groupby("chr").size().reindex(main_chroms, fill_value=0)

# ── Parse GTF for protein-coding genes per chromosome ────────────────────────
gene_counts = {c: 0 for c in main_chroms}
seen_genes = set()

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
        if chrom not in main_chroms:
            continue
        m_gene_id = re.search(r'gene_id "([^"]+)"', attrs)
        m_gene_type = re.search(r'gene_type "([^"]+)"', attrs)
        if not m_gene_id or not m_gene_type:
            continue
        gene_id = m_gene_id.group(1)
        gene_type = m_gene_type.group(1)
        if gene_type != "protein_coding":
            continue
        if gene_id in seen_genes:
            continue
        seen_genes.add(gene_id)
        gene_counts[chrom] += 1

# ── Build dataframe of densities ─────────────────────────────────────────────
records = []
for chrom in main_chroms:
    length_mb = chrom_lengths[chrom] / 1e6
    records.append({
        "chrom": chrom,
        "length_bp": chrom_lengths[chrom],
        "length_mb": length_mb,
        "mod_count": int(mod_counts[chrom]),
        "gene_count": gene_counts[chrom],
        "mods_per_mb": mod_counts[chrom] / length_mb,
        "genes_per_mb": gene_counts[chrom] / length_mb,
    })

df = pd.DataFrame(records)
print(df.to_string(index=False))

# Exclude chrM and chrY from regression/plot:
#  - chrM is an extreme outlier (~17 kb, very high density)
#  - chrY has zero polyA modifications in this dataset
df_fit = df[~df["chrom"].isin(["chrM", "chrY"])].copy()

# ── Linear regression ────────────────────────────────────────────────────────
slope, intercept, r_value, p_value, std_err = stats.linregress(
    df_fit["genes_per_mb"], df_fit["mods_per_mb"]
)
print(f"\nLinear fit (chrM and chrY excluded): y = {slope:.3f}x + {intercept:.3f}")
print(f"R = {r_value:.3f}, R² = {r_value**2:.3f}, p = {p_value:.4g}")

# ── Plotting with Matplotlib ──────────────────────────────────────────────────
# Sized for consistency with panel c and d: figsize=(1.2, 1.2)
fig, ax = plt.subplots(figsize=(1.2, 1.2))

# Plot linear fit line
x_min = df_fit["genes_per_mb"].min()
x_max = df_fit["genes_per_mb"].max()
x_vals = np.linspace(1.2, x_max * 1.25, 100)
y_vals = slope * x_vals + intercept
ax.plot(x_vals, y_vals, color="#C0392B", linestyle="--", linewidth=0.5, zorder=1)

# Scatter points
ax.scatter(
    df_fit["genes_per_mb"],
    df_fit["mods_per_mb"],
    facecolor="#4F7CAC",
    edgecolor="#0D3B6E",
    linewidths=0.3,
    s=3.0,
    alpha=0.9,
    zorder=2
)

# Label top 10 chromosomes by mods_per_mb
top10_idx = df_fit.sort_values("mods_per_mb", ascending=False).head(10).index
# Tailored offsets for the labels to avoid overlapping with points
offsets = {
    "chr19": (0, 18, "center", "bottom"),
    "chr17": (-18, 12, "right", "bottom"),
    "chr16": (18, -4, "left", "center"),
    "chr22": (12, -18, "left", "top"),
    "chr1":  (-20, 0, "right", "center"),
    "chr20": (-15, -15, "right", "top"),
    "chr12": (0, -22, "center", "top"),
    "chr6":  (20, 15, "left", "bottom"),
    "chr11": (-8, 20, "right", "bottom"),
    "chr14": (22, -8, "left", "top"),
}

for idx, row in df_fit.iterrows():
    if idx in top10_idx:
        chrom = row["chrom"]
        label = chrom.replace("chr", "")
        x_off, y_off, ha, va = offsets.get(chrom, (10, 10, "left", "bottom"))
        ax.annotate(
            label,
            xy=(row["genes_per_mb"], row["mods_per_mb"]),
            xytext=(x_off, y_off),
            textcoords="offset points",
            fontsize=1.0,
            fontweight="bold",
            color="black",
            ha=ha,
            va=va,
            arrowprops=dict(
                arrowstyle="-|>,head_width=0.2,head_length=0.3",
                color="black",
                lw=0.4,
                shrinkA=1,
                shrinkB=1,
                connectionstyle="arc3,rad=0.2",
            ),
            zorder=6
        )

# Labels
ax.set_xlabel("Protein-coding genes per Mb", fontweight="bold", labelpad=1)
ax.set_ylabel("PolyA modification sites per Mb", fontweight="bold", labelpad=1)

# Axis limits and spines
ax.set_xlim(1.2, x_max * 1.25)
ax.set_ylim(1.5, df_fit["mods_per_mb"].max() * 1.22)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.3)
ax.spines["bottom"].set_linewidth(0.3)
ax.tick_params(axis="both", which="major", width=0.3, pad=1)

plt.subplots_adjust(left=0.25, right=0.95, bottom=0.25, top=0.95)

pdf_out = OUTDIR / "panel_c_polyA_perMb_vs_genes_perMb.pdf"
png_out = OUTDIR / "panel_c_polyA_perMb_vs_genes_perMb.png"
plt.savefig(pdf_out, format="pdf", dpi=300)
plt.savefig(png_out, format="png", dpi=300)
plt.close()

print(f"Saved: {pdf_out}")
print(f"Saved: {png_out}")

