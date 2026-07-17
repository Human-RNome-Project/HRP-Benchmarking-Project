#!/usr/bin/env python3
"""
Observed vs expected polyA modification load per gene.

Expected modification counts are obtained from a negative-binomial GLM
regressing observed counts on gene length (longest transcript) and
polyA transcript abundance (TPM). Genes above the solid y=x line carry
more modifications than predicted from length and expression alone.
Immune receptors and G-protein-coupled receptors (GPCRs) are highlighted.

Inputs:
  - draft_reference/tiered_lists/tiered_polyA.tsv
  - ~/ref/gencode.v49.primary_assembly.annotation.gtf.gz
  - data_analysis/expression_quant/OUT/OUT.gene_tpm.tsv

Outputs:
  - figures/panel_e_observed_vs_expected_mod_load/panel_e_observed_vs_expected_mod_load.pdf
  - figures/panel_e_observed_vs_expected_mod_load/panel_e_observed_vs_expected_mod_load.png
"""

import gzip
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ────────────────────────────────────────────────────────────────────
GITHUB_ROOT = Path(__file__).resolve().parents[2]
OUTDIR = GITHUB_ROOT / "figures" / "panel_e_observed_vs_expected_mod_load"
OUTDIR.mkdir(parents=True, exist_ok=True)
POLYA = GITHUB_ROOT / "outputs" / "tiered_lists" / "tiered_polyA.tsv"

# Locate the references: prefer a local copy, fall back to the reference area
LOCAL_GTF = GITHUB_ROOT / "inputs" / "gencode.v49.primary_assembly.annotation.gtf.gz"
FALLBACK_GTF = Path.home() / "ref" / "gencode.v49.primary_assembly.annotation.gtf.gz"
GTF = LOCAL_GTF if LOCAL_GTF.exists() else FALLBACK_GTF

LOCAL_TPM = GITHUB_ROOT / "inputs" / "OUT.gene_tpm.tsv"
FALLBACK_TPM = Path("/home/phoenix/projects/HRP_benchmarking_project/data_analysis/expression_quant/OUT/OUT.gene_tpm.tsv")
TPM = LOCAL_TPM if LOCAL_TPM.exists() else FALLBACK_TPM

if not GTF.exists():
    raise FileNotFoundError(f"GENCODE GTF not found at {LOCAL_GTF} or {FALLBACK_GTF}.")
if not TPM.exists():
    raise FileNotFoundError(f"TPM data not found at {LOCAL_TPM} or {FALLBACK_TPM}.")

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
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["text.color"] = "black"
plt.rcParams["axes.labelcolor"] = "black"
plt.rcParams["xtick.color"] = "black"
plt.rcParams["ytick.color"] = "black"

# ── Parse GTF for genes, transcripts and exons ───────────────────────────────
genes = {}
transcripts = {}
exons = []

print("Parsing GTF...")
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = p
        start = int(start)
        end = int(end)
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
            genes[gene_id] = {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_type": gene_type,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
            }
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
                "start": start,
                "end": end,
            })

exons_df = pd.DataFrame(exons)
print(f"Genes: {len(genes):,}, transcripts: {len(transcripts):,}, exons: {len(exons_df):,}")

# ── Compute representative transcript length per gene ────────────────────────
# Use the longest cumulative exon length among transcripts of each gene.
transcript_lengths = exons_df.groupby("transcript_id").apply(
    lambda x: (x["end"] - x["start"]).sum(), include_groups=False
).reset_index(name="transcript_length")
transcript_lengths["gene_id"] = transcript_lengths["transcript_id"].map(
    {t["transcript_id"]: t["gene_id"] for t in transcripts.values()}
)
rep_lengths = transcript_lengths.sort_values("transcript_length", ascending=False).drop_duplicates("gene_id", keep="first")
rep_lengths = rep_lengths.set_index("gene_id")["transcript_length"]

# ── Load TPM ─────────────────────────────────────────────────────────────────
print("Loading TPM...")
tpm = pd.read_csv(TPM, sep="\t")
tpm.columns = ["gene_id", "TPM"]
# GENCODE gene IDs in TPM include version (e.g. ENSG... .15); strip version for matching
tpm["gene_id_base"] = tpm["gene_id"].str.replace(r"\.\d+$", "", regex=True)
tpm = tpm.groupby("gene_id_base", as_index=False)["TPM"].mean()

# ── Build gene table ─────────────────────────────────────────────────────────
print("Building gene table...")
gene_records = []
for gene_id, g in genes.items():
    if g["gene_type"] != "protein_coding":
        continue
    gene_id_base = re.sub(r"\.\d+$", "", gene_id)
    length = rep_lengths.get(gene_id, np.nan)
    tpm_val = tpm.loc[tpm["gene_id_base"] == gene_id_base, "TPM"]
    tpm_val = tpm_val.values[0] if len(tpm_val) > 0 else np.nan
    gene_records.append({
        "gene_id": gene_id,
        "gene_id_base": gene_id_base,
        "gene_name": g["gene_name"],
        "chrom": g["chrom"],
        "start": g["start"],
        "end": g["end"],
        "strand": g["strand"],
        "length": int(length) if not pd.isna(length) else np.nan,
        "TPM": float(tpm_val) if not pd.isna(tpm_val) else 0.0,
    })

gene_df = pd.DataFrame(gene_records)

# ── Load polyA modifications and count per gene ──────────────────────────────
print("Loading polyA modifications...")
polya = pd.read_csv(POLYA, sep="\t")
print(f"Total polyA modifications: {len(polya):,}")

# Strand-matched overlap of mods with gene bodies
mod_counts = {gid: 0 for gid in gene_df["gene_id"]}
gene_bed = gene_df[["chrom", "start", "end", "strand", "gene_id"]].copy()

# Use interval-based overlap for speed
from collections import defaultdict
chrom_strand_genes = defaultdict(list)
for _, g in gene_bed.iterrows():
    chrom_strand_genes[(g["chrom"], g["strand"])].append((g["start"], g["end"], g["gene_id"]))

for _, m in polya.iterrows():
    key = (m["chr"], m["strand"])
    if key not in chrom_strand_genes:
        continue
    for gstart, gend, gid in chrom_strand_genes[key]:
        if gstart <= m["start"] <= gend:
            mod_counts[gid] += 1
            break  # assign to first overlapping gene (should be rare to overlap multiple)

# Resolve any mods that might overlap multiple genes by counting all overlaps
# (re-do more carefully)
mod_counts = {gid: 0 for gid in gene_df["gene_id"]}
for _, m in polya.iterrows():
    key = (m["chr"], m["strand"])
    if key not in chrom_strand_genes:
        continue
    for gstart, gend, gid in chrom_strand_genes[key]:
        if gstart <= m["start"] <= gend:
            mod_counts[gid] += 1

gene_df["observed"] = gene_df["gene_id"].map(mod_counts)

# ── Filter and prepare for NB regression ─────────────────────────────────────
# Keep genes with length > 0 and TPM > 0
model_df = gene_df[(gene_df["length"] > 0) & (gene_df["TPM"] > 0)].copy()
model_df["log10_length"] = np.log10(model_df["length"])
model_df["log10_TPM"] = np.log10(model_df["TPM"] + 1e-6)
model_df = model_df.reset_index(drop=True)
print(f"Genes in model: {len(model_df):,}")
print(f"Genes with ≥1 mod: {(model_df['observed'] > 0).sum():,}")

# ── Fit negative-binomial GLM ────────────────────────────────────────────────
print("Fitting negative-binomial GLM...")
formula = "observed ~ log10_length + log10_TPM"

# First fit a Poisson model to estimate overdispersion and derive NB alpha
pois_model = smf.glm(formula=formula, data=model_df, family=sm.families.Poisson()).fit()
pearson_chi2 = pois_model.pearson_chi2
df_resid = pois_model.df_resid
# NB2 alpha estimate: (Pearson chi2/df - 1) / mean(y) when > 0
mean_y = model_df["observed"].mean()
alpha_est = max(0.01, (pearson_chi2 / df_resid - 1) / mean_y) if pearson_chi2 > df_resid else 0.5
print(f"Poisson Pearson chi2/df = {pearson_chi2/df_resid:.3f}; estimated NB alpha = {alpha_est:.4f}")

# Fit NB GLM with estimated alpha
nb_model = smf.glm(
    formula=formula,
    data=model_df,
    family=sm.families.NegativeBinomial(alpha=alpha_est),
).fit()
print(nb_model.summary())

model_df["expected"] = nb_model.predict(model_df)
model_df["expected"] = model_df["expected"].clip(lower=1e-6)
model_df["fold"] = model_df["observed"] / model_df["expected"]
model_df["residual"] = model_df["observed"] - model_df["expected"]
print(f"NB GLM deviance/df = {nb_model.deviance/nb_model.df_resid:.3f}")

# ── Classify immune receptors and GPCRs by gene symbol ───────────────────────
def classify_gene_symbol(name):
    """Heuristic classification of immune receptors and GPCRs by HGNC symbol."""
    name = str(name)
    # GPCR families (partial list of common patterns)
    gpcr_prefixes = (
        "OR", "GRM", "GABBR", "ADRA", "ADRB", "CHRM", "HRH", "DRD", "HTR",
        "OPR", "TAS", "LGR", "FPR", "CCR", "CXCR", "CX3CR", "XCR", "GPR",
        "P2RY", "S1PR", "LPAR", "AGTR", "EDNRA", "EDNRB", "GHSR", "MTNR",
        "NPY", "VIPR", "CALCR", "PTH1R", "PTH2R", "GLP1R", "GIPR", "GCGR",
        "FFAR", "HCAR", "TAAR", "MAS", "MRGPR", "RXFP", "RHO", "OPS",
    )
    # Immune receptor families
    immune_prefixes = (
        "CD", "IL", "TLR", "TNFRSF", "TNF", "KLRC", "KLRK", "KLRD", "KLRB",
        "FCGR", "FCRL", "FCER", "NCR", "NKG", "SLAMF", "SIGLEC", "LY9",
        "CD48", "CD244", "CD28", "CTLA4", "PDCD1", "ICOS", "BTLA", "LAG3",
        "HAVCR", "TRBC", "TRAC", "TRDC", "TRGC", "TRBV", "TRAV", "TRDV",
        "IGH", "IGL", "IGK", "FCGRT", "FCMR", "PILR", "CEACAM", "CD300",
        "LILR", "ILT", "KIR", "NCR", "DECTIN", "CLEC", "LY96", "CD14",
    )
    is_gpcr = name.startswith(gpcr_prefixes) or name in (
        "PTGDR", "PTGER", "PTGFR", "TBXA2R", "LPAR", "S1PR"
    )
    is_immune = name.startswith(immune_prefixes) or name in (
        "CD48", "CD244", "LY9", "SLAMF1", "SLAMF6", "SLAMF7", "CTLA4", "PDCD1",
        "ICOS", "BTLA", "LAG3", "HAVCR2", "TIGIT", "CD274", "PDCD1LG2"
    )
    return is_immune, is_gpcr

model_df[["is_immune", "is_gpcr"]] = model_df["gene_name"].apply(
    lambda x: pd.Series(classify_gene_symbol(x))
)

n_immune = model_df["is_immune"].sum()
n_gpcr = model_df["is_gpcr"].sum()
both = model_df["is_immune"] & model_df["is_gpcr"]
n_both = both.sum()
print(f"\nImmune receptor genes: {n_immune:,}")
print(f"GPCR genes: {n_gpcr:,}")
print(f"Both (e.g. chemokine receptors): {n_both:,}")

# ── Enrichment among over-modified genes ─────────────────────────────────────
# Define over-modified as observed > 2 * expected
over_mod = model_df[model_df["observed"] > 2 * model_df["expected"]].copy()
print(f"\nOver-modified genes (>2-fold above expected): {len(over_mod):,}")

from scipy.stats import fisher_exact
def enrich(category_col, label):
    a = ((model_df[category_col]) & (model_df["observed"] > 2 * model_df["expected"])).sum()
    b = (model_df[category_col] & ~(model_df["observed"] > 2 * model_df["expected"])).sum()
    c = (~model_df[category_col] & (model_df["observed"] > 2 * model_df["expected"])).sum()
    d = (~model_df[category_col] & ~(model_df["observed"] > 2 * model_df["expected"])).sum()
    table = [[a, b], [c, d]]
    oddsratio, pval = fisher_exact(table, alternative="greater")
    print(f"{label}: {a}/{b+a} over-modified vs {c}/{c+d} others; OR={oddsratio:.2f}, p={pval:.3g}")
    return oddsratio, pval, a, b+a

print("\nEnrichment among >2-fold over-modified genes:")
immune_or, immune_p, immune_n_over, immune_n = enrich("is_immune", "Immune receptors")
gpcr_or, gpcr_p, gpcr_n_over, gpcr_n = enrich("is_gpcr", "GPCRs")

# ── Plot ─────────────────────────────────────────────────────────────────────
print("\nGenerating plot...")
fig, ax = plt.subplots(figsize=(1.2, 1.2))

# Use pseudocounts for log-scale plotting
model_df["observed_plot"] = model_df["observed"] + 0.1
model_df["expected_plot"] = model_df["expected"] + 0.1

# Jitter points slightly for visualization
np.random.seed(42)
model_df["expected_plot"] *= np.random.uniform(0.90, 1.11, size=len(model_df))
model_df["observed_plot"] *= np.random.uniform(0.90, 1.11, size=len(model_df))

# Reference lines (log scale so 2-fold lines are parallel to y=x)
# y = x (solid)
line_min = min(model_df["expected_plot"].min(), model_df["observed_plot"].min())
line_max = max(model_df["expected_plot"].max(), model_df["observed_plot"].max())
line_lims = [line_min, line_max]
ax.plot(line_lims, line_lims, color="#2C3E50", linewidth=0.5, linestyle="-", zorder=2)
# y = 2x and y = x/2 (dotted) — parallel on log-log axes
ax.plot(line_lims, [2 * x for x in line_lims], color="#7F8C8D", linewidth=0.3, linestyle=":", zorder=2)
ax.plot(line_lims, [0.5 * x for x in line_lims], color="#7F8C8D", linewidth=0.3, linestyle=":", zorder=2)

# Helper to plot subsets based on residual and classification
def plot_subset(subset, marker, label_prefix):
    # Residual <= 4 (gray)
    sub_gray = subset[subset["residual"] <= 4]
    if not sub_gray.empty:
        ax.scatter(
            sub_gray["expected_plot"],
            sub_gray["observed_plot"],
            marker=marker,
            s=1.5,
            c="#B0BEC5",
            edgecolors="none",
            alpha=0.35,
            zorder=3,
        )
    # Residual > 4 (red)
    sub_red = subset[subset["residual"] > 4]
    if not sub_red.empty:
        ax.scatter(
            sub_red["expected_plot"],
            sub_red["observed_plot"],
            marker=marker,
            s=2.5,
            c="#E74C3C",
            edgecolors="#922B21",
            linewidths=0.2,
            alpha=0.65,
            zorder=4,
        )

# Group and plot
other_genes = model_df[~(model_df["is_immune"] | model_df["is_gpcr"])]
immune_genes = model_df[model_df["is_immune"]]
gpcr_genes = model_df[model_df["is_gpcr"]]

plot_subset(other_genes, "o", "Other")
plot_subset(immune_genes, "^", "Immune")
plot_subset(gpcr_genes, "s", "GPCR")

ax.set_xscale("log")
ax.set_yscale("log")
# Independent padded limits to give top-10 labels room
x_max = model_df["expected_plot"].max()
y_max = model_df["observed_plot"].max()
ax.set_xlim(0.08, x_max * 2.5)
ax.set_ylim(0.08, y_max * 1.5)

# Format axis ticks with actual numbers instead of powers of 10
from matplotlib.ticker import FuncFormatter, NullLocator
def log_formatter(x, pos):
    if x >= 1:
        return f"{int(x)}"
    elif x >= 0.1:
        return f"{x:.1f}"
    else:
        return f"{x}"
ax.xaxis.set_major_formatter(FuncFormatter(log_formatter))
ax.yaxis.set_major_formatter(FuncFormatter(log_formatter))
ax.set_xticks([0.1, 1, 10])
ax.set_yticks([0.1, 1, 10, 100])
ax.xaxis.set_minor_locator(NullLocator())
ax.yaxis.set_minor_locator(NullLocator())

ax.set_xlabel("Expected modification count", fontweight="bold", labelpad=1)
ax.set_ylabel("Observed modification count", fontweight="bold", labelpad=1)
# No title for publication-style figure

# Custom legend explaining colors (grey and red) - placed at upper left to avoid dots at 0.1
legend_elements = [
    plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="#B0BEC5", markeredgecolor="none", markersize=2, alpha=0.35, label=r"Residual $\leq$ 4"),
    plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="#E74C3C", markeredgecolor="#922B21", markersize=2, alpha=0.65, label="Residual > 4"),
]

ax.legend(
    handles=legend_elements,
    loc="lower right",
    bbox_to_anchor=(1.0, 0.15),
    fontsize=3,
    frameon=False,
    handletextpad=0.1,
    borderaxespad=0.1
)

# Annotate top 10 genes by observed modification count
# Place each gene name next to its point with tailored offsets to avoid overlap
top10 = model_df.nlargest(10, "observed").copy()
top10_offsets = {
    "BCL2":    (0,   18, "center", "bottom"),
    "TLR10":   (-18, 8,  "right",  "center"),
    "SLAMF1":  (0,   -18, "center", "top"),
    "GPR183":  (15,  16, "left",   "bottom"),
    "GPR15":   (-16, -22, "right",  "top"),
    "APOL1":   (15,  -15, "left",   "top"),
    "CCR7":    (18,  0,  "left",   "center"),
    "BCHE":    (-20, 0,  "right",  "center"),
    "RHOH":    (-15, -16, "right",  "top"),
    "QPRT":    (15,  -16, "left",   "top"),
}
for _, row in top10.iterrows():
    name = row["gene_name"]
    x_off, y_off, ha, va = top10_offsets.get(
        name, (10, 10, "left", "bottom")
    )
    # Scale offsets slightly to make room for the arrow
    ax.annotate(
        name,
        xy=(row["expected_plot"], row["observed_plot"]),
        xytext=(x_off * 0.5, y_off * 0.5),
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
            connectionstyle="arc3,rad=0.08",
        ),
        zorder=6,
    )

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.3)
ax.spines["bottom"].set_linewidth(0.3)
ax.tick_params(axis="both", which="major", width=0.3, pad=1)

plt.subplots_adjust(left=0.25, right=0.95, bottom=0.25, top=0.95)
plt.savefig(OUTDIR / "panel_e_observed_vs_expected_mod_load.pdf", format="pdf", dpi=300)
plt.savefig(OUTDIR / "panel_e_observed_vs_expected_mod_load.png", format="png", dpi=300)
plt.savefig(OUTDIR / "panel_e_observed_vs_expected_mod_load_jitter.pdf", format="pdf", dpi=300)
plt.savefig(OUTDIR / "panel_e_observed_vs_expected_mod_load_jitter.png", format="png", dpi=300)
plt.close()

print(f"\nSaved: {OUTDIR / 'panel_e_observed_vs_expected_mod_load.pdf'}")
print(f"Saved: {OUTDIR / 'panel_e_observed_vs_expected_mod_load.png'}")