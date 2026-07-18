#%%
# ============================================================================
# STEP 0 — Imports
# ----------------------------------------------------------------------------
# Two dataframe engines are used side by side on purpose:
#   * pandas  -> needed for PyRanges and seaborn interop (they speak pandas)
#   * polars  -> used for the fast / memory-light filtering + joins on the large
#                per-site BED tables (millions of rows)
# The rest of the stack:
#   * matplotlib / seaborn  -> plotting
#   * pyranges              -> genomic interval overlap (join sites <-> gene models)
#   * scipy.stats           -> hypergeometric / Fisher test for native-vs-IVT calls
#   * statsmodels           -> Benjamini-Hochberg multiple-testing correction
#   * upsetplot             -> set-overlap visualisation
#   * tqdm                  -> progress bars for the slow per-row annotation loops
#   * dmode                 -> the project's in-house helper package (metagene utils etc.)
# ============================================================================
import pandas as pd
import polars as pl
import os
import dmode 
import numpy as np
import matplotlib.pyplot as plt
import pyranges as pr
import seaborn as sns
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests
from upsetplot import UpSet, from_contents
from tqdm import tqdm

#%%
# ============================================================================
# STEP 1 — Central plotting theme ("RNome" / Nature house style)
# ----------------------------------------------------------------------------
# This block is effectively a self-contained module. It defines:
#   * a fixed colour vocabulary for every RNA modification (grouped by base),
#   * platform colours (LRS / NGS / MS-seq / Combined),
#   * Nature figure sizing helpers,
#   * a set_rnome_theme() function that pushes everything into matplotlib rcParams,
#   * small helpers for panel labels and vector-only saving.
# Centralising this guarantees every figure in the paper looks identical.
# ============================================================================
"""
RNome paper visualisation theme.
Encodes Nature figure guidelines (updated 2026-06-03) as a seaborn theme.

Usage:
    from rnome_theme import set_rnome_theme, MODIFICATION_COLORS, PLATFORM_COLORS, \
        PLATFORM_ALPHAS, FAMILY, MODIFICATION_ORDER, mm2in, NATURE_WIDTHS, NATURE_MAX_HEIGHT
    set_rnome_theme()
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------------------------------------------------------------
# Color codes
# ----------------------------------------------------------------------
# Logic: one hue family per canonical base (A/C/G/U). Within a family the shade
# goes dark -> light so that related modifications read as "the same colour group"
# in a legend while still being distinguishable. "reserved_*" and "*?" entries are
# placeholders/ambiguous calls kept so the palette length is stable.
MODIFICATION_COLORS = {
    # A — m6A family (deep crimson → light blush)
    "m6A":   "#721817",
    "Am":    "#D44F3E",
    "m1A":   "#A52020",
    "mA?":   "#E8907A",
    "m6,6A": "#F5C4B8",
    "I":     "#6E1E1F",
    # C — m5C family (near-black navy → pale sky)
    "m5C":   "#001427",
    "Cm":    "#0D3B6E",
    "mC?":   "#1E6EB5",
    "ac4c":  "#6AAED6",
    "reserved_C_3": "#B8D9F0",
    # G — Inosine family (dark forest → pale mint)

    "Gm":    "#74B354",
    "mG?":   "#4A8532",
    "m2,2,7G": "#A8D48A",
    "reserved_G_3": "#D4EEC4",
    # U — Psi family (deep amber → pale gold)
    "Y":     "#F0A202",
    "Um":    "#C47A02",
    "mU?":   "#F5BE45",
    "reserved_U_2": "#F9D47E",
    "reserved_U_3": "#FDE9B8",
}

# Colours for the three sequencing/measurement platforms plus their merged set.
PLATFORM_COLORS = {
    "LRS":    "#6EA359",
    "NGS":    "#2261AA",
    "MS-seq": "#E2813C",
    "Combined": "#515963"
}

# Semi-transparency used when platform layers are overplotted (so overlaps show).
PLATFORM_ALPHAS = {"LRS": 0.5, "NGS": 0.5, "MS-seq": 0.5}

# Reverse index: which modifications belong to which base. Used to build the
# canonical ordering below and to group things in plots.
FAMILY = {
    "A": ["m6A", "Am", "m1A", "mA?", "m6,6A"],
    "C": ["m5C", "Cm", "mC?", "ac4c", "reserved_C_3"],
    "G": ["I", "Gm", "mG?", "m2,2,7G", "reserved_G_3"],
    "U": ["Y", "Um", "mU?", "reserved_U_2", "reserved_U_3"],
}

# Canonical modification order: by base (A,C,G,U) then lexicographic within family.
# This single ordering is reused as the default categorical/palette order so bars,
# heatmap rows, legends etc. are always sorted the same way across every figure.
MODIFICATION_ORDER = [m for base in ("A", "C", "G", "U")
                      for m in sorted(FAMILY[base], key=str.lower)]

# ----------------------------------------------------------------------
# Sizing helpers (Nature)
# ----------------------------------------------------------------------
# Nature specifies column widths in millimetres; matplotlib wants inches.
def mm2in(mm):
    return mm / 25.4

NATURE_WIDTHS = {"single": mm2in(89), "double": mm2in(183)}  # inches
NATURE_MAX_HEIGHT = mm2in(170)  # inches

# ----------------------------------------------------------------------
# Theme
# ----------------------------------------------------------------------
def set_rnome_theme(context="paper"):
    """Apply the RNome/Nature seaborn theme. Fonts 5–7 pt, embedded TrueType."""
    # Start from a clean seaborn base (ticks style, Helvetica) then override rc.
    sns.set_theme(context=context, style="ticks",
                  font="Helvetica", rc=None)
    rc = {
        # Fonts — sans-serif, Helvetica/Arial; embed TrueType 42 so text stays
        # editable/selectable in the vector PDF (fonttype 42 = TrueType embedding).
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",

        # Text sizes: max 7 pt, min 5 pt (panel labels 8 pt bold set manually)
        # NOTE: the docstring/Nature guideline says 5–7 pt, but the values below
        # are ~20–28 pt. These are the *actual* sizes applied. If true Nature sizing
        # is wanted these need lowering; left as-is per "keep everything as it is".
        "font.size": 20,
        "axes.labelsize": 28,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "legend.fontsize": 20,
        "legend.title_fontsize": 22,
        "figure.titlesize": 24,

        # No coloured text — Nature wants all text/axes black.
        "text.color": "black",
        "axes.labelcolor": "black",
        "axes.titlecolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.edgecolor": "black",

        # Thin lines suited to small print (hairlines read cleanly when shrunk).
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "lines.linewidth": 0.75,
        "grid.linewidth": 0.4,
        "patch.linewidth": 0.5,

        # Default to single-column figure at print resolution.
        "figure.figsize": (NATURE_WIDTHS["single"], mm2in(60)),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.format": "pdf",
        "savefig.bbox": "tight",
        "savefig.transparent": False,

        # No grid, and drop the top/right spines (cleaner "ticks" look).
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
    mpl.rcParams.update(rc)
    # Default categorical palette = canonical modification order, so any plot that
    # doesn't pass an explicit palette still gets the house colours in order.
    sns.set_palette([MODIFICATION_COLORS[m] for m in MODIFICATION_ORDER])

def panel_label(ax, letter, x=-0.15, y=1.05):
    """Add a Nature panel label: 8-pt bold, upright, lowercase."""
    # Placed in axes-fraction coordinates so it sits just outside the top-left corner
    # regardless of data ranges — this is the "a/b/c" tag on multi-panel figures.
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=8, fontweight="bold", fontstyle="normal",
            va="bottom", ha="right", color="black")

def save_nature(fig, path):
    """Save as vector (.pdf/.svg/.eps/.ai). Rejects raster formats."""
    # Guardrail: Nature only accepts vector figures, so refuse raster extensions
    # rather than silently producing a non-compliant file.
    bad = (".png", ".jpg", ".jpeg", ".tiff", ".tif")
    if path.lower().endswith(bad):
        raise ValueError("Nature requires vector formats (.pdf/.svg/.eps/.ai); "
                         "raster formats are not accepted.")
    fig.savefig(path)

# Apply the theme immediately so every figure produced later inherits it.
set_rnome_theme()

#%%
# ============================================================================
# STEP 2 — Load the genome annotation (GTF) and build gene-body coordinates
# ----------------------------------------------------------------------------
# prepare_gene_body_coverage() (from the in-house dmode package) parses the
# GENCODE GTF into:
#   * gtf_df       -> a tidy table of all annotation records (genes/exons/CDS/UTR...)
#   * gene_body_df -> collapsed gene-body intervals, used later to assign each
#                     modification site to a gene / metagene position.
# ============================================================================
gtf_file = "/home/stefan/Synology/Data_jamboree/final_bedRmods/gencode.v49.primary_assembly.annotation.gtf"
gtf_df, gene_body_df = dmode.metagene_plot.prepare_gene_body_coverage(gtf_file)

# %%
# ============================================================================
# STEP 3 — Load the combined ONT (nanopore) per-site modification BED
# ----------------------------------------------------------------------------
# This is the bedRMod/bedMethyl-style table: one row per detected modified site,
# with coverage, modification frequency, the single-letter mod code, etc.
# chrY is dropped up front (GM12878 is female; chrY calls would be spurious).
# The GTF chromosome set is captured so we can later tell "genome" contigs apart
# from the spiked-in ribosomal reference contigs.
# ============================================================================
bedfile = pd.read_csv("/global/cfs/cdirs/m5243/final_bedRmods/ont/ONT_polyARNA_rRNA_combined.filtered.bed",sep="\t", header=None,names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","coverage","frequency","single_letter_code","mod_id"], comment="#")
chromosomes_in_gtf = set(gtf_df['seqname'].unique())
bedfile = bedfile[bedfile["chrom"] != "chrY"]
#%%
bedfile.shape
# ----------------------------------------------------------------------------
# Split the BED into two logical sub-tables (polars used here for speed):
#   * ribosomal_bedfile -> contigs NOT in the GTF whose name starts with "hs_"
#                          (these are the custom human rRNA reference sequences)
#   * genome_bedfile    -> contigs that ARE in the GTF (the real genome)
# ----------------------------------------------------------------------------
ribosomal_bedfile = (
    pl.from_pandas(bedfile)
    .with_columns(pl.col("chrom").cast(pl.String))
    .filter(
        ~pl.col("chrom").is_in(chromosomes_in_gtf) &
        pl.col("chrom").str.starts_with("hs_")
    )
)

genome_bedfile = (
    pl.from_pandas(bedfile)
    .with_columns(pl.col("chrom").cast(pl.String))
    .filter(
        pl.col("chrom").is_in(chromosomes_in_gtf)
    )
)

#%%
# ============================================================================
# STEP 4 — Overlap every modification site with gene bodies (PyRanges join)
# ----------------------------------------------------------------------------
# PyRanges needs the canonical column names Chromosome/Start/End/Strand.
# A LEFT join keeps every modification site even if it falls outside any gene
# body (those rows get feature == -1, handled in the next step).
# ============================================================================
# Create PyRanges objects
bed_ranges = pr.PyRanges(bedfile.rename(columns={
    "chrom": "Chromosome",
    "chromStart": "Start",
    "chromEnd": "End",
    "strand": "Strand",
}))



gene_body_ranges = pr.PyRanges(gene_body_df.to_pandas().rename(columns={
    "seqname": "Chromosome",
    "start": "Start",
    "end": "End",
    "strand": "Strand",
}))
# Find overlaps
overlaps = bed_ranges.join(gene_body_ranges, how="left")
# Convert overlaps to DataFrame
overlaps_df = overlaps.df

# %%
# ============================================================================
# STEP 5 — Classify each site into a genomic feature category
# ----------------------------------------------------------------------------
# This is the core annotation decision tree. For every overlap row:
#   * If the gene-body join produced a real feature -> keep it as-is.
#   * If feature == "-1" (no gene-body overlap), decide what the site is:
#       - rRNA sub-species detected by substring in the contig name (18S/28S/5S/5.8S)
#       - otherwise look up a *containing* gene in the full GTF (start<=site<=end,
#         same strand). If none -> "intergenic".
#       - if a containing gene exists and is protein_coding, refine into
#         CDS / 5'UTR / 3'UTR using the record's feature + tags; anything else
#         protein_coding but unresolved -> "undefined intragenic".
#       - non-coding genes -> coarse bucket "lncRNA/rRNA/snoRNA/other", while the
#         *detailed* column records the exact gene_type.
# Two parallel lists are built: a coarse label and a fine-grained ("detailed") one.
# NOTE: candidate_gene[...][-1] takes the LAST matching GTF record; when several
# annotations overlap, the last one wins (a deliberate but arbitrary tie-break).
# ============================================================================
new_features = []
new_features_detailed = []
for feature,chrom,start,end,strand in tqdm(zip(overlaps_df["feature"], 
                                          overlaps_df["Chromosome"], 
                                          overlaps_df["Start"], 
                                          overlaps_df["End"], 
                                          overlaps_df["Strand"]
                                          ), total=len(overlaps_df)):
    if feature == "-1":
        # No gene-body overlap: first try to name the rRNA species from the contig.
        if chrom.__contains__("18S"):
            new_features.append("18S")
            new_features_detailed.append("18S")
        elif chrom.__contains__("28S"):
            new_features.append("28S")
            new_features_detailed.append("28S")
        elif chrom.__contains__("5S"):
            new_features.append("5S")
            new_features_detailed.append("5S")
        elif chrom.__contains__("5.8S"):
            new_features.append("5.8S")
            new_features_detailed.append("5.8S")
        else:
            # Not rRNA: search the full GTF for a gene that fully contains this site.
            candidate_gene = gtf_df.filter(
                (pl.col("seqname") == chrom) &
                (pl.col("start") <= start) &
                (pl.col("end") >= end) &
                (pl.col("strand") == strand)
            )
            if candidate_gene.is_empty():
                new_features.append("intergenic")
                new_features_detailed.append("intergenic")
            else:
                # Protein-coding + a UTR/CDS record -> resolve exact coding sub-region.
                if candidate_gene["gene_type"][-1] == "protein_coding" and candidate_gene["feature"][-1] in ["UTR","CDS"]:
                    if candidate_gene["feature"][-1] == "CDS":
                        new_features.append("CDS")
                        new_features_detailed.append("CDS")
                    elif candidate_gene["tags"][-1].__contains__("CCDS"):
                        # CCDS tag implies coding sequence even if labelled generically.
                        new_features.append("CDS")
                        new_features_detailed.append("CDS")
                    elif candidate_gene["tags"][-1].__contains__("5_UTR"):
                        new_features.append("5'UTR")
                        new_features_detailed.append("5'UTR")
                    elif candidate_gene["tags"][-1].__contains__("3_UTR"):
                        new_features.append("3'UTR")
                        new_features_detailed.append("3'UTR")
                    else:
                        new_features.append("undefined intragenic")
                        new_features_detailed.append("undefined intragenic")
                elif candidate_gene["gene_type"][-1] == "protein_coding":
                    # Coding gene but the record isn't a clean UTR/CDS -> unresolved.
                    new_features.append("undefined intragenic")
                    new_features_detailed.append("undefined intragenic")
                else:
                    # Any non-coding gene: coarse bucket + keep the true type in detail.
                    new_features.append("lncRNA/rRNA/snoRNA/other")
                    new_features_detailed.append(candidate_gene["gene_type"][-1])
    else:
        # The gene-body join already gave a usable feature -> pass it through.
        new_features.append(feature)
        new_features_detailed.append(feature)
overlaps_df["feature"] = new_features
overlaps_df["feature_detailed"] = new_features_detailed

# %%
# ----------------------------------------------------------------------------
# Re-split the *annotated* overlaps the same way as before (rRNA vs genome),
# now that every row carries a feature label.
# ----------------------------------------------------------------------------
ribosomal_df = (
    pl.from_pandas(overlaps_df)
    .with_columns(pl.col("Chromosome").cast(pl.String))
    .filter(
        ~pl.col("Chromosome").is_in(chromosomes_in_gtf) &
        pl.col("Chromosome").str.starts_with("hs_")
    )
)

genome_df = (
    pl.from_pandas(overlaps_df)
    .with_columns(pl.col("Chromosome").cast(pl.String))
    .filter(
        pl.col("Chromosome").is_in(chromosomes_in_gtf)
    )
)

#%%
# ----------------------------------------------------------------------------
# Convenience views with IVT (in-vitro-transcribed, unmodified control) removed.
# score == 0 marks the non-IVT (native) rows here.
# ----------------------------------------------------------------------------
ribosomal_without_ivt_df = ribosomal_df.filter(pl.col("score") == 0)
genome_without_ivt_df = genome_df.filter(pl.col("score") == 0)

#%%
# ----------------------------------------------------------------------------
# Define the canonical left-to-right ordering of feature categories for plots,
# and derive a sorted feature list that puts the "interesting" transcript regions
# (5'UTR, CDS, 3'UTR, intronic) first and everything else alphabetically after.
# ----------------------------------------------------------------------------
category_order = ["5'UTR", "CDS", "3'UTR", "undefined intragenic","lncRNA/rRNA/snoRNA/other", "intergenic"]
genome_df_features = list(set(genome_df["feature"]))
priority = ["5'UTR", "CDS", "3'UTR", "intronic"]  # entries you want in front
genome_df_features = list(set(genome_df["feature"]))
# Sort key: (rank-in-priority if present else after all priorities, then name).
genome_df_features.sort(key=lambda x: (priority.index(x) if x in priority else len(priority), x))

#%%
# ============================================================================
# STEP 6 — Load the 15 individual lab datasets (ONT direct-RNA, polyA+)
# ----------------------------------------------------------------------------
# Each lab contributed a re-basecalled, frequency-filtered per-site BED. They all
# share the same schema (note this schema differs from STEP 3: it has the split
# per-call counts). For each file we:
#   * tag a fixed detection threshold (score_nmod_thresholds = 14),
#   * tag the lab ID,
#   * convert to polars for the later concat/join.
# NOTE: the trailing "#HRP_..." comments are the author's dataset codes. A couple
# don't match the assigned `lab` string (Lab12 comment says A016 but lab=A_011_1;
# Lab15 comment says A_015 but lab=A_003_1) — worth double-checking, left as-is.
# ============================================================================
#HRP_A_007_1
Lab1_df = pd.read_csv("polyA_SUP_CHEUNG_V_pod5_ONT_HRP_GM12878_TOT_RNA_DRS_20260112_1.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab1_df["score_nmod_thresholds"] = 14
Lab1_df["lab"] = "HRP_A_007_1"
Lab1_df = pl.from_pandas(Lab1_df)

#HRP_A_013_1
Lab2_df = pd.read_csv("polyA_SUP_CONTICELLO_S_ONT_totalRNA_polyAenriched.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab2_df["score_nmod_thresholds"] = 14
Lab2_df["lab"] = "HRP_A_013_1"
Lab2_df = pl.from_pandas(Lab2_df)

#HRP_A_001_1
Lab3_df = pd.read_csv("polyA_SUP_DIETERICH_C_GM12878_native_RTA.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab3_df["score_nmod_thresholds"] = 14
Lab3_df["lab"] = "HRP_A_001_1"
Lab3_df = pl.from_pandas(Lab3_df)

#HRP_A_002_1
Lab4_df = pd.read_csv("polyA_SUP_GALLO_A_AG_poly_a_GM12878.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab4_df["score_nmod_thresholds"] = 14
Lab4_df["lab"] = "HRP_A_002_1"
Lab4_df = pl.from_pandas(Lab4_df)


#HRP_A_006_4
Lab5_df = pd.read_csv("polyA_SUP_GERBER_S_pod5_ONT_Gerber_native_mRNA_251209_r1_v1.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab5_df["score_nmod_thresholds"] = 14
Lab5_df["lab"] = "HRP_A_006_4"
Lab5_df = pl.from_pandas(Lab5_df)

#HRP_A_005_1
Lab6_df = pd.read_csv("polyA_SUP_GOEKE_J_ONT_1_SQKRNA004_20251125_rep1.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab6_df["score_nmod_thresholds"] = 14
Lab6_df["lab"] = "HRP_A_005_1"
Lab6_df = pl.from_pandas(Lab6_df)

#HRP_A_014_1
Lab7_df = pd.read_csv("polyA_SUP_Ohtan_Wang_pod5_ONT_Wang_totalRNA_dRNA_Jan31_rep1_version1.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab7_df["score_nmod_thresholds"] = 14
Lab7_df["lab"] = "HRP_A_014_1"
Lab7_df = pl.from_pandas(Lab7_df)

#HRP_A_009_1
Lab8_df = pd.read_csv("polyA_SUP_ONT_ont_polyA_rna_direct_rna_20260223_rep1.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab8_df["score_nmod_thresholds"] = 14
Lab8_df["lab"] = "HRP_A_009_1"
Lab8_df = pl.from_pandas(Lab8_df)

#HRP_A_009_2
Lab9_df = pd.read_csv("polyA_SUP_ONT_ont_polyA_rna_direct_rna_20260223_rep2.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab9_df["score_nmod_thresholds"] = 14
Lab9_df["lab"] = "HRP_A_009_2"
Lab9_df = pl.from_pandas(Lab9_df)

#HRP_A010_1
Lab10_df = pd.read_csv("polyA_SUP_PICARDI_promethion_RNome-uniba_dRNAtotal_200226_1_0.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab10_df["score_nmod_thresholds"] = 14
Lab10_df["lab"] = "HRP_A_010_1"
Lab10_df = pl.from_pandas(Lab10_df)


#HRP_A010_2
# NOTE: this file uses the "greater_null_polyA" prefix (no "_freq_") unlike the
# datasets above — just a naming difference in the upstream pipeline output.
Lab11_df = pd.read_csv("polyA_SUP_PICARDI_promethion_RNome-uniba_dRNAtotal_200226_2_0.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab11_df["score_nmod_thresholds"] = 14
Lab11_df["lab"] = "HRP_A_010_2"
Lab11_df = pl.from_pandas(Lab11_df)

#HRP_A016_1
Lab12_df = pd.read_csv("polyA_SUP_RAGOUSSIS_J_ONT_GM12878_directRNA_01022026_01.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab12_df["score_nmod_thresholds"] = 14
Lab12_df["lab"] = "HRP_A_011_1"
Lab12_df = pl.from_pandas(Lab12_df)

# Same lab/donor as Lab12, replicate 2 (RAGOUSSIS 10-02-2026 run).
Lab13_df = pd.read_csv("polyA_SUP_RAGOUSSIS_J_ONT_GM12878_directRNA_10022026_02.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab13_df["score_nmod_thresholds"] = 14
Lab13_df["lab"] = "HRP_A_011_2"
Lab13_df = pl.from_pandas(Lab13_df)


# Same lab/donor as Lab12, replicate 3 (RAGOUSSIS 12-02-2026 run).
Lab14_df = pd.read_csv("polyA_SUP_RAGOUSSIS_J_ONT_GM12878_directRNA_12022026_03.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab14_df["score_nmod_thresholds"] = 14
Lab14_df["lab"] = "HRP_A_011_3"
Lab14_df = pl.from_pandas(Lab14_df)


#HRP_A_015_1
Lab15_df = pd.read_csv("HRP_A_003_native_polyA_RNA_001.merged.aligned.sorted.genome.bed", sep="\t", comment="#", header=None, names=["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","color","coverage","frequency","n_mod","count_canonical","count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
Lab15_df["score_nmod_thresholds"] = 14
Lab15_df["lab"] = "HRP_A_003_1"
Lab15_df = pl.from_pandas(Lab15_df)




#%%
# ----------------------------------------------------------------------------
# Stack all 15 labs into one long table (vertical concat) and drop chrY again
# (the per-lab files may still contain it).
# ----------------------------------------------------------------------------
combined_labs_df = pl.concat([Lab1_df, Lab2_df, Lab3_df, Lab4_df, Lab5_df, Lab6_df, Lab7_df, Lab8_df, Lab9_df, Lab10_df, Lab11_df, Lab12_df, Lab13_df, Lab14_df, Lab15_df], how="vertical")
combined_labs_df = combined_labs_df.filter(pl.col("chrom") != "chrY")
# %%
# ============================================================================
# STEP 7 — Composite "confidence" score per site (weighted_sum)
# ----------------------------------------------------------------------------
# Idea: combine three normalised signals into a single 0–1 score so sites can be
# ranked/filtered on one number. Each input is min-max scaled to [0,1]:
#   * coverage  -> log2(x+1) first (heavy-tailed) then min-max
#   * n_mod     -> same log2 treatment (computed but not used in the final sum)
#   * score_nmod_thresholds -> divided by its max possible value (15)
#   * frequency -> percent -> fraction (/100)
# The final weighted_sum is an equal-weight (1/3 each) average of coverage,
# frequency, and the threshold score. Higher = more trustworthy detection.
# ============================================================================
def compute_weighted_sums(df):
    def _minmax_expr(col: str) -> pl.Expr:
        # log-transform the highly skewed columns before scaling; others as-is.
        if col == "coverage" or col == "n_mod":
            c = (pl.col(col) + 1).log(base=2)
        else:
            c = pl.col(col)
        return ((c - c.min()) / (c.max() - c.min())).alias(f"std_{col}")
    df = df.with_columns([
            _minmax_expr("coverage").alias("std_coverage"),
            _minmax_expr("n_mod").alias("std_n_mod"),
            (pl.col("score_nmod_thresholds") / 15).alias("std_score_nmod_thresholds"),
            (pl.col("frequency") / 100).alias("std_frequency"),
        ])
    # Equal 1/3 weights on coverage, frequency, threshold-score (n_mod is omitted).
    df = df.with_columns((pl.col("std_coverage") * 0.3333 + 
                        pl.col("std_frequency") * 0.3333 + 
                        pl.col("std_score_nmod_thresholds") * 0.3333).alias("weighted_sum"))
    return df

combined_labs_df = compute_weighted_sums(combined_labs_df)


# ============================================================================
# STEP 8 — Modification-name harmonisation dictionary
# ----------------------------------------------------------------------------
# Different tools / ChEBI IDs / aliases all refer to the same modification. This
# map collapses every synonym (short code, ChEBI number, arrow-style substitution)
# onto one canonical short name so tables from different sources can be joined and
# plotted consistently. Grouped by chemistry for readability; "?" entries are
# ambiguous calls, "999xx" are internal placeholders for mods without a public
# ChEBI ID, and the "X->Y" entries are base substitutions (kept, not mods).
# ============================================================================
_MOD_MAP: dict[str, str] = {
    # --- Unmodified bases ---
    "A": "A",
    "16335": "A",                 # CHEBI:16335 — adenosine
    "G": "G",
    "16750": "G",                 # CHEBI:16750 — guanosine
    "C": "C",
    "17562": "C",                 # CHEBI:17562 — cytidine
    "U": "U",
    "16704": "U",                 # CHEBI:16704 — uridine
    "T": "T",

    # --- 2'-O-methyl ---
    "Am": "Am",
    "69426": "Am",                # CHEBI:69426 — 2'-O-methyladenosine
    "Cm": "Cm",
    "19228": "Cm",                # CHEBI:19228 — 2'-O-methylcytidine
    "Gm": "Gm",
    "19229": "Gm",                # CHEBI:19229 — 2'-O-methylguanosine
    "Um": "Um",
    "19227": "Um",                # CHEBI:19227 — 2'-O-methyluridine
    "Ym": "Ym",
    "m5Cm": "m5Cm",
    "184012": "m5Cm",             # CHEBI:184012 — 5,2'-O-dimethylcytidine
    "hm5Cm": "hm5Cm",
    "99997": "hm5Cm",             # placeholder (no public ChEBI ID)
    "f5Cm": "f5Cm",
    "99996": "f5Cm",              # placeholder (no public ChEBI ID)

    # --- Inosine ---
    "I": "I",
    "Ino": "I",
    "17596": "I",                 # CHEBI:17596 — inosine
    "A->G": "I",
    "m1I": "m1I",
    "19065": "m1I",               # CHEBI:19065 — 1-methylinosine

    # --- Pseudouridine ---
    "Y": "Y",
    "psU": "Y",
    "Psi": "Y",
    "17802": "Y",                 # CHEBI:17802 — pseudouridine
    "m1acp3Y": "m1acp3Y",
    "m1ap3U": "m1acp3Y",

    # --- Methyl-A ---
    "m1A": "m1A",
    "m1A ": "m1A",
    "16020": "m1A",               # CHEBI:16020 — 1-methyladenosine
    "m6A": "m6A",
    "a": "m6A",
    "x6A": "m6A",
    "21891": "m6A",               # CHEBI:21891 — N6-methyladenosine
    "m6,6A": "m6,6A",
    "m66A": "m6,6A",
    "28284": "m6,6A",             # CHEBI:28284 — N6,N6-dimethyladenosine
    "m7A": "m7A",

    # --- Methyl-C ---
    "m3C": "m3C",
    "m4C": "m4C",
    "m5C": "m5C",
    "m": "m5C",
    "20607": "m5C",               # CHEBI:20607 — 5-methylcytidine

    # --- Methyl-G ---
    "m1G": "m1G",
    "19062": "m1G",               # CHEBI:19062 — 1-methylguanosine
    "m2G": "m2G",
    "19702": "m2G",               # CHEBI:19702 — N2-methylguanosine
    "m2,2G": "m2,2G",
    "m22G": "m2,2G",
    "19289": "m2,2G",             # CHEBI:19289 — N2,N2-dimethylguanosine
    "m2,2,7G": "m2,2,7G",
    "143283": "m2,2,7G",          # CHEBI:143283 — N2,N2,7-trimethylguanosine
    "m6G": "m6G",
    "m7G": "m7G",
    "20794": "m7G",               # CHEBI:20794 — 7-methylguanosine

    # --- Methyl-U / dihydrouridine ---
    "m3U": "m3U",
    "89487": "m3U",               # CHEBI:89487 — 3-methyluridine
    "m5U": "m5U",
    "45996": "m5U",               # CHEBI:45996 — ribothymidine (5-methyluridine)
    "D": "D",
    "23774": "D",                 # CHEBI:23774 — dihydrouridine

    # --- Acetyl / acyl / threonyl ---
    "ac4C": "ac4C",
    "70989": "ac4C",              # CHEBI:70989 — N4-acetylcytidine
    "ac7G": "ac7G",
    "t6A": "t6A",
    "21440": "t6A",               # CHEBI:21440 — N6-threonylcarbamoyladenosine
    "m6t6A": "m6t6A",
    "133071": "m6t6A",            # CHEBI:133071 — N6-methyl-N6-threonylcarbamoyladenosine

    # --- C5-substituted U family ---
    "cm5U": "cm5U",
    "75654": "cm5U",              # CHEBI:75654 — 5-carboxymethyluridine
    "mcm5U": "mcm5U",
    "20598": "mcm5U",             # CHEBI:20598 — 5-methoxycarbonylmethyluridine
    "mcmo5U": "mcmo5U",
    "27241": "mcmo5U",            # CHEBI:27241 — 5-methoxycarbonylmethoxyuridine
    "mchm5U": "mchm5U",
    "99998": "mchm5U",            # placeholder (no public ChEBI ID)
    "0": "mchm5U",                # legacy placeholder
    "ncm5U": "ncm5U",
    "62005": "ncm5U",             # CHEBI:62005 — 5-carbamoylmethyluridine
    "ncm5Um": "ncm5Um",
    "99994": "ncm5Um",            # placeholder (no public ChEBI ID)
    "ncm5s2U": "ncm5s2U",
    "99995": "ncm5s2U",           # placeholder (no public ChEBI ID)

    # --- Hydroxymethyl / formyl C ---
    "hm5C": "hm5C",
    "191041": "hm5C",             # CHEBI:191041 — 5-hydroxymethylcytidine
    "f5C": "f5C",
    "234279": "f5C",              # CHEBI:234279 — 5-formylcytidine

    # --- Isopentenyl / thio A ---
    "i6A": "i6A",
    "62881": "i6A",               # CHEBI:62881 — N6-(Δ2-isopentenyl)adenosine
    "ms2i6A": "ms2i6A",
    "62875": "ms2i6A",            # CHEBI:62875 — 2-methylthio-N6-(Δ2-isopentenyl)adenosine
    "io6A": "io6A",
    "71693": "io6A",              # CHEBI:71693 — N6-(cis-hydroxyisopentenyl)adenosine
    "ms2io6A": "ms2io6A",
    "62879": "ms2io6A",           # CHEBI:62879 — 2-methylthio-N6-(cis-hydroxyisopentenyl)adenosine

    # --- 3-amino-3-carboxypropyl ---
    "acp3U": "acp3U",
    "acp3D": "acp3D",
    "71588": "acp3D",             # CHEBI:71588 — 3-(3-amino-3-carboxypropyl)dihydrouridine

    # --- Queuosine / 2-thio ---
    "Q": "Q",
    "60193": "Q",                 # CHEBI:60193 — queuosine
    "s2U": "s2U",
    "60731": "s2U",               # CHEBI:60731 — 2-thiouridine

    # --- Ambiguous / generic placeholders ---
    "mA?": "mA?",
    "99990": "mA?",               # placeholder
    "mC?": "mC?",
    "99991": "mC?",               # placeholder
    "mG?": "mG?",
    "99992": "mG?",               # placeholder
    "mU?": "mU?",
    "99993": "mU?",               # placeholder
    "U*": "U*",

    # --- Non-canonical / experimental placeholders ---
    "xp3Cm": "xp3Cm",
    "xp4U": "xp4U",
    "xp6G": "xp6G",
    "xp7G": "xp7G",
    "m2xp7G": "m2xp7G",

    # --- Base substitutions (not modifications) ---
    "G->T": "G->T",
    "G->A": "G->A",
    "G->C": "G->C",
    "C->A": "C->A",
    "C->T": "C->T",
    "C->G": "C->G",
    "A->C": "A->C",
    "A->T": "A->T",
    "T->C": "T->C",
    "T->G": "T->G",
    "T->A": "T->A",
}

def align_modification_names_polars(df: pl.DataFrame) -> pl.DataFrame:
    """
    Maps the 'name' column to standardised modification names using _MOD_MAP.

    Unmapped values are kept as-is. The 'name' column is cast to String before mapping.

    Parameters:
        df (pl.DataFrame): Input DataFrame containing a 'name' column with modification codes.

    Returns:
        pl.DataFrame: DataFrame with the 'name' column replaced by standardised modification names.
    """
    # replace(..., default=col) => leave anything not in the map untouched.
    return df.with_columns(
        pl.col("name")
        .cast(pl.String)
        .replace(_MOD_MAP, default=pl.col("name"))
        .alias("name")
    )


# %%
# ----------------------------------------------------------------------------
# Recover the composite score for the STEP 3 tables from the RGBA colour string.
# The upstream pipeline encoded the score into the itemRgb field; here the first
# numeric component after "rgba(" (the red channel, 0–255) is pulled out and
# rescaled to 0–1 to reconstruct weighted_sum for genome_df / ribosomal_df.
# ----------------------------------------------------------------------------
genome_df = genome_df.with_columns(
    ((pl.col("itemRgb").str.extract(r"rgba\((\d+)", 1).cast(pl.Float64)) / 255)
    .alias("weighted_sum")
)

ribosomal_df = ribosomal_df.with_columns(
    ((pl.col("itemRgb").str.extract(r"rgba\((\d+)", 1).cast(pl.Float64)) / 255)
    .alias("weighted_sum")
)

#%%
# ============================================================================
# STEP 9 — Load the native-vs-IVT Fisher table (for the unmodified control)
# ----------------------------------------------------------------------------
# The IVT sample is transcribed in vitro so it carries NO real modifications; it
# is the false-positive baseline. This large file is read lazily (scan_csv) and
# immediately narrowed to just the columns we need (position key + IVT freq/cov),
# since only those feed the enrichment test later.
# ============================================================================
ivt_fisher_bed = pl.scan_csv("/global/cfs/cdirs/m5243/final_bedRmods/ont/raw/polyA.native_vs_ivt_fisher.bed", 
                            has_header=False, 
                            comment_prefix="#",
                            new_columns=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 
                            'thickEnd', 'itemRgb', 'coverage', 'frequency', 'n_mod', 'count_canonical', 'count_other_mod', 
                            'count_delete', 'count_fail', 'count_diff', 'count_nocall', 'score_nmod_thresholds', 
                            'std_coverage', 'std_n_mod', 'std_score_nmod_thresholds', 'std_frequency', 
                            'pvalue_shuffle_axis_1', 'pvalue_adj_shuffle_axis_1', 'pvalue_shuffle_axis_0', 
                            'pvalue_ranks', 'pvalue_shuffle_axis_1_IVT', 'ivt_frequency', 'ivt_coverage', 'ivt_n_mod', 'pvalue', 'padj'
                            ],
                            separator="\t"
                            ).select(['chrom', 'chromStart', 
                                      'chromEnd', 'name', 'strand', 
                                      'ivt_frequency', 'ivt_coverage'
                                      ]).collect()
# %%
# ----------------------------------------------------------------------------
# Clean the IVT numeric columns: literal "NA" strings -> NaN, then cast to float
# so arithmetic works. (String "NA" would otherwise poison a numeric cast.)
# ----------------------------------------------------------------------------
ivt_fisher_bed = ivt_fisher_bed.with_columns(
    pl.col("ivt_frequency")
    .str.replace("^NA$", float("nan"))   # swap "NA" strings first
    .cast(pl.Float64)
    .alias("ivt_frequency")
)


ivt_fisher_bed = ivt_fisher_bed.with_columns(
    pl.col("ivt_coverage")
    .str.replace("^NA$", float("nan"))   # swap "NA" strings first
    .cast(pl.Float64)
    .alias("ivt_coverage")
)



#%% 
# ----------------------------------------------------------------------------
# Harmonise modification names on BOTH sides before joining, otherwise sites with
# different-but-equivalent name codes wouldn't match on the join key.
# ----------------------------------------------------------------------------
# Align modification names in combined_labs_df to match those in ivt_fisher_bed before the join
ivt_fisher_bed = align_modification_names_polars(ivt_fisher_bed)
combined_labs_df = align_modification_names_polars(combined_labs_df)

#%%
# ----------------------------------------------------------------------------
# Attach the IVT frequency/coverage to every native site via an INNER join on the
# full position key. Inner => only keep sites that also exist in the IVT table
# (i.e. sites we can actually test against the control).
# ----------------------------------------------------------------------------
combined_labs_df = combined_labs_df.join(ivt_fisher_bed, on=["chrom", "chromStart", "chromEnd", "name", "strand"], how="inner")

#%%
# ============================================================================
# STEP 10 — Native-vs-IVT enrichment test (vectorised hypergeometric / Fisher)
# ----------------------------------------------------------------------------
# For each site we build a 2x2 contingency table of modified vs canonical calls in
# native vs IVT, and ask: is the modification significantly MORE frequent in native
# than in the unmodified IVT control? A one-sided hypergeometric survival function
# gives that p-value for the whole column at once (no Python loop).
#   N = total calls across both conditions
#   K = total modified calls (native + IVT)         -> "successes" in the urn
#   n = total native calls (mod + canonical)        -> draws
#   k = native modified calls                        -> observed successes
#   p = P(X >= k) = hypergeom.sf(k-1, N, K, n)
# Then Benjamini-Hochberg FDR-correct across all valid sites.
# ============================================================================

# ── 1. Pull columns to NumPy in one shot ──────────────────────────────────────
n_mod        = combined_labs_df["n_mod"].to_numpy(allow_copy=True)
coverage     = combined_labs_df["coverage"].to_numpy(allow_copy=True)
ivt_freq     = combined_labs_df["ivt_frequency"].to_numpy(allow_copy=True)
ivt_coverage = combined_labs_df["ivt_coverage"].to_numpy(allow_copy=True)

# ── 2. Build valid mask (handles NaN in one vectorized op) ────────────────────
# Any site missing native or IVT numbers can't be tested -> excluded, p stays NaN.
valid = ~(np.isnan(n_mod) | np.isnan(coverage) |
          np.isnan(ivt_freq) | np.isnan(ivt_coverage))

pvalues = np.full(len(n_mod), np.nan, dtype=np.float64)

# ── 3. Fisher / hypergeometric SF — fully vectorized over valid rows ──────────
# Reconstruct raw counts: IVT modified count = freq(%) * coverage (rounded to int).
n_mod_nat = n_mod[valid].astype(np.int64)                          # native modified
n_can_nat = coverage[valid].astype(np.int64)                       # native "canonical"/total
n_mod_ivt = (ivt_freq[valid] * ivt_coverage[valid]).astype(np.int64)  # IVT modified
n_can_ivt = ivt_coverage[valid].astype(np.int64)                   # IVT total

N = n_mod_nat + n_can_nat + n_mod_ivt + n_can_ivt
K = n_mod_nat + n_mod_ivt
n = n_mod_nat + n_can_nat

pvalues[valid] = hypergeom.sf(n_mod_nat - 1, N, K, n)  # vectorized; no Python loop

# ── 4. BH correction — only on non-NaN p-values ───────────────────────────────
padj = np.full(len(pvalues), np.nan, dtype=np.float64)
valid_p = ~np.isnan(pvalues)
if valid_p.any():
    _, padj_valid, _, _ = multipletests(pvalues[valid_p], method="fdr_bh")
    padj[valid_p] = padj_valid

# ── 5. Write both columns back in one pass ────────────────────────────────────
combined_labs_df = combined_labs_df.with_columns([
    pl.Series("pvalue", pvalues),
    pl.Series("padj",   padj),
])



# %%
# ----------------------------------------------------------------------------
# Checkpoint: persist the fully-annotated + tested table to disk, then read it
# straight back. This makes the expensive STEP 6–10 work re-runnable without
# recomputing (comment out the write / keep the read to resume from here).
# ----------------------------------------------------------------------------
combined_labs_df.write_csv("combined_labs_with_ivt_fisher.csv", separator="\t")

#%%
combined_labs_df = pl.read_csv("combined_labs_with_ivt_fisher.csv", separator="\t", has_header=True)  

#%%
# Sanity check: how many rows before/after dropping chrY (again, defensively).
print(combined_labs_df.shape)
combined_labs_df = combined_labs_df.filter(pl.col("chrom") != "chrY")
print(combined_labs_df.shape)
combined_labs_df.shape
#%%
# ============================================================================
# STEP 11 — Final high-confidence site filter
# ----------------------------------------------------------------------------
# First force explicit dtypes (a CSV round-trip can leave columns as strings),
# then keep a site if EITHER of these "is it real?" criteria passes:
#   (a) significantly enriched vs IVT: padj <= 0.05 AND the IVT frequency is at
#       most half the native frequency (ivt_freq/freq <= 0.5), i.e. the control is
#       clearly cleaner than the native sample; OR
#   (b) very high composite confidence: weighted_sum >= 0.9.
# AND, regardless of which branch, require solid support:
#       coverage >= 30 reads AND frequency >= 3%.
# Result is materialised back to pandas for the PyRanges annotation that follows.
# ============================================================================
combined_labs_df = combined_labs_df.with_columns([
    pl.col("chrom").cast(pl.String),
    pl.col("chromStart").cast(pl.Int64),
    pl.col("chromEnd").cast(pl.Int64),
    pl.col("strand").cast(pl.String),
    pl.col("name").cast(pl.String),
    pl.col("lab").cast(pl.String),
    pl.col("pvalue").cast(pl.Float64),
    pl.col("padj").cast(pl.Float64),
    pl.col("coverage").cast(pl.Int64),
    pl.col("frequency").cast(pl.Float64),
    pl.col("ivt_frequency").cast(pl.Float64),
    pl.col("ivt_coverage").cast(pl.Float64),
    pl.col("weighted_sum").cast(pl.Float64),
]).filter(((pl.col("padj") <= 0.05) & (pl.col("ivt_frequency") / pl.col("frequency") <= 0.5)) | (pl.col("weighted_sum") >= 0.9)).filter((pl.col("coverage") >= 30) & (pl.col("frequency") >= 3)).to_pandas()

# %%
# ----------------------------------------------------------------------------
# Overlap the filtered high-confidence sites with gene bodies (same PyRanges
# join as STEP 4), to prepare per-feature annotation of the merged call-set.
# ----------------------------------------------------------------------------
combined_bed_ranges = pr.PyRanges(combined_labs_df.rename(columns={
    "chrom": "Chromosome",
    "chromStart": "Start",
    "chromEnd": "End",
    "strand": "Strand",
}))


gene_body_ranges = pr.PyRanges(gene_body_df.to_pandas().rename(columns={
    "seqname": "Chromosome",
    "start": "Start",
    "end": "End",
    "strand": "Strand",
}))


# Find overlaps
combined_overlaps = combined_bed_ranges.join(gene_body_ranges, how="left")
# Convert overlaps to DataFrame
combined_overlaps_df = combined_overlaps.df


# %%
# ----------------------------------------------------------------------------
# Identical feature-classification decision tree as STEP 5, now applied to the
# filtered merged call-set (see STEP 5 for the full logic explanation).
# ----------------------------------------------------------------------------
new_features = []
new_features_detailed = []
for feature,chrom,start,end,strand in tqdm(zip(combined_overlaps_df["feature"], 
                                          combined_overlaps_df["Chromosome"], 
                                          combined_overlaps_df["Start"], 
                                          combined_overlaps_df["End"], 
                                          combined_overlaps_df["Strand"]
                                          ), total=len(combined_overlaps_df)):
    if feature == "-1":
        if chrom.__contains__("18S"):
            new_features.append("18S")
            new_features_detailed.append("18S")
        elif chrom.__contains__("28S"):
            new_features.append("28S")
            new_features_detailed.append("28S")
        elif chrom.__contains__("5S"):
            new_features.append("5S")
            new_features_detailed.append("5S")
        elif chrom.__contains__("5.8S"):
            new_features.append("5.8S")
            new_features_detailed.append("5.8S")
        else:
            candidate_gene = gtf_df.filter(
                (pl.col("seqname") == chrom) &
                (pl.col("start") <= start) &
                (pl.col("end") >= end) &
                (pl.col("strand") == strand)
            )
            if candidate_gene.is_empty():
                new_features.append("intergenic")
                new_features_detailed.append("intergenic")
            else:
                if candidate_gene["gene_type"][-1] == "protein_coding" and candidate_gene["feature"][-1] in ["UTR","CDS"]:
                    if candidate_gene["feature"][-1] == "CDS":
                        new_features.append("CDS")
                        new_features_detailed.append("CDS")
                    elif candidate_gene["tags"][-1].__contains__("CCDS"):
                        new_features.append("CDS")
                        new_features_detailed.append("CDS")
                    elif candidate_gene["tags"][-1].__contains__("5_UTR"):
                        new_features.append("5'UTR")
                        new_features_detailed.append("5'UTR")
                    elif candidate_gene["tags"][-1].__contains__("3_UTR"):
                        new_features.append("3'UTR")
                        new_features_detailed.append("3'UTR")
                    else:
                        new_features.append("undefined intragenic")
                        new_features_detailed.append("undefined intragenic")
                elif candidate_gene["gene_type"][-1] == "protein_coding":
                    new_features.append("undefined intragenic")
                    new_features_detailed.append("undefined intragenic")
                else:
                    new_features.append("lncRNA/rRNA/snoRNA/other")
                    new_features_detailed.append(candidate_gene["gene_type"][-1])
    else:
        new_features.append(feature)
        new_features_detailed.append(feature)
combined_overlaps_df["feature"] = new_features
combined_overlaps_df["feature_detailed"] = new_features_detailed

# %%
# ----------------------------------------------------------------------------
# Prep the per-lab genome table (from STEP 5) as a comparison group and relabel
# every row's lab as "combined" so it can be plotted as one aggregate track.
# ----------------------------------------------------------------------------
genome_df = genome_df.to_pandas()
genome_df["lab"] = "combined"
print(genome_df.columns)

#%%
# ----------------------------------------------------------------------------
# Reduce both tables to a common set of columns and stack them, giving a single
# tidy frame of (position, modification, lab, feature, coverage) for downstream
# feature/modification summaries.
# ----------------------------------------------------------------------------
genome_df = genome_df[["Chromosome", "Start", "End", "Strand", "name", "lab", "feature", "coverage"]]
combined_overlaps_df = combined_overlaps_df[["Chromosome", "Start", "End", "Strand", "name", "lab", "feature", "coverage"]]


all_overlaps_df = pd.concat([genome_df, combined_overlaps_df], ignore_index=True)



#%%
# ============================================================================
# STEP 12 — Summarise modification-by-feature distribution (for the heatmap)
# ----------------------------------------------------------------------------
# Restrict to the "combined" track and compute, per modification:
#   * value_counts_mods     -> raw (feature, name) counts, descending
#   * rel                   -> proportion of each feature WITHIN each modification
#                              (normalised so each modification's features sum to 1)
#   * abs_                   -> the matching absolute counts
# rel + abs_ are merged so the heatmap can colour by proportion but annotate with
# the raw count. A "modification" alias column is added for clean plotting labels.
# ============================================================================
value_counts_mods = (
    all_overlaps_df[all_overlaps_df["lab"] =="combined"]
    .value_counts(["feature", "name"])
    .sort_values(ascending=False)
).reset_index()

# Relative proportions (normalized within each feature group)
# Relative proportions (normalized within each lab)
rel = (
    all_overlaps_df[all_overlaps_df["lab"] =="combined"]
    .groupby("name")["feature"]
    .value_counts(normalize=True)
    .rename("proportion")
)

# Absolute counts
abs_ = (
    all_overlaps_df[all_overlaps_df["lab"] =="combined"]
    .groupby("name")["feature"]
    .value_counts(normalize=False)
    .rename("count")
)

value_counts_mods_rel = (
    pd.concat([rel, abs_], axis=1)
    .reset_index()
    .sort_values(["name", "feature"], ascending=False)
    .reset_index(drop=True)
)

value_counts_mods["modification"] = value_counts_mods["name"]
value_counts_mods_rel["modification"] = value_counts_mods_rel["name"]

# %%
# ----------------------------------------------------------------------------
# Impose the canonical feature ordering on the summary so the heatmap columns come
# out in biological 5'->3' order rather than alphabetical.
# ----------------------------------------------------------------------------
category_order = ["5'UTR", "CDS", "3'UTR", "undefined intragenic","lncRNA/rRNA/snoRNA/other", "intergenic"]

value_counts_mods_rel["feature"] = pd.Categorical(
    value_counts_mods_rel["feature"],
    categories=category_order,
    ordered=True
)

# %%
# ============================================================================
# STEP 13 — figure 1: modification × feature heatmap
# ----------------------------------------------------------------------------
# Pivot the long summary into two matrices (rows = modification, cols = feature):
#   * prop_pivot  -> proportion, drives the colour (0..1, coolwarm)
#   * count_pivot -> raw count, printed as the cell annotation
# So colour answers "what fraction of this mod's sites are in this region?" while
# the number gives the absolute support behind that fraction.
# ============================================================================
fig, ax = plt.subplots(figsize=(16, 9))

prop_pivot = value_counts_mods_rel.pivot(index="modification", columns="feature", values="proportion")
count_pivot = value_counts_mods_rel.pivot(index="modification", columns="feature", values="count")

sns.heatmap(
    prop_pivot,
    annot=count_pivot,
    fmt=".0f",
    cmap="coolwarm",
    vmin=0,
    vmax=1,
    ax=ax,
)
ax.set_xlabel("Genomic Feature")
ax.set_ylabel("Modification type")
ax.set_xticklabels(ax.get_xticklabels(),rotation=45)


fig.savefig("Figure1_modification_types_on_features_heatmap.pdf", format="pdf")



# %%
# ============================================================================
# STEP 14 — PolyA tail-length data: load per-lab parquet tables
# ----------------------------------------------------------------------------
# Nanopore also yields a poly(A) tail-length estimate per read. These parquet
# tables (one per lab, plus a merged "combined") are read and tagged with a lab
# label for a distribution comparison. Parquet is used because it's compact/fast.
# ============================================================================
def load_lab(path, lab, add_threshold=True):
    df = pl.read_parquet(path)
    return df.with_columns(pl.lit(lab).alias("lab"))

base = "/global/cfs/cdirs/m5243/final_bedRmods/ont/taillength_tables/"

Lab1_df  = load_lab(base + "HRP_A_007_native_polyA_RNA_001.merged.taillengths.parquet", "HRP_A_007_1")
Lab2_df  = load_lab(base + "HRP_A_013_native_polyA_RNA_001.taillengths.parquet",        "HRP_A_013_1")
Lab3_df  = load_lab(base + "HRP_A_001_native_polyA_RNA_001.taillengths.parquet",        "HRP_A_001_1")
Lab4_df  = load_lab(base + "HRP_A_002_native_polyA_RNA_001.taillengths.parquet",        "HRP_A_002_1")
Lab15_df = load_lab(base + "HRP_A_003_native_polyA_RNA_001.taillengths.parquet",        "HRP_A_003_1")
Lab5_df  = load_lab(base + "HRP_A_006_native_polyA_RNA_004.taillengths.parquet",        "HRP_A_006_4")
Lab6_df  = load_lab(base + "HRP_A_005_native_polyA_RNA_001.taillengths.parquet",        "HRP_A_005_1")
Lab7_df  = load_lab(base + "HRP_A_014_native_polyA_RNA_001.taillengths.parquet",        "HRP_A_014_1")
Lab8_df  = load_lab(base + "HRP_A_009_native_polyA_RNA_001.taillengths.parquet",        "HRP_A_009_1")
Lab9_df  = load_lab(base + "HRP_A_009_native_polyA_RNA_002.taillengths.parquet",        "HRP_A_009_2")
Lab10_df = load_lab(base + "HRP_A_010_native_polyA_RNA_001.taillengths.parquet",        "HRP_A_010_1")
Lab11_df = load_lab(base + "HRP_A_010_native_polyA_RNA_002.taillengths.parquet",        "HRP_A_010_2")
Lab12_df = load_lab(base + "HRP_A_011_native_polyA_RNA_001.taillengths.parquet",        "HRP_A_011_1")
Lab13_df = load_lab(base + "HRP_A_011_native_polyA_RNA_002.taillengths.parquet",        "HRP_A_011_2")
Lab14_df = load_lab(base + "HRP_A_011_native_polyA_RNA_003.taillengths.parquet",        "HRP_A_011_3")

combined_df = load_lab(base + "merged_dataset.taillengths.parquet",                     "combined")
# %%
# ----------------------------------------------------------------------------
# Stack all tail-length tables (individual labs + merged) into one long frame,
# then checkpoint to a single parquet and read it back (resume point).
# ----------------------------------------------------------------------------
frames = [Lab1_df, Lab2_df, Lab3_df, Lab4_df, Lab5_df, Lab6_df, Lab7_df, Lab8_df, Lab9_df, Lab10_df, Lab11_df, Lab12_df, Lab13_df, Lab14_df, Lab15_df, combined_df]
concat_taillengths = pl.concat(frames, how="vertical")

#%%
concat_taillengths.write_parquet("/global/cfs/cdirs/m5243/final_bedRmods/ont/taillength_tables/single_datasets_and_merged_taillengths.parquet")
#%%

concat_taillengths = pl.read_parquet("/global/cfs/cdirs/m5243/final_bedRmods/ont/taillength_tables/single_datasets_and_merged_taillengths.parquet")


#%%
# ----------------------------------------------------------------------------
# Ordering for the x-axis: put the aggregate "combined" violin first, then the
# individual HRP labs in sorted order so the plot reads left-to-right sensibly.
# ----------------------------------------------------------------------------
# Build the desired category order: "combined" first, HRPs sorted after
hrp_labels = sorted(l for l in concat_taillengths["lab"].unique() if l != "combined")
lab_order = ["combined"] + hrp_labels

#%%
# ============================================================================
# STEP 15 — Tail-length distribution violin plot (per lab)
# ----------------------------------------------------------------------------
# One violin per lab showing the full poly(A) tail-length density. inner="quartile"
# draws median + quartile lines; cut=0 stops the kernel density from extending past
# the observed data. hue duplicates x purely for colour, so the legend is dropped.
# Saved as SVG (vector) for publication.
# ============================================================================
fig, ax = plt.subplots(nrows=1, ncols=1)

sns.violinplot(
    data=concat_taillengths,
    x="lab", y="taillength", hue="lab",
    order=lab_order, hue_order=lab_order,
    inner="quartile",        # show median + quartile lines inside each violin
    cut=0,                   # don't extend density past observed data range
    linewidth=1.2,
    saturation=0.85,
    legend=False,            # hue duplicates the x axis, so drop the legend
    ax=ax,
)

# Labels and title
ax.set_xlabel("Lab", fontsize=12, labelpad=8)
ax.set_ylabel("PolyA taillength", fontsize=12, labelpad=8)   # adjust unit
ax.set_title("Tail length distribution by lab", fontsize=13, pad=12)

# Clean up the frame
sns.despine(ax=ax, trim=True)
ax.grid(axis="y", linestyle="--", alpha=0.3)
ax.tick_params(labelsize=10)

fig.tight_layout()
fig.savefig("/global/cfs/cdirs/m5243/final_bedRmods/ont/taillength_tables/taillengths_violin.svg", bbox_inches="tight")   # vector for publication
plt.show()

# %%
# ============================================================================
# STEP 16 — (Re-)imports for the PCA section + PCA helper
# ----------------------------------------------------------------------------
# Because this is a cell-based (#%%) script, this block re-imports what the PCA
# work needs so the section can run independently. It then defines pca_plot(),
# which turns a list of BED tables into a samples × positions matrix and runs PCA
# to compare whole datasets against each other on a chosen per-site metric.
# ============================================================================
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from upsetplot import UpSet,from_contents,plot
from matplotlib import pyplot
import re
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.metrics import jaccard_score
from scipy.stats import norm
from tqdm import tqdm
import polars as pl
import math
from statsmodels.stats.multitest import multipletests
import io
import pyarrow
from joblib import Parallel, delayed


def pca_plot(
    bed_df_list: list[pd.DataFrame],
    metric: str,
    sample_labels: list[str] | None = None,
    n_components: int = 2,
    figsize: tuple = (8, 6),
) -> tuple[PCA, pd.DataFrame]:
    """
    Perform and plot PCA across multiple BED DataFrames using a chosen numeric metric.

    Each DataFrame becomes one sample (row) in the PCA matrix.
    Genomic positions are the features (columns); the chosen metric supplies the values.
    Positions absent in a given sample are filled with 0.

    Parameters
    ----------
    bed_df_list   : List of BED DataFrames. Each must contain the columns
                    chrom, chromStart, chromEnd, strand, name, and `metric`.
    metric        : Column name of the numeric value to compare (e.g. "score",
                    "coverage", "signal", "pValue").
    sample_labels : Display names for each DataFrame. Defaults to
                    ["Sample_0", "Sample_1", ...].
    n_components  : Number of principal components to compute (≥ 2).
    figsize       : Figure size passed to matplotlib.

    Returns
    -------
    pca_model     : Fitted sklearn PCA object.
    scores_df     : DataFrame of PC scores (rows = samples).
    """
    # Default sample names if none supplied.
    if sample_labels is None:
        sample_labels = [f"Sample_{i}" for i in range(len(bed_df_list))]

    if len(sample_labels) != len(bed_df_list):
        raise ValueError("`sample_labels` length must match `bed_df_list` length.")

    # ------------------------------------------------------------------ #
    # 1. Build a position token for every DataFrame and pivot to wide form #
    #    Each site becomes a unique column key; each dataset becomes a row. #
    # ------------------------------------------------------------------ #
    wide_frames = []

    for label, df in tqdm(zip(sample_labels, bed_df_list)):
        # Guard: every input must carry the key columns + the requested metric.
        required = {"chrom", "chromStart", "chromEnd", "strand", "name", metric}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"[{label}] Missing columns: {missing}")

        df = df.copy()
        # A single string uniquely identifying a modified position (incl. mod name).
        df["position_token"] = [
            f"{chrom}:{start}:{end}:{strand}:{name}"
            for chrom, start, end, strand, name in zip(
                df["chrom"],
                df["chromStart"],
                df["chromEnd"],
                df["strand"],
                df["name"],
            )
        ]

        # One row per sample → columns are position tokens, values are the metric.
        # .mean() collapses any duplicate rows for the same site within a dataset.
        wide = (
            df.groupby("position_token")[metric]
            .mean()                        # collapse duplicate positions
            .rename(label)
            .to_frame()
            .T
        )
        wide_frames.append(wide)

    # Align all samples on the union of positions; fill missing with 0
    # (a site not called in a sample contributes 0 for that feature).
    combined = pd.concat(wide_frames, axis=0).fillna(0)
    combined.index = sample_labels

    if combined.shape[0] < 2:
        raise ValueError("PCA requires at least 2 samples.")

    # Can't ask for more components than samples or features.
    n_components = min(n_components, combined.shape[0], combined.shape[1])

    # ------------------------------------------------------------------ #
    # 2. Scale and run PCA                                                #
    #    Standardise each position (feature) before PCA so high-coverage  #
    #    sites don't dominate purely by magnitude.                        #
    # ------------------------------------------------------------------ #
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(combined.values)

    pca_model = PCA(n_components=n_components)
    scores = pca_model.fit_transform(X_scaled)

    pc_cols = [f"PC{i + 1}" for i in range(n_components)]
    scores_df = pd.DataFrame(scores, index=sample_labels, columns=pc_cols)

    # ------------------------------------------------------------------ #
    # 3. Plot PC1 vs PC2, one labelled point per dataset.                #
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(figsize=figsize)

    for label, row in scores_df.iterrows():
        ax.scatter(row["PC1"], row["PC2"], s=80, zorder=3)
        ax.annotate(label, (row["PC1"], row["PC2"]),
                    textcoords="offset points", xytext=(6, 4), fontsize=9)

    # Axis labels report the % variance each PC explains.
    var = pca_model.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({var[0]:.1%} variance explained)")
    ax.set_ylabel(f"PC2 ({var[1]:.1%} variance explained)")
    ax.set_title(f"PCA — metric: modification level")
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    return fig, pca_model, scores_df

#%%
# ============================================================================
# STEP 17 — Figure 1S: PCA of the 15 polyA-RNA labs (metric = mod frequency)
# ----------------------------------------------------------------------------
# Reuse the per-lab polars frames from STEP 6 (still in memory), pull each lab's
# name, and run pca_plot on modification frequency. This shows how similar the
# labs are to each other in their per-site mod-frequency profiles.
# ============================================================================
bed_files_list = [Lab1_df, Lab2_df, Lab3_df, Lab4_df, Lab5_df, Lab6_df, Lab7_df, Lab8_df, Lab9_df, Lab10_df, Lab11_df, Lab12_df, Lab13_df, Lab14_df, Lab15_df]
names_list = [i["lab"][0] for i in bed_files_list]

print(len(bed_files_list))
print(len(names_list))

# n_components=15 = number of samples, so all PCs are available (plot uses PC1/PC2).
fig, pca_model, scores_df = pca_plot(bed_df_list=[i.to_pandas() for i in bed_files_list],sample_labels=names_list,metric="frequency", n_components=15)

fig.savefig("/global/cfs/cdirs/m5243/final_bedRmods/ont/Figure1S_polyARNA_PCA_mod_frequency.pdf", bbox_inches="tight") 


# %%
# ============================================================================
# STEP 18 — Figure 1S: PCA of the rRNA datasets
# ----------------------------------------------------------------------------
# Load 7 rRNA per-site BEDs (transcriptome-aligned; note the schema uses
# count_modified rather than n_mod), tag each with its lab, and run the same
# frequency PCA to compare rRNA modification profiles across labs/replicates.
# ============================================================================
rRNA1_df = pd.read_csv("/home/stefan/Synology/Data_jamboree/rRNA_session/ONT/Gerber_native_rebasecalled/native_rRNA_SUP_Gerber_native_rRNA_251118_r1_v1.merged.transcriptome.aligned.sorted.transcriptome.bed", 
                       sep="\t", comment="#", header=None, 
                       names=["chrom","chromStart","chromEnd","name",
                              "score","strand","thickStart","thickEnd","color",
                              "coverage","frequency","count_modified","count_canonical",
                              "count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
rRNA1_df["lab"] = "HRP_A_006_1"


rRNA2_df = pd.read_csv("/home/stefan/Synology/Data_jamboree/rRNA_session/ONT/Novoa_rRNA_Genterga_rebasecalled/native_rRNA_Novoa_Eppendorf_Rep1.merged.transcriptome.aligned.sorted.transcriptome.bed", 
                       sep="\t", comment="#", header=None, 
                       names=["chrom","chromStart","chromEnd","name",
                              "score","strand","thickStart","thickEnd","color",
                              "coverage","frequency","count_modified","count_canonical",
                              "count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
rRNA2_df["lab"] = "HRP_A_018_1"


rRNA3_df = pd.read_csv("/home/stefan/Synology/Data_jamboree/rRNA_session/ONT/Novoa_rRNA_Genterga_rebasecalled/native_rRNA_Novoa_Eppendorf_Rep2.merged.transcriptome.aligned.sorted.transcriptome.bed", 
                       sep="\t", comment="#", header=None, 
                       names=["chrom","chromStart","chromEnd","name",
                              "score","strand","thickStart","thickEnd","color",
                              "coverage","frequency","count_modified","count_canonical",
                              "count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
rRNA3_df["lab"] = "HRP_A_018_2"

rRNA4_df = pd.read_csv("/home/stefan/Synology/Data_jamboree/rRNA_session/ONT/Novoa_rRNA_Genterga_rebasecalled/native_rRNA_Novoa_Eppendorf_Rep3.merged.transcriptome.aligned.sorted.transcriptome.bed", 
                       sep="\t", comment="#", header=None, 
                       names=["chrom","chromStart","chromEnd","name",
                              "score","strand","thickStart","thickEnd","color",
                              "coverage","frequency","count_modified","count_canonical",
                              "count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
rRNA4_df["lab"] = "HRP_A_018_3"



rRNA5_df = pd.read_csv("/home/stefan/Synology/Data_jamboree/rRNA_session/ONT/Novoa_rRNA_MaximaH_SSII_rebasecalled/native_rRNA_Novoa_rRNA_native_rep1.merged.transcriptome.aligned.sorted.transcriptome.bed", 
                       sep="\t", comment="#", header=None, 
                       names=["chrom","chromStart","chromEnd","name",
                              "score","strand","thickStart","thickEnd","color",
                              "coverage","frequency","count_modified","count_canonical",
                              "count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
rRNA5_df["lab"] = "HRP_A_017_1"


rRNA6_df = pd.read_csv("/home/stefan/Synology/Data_jamboree/rRNA_session/ONT/Novoa_rRNA_MaximaH_SSII_rebasecalled/native_rRNA_Novoa_rRNA_native_rep2.merged.transcriptome.aligned.sorted.transcriptome.bed", 
                       sep="\t", comment="#", header=None, 
                       names=["chrom","chromStart","chromEnd","name",
                              "score","strand","thickStart","thickEnd","color",
                              "coverage","frequency","count_modified","count_canonical",
                              "count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
rRNA6_df["lab"] = "HRP_A_017_2"


rRNA7_df = pd.read_csv("/home/stefan/Synology/Data_jamboree/rRNA_session/ONT/Novoa_rRNA_MaximaH_SSII_rebasecalled/native_rRNA_Novoa_rRNA_native_rep3.merged.transcriptome.aligned.sorted.transcriptome.bed", 
                       sep="\t", comment="#", header=None, 
                       names=["chrom","chromStart","chromEnd","name",
                              "score","strand","thickStart","thickEnd","color",
                              "coverage","frequency","count_modified","count_canonical",
                              "count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
rRNA7_df["lab"] = "HRP_A_017_3"

bed_files_list = [rRNA1_df,rRNA2_df,rRNA3_df,rRNA4_df,rRNA5_df,rRNA6_df,rRNA7_df]
names_list = [i["lab"][0] for i in bed_files_list]

fig, pca_model, scores_df = pca_plot(bed_df_list=[i for i in bed_files_list],sample_labels=names_list,metric="frequency", n_components=8)

fig.savefig("/global/cfs/cdirs/m5243/final_bedRmods/ont/Figure1S_rRNA_PCA_mod_frequency.pdf", bbox_inches="tight") 


#%%
# ============================================================================
# STEP 19 — Figure 1S: PCA of the tRNA datasets
# ----------------------------------------------------------------------------
# load_tRNA reads a threshold-0.98 tRNA BED, drops zero-frequency sites (only
# keep positions with a detected modification), and tags the lab. All 15 tRNA
# files are loaded into a dict keyed by lab, then fed to the same frequency PCA.
# NOTE: unlike STEP 17/18, bed_files_list and names_list are wrapped in an extra
# list ([ ... ]) here, so pca_plot receives a single-element list-of-lists. If the
# PCA errors or shows one point, that nesting is the likely cause. Left as-is.
# ============================================================================
import os
import pandas as pd

def load_tRNA(path, lab):
    df = pd.read_csv(path,
                    sep="\t", comment="#", header=None,
                    names=["chrom","chromStart","chromEnd","name",
                            "score","strand","thickStart","thickEnd","color",
                            "coverage","frequency","count_modified","count_canonical",
                            "count_other_mod","count_delete","count_fail","count_diff","count_nocall"])
    df = df[df["frequency"] > 0]   # keep only positions with a called modification
    df["lab"] = lab
    return df

folder = "/home/stefan/Synology/Data_jamboree/tRNA_sessions/single_tRNA"

# (filename, lab label) pairs: three donors x three replicates plus three
# single-run contributor datasets.
files_labels = [
    ("HRP_A_020_1_native_tRNA_001_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_020_1_1"),
    ("HRP_A_020_1_native_tRNA_002_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_020_1_2"),
    ("HRP_A_020_1_native_tRNA_003_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_020_1_3"),
    ("HRP_A_020_2_native_tRNA_001_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_020_2_1"),
    ("HRP_A_020_2_native_tRNA_002_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_020_2_2"),
    ("HRP_A_020_2_native_tRNA_003_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_020_2_3"),
    ("HRP_A_020_3_native_tRNA_001_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_020_3_1"),
    ("HRP_A_020_3_native_tRNA_002_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_020_3_2"),
    ("HRP_A_020_3_native_tRNA_003_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_020_3_3"),
    ("HRP_A_021_1_native_tRNA_001_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_021_1_1"),
    ("HRP_A_021_1_native_tRNA_002_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_021_1_2"),
    ("HRP_A_021_1_native_tRNA_003_SUP.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_021_1_3"),
    ("native_tRNA_SUP_ADAMCZYK_M_NATIVE_tRNA.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_004_1"),
    ("native_tRNA_SUP_GERBER_ONT_native_tRNA.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_006_1"),
    ("native_tRNA_SUP_Soares_ONT_native_tRNA.merged.transcriptome.aligned.sorted.threshold_0.98.tRNA.bed", "HRP_A_012_1"),
]

# Load every file into {lab: dataframe}.
dfs = {lab: load_tRNA(os.path.join(folder, fname), lab) for fname, lab in files_labels}


bed_files_list =  [list(dfs.values())]   # NOTE: extra [] nesting (see step header)
names_list = [list(dfs.keys())]

fig, pca_model, scores_df = pca_plot(bed_df_list=bed_files_list,sample_labels=names_list,metric="frequency", n_components=15)

fig.savefig("/global/cfs/cdirs/m5243/final_bedRmods/ont/Figure1S_tRNA_PCA_mod_frequency.pdf", bbox_inches="tight") 
