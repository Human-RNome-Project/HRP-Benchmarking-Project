import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import numpy as np
from logomaker import Logo, alignment_to_matrix

matplotlib.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype":     42,   # TrueType — required by Nature
    "ps.fonttype":      42,
    "svg.fonttype":     "none",
    "font.size":        7,    # max body text
    "axes.titlesize":   7,
    "axes.labelsize":   7,
    "xtick.labelsize":  6,
    "ytick.labelsize":  6,
    "axes.linewidth":   0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

MM = 1 / 25.4   # millimetres → inches, for Nature column widths
W_SINGLE = 89  * MM   # single-column width
W_DOUBLE = 183 * MM   # double-column width
H_MAX    = 170 * MM   # maximum figure height

PT_PANEL_LETTER = 8   # panel labels: 8-pt bold lowercase (Nature spec)
PT_BODY = 7           # all other text ≤ 7 pt
PT_TICK = 6

# ── Configuration ──────────────────────────────────────────────────────────────

HERE   = os.path.dirname(os.path.abspath(__file__))

GENOME  = os.path.join(HERE, "hg38.fa")
BED_DIR = os.path.expanduser("~/RNome/FINAL_BEDRMOD_FINAL")


ANALYSES = {
    "ONT":      os.path.join(BED_DIR, "ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed"),
    "Illumina": os.path.join(BED_DIR, "Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed"),
}

# Width of the sequence context window (odd number, centred on the modified base)
KMER = 9
FLANK = KMER // 2
# x-axis tick labels: -FLANK … 0 … +FLANK (the modified base sits at 0)
POS_LABELS = [("0" if d == 0 else f"{d:+d}") for d in range(-FLANK, FLANK + 1)]

# Biotype assigned from the reference-contig name. The custom rRNA/tRNA contigs
# (cytoplasmic + mitochondrial) carry the hs_ prefix; everything else is mRNA.
def _is_rRNA(c): return c.startswith("hs_rRNA") or c.startswith("hs_mt-rRNA")
def _is_tRNA(c): return c.startswith("hs_tRNA") or c.startswith("hs_mttRNA")
BIOTYPES = {
    "mRNA": lambda c: not (_is_rRNA(c) or _is_tRNA(c)),
    "rRNA": _is_rRNA,
    "tRNA": _is_tRNA,
}

MODIFICATION_COLORS = {
    # A — m6A family (deep crimson → light blush)
    "Am": "#D44F3E", "i6A": "#C0392B", "m1A": "#A52020", "m6A": "#721817",
    "m6,6A": "#F5C4B8", "m6t6A": "#F0A898", "ms2i6A": "#8B1A1A",
    "mxA": "#E8907A", "t6A": "#E06050",
    # C — m5C family (near-black navy → pale sky)
    "ac4C": "#6AAED6", "Cm": "#0D3B6E", "f5C": "#B8D9F0", "f5Cm": "#D4E4F0",
    "hm5C": "#2A7FC4", "hm5Cm": "#559CC0", "m3C": "#B8D9F0", "m5C": "#001427",
    "m5Cm": "#3A8FD4", "mxC": "#1E6EB5",
    # G — Inosine family (dark forest → pale mint)
    "Gm": "#74B354", "I": "#2D6E1E", "m1G": "#3A7A28", "m1I": "#BAE09E",
    "m2,2,7G": "#A8D48A", "m2,2G": "#8DC46A", "m2G": "#5C9840", "m7G": "#D4EEC4",
    "mxG": "#4A8532", "Q": "#ACC49B",
    # U — Psi family (deep amber → pale gold)
    "acp3D": "#FDE9B8", "acp3U": "#A86A00", "cm5U": "#F6C830", "D": "#F9D47E",
    "m5U": "#FBDD80", "mchm5U": "#E8950A", "mcm5U": "#FAD55A", "mcmo5U": "#FBD96A",
    "mxU": "#F5BE45", "ncm5U": "#FCE08A", "s2U": "#E0A858", "U*": "#FBE8AE",
    "Um": "#C47A02", "Y": "#F0A202",
}


NT_COLORS = {
    "A": MODIFICATION_COLORS["m6A"],   # #721817  A — m6A family
    "C": MODIFICATION_COLORS["m5C"],   # #001427  C — m5C family
    "G": MODIFICATION_COLORS["I"],     # #2D6E1E  G — Inosine family
    "T": MODIFICATION_COLORS["Y"],     # #F0A202  U — Psi family
    "U": MODIFICATION_COLORS["Y"],
}

_MODS = {
    "m6A": ("m6A",  "m6A", "A"),
    "Am":  ("Am",   "Am",  "A"),
    "I":   ("I",    "I",   "G"),
    "Y":   ("Y",    "Y",   "U"),
    "Um":  ("Um",   "Um",  "U"),
    "Gm":  ("Gm",   "Gm",  "G"),
    "Cm":  ("Cm",   "Cm",  "C"),
    "m5C": ("m5C",  "m5C", "C"),
}

# Order: (1) canonical base A,C,G,U, then (2) lexicographically by name.
_BASE_RANK = {"A": 0, "C": 1, "G": 2, "U": 3}
MOD_ORDER = sorted(_MODS, key=lambda k: (_BASE_RANK[_MODS[k][2]], _MODS[k][1]))

# Mapping kept for the plotting loops: name → (label, modification color)
MOD_CONFIG = {
    k: (_MODS[k][0], MODIFICATION_COLORS[_MODS[k][1]]) for k in MOD_ORDER
}

# ── Helpers ────────────────────────────────────────────────────────────────────

def get_sequence(fasta, chrom, pos, strand):
    """Extract KMER-wide context centred on pos (0-based half-open). Returns DNA string."""
    start = max(0, pos - FLANK)
    end   = pos + FLANK + 1          # KMER bases total
    seq   = fasta.fetch(chrom, start, end).upper()
    if strand == "-":
        seq = seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
    return seq

# ── Load bedfile ───────────────────────────────────────────────────────────────

BED_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "mod_start", "mod_end", "color", "n_valid_cov", "frac_mod",
]


def info_matrix_RNA(kmers):
    """Information-content (bits) matrix with the DNA T column relabelled to RNA U.

    Per-position height = 2 - H (Shannon entropy), so a position only appears
    large when the sequence is genuinely conserved (max 2 bits for 4 bases).
    Uninformative flanks collapse toward 0, unlike a probability logo.
    """
    info_matrix = alignment_to_matrix(kmers, to_type="information")
    # RNA: thymine → uracil
    info_matrix = info_matrix.rename(columns={"T": "U"})
    return info_matrix


def load_kmers(bed_path, fasta):
    """Read a bedRMod file and return a frame of valid KMER contexts."""
    df = pd.read_csv(
        bed_path, sep="\t", comment="#", header=None,
        usecols=range(11), names=BED_COLS,
    )
    df = df[df["name"].isin(MOD_CONFIG)]
    print(f"  {len(df):,} sites across {df['name'].nunique()} modification types")

    fasta_chroms = set(fasta.references)
    n_before = len(df)
    df = df[df["chrom"].isin(fasta_chroms)]
    dropped = n_before - len(df)
    if dropped:
        print(f"  {dropped:,} sites on contigs absent from the reference (dropped)")

    print(f"  Extracting {KMER}-mer context …")
    df["kmer"] = df.apply(
        lambda r: get_sequence(fasta, r["chrom"], r["start"], r["strand"]), axis=1
    )
    df = df[df["kmer"].str.len() == KMER]   # drop edge cases near contig boundaries
    return df


def report_and_drop_N(df, outdir):
    """Report k-mers with ambiguous N bases, save the breakdown, drop them."""
    n_mask = df["kmer"].str.contains("N")
    if n_mask.any():
        print(f"  Sites excluded due to N in {KMER}-mer: {int(n_mask.sum()):,}")
        stats = (
            df[n_mask].groupby(["chrom", "name"]).size()
            .reset_index(name="n_sites").sort_values("n_sites", ascending=False)
        )
        stats.to_csv(f"{outdir}N_kmer_failures.tsv", sep="\t", index=False)
    df = df[~n_mask]
    print(f"  {len(df):,} sites retained after N-filter")
    return df


# ── Per-analysis pipeline ───────────────────────────────────────────────────────

N_COLS = 2
PANEL_LETTERS = "abcdefghijklmnop"


def make_logo_axis(ax, info_matrix):
    """Draw an information-content (bits) logo with the shared Nature-styled axis."""
    Logo(info_matrix, ax=ax, color_scheme=NT_COLORS)
    ax.set_xticks(range(KMER))
    ax.set_xticklabels(POS_LABELS, fontsize=PT_TICK)
    ax.set_ylim(0, 2)
    ax.set_yticks([0, 1, 2])
    ax.tick_params(labelsize=PT_TICK)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def run_analysis(name, bed_path, fasta):
    """Full motif analysis for one platform → HERE/<name>/."""
    outdir = os.path.join(HERE, name) + os.sep
    os.makedirs(outdir, exist_ok=True)
    print(f"\n========== {name} ==========")
    print(f"Loading {os.path.basename(bed_path)} …")

    df = load_kmers(bed_path, fasta)
    df = report_and_drop_N(df, outdir)

    # results[biotype][mod_name] = (info_matrix, n)  — reused for the panel
    results = {}
    for biotype, chrom_filter in BIOTYPES.items():
        svg_dir = f"{outdir}{biotype}/svg/"
        pdf_dir = f"{outdir}{biotype}/pdf/"
        os.makedirs(svg_dir, exist_ok=True)
        os.makedirs(pdf_dir, exist_ok=True)
        df_bio = df[df["chrom"].apply(chrom_filter)]
        print(f"\n── {name} / {biotype}  ({len(df_bio):,} sites) ──")
        results[biotype] = {}

        for mod_name in MOD_ORDER:
            label = MOD_CONFIG[mod_name][0]
            sub = df_bio[df_bio["name"] == mod_name]["kmer"]
            n   = len(sub)
            if n == 0:
                print(f"  Skipping {label}: no sites")
                continue
            print(f"  Plotting {label} (n={n:,}) …")

            info_matrix = info_matrix_RNA(sub.tolist())
            results[biotype][mod_name] = (info_matrix, n)

            fig, ax = plt.subplots(figsize=(W_SINGLE, 38 * MM))
            make_logo_axis(ax, info_matrix)
            ax.set_xlabel("Position relative to modification", fontsize=PT_BODY)
            ax.set_ylabel("Information (bits)", fontsize=PT_BODY)
            ax.set_title(f"{label}  [{name} {biotype}]  (n={n:,})",
                         fontsize=PT_BODY, color="black")
            fig.tight_layout()
            stem = f"motif_{label}"
            fig.savefig(f"{svg_dir}{stem}.svg", transparent=True)
            fig.savefig(f"{pdf_dir}{stem}.pdf", transparent=True)
            plt.close(fig)

    # ── one panel per biotype: 8 mods, 4 rows (A,C,G,U) × 2 cols ──────────────
    n_rows = -(-len(MOD_ORDER) // N_COLS)   # ceil
    fig_h = min(H_MAX, 36 * MM * n_rows)
    for biotype in BIOTYPES:
        print(f"Building {name} / {biotype} panel …")
        fig, axes = plt.subplots(n_rows, N_COLS, figsize=(W_DOUBLE, fig_h),
                                 squeeze=False)
        for idx, mod_name in enumerate(MOD_ORDER):
            r, c = divmod(idx, N_COLS)
            ax = axes[r][c]
            label = MOD_CONFIG[mod_name][0]
            entry = results[biotype].get(mod_name)

            if entry is None:
                ax.text(0.5, 0.5, "no sites", ha="center", va="center",
                        fontsize=PT_BODY, color="0.4", transform=ax.transAxes)
                ax.set_xticks([]); ax.set_yticks([])
                for s in ax.spines.values():
                    s.set_visible(False)
                ax.set_title(label, fontsize=PT_BODY, color="black")
            else:
                info_matrix, n = entry
                make_logo_axis(ax, info_matrix)
                ax.set_title(f"{label}  (n={n:,})", fontsize=PT_BODY, color="black")

            if c == 0:
                ax.set_ylabel("Information (bits)", fontsize=PT_BODY)
            if r == n_rows - 1:
                ax.set_xlabel("Position relative to modification", fontsize=PT_BODY)

        fig.tight_layout(h_pad=1.4, w_pad=2.0)
        stem = f"motif_panel_{biotype}"
        # vector formats only (Nature does not accept raster for main figures)
        fig.savefig(f"{outdir}{stem}.svg", transparent=True)
        fig.savefig(f"{outdir}{stem}.pdf", transparent=True)
        plt.close(fig)
        print(f"  saved {name}/{stem}.svg / .pdf")


# ── Run both analyses, kept entirely apart ──────────────────────────────────────

print("Opening reference genome …")
fasta = pysam.FastaFile(GENOME)
for name, bed_path in ANALYSES.items():
    run_analysis(name, bed_path, fasta)
fasta.close()

print("\nDone.")
