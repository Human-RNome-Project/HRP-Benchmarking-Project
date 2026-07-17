#!/usr/bin/env python3
"""
calculate_native_mismatch_rates_fast.py

Streamlined native DRS mismatch rate calculation with two execution modes:

  --mode chrom  (default)
      Process a single per-chromosome pysamstats TSV.
      Outputs per_position_mismatch.tsv only — no aggregation.
      Run in parallel across all chroms × all samples in WDL scatter.

  --mode merge
      Takes a list of per-chrom per_position_mismatch.tsv files via
      --position-tsvs, concatenates them, and runs full aggregation:
      error rates, substitution spectrum, modkit stratification,
      native-IVT delta, masking report, save + print.

WDL execution model:
  scatter(sample) {
    scatter(chrom) {
      call ChromTask   → per_position_mismatch.tsv per chrom
    }
    call MergeTask     → all summary outputs
  }

Requirements:
  bedtools >= 2.30  (chrom mode only)
  pandas, numpy, scipy

Pre-build shared masks once:
    # variants.bed
    bcftools query -r chr1,...,chrX -R giab.bed \\
        -f '%CHROM\\t%POS0\\t%END\\n' giab.vcf.gz \\
        | sort -k1,1 -k2,2n | bedtools merge -i stdin > shared_masks/variants.bed

    # junctions.bed
    awk '$3=="exon" && $1~/^chr([0-9]+|X)$/ {
        s=$4-1; e=$5;
        print $1"\\t"(s-90<0?0:s-90)"\\t"(s+91);
        print $1"\\t"(e-91<0?0:e-91)"\\t"(e+90);
    }' gencode.v49.annotation.gtf \\
        | sort -k1,1 -k2,2n | bedtools merge -i stdin > shared_masks/junctions.bed

    # transcript_strand.bed
    awk '$3=="exon" {
        print $1"\\t"($4-1)"\\t"$5"\\t.\\t.\\t"$7
    }' gencode.v49.annotation.gtf \\
        | sort -k1,1 -k2,2n \\
        | bedtools merge -s -c 6 -o distinct \\
        > shared_masks/transcript_strand.bed
    # Note: positions with both strands (e.g. "+,-") are excluded as ambiguous

    # analysis_regions.bed
    grep -E '^chr([0-9]+|X)\\s' giab.bed | sort -k1,1 -k2,2n > shared_masks/giab_filtered.bed
    bedtools subtract -a shared_masks/giab_filtered.bed -b shared_masks/variants.bed \\
        | sort -k1,1 -k2,2n > shared_masks/giab_minus_variants.bed
    bedtools subtract -a shared_masks/giab_minus_variants.bed -b shared_masks/junctions.bed \\
        | sort -k1,1 -k2,2n > shared_masks/analysis_regions.bed
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

# =============================================================================
# Constants
# =============================================================================

AUTOSOMES_CHRX = [f"chr{i}" for i in range(1, 23)]
CHUNK_SIZE     = 500_000

# Complement map for minus-strand orientation correction
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


# =============================================================================
# Utilities
# =============================================================================

def log(msg):  print(f"[INFO]  {msg}", flush=True)
def warn(msg): print(f"[WARN]  {msg}", flush=True)
def err(msg):
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


def require_file(path, label):
    if not path or not os.path.exists(path):
        err(f"{label} not found: {path}")


def require_tool(name):
    if shutil.which(name) is None:
        err(f"Required tool not found in PATH: {name}\n"
            f"  Install with: conda install -c bioconda {name}")


def run(cmd, desc="", check=True, capture=False):
    if desc:
        log(f"Running: {desc}")
    log(f"  CMD: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    result = subprocess.run(
        cmd, shell=isinstance(cmd, str),
        capture_output=capture, text=True,
        check=False
    )
    if check and result.returncode != 0:
        stderr = result.stderr if capture else "(run with capture=True to see stderr)"
        err(f"Command failed (exit {result.returncode}):\n"
            f"  {cmd}\n  stderr: {stderr}")
    return result


def wilson_ci(successes, total, z=1.96):
    successes = float(successes)
    total     = float(total)
    if total == 0.0:
        return 0.0, 0.0
    p      = successes / total
    denom  = 1 + z**2 / total
    ctr    = (p + z**2 / (2 * total)) / denom
    margin = (z * np.sqrt(p * (1 - p) / total +
                           z**2 / (4 * total**2))) / denom
    return max(0.0, ctr - margin), ctr + margin


# =============================================================================
# CHROM MODE — Step 1: Stream pysamstats TSV → filtered BED
# =============================================================================

def pysamstats_to_filtered_bed(pysamstats_path, out_bed, chroms, min_coverage=1):
    """
    Stream a per-chrom pysamstats TSV to a filtered BED.
    Input is pre-sorted; no sort subprocess needed.
    Returns list of extra column names (all columns after chrom/pos).
    """
    log(f"Streaming pysamstats TSV → filtered BED")
    log(f"  Input: {pysamstats_path}")

    chroms_set = set(chroms)
    opener = gzip.open if pysamstats_path.endswith(".gz") else open
    with opener(pysamstats_path, "rt") as f:
        header = f.readline().strip().split("\t")

    required = ["chrom", "pos", "reads_all", "mismatches", "insertions", "deletions"]
    missing  = [c for c in required if c not in header]
    if missing:
        err(f"Missing columns in pysamstats output: {missing}\nAvailable: {header}")

    extra_cols = [c for c in header if c not in ("chrom", "pos")]

    n_written = n_total = 0
    with open(out_bed, "w") as out:
        for chunk in pd.read_csv(pysamstats_path, sep="\t",
                                  chunksize=CHUNK_SIZE, low_memory=False):
            n_total += len(chunk)
            chunk = chunk[
                chunk["chrom"].isin(chroms_set) &
                (chunk["reads_all"] >= min_coverage)
            ]
            if chunk.empty:
                continue
            chunk = chunk.copy()
            chunk["_start"] = chunk["pos"].astype(int)
            chunk["_end"]   = chunk["_start"] + 1
            chunk[["chrom", "_start", "_end"] + extra_cols].to_csv(
                out, sep="\t", header=False, index=False
            )
            n_written += len(chunk)

    log(f"  Total: {n_total:,}  After filter: {n_written:,}")
    if n_written == 0:
        err("No positions passed chromosome + coverage filter.")
    return extra_cols


# =============================================================================
# CHROM MODE — Step 2: Intersect with pre-built analysis regions
# =============================================================================

def intersect_with_regions(pysamstats_bed, analysis_bed, strand_bed,
                            out_bed, extra_cols, min_coverage):
    """
    Step 1: intersect pysamstats BED with analysis regions.
    Step 2: intersect result with transcript_strand.bed to annotate strand.
            Positions with ambiguous strand (+,-) are excluded.
    """
    log("Intersecting with analysis regions...")
    tmp_regions = out_bed + ".regions.tmp"
    run(f"bedtools intersect -a {pysamstats_bed} -b {analysis_bed} -u > {tmp_regions}",
        desc="bedtools intersect analysis regions")

    n_before = int(run(f"wc -l < {tmp_regions}", capture=True).stdout.strip())
    log(f"  Positions in analysis regions: {n_before:,}")
    if n_before == 0:
        err("No positions intersect analysis regions. "
            "Check chrom name consistency (chr1 vs 1).")

    if "reads_all" in extra_cols:
        cov_col_awk = extra_cols.index("reads_all") + 4
        tmp_cov = out_bed + ".cov.tmp"
        run(f"awk '${cov_col_awk} >= {min_coverage}' {tmp_regions} > {tmp_cov}")
        os.remove(tmp_regions)
    else:
        warn("reads_all column not found — skipping coverage refilter")
        tmp_cov = tmp_regions

    n_after_cov = int(run(f"wc -l < {tmp_cov}", capture=True).stdout.strip())
    log(f"  After min coverage (≥{min_coverage}x): {n_after_cov:,}")
    if n_after_cov == 0:
        err(f"No positions remain after min_coverage={min_coverage}.")

    # Annotate strand from transcript_strand.bed
    # -wa -wb: keep all pysamstats columns, append strand BED columns
    # Last column of strand BED (col 6) is the strand (+, -, or +,-)
    if strand_bed and os.path.exists(strand_bed):
        log("  Annotating strand from transcript_strand.bed...")
        tmp_strand = out_bed + ".strand.tmp"
        # -loj: left outer join — unannotated positions get "." for strand
        run(
            f"bedtools intersect -a {tmp_cov} -b {strand_bed} "
            f"-wa -wb -loj "
            f"| awk '{{print $0}}' > {tmp_strand}",
            desc="bedtools intersect strand"
        )
        # transcript_strand.bed has 6 columns; -wb appends them after original cols
        # strand is always the last column ($NF) of the joined output
        # original columns are 1 to NF-6
        # Keep +, -, AND . (unmatched/intergenic) — only drop ambiguous +,-
        # Deduplicate on chrom+start+end: a position overlapping two exon records
        # (e.g. isoforms) produces multiple rows — keep first (highest-priority strand)
        tmp_annotated = out_bed + ".annotated.tmp"
        run(
            f"awk -v OFS='\\t' '{{"
            f"strand=$NF; "
            f"if (strand==\"+\" || strand==\"-\" || strand==\".\") {{"
            f"n=NF-6; "
            f"for(i=1;i<=n;i++) printf $i OFS; "
            f"print strand"
            f"}}}}' {tmp_strand} "
            f"| awk -v OFS='\\t' '!seen[$1,$2,$3]++' "
            f"> {tmp_annotated}",
            desc="Extract strand-annotated positions (deduplicated)"
        )
        os.remove(tmp_cov)
        os.remove(tmp_strand)
        n_strand = int(run(f"wc -l < {tmp_annotated}", capture=True).stdout.strip())
        log(f"  After strand annotation: {n_strand:,} positions "
            f"(+/- flipped, intergenic/intronic kept as '.')")
        if n_strand == 0:
            err("No positions remain after strand annotation. "
                "Check that transcript_strand.bed uses same chrom names.")
        os.rename(tmp_annotated, out_bed)
        return n_strand, True   # True = strand column present
    else:
        warn("transcript_strand.bed not found — skipping strand correction")
        os.rename(tmp_cov, out_bed)
        n_final = int(run(f"wc -l < {out_bed}", capture=True).stdout.strip())
        return n_final, False   # False = no strand column


# =============================================================================
# CHROM MODE — Step 3: Load BED → per-position TSV
# Output is a minimal TSV: chrom, pos, reads_all, mismatches,
#                          insertions, deletions, ref, A, C, G, T
# This is the only output of chrom mode — merge mode reads it back.
# =============================================================================

def load_and_save_positions(filtered_bed, extra_cols, out_tsv, has_strand):
    """
    Load filtered BED, apply strand correction for minus-strand positions,
    compute per-position rates, write per_position_mismatch.tsv.
    Returns (df_with_rates, n_positions).

    Strand correction (RNA orientation):
      For minus-strand positions, flip ref base and base counts to complement.
      e.g. a minus-strand site with ref=A, G-count=5 becomes ref=T, C-count=5
      This puts all mismatch calls in RNA space regardless of genomic strand.
    """
    log("Loading filtered positions...")
    col_names = ["chrom", "start", "end"] + extra_cols
    if has_strand:
        # strand is always appended as the last column by the awk filter
        # regardless of how many extra_cols there are
        col_names = col_names + ["strand"]

    df = pd.read_csv(filtered_bed, sep="\t", header=None,
                     names=col_names, low_memory=False)
    df["pos"] = df["start"].astype(int)

    for c in ["reads_all", "mismatches", "insertions", "deletions"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["reads_all", "mismatches", "insertions", "deletions"])
    df = df[df["reads_all"] > 0]

    # ── Strand correction ──────────────────────────────────────────────────────
    if has_strand and "strand" in df.columns:
        base_cols_present = [c for c in ["A", "C", "G", "T"] if c in df.columns]
        ref_col_present   = "ref" in df.columns

        minus = df["strand"] == "-"
        n_minus = minus.sum()
        log(f"  Strand correction: {n_minus:,} minus-strand positions "
            f"({n_minus/len(df)*100:.1f}%)")

        if n_minus > 0:
            # Flip ref base to complement
            if ref_col_present:
                df.loc[minus, "ref"] = df.loc[minus, "ref"].map(
                    lambda b: COMP.get(b, b)
                )

            # Flip base counts: swap A↔T and C↔G
            if all(c in df.columns for c in ["A", "T", "C", "G"]):
                # Temp copies to avoid overwriting during swap
                a_vals = df.loc[minus, "A"].copy()
                t_vals = df.loc[minus, "T"].copy()
                c_vals = df.loc[minus, "C"].copy()
                g_vals = df.loc[minus, "G"].copy()
                df.loc[minus, "A"] = t_vals
                df.loc[minus, "T"] = a_vals
                df.loc[minus, "C"] = g_vals
                df.loc[minus, "G"] = c_vals

        log(f"  Strand correction applied — mismatches now in RNA space")
    elif not has_strand:
        warn("No strand column — mismatch attribution is in genomic space "
             "(minus-strand genes will have flipped substitution signatures)")

    # ── Compute per-position rates ─────────────────────────────────────────────
    df["sub_rate"]   = df["mismatches"] / df["reads_all"]
    df["ins_rate"]   = df["insertions"] / df["reads_all"]
    df["del_rate"]   = df["deletions"]  / df["reads_all"]
    df["total_rate"] = (df["mismatches"] + df["insertions"] +
                        df["deletions"]) / df["reads_all"]

    save_cols = ["chrom", "pos", "reads_all", "mismatches", "insertions",
                 "deletions", "sub_rate", "ins_rate", "del_rate", "total_rate"]
    # Include strand, ref, base cols for spectrum and downstream use
    for c in ["strand", "ref", "A", "C", "G", "T"]:
        if c in df.columns:
            save_cols.append(c)

    available = [c for c in save_cols if c in df.columns]
    df[available].to_csv(out_tsv, sep="\t", index=False)
    log(f"  Per-position TSV written: {len(df):,} positions → {out_tsv}")
    return df, len(df)


# =============================================================================
# MERGE MODE — aggregation functions (operate on concatenated per-position TSV)
# =============================================================================

def compute_error_rates(df, label=""):
    for c in ["reads_all", "mismatches", "insertions", "deletions"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["reads_all", "mismatches", "insertions", "deletions"])
    df = df[df["reads_all"] > 0].copy()

    total_bases = df["reads_all"].sum()
    total_sub   = df["mismatches"].sum()
    total_ins   = df["insertions"].sum()
    total_del   = df["deletions"].sum()
    total_err   = total_sub + total_ins + total_del

    sub_rate = total_sub  / total_bases
    ins_rate = total_ins  / total_bases
    del_rate = total_del  / total_bases
    tot_rate = total_err  / total_bases

    sub_lo,  sub_hi  = wilson_ci(total_sub, total_bases)
    ins_lo,  ins_hi  = wilson_ci(total_ins, total_bases)
    del_lo,  del_hi  = wilson_ci(total_del, total_bases)
    tot_lo,  tot_hi  = wilson_ci(total_err, total_bases)

    # Recompute per-position rates if not already present
    if "sub_rate" not in df.columns:
        df["sub_rate"]   = df["mismatches"] / df["reads_all"]
        df["ins_rate"]   = df["insertions"] / df["reads_all"]
        df["del_rate"]   = df["deletions"]  / df["reads_all"]
        df["total_rate"] = (df["mismatches"] + df["insertions"] + df["deletions"]) / df["reads_all"]

    agg = {
        "label":                    label,
        "n_positions":              len(df),
        "total_bases_analysed":     int(total_bases),
        "substitution_rate":        sub_rate,
        "substitution_rate_lo95":   sub_lo,
        "substitution_rate_hi95":   sub_hi,
        "insertion_rate":           ins_rate,
        "insertion_rate_lo95":      ins_lo,
        "insertion_rate_hi95":      ins_hi,
        "deletion_rate":            del_rate,
        "deletion_rate_lo95":       del_lo,
        "deletion_rate_hi95":       del_hi,
        "total_mismatch_rate":      tot_rate,
        "total_mismatch_rate_lo95": tot_lo,
        "total_mismatch_rate_hi95": tot_hi,
    }
    return agg, df


def compute_substitution_spectrum(df):
    if "ref" not in df.columns or not all(c in df.columns for c in ["A","C","G","T"]):
        warn("ref/base columns missing — skipping substitution spectrum")
        return None
    rows = []
    for ref_base in ["A", "C", "G", "T"]:
        ref_rows  = df[df["ref"] == ref_base]
        if ref_rows.empty:
            continue
        ref_total = ref_rows["reads_all"].sum()
        for alt_base in ["A", "C", "G", "T"]:
            if alt_base == ref_base:
                continue
            rate = ref_rows[alt_base].sum() / ref_total if ref_total > 0 else 0
            rows.append({
                "substitution": f"{ref_base}>{alt_base}",
                "ref_base":     ref_base,
                "alt_base":     alt_base,
                "n_ref_bases":  int(ref_total),
                "n_alt_calls":  int(ref_rows[alt_base].sum()),
                "rate":         rate,
            })
    return pd.DataFrame(rows).sort_values("rate", ascending=False)


def load_modkit(path, chroms, min_fraction, min_coverage):
    if not path or not os.path.exists(path):
        warn("modkit BED not found or not provided — skipping")
        return None
    log(f"Loading modkit: {path}")
    try:
        df = pd.read_csv(
            path, sep="\t", header=None, comment="#",
            names=["chrom","start","end","mod_code","score","strand",
                   "thick_start","thick_end","color",
                   "n_valid_cov","frac_modified","n_mod","n_canon",
                   "n_other_mod","n_delete","n_fail","n_diff","n_nocall"]
        )
    except Exception as e:
        warn(f"Could not parse modkit BED: {e}")
        return None
    df = df[df["chrom"].isin(set(chroms))].copy()
    before = len(df)
    df = df[(df["n_valid_cov"] >= min_coverage) &
            (df["frac_modified"] >= min_fraction)]
    log(f"  After filters: {len(df):,} sites (removed {before-len(df):,})")
    return df


def stratify_by_modkit(df, modkit_df, ivt_summary_df=None):
    if modkit_df is None or modkit_df.empty:
        return None, None
    log("Stratifying by modkit modification sites...")

    mod_positions   = {}
    all_mod_pos_set = set()
    for mod_type, grp in modkit_df.groupby("mod_code"):
        s = set(zip(grp["chrom"], grp["start"].astype(int)))
        mod_positions[mod_type] = s
        all_mod_pos_set |= s

    df = df.copy()
    df["is_modified"] = [
        (r_chrom, r_pos) in all_mod_pos_set
        for r_chrom, r_pos in zip(df["chrom"], df["pos"])
    ]

    results = {}
    for strat_label, subset in [
        ("all_sites",        df),
        ("unmodified_sites", df[~df["is_modified"]]),
        ("modified_sites",   df[df["is_modified"]]),
    ]:
        if subset.empty:
            continue
        agg, _ = compute_error_rates(subset.copy(), label=strat_label)
        results[strat_label] = agg

    for mod_type, mod_set in mod_positions.items():
        mask   = [(c, p) in mod_set for c, p in zip(df["chrom"], df["pos"])]
        subset = df[mask]
        if len(subset) < 10:
            continue
        agg, _ = compute_error_rates(subset.copy(), label=f"mod_{mod_type}")
        results[f"mod_{mod_type}"] = agg

    return results, pd.DataFrame(list(results.values()))


def _get_ivt_rate(df, key):
    if df is None:
        return None
    if "metric" in df.columns:
        row = df[df["metric"] == key]
        return float(row["value"].iloc[0]) if not row.empty else None
    return float(df[key].iloc[0]) if key in df.columns else None


def compute_delta(native_agg, ivt_summary_df):
    if ivt_summary_df is None:
        warn("IVT summary not provided — skipping native−IVT delta")
        return None
    pairs = [
        ("substitution_rate",   "substitution_rate"),
        ("insertion_rate",      "insertion_rate"),
        ("deletion_rate",       "deletion_rate"),
        ("total_mismatch_rate", "total_error_rate"),
    ]
    rows = []
    for nat_key, ivt_key in pairs:
        nat_val = native_agg.get(nat_key, np.nan)
        ivt_val = _get_ivt_rate(ivt_summary_df, ivt_key)
        if ivt_val is None:
            continue
        rows.append({
            "metric":      nat_key,
            "ivt_rate":    ivt_val,
            "native_rate": nat_val,
            "delta":       nat_val - ivt_val,
            "fold_change": nat_val / ivt_val if ivt_val > 0 else np.nan,
        })
    delta_df = pd.DataFrame(rows)
    print(f"\n{'='*68}")
    print(f"  NATIVE − IVT DELTA")
    print(f"{'='*68}")
    for _, r in delta_df.iterrows():
        print(f"  {r['metric']:<28} "
              f"IVT={r['ivt_rate']*100:.5f}%  "
              f"Native={r['native_rate']*100:.5f}%  "
              f"Δ={r['delta']*100:+.5f}%  "
              f"({r['fold_change']:.2f}×)")
    print(f"{'='*68}")
    return delta_df


def build_masking_report(masks_dir, giab_bed, analysis_bed,
                          final_n_positions, min_coverage):
    log("Building masking report...")

    def count_bases(bed):
        if not bed or not os.path.exists(bed):
            return "N/A"
        r = run(f"awk '{{sum+=$3-$2}} END{{print sum+0}}' {bed}", capture=True)
        return int(r.stdout.strip() or "0")

    giab_minus_variants = os.path.join(masks_dir, "giab_minus_variants.bed")
    minus_junctions     = os.path.join(masks_dir, "minus_junctions.bed")

    return pd.DataFrame([
        {"step": "1. GIAB high-conf regions",
         "region_bases": count_bases(giab_bed)},
        {"step": "2. Minus variant sites (HET+HOM)",
         "region_bases": count_bases(giab_minus_variants)},
        {"step": "3. Minus junction exclusion zones",
         "region_bases": count_bases(minus_junctions)},
        {"step": "4. Minus editing sites (SKIPPED — by design)",
         "region_bases": count_bases(analysis_bed)},
        {"step": f"5. Positions with ≥{min_coverage}x coverage",
         "region_bases": f"{final_n_positions:,} positions"},
    ])


def save_outputs(outdir, label, agg, mask_report, df_with_rates,
                 spectrum_df, modkit_strat_df, delta_df,
                 editing_applied, min_coverage, junction_bp):
    Path(outdir).mkdir(parents=True, exist_ok=True)

    agg["junction_exclusion_bp"] = junction_bp
    agg["min_coverage"]          = min_coverage
    agg["editing_mask_applied"]  = editing_applied
    agg["note"] = ("Reference mismatch rate — not true platform error rate. "
                   "Compare to IVT for platform baseline.")

    pd.DataFrame([agg]).to_csv(
        f"{outdir}/{label}_mismatch_rate_summary.tsv", sep="\t", index=False)
    mask_report.to_csv(
        f"{outdir}/{label}_masking_report.tsv", sep="\t", index=False)

    cols = ["chrom", "pos", "reads_all", "mismatches", "insertions",
            "deletions", "sub_rate", "ins_rate", "del_rate", "total_rate"]
    available = [c for c in cols if c in df_with_rates.columns]
    df_with_rates[available].to_csv(
        f"{outdir}/{label}_per_position_mismatch.tsv", sep="\t", index=False)

    if spectrum_df is not None:
        spectrum_df.to_csv(
            f"{outdir}/{label}_substitution_spectrum.tsv", sep="\t", index=False)
    if modkit_strat_df is not None:
        modkit_strat_df.to_csv(
            f"{outdir}/{label}_modkit_stratified_rates.tsv", sep="\t", index=False)
    if delta_df is not None:
        delta_df.to_csv(
            f"{outdir}/{label}_native_ivt_delta.tsv", sep="\t", index=False)

    log(f"Outputs written to: {outdir}/")


def print_summary(agg, mask_report, junction_bp, editing_applied, label):
    print(f"\n{'='*68}")
    print(f"  NATIVE DRS MISMATCH RATE — {label}")
    print(f"{'='*68}")
    print(f"  NOTE: Reference mismatch rate — not true platform error rate.")
    print(f"  Editing mask: {'Applied' if editing_applied else 'NOT applied — A>G inflated'}")
    print(f"  Junction exclusion: {junction_bp}bp")
    print()
    for display, key in [
        ("Substitution",   "substitution_rate"),
        ("Insertion",      "insertion_rate"),
        ("Deletion",       "deletion_rate"),
        ("Total mismatch", "total_mismatch_rate"),
    ]:
        val = agg[key]
        lo  = agg.get(f"{key}_lo95", np.nan)
        hi  = agg.get(f"{key}_hi95", np.nan)
        ci  = f"[{lo*100:.5f}%, {hi*100:.5f}%]" if not np.isnan(lo) else ""
        print(f"  {display:<18} {val*100:.5f}%  {ci}")
    print(f"\n  Masking summary:")
    print(f"  {mask_report.to_string(index=False)}")
    print(f"{'='*68}")


# =============================================================================
# Argument parsing
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Native DRS mismatch rate — chrom or merge mode"
    )
    p.add_argument("--mode", choices=["chrom", "merge"], default="chrom",
                   help="chrom: process one per-chrom pysamstats TSV. "
                        "merge: aggregate per-chrom per_position TSVs into final summary.")

    # ── chrom mode args ───────────────────────────────────────────────────────
    chrom_grp = p.add_argument_group("chrom mode")
    chrom_grp.add_argument("--pysamstats",
                           help="Per-chrom pre-sorted pysamstats TSV (or .gz)")
    chrom_grp.add_argument("--masks-dir",
                           help="Directory with analysis_regions.bed, variants.bed, junctions.bed")
    chrom_grp.add_argument("--chrom",
                           help="Chromosome being processed (e.g. chr1). "
                                "Used to filter TSV and name output.")
    chrom_grp.add_argument("--strand-bed", default=None,
                           help="transcript_strand.bed (in --masks-dir). "
                                "Flips ref+base counts to RNA orientation for "
                                "minus-strand positions. If absent, mismatches "
                                "are reported in genomic space (not recommended).")

    # ── merge mode args ───────────────────────────────────────────────────────
    merge_grp = p.add_argument_group("merge mode")
    merge_grp.add_argument("--position-tsvs", nargs="+",
                           help="Per-chrom per_position_mismatch.tsv files to merge")
    merge_grp.add_argument("--masks-dir-merge",
                           help="Masks dir for masking report (same shared_masks/)")
    merge_grp.add_argument("--giab-bed",
                           help="GIAB BED for masking report")
    merge_grp.add_argument("--modkit",       default=None)
    merge_grp.add_argument("--ivt-summary",  default=None)

    # ── shared args ───────────────────────────────────────────────────────────
    p.add_argument("--outdir",           default="results/native")
    p.add_argument("--sample-label",     default="GM12878_native")
    p.add_argument("--min-coverage",     type=int,   default=20)
    p.add_argument("--junction-bp",      type=int,   default=90)
    p.add_argument("--min-mod-fraction", type=float, default=0.5)
    p.add_argument("--min-mod-coverage", type=int,   default=20)
    p.add_argument("--chroms",           default=None)
    p.add_argument("--keep-tmp",         action="store_true")
    return p.parse_args()


# =============================================================================
# Main
# =============================================================================

def main():
    args   = parse_args()
    chroms = args.chroms.split(",") if args.chroms else AUTOSOMES_CHRX
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # CHROM MODE
    # =========================================================================
    if args.mode == "chrom":
        require_tool("bedtools")
        if not args.pysamstats:
            err("--pysamstats required in chrom mode")
        if not args.masks_dir:
            err("--masks-dir required in chrom mode")

        require_file(args.pysamstats, "pysamstats TSV")
        masks_dir    = args.masks_dir.rstrip("/")
        analysis_bed = os.path.join(masks_dir, "analysis_regions.bed")
        require_file(analysis_bed, "analysis_regions.bed")

        # Strand BED — look in masks_dir if not explicitly passed
        strand_bed = args.strand_bed or os.path.join(masks_dir, "transcript_strand.bed")
        if not os.path.exists(strand_bed):
            warn("transcript_strand.bed not found — proceeding without strand correction. "
                 "Substitution spectrum will be in genomic not RNA space.")
            strand_bed = None

        # Restrict to single chrom if provided
        chrom_filter = [args.chrom] if args.chrom else chroms

        tmpdir = tempfile.mkdtemp(prefix="drs_chrom_", dir=args.outdir)
        log(f"Temporary directory: {tmpdir}")
        log(f"Processing chrom(s): {chrom_filter}")

        try:
            pysamstats_bed = f"{tmpdir}/pysamstats.bed"
            extra_cols = pysamstats_to_filtered_bed(
                args.pysamstats, pysamstats_bed,
                chrom_filter, min_coverage=max(1, args.min_coverage)
            )

            filtered_bed = f"{tmpdir}/filtered_positions.bed"
            n_positions, has_strand = intersect_with_regions(
                pysamstats_bed, analysis_bed, strand_bed,
                filtered_bed, extra_cols, args.min_coverage
            )

            chrom_tag = args.chrom if args.chrom else "all"
            out_tsv   = os.path.join(
                args.outdir,
                f"{args.sample_label}_{chrom_tag}_per_position_mismatch.tsv"
            )
            load_and_save_positions(filtered_bed, extra_cols, out_tsv, has_strand)

        finally:
            if not args.keep_tmp:
                shutil.rmtree(tmpdir, ignore_errors=True)

        log(f"Chrom mode done → {out_tsv}")

    # =========================================================================
    # MERGE MODE — streaming accumulator, one TSV at a time
    # Never holds more than one chrom TSV in RAM simultaneously.
    # Accumulates raw count sums for rates/spectrum/modkit/delta.
    # Writes per-position TSV in append mode chrom by chrom.
    # =========================================================================
    elif args.mode == "merge":
        if not args.position_tsvs:
            err("--position-tsvs required in merge mode")

        ivt_summary = (pd.read_csv(args.ivt_summary, sep="\t")
                       if args.ivt_summary and os.path.exists(args.ivt_summary)
                       else None)
        modkit_df = load_modkit(args.modkit, chroms,
                                args.min_mod_fraction, args.min_mod_coverage)

        # Build modkit position set once (small, keep in RAM)
        all_mod_pos_set  = set()
        mod_positions    = {}
        if modkit_df is not None and not modkit_df.empty:
            for mod_type, grp in modkit_df.groupby("mod_code"):
                s = set(zip(grp["chrom"], grp["start"].astype(int)))
                mod_positions[mod_type] = s
                all_mod_pos_set |= s
            log(f"modkit positions loaded into RAM: {len(all_mod_pos_set):,}")

        # Running accumulators
        total_reads_all  = 0
        total_mismatches = 0
        total_insertions = 0
        total_deletions  = 0
        n_positions      = 0

        # Spectrum accumulators: ref_base → {alt_base → count}
        ref_base_list = ["A", "C", "G", "T"]
        spec_counts   = {r: {a: 0 for a in ref_base_list if a != r}
                         for r in ref_base_list}
        ref_totals    = {r: 0 for r in ref_base_list}

        # modkit accumulators: strat_label → {reads_all, mismatches, ins, del}
        strat_keys = ["all_sites", "unmodified_sites", "modified_sites"] + \
                     [f"mod_{k}" for k in mod_positions.keys()]
        strat_acc  = {k: {"reads_all":0,"mismatches":0,
                          "insertions":0,"deletions":0}
                      for k in strat_keys}

        # Output per-position TSV (written in append mode)
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        out_pos_tsv = os.path.join(
            args.outdir,
            f"{args.sample_label}_per_position_mismatch.tsv"
        )
        pos_header_written = False

        log(f"Merge mode: streaming {len(args.position_tsvs)} per-chrom TSVs")

        for tsv_path in args.position_tsvs:
            if not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0:
                warn(f"Skipping missing/empty TSV: {tsv_path}")
                continue

            log(f"  Processing: {os.path.basename(tsv_path)}")
            df = pd.read_csv(tsv_path, sep="\t", low_memory=False)

            for c in ["reads_all", "mismatches", "insertions", "deletions"]:
                df[c] = pd.to_numeric(df[c], errors="coerce")
            df = df.dropna(subset=["reads_all","mismatches",
                                    "insertions","deletions"])
            df = df[df["reads_all"] > 0]
            if df.empty:
                continue

            # Recompute per-position rates
            df = df.copy()
            df["sub_rate"]   = df["mismatches"] / df["reads_all"]
            df["ins_rate"]   = df["insertions"] / df["reads_all"]
            df["del_rate"]   = df["deletions"]  / df["reads_all"]
            df["total_rate"] = (df["mismatches"] + df["insertions"] +
                                df["deletions"]) / df["reads_all"]

            # Accumulate raw counts
            total_reads_all  += int(df["reads_all"].sum())
            total_mismatches += int(df["mismatches"].sum())
            total_insertions += int(df["insertions"].sum())
            total_deletions  += int(df["deletions"].sum())
            n_positions      += len(df)

            # Substitution spectrum
            if all(c in df.columns for c in ["ref", "A", "C", "G", "T"]):
                for ref_base in ref_base_list:
                    sub = df[df["ref"] == ref_base]
                    if sub.empty:
                        continue
                    ref_totals[ref_base] += int(sub["reads_all"].sum())
                    for alt_base in ref_base_list:
                        if alt_base == ref_base:
                            continue
                        spec_counts[ref_base][alt_base] += int(sub[alt_base].sum())

            # modkit stratification accumulators
            if all_mod_pos_set:
                pos_set = set(zip(df["chrom"], df["pos"].astype(int)))
                is_mod  = np.array([(c, p) in all_mod_pos_set
                                    for c, p in zip(df["chrom"],
                                                    df["pos"].astype(int))])
                for strat_label, mask in [
                    ("all_sites",        np.ones(len(df), dtype=bool)),
                    ("unmodified_sites", ~is_mod),
                    ("modified_sites",   is_mod),
                ]:
                    sub = df[mask]
                    if sub.empty:
                        continue
                    strat_acc[strat_label]["reads_all"]  += int(sub["reads_all"].sum())
                    strat_acc[strat_label]["mismatches"] += int(sub["mismatches"].sum())
                    strat_acc[strat_label]["insertions"] += int(sub["insertions"].sum())
                    strat_acc[strat_label]["deletions"]  += int(sub["deletions"].sum())

                for mod_type, mod_set in mod_positions.items():
                    mod_mask = np.array([(c, p) in mod_set
                                         for c, p in zip(df["chrom"],
                                                         df["pos"].astype(int))])
                    sub = df[mod_mask]
                    if sub.empty:
                        continue
                    key = f"mod_{mod_type}"
                    strat_acc[key]["reads_all"]  += int(sub["reads_all"].sum())
                    strat_acc[key]["mismatches"] += int(sub["mismatches"].sum())
                    strat_acc[key]["insertions"] += int(sub["insertions"].sum())
                    strat_acc[key]["deletions"]  += int(sub["deletions"].sum())

            # Write per-position rows (append mode)
            pos_cols  = ["chrom","pos","reads_all","mismatches","insertions",
                         "deletions","sub_rate","ins_rate","del_rate","total_rate"]
            available = [c for c in pos_cols if c in df.columns]
            df[available].to_csv(
                out_pos_tsv, sep="\t", index=False,
                header=not pos_header_written,
                mode="a" if pos_header_written else "w"
            )
            pos_header_written = True

            del df   # free RAM immediately

        if n_positions == 0:
            err("No positions found across all per-chrom TSVs.")

        log(f"  Total positions: {n_positions:,}")

        # ── Compute aggregate rates from sums ──────────────────────────────────
        def _rates_from_counts(ra, mm, ins, dl, label):
            if ra == 0:
                return {}
            sub_lo,  sub_hi  = wilson_ci(mm,       ra)
            ins_lo,  ins_hi  = wilson_ci(ins,       ra)
            del_lo,  del_hi  = wilson_ci(dl,        ra)
            tot_lo,  tot_hi  = wilson_ci(mm+ins+dl, ra)
            return {
                "label":                    label,
                "n_positions":              n_positions,
                "total_bases_analysed":     ra,
                "substitution_rate":        mm / ra,
                "substitution_rate_lo95":   sub_lo,
                "substitution_rate_hi95":   sub_hi,
                "insertion_rate":           ins / ra,
                "insertion_rate_lo95":      ins_lo,
                "insertion_rate_hi95":      ins_hi,
                "deletion_rate":            dl / ra,
                "deletion_rate_lo95":       del_lo,
                "deletion_rate_hi95":       del_hi,
                "total_mismatch_rate":      (mm+ins+dl) / ra,
                "total_mismatch_rate_lo95": tot_lo,
                "total_mismatch_rate_hi95": tot_hi,
            }

        agg = _rates_from_counts(
            total_reads_all, total_mismatches,
            total_insertions, total_deletions,
            args.sample_label
        )

        # ── Substitution spectrum from accumulated counts ───────────────────────
        spectrum_rows = []
        for ref_base in ref_base_list:
            denom = ref_totals[ref_base]
            for alt_base in ref_base_list:
                if alt_base == ref_base:
                    continue
                cnt = spec_counts[ref_base][alt_base]
                spectrum_rows.append({
                    "substitution": f"{ref_base}>{alt_base}",
                    "ref_base":     ref_base,
                    "alt_base":     alt_base,
                    "n_ref_bases":  denom,
                    "n_alt_calls":  cnt,
                    "rate":         cnt / denom if denom > 0 else 0,
                })
        spectrum_df = pd.DataFrame(spectrum_rows).sort_values(
            "rate", ascending=False) if spectrum_rows else None

        # ── modkit stratified rates from accumulated counts ────────────────────
        modkit_strat_rows = []
        for strat_label, acc in strat_acc.items():
            ra = acc["reads_all"]
            if ra == 0:
                continue
            mm  = acc["mismatches"]
            ins = acc["insertions"]
            dl  = acc["deletions"]
            modkit_strat_rows.append({
                "label":               strat_label,
                "n_positions":         ra,
                "substitution_rate":   mm / ra,
                "insertion_rate":      ins / ra,
                "deletion_rate":       dl / ra,
                "total_mismatch_rate": (mm+ins+dl) / ra,
            })
        modkit_strat_df = pd.DataFrame(modkit_strat_rows) \
                          if modkit_strat_rows else None

        # ── Native − IVT delta ────────────────────────────────────────────────
        delta_df = compute_delta(agg, ivt_summary)

        # ── Masking report ────────────────────────────────────────────────────
        masks_dir_merge = (args.masks_dir_merge or "").rstrip("/")
        analysis_bed    = os.path.join(masks_dir_merge, "analysis_regions.bed") \
                          if masks_dir_merge else ""
        mask_report = build_masking_report(
            masks_dir_merge, args.giab_bed, analysis_bed,
            n_positions, args.min_coverage
        )

        # ── Print and save ─────────────────────────────────────────────────────
        editing_applied = False
        print_summary(agg, mask_report, args.junction_bp,
                      editing_applied, args.sample_label)

        # save_outputs would overwrite per_position_mismatch.tsv — skip that
        # and write remaining outputs directly
        agg["junction_exclusion_bp"] = args.junction_bp
        agg["min_coverage"]          = args.min_coverage
        agg["editing_mask_applied"]  = editing_applied
        agg["note"] = ("Reference mismatch rate — not true platform error rate. "
                       "Compare to IVT for platform baseline.")

        pd.DataFrame([agg]).to_csv(
            f"{args.outdir}/{args.sample_label}_mismatch_rate_summary.tsv",
            sep="\t", index=False)
        mask_report.to_csv(
            f"{args.outdir}/{args.sample_label}_masking_report.tsv",
            sep="\t", index=False)
        if spectrum_df is not None:
            spectrum_df.to_csv(
                f"{args.outdir}/{args.sample_label}_substitution_spectrum.tsv",
                sep="\t", index=False)
        if modkit_strat_df is not None:
            modkit_strat_df.to_csv(
                f"{args.outdir}/{args.sample_label}_modkit_stratified_rates.tsv",
                sep="\t", index=False)
        if delta_df is not None:
            delta_df.to_csv(
                f"{args.outdir}/{args.sample_label}_native_ivt_delta.tsv",
                sep="\t", index=False)

        log(f"Outputs written to: {args.outdir}/")
        log("Merge mode done.")


if __name__ == "__main__":
    main()
