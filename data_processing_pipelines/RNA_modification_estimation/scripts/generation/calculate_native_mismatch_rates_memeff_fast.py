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

AUTOSOMES_CHRX = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CHUNK_SIZE     = 500_000


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

def intersect_with_regions(pysamstats_bed, analysis_bed, out_bed,
                            extra_cols, min_coverage):
    log("Intersecting with analysis regions...")
    tmp = out_bed + ".tmp"
    run(f"bedtools intersect -a {pysamstats_bed} -b {analysis_bed} -u > {tmp}",
        desc="bedtools intersect")

    n_before = int(run(f"wc -l < {tmp}", capture=True).stdout.strip())
    log(f"  Positions in analysis regions: {n_before:,}")
    if n_before == 0:
        err("No positions intersect analysis regions. "
            "Check chrom name consistency (chr1 vs 1).")

    if "reads_all" in extra_cols:
        cov_col_awk = extra_cols.index("reads_all") + 4
        run(f"awk '${cov_col_awk} >= {min_coverage}' {tmp} > {out_bed}")
    else:
        warn("reads_all column not found — skipping coverage refilter")
        shutil.copy(tmp, out_bed)
    os.remove(tmp)

    n_after = int(run(f"wc -l < {out_bed}", capture=True).stdout.strip())
    log(f"  After min coverage (≥{min_coverage}x): {n_after:,}")
    if n_after == 0:
        err(f"No positions remain after min_coverage={min_coverage}.")
    return n_after


# =============================================================================
# CHROM MODE — Step 3: Load BED → per-position TSV
# Output is a minimal TSV: chrom, pos, reads_all, mismatches,
#                          insertions, deletions, ref, A, C, G, T
# This is the only output of chrom mode — merge mode reads it back.
# =============================================================================

def load_and_save_positions(filtered_bed, extra_cols, out_tsv):
    """
    Load filtered BED, compute per-position rates, write per_position_mismatch.tsv.
    Returns (df_with_rates, n_positions).
    """
    log("Loading filtered positions...")
    col_names = ["chrom", "start", "end"] + extra_cols
    df = pd.read_csv(filtered_bed, sep="\t", header=None,
                     names=col_names, low_memory=False)
    df["pos"] = df["start"].astype(int)

    for c in ["reads_all", "mismatches", "insertions", "deletions"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["reads_all", "mismatches", "insertions", "deletions"])
    df = df[df["reads_all"] > 0]

    df["sub_rate"]   = df["mismatches"] / df["reads_all"]
    df["ins_rate"]   = df["insertions"] / df["reads_all"]
    df["del_rate"]   = df["deletions"]  / df["reads_all"]
    df["total_rate"] = (df["mismatches"] + df["insertions"] + df["deletions"]) / df["reads_all"]

    save_cols = ["chrom", "pos", "reads_all", "mismatches", "insertions",
                 "deletions", "sub_rate", "ins_rate", "del_rate", "total_rate"]
    # include ref + base cols for spectrum if present
    for c in ["ref", "A", "C", "G", "T"]:
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
            n_positions  = intersect_with_regions(
                pysamstats_bed, analysis_bed, filtered_bed,
                extra_cols, args.min_coverage
            )

            chrom_tag = args.chrom if args.chrom else "all"
            out_tsv   = os.path.join(
                args.outdir,
                f"{args.sample_label}_{chrom_tag}_per_position_mismatch.tsv"
            )
            load_and_save_positions(filtered_bed, extra_cols, out_tsv)

        finally:
            if not args.keep_tmp:
                shutil.rmtree(tmpdir, ignore_errors=True)

        log(f"Chrom mode done → {out_tsv}")

    # =========================================================================
    # MERGE MODE
    # =========================================================================
    elif args.mode == "merge":
        if not args.position_tsvs:
            err("--position-tsvs required in merge mode")

        log(f"Merge mode: concatenating {len(args.position_tsvs)} per-chrom TSVs")
        dfs = []
        for f in args.position_tsvs:
            if not os.path.exists(f):
                warn(f"Missing per-chrom TSV (skipping): {f}")
                continue
            dfs.append(pd.read_csv(f, sep="\t", low_memory=False))

        if not dfs:
            err("No per-chrom TSVs could be loaded.")

        df_all = pd.concat(dfs, ignore_index=True)
        log(f"  Total positions after merge: {len(df_all):,}")

        # ── Aggregate rates ────────────────────────────────────────────────────
        agg, df_with_rates = compute_error_rates(df_all, label=args.sample_label)

        # ── Substitution spectrum ──────────────────────────────────────────────
        spectrum_df = compute_substitution_spectrum(df_with_rates)
        if spectrum_df is not None:
            log(f"Top substitutions:")
            log(spectrum_df.head(6).to_string(index=False))

        # ── modkit stratification ──────────────────────────────────────────────
        ivt_summary = (pd.read_csv(args.ivt_summary, sep="\t")
                       if args.ivt_summary and os.path.exists(args.ivt_summary)
                       else None)
        modkit_df = load_modkit(args.modkit, chroms,
                                args.min_mod_fraction, args.min_mod_coverage)
        _, modkit_strat_df = stratify_by_modkit(df_with_rates, modkit_df, ivt_summary)

        # ── Native − IVT delta ────────────────────────────────────────────────
        delta_df = compute_delta(agg, ivt_summary)

        # ── Masking report ────────────────────────────────────────────────────
        masks_dir_merge = (args.masks_dir_merge or "").rstrip("/")
        analysis_bed    = os.path.join(masks_dir_merge, "analysis_regions.bed") \
                          if masks_dir_merge else ""
        mask_report = build_masking_report(
            masks_dir_merge, args.giab_bed, analysis_bed,
            len(df_with_rates), args.min_coverage
        )

        # ── Print and save ─────────────────────────────────────────────────────
        editing_applied = False   # REDIportal masking omitted by design
        print_summary(agg, mask_report, args.junction_bp,
                      editing_applied, args.sample_label)
        save_outputs(
            args.outdir, args.sample_label,
            agg, mask_report, df_with_rates,
            spectrum_df, modkit_strat_df, delta_df,
            editing_applied, args.min_coverage, args.junction_bp
        )

        log("Merge mode done.")


if __name__ == "__main__":
    main()
