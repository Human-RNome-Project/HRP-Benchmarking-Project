#!/usr/bin/env python3
"""
calculate_coverage_error_rate.py

Calculates per-position error rates binned by coverage depth.

Answers the question: does sequencing error rate depend on read depth?

For a platform where errors are random and independent across reads:
  - Error rate should be CONSTANT across coverage bins
  - High coverage positions should have the same error rate as low coverage

Deviations from this expectation reveal:
  - Systematic errors (modification signals, alignment artefacts) that
    are diluted at higher coverage but present at all depths
  - Coverage-dependent calling biases in the basecaller
  - Regions with both high coverage and high error (e.g. repetitive regions
    that attract spurious alignments)

Two input modes:
  Mode A (recommended): takes pysamstats per-position TSV directly
  Mode B:               takes a filtered_positions BED from the pipeline tmpdir

Usage:
    python calculate_coverage_error_rate.py \
        --pysamstats    native_variation_per_position.tsv \
        --outdir        results/coverage_error \
        --sample-label  GM12878_native \
        --giab-bed      HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
        --chroms        chr19
"""

import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Constants
# =============================================================================

AUTOSOMES_CHRX = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CHUNK_SIZE     = 500_000

# Coverage bins used throughout — chosen to give roughly equal log-spacing
# and enough positions in each bin for reliable estimates
COV_BINS   = [1, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 1000, np.inf]
COV_LABELS = ["1–5", "5–10", "10–20", "20–30", "30–50", "50–75",
              "75–100", "100–150", "150–200", "200–300", "300–500",
              "500–1k", "1k+"]


# =============================================================================
# Utilities
# =============================================================================

def log(msg):  print(f"[INFO]  {msg}", flush=True)
def warn(msg): print(f"[WARN]  {msg}", flush=True)
def err(msg):
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


# =============================================================================
# Step 1 — Load pysamstats TSV in chunks, applying optional GIAB filter
# =============================================================================

def load_pysamstats_chunked(pysamstats_path, chroms, giab_set=None):
    """
    Stream pysamstats TSV in chunks and compute per-position error rates.
    Optionally restrict to GIAB high-confidence regions.
    Returns a DataFrame with essential columns only to keep RAM low.
    """
    log(f"Loading pysamstats: {pysamstats_path}")
    chroms_set = set(chroms)

    opener = gzip.open if pysamstats_path.endswith(".gz") else open
    with opener(pysamstats_path, "rt") as f:
        header = f.readline().strip().split("\t")

    required = ["chrom", "pos", "reads_all", "mismatches",
                "insertions", "deletions"]
    missing  = [c for c in required if c not in header]
    if missing:
        err(f"Missing columns: {missing}\nAvailable: {header}")

    chunks_out = []
    n_total    = 0
    n_kept     = 0

    for chunk in pd.read_csv(
        pysamstats_path, sep="\t", chunksize=CHUNK_SIZE, low_memory=False
    ):
        n_total += len(chunk)

        # Chromosome filter
        chunk = chunk[chunk["chrom"].isin(chroms_set)].copy()

        # Zero-coverage filter
        chunk = chunk[chunk["reads_all"] > 0]

        # GIAB filter
        if giab_set is not None:
            chunk = chunk[
                chunk.apply(
                    lambda r: (r["chrom"], r["pos"]) in giab_set, axis=1
                )
            ]

        if chunk.empty:
            continue

        # Compute per-position rates
        chunk["reads_all"]  = pd.to_numeric(chunk["reads_all"],  errors="coerce")
        chunk["mismatches"] = pd.to_numeric(chunk["mismatches"], errors="coerce")
        chunk["insertions"] = pd.to_numeric(chunk["insertions"], errors="coerce")
        chunk["deletions"]  = pd.to_numeric(chunk["deletions"],  errors="coerce")
        chunk = chunk.dropna(subset=required)
        chunk = chunk[chunk["reads_all"] > 0]

        chunk["sub_rate"]   = chunk["mismatches"] / chunk["reads_all"]
        chunk["ins_rate"]   = chunk["insertions"] / chunk["reads_all"]
        chunk["del_rate"]   = chunk["deletions"]  / chunk["reads_all"]
        chunk["total_rate"] = (chunk["mismatches"] +
                               chunk["insertions"] +
                               chunk["deletions"]) / chunk["reads_all"]

        # Keep only essential columns
        keep = ["chrom", "pos", "reads_all",
                "mismatches", "insertions", "deletions",
                "sub_rate", "ins_rate", "del_rate", "total_rate"]
        if "ref" in chunk.columns:
            keep.append("ref")

        chunks_out.append(chunk[keep])
        n_kept += len(chunk)

    log(f"  Rows loaded: {n_total:,}  |  After filters: {n_kept:,}")

    if n_kept == 0:
        err("No positions loaded. Check chromosome naming and file path.")

    return pd.concat(chunks_out, ignore_index=True)


# =============================================================================
# Step 2 — Load GIAB BED as position set (optional, for validation runs)
# Note: for genome-wide runs this is memory-intensive. Use --chroms to restrict.
# =============================================================================

def load_giab_bed_set(bed_path, chroms):
    """
    Load GIAB high-confidence regions into a set of (chrom, pos) tuples.
    Only use this when restricting to specific chromosomes — genome-wide
    will OOM for the same reason as the main pipeline.
    """
    if bed_path is None:
        return None
    log(f"Loading GIAB BED: {bed_path}")
    chroms_set = set(chroms)
    pos_set    = set()
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if chrom not in chroms_set:
                continue
            for pos in range(start, end):
                pos_set.add((chrom, pos))
    log(f"  GIAB positions loaded: {len(pos_set):,}")
    return pos_set


# =============================================================================
# Step 3 — Bin positions by coverage and compute aggregate error rates
# =============================================================================

def bin_by_coverage(df, label):
    """
    For each coverage bin compute:
      - n_sites: number of positions in this bin
      - total_bases: sum of reads_all across positions
      - weighted aggregate error rate for each error type
        (sum of counts / sum of bases — not mean of rates)
      - median per-position rate
      - 25th and 75th percentile per-position rate

    Weighted aggregate rate is the primary metric — it is consistent
    with how the main pipeline computes error rates and is not biased
    by the distribution of coverage within the bin.
    """
    log("Binning positions by coverage depth...")

    df["cov_bin"] = pd.cut(
        df["reads_all"], bins=COV_BINS, labels=COV_LABELS, right=False
    )

    rows = []
    for bin_label in COV_LABELS:
        subset = df[df["cov_bin"] == bin_label]
        if subset.empty:
            continue

        total_bases = float(subset["reads_all"].sum())
        total_sub   = float(subset["mismatches"].sum())
        total_ins   = float(subset["insertions"].sum())
        total_del   = float(subset["deletions"].sum())
        total_err   = total_sub + total_ins + total_del

        row = {
            "label":          label,
            "cov_bin":        bin_label,
            "n_sites":        len(subset),
            "total_bases":    int(total_bases),
            "cov_min":        subset["reads_all"].min(),
            "cov_max":        subset["reads_all"].max(),
            "cov_median":     subset["reads_all"].median(),
        }

        for rate_name, count in [
            ("sub",   total_sub),
            ("ins",   total_ins),
            ("del",   total_del),
            ("total", total_err),
        ]:
            col = f"{rate_name}_rate" if rate_name != "total" else "total_rate"
            row[f"{rate_name}_weighted_rate"] = count / total_bases if total_bases > 0 else np.nan
            row[f"{rate_name}_median_rate"]   = float(subset[col].median()) if col in subset.columns else np.nan
            row[f"{rate_name}_p25_rate"]      = float(subset[col].quantile(0.25)) if col in subset.columns else np.nan
            row[f"{rate_name}_p75_rate"]      = float(subset[col].quantile(0.75)) if col in subset.columns else np.nan

        rows.append(row)

    result = pd.DataFrame(rows)
    log(f"  Coverage bins with data: {len(result)}")

    # Print summary table
    print(f"\n{'='*80}")
    print(f"  COVERAGE vs ERROR RATE — {label}")
    print(f"{'='*80}")
    print(f"  {'Bin':<10} {'N sites':>10} {'Total bases':>14} "
          f"{'Sub%':>8} {'Del%':>8} {'Ins%':>8} {'Total%':>8}")
    print(f"  {'-'*70}")
    for _, row in result.iterrows():
        print(f"  {row['cov_bin']:<10} {row['n_sites']:>10,} "
              f"{row['total_bases']:>14,} "
              f"{row['sub_weighted_rate']*100:>7.4f}% "
              f"{row['del_weighted_rate']*100:>7.4f}% "
              f"{row['ins_weighted_rate']*100:>7.4f}% "
              f"{row['total_weighted_rate']*100:>7.4f}%")
    print(f"{'='*80}\n")

    return result


# =============================================================================
# Step 4 — Per-base-quality stratification (if base quality available)
# Nanopore per-base quality scores are noisy but provide a useful stratifier
# =============================================================================

def stratify_by_ref_base(df, label):
    """
    Compute error rates stratified by reference base.
    Reveals base-specific error biases — important for DRS where
    C>T and G>A mismatches are dominant (as seen in Fig 3).
    Only possible if 'ref' column is present in pysamstats output.
    """
    if "ref" not in df.columns:
        warn("No 'ref' column — skipping per-base stratification")
        return None

    log("Computing per-reference-base error rates...")

    rows = []
    for ref_base in ["A", "C", "G", "T"]:
        subset = df[df["ref"] == ref_base]
        if subset.empty:
            continue
        total_bases = float(subset["reads_all"].sum())
        total_sub   = float(subset["mismatches"].sum())
        total_ins   = float(subset["insertions"].sum())
        total_del   = float(subset["deletions"].sum())
        rows.append({
            "label":         label,
            "ref_base":      ref_base,
            "n_sites":       len(subset),
            "total_bases":   int(total_bases),
            "sub_rate":      total_sub / total_bases,
            "ins_rate":      total_ins / total_bases,
            "del_rate":      total_del / total_bases,
            "total_rate":    (total_sub + total_ins + total_del) / total_bases,
        })

    return pd.DataFrame(rows)


# =============================================================================
# Step 5 — Decile-based coverage profile
# Bin sites into deciles by coverage and compute mean error rate per decile.
# This gives a smooth coverage vs error curve without fixed bin boundaries.
# =============================================================================

def coverage_decile_profile(df, label, n_deciles=20):
    """
    Divide positions into n_deciles equal-count bins by coverage depth
    and compute weighted aggregate error rate per bin.
    This smooths over the irregular distribution of coverage values and
    gives a cleaner curve for plotting.
    """
    log(f"Computing coverage decile profile ({n_deciles} bins)...")

    df = df[df["reads_all"] > 0].copy()
    df["cov_decile"] = pd.qcut(
        df["reads_all"], q=n_deciles, duplicates="drop"
    )

    rows = []
    for interval, grp in df.groupby("cov_decile", observed=True):
        total_bases = float(grp["reads_all"].sum())
        if total_bases == 0:
            continue
        rows.append({
            "label":        label,
            "cov_mid":      float(interval.mid),
            "cov_min":      float(interval.left),
            "cov_max":      float(interval.right),
            "n_sites":      len(grp),
            "total_bases":  int(total_bases),
            "sub_rate":     grp["mismatches"].sum() / total_bases,
            "ins_rate":     grp["insertions"].sum() / total_bases,
            "del_rate":     grp["deletions"].sum()  / total_bases,
            "total_rate":   (grp["mismatches"].sum() +
                             grp["insertions"].sum() +
                             grp["deletions"].sum()) / total_bases,
            # Per-position IQR within this decile bin — used for error bars
            # in the smooth decile panel (Fig 1C).
            # Computed on per-position rates not weighted rates, so they
            # reflect the actual spread of positions in the bin.
            "sub_p25":   float(grp["sub_rate"].quantile(0.25))   if "sub_rate"   in grp.columns else np.nan,
            "sub_p75":   float(grp["sub_rate"].quantile(0.75))   if "sub_rate"   in grp.columns else np.nan,
            "del_p25":   float(grp["del_rate"].quantile(0.25))   if "del_rate"   in grp.columns else np.nan,
            "del_p75":   float(grp["del_rate"].quantile(0.75))   if "del_rate"   in grp.columns else np.nan,
            "total_p25": float(grp["total_rate"].quantile(0.25)) if "total_rate" in grp.columns else np.nan,
            "total_p75": float(grp["total_rate"].quantile(0.75)) if "total_rate" in grp.columns else np.nan,
        })

    return pd.DataFrame(rows)


# =============================================================================
# Save
# =============================================================================

def save_outputs(outdir, label, binned_df, decile_df, base_df):
    Path(outdir).mkdir(parents=True, exist_ok=True)

    binned_df.to_csv(
        f"{outdir}/{label}_coverage_binned_error_rates.tsv",
        sep="\t", index=False)
    log(f"  Saved: {label}_coverage_binned_error_rates.tsv")

    if decile_df is not None and not decile_df.empty:
        decile_df.to_csv(
            f"{outdir}/{label}_coverage_decile_error_rates.tsv",
            sep="\t", index=False)
        log(f"  Saved: {label}_coverage_decile_error_rates.tsv")

    if base_df is not None and not base_df.empty:
        base_df.to_csv(
            f"{outdir}/{label}_per_base_error_rates.tsv",
            sep="\t", index=False)
        log(f"  Saved: {label}_per_base_error_rates.tsv")


# =============================================================================
# Argument parsing
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Coverage vs error rate analysis from pysamstats output"
    )
    p.add_argument("--pysamstats",   required=True,
                   help="pysamstats --type variation TSV output")
    p.add_argument("--giab-bed",     default=None,
                   help="GIAB high-confidence BED for position filtering. "
                        "Only use with --chroms to avoid OOM. "
                        "If not provided, all positions are used.")
    p.add_argument("--outdir",       default="results/coverage_error")
    p.add_argument("--sample-label", default="sample")
    p.add_argument("--chroms",       default=None,
                   help="Comma-separated chromosomes to analyse. "
                        "Strongly recommended when using --giab-bed. "
                        "Default: all autosomes + chrX")
    p.add_argument("--n-deciles",    type=int, default=20,
                   help="Number of coverage decile bins for smooth curve "
                        "(default: 20)")
    return p.parse_args()


# =============================================================================
# Main
# =============================================================================

def main():
    args   = parse_args()
    chroms = args.chroms.split(",") if args.chroms else AUTOSOMES_CHRX

    if not os.path.exists(args.pysamstats):
        err(f"pysamstats file not found: {args.pysamstats}")

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # Load GIAB set (only if chromosomes are restricted)
    giab_set = None
    if args.giab_bed:
        if args.chroms is None:
            warn("--giab-bed provided without --chroms. This will load the "
                 "entire GIAB position set into RAM and may OOM.")
            warn("Recommend adding --chroms chr19 or similar.")
        giab_set = load_giab_bed_set(args.giab_bed, chroms)

    # Load data
    df = load_pysamstats_chunked(args.pysamstats, chroms, giab_set)

    # Analysis
    binned_df = bin_by_coverage(df, args.sample_label)
    decile_df = coverage_decile_profile(df, args.sample_label, args.n_deciles)
    base_df   = stratify_by_ref_base(df, args.sample_label)

    # Save
    save_outputs(args.outdir, args.sample_label, binned_df, decile_df, base_df)
    log("Done.")


if __name__ == "__main__":
    main()
