#!/usr/bin/env python3
"""
calculate_read_error_rates_parallel.py

Per-read error rate calculation from BAM with two execution modes for
parallel WDL execution:

  --mode chrom  (default)
      Single-pass reservoir sampling over one chromosome of a BAM.
      Outputs a per-read TSV and raw spectrum/positional bin arrays
      as numpy .npz files for merging.
      sample_n reads are drawn from this chromosome proportionally.

  --mode merge
      Takes per-chrom per-read TSVs and .npz files, concatenates them,
      resamples down to --sample-n if needed, and produces final outputs:
        - per_read_error_rates.tsv
        - read_error_rate_summary.tsv
        - read_substitution_spectrum.tsv
        - read_positional_bias.tsv

WDL execution model:
  scatter(chrom) {
    call ChromTask   → per_read.tsv + counts.npz
  }
  call MergeTask   → all final outputs

T is displayed as U throughout (RNA convention).

Usage — chrom mode:
    python calculate_read_error_rates_parallel.py \\
        --mode         chrom \\
        --bam          merged_native.bam \\
        --chrom        chr1 \\
        --sample-n     2274 \\
        --min-mapq     20 \\
        --min-length   200 \\
        --sample-label GM12878_native \\
        --seed         42 \\
        --outdir       results/chr1/

Usage — merge mode:
    python calculate_read_error_rates_parallel.py \\
        --mode         merge \\
        --per-read-tsvs chr1/GM12878_native_chr1_per_read.tsv ... \\
        --npz-files     chr1/GM12878_native_chr1_counts.npz ... \\
        --sample-n      50000 \\
        --sample-label  GM12878_native \\
        --seed          42 \\
        --outdir        results/merged/
"""

import argparse
import os
import random
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import pysam
except ImportError:
    print("[ERROR] pysam not installed. Run: pip install pysam")
    sys.exit(1)


# =============================================================================
# Constants
# =============================================================================

AUTOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

REF_BASES   = ["A", "C", "G", "T"]
SUB_TYPES   = [f"{r}>{a}" for r in REF_BASES for a in REF_BASES if r != a]
SUB_DISPLAY = {s: s.replace("T>", "U>").replace(">T", ">U") for s in SUB_TYPES}

N_POS_BINS    = 20
POS_BIN_EDGES = np.linspace(0, 1, N_POS_BINS + 1)
POS_BIN_MIDS  = (POS_BIN_EDGES[:-1] + POS_BIN_EDGES[1:]) / 2


# =============================================================================
# Utilities
# =============================================================================

def log(msg):  print(f"[INFO]  {msg}", flush=True)
def warn(msg): print(f"[WARN]  {msg}", flush=True)
def err(msg):
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


# =============================================================================
# MD tag and CIGAR parsers
# =============================================================================

def get_mismatch_read_positions(read):
    if read.cigartuples is None:
        return []
    try:
        md = read.get_tag("MD")
    except KeyError:
        return []
    seq = read.query_sequence
    if not seq:
        return []

    subs_at_offset = {}
    offset = 0
    for token in re.findall(r'(\d+|\^[ACGT]+|[ACGT])', md):
        if token.isdigit():
            offset += int(token)
        elif token.startswith("^"):
            pass
        else:
            subs_at_offset[offset] = token.upper()
            offset += 1

    sub_records = []
    cigar_qpos  = 0
    md_offset   = 0
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            for _ in range(length):
                if md_offset in subs_at_offset:
                    ref_b  = subs_at_offset[md_offset]
                    read_b = seq[cigar_qpos].upper()
                    if read_b != ref_b:
                        sub_records.append((cigar_qpos, ref_b, read_b))
                cigar_qpos += 1
                md_offset  += 1
        elif op == 1:  cigar_qpos += length
        elif op == 2:  md_offset  += length
        elif op == 4:  cigar_qpos += length
    return sub_records


def get_indel_read_positions(cigartuples, qlen):
    if cigartuples is None or qlen == 0:
        return [], []
    ins_pos, del_pos, qpos = [], [], 0
    for op, length in cigartuples:
        if op in (0, 7, 8):  qpos += length
        elif op == 1:
            frac = qpos / qlen
            ins_pos.extend([frac] * length)
            qpos += length
        elif op == 2:
            frac = qpos / qlen
            del_pos.extend([frac] * length)
        elif op == 4:
            qpos += length
    return ins_pos, del_pos


# =============================================================================
# Extract stats from one read
# =============================================================================

def extract_one_read(read):
    if read.cigartuples is None:
        return None
    qlen = read.query_length
    if qlen is None or qlen == 0:
        return None

    n_mismatches = n_insertions = n_deletions = n_soft_clip = aligned_len = 0
    for op, length in read.cigartuples:
        if   op == 0:  aligned_len  += length
        elif op == 1:  n_insertions += length
        elif op == 2:  n_deletions  += length; aligned_len += length
        elif op == 4:  n_soft_clip  += length
        elif op == 7:  aligned_len  += length
        elif op == 8:  n_mismatches += length; aligned_len += length

    has_xeq = any(op in (7, 8) for op, _ in read.cigartuples)
    if not has_xeq:
        try:
            nm = read.get_tag("NM")
            n_mismatches = max(0, nm - n_insertions - n_deletions)
        except KeyError:
            return None

    if aligned_len == 0:
        return None

    sub_rate   = n_mismatches / aligned_len
    ins_rate   = n_insertions / aligned_len
    del_rate   = n_deletions  / aligned_len
    total_rate = (n_mismatches + n_insertions + n_deletions) / aligned_len

    bq      = read.query_qualities
    mean_bq = float(np.mean(bq)) if bq is not None else np.nan
    seq     = read.query_sequence or ""
    gc      = (seq.count("G") + seq.count("C")) / len(seq) if seq else np.nan

    record = {
        "read_name":         read.query_name,
        "chrom":             read.reference_name,
        "start":             read.reference_start,
        "end":               read.reference_end,
        "strand":            "-" if read.is_reverse else "+",
        "mapq":              read.mapping_quality,
        "query_length":      qlen,
        "aligned_length":    aligned_len,
        "n_mismatches":      n_mismatches,
        "n_insertions":      n_insertions,
        "n_deletions":       n_deletions,
        "n_soft_clip":       n_soft_clip,
        "sub_rate":          sub_rate,
        "ins_rate":          ins_rate,
        "del_rate":          del_rate,
        "total_rate":        total_rate,
        "mean_base_quality": mean_bq,
        "gc_content":        gc,
    }

    sub_records          = get_mismatch_read_positions(read)
    ins_fracs, del_fracs = get_indel_read_positions(read.cigartuples, qlen)
    return record, sub_records, ins_fracs, del_fracs


# =============================================================================
# Single-pass reservoir sampling over one chromosome
# =============================================================================

def sample_chrom(bam_path, chrom, sample_n, min_mapq, min_length, seed):
    log(f"Sampling {chrom}: n={sample_n:,}, MAPQ>={min_mapq}, len>={min_length}")
    random.seed(seed)

    reservoir  = []
    n_seen     = 0
    n_filtered = 0
    n_no_md    = 0

    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch(chrom):
        if (read.is_secondary or read.is_supplementary or
                read.is_unmapped or
                read.mapping_quality < min_mapq or
                read.query_length is None or
                read.query_length < min_length):
            n_filtered += 1
            continue

        result = extract_one_read(read)
        if result is None:
            n_filtered += 1
            continue

        if not result[1]:
            n_no_md += 1

        n_seen += 1
        if len(reservoir) < sample_n:
            reservoir.append(result)
        else:
            j = random.randint(0, n_seen - 1)
            if j < sample_n:
                reservoir[j] = result

    bam.close()
    log(f"  {chrom}: {n_seen:,} qualifying reads, {len(reservoir):,} sampled")
    if n_no_md > 0:
        warn(f"  {n_no_md:,} reads had no MD tag on {chrom}")
    return reservoir


# =============================================================================
# Build per-chrom outputs from reservoir
# Returns (df, sub_counts dict, ref_counts dict, pos arrays)
# =============================================================================

def build_chrom_outputs(reservoir):
    records       = []
    sub_counts    = {s: 0 for s in SUB_TYPES}
    ref_counts    = {b: 0 for b in REF_BASES}
    pos_sub_bins  = np.zeros(N_POS_BINS)
    pos_ins_bins  = np.zeros(N_POS_BINS)
    pos_del_bins  = np.zeros(N_POS_BINS)
    pos_base_bins = np.zeros(N_POS_BINS)

    for record, sub_records, ins_fracs, del_fracs in reservoir:
        records.append(record)
        qlen        = record["query_length"]
        aligned_len = record["aligned_length"]

        for _, ref_base, read_base in sub_records:
            key = f"{ref_base}>{read_base}"
            if key in sub_counts:
                sub_counts[key] += 1
            if ref_base in ref_counts:
                ref_counts[ref_base] += 1

        for qpos, _, _ in sub_records:
            if qlen > 0:
                bin_idx = min(int(qpos / qlen * N_POS_BINS), N_POS_BINS - 1)
                pos_sub_bins[bin_idx] += 1

        for frac in ins_fracs:
            bin_idx = min(int(frac * N_POS_BINS), N_POS_BINS - 1)
            pos_ins_bins[bin_idx] += 1

        for frac in del_fracs:
            bin_idx = min(int(frac * N_POS_BINS), N_POS_BINS - 1)
            pos_del_bins[bin_idx] += 1

        pos_base_bins += aligned_len / N_POS_BINS

    df = pd.DataFrame(records)
    return df, sub_counts, ref_counts, pos_sub_bins, pos_ins_bins, pos_del_bins, pos_base_bins


# =============================================================================
# Summary statistics
# =============================================================================

def compute_summary(df, label):
    if df.empty:
        return {}
    summary = {"label": label, "n_reads": len(df)}
    for rate_col in ["sub_rate", "ins_rate", "del_rate", "total_rate"]:
        vals   = df[rate_col].dropna()
        prefix = rate_col.replace("_rate", "")
        if vals.empty:
            continue
        summary.update({
            f"{prefix}_mean":   float(vals.mean()),
            f"{prefix}_median": float(vals.median()),
            f"{prefix}_std":    float(vals.std()),
            f"{prefix}_p5":     float(vals.quantile(0.05)),
            f"{prefix}_p25":    float(vals.quantile(0.25)),
            f"{prefix}_p75":    float(vals.quantile(0.75)),
            f"{prefix}_p95":    float(vals.quantile(0.95)),
            f"{prefix}_p99":    float(vals.quantile(0.99)),
        })
    lens = df["query_length"].dropna()
    if not lens.empty:
        sl  = np.sort(lens.values)[::-1]
        cs  = np.cumsum(sl)
        idx = np.searchsorted(cs, cs[-1] / 2)
        summary.update({
            "read_length_mean":   float(lens.mean()),
            "read_length_median": float(lens.median()),
            "read_length_n50":    float(sl[idx]),
        })
    bq = df["mean_base_quality"].dropna()
    if not bq.empty:
        summary["mean_base_quality_mean"]   = float(bq.mean())
        summary["mean_base_quality_median"] = float(bq.median())
    return summary


def print_summary(summary, spec_df=None):
    label = summary.get("label", "")
    print(f"\n{'='*65}")
    print(f"  PER-READ ERROR RATE SUMMARY — {label}")
    print(f"{'='*65}")
    print(f"  Reads analysed: {summary.get('n_reads', 0):,}")
    print(f"\n  {'Metric':<28} {'Mean':>9}  {'Median':>9}  {'P5':>7}  {'P95':>7}")
    print(f"  {'-'*60}")
    for prefix, display in [("sub","Substitution rate"),("ins","Insertion rate"),
                             ("del","Deletion rate"),("total","Total error rate")]:
        m   = summary.get(f"{prefix}_mean",   np.nan)
        md  = summary.get(f"{prefix}_median", np.nan)
        p5  = summary.get(f"{prefix}_p5",     np.nan)
        p95 = summary.get(f"{prefix}_p95",    np.nan)
        print(f"  {display:<28} {m*100:>8.4f}%  {md*100:>8.4f}%  "
              f"{p5*100:>6.4f}%  {p95*100:>6.4f}%")
    print(f"\n  Read length N50: {summary.get('read_length_n50',0):.0f}bp  "
          f"median: {summary.get('read_length_median',0):.0f}bp")
    if spec_df is not None and not spec_df.empty:
        print(f"\n  Top substitutions (T displayed as U):")
        for _, row in spec_df.head(6).iterrows():
            print(f"    {row['display']:<6}  count={int(row['count']):>8,}  "
                  f"rate={row['rate']*100:.5f}%")
    print(f"{'='*65}")


def save_outputs(outdir, label, df, summary, spec_df, pos_df):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    df.to_csv(f"{outdir}/{label}_per_read_error_rates.tsv",    sep="\t", index=False)
    pd.DataFrame([summary]).to_csv(
        f"{outdir}/{label}_read_error_rate_summary.tsv",       sep="\t", index=False)
    if spec_df is not None and not spec_df.empty:
        spec_df.to_csv(
            f"{outdir}/{label}_read_substitution_spectrum.tsv", sep="\t", index=False)
    if pos_df is not None and not pos_df.empty:
        pos_df.to_csv(
            f"{outdir}/{label}_read_positional_bias.tsv",       sep="\t", index=False)
    log(f"Outputs written to {outdir}/")


def build_spec_df(sub_counts, ref_counts):
    rows = []
    for sub in SUB_TYPES:
        rb    = sub[0]
        denom = ref_counts.get(rb, 0)
        cnt   = sub_counts[sub]
        rows.append({
            "substitution":   sub,
            "display":        SUB_DISPLAY[sub],
            "ref_base":       rb,
            "alt_base":       sub[2],
            "count":          cnt,
            "ref_base_total": denom,
            "rate":           cnt / denom if denom > 0 else 0,
        })
    return pd.DataFrame(rows).sort_values("rate", ascending=False)


def build_pos_df(pos_sub_bins, pos_ins_bins, pos_del_bins, pos_base_bins):
    rows = []
    for i in range(N_POS_BINS):
        denom = pos_base_bins[i]
        rows.append({
            "bin_idx":       i,
            "bin_mid":       float(POS_BIN_MIDS[i]),
            "bin_start":     float(POS_BIN_EDGES[i]),
            "bin_end":       float(POS_BIN_EDGES[i + 1]),
            "sub_count":     int(pos_sub_bins[i]),
            "ins_count":     int(pos_ins_bins[i]),
            "del_count":     int(pos_del_bins[i]),
            "total_count":   int(pos_sub_bins[i] + pos_ins_bins[i] + pos_del_bins[i]),
            "aligned_bases": float(denom),
            "sub_rate":      pos_sub_bins[i] / denom if denom > 0 else 0,
            "ins_rate":      pos_ins_bins[i] / denom if denom > 0 else 0,
            "del_rate":      pos_del_bins[i] / denom if denom > 0 else 0,
            "total_rate":    (pos_sub_bins[i] + pos_ins_bins[i] +
                              pos_del_bins[i]) / denom if denom > 0 else 0,
        })
    return pd.DataFrame(rows)


# =============================================================================
# Argument parsing
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Per-read DRS error rates — chrom or merge mode"
    )
    p.add_argument("--mode", choices=["chrom", "merge"], default="chrom")

    # chrom mode
    cg = p.add_argument_group("chrom mode")
    cg.add_argument("--bam",    help="BAM file (must be indexed)")
    cg.add_argument("--chrom",  help="Single chromosome to process (e.g. chr1)")

    # merge mode
    mg = p.add_argument_group("merge mode")
    mg.add_argument("--per-read-tsvs", nargs="+",
                    help="Per-chrom per_read_error_rates TSV files")
    mg.add_argument("--npz-files", nargs="+",
                    help="Per-chrom counts .npz files")

    # shared
    p.add_argument("--outdir",       default="results/read_error_rates")
    p.add_argument("--sample-n",     type=int, default=50000,
                   help="Total reads to sample. In chrom mode: reads per chrom. "
                        "In merge mode: downsample merged df to this if larger.")
    p.add_argument("--min-mapq",     type=int, default=20)
    p.add_argument("--min-length",   type=int, default=200)
    p.add_argument("--sample-label", default="sample")
    p.add_argument("--seed",         type=int, default=42)
    return p.parse_args()


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # CHROM MODE
    # =========================================================================
    if args.mode == "chrom":
        if not args.bam:
            err("--bam required in chrom mode")
        if not args.chrom:
            err("--chrom required in chrom mode")
        if not os.path.exists(args.bam):
            err(f"BAM not found: {args.bam}")
        if (not os.path.exists(args.bam + ".bai") and
                not os.path.exists(args.bam + ".csi")):
            err(f"BAM index not found: {args.bam}.bai")

        reservoir = sample_chrom(
            args.bam, args.chrom, args.sample_n,
            args.min_mapq, args.min_length, args.seed
        )

        if not reservoir:
            warn(f"No reads sampled from {args.chrom} — writing empty outputs")
            # Write empty files so WDL output declarations are satisfied
            tag = f"{args.sample_label}_{args.chrom}"
            open(f"{args.outdir}/{tag}_per_read.tsv",    "w").close()
            np.savez(f"{args.outdir}/{tag}_counts.npz",
                     sub_counts=np.zeros(len(SUB_TYPES)),
                     ref_counts=np.zeros(len(REF_BASES)),
                     pos_sub_bins=np.zeros(N_POS_BINS),
                     pos_ins_bins=np.zeros(N_POS_BINS),
                     pos_del_bins=np.zeros(N_POS_BINS),
                     pos_base_bins=np.zeros(N_POS_BINS),
                     sub_type_keys=np.array(SUB_TYPES),
                     ref_base_keys=np.array(REF_BASES))
            log("Done (empty).")
            return

        df, sub_counts, ref_counts, \
            pos_sub_bins, pos_ins_bins, pos_del_bins, pos_base_bins = \
            build_chrom_outputs(reservoir)

        tag = f"{args.sample_label}_{args.chrom}"

        # Save per-read TSV
        df.to_csv(f"{args.outdir}/{tag}_per_read.tsv", sep="\t", index=False)
        log(f"Per-read TSV: {len(df):,} rows → {args.outdir}/{tag}_per_read.tsv")

        # Save raw counts as npz for merge mode
        npz_path = f"{args.outdir}/{tag}_counts.npz"
        np.savez(npz_path,
                 sub_counts   = np.array([sub_counts[s] for s in SUB_TYPES]),
                 ref_counts   = np.array([ref_counts[b] for b in REF_BASES]),
                 pos_sub_bins = pos_sub_bins,
                 pos_ins_bins = pos_ins_bins,
                 pos_del_bins = pos_del_bins,
                 pos_base_bins= pos_base_bins,
                 sub_type_keys= np.array(SUB_TYPES),
                 ref_base_keys= np.array(REF_BASES))
        log(f"Count arrays: {npz_path}")
        log(f"Chrom mode done: {args.chrom}")

    # =========================================================================
    # MERGE MODE
    # =========================================================================
    elif args.mode == "merge":
        if not args.per_read_tsvs:
            err("--per-read-tsvs required in merge mode")
        if not args.npz_files:
            err("--npz-files required in merge mode")

        # ── Concatenate per-read TSVs ─────────────────────────────────────────
        log(f"Merging {len(args.per_read_tsvs)} per-chrom TSVs...")
        dfs = []
        for f in args.per_read_tsvs:
            if not os.path.exists(f) or os.path.getsize(f) == 0:
                warn(f"Skipping empty/missing TSV: {f}")
                continue
            dfs.append(pd.read_csv(f, sep="\t", low_memory=False))

        if not dfs:
            err("No valid per-chrom TSVs found.")

        df_all = pd.concat(dfs, ignore_index=True)
        log(f"  Total reads after merge: {len(df_all):,}")

        # Downsample if more reads than sample_n
        if len(df_all) > args.sample_n:
            df_all = df_all.sample(n=args.sample_n,
                                   random_state=args.seed).reset_index(drop=True)
            log(f"  Downsampled to {args.sample_n:,} reads")

        # ── Aggregate npz count arrays ────────────────────────────────────────
        log(f"Merging {len(args.npz_files)} count .npz files...")
        sub_counts_agg    = {s: 0 for s in SUB_TYPES}
        ref_counts_agg    = {b: 0 for b in REF_BASES}
        pos_sub_bins_agg  = np.zeros(N_POS_BINS)
        pos_ins_bins_agg  = np.zeros(N_POS_BINS)
        pos_del_bins_agg  = np.zeros(N_POS_BINS)
        pos_base_bins_agg = np.zeros(N_POS_BINS)

        for f in args.npz_files:
            if not os.path.exists(f):
                warn(f"Skipping missing npz: {f}")
                continue
            npz = np.load(f, allow_pickle=False)
            keys = list(npz["sub_type_keys"])
            for i, s in enumerate(keys):
                if s in sub_counts_agg:
                    sub_counts_agg[s] += int(npz["sub_counts"][i])
            ref_keys = list(npz["ref_base_keys"])
            for i, b in enumerate(ref_keys):
                if b in ref_counts_agg:
                    ref_counts_agg[b] += int(npz["ref_counts"][i])
            pos_sub_bins_agg  += npz["pos_sub_bins"]
            pos_ins_bins_agg  += npz["pos_ins_bins"]
            pos_del_bins_agg  += npz["pos_del_bins"]
            pos_base_bins_agg += npz["pos_base_bins"]

        # ── Build final outputs ───────────────────────────────────────────────
        spec_df = build_spec_df(sub_counts_agg, ref_counts_agg)
        pos_df  = build_pos_df(pos_sub_bins_agg, pos_ins_bins_agg,
                                pos_del_bins_agg, pos_base_bins_agg)
        summary = compute_summary(df_all, args.sample_label)
        print_summary(summary, spec_df)
        save_outputs(args.outdir, args.sample_label,
                     df_all, summary, spec_df, pos_df)
        log("Merge mode done.")


if __name__ == "__main__":
    main()
