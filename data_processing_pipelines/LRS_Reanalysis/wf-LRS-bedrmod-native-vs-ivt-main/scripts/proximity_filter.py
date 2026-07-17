#!/usr/bin/env python3
"""
Remove modification calls too close to a higher-scored neighbour.

A site is discarded when another site on the same (chrom, strand) at a
DIFFERENT chromStart within --max-proximity nt has a strictly higher score.
Sites at the same chromStart (different modification types at the same
position) are never compared against each other.

The coordinate order of the input is preserved in the output.
The input must be sorted by (chrom, chromStart).
"""

import ast
import gzip
import sys
import argparse
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",  help="Input BEDRMod file (plain or .gz)")
    p.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    p.add_argument("--max-proximity", type=int, required=True, metavar="N",
                   help="Remove sites within N nt of a higher-scored site at a "
                        "different coordinate on the same (chrom, strand)")
    p.add_argument("--output-counts", default=None, metavar="FILE",
                   help="Write per-modification before/after dedup counts to this TSV")
    return p.parse_args()


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def detect_header(path):
    col_names = None
    n_skip = 0
    with _open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            n_skip += 1
            stripped = line.lstrip("#").strip()
            if stripped.startswith("["):
                try:
                    col_names = ast.literal_eval(stripped)
                except (ValueError, SyntaxError):
                    pass
    if col_names is None:
        raise ValueError(f"No column header found in {path}")
    return col_names, n_skip


def proximity_mask(df, max_dist):
    """
    Return a boolean numpy array; True = keep the row.

    A row is removed if any other row on the same (chrom, strand) at a
    DIFFERENT chromStart within max_dist has a strictly higher score.
    Rows at the same chromStart are never compared.

    df must be sorted by (chrom, chromStart).
    """
    chrom  = df["chrom"].to_numpy(dtype=str)
    start  = df["chromStart"].astype(int).to_numpy()
    strand = df["strand"].to_numpy(dtype=str)
    score  = pd.to_numeric(df["score"], errors="coerce").to_numpy(dtype=float)

    n    = len(df)
    keep = np.ones(n, dtype=bool)

    for i in range(n):
        if not keep[i]:
            continue
        c_i, p_i, s_i, sc_i = chrom[i], start[i], strand[i], score[i]

        # Scan forward (larger coordinates, same chrom)
        for j in range(i + 1, n):
            if chrom[j] != c_i:
                break                          # new chrom — stop
            if start[j] - p_i > max_dist:
                break                          # beyond window — stop
            if strand[j] != s_i or start[j] == p_i:
                continue                       # different strand or same pos — skip
            if score[j] > sc_i:
                keep[i] = False
                break

        if not keep[i]:
            continue

        # Scan backward (smaller coordinates, same chrom)
        for j in range(i - 1, -1, -1):
            if chrom[j] != c_i:
                break
            if p_i - start[j] > max_dist:
                break
            if strand[j] != s_i or start[j] == p_i:
                continue
            if score[j] > sc_i:
                keep[i] = False
                break

    return keep


def main():
    args = parse_args()

    col_names, n_skip = detect_header(args.input)
    for col in ("chrom", "chromStart", "strand", "score"):
        if col not in col_names:
            raise ValueError(f"Column '{col}' not found in {args.input}.")

    print(f"[info] Reading {args.input} ...", file=sys.stderr)
    df = pd.read_csv(
        args.input, sep="\t", skiprows=n_skip, header=None,
        names=col_names, dtype=str,
    )
    print(f"[info] {len(df):,} sites loaded", file=sys.stderr)

    mask    = proximity_mask(df, args.max_proximity)
    n_kept  = int(mask.sum())
    print(
        f"[info] After proximity filter (≤{args.max_proximity} nt): "
        f"{n_kept:,} kept, {len(df) - n_kept:,} removed",
        file=sys.stderr,
    )

    out = open(args.output, "w") if args.output != "-" else sys.stdout
    try:
        with _open(args.input) as fh:
            for line in fh:
                if line.startswith("#"):
                    out.write(line)
                else:
                    break
        df[mask].to_csv(out, sep="\t", header=False, index=False)
    finally:
        if args.output != "-":
            out.close()

    if args.output_counts is not None:
        if "name" not in col_names:
            raise ValueError(f"Column 'name' not found in {args.input}; required for --output-counts")
        mod = df["name"].astype(str)
        mods = sorted(mod.unique())
        with open(args.output_counts, "w") as fh:
            fh.write(f"#max_proximity={args.max_proximity}\n")
            fh.write("modification\tbefore_dedup\tafter_dedup\n")
            for m in mods:
                m_mask = mod == m
                before = int(m_mask.sum())
                after  = int((m_mask & mask).sum())
                fh.write(f"{m}\t{before}\t{after}\n")
        print(f"[info] Counts written → {args.output_counts}", file=sys.stderr)

    print("[info] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
