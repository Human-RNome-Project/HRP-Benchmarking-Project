#!/usr/bin/env python3
"""Call m5C sites from UBS-seq per-position conversion tables.

Reproduces `05_m5C_calling.ipynb`.

At every cytosine with depth *d* and *k* unconverted reads the binomial tail
P(X >= k) ~ Binomial(d, e) is evaluated, where *e* is the per-sample background
non-conversion rate (Legrand et al. 2017; Dai et al. 2024).  Sites passing the
depth, stoichiometry and p-value thresholds are written out.

Two modes:

  calling     m5C_calling.py --input conv_unconv3n_filter/YSL-5.tsv.gz \
                             --sample YSL-5 --error-rate 0.000731 \
                             --min-depth 3 --min-ratio 0.0 --max-p-value 1e-12 \
                             --output m5C_sites_YSL-5_dpth_3_p_val_1e-12_min_rat_0.0.tsv

  intersect   m5C_calling.py --intersect --input A.tsv B.tsv --output shared.tsv

Requires pandas, numpy and scipy.
"""

import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import binom

OUT_COLS = ["Sample", "Chrom", "Pos", "Strand", "Motif", "Converted",
            "Unconverted", "Depth", "Ratio_conv", "Ratio_unconv", "p_val", "ID"]

DTYPES = {"Sample": "string", "Chrom": "string", "Pos": "int64",
          "Strand": "string", "Converted": "int64", "Unconverted": "int64",
          "Depth": "int64", "Motif": "string"}


def estimate_error_rate(path, chunksize):
    """Background non-conversion rate = sum(Unconverted) / sum(Depth)."""
    unconv = depth = 0
    for chunk in pd.read_csv(path, sep="\t", chunksize=chunksize,
                             usecols=["Unconverted", "Depth"]):
        unconv += int(chunk["Unconverted"].sum())
        depth += int(chunk["Depth"].sum())
    if depth == 0:
        sys.exit(f"error: no covered positions in {path}")
    return unconv / depth


def call_sites(args):
    error_rate = args.error_rate
    if error_rate is None:
        source = args.error_rate_from or args.input[0]
        error_rate = estimate_error_rate(source, args.chunksize)
        sys.stderr.write(f"estimated background non-conversion rate from "
                         f"{source}: {error_rate:.6g}\n")
    else:
        sys.stderr.write(f"background non-conversion rate: {error_rate:.6g}\n")

    kept = []
    total = 0
    for chunk in pd.read_csv(args.input[0], sep="\t", chunksize=args.chunksize,
                             dtype={k: v for k, v in DTYPES.items()}):
        total += len(chunk)
        if "Motif" not in chunk.columns:
            chunk["Motif"] = pd.NA
        sub = chunk[(chunk["Depth"] >= args.min_depth) &
                    (chunk["Unconverted"] > 0)]
        if args.min_ratio > 0:
            sub = sub[sub["Unconverted"] / sub["Depth"] >= args.min_ratio]
        if sub.empty:
            continue
        sub = sub.copy()
        sub["p_val"] = binom.sf(sub["Unconverted"].to_numpy() - 1,
                                sub["Depth"].to_numpy(), error_rate)
        sub = sub[sub["p_val"] <= args.max_p_value]
        if not sub.empty:
            kept.append(sub)

    if kept:
        df = pd.concat(kept, ignore_index=True)
    else:
        df = pd.DataFrame(columns=OUT_COLS)

    if len(df):
        if args.sample:
            df["Sample"] = args.sample
        df["Ratio_conv"] = df["Converted"] / df["Depth"]
        df["Ratio_unconv"] = df["Unconverted"] / df["Depth"]
        df["ID"] = df["Chrom"].astype(str) + "_" + df["Pos"].astype(str)
        df = df.sort_values(["Chrom", "Pos", "Strand"], kind="mergesort")
        df = df[OUT_COLS]

    df.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
    sys.stderr.write(f"{len(df)} of {total} positions called; "
                     f"written to {args.output}\n")


def intersect(args):
    if len(args.input) < 2:
        sys.exit("error: --intersect needs at least two input files")
    frames = [pd.read_csv(f, sep="\t") for f in args.input]
    shared = set(frames[0]["ID"])
    for fr in frames[1:]:
        shared &= set(fr["ID"])
    out = pd.concat([fr[fr["ID"].isin(shared)] for fr in frames],
                    ignore_index=True)
    out = out.sort_values(["Chrom", "Pos", "Strand", "Sample"], kind="mergesort")
    out.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
    sys.stderr.write(f"{len(shared)} shared sites "
                     f"({len(out)} site x sample rows) written to {args.output}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", nargs="+", required=True,
                   help="per-position table (calling) or called-site tables (--intersect)")
    p.add_argument("--output", required=True)
    p.add_argument("--sample", help="sample ID written into the Sample column")
    p.add_argument("--error-rate", type=float,
                   help="background non-conversion rate; if omitted it is "
                        "estimated from --error-rate-from or from --input")
    p.add_argument("--error-rate-from",
                   help="unfiltered per-position table used to estimate the "
                        "background non-conversion rate (recommended)")
    p.add_argument("--min-depth", type=int, default=3)
    p.add_argument("--min-ratio", type=float, default=0.0)
    p.add_argument("--max-p-value", type=float, default=1e-12)
    p.add_argument("--chunksize", type=int, default=2_000_000)
    p.add_argument("--intersect", action="store_true",
                   help="intersect called-site tables instead of calling")
    args = p.parse_args()

    if args.intersect:
        intersect(args)
    else:
        call_sites(args)


if __name__ == "__main__":
    main()
