#!/usr/bin/env python3
"""
Apply final site-level filters to a native_vs_ivt_fisher BEDRMod file.

Two classes of sites are retained:

  Testable sites (padj != NA):
    (padj <= --max-padj AND log2(native_freq / ivt_freq) > --min-log2fc)
    OR score > --min-ivt-absent-score
    AND native_freq > --min-native-freq

  IVT-absent sites (padj == NA):
    score > --min-ivt-absent-score
    native_freq > --min-native-freq

All header lines are preserved in the output.
Streams the input in chunks — constant memory regardless of file size.
"""

import ast
import gzip
import sys
import argparse
import numpy as np
import pandas as pd


CHUNK_SIZE = 500_000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",  help="Input BEDRMod file (plain or .gz)")
    p.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    p.add_argument("--min-log2fc",           type=float, required=True,
                   metavar="N", help="Minimum log2(native/IVT freq) for testable sites")
    p.add_argument("--max-padj",             type=float, required=True,
                   metavar="N", help="Maximum adjusted p-value for testable sites")
    p.add_argument("--min-ivt-absent-score", type=float, required=True,
                   metavar="N", help="Minimum score for IVT-absent sites")
    p.add_argument("--min-native-freq",      type=float, required=True,
                   metavar="N", help="Minimum native frequency (same units as column)")
    p.add_argument("--output-counts", default=None, metavar="FILE",
                   help="Write per-modification filter step counts to this TSV")
    return p.parse_args()


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


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
        raise ValueError(f"Could not detect column names from header of {path}")
    return col_names, n_skip


def _filter_chunk(chunk, args, count_mods=False):
    freq     = pd.to_numeric(chunk["frequency"],     errors="coerce")
    ivt_freq = pd.to_numeric(chunk["ivt_frequency"], errors="coerce")
    padj     = pd.to_numeric(chunk["padj"],          errors="coerce")
    score    = pd.to_numeric(chunk["score"],         errors="coerce")

    is_testable = padj.notna()

    with np.errstate(divide="ignore", invalid="ignore"):
        log2fc = np.where(
            is_testable & (freq > 0) & (ivt_freq > 0),
            np.log2(freq.values / ivt_freq.values),
            np.nan,
        )

    padj_ok    = is_testable & (padj <= args.max_padj)
    padj_fc_ok = padj_ok & (log2fc > args.min_log2fc)
    or_ok      = is_testable & (padj_fc_ok | (score > args.min_ivt_absent_score))
    abs_sc_ok  = ~is_testable & (score > args.min_ivt_absent_score)

    pass_testable = or_ok    & (freq > args.min_native_freq)
    pass_absent   = abs_sc_ok & (freq > args.min_native_freq)
    keep          = pass_testable | pass_absent

    mod_counts = None
    if count_mods:
        mod      = chunk["name"].astype(str)
        mod_counts = {}
        for m in mod.unique():
            mask = mod == m
            mod_counts[m] = {
                "total":            int(mask.sum()),
                "ivt_testable":     int((mask & is_testable).sum()),
                "ivt_absent":       int((mask & ~is_testable).sum()),
                "padj":             int((mask & padj_ok).sum()),
                "padj_fc":          int((mask & padj_fc_ok).sum()),
                "padj_fc_freq":     int((mask & padj_fc_ok & (freq > args.min_native_freq)).sum()),
                "or_condition":     int((mask & or_ok).sum()),
                "absent_score":     int((mask & abs_sc_ok).sum()),
                "or_condition_freq":  int((mask & pass_testable).sum()),
                "absent_score_freq":  int((mask & pass_absent).sum()),
                "after_filter":     int((mask & keep).sum()),
            }

    return (
        chunk[keep],
        mod_counts,
        int(len(chunk)),
        int(is_testable.sum()),
        int(pass_testable.sum()),
        int(pass_absent.sum()),
    )


def _write_counts_tsv(all_mod_counts, args):
    mp = args.max_padj
    lf = args.min_log2fc
    sc = args.min_ivt_absent_score

    fr = args.min_native_freq
    col_map = [
        ("total",              "total"),
        ("ivt_testable",       "ivt_testable"),
        ("ivt_absent",         "ivt_absent"),
        ("padj",               f"padj_le_{mp}"),
        ("padj_fc",            f"padj_le_{mp}_and_log2fc_gt_{lf}"),
        ("padj_fc_freq",       f"padj_le_{mp}_and_log2fc_gt_{lf}_freq_gt_{fr}"),
        ("or_condition",       f"padj_le_{mp}_log2fc_gt_{lf}_or_score_gt_{sc}"),
        ("absent_score",       f"absent_and_score_gt_{sc}"),
        ("or_condition_freq",  f"padj_le_{mp}_log2fc_gt_{lf}_or_score_gt_{sc}_freq_gt_{fr}"),
        ("absent_score_freq",  f"absent_and_score_gt_{sc}_freq_gt_{fr}"),
        ("after_filter",       "after_filter"),
    ]

    with open(args.output_counts, "w") as fh:
        fh.write(
            f"#max_padj={mp}\tmin_log2fc={lf}\t"
            f"min_ivt_absent_score={sc}\tmin_native_freq={args.min_native_freq}\n"
        )
        headers = ["modification"] + [label for _, label in col_map]
        fh.write("\t".join(headers) + "\n")
        for mod in sorted(all_mod_counts):
            row = [mod] + [
                str(all_mod_counts[mod].get(key, 0)) for key, _ in col_map
            ]
            fh.write("\t".join(row) + "\n")


def main():
    args = parse_args()

    col_names, n_skip = detect_header(args.input)
    for col in ("name", "frequency", "ivt_frequency", "score", "padj"):
        if col not in col_names:
            raise ValueError(f"Column '{col}' not found in {args.input}.")

    print(f"[info] Filtering {args.input} ...", file=sys.stderr)

    out = open(args.output, "w") if args.output != "-" else sys.stdout
    try:
        with _open(args.input) as fh:
            for line in fh:
                if line.startswith("#"):
                    out.write(line)
                else:
                    break

        n_total = n_test = n_kt = n_ka = 0
        count_mods  = args.output_counts is not None
        all_mod_counts: dict = {}

        for chunk in pd.read_csv(
            args.input, sep="\t", skiprows=n_skip, header=None,
            names=col_names, dtype=str, chunksize=CHUNK_SIZE,
        ):
            filtered, mod_counts, nt, ns, kt, ka = _filter_chunk(
                chunk, args, count_mods=count_mods,
            )
            n_total += nt
            n_test  += ns
            n_kt    += kt
            n_ka    += ka
            filtered.to_csv(out, sep="\t", header=False, index=False)

            if mod_counts:
                for m, c in mod_counts.items():
                    if m not in all_mod_counts:
                        all_mod_counts[m] = {k: 0 for k in c}
                    for k, v in c.items():
                        all_mod_counts[m][k] += v

    finally:
        if args.output != "-":
            out.close()

    if count_mods:
        _write_counts_tsv(all_mod_counts, args)
        print(f"[info] Counts written → {args.output_counts}", file=sys.stderr)

    print(
        f"[info] Total: {n_total:,} | Testable: {n_test:,} "
        f"| IVT-absent: {n_total - n_test:,}",
        file=sys.stderr,
    )
    print(
        f"[info] Kept: {n_kt + n_ka:,} "
        f"({n_kt:,} testable + {n_ka:,} IVT-absent)",
        file=sys.stderr,
    )
    print("[info] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
