#!/usr/bin/env python3
"""
Annotate native BEDRMod sites with IVT comparison statistics.

For each site in the native file, runs a one-sided Fisher's exact test
(native > IVT) using n_mod vs count_canonical read counts. The test is
vectorized via the hypergeometric distribution (equivalent to Fisher's exact).
Sites absent from IVT (or below IVT coverage threshold) are retained with
NA p-values. BH FDR correction is applied over all testable sites.

Coverage filters:
  - Native sites with coverage <= --min-native-coverage are included in output
    but receive NA for pvalue, padj, and IVT fields (not tested).
  - IVT sites with coverage <= --min-ivt-coverage are treated as absent.

Memory strategy:
  - Load only the columns needed for the test with efficient dtypes.
  - After BH correction, stream the native file line-by-line for output.
"""

import ast
import gzip
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


KEY_COLS  = ["chrom", "chromStart", "chromEnd", "strand", "name"]
TEST_COLS = KEY_COLS + ["coverage", "n_mod", "count_canonical"]
IVT_COLS  = KEY_COLS + ["coverage", "frequency", "n_mod", "count_canonical"]

TEST_DTYPES = {
    "chrom":           "category",
    "chromStart":      "int32",
    "chromEnd":        "int32",
    "strand":          "category",
    "name":            "category",
    "coverage":        "int32",
    "n_mod":           "int32",
    "count_canonical": "int32",
}
IVT_DTYPES = {
    **TEST_DTYPES,
    "frequency": "float32",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("native", help="Native BEDRMod file")
    p.add_argument("ivt",    help="IVT BEDRMod file")
    p.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
    p.add_argument(
        "--min-native-coverage", type=int, default=20, metavar="N",
        help="Minimum coverage for native sites (default: 20)",
    )
    p.add_argument(
        "--min-ivt-coverage", type=int, default=10, metavar="N",
        help="Minimum coverage for IVT sites (default: 10)",
    )
    return p.parse_args()


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def detect_header(path):
    """Return (col_names, n_header_lines) from the BEDRMod #[...] line."""
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


def read_cols(path, usecols, col_names, n_skip, dtypes):
    return pd.read_csv(
        path,
        sep="\t",
        skiprows=n_skip,
        header=None,
        names=col_names,
        usecols=usecols,
        dtype=dtypes,
    )


def vectorized_fisher_greater(n_mod_nat, n_can_nat, n_mod_ivt, n_can_ivt):
    """
    One-sided Fisher's exact test (alternative='greater') via hypergeometric sf.
    p = P(X >= n_mod_nat), X ~ Hypergeom(N, K, n).
    """
    N = n_mod_nat + n_can_nat + n_mod_ivt + n_can_ivt
    K = n_mod_nat + n_mod_ivt
    n = n_mod_nat + n_can_nat
    return hypergeom.sf(n_mod_nat - 1, N, K, n)


def main():
    args = parse_args()
    min_nat_cov = args.min_native_coverage
    min_ivt_cov = args.min_ivt_coverage

    nat_cols, nat_skip = detect_header(args.native)
    ivt_cols, ivt_skip = detect_header(args.ivt)

    # --- Load IVT (minimal columns, efficient dtypes) ---
    print("[info] Loading IVT file...", file=sys.stderr)
    ivt_df = read_cols(args.ivt, IVT_COLS, ivt_cols, ivt_skip, IVT_DTYPES)
    ivt_df = ivt_df.drop_duplicates(subset=KEY_COLS)

    # Apply IVT coverage filter (below threshold → treated as absent)
    n_ivt_total = len(ivt_df)
    ivt_df = ivt_df[ivt_df["coverage"] > min_ivt_cov]
    print(
        f"[info] IVT sites: {n_ivt_total} → {len(ivt_df)} after coverage > {min_ivt_cov}",
        file=sys.stderr,
    )

    # IVT lookup dict for the streaming output pass
    ivt_lookup = {}
    for row in ivt_df.itertuples(index=False):
        k = (str(row.chrom), str(row.chromStart), str(row.chromEnd), str(row.strand), str(row.name))
        ivt_lookup[k] = (f"{row.frequency:.4g}", str(row.coverage), str(row.n_mod))

    ivt_test = ivt_df[KEY_COLS + ["n_mod", "count_canonical"]].rename(columns={
        "n_mod": "n_mod_ivt",
        "count_canonical": "count_canonical_ivt",
    })
    del ivt_df

    # --- Load native (minimal columns, efficient dtypes) ---
    print("[info] Loading native file (minimal columns)...", file=sys.stderr)
    nat_test = read_cols(args.native, TEST_COLS, nat_cols, nat_skip, TEST_DTYPES)
    n_nat_total = len(nat_test)

    # Mark which rows pass coverage (determines testability, not output inclusion)
    nat_pass = (nat_test["coverage"] > min_nat_cov).values
    n_passing = int(nat_pass.sum())
    print(
        f"[info] Native sites: {n_nat_total} total | {n_passing} pass coverage > {min_nat_cov} "
        f"| {n_nat_total - n_passing} below threshold (written with NA)",
        file=sys.stderr,
    )

    # --- Merge only coverage-passing rows with IVT ---
    print("[info] Merging...", file=sys.stderr)
    nat_passing = nat_test[nat_pass].drop(columns=["coverage"])
    del nat_test
    merged = nat_passing.merge(ivt_test, on=KEY_COLS, how="left")
    del nat_passing, ivt_test

    has_ivt    = merged["n_mod_ivt"].notna()
    n_testable = int(has_ivt.sum())
    print(
        f"[info] Coverage-passing sites: {len(merged)} | In IVT: {n_testable} | "
        f"IVT-absent (NA): {(~has_ivt).sum()}",
        file=sys.stderr,
    )

    # --- Vectorized Fisher's exact test ---
    # pvalues/padj span ALL native rows; only testable (passing + in IVT) rows are filled.
    print("[info] Running Fisher's exact tests...", file=sys.stderr)
    pvalues = np.full(n_nat_total, np.nan, dtype=np.float64)
    padj    = np.full(n_nat_total, np.nan, dtype=np.float64)

    if n_testable > 0:
        # Map from merged-row index to global (file-row) index
        global_pass_idx      = np.where(nat_pass)[0]
        testable_global_idx  = global_pass_idx[has_ivt.values]
        t = merged[has_ivt]
        pvals = vectorized_fisher_greater(
            t["n_mod"].values.astype(np.int64),
            t["count_canonical"].values.astype(np.int64),
            t["n_mod_ivt"].values.astype(np.int64),
            t["count_canonical_ivt"].values.astype(np.int64),
        )
        pvalues[testable_global_idx] = pvals
        _, padj_vals, _, _ = multipletests(pvals, method="fdr_bh")
        padj[testable_global_idx] = padj_vals

    del merged

    # --- Stream native file line-by-line for output ---
    print("[info] Writing output...", file=sys.stderr)
    out_cols    = nat_cols + ["ivt_frequency", "ivt_coverage", "ivt_n_mod", "pvalue", "padj"]
    out = open(args.output, "w") if args.output != "-" else sys.stdout

    try:
        orig_row_idx = 0
        written = 0
        col_written = False

        with open(args.native) as fh:
            for line in fh:
                if line.startswith("#"):
                    stripped = line.lstrip("#").strip()
                    is_col = False
                    if stripped.startswith("["):
                        try:
                            parsed = ast.literal_eval(stripped)
                            is_col = isinstance(parsed, list) and "chrom" in parsed
                        except (ValueError, SyntaxError):
                            pass
                    if is_col:
                        # Emit the column-definition line exactly once; other
                        # bracketed metadata lines pass through unchanged.
                        if not col_written:
                            out.write(f"#{out_cols}\n")
                            col_written = True
                    else:
                        out.write(line)
                    continue

                fields = line.rstrip("\n").split("\t")

                if nat_pass[orig_row_idx]:
                    key = (fields[0], fields[1], fields[2], fields[5], fields[3])
                    ivt_freq, ivt_cov, ivt_nmod = ivt_lookup.get(key, ("NA", "NA", "NA"))
                else:
                    ivt_freq = ivt_cov = ivt_nmod = "NA"

                pv  = f"{pvalues[orig_row_idx]:.6e}" \
                      if not np.isnan(pvalues[orig_row_idx]) else "NA"
                paj = f"{padj[orig_row_idx]:.6e}" \
                      if not np.isnan(padj[orig_row_idx]) else "NA"

                out.write(
                    line.rstrip("\n")
                    + f"\t{ivt_freq}\t{ivt_cov}\t{ivt_nmod}\t{pv}\t{paj}\n"
                )
                orig_row_idx += 1
                written += 1

    finally:
        if args.output != "-":
            out.close()

    print(f"[info] Done. Wrote {written} rows.", file=sys.stderr)


if __name__ == "__main__":
    main()
