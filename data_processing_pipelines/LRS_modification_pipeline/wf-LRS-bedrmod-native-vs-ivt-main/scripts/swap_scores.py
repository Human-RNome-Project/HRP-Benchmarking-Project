#!/usr/bin/env python3
"""
Produce a "scores-swapped" copy of a filtered BEDRMod for visualisation.

Two columns are rewritten; every row, coordinate and other column is preserved:

  score  ->  statistical significance  -log10(padj)
             * padj is NA (non-testable)  -> -log10(1) = 0
             * padj == 0 (most significant) -> a fixed constant (--zero-score),
               chosen above the ~323 ceiling of -log10(padj) so the best sites
               always get the top code, consistently across biotypes
  color  ->  original score (0-1) encoded as red intensity  "R,0,0"  (R = score*255)

The original score is also retained verbatim in a new trailing column,
`composite_score`, so no information is lost.

The original score therefore remains readable as redness, while the score
column now drives significance-based shading in a genome browser.

Column names are auto-detected from the #[...] header line, which is rewritten
with the appended composite_score column.
"""

import ast
import gzip
import sys
import argparse
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", help="Filtered BEDRMod file (plain or .gz)")
    p.add_argument("-o", "--output", required=True,
                   help="Output BEDRMod path (.gz to compress)")
    p.add_argument(
        "--zero-score", type=int, default=350,
        help="fixed score assigned to padj==0 sites (default: 350). Set above "
             "the ~323 maximum of -log10(padj) so the best sites always rank "
             "top, identically across biotypes.",
    )
    return p.parse_args()


def _open(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def read_bedrmod(path):
    """Return (DataFrame[str], comment_lines, col_names)."""
    comments, col_names = [], None
    with _open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            comments.append(line.rstrip("\n"))
            stripped = line.lstrip("#").strip()
            if stripped.startswith("["):
                try:
                    parsed = ast.literal_eval(stripped)
                except (ValueError, SyntaxError):
                    continue
                if isinstance(parsed, list) and "chrom" in parsed:
                    col_names = parsed
    if col_names is None:
        raise ValueError(f"No '#[...]' column header found in {path}")

    df = pd.read_csv(
        path, sep="\t", skiprows=len(comments), header=None,
        names=col_names, dtype=str,
    )
    return df, comments, col_names


def main():
    args = parse_args()
    df, comments, cols = read_bedrmod(args.input)
    for c in ("score", "padj", "color"):
        if c not in cols:
            raise ValueError(f"Column '{c}' missing from {args.input}")

    composite = df["score"].copy()   # keep the original score values verbatim
    score = pd.to_numeric(df["score"], errors="coerce").to_numpy(dtype=float)
    padj  = pd.to_numeric(df["padj"],  errors="coerce").to_numpy(dtype=float)

    # ── color: original score -> red intensity "R,0,0" ────────────────────
    red = np.clip(np.round(np.nan_to_num(score, nan=0.0) * 255), 0, 255).astype(int)
    df["color"] = [f"{r},0,0" for r in red]

    # ── score: -log10(padj), rounded to an integer (BED spec: int in 0-1000) ─
    neglog = np.zeros(len(padj))                   # NA (non-testable) -> 0
    na   = np.isnan(padj)
    zero = ~na & (padj <= 0)
    pos  = ~na & (padj > 0)

    # -log10(padj) tops out around 323 for the smallest double, so it always
    # fits the 0-1000 range; round to the nearest integer.
    neglog[pos]  = np.round(-np.log10(padj[pos]))
    neglog[zero] = args.zero_score                 # fixed top code, all biotypes

    max_finite = neglog[pos].max() if pos.any() else 0.0
    score = np.clip(np.round(neglog), 0, 1000).astype(int)
    df["score"] = score

    # ── retain original score in a new trailing column ────────────────────
    df["composite_score"] = composite
    out_cols = cols + ["composite_score"]

    print(
        f"[info] {len(df):,} sites | testable(pos): {int(pos.sum()):,} | "
        f"padj==0: {int(zero.sum()):,} (score={args.zero_score}) | "
        f"NA: {int(na.sum()):,} (score=0) | max finite -log10(padj)={int(max_finite)}",
        file=sys.stderr,
    )

    with _open(args.output, "wt") as out:
        for c in comments:
            stripped = c.lstrip("#").strip()
            if stripped.startswith("[") and "chrom" in stripped:
                out.write("#" + repr(out_cols) + "\n")   # header + composite_score
            else:
                out.write(c + "\n")
        df.to_csv(out, sep="\t", header=False, index=False)

    print(f"[info] Written -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
