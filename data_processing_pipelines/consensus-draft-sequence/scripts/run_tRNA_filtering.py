#!/usr/bin/env python3
"""
Apply optimal tRNA grid-search thresholds to platform BED files
to create filtered subset BEDs for tRNA.

Uses tRNA_grid_search_best.tsv for 2-way mods and tRNA_3way_best.tsv for 3-way mods.
"""

import pandas as pd
from pathlib import Path
import numpy as np
from paths import ROOT, INPUTS, OUTPUTS, THRESHOLDS, FILTERED, STATE, DRAFT, TIERED, TIERED_TRNA

OUTDIR = FILTERED
OUTDIR.mkdir(exist_ok=True, parents=True)

TRNA_BEST = THRESHOLDS / "tRNA_grid_search_best.tsv"
TRNA_3WAY_BEST = THRESHOLDS / "tRNA_3way_best.tsv"

ILLUMINA_BED = INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed"
ONT_BED = INPUTS / "ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed"
MS_BED = OUTPUTS / "MS_rRNA_tRNA_harmonized.bed"

def load_best_thresholds(tsv_path):
    """Load grid-search best results; select optimal thresholds per mod type."""
    df = pd.read_csv(tsv_path, sep="\t")
    platform_labels = ["Illumina", "ONT", "MassSpec"]
    thresholds = {}

    for mod in df["mod_type"].unique():
        sub = df[df["mod_type"] == mod].copy()
        jaccard_col = "best_jaccard" if "best_jaccard" in sub.columns else "jaccard"
        best_row = sub.loc[sub[jaccard_col].idxmax()]
        mod_thresh = {"jaccard": best_row[jaccard_col], "comparison": best_row["comparison"]}

        for plat in platform_labels:
            score_col = f"{plat}_score"
            cov_col = f"{plat}_coverage"
            lvl_col = f"{plat}_level"
            if score_col in best_row and pd.notna(best_row[score_col]):
                mod_thresh[plat] = {
                    "score": float(best_row[score_col]),
                    "coverage": float(best_row[cov_col]),
                    "level": float(best_row[lvl_col]),
                }
        thresholds[mod] = mod_thresh

    return thresholds

def load_3way_thresholds(tsv_path):
    """Load three-way grid-search best results. Each row is already the best combo."""
    df = pd.read_csv(tsv_path, sep="\t")
    platform_labels = ["Illumina", "ONT", "MassSpec"]
    thresholds = {}

    for _, row in df.iterrows():
        mod = row["mod_type"]
        mod_thresh = {"jaccard": row["best_jaccard"], "comparison": "3-way"}

        for plat in platform_labels:
            score_col = f"{plat}_score"
            cov_col = f"{plat}_coverage"
            lvl_col = f"{plat}_level"
            if score_col in row and pd.notna(row[score_col]):
                mod_thresh[plat] = {
                    "score": float(row[score_col]),
                    "coverage": float(row[cov_col]),
                    "level": float(row[lvl_col]),
                }
        thresholds[mod] = mod_thresh

    return thresholds

def filter_bed(in_path, out_path, thresholds_dict, platform_key):
    """Read an original BED line-by-line, apply tRNA thresholds, write subset."""
    n_in = 0
    n_out = 0
    mod_counts = {}

    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue

            chrom = p[0]
            if not (chrom.startswith("hs_tRNA") or chrom.startswith("hs_mttRNA")):
                continue

            n_in += 1

            raw_name = p[3]
            score = float(p[4]) if p[4] not in (".", "NA", "") else 0.0
            coverage = float(p[9]) if len(p) > 9 and p[9] not in (".", "NA", "") else 0.0
            frequency = float(p[10]) if len(p) > 10 and p[10] not in (".", "NA", "") else 0.0

            if raw_name in thresholds_dict and platform_key in thresholds_dict[raw_name]:
                t = thresholds_dict[raw_name][platform_key]
                if not (score >= t["score"] and coverage >= t["coverage"] and frequency >= t["level"]):
                    continue

            fout.write(line)
            n_out += 1
            mod_counts[raw_name] = mod_counts.get(raw_name, 0) + 1

    return n_in, n_out, mod_counts

def main():
    print("=== Loading tRNA thresholds ===")
    trna_thresholds = load_best_thresholds(TRNA_BEST)

    if TRNA_3WAY_BEST.exists():
        three_way = load_3way_thresholds(TRNA_3WAY_BEST)
        for mod in three_way:
            trna_thresholds[mod] = three_way[mod]
            print(f"  [3-way] {mod}: Jaccard={three_way[mod]['jaccard']:.4f}")

    print(f"  tRNA: {len(trna_thresholds)} mod types")

    files = [
        (ILLUMINA_BED, "Illumina", "Illumina"),
        (ONT_BED, "ONT", "ONT"),
        (MS_BED, "MassSpec", "MassSpec"),
    ]

    summary = []

    for bed_path, plat_key, plat_label in files:
        if not bed_path.exists():
            print(f"\n[SKIP] {plat_label}: file not found: {bed_path}")
            continue

        plat_fname = plat_key.replace(" ", "_").lower()
        out_path = OUTDIR / f"tRNA_{plat_fname}_filtered.bed"

        n_in, n_out, mod_counts = filter_bed(
            bed_path, out_path, trna_thresholds, plat_key
        )

        print(f"\n{plat_label} / tRNA: {n_in:,} → {n_out:,} rows → {out_path.name}")
        for mod, cnt in sorted(mod_counts.items()):
            print(f"  {mod}: {cnt}")

        summary.append({
            "biotype": "tRNA",
            "platform": plat_label,
            "n_in": n_in,
            "n_out": n_out,
            "pct": round(100 * n_out / n_in, 2) if n_in > 0 else 0.0,
        })

    # Save summary
    summary_path = OUTDIR / "tRNA_filtering_summary.tsv"
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\nSummary saved to: {summary_path}")

    print("\n=== tRNA Filtering Summary ===")
    for _, row in summary_df.iterrows():
        print(f"  {row['biotype']}/{row['platform']}: {row['n_in']:,} → {row['n_out']:,} ({row['pct']:.1f}%)")

if __name__ == "__main__":
    main()
