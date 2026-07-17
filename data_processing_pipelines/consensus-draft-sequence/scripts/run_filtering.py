#!/usr/bin/env python3
"""
Apply optimal grid-search thresholds to platform BED files
to create filtered subset BEDs.

Preserves all 13 original columns exactly, except that for rRNA MassSpec
rows the legacy ambiguous symbols mxA/mxC/mxG/mxU in the last two columns
are normalized to mA?/mC?/mG?/mU? to match the original filtered BEDs.

Uses standardized names only for threshold matching;
original names (Y, Psi, etc.) are preserved in output.
"""

import pandas as pd
from pathlib import Path
from paths import ROOT, INPUTS, OUTPUTS, THRESHOLDS, FILTERED, STATE, DRAFT, TIERED, TIERED_TRNA

# Legacy ambiguous symbol normalization for rRNA MassSpec raw columns.
MX_TO_MQ = {"mxA": "mA?", "mxC": "mC?", "mxG": "mG?", "mxU": "mU?"}

OUTDIR = FILTERED
OUTDIR.mkdir(exist_ok=True, parents=True)

RRNA_BEST = THRESHOLDS / "rRNA_grid_search_best.tsv"
RRNA_3WAY_BEST = THRESHOLDS / "rRNA_3way_Nm_best.tsv"
POLYA_BEST = THRESHOLDS / "polyA_grid_search_best.tsv"

ILLUMINA_BED = INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed"
ONT_BED = INPUTS / "ONT_polyARNA_rRNA_combined.filtered.bed"
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

def filter_bed(in_path, out_path, thresholds_dict, platform_key, biotype):
    """
    Read an original BED file line-by-line, apply thresholds,
    and write an exact subset preserving all 13 columns and headers.
    """
    n_in = 0
    n_out = 0
    mod_counts = {}

    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            # Preserve comment/header lines unchanged
            if line.startswith("#"):
                fout.write(line)
                continue

            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue

            chrom = p[0]

            # Biotype filter by chromosome prefix
            if biotype == "rRNA":
                if not chrom.startswith("hs_rRNA"):
                    continue
            else:  # polyA
                if chrom.startswith("hs_rRNA") or chrom.startswith("hs_tRNA"):
                    continue

            n_in += 1

            # Extract fields
            raw_name = p[3]
            score = float(p[4]) if p[4] not in (".", "NA", "") else 0.0
            coverage = float(p[9]) if len(p) > 9 and p[9] not in (".", "NA", "") else 0.0
            frequency = float(p[10]) if len(p) > 10 and p[10] not in (".", "NA", "") else 0.0

            # Apply thresholds if they exist for this mod+platform
            if raw_name in thresholds_dict and platform_key in thresholds_dict[raw_name]:
                t = thresholds_dict[raw_name][platform_key]
                if not (score >= t["score"] and coverage >= t["coverage"] and frequency >= t["level"]):
                    continue

            # Row passes — write it back preserving columns.
            # For rRNA rows, normalize any remaining mx* symbols in the last
            # two raw-name columns (single_letter_code, mod_id) to legacy m? form.
            if biotype == "rRNA" and len(p) >= 13:
                p[11] = MX_TO_MQ.get(p[11], p[11])
                p[12] = MX_TO_MQ.get(p[12], p[12])

            fout.write("\t".join(p) + "\n")
            n_out += 1
            mod_counts[raw_name] = mod_counts.get(raw_name, 0) + 1

    return n_in, n_out, mod_counts

def main():
    print("=== Loading thresholds ===")
    rrna_thresholds = load_best_thresholds(RRNA_BEST)
    polya_thresholds = load_best_thresholds(POLYA_BEST)
    
    # Override Nm mods with 3-way thresholds
    if RRNA_3WAY_BEST.exists():
        three_way = load_3way_thresholds(RRNA_3WAY_BEST)
        for mod in ["Am", "Cm", "Gm", "Um"]:
            if mod in three_way:
                rrna_thresholds[mod] = three_way[mod]
                print(f"  [3-way] {mod}: Jaccard={three_way[mod]['jaccard']:.4f}")
    
    print(f"  rRNA: {len(rrna_thresholds)} mod types")
    print(f"  polyA: {len(polya_thresholds)} mod types")

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

        for biotype, thresholds in [("rRNA", rrna_thresholds), ("polyA", polya_thresholds)]:
            plat_fname = plat_key.replace(" ", "_").lower()
            out_path = OUTDIR / f"{biotype}_{plat_fname}_filtered.bed"

            n_in, n_out, mod_counts = filter_bed(
                bed_path, out_path, thresholds, plat_key, biotype
            )

            if n_in == 0 and biotype == "polyA" and plat_label == "MassSpec":
                # MS has no polyA data — expected, remove empty file
                out_path.unlink(missing_ok=True)
                print(f"\n{plat_label} / {biotype}: no data (expected for MassSpec)")
                continue

            print(f"\n{plat_label} / {biotype}: {n_in:,} → {n_out:,} rows → {out_path.name}")
            for mod, cnt in sorted(mod_counts.items()):
                print(f"  {mod}: {cnt}")

            summary.append({
                "biotype": biotype,
                "platform": plat_label,
                "n_in": n_in,
                "n_out": n_out,
                "pct": round(100 * n_out / n_in, 2) if n_in > 0 else 0.0,
            })

    # Save summary
    summary_df = pd.DataFrame(summary)
    summary_path = OUTDIR / "filtering_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\nSummary saved to: {summary_path}")

    print("\n=== Filtering Summary ===")
    for _, row in summary_df.iterrows():
        print(f"  {row['biotype']}/{row['platform']}: {row['n_in']:,} → {row['n_out']:,} ({row['pct']:.1f}%)")

if __name__ == "__main__":
    main()
