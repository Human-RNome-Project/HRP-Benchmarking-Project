#!/usr/bin/env python3
"""
rRNA 3-way grid search for m5C using quantile grids from parameter ranges TSV.

Reads grids from rRNA_grid_search_parameter_ranges.tsv.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import product

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import parse_bed
from paths import ROOT, INPUTS, OUTPUTS, THRESHOLDS, FILTERED, STATE, DRAFT, TIERED, TIERED_TRNA

OUTDIR = THRESHOLDS

RRNA_FILES = {
    "Illumina": INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed",
    "ONT": INPUTS / "ONT_polyARNA_rRNA_combined.filtered.bed",
    "MassSpec": OUTPUTS / "MS_rRNA_tRNA_harmonized.bed",
}

RANGES_TSV = THRESHOLDS / "rRNA_grid_search_parameter_ranges.tsv"
TARGET_MODS = ["m6A"]

# ── Load BED data ───────────────────────────────────────────────────────────
def load_and_filter(files):
    data = {}
    for plat, path in files.items():
        df = parse_bed(path, plat)
        df = df[df["chrom"].astype(str).str.startswith("hs_rRNA")].copy()
        data[plat] = df
    return data

print("Loading rRNA BEDs...")
rrna_data = load_and_filter(RRNA_FILES)
for plat, df in rrna_data.items():
    print(f"  [{plat}] {len(df):,} rows, mods: {sorted(df['name'].unique())}")

# ── Load parameter ranges ───────────────────────────────────────────────────
print(f"\nLoading parameter ranges from {RANGES_TSV}...")
ranges_df = pd.read_csv(RANGES_TSV, sep="\t")

grid_lookup = {}
for _, row in ranges_df.iterrows():
    key = (row["platform"], row["mod_type"], row["parameter"])
    thresholds = [float(x.strip()) for x in str(row["grid_thresholds"]).split(",")]
    grid_lookup[key] = thresholds

# ── Helper functions ────────────────────────────────────────────────────────
def classify_positions(df, mod_type, score_t, cov_t, level_t):
    sub = df[df["name"] == mod_type].copy()
    mask = (sub["score"] >= score_t) & (sub["coverage"] >= cov_t) & (sub["level"] >= level_t)
    return set(zip(sub.loc[mask, "chrom"], sub.loc[mask, "start"]))

def get_grid(plat, mod, param):
    return grid_lookup.get((plat, mod, param), [])

def grid_search_3way_from_ranges(mod_type):
    g_ill = (get_grid("Illumina", mod_type, "score"),
             get_grid("Illumina", mod_type, "coverage"),
             get_grid("Illumina", mod_type, "level"))
    g_ont = (get_grid("ONT", mod_type, "score"),
             get_grid("ONT", mod_type, "coverage"),
             get_grid("ONT", mod_type, "level"))
    g_ms = (get_grid("MassSpec", mod_type, "score"),
            get_grid("MassSpec", mod_type, "level"))  # no coverage for MassSpec

    if not all(g_ill) or not all(g_ont) or not all(g_ms):
        return None, 0

    # Pre-compute sets
    sets_ill = []
    for s, c, l in product(*g_ill):
        sets_ill.append((classify_positions(rrna_data["Illumina"], mod_type, s, c, l), s, c, l))

    sets_ont = []
    for s, c, l in product(*g_ont):
        sets_ont.append((classify_positions(rrna_data["ONT"], mod_type, s, c, l), s, c, l))

    sets_ms = []
    for s, l in product(*g_ms):
        sets_ms.append((classify_positions(rrna_data["MassSpec"], mod_type, s, 0.0, l), s, 0.0, l))

    total = len(sets_ill) * len(sets_ont) * len(sets_ms)
    best_jaccard = -1
    best_result = None
    count = 0

    for set1, s1, c1, l1 in sets_ill:
        for set2, s2, c2, l2 in sets_ont:
            for set3, s3, c3, l3 in sets_ms:
                inter = set1 & set2 & set3
                union = set1 | set2 | set3
                jaccard = len(inter) / len(union) if union else 0.0

                if jaccard > best_jaccard:
                    best_jaccard = jaccard
                    best_result = {
                        "mod_type": mod_type,
                        "best_jaccard": jaccard,
                        "n_intersection": len(inter),
                        "n_union": len(union),
                        "n_Illumina": len(set1),
                        "n_ONT": len(set2),
                        "n_MassSpec": len(set3),
                        "Illumina_score": s1,
                        "Illumina_coverage": c1,
                        "Illumina_level": l1,
                        "ONT_score": s2,
                        "ONT_coverage": c2,
                        "ONT_level": l2,
                        "MassSpec_score": s3,
                        "MassSpec_coverage": c3,
                        "MassSpec_level": l3,
                    }
                count += 1

    return best_result, total

# ── Run 3-way searches ──────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("3-WAY GRID SEARCH for m5C (from parameter ranges TSV)")
print("=" * 70)

best_results = []
for mod in TARGET_MODS:
    result, total = grid_search_3way_from_ranges(mod)
    if result is None:
        print(f"\n{mod}: skipped (missing grids)")
        continue
    print(f"\n{mod}: {total:,} combos, best J={result['best_jaccard']:.4f}")
    best_results.append(result)

if best_results:
    df = pd.DataFrame(best_results)
    df.to_csv(OUTDIR / "rRNA_3way_m6A_best.tsv", sep="\t", index=False)
    print(f"\nSaved rRNA_3way_m6A_best.tsv ({len(df)} rows)")
else:
    print("\nNo results generated.")

print("\nDone!")
