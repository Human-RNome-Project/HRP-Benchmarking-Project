#!/usr/bin/env python3
"""
rRNA 2-way pairwise grid search using quantile grids from parameter ranges TSV.

Reads grids from rRNA_grid_search_parameter_ranges.tsv instead of unique data triples.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import product
from collections import defaultdict

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

# Build grid lookup: (platform, mod_type, parameter) -> list of thresholds
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

def grid_search_pairwise_from_ranges(df1, plat1, df2, plat2, mod_type):
    g1_score = get_grid(plat1, mod_type, "score")
    g1_cov = get_grid(plat1, mod_type, "coverage")
    g1_level = get_grid(plat1, mod_type, "level")
    g2_score = get_grid(plat2, mod_type, "score")
    g2_cov = get_grid(plat2, mod_type, "coverage")
    g2_level = get_grid(plat2, mod_type, "level")

    # MassSpec has no coverage grid — build parameter lists accordingly
    params1 = [g1_score, g1_level] if plat1 == "MassSpec" else [g1_score, g1_cov, g1_level]
    params2 = [g2_score, g2_level] if plat2 == "MassSpec" else [g2_score, g2_cov, g2_level]

    if not all(params1) or not all(params2):
        return pd.DataFrame(), None

    # Pre-compute all classified sets for platform 1
    sets1 = []
    for vals1 in product(*params1):
        if plat1 == "MassSpec":
            s, l = vals1
            c = 0.0  # coverage is 0 for MassSpec, threshold doesn't matter
        else:
            s, c, l = vals1
        classified = classify_positions(df1, mod_type, s, c, l)
        sets1.append((classified, s, c, l))

    # Pre-compute all classified sets for platform 2
    sets2 = []
    for vals2 in product(*params2):
        if plat2 == "MassSpec":
            s, l = vals2
            c = 0.0
        else:
            s, c, l = vals2
        classified = classify_positions(df2, mod_type, s, c, l)
        sets2.append((classified, s, c, l))

    results = []
    best_jaccard = -1
    best_result = None

    for set1, s1, c1, l1 in sets1:
        for set2, s2, c2, l2 in sets2:
            intersection = set1 & set2
            union = set1 | set2
            jaccard = len(intersection) / len(union) if union else 0.0

            result = {
                "mod_type": mod_type,
                "comparison": f"{plat1}-{plat2}",
                "jaccard": jaccard,
                "n_intersection": len(intersection),
                "n_union": len(union),
                "n_plat1": len(set1),
                "n_plat2": len(set2),
                f"{plat1}_score": s1,
                f"{plat1}_coverage": c1,
                f"{plat1}_level": l1,
                f"{plat2}_score": s2,
                f"{plat2}_coverage": c2,
                f"{plat2}_level": l2,
            }
            results.append(result)
            if jaccard > best_jaccard:
                best_jaccard = jaccard
                best_result = result.copy()

    return pd.DataFrame(results), best_result

# ── Determine which mods to search ──────────────────────────────────────────
# Use mods from ranges TSV that appear in at least 2 platforms
plat_mods = defaultdict(set)
for _, row in ranges_df.iterrows():
    if row["biotype"] == "rRNA":
        plat_mods[row["mod_type"]].add(row["platform"])

target_mods = sorted(m for m, plats in plat_mods.items() if len(plats) >= 2)
print(f"\nTarget mods (≥2 platforms in ranges): {target_mods}")

# Also need platform pairs per mod
all_results = []
best_results = []

print("\n" + "=" * 70)
print("PAIRWISE GRID SEARCH (from parameter ranges TSV)")
print("=" * 70)

for mod_type in target_mods:
    platforms = sorted(plat_mods[mod_type])
    if len(platforms) < 2:
        continue

    print(f"\n--- {mod_type} (platforms: {', '.join(platforms)}) ---")

    for i, p1 in enumerate(platforms):
        for p2 in platforms[i + 1:]:
            # Check both platforms actually have data for this mod
            if mod_type not in rrna_data[p1]["name"].unique() or mod_type not in rrna_data[p2]["name"].unique():
                print(f"  {p1}-{p2}: skipped (no data in BED)")
                continue

            res_df, best = grid_search_pairwise_from_ranges(
                rrna_data[p1], p1, rrna_data[p2], p2, mod_type
            )
            if res_df.empty:
                print(f"  {p1}-{p2}: no results")
                continue

            print(f"  {p1}-{p2}: {len(res_df):,} combos, best J={best['jaccard']:.4f}")
            all_results.append(res_df)
            best_results.append(best)

# ── Save results ────────────────────────────────────────────────────────────
if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUTDIR / "rRNA_grid_search_all_results.tsv", sep="\t", index=False, float_format="%.6f")
    print(f"\nSaved rRNA_grid_search_all_results.tsv ({len(combined):,} rows)")

    best_df = pd.DataFrame(best_results)
    best_df = best_df.loc[best_df.groupby(["mod_type", "comparison"])["jaccard"].idxmax()].reset_index(drop=True)
    best_df.to_csv(OUTDIR / "rRNA_grid_search_best.tsv", sep="\t", index=False, float_format="%.6f")
    print(f"Saved rRNA_grid_search_best.tsv ({len(best_df)} rows)")
    for _, row in best_df.iterrows():
        print(f"  {row['mod_type']:8s} | {row['comparison']:25s} | J={row['jaccard']:.4f}")
else:
    print("\nNo results generated.")

print("\nDone!")
