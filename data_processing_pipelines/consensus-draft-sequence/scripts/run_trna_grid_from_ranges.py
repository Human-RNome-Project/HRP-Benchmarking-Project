#!/usr/bin/env python3
"""
tRNA grid search using quantile grids from tRNA_grid_search_parameter_ranges.tsv.

Runs 3-way exhaustive search for mods present in all three platforms,
and pairwise searches for mods present in exactly two platforms.
MassSpec has no coverage, so its parameter lists are score+level only.

Inputs:
  - Illumina_combined_polyARNA_tRNA_rRNA.bed
  - ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed
  - MS_rRNA_tRNA_harmonized.bed

Outputs:
  - tRNA_3way_best.tsv
  - tRNA_grid_search_all_results.tsv
  - tRNA_grid_search_best.tsv
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

TRNA_FILES = {
    "Illumina": INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed",
    "ONT": INPUTS / "ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed",
    "MassSpec": OUTPUTS / "MS_rRNA_tRNA_harmonized.bed",
}

RANGES_TSV = THRESHOLDS / "tRNA_grid_search_parameter_ranges.tsv"

def load_and_filter(files):
    data = {}
    for plat, path in files.items():
        df = parse_bed(path, plat)
        df = df[df["chrom"].astype(str).str.startswith("hs_tRNA") | df["chrom"].astype(str).str.startswith("hs_mttRNA")].copy()
        data[plat] = df
    return data

print("Loading tRNA BEDs...")
trna_data = load_and_filter(TRNA_FILES)
for plat, df in trna_data.items():
    print(f"  [{plat}] {len(df):,} rows, mods: {sorted(df['name'].unique())}")

print(f"\nLoading parameter ranges from {RANGES_TSV}...")
ranges_df = pd.read_csv(RANGES_TSV, sep="\t")

grid_lookup = {}
for _, row in ranges_df.iterrows():
    key = (row["platform"], row["mod_type"], row["parameter"])
    thresholds = [float(x.strip()) for x in str(row["grid_thresholds"]).split(",")]
    grid_lookup[key] = thresholds

def classify_positions(df, mod_type, score_t, cov_t, level_t):
    sub = df[df["name"] == mod_type].copy()
    mask = (sub["score"] >= score_t) & (sub["coverage"] >= cov_t) & (sub["level"] >= level_t)
    return set(zip(sub.loc[mask, "chrom"], sub.loc[mask, "start"]))

def get_grid(plat, mod, param):
    return grid_lookup.get((plat, mod, param), [])

def make_params(plat, mod):
    score = get_grid(plat, mod, "score")
    level = get_grid(plat, mod, "level")
    if plat == "MassSpec":
        return [score, level]
    cov = get_grid(plat, mod, "coverage")
    return [score, cov, level]

def unpack_params(plat, vals):
    if plat == "MassSpec":
        s, l = vals
        return s, 0.0, l
    return vals

def grid_search_pairwise(df1, plat1, df2, plat2, mod_type):
    params1 = make_params(plat1, mod_type)
    params2 = make_params(plat2, mod_type)

    if not all(params1) or not all(params2):
        return pd.DataFrame(), None

    sets1 = []
    for vals1 in product(*params1):
        s, c, l = unpack_params(plat1, vals1)
        sets1.append((classify_positions(df1, mod_type, s, c, l), s, c, l))

    sets2 = []
    for vals2 in product(*params2):
        s, c, l = unpack_params(plat2, vals2)
        sets2.append((classify_positions(df2, mod_type, s, c, l), s, c, l))

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

def grid_search_3way(mod_type):
    params_ill = make_params("Illumina", mod_type)
    params_ont = make_params("ONT", mod_type)
    params_ms = make_params("MassSpec", mod_type)

    if not all(params_ill) or not all(params_ont) or not all(params_ms):
        return None, 0

    sets_ill = []
    for vals in product(*params_ill):
        s, c, l = unpack_params("Illumina", vals)
        sets_ill.append((classify_positions(trna_data["Illumina"], mod_type, s, c, l), s, c, l))

    sets_ont = []
    for vals in product(*params_ont):
        s, c, l = unpack_params("ONT", vals)
        sets_ont.append((classify_positions(trna_data["ONT"], mod_type, s, c, l), s, c, l))

    sets_ms = []
    for vals in product(*params_ms):
        s, c, l = unpack_params("MassSpec", vals)
        sets_ms.append((classify_positions(trna_data["MassSpec"], mod_type, s, c, l), s, c, l))

    total = len(sets_ill) * len(sets_ont) * len(sets_ms)
    best_jaccard = -1
    best_result = None

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

    return best_result, total

plat_mods = defaultdict(set)
for _, row in ranges_df.iterrows():
    if row["biotype"] == "tRNA":
        plat_mods[row["mod_type"]].add(row["platform"])

mods_3way = sorted(m for m, plats in plat_mods.items() if len(plats) == 3)
mods_2way = sorted(m for m, plats in plat_mods.items() if len(plats) == 2)
print(f"\n3-way mods (all platforms): {mods_3way}")
print(f"2-way mods: {mods_2way}")

all_results = []
best_results = []

print("\n" + "=" * 70)
print("tRNA 3-WAY GRID SEARCH (from parameter ranges TSV)")
print("=" * 70)

for mod in mods_3way:
    result, total = grid_search_3way(mod)
    if result is None:
        print(f"\n{mod}: skipped (missing grids)")
        continue
    print(f"\n{mod}: {total:,} combos, best J={result['best_jaccard']:.4f}")
    best_results.append(result)

if best_results:
    df_3way = pd.DataFrame(best_results)
    df_3way.to_csv(OUTDIR / "tRNA_3way_best.tsv", sep="\t", index=False)
    print(f"\nSaved tRNA_3way_best.tsv ({len(df_3way)} rows)")
else:
    print("\nNo 3-way results generated.")

print("\n" + "=" * 70)
print("tRNA PAIRWISE GRID SEARCH (from parameter ranges TSV)")
print("=" * 70)

best_results = []  # reset for pairwise results only
for mod_type in mods_2way:
    platforms = sorted(plat_mods[mod_type])
    if len(platforms) < 2:
        continue

    print(f"\n--- {mod_type} (platforms: {', '.join(platforms)}) ---")

    for i, p1 in enumerate(platforms):
        for p2 in platforms[i + 1:]:
            if mod_type not in trna_data[p1]["name"].unique() or mod_type not in trna_data[p2]["name"].unique():
                print(f"  {p1}-{p2}: skipped (no data in BED)")
                continue

            res_df, best = grid_search_pairwise(trna_data[p1], p1, trna_data[p2], p2, mod_type)
            if res_df.empty:
                print(f"  {p1}-{p2}: no results")
                continue

            print(f"  {p1}-{p2}: {len(res_df):,} combos, best J={best['jaccard']:.4f}")
            all_results.append(res_df)
            best_results.append(best)

if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUTDIR / "tRNA_grid_search_all_results.tsv", sep="\t", index=False, float_format="%.6f")
    print(f"\nSaved tRNA_grid_search_all_results.tsv ({len(combined):,} rows)")

    best_df = pd.DataFrame(best_results)
    best_df = best_df.loc[best_df.groupby(["mod_type", "comparison"])["jaccard"].idxmax()].reset_index(drop=True)
    best_df.to_csv(OUTDIR / "tRNA_grid_search_best.tsv", sep="\t", index=False, float_format="%.6f")
    print(f"Saved tRNA_grid_search_best.tsv ({len(best_df)} rows)")
    for _, row in best_df.iterrows():
        print(f"  {row['mod_type']:8s} | {row['comparison']:25s} | J={row['jaccard']:.4f}")
else:
    print("\nNo pairwise results generated.")

print("\nDone!")
