#!/usr/bin/env python3
"""
polyA pairwise grid search using numpy boolean arrays for speed.

Reads grids from polyA_grid_search_parameter_ranges.tsv.
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

POLYA_FILES = {
    "Illumina": INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed",
    "ONT": INPUTS / "ONT_polyARNA_rRNA_combined.filtered.bed",
}

RANGES_TSV = THRESHOLDS / "polyA_grid_search_parameter_ranges.tsv"

# ── Load BED data ───────────────────────────────────────────────────────────
def load_and_filter(files):
    data = {}
    for plat, path in files.items():
        df = parse_bed(path, plat)
        df = df[~(df["chrom"].astype(str).str.startswith("hs_rRNA") |
                  df["chrom"].astype(str).str.startswith("hs_tRNA"))].copy()
        data[plat] = df
    return data

print("Loading polyA BEDs...")
polya_data = load_and_filter(POLYA_FILES)
for plat, df in polya_data.items():
    print(f"  [{plat}] {len(df):,} rows, mods: {sorted(df['name'].unique())}")

# ── Load parameter ranges ───────────────────────────────────────────────────
print(f"\nLoading parameter ranges from {RANGES_TSV}...")
ranges_df = pd.read_csv(RANGES_TSV, sep="\t")

grid_lookup = {}
for _, row in ranges_df.iterrows():
    key = (row["platform"], row["mod_type"], row["parameter"])
    thresholds = [float(x.strip()) for x in str(row["grid_thresholds"]).split(",")]
    grid_lookup[key] = thresholds

def get_grid(plat, mod, param):
    return grid_lookup.get((plat, mod, param), [])

# ── Fast grid search with boolean arrays ────────────────────────────────────
def grid_search_pairwise_fast(df1, plat1, df2, plat2, mod_type):
    g1_score = get_grid(plat1, mod_type, "score")
    g1_cov = get_grid(plat1, mod_type, "coverage")
    g1_level = get_grid(plat1, mod_type, "level")
    g2_score = get_grid(plat2, mod_type, "score")
    g2_cov = get_grid(plat2, mod_type, "coverage")
    g2_level = get_grid(plat2, mod_type, "level")

    if not all([g1_score, g1_cov, g1_level, g2_score, g2_cov, g2_level]):
        return pd.DataFrame(), None

    # Build position index mapping
    sub1 = df1[df1["name"] == mod_type]
    sub2 = df2[df2["name"] == mod_type]
    all_pos = sorted(set(zip(sub1["chrom"].values, sub1["start"].values)) |
                     set(zip(sub2["chrom"].values, sub2["start"].values)))
    pos_to_idx = {p: i for i, p in enumerate(all_pos)}
    N = len(all_pos)

    # Pre-compute boolean masks for platform 1
    score1 = sub1["score"].values
    cov1 = sub1["coverage"].values
    level1 = sub1["level"].values
    idx1 = np.array([pos_to_idx[(c, s)] for c, s in zip(sub1["chrom"].values, sub1["start"].values)])

    m1s = {s: np.zeros(N, dtype=bool) for s in g1_score}
    m1c = {c: np.zeros(N, dtype=bool) for c in g1_cov}
    m1l = {l: np.zeros(N, dtype=bool) for l in g1_level}
    for i in range(len(sub1)):
        for s in g1_score:
            if score1[i] >= s: m1s[s][idx1[i]] = True
        for c in g1_cov:
            if cov1[i] >= c: m1c[c][idx1[i]] = True
        for l in g1_level:
            if level1[i] >= l: m1l[l][idx1[i]] = True

    # Pre-compute boolean masks for platform 2
    score2 = sub2["score"].values
    cov2 = sub2["coverage"].values
    level2 = sub2["level"].values
    idx2 = np.array([pos_to_idx[(c, s)] for c, s in zip(sub2["chrom"].values, sub2["start"].values)])

    m2s = {s: np.zeros(N, dtype=bool) for s in g2_score}
    m2c = {c: np.zeros(N, dtype=bool) for c in g2_cov}
    m2l = {l: np.zeros(N, dtype=bool) for l in g2_level}
    for i in range(len(sub2)):
        for s in g2_score:
            if score2[i] >= s: m2s[s][idx2[i]] = True
        for c in g2_cov:
            if cov2[i] >= c: m2c[c][idx2[i]] = True
        for l in g2_level:
            if level2[i] >= l: m2l[l][idx2[i]] = True

    # Grid search with vectorized ops
    results = []
    best_jaccard = -1.0
    best_result = None

    for s1, c1, l1 in product(g1_score, g1_cov, g1_level):
        mask1 = m1s[s1] & m1c[c1] & m1l[l1]
        n1 = int(mask1.sum())
        for s2, c2, l2 in product(g2_score, g2_cov, g2_level):
            mask2 = m2s[s2] & m2c[c2] & m2l[l2]
            n2 = int(mask2.sum())
            n_inter = int((mask1 & mask2).sum())
            n_union = int((mask1 | mask2).sum())
            jaccard = n_inter / n_union if n_union > 0 else 0.0

            result = {
                "mod_type": mod_type,
                "comparison": f"{plat1}-{plat2}",
                "jaccard": jaccard,
                "n_intersection": n_inter,
                "n_union": n_union,
                "n_plat1": n1,
                "n_plat2": n2,
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

# ── Determine mods to search ────────────────────────────────────────────────
plat_mods = {}
for _, row in ranges_df.iterrows():
    if row["biotype"] == "polyA":
        mod = row["mod_type"]
        if mod not in plat_mods:
            plat_mods[mod] = set()
        plat_mods[mod].add(row["platform"])

target_mods = sorted(m for m, plats in plat_mods.items() if len(plats) >= 2)
print(f"\nTarget mods: {target_mods}")

# ── Run searches ────────────────────────────────────────────────────────────
all_results = []
best_results = []

print("\n" + "=" * 70)
print("polyA PAIRWISE GRID SEARCH (numpy boolean arrays)")
print("=" * 70)

for mod_type in target_mods:
    platforms = sorted(plat_mods[mod_type])
    if len(platforms) < 2:
        continue

    print(f"\n--- {mod_type} ---")
    for i, p1 in enumerate(platforms):
        for p2 in platforms[i + 1:]:
            if mod_type not in polya_data[p1]["name"].unique() or mod_type not in polya_data[p2]["name"].unique():
                print(f"  {p1}-{p2}: skipped")
                continue

            res_df, best = grid_search_pairwise_fast(polya_data[p1], p1, polya_data[p2], p2, mod_type)
            if res_df.empty:
                print(f"  {p1}-{p2}: no results")
                continue

            print(f"  {p1}-{p2}: {len(res_df):,} combos, best J={best['jaccard']:.6f}")
            all_results.append(res_df)
            best_results.append(best)

# ── Save results ────────────────────────────────────────────────────────────
if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(OUTDIR / "polyA_grid_search_all_results.tsv", sep="\t", index=False, float_format="%.6f")
    print(f"\nSaved polyA_grid_search_all_results.tsv ({len(combined):,} rows)")

    best_df = pd.DataFrame(best_results)
    best_df = best_df.loc[best_df.groupby(["mod_type", "comparison"])["jaccard"].idxmax()].reset_index(drop=True)
    best_df.to_csv(OUTDIR / "polyA_grid_search_best.tsv", sep="\t", index=False, float_format="%.6f")
    print(f"Saved polyA_grid_search_best.tsv ({len(best_df)} rows)")
    for _, row in best_df.iterrows():
        print(f"  {row['mod_type']:8s} | {row['comparison']:25s} | J={row['jaccard']:.6f}")
else:
    print("\nNo results generated.")

print("\nDone!")
