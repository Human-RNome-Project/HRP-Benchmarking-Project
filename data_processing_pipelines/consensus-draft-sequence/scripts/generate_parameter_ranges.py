#!/usr/bin/env python3
"""
Draft Reference Pipeline — Step 1: Grid Search Parameter Ranges

Generates quantile-based parameter ranges (score, coverage, level) per platform
and mod type for both rRNA and polyA biotypes.

Input: Local BED files in draft_reference/
Output: TSVs saved to outputs/thresholds/
"""

import pandas as pd
import numpy as np
from pathlib import Path

from utils import parse_bed, make_grid
from paths import ROOT, INPUTS, OUTPUTS, THRESHOLDS, FILTERED, STATE, DRAFT, TIERED, TIERED_TRNA

# ── Paths ────────────────────────────────────────────────────────────────────

OUTDIR = THRESHOLDS
OUTDIR.mkdir(exist_ok=True)

ILLUMINA_BED = INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed"
ONT_BED = INPUTS / "ONT_polyARNA_rRNA_combined.filtered.bed"
MS_BED = OUTPUTS / "MS_rRNA_tRNA_harmonized.bed"

for label, path in [("Illumina", ILLUMINA_BED), ("ONT", ONT_BED), ("MassSpec", MS_BED)]:
    exists = "✓" if path.exists() else "✗ MISSING"
    print(f"{exists} {label}: {path}")

# ── Load and split by biotype ───────────────────────────────────────────────
def load_and_split(bed_path, platform_label):
    """Load BED and split into rRNA and polyA DataFrames by chromosome."""
    df = parse_bed(bed_path, platform_label)
    chrom = df["chrom"].astype(str)
    rrna = df[chrom.str.startswith("hs_rRNA")].copy()
    polya = df[~(chrom.str.startswith("hs_rRNA") | chrom.str.startswith("hs_tRNA"))].copy()
    return rrna, polya

print("\n=== Illumina ===")
ill_rrna, ill_polya = load_and_split(ILLUMINA_BED, "Illumina")
print(f"  rRNA: {len(ill_rrna):,} rows, mods: {sorted(ill_rrna['name'].unique())}")
print(f"  polyA: {len(ill_polya):,} rows, mods: {sorted(ill_polya['name'].unique())}")

print("\n=== ONT ===")
ont_rrna, ont_polya = load_and_split(ONT_BED, "ONT")
print(f"  rRNA: {len(ont_rrna):,} rows, mods: {sorted(ont_rrna['name'].unique())}")
print(f"  polyA: {len(ont_polya):,} rows, mods: {sorted(ont_polya['name'].unique())}")

print("\n=== MassSpec ===")
ms_rrna, ms_polya = load_and_split(MS_BED, "MassSpec")
print(f"  rRNA: {len(ms_rrna):,} rows, mods: {sorted(ms_rrna['name'].unique())}")
print(f"  polyA: {len(ms_polya):,} rows (expected 0)")

rrna_data = {
    "Illumina": ill_rrna,
    "ONT": ont_rrna,
    "MassSpec": ms_rrna,
}

polya_data = {
    "Illumina": ill_polya,
    "ONT": ont_polya,
    "MassSpec": ms_polya,
}

# ── Compute parameter ranges ────────────────────────────────────────────────
def compute_param_ranges(df, platform, biotype):
    """Compute parameter ranges per mod type for a given platform DataFrame.
    
    Skip coverage for MassSpec — it is uniformly 0.0 and has no discriminative power.
    """
    records = []
    for mod in sorted(df["name"].unique()):
        sub = df[df["name"] == mod]
        params = ["score", "coverage", "level"]
        if platform == "MassSpec":
            params = ["score", "level"]  # exclude coverage
        for param in params:
            vals = sub[param].dropna()
            if len(vals) == 0:
                continue
            grid = make_grid(vals, n_points=7)
            records.append({
                "biotype": biotype,
                "platform": platform,
                "mod_type": mod,
                "parameter": param,
                "n_values": len(vals),
                "min": round(vals.min(), 4),
                "max": round(vals.max(), 4),
                "mean": round(vals.mean(), 4),
                "median": round(vals.median(), 4),
                "n_grid_points": len(grid),
                "grid_thresholds": ", ".join(str(round(g, 4)) for g in grid),
            })
    return pd.DataFrame(records)

# rRNA
print("\n=== Computing rRNA parameter ranges ===")
rrna_ranges = pd.concat([
    compute_param_ranges(d, p, "rRNA") for p, d in rrna_data.items() if not d.empty
], ignore_index=True)
rrna_ranges.sort_values(["mod_type", "platform", "parameter"]).to_csv(
    OUTDIR / "rRNA_grid_search_parameter_ranges.tsv", sep="\t", index=False
)
print(f"Saved rRNA_grid_search_parameter_ranges.tsv ({len(rrna_ranges)} rows)")

# polyA
print("\n=== Computing polyA parameter ranges ===")
polya_ranges = pd.concat([
    compute_param_ranges(d, p, "polyA") for p, d in polya_data.items() if not d.empty
], ignore_index=True)
polya_ranges.sort_values(["mod_type", "platform", "parameter"]).to_csv(
    OUTDIR / "polyA_grid_search_parameter_ranges.tsv", sep="\t", index=False
)
print(f"Saved polyA_grid_search_parameter_ranges.tsv ({len(polya_ranges)} rows)")

# ── Search space size summary ───────────────────────────────────────────────
def compute_search_space_2way(ranges_df, biotype):
    results = []
    for mod in sorted(ranges_df[ranges_df["biotype"] == biotype]["mod_type"].unique()):
        sub = ranges_df[(ranges_df["mod_type"] == mod) & (ranges_df["biotype"] == biotype)]
        plat_params = {}
        for _, row in sub.iterrows():
            plat = row["platform"]
            param = row["parameter"]
            n_grid = row["n_grid_points"]
            if plat not in plat_params:
                plat_params[plat] = {}
            plat_params[plat][param] = n_grid
        platforms = list(plat_params.keys())
        for i in range(len(platforms)):
            for j in range(i + 1, len(platforms)):
                p1, p2 = platforms[i], platforms[j]
                n1 = int(np.prod(list(plat_params[p1].values())))
                n2 = int(np.prod(list(plat_params[p2].values())))
                results.append({
                    "biotype": biotype,
                    "mod_type": mod,
                    "comparison": f"{p1} × {p2}",
                    "p1_combos": n1,
                    "p2_combos": n2,
                    "total_combos": n1 * n2,
                    "search_type": "2-way",
                })
    return pd.DataFrame(results)

def compute_search_space_3way_from_ranges(ranges_df, mod_types=["Am", "Cm", "Gm", "Um"]):
    """Compute 3-way search space for Nm mods using quantile grid sizes from ranges TSV."""
    results = []
    for mod in mod_types:
        sub = ranges_df[(ranges_df["mod_type"] == mod) & (ranges_df["biotype"] == "rRNA")]
        plat_params = {}
        for _, row in sub.iterrows():
            plat = row["platform"]
            param = row["parameter"]
            n_grid = row["n_grid_points"]
            if plat not in plat_params:
                plat_params[plat] = {}
            plat_params[plat][param] = n_grid
        
        # Only include if all 3 platforms have grids
        if len(plat_params) == 3:
            p1, p2, p3 = sorted(plat_params.keys())
            n1 = int(np.prod(list(plat_params[p1].values())))
            n2 = int(np.prod(list(plat_params[p2].values())))
            n3 = int(np.prod(list(plat_params[p3].values())))
            results.append({
                "biotype": "rRNA",
                "mod_type": mod,
                "comparison": f"{p1} × {p2} × {p3}",
                "p1_combos": n1,
                "p2_combos": n2,
                "p3_combos": n3,
                "total_combos": n1 * n2 * n3,
                "search_type": "3-way",
            })
    return pd.DataFrame(results)

for bio in ["rRNA", "polyA"]:
    df = rrna_ranges if bio == "rRNA" else polya_ranges
    space2 = compute_search_space_2way(df, bio)
    if not space2.empty:
        print(f"\n=== {bio} 2-Way Search Space ===")
        print(space2.to_string(index=False))

# 3-way only for rRNA Nm mods
space3 = compute_search_space_3way_from_ranges(rrna_ranges)
if not space3.empty:
    print(f"\n=== rRNA 3-Way Search Space (Nm mods) ===")
    print(space3.to_string(index=False))

# Combine and save
all_space = pd.concat([compute_search_space_2way(rrna_ranges, "rRNA"),
                       compute_search_space_2way(polya_ranges, "polyA"),
                       space3], ignore_index=True)
if not all_space.empty:
    all_space.to_csv(OUTDIR / "grid_search_space_sizes.tsv", sep="\t", index=False)
    print(f"\nSaved combined grid_search_space_sizes.tsv ({len(all_space)} rows)")
