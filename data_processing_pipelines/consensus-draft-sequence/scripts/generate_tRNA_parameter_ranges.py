#!/usr/bin/env python3
"""
Generate quantile-based parameter ranges for tRNA grid search.

Input BEDs:
  - Illumina_combined_polyARNA_tRNA_rRNA.bed
  - ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed
  - MS_rRNA_tRNA_harmonized.bed

Output:
  - tRNA_grid_search_parameter_ranges.tsv
  - tRNA_grid_search_space_sizes.tsv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import parse_bed, make_grid
from paths import ROOT, INPUTS, OUTPUTS, THRESHOLDS, FILTERED, STATE, DRAFT, TIERED, TIERED_TRNA

BASE = Path(__file__).parent.parent
OUTDIR = Path(__file__).parent

ILLUMINA_BED = INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed"
ONT_BED = INPUTS / "ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed"
MS_BED = OUTPUTS / "MS_rRNA_tRNA_harmonized.bed"

for label, path in [("Illumina", ILLUMINA_BED), ("ONT", ONT_BED), ("MassSpec", MS_BED)]:
    exists = "✓" if path.exists() else "✗ MISSING"
    print(f"{exists} {label}: {path}")

def load_tRNA(bed_path, platform_label):
    """Load BED and return tRNA rows only."""
    df = parse_bed(bed_path, platform_label)
    chrom = df["chrom"].astype(str)
    trna = df[chrom.str.startswith("hs_tRNA") | chrom.str.startswith("hs_mttRNA")].copy()
    return trna

print("\n=== Loading tRNA BEDs ===")
trna_data = {
    "Illumina": load_tRNA(ILLUMINA_BED, "Illumina"),
    "ONT": load_tRNA(ONT_BED, "ONT"),
    "MassSpec": load_tRNA(MS_BED, "MassSpec"),
}
for plat, df in trna_data.items():
    print(f"  [{plat}] {len(df):,} tRNA rows, mods: {sorted(df['name'].unique())}")

def compute_param_ranges(df, platform, biotype):
    """Compute parameter ranges per mod type. MassSpec has no coverage."""
    records = []
    for mod in sorted(df["name"].unique()):
        sub = df[df["name"] == mod]
        params = ["score", "coverage", "level"]
        if platform == "MassSpec":
            params = ["score", "level"]
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

print("\n=== Computing tRNA parameter ranges ===")
trna_ranges = pd.concat([
    compute_param_ranges(d, p, "tRNA") for p, d in trna_data.items() if not d.empty
], ignore_index=True)
trna_ranges.sort_values(["mod_type", "platform", "parameter"]).to_csv(
    OUTDIR / "tRNA_grid_search_parameter_ranges.tsv", sep="\t", index=False
)
print(f"Saved tRNA_grid_search_parameter_ranges.tsv ({len(trna_ranges)} rows)")

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

def compute_search_space_3way_from_ranges(ranges_df, biotype):
    """Compute 3-way search space for mods present in all 3 platforms."""
    results = []
    mod_types = sorted(ranges_df[ranges_df["biotype"] == biotype]["mod_type"].unique())
    for mod in mod_types:
        sub = ranges_df[(ranges_df["mod_type"] == mod) & (ranges_df["biotype"] == biotype)]
        plat_params = {}
        for _, row in sub.iterrows():
            plat = row["platform"]
            param = row["parameter"]
            n_grid = row["n_grid_points"]
            if plat not in plat_params:
                plat_params[plat] = {}
            plat_params[plat][param] = n_grid

        if len(plat_params) == 3:
            p1, p2, p3 = sorted(plat_params.keys())
            n1 = int(np.prod(list(plat_params[p1].values())))
            n2 = int(np.prod(list(plat_params[p2].values())))
            n3 = int(np.prod(list(plat_params[p3].values())))
            results.append({
                "biotype": biotype,
                "mod_type": mod,
                "comparison": f"{p1} × {p2} × {p3}",
                "p1_combos": n1,
                "p2_combos": n2,
                "p3_combos": n3,
                "total_combos": n1 * n2 * n3,
                "search_type": "3-way",
            })
    return pd.DataFrame(results)

space2 = compute_search_space_2way(trna_ranges, "tRNA")
if not space2.empty:
    print("\n=== tRNA 2-Way Search Space ===")
    print(space2.to_string(index=False))

space3 = compute_search_space_3way_from_ranges(trna_ranges, "tRNA")
if not space3.empty:
    print("\n=== tRNA 3-Way Search Space ===")
    print(space3.to_string(index=False))

all_space = pd.concat([space2, space3], ignore_index=True)
if not all_space.empty:
    all_space.to_csv(OUTDIR / "tRNA_grid_search_space_sizes.tsv", sep="\t", index=False)
    print(f"\nSaved tRNA_grid_search_space_sizes.tsv ({len(all_space)} rows)")

print("\nDone.")
