#!/usr/bin/env python3
"""
Create tiered tRNA site lists.

Tier1 uses the grid-search-optimized comparison intersections:
  - tRNA 3-way mods (Cm, Gm, Um, m5C): 3-way intersections (from tRNA_3way_best.tsv)
  - tRNA 2-way mods: all 2-way intersections (from tRNA_grid_search_best.tsv)
  tier1 = union of all grid-search-optimized comparison intersections

Tier2: site appears in >=2 raw combined BEDs, but NOT in tier1.

Inputs:
  - tRNA_3way_best.tsv
  - tRNA_grid_search_best.tsv
  - Illumina_combined_polyARNA_tRNA_rRNA.bed
  - ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed
  - MS_rRNA_tRNA_harmonized.bed

Outputs:
  - tiered_tRNA/tiered_tRNA.tsv
  - tiered_tRNA/tiered_tRNA_counts.tsv
"""

import pandas as pd
from pathlib import Path
from collections import Counter
from paths import ROOT, INPUTS, OUTPUTS, THRESHOLDS, FILTERED, STATE, DRAFT, TIERED, TIERED_TRNA

BASE = ROOT
OUTDIR = TIERED_TRNA
NEW_RUN = THRESHOLDS
OUTDIR.mkdir(exist_ok=True, parents=True)

# tRNA mods that use 3-way grid-search intersections for tier1
THREE_WAY_MODS = {"Cm", "Gm", "Um", "m5C"}

# ── Parse grid-search optimal thresholds ─────────────────────────────────────

def parse_tRNA_2way(path):
    """Parse tRNA_grid_search_best.tsv → dict mod → list of comparison specs."""
    comps = {}
    with open(path) as fh:
        header = next(fh).strip().split("\t")
        for line in fh:
            p = line.rstrip("\n").rstrip("\r").split("\t")
            if len(p) < 16:
                continue
            mod = p[0]
            comparison = p[1]
            plat_thresh = {}
            # Illumina
            if p[7] not in ("", "NA"):
                plat_thresh["illumina"] = (float(p[7]), float(p[8]), float(p[9]))
            # ONT
            if p[10] not in ("", "NA"):
                plat_thresh["ont"] = (float(p[10]), float(p[11]), float(p[12]))
            # MassSpec
            if p[13] not in ("", "NA"):
                plat_thresh["massspec"] = (float(p[13]), float(p[14]), float(p[15]))
            jaccard = float(p[2]) if p[2] not in ("", "NA") else "NA"
            if mod not in comps:
                comps[mod] = []
            comps[mod].append({"comparison": comparison, "platforms": plat_thresh, "jaccard": jaccard})
    return comps

def parse_tRNA_3way(path):
    """Parse tRNA_3way_best.tsv → dict mod → list with single 3-way spec."""
    comps = {}
    with open(path) as fh:
        header = next(fh).strip().split("\t")
        for line in fh:
            p = line.rstrip("\n").rstrip("\r").split("\t")
            if len(p) < 16:
                continue
            mod = p[0]
            plat_thresh = {
                "illumina": (float(p[7]), float(p[8]), float(p[9])),
                "ont": (float(p[10]), float(p[11]), float(p[12])),
                "massspec": (float(p[13]), float(p[14]), float(p[15])),
            }
            jaccard = float(p[1]) if p[1] not in ("", "NA") else "NA"
            comps[mod] = [{"comparison": "Illumina-ONT-MassSpec", "platforms": plat_thresh, "best_jaccard": jaccard}]
    return comps

print("Parsing grid-search thresholds...")
tRNA_2way = parse_tRNA_2way(NEW_RUN / "tRNA_grid_search_best.tsv")
tRNA_3way = parse_tRNA_3way(NEW_RUN / "tRNA_3way_best.tsv")

# Merge: 3-way specs override 2-way entries for their mods
tRNA_thresholds = dict(tRNA_2way)
for mod, specs in tRNA_3way.items():
    tRNA_thresholds[mod] = specs

# Also keep jaccard lookup for each mod/comparison
jaccard_lookup = {}
for mod, specs in tRNA_2way.items():
    for spec in specs:
        jaccard_lookup[(mod, spec["comparison"])] = spec.get("jaccard", "NA")
for mod, specs in tRNA_3way.items():
    for spec in specs:
        jaccard_lookup[(mod, spec["comparison"])] = spec.get("best_jaccard", "NA")

print(f"  tRNA mods with thresholds: {sorted(tRNA_thresholds.keys())}")

# ── Helper: load BED rows with raw values ────────────────────────────────────
def load_bed_rows(path, chrom_filter=None):
    """Load BED file, return list of (chrom, start, end, name, strand, score, coverage, level)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 11:
                continue
            chrom, start, end, name = p[0], int(p[1]), int(p[2]), p[3]
            strand = p[5] if len(p) > 5 else "+"
            try:
                score = float(p[4]) if p[4] not in (".", "NA", "") else 0.0
                cov = float(p[9]) if p[9] not in (".", "NA", "") else 0.0
                lvl = float(p[10]) if p[10] not in (".", "NA", "") else 0.0
            except ValueError:
                continue
            if chrom_filter and not chrom_filter(chrom):
                continue
            rows.append((chrom, start, end, name, strand, score, cov, lvl))
    return rows

def is_trna(chrom):
    return str(chrom).startswith("hs_tRNA") or str(chrom).startswith("hs_mttRNA")

# ── Apply thresholds for a specific comparison ───────────────────────────────
def apply_comparison_thresholds(rows, plat_key, mod_type, thresh_spec):
    """Filter rows by grid-search optimal thresholds for one platform/mod/comparison."""
    if plat_key not in thresh_spec:
        return []
    score_th, cov_th, lvl_th = thresh_spec[plat_key]
    result = []
    for chrom, start, end, name, strand, score, cov, lvl in rows:
        if name != mod_type:
            continue
        if score < score_th:
            continue
        if cov < cov_th:
            continue
        if lvl < lvl_th:
            continue
        result.append((chrom, start, end, name, strand))
    return result

# ── Tier computation ─────────────────────────────────────────────────────────
def compute_tiers(raw_beds, threshold_map):
    """
    Compute tier1/tier2 for tRNA.
    raw_beds: list of (label, path)
    threshold_map: dict mod → list of comparison specs
    """
    # Load all raw rows per platform
    raw_rows = {}
    raw_strands = {}
    representative_row = {}
    for label, path in raw_beds:
        rows = load_bed_rows(path, chrom_filter=is_trna)
        raw_rows[label] = rows
        for chrom, start, end, name, strand, *_ in rows:
            key = (chrom, start, name)
            raw_strands[key] = strand
            if key not in representative_row:
                representative_row[key] = (chrom, start, end, name, strand)
        print(f"  Raw {label}: {len(rows):,} tRNA rows")

    # Collect all mod types seen
    all_mods = set()
    for rows in raw_rows.values():
        for row in rows:
            all_mods.add(row[3])

    # tier1/tier2 metadata dicts: key (chrom, start, name) -> metadata dict
    tier1_meta = {}
    tier2_meta = {}

    for mod in sorted(all_mods):
        # Build per-platform raw position sets for this mod (chrom, start, name)
        plat_raw_sets = {}
        for label, rows in raw_rows.items():
            s = set((r[0], r[1], r[3]) for r in rows if r[3] == mod)
            plat_raw_sets[label] = s

        # Tier1: grid-search comparison intersections
        if mod in threshold_map:
            for comp_spec in threshold_map[mod]:
                comparison = comp_spec["comparison"]
                plat_thresh = comp_spec["platforms"]
                # Apply thresholds per platform for this comparison
                plat_filtered_sets = {}
                for label, rows in raw_rows.items():
                    if label in plat_thresh:
                        filtered = apply_comparison_thresholds(rows, label, mod, plat_thresh)
                        plat_filtered_sets[label] = set((r[0], r[1], r[3]) for r in filtered)
                # Compute intersection across all platforms in this comparison
                if plat_filtered_sets:
                    inter = set.intersection(*plat_filtered_sets.values())
                    jaccard = jaccard_lookup.get((mod, comparison), "NA")
                    for key in inter:
                        if key not in tier1_meta:
                            tier1_meta[key] = {
                                "comparison": comparison,
                                "jaccard": jaccard,
                                "platforms": sorted(plat_filtered_sets.keys()),
                            }
                    if len(inter) > 0:
                        print(f"    {mod} {comparison}: {len(inter)} sites (J={jaccard})")

        # Tier2: raw overlap in 2+ platforms, not tier1
        all_raw = set()
        for s in plat_raw_sets.values():
            all_raw.update(s)
        for key in all_raw:
            if key in tier1_meta:
                continue
            supporting = sorted([label for label, s in plat_raw_sets.items() if key in s])
            n = len(supporting)
            if n >= 2 and key not in tier2_meta:
                tier2_meta[key] = {
                    "comparison": "raw_overlap",
                    "jaccard": "NA",
                    "platforms": supporting,
                }

    print(f"  Tier1 (grid-search comparisons): {len(tier1_meta):,}")
    print(f"  Tier2 (raw 2+ platform, not tier1): {len(tier2_meta):,}")

    # Build output DataFrame
    seen = set()
    results = []
    def add(meta_dict, tier):
        for key, meta in meta_dict.items():
            if key in seen:
                continue
            seen.add(key)
            chrom, start, name = key
            rep = representative_row.get(key)
            if rep is None:
                end = start + 1
                strand = raw_strands.get(key, "+")
            else:
                chrom, start, end, name, strand = rep
            results.append({
                "chr": chrom, "start": start, "end": end, "name": name,
                "tier": tier, "strand": strand,
                "platforms": ",".join(meta["platforms"]),
                "comparison": meta["comparison"],
                "jaccard": meta["jaccard"],
            })

    add(tier1_meta, "tier1")
    add(tier2_meta, "tier2")

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(["chr", "start", "name", "tier"])
        print(f"  Total tiered sites: {len(df):,}")
        for t, c in df["tier"].value_counts().sort_index().items():
            print(f"    {t}: {c}")
    return df

# ── tRNA raw beds ────────────────────────────────────────────────────────────
print("\n=== Processing tRNA ===")

tRNA_raw = [
    ("illumina", INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed"),
    ("ont", INPUTS / "ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed"),
    ("massspec", OUTPUTS / "MS_rRNA_tRNA_harmonized.bed"),
]

df_tRNA = compute_tiers(tRNA_raw, threshold_map=tRNA_thresholds)

out_tRNA = OUTDIR / "tiered_tRNA.tsv"
df_tRNA.to_csv(out_tRNA, sep="\t", index=False)
print(f"\nSaved: {out_tRNA}")

# Per-mod counts
counts = []
all_mods = sorted(df_tRNA["name"].unique())
for mod in all_mods:
    sub = df_tRNA[df_tRNA["name"] == mod]
    tier1 = len(sub[sub["tier"] == "tier1"])
    tier2 = len(sub[sub["tier"] == "tier2"])
    counts.append({
        "mod": mod,
        "tier1": tier1,
        "tier2": tier2,
        "total": tier1 + tier2,
    })
counts_df = pd.DataFrame(counts)
counts_path = OUTDIR / "tiered_tRNA_counts.tsv"
counts_df.to_csv(counts_path, sep="\t", index=False)
print(f"Saved: {counts_path}")
print(counts_df.to_string(index=False))

print("\nDone.")
