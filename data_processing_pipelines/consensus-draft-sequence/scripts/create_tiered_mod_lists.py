#!/usr/bin/env python3
"""
Create tiered modification site lists for rRNA and polyA.

Tier1 uses the specific comparison intersections from grid search:
  - rRNA Nm mods: 3-way intersections (from rRNA_3way_Nm_best.tsv)
  - rRNA m5C: 3-way intersections (from rRNA_3way_m5C_best.tsv)
  - rRNA m6A: 3-way intersections (from rRNA_3way_m6A_best.tsv)
  - rRNA other non-Nm (e.g. Y): all 2-way intersections (from rRNA_grid_search_best.tsv)
  - polyA: 2-way Illumina-ONT intersections (from polyA_grid_search_best.tsv)
  tier1 = union of all grid-search-optimized comparison intersections

Tier2: site appears in >=2 raw combined BEDs, but NOT in tier1

Outputs:
  1. tiered_rRNA_shared_mod_types.tsv
  2. tiered_rRNA_only.tsv (+ .bed BED6 version)
  3. tiered_polyA.tsv

Table format: chr start end name tier strand in_ref ref_mod_type

v6 fixes: rRNA 2-way best-file parsing now preserves trailing empty columns
so Illumina-MassSpec comparisons are no longer skipped.
v7 change: m5C and m6A tier1 are now derived from 3-way grid search (same rule as Nm mods).
"""

import pandas as pd
from pathlib import Path
from collections import Counter
from paths import ROOT, INPUTS, OUTPUTS, THRESHOLDS, FILTERED, STATE, DRAFT, TIERED, TIERED_TRNA

BASE = ROOT
OUTDIR = TIERED
OUTDIR.mkdir(exist_ok=True)

# Mods that use 3-way grid-search intersections for tier1 (Nm + m5C + m6A)
THREE_WAY_MODS = {"Am", "Cm", "Gm", "Um", "m5C", "m6A"}

# ── Parse grid-search optimal thresholds ─────────────────────────────────────

def parse_rRNA_2way(path):
    """Parse rRNA_grid_search_best.tsv → dict mod → list of comparison specs."""
    comps = {}
    with open(path) as fh:
        header = next(fh).strip().split("\t")
        for line in fh:
            # Use rstrip('\n') only; strip() would drop trailing empty ONT columns
            # for Illumina-MassSpec comparisons.
            p = line.rstrip("\n").rstrip("\r").split("\t")
            if len(p) < 15:
                continue
            mod = p[0]
            comparison = p[1]
            plat_thresh = {}
            # Illumina
            if p[7] not in ("", "NA"):
                plat_thresh["illumina"] = (float(p[7]), float(p[8]), float(p[9]))
            # MassSpec
            if p[10] not in ("", "NA"):
                plat_thresh["massspec"] = (float(p[10]), float(p[11]), float(p[12]))
            # ONT
            if p[13] not in ("", "NA"):
                plat_thresh["ont"] = (float(p[13]), float(p[14]), float(p[15]))
            if mod not in comps:
                comps[mod] = []
            comps[mod].append({"comparison": comparison, "platforms": plat_thresh})
    return comps

def parse_rRNA_3way(path):
    """Parse a 3-way best TSV → dict mod → list with single 3-way spec."""
    comps = {}
    with open(path) as fh:
        header = next(fh).strip().split("\t")
        for line in fh:
            p = line.rstrip("\n").rstrip("\r").split("\t")
            if len(p) < 15:
                continue
            mod = p[0]
            plat_thresh = {
                "illumina": (float(p[7]), float(p[8]), float(p[9])),
                "ont": (float(p[10]), float(p[11]), float(p[12])),
                "massspec": (float(p[13]), float(p[14]), float(p[15])),
            }
            comps[mod] = [{"comparison": "Illumina-ONT-MassSpec", "platforms": plat_thresh}]
    return comps

def parse_polyA(path):
    """Parse polyA_grid_search_best.tsv → dict mod → list with single 2-way spec."""
    comps = {}
    with open(path) as fh:
        header = next(fh).strip().split("\t")
        for line in fh:
            p = line.rstrip("\n").rstrip("\r").split("\t")
            if len(p) < 12:
                continue
            mod = p[0]
            plat_thresh = {
                "illumina": (float(p[7]), float(p[8]), float(p[9])),
                "ont": (float(p[10]), float(p[11]), float(p[12])),
            }
            comps[mod] = [{"comparison": p[1], "platforms": plat_thresh}]
    return comps

print("Parsing grid-search thresholds...")
rRNA_2way = parse_rRNA_2way(THRESHOLDS / "rRNA_grid_search_best.tsv")
rRNA_3way_nm = parse_rRNA_3way(THRESHOLDS / "rRNA_3way_Nm_best.tsv")
rRNA_3way_m5c = parse_rRNA_3way(THRESHOLDS / "rRNA_3way_m5C_best.tsv")
rRNA_3way_m6a = parse_rRNA_3way(THRESHOLDS / "rRNA_3way_m6A_best.tsv")
polyA_gs = parse_polyA(THRESHOLDS / "polyA_grid_search_best.tsv")

# Merge: 3-way specs override 2-way entries for their mods
rRNA_thresholds = dict(rRNA_2way)
for mod, specs in rRNA_3way_nm.items():
    rRNA_thresholds[mod] = specs
for mod, specs in rRNA_3way_m5c.items():
    rRNA_thresholds[mod] = specs
for mod, specs in rRNA_3way_m6a.items():
    rRNA_thresholds[mod] = specs

print(f"  rRNA mods with thresholds: {sorted(rRNA_thresholds.keys())}")
print(f"  polyA mods with thresholds: {sorted(polyA_gs.keys())}")

# ── Load reference BED ───────────────────────────────────────────────────────
def load_ref_bed(path):
    ref = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            p = line.rstrip("\n").rstrip("\r").split("\t")
            if len(p) < 4:
                continue
            chrom, start, name = p[0], int(p[1]), p[3]
            ref[(chrom, start)] = name
    return ref

print("\nLoading reference BED...")
ref_map = load_ref_bed(INPUTS / "H.sapiens_rRNA_ref_mods.bed")
print(f"  Reference positions: {len(ref_map):,}")

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

# tRNA mod types for shared-mod-type table
ms_rows_all = load_bed_rows(OUTPUTS / "MS_rRNA_tRNA_harmonized.bed")
trna_mods = set()
for chrom, start, end, name, *_ in ms_rows_all:
    if chrom.startswith("hs_tRNA"):
        trna_mods.add(name)
print(f"  tRNA mod types in MS: {sorted(trna_mods)}")

# ── Tier computation ─────────────────────────────────────────────────────────
def compute_tiers(raw_beds, biotype_name, target_mods=None, threshold_map=None):
    """
    Compute tier1/tier2/tier3.
    raw_beds: list of (label, path, chrom_filter_func)
    threshold_map: dict mod → list of comparison specs
    """
    # Load all raw rows per platform
    raw_rows = {}
    raw_strands = {}
    for label, path, cfilter in raw_beds:
        rows = load_bed_rows(path, chrom_filter=cfilter)
        raw_rows[label] = rows
        for chrom, start, end, name, strand, *_ in rows:
            raw_strands[(chrom, start, end, name)] = strand
        print(f"  [{biotype_name}] Raw {label}: {len(rows):,} rows")

    # Collect all mod types seen
    all_mods = set()
    for rows in raw_rows.values():
        for row in rows:
            name = row[3]
            if target_mods is None or name in target_mods:
                all_mods.add(name)

    tier1_keys = set()
    tier2_keys = set()

    for mod in sorted(all_mods):
        # Build per-platform raw sets for this mod
        plat_raw_sets = {}
        for label, rows in raw_rows.items():
            s = set((r[0], r[1], r[2], r[3]) for r in rows if r[3] == mod)
            plat_raw_sets[label] = s

        # Tier1: grid-search comparison intersections
        if threshold_map and mod in threshold_map:
            for comp_spec in threshold_map[mod]:
                comparison = comp_spec["comparison"]
                plat_thresh = comp_spec["platforms"]
                # Apply thresholds per platform for this comparison
                plat_filtered_sets = {}
                for label, rows in raw_rows.items():
                    if label in plat_thresh:
                        filtered = apply_comparison_thresholds(rows, label, mod, plat_thresh)
                        plat_filtered_sets[label] = set((r[0], r[1], r[2], r[3]) for r in filtered)
                # Compute intersection across all platforms in this comparison
                if plat_filtered_sets:
                    inter = set.intersection(*plat_filtered_sets.values())
                    tier1_keys |= inter
                    if len(inter) > 0:
                        print(f"    {mod} {comparison}: {len(inter)} sites")

        # Tier2: raw overlap in 2+ platforms, not tier1
        all_raw = set()
        for s in plat_raw_sets.values():
            all_raw.update(s)
        for key in all_raw:
            n = sum(key in s for s in plat_raw_sets.values())
            if n >= 2 and key not in tier1_keys:
                tier2_keys.add(key)

    print(f"  [{biotype_name}] Tier1 (grid-search comparisons): {len(tier1_keys):,}")
    print(f"  [{biotype_name}] Tier2 (raw 2+ platform, not tier1): {len(tier2_keys):,}")

    # Combine into DataFrame
    seen = set()
    results = []
    def add(keys, tier):
        for chrom, start, end, name in keys:
            # exclude chrY
            if str(chrom).startswith("chrY"):
                continue
            key = (chrom, start, end, name)
            if key in seen:
                continue
            seen.add(key)
            strand = raw_strands.get(key, "+")
            in_ref = "TRUE" if (chrom, start) in ref_map else "FALSE"
            ref_mod_type = ref_map.get((chrom, start), "NA")
            results.append({
                "chr": chrom, "start": start, "end": end, "name": name,
                "tier": tier, "strand": strand, "in_ref": in_ref,
                "ref_mod_type": ref_mod_type,
            })

    add(tier1_keys, "tier1")
    add(tier2_keys, "tier2")

    df = pd.DataFrame(results)
    if not df.empty:
        print(f"  [{biotype_name}] Total tiered sites: {len(df):,}")
        for t, c in df["tier"].value_counts().sort_index().items():
            print(f"    {t}: {c}")
    return df

# ── Filters ──────────────────────────────────────────────────────────────────
def is_rrna(chrom):
    return chrom.startswith("hs_rRNA")

def is_polya(chrom):
    return not chrom.startswith("hs_rRNA") and not chrom.startswith("hs_tRNA")

# ── rRNA ─────────────────────────────────────────────────────────────────────
print("\n=== Processing rRNA ===")

rRNA_raw = [
    ("illumina", INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed", is_rrna),
    ("massspec", OUTPUTS / "MS_rRNA_tRNA_harmonized.bed", is_rrna),
    ("ont", INPUTS / "ONT_polyARNA_rRNA_combined.filtered.bed", is_rrna),
]

df_rRNA_all = compute_tiers(rRNA_raw, "rRNA", threshold_map=rRNA_thresholds)

out_rRNA_only = OUTDIR / "tiered_rRNA_only.tsv"
df_rRNA_all.to_csv(out_rRNA_only, sep="\t", index=False)
print(f"  Saved: {out_rRNA_only}")

# Also write a BED6 version: name = "mod|tier"
bed_rows = []
for _, row in df_rRNA_all.iterrows():
    bed_rows.append({
        "chrom": row["chr"],
        "start": row["start"],
        "end": row["end"],
        "name": f"{row['name']}|{row['tier']}",
        "score": 0,
        "strand": row["strand"],
    })
df_bed = pd.DataFrame(bed_rows)
out_rRNA_bed = OUTDIR / "tiered_rRNA_only.bed"
df_bed.to_csv(out_rRNA_bed, sep="\t", index=False, header=False)
print(f"  Saved: {out_rRNA_bed}")

# Shared mod types
rRNA_mods = set(df_rRNA_all["name"].unique())
shared_mods = sorted(rRNA_mods & trna_mods)
print(f"  Shared mod types (rRNA ∩ tRNA MS): {shared_mods}")

df_rRNA_shared = df_rRNA_all[df_rRNA_all["name"].isin(shared_mods)].copy()
out_rRNA_shared = OUTDIR / "tiered_rRNA_shared_mod_types.tsv"
df_rRNA_shared.to_csv(out_rRNA_shared, sep="\t", index=False)
print(f"  Shared-mod-type tiered sites: {len(df_rRNA_shared):,}")
print(f"  Saved: {out_rRNA_shared}")

# ── polyA ────────────────────────────────────────────────────────────────────
print("\n=== Processing polyA ===")

target_polya_mods = {"m6A", "m5C", "Y", "I"}

polyA_raw = [
    ("illumina", INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed", is_polya),
    ("ont", INPUTS / "ONT_polyARNA_rRNA_combined.filtered.bed", is_polya),
]

df_polyA = compute_tiers(polyA_raw, "polyA", target_mods=target_polya_mods, threshold_map=polyA_gs)

# Drop reference columns for polyA (no rRNA reference)
df_polyA = df_polyA.drop(columns=["in_ref", "ref_mod_type"], errors="ignore")

out_polyA = OUTDIR / "tiered_polyA.tsv"
df_polyA.to_csv(out_polyA, sep="\t", index=False)
print(f"  Saved: {out_polyA}")

# Per-mod polyA files
for mod in sorted(target_polya_mods):
    df_mod = df_polyA[df_polyA["name"] == mod].copy()
    out_mod = OUTDIR / f"tiered_polyA_{mod}.tsv"
    df_mod.to_csv(out_mod, sep="\t", index=False)
    print(f"  Saved: {out_mod} ({len(df_mod):,} rows)")

print("\n=== All done ===")
