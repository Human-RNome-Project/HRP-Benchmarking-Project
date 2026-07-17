#!/usr/bin/env python3
"""
MassSpec BED harmonization script for tRNA + rRNA.

Resolves ambiguous MS modification calls (mA?/mxA, mC?/mxC, mG?/mxG, mU?/mxU)
to specific modifications by cross-referencing with:
  1. MS-self (other non-ambiguous mods at same position)
  2. Reference BED (H.sapiens_rRNA_ref_mods.bed) — rRNA only
  3. Illumina combined BED — rRNA and tRNA
  4. ONT combined BED — rRNA and tRNA

Priority for ambiguous resolution: MS-self > Reference (rRNA) > Illumina > ONT.
Critical fix for mU?/mxU: only blocks resolution to pseU/Y/psU (pseudouridine);
allows resolution to Um (2'-O-methyluridine).

Input:  inputs/MS_rRNA_tRNA.bed (rmchrY version, mxA/mxC/mxG/mxU ambiguous names)
Output: outputs/MS_rRNA_tRNA_harmonized.bed
"""

import pandas as pd
from pathlib import Path
from collections import Counter
from paths import INPUTS, OUTPUTS

# ── Paths ────────────────────────────────────────────────────────────────────
REF_BED = INPUTS / "H.sapiens_rRNA_ref_mods.bed"
ILLUMINA_BED = INPUTS / "Illumina_combined_polyARNA_tRNA_rRNA.bed"
ONT_BED = INPUTS / "ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed"
MS_INPUT = INPUTS / "MS_rRNA_tRNA.bed"
MS_OUTPUT = OUTPUTS / "MS_rRNA_tRNA_harmonized.bed"

# ── Helper functions ─────────────────────────────────────────────────────────

AMBIGUOUS_NAMES = {"ma?", "mc?", "mg?", "mu?", "mxa", "mxc", "mxg", "mxu"}

# Canonical ambiguous symbols to write for unresolved calls in column 4.
# The raw last-two columns are preserved unchanged from the input BED.
AMBIGUOUS_OLD_TO_NEW = {
    "ma?": "mA?",
    "mc?": "mC?",
    "mg?": "mG?",
    "mu?": "mU?",
    "mxa": "mA?",
    "mxc": "mC?",
    "mxg": "mG?",
    "mxu": "mU?",
}


def is_ambiguous(name):
    return str(name).strip().lower() in AMBIGUOUS_NAMES


def canonical_ambiguous(name):
    """Return the canonical symbol for an ambiguous name, or the name itself."""
    return AMBIGUOUS_OLD_TO_NEW.get(str(name).strip().lower(), name)


def parse_bed(path, platform_label):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            cov, freq = 0.0, 0.0
            if len(p) >= 11:
                try:
                    cov = float(p[9]) if p[9] not in (".", "NA", "") else 0.0
                except Exception:
                    pass
                try:
                    freq = float(p[10]) if p[10] not in (".", "NA", "") else 0.0
                except Exception:
                    pass
            rows.append({
                "chrom": p[0],
                "start": int(p[1]),
                "end": int(p[2]),
                "raw_name": p[3],
                "name": p[3],
                "score": float(p[4]) if p[4] not in (".", "NA", "") else 0.0,
                "strand": p[5],
                "coverage": cov,
                "frequency": freq,
                "platform": platform_label,
            })
    return pd.DataFrame(rows)


def build_position_map(df, chrom_prefix):
    """Build (chrom, start) -> name map for rows whose chrom starts with prefix.

    Skips ambiguous mod names so that ambiguous calls cannot resolve other
    ambiguous calls.
    """
    sub = df[df["chrom"].astype(str).str.startswith(chrom_prefix)].copy()
    pos_map = {}
    for _, row in sub.iterrows():
        if is_ambiguous(row["name"]):
            continue
        key = (row["chrom"], row["start"])
        if key not in pos_map:
            pos_map[key] = row["name"]
    return pos_map


def build_self_map(df, chrom_prefix):
    """Build MS self-map of specific (non-ambiguous) mods for given chrom prefix."""
    self_map = {}
    sub = df[df["chrom"].astype(str).str.startswith(chrom_prefix)].copy()
    for _, row in sub.iterrows():
        if is_ambiguous(row["raw_name"]):
            continue
        key = (row["chrom"], row["start"])
        if key not in self_map:
            self_map[key] = row["name"]
    return self_map


# ── Load cross-reference data ────────────────────────────────────────────────
print("Loading reference BED...")
ref_df = parse_bed(REF_BED, "Reference")
ref_map = build_position_map(ref_df, "hs_rRNA")
print(f"  Reference rRNA: {len(ref_map):,} unique positions")

print(f"Loading Illumina BED: {ILLUMINA_BED}")
ill_df = parse_bed(ILLUMINA_BED, "Illumina")
ill_rrna_map = build_position_map(ill_df, "hs_rRNA")
ill_trna_map = build_position_map(ill_df, "hs_tRNA")
print(f"  Illumina rRNA: {len(ill_rrna_map):,} positions, tRNA: {len(ill_trna_map):,} positions")

print(f"Loading ONT BED: {ONT_BED}")
ont_df = parse_bed(ONT_BED, "ONT")
ont_rrna_map = build_position_map(ont_df, "hs_rRNA")
ont_trna_map = build_position_map(ont_df, "hs_tRNA")
print(f"  ONT rRNA: {len(ont_rrna_map):,} positions, tRNA: {len(ont_trna_map):,} positions")

print(f"Loading MassSpec BED: {MS_INPUT}")
ms_df = parse_bed(MS_INPUT, "MassSpec")
ms_rrna = ms_df[ms_df["chrom"].astype(str).str.startswith("hs_rRNA")].copy()
ms_trna = ms_df[ms_df["chrom"].astype(str).str.startswith("hs_tRNA")].copy()
print(f"  MassSpec total: {len(ms_df)} rows, rRNA: {len(ms_rrna)} rows, tRNA: {len(ms_trna)} rows")
print(f"  rRNA mod types: {sorted(ms_rrna['name'].unique())}")
print(f"  tRNA mod types: {sorted(ms_trna['name'].unique())}")

# Build MS self-maps
ms_rrna_self_map = build_self_map(ms_df, "hs_rRNA")
ms_trna_self_map = build_self_map(ms_df, "hs_tRNA")
print(f"  MS rRNA self-map: {len(ms_rrna_self_map)} positions with specific mods")
print(f"  MS tRNA self-map: {len(ms_trna_self_map)} positions with specific mods")


# ── Resolve ambiguous mods ───────────────────────────────────────────────────
def resolve_ambiguous(ms_subset, self_map, ref_map, ill_map, ont_map, biotype):
    """Resolve ambiguous names in an MS subset and return a copy with 'name' updated."""
    resolution_stats = Counter()
    resolution_detail = []
    ms_harmonized = ms_subset.copy()

    for idx in ms_harmonized.index:
        raw_name = ms_harmonized.at[idx, "raw_name"]
        if not is_ambiguous(raw_name):
            continue

        chrom = ms_harmonized.at[idx, "chrom"]
        start = ms_harmonized.at[idx, "start"]
        key = (chrom, start)

        resolved = False
        resolved_by = None
        new_name = None

        # Priority: MS-self > Reference (rRNA only) > Illumina > ONT
        if key in self_map and self_map[key]:
            new_name = self_map[key]
            resolved_by = "MS-self"
            resolved = True
        elif biotype == "rRNA" and key in ref_map and ref_map[key]:
            new_name = ref_map[key]
            resolved_by = "Reference"
            resolved = True
        elif key in ill_map and ill_map[key]:
            new_name = ill_map[key]
            resolved_by = "Illumina"
            resolved = True
        elif key in ont_map and ont_map[key]:
            new_name = ont_map[key]
            resolved_by = "ONT"
            resolved = True

        # Critical fix: only block mU?/mxU -> pseU/Y/psU, allow -> Um
        if resolved and raw_name.lower() in ("mu?", "mxu") and new_name in ("Y", "pseU", "psU"):
            resolved = False
            resolved_by = "Unresolved"

        if resolved:
            ms_harmonized.at[idx, "name"] = new_name
            resolution_stats[resolved_by] += 1
        else:
            resolution_stats["Unresolved"] += 1
            resolved_by = "Unresolved"
            # Write unresolved ambiguous calls with the canonical symbol
            new_name = canonical_ambiguous(raw_name)
            ms_harmonized.at[idx, "name"] = new_name

        resolution_detail.append({
            "chrom": chrom,
            "start": start,
            "original": raw_name,
            "resolved_to": new_name,
            "resolved_by": resolved_by,
        })

    return ms_harmonized, resolution_stats, resolution_detail


print("\n=== Resolving rRNA ambiguous mods ===")
ms_rrna_harm, rrna_stats, rrna_detail = resolve_ambiguous(
    ms_rrna, ms_rrna_self_map, ref_map, ill_rrna_map, ont_rrna_map, "rRNA"
)
for source, count in rrna_stats.most_common():
    print(f"  {source}: {count}")

print("\n=== Resolving tRNA ambiguous mods ===")
ms_trna_harm, trna_stats, trna_detail = resolve_ambiguous(
    ms_trna, ms_trna_self_map, None, ill_trna_map, ont_trna_map, "tRNA"
)
for source, count in trna_stats.most_common():
    print(f"  {source}: {count}")

print(f"\nBefore harmonization rRNA mod types: {sorted(ms_rrna['name'].unique())}")
print(f"After  harmonization rRNA mod types: {sorted(ms_rrna_harm['name'].unique())}")
print(f"\nBefore harmonization tRNA mod types: {sorted(ms_trna['name'].unique())}")
print(f"After  harmonization tRNA mod types: {sorted(ms_trna_harm['name'].unique())}")

# Combine lookup tables
harmonize_lookup = {}
for _, row in ms_rrna_harm.iterrows():
    harmonize_lookup[(row["chrom"], row["start"], row["raw_name"])] = row["name"]
for _, row in ms_trna_harm.iterrows():
    harmonize_lookup[(row["chrom"], row["start"], row["raw_name"])] = row["name"]

# ── Write output BED ─────────────────────────────────────────────────────────
OUTPUTS.mkdir(parents=True, exist_ok=True)
written = 0
changed = 0
with open(MS_INPUT) as fin, open(MS_OUTPUT, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            if line.startswith("#sequencing_plattform="):
                fout.write("#sequencing_plattform=MS\n")
            else:
                fout.write(line)
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        chrom, start, raw_name = p[0], int(p[1]), p[3]
        key = (chrom, start, raw_name)
        if key in harmonize_lookup and harmonize_lookup[key] != raw_name:
            p[3] = harmonize_lookup[key]
            changed += 1

        fout.write("\t".join(str(x) for x in p) + "\n")
        written += 1

print(f"\nSaved harmonized BED: {MS_OUTPUT}")
print(f"  Total rows written: {written}")
print(f"  Names changed: {changed}")

# ── Verification ─────────────────────────────────────────────────────────────
verify_df = parse_bed(MS_OUTPUT, "MassSpec")
verify_rrna = verify_df[verify_df["chrom"].astype(str).str.startswith("hs_rRNA")].copy()
verify_trna = verify_df[verify_df["chrom"].astype(str).str.startswith("hs_tRNA")].copy()

print("\n=== Original vs Harmonized (rRNA only) ===")
orig_counts = ms_rrna["name"].value_counts().sort_index()
harm_counts = verify_rrna["name"].value_counts().sort_index()
compare = pd.DataFrame({"original": orig_counts, "harmonized": harm_counts}).fillna(0).astype(int)
print(compare.to_string())

print("\n=== Original vs Harmonized (tRNA only) ===")
orig_counts = ms_trna["name"].value_counts().sort_index()
harm_counts = verify_trna["name"].value_counts().sort_index()
compare = pd.DataFrame({"original": orig_counts, "harmonized": harm_counts}).fillna(0).astype(int)
print(compare.to_string())
