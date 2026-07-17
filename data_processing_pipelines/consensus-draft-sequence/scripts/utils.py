"""
Shared utilities for draft reference pipeline.

Key conventions:
  1. Nm is resolved to Am, Cm, Gm, or Um based on canonical base at position.
  2. "frequency" is relabeled as "level" in all outputs.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from paths import ROOT, INPUTS, OUTPUTS, THRESHOLDS, FILTERED, STATE, DRAFT, TIERED, TIERED_TRNA

# ── BED parsing ──────────────────────────────────────────────────────────────

def parse_bed(path, platform_label):
    """
    Generic BED parser.
    Returns DataFrame with columns:
      chrom, start, end, name, score, strand, coverage, level, platform
    Note: 'level' is used instead of 'frequency'.
    """
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            cov, lvl = 0.0, 0.0
            if len(p) >= 11:
                try:
                    cov = float(p[9]) if p[9] not in (".", "NA", "") else 0.0
                except:
                    pass
                try:
                    lvl = float(p[10]) if p[10] not in (".", "NA", "") else 0.0
                except:
                    pass
            rows.append({
                "chrom": p[0],
                "start": int(p[1]),
                "end": int(p[2]),
                "name": p[3],
                "score": float(p[4]) if p[4] not in (".", "NA", "") else 0.0,
                "strand": p[5],
                "coverage": cov,
                "level": lvl,  # renamed from frequency
                "platform": platform_label,
            })
    return pd.DataFrame(rows)

# ── Grid helper ──────────────────────────────────────────────────────────────

def make_grid(series, n_points=7):
    """Generate up to n_points quantile-based thresholds from a data series.
    
    Each threshold is snapped to the nearest actual value from the original data,
    ensuring all grid points exist in the input BED file.
    """
    vals = series.dropna().sort_values()
    if len(vals) == 0:
        return [0.0]
    unique_vals = np.unique(vals.values)
    if len(unique_vals) <= n_points:
        return [float(v) for v in unique_vals]
    
    q = np.linspace(0, 1, n_points)
    grid = np.quantile(vals, q)
    
    # Snap each quantile to the nearest actual data value
    sorted_vals = np.sort(vals.values)
    snapped = []
    for g in grid:
        idx = np.searchsorted(sorted_vals, g)
        candidates = []
        if idx > 0:
            candidates.append(sorted_vals[idx - 1])
        if idx < len(sorted_vals):
            candidates.append(sorted_vals[idx])
        if not candidates:
            candidates = [sorted_vals[0]]
        # Pick closest actual value
        best = min(candidates, key=lambda x: abs(x - g))
        snapped.append(best)
    
    # Deduplicate while preserving order
    seen = set()
    result = []
    for v in snapped:
        if v not in seen:
            seen.add(v)
            result.append(float(v))
    return result

# ── Level-aware column renaming ──────────────────────────────────────────────

LEVEL_COLUMNS = ["score", "coverage", "level"]

