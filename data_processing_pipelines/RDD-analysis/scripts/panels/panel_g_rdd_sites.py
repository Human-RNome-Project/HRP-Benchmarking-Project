import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import os

GITHUB_ROOT = Path(__file__).resolve().parents[2]
JSON_PATH = GITHUB_ROOT / "inputs" / "genomewide_rdd.json"
OUT_DIR = GITHUB_ROOT / "figures"
os.makedirs(OUT_DIR, exist_ok=True)

CONSENSUS_MOD_COUNT = 17565

# Nature / Human RNome guidelines
plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

MM = 1 / 25.4
SINGLE_COL_W = 89 * MM
FIG_H = 75 * MM

with open(JSON_PATH) as fh:
    rdd = json.load(fh)

thresholds_str = sorted(rdd["thresh"].keys(), key=lambda k: float(k))
counts = np.array([rdd["thresh"][t] for t in thresholds_str])
thresholds = np.array([float(t) for t in thresholds_str])
thresholds_pct = thresholds * 100

fig, ax = plt.subplots(figsize=(SINGLE_COL_W, FIG_H))

ax.plot(thresholds_pct, counts, "o-", color="#1E88E5",
        markersize=4, linewidth=1.2, markeredgecolor="white", markeredgewidth=0.4,
        zorder=3)

# Horizontal dashed line at consensus modification count
ax.axhline(CONSENSUS_MOD_COUNT, color="#000000", linestyle="--", linewidth=0.8, zorder=2)

# Add text above the dashed line (moved to the left)
ax.text(
    x=min(thresholds_pct),
    y=CONSENSUS_MOD_COUNT,
    s=f"{CONSENSUS_MOD_COUNT:,} sites",
    va="bottom",
    ha="left",
    fontsize=7,
    color="#000000",
    zorder=4
)

# Find 15% point
if 15.0 in thresholds_pct:
    idx_15 = list(thresholds_pct).index(15.0)
    count_15 = counts[idx_15]

    # Highlight 15% threshold with a red circle around the blue dot
    ax.scatter([15.0], [count_15], edgecolors="red", facecolors="none", s=60, linewidths=1.2, zorder=5)
    ax.annotate(
        f"{count_15:,} sites",
        xy=(15.0, count_15),
        xytext=(8, 8),
        textcoords="offset points",
        fontsize=7,
        color="red",
        zorder=6
    )

ax.set_xlabel("RDD level threshold (%)")
ax.set_ylabel("Number of RDD sites")
ax.set_yscale("log")
ax.grid(axis="both", alpha=0.25, linewidth=0.5)

plt.tight_layout()
fig.savefig(f"{OUT_DIR}/panel_g_rdd_sites.pdf", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT_DIR}/panel_g_rdd_sites.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved to {OUT_DIR}/panel_g_rdd_sites.png")
