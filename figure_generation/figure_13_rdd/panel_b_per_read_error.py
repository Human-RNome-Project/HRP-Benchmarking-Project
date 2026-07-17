#!/usr/bin/env python3
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)

GITHUB_ROOT = Path(__file__).resolve().parents[2]
READ_DIR = GITHUB_ROOT / "inputs"
OUTDIR = GITHUB_ROOT / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

EXCLUDE_NATIVE = {"HRP_A_013", "HRP_A_014"}

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.labelsize"] = 8
plt.rcParams["axes.titlesize"] = 8
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 7
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
mm = 1 / 25.4

all_files = sorted(READ_DIR.glob("*_per_read_error_rates.tsv"))
valid_samples = []
for f in all_files:
    sample = f.stem.replace("_per_read_error_rates", "")
    sample_id = "_".join(sample.split("_")[:3])
    if "IVT" not in sample and any(ex in sample_id for ex in EXCLUDE_NATIVE):
        continue
    valid_samples.append(sample)

native_samples = sorted([s for s in valid_samples if "IVT" not in s])
ivt_samples = sorted([s for s in valid_samples if "IVT" in s])
ordered_samples = native_samples + ivt_samples

sample_data = {}
for s in ordered_samples:
    f = READ_DIR / f"{s}_per_read_error_rates.tsv"
    df = pd.read_csv(f, sep="\t", usecols=["chrom", "total_rate"])
    df = df[df["chrom"] != "chrY"]
    sample_data[s] = df["total_rate"].dropna().values.astype(np.float64)

values = [sample_data[s] for s in ordered_samples]
positions = np.arange(1, len(ordered_samples) + 1)

fig_width = max(130 * mm, len(ordered_samples) * 8 * mm)
fig, ax = plt.subplots(figsize=(fig_width, 85 * mm))

bp = ax.boxplot(
    values,
    positions=positions,
    widths=0.6,
    patch_artist=True,
    showfliers=False,
    boxprops=dict(linewidth=1.0, color="#2C3E50"),
    medianprops=dict(linewidth=2.0, color="#2C3E50"),
    whiskerprops=dict(linewidth=1.0, color="#2C3E50"),
    capprops=dict(linewidth=1.0, color="#2C3E50"),
    zorder=2,
)

for patch, s in zip(bp["boxes"], ordered_samples):
    patch.set_facecolor("#bebcb7" if "IVT" in s else "#70AB5B")
    patch.set_alpha(0.8)

native_center = np.mean(np.arange(1, len(native_samples) + 1))
ivt_center = np.mean(np.arange(len(native_samples) + 1, len(ordered_samples) + 1))

ax.set_xticks([native_center, ivt_center])
ax.set_xticklabels(["Native", "IVT"], fontweight="bold", fontsize=10)
ax.tick_params(axis="x", length=0)
ax.set_ylabel("Per-read total error rate (%)", fontweight="bold")
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x*100:.1f}"))

# CRUCIAL FIX: Drop the Y-axis slightly below 0 so the 0% whisker floats above the spine
ax.set_ylim(bottom=-0.002, top=0.06) # Capping top at 6% to force separation from 0

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout(pad=1.5)
plt.savefig(OUTDIR / "panel_b_per_read_error.png", bbox_inches="tight", dpi=300)
plt.savefig(OUTDIR / "panel_b_per_read_error.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.close()
