#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

GITHUB_ROOT = Path(__file__).resolve().parents[2]
OUTDIR = GITHUB_ROOT / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(GITHUB_ROOT / "inputs" / "no_chrY_sub_summary.tsv", sep="\t")
# Use only the 2 specific IVT samples, drop the single 'IVT_noSupp' entry
df = df[df["sample"] != "IVT_noSupp"]

native_rates = df[df["group"] == "Native"]["total_mismatch_rate"].values * 100
ivt_rates = df[df["group"] == "IVT"]["total_mismatch_rate"].values * 100

print(f"Native: Mean={np.mean(native_rates):.4f}%, Median={np.median(native_rates):.4f}%")
print(f"IVT: Mean={np.mean(ivt_rates):.4f}%, Median={np.median(ivt_rates):.4f}%")

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.labelsize"] = 8
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 7
plt.rcParams["axes.labelweight"] = "bold"

mm = 1 / 25.4
fig, ax = plt.subplots(figsize=(40 * mm, 50 * mm))

bp = ax.boxplot(
    [native_rates, ivt_rates],
    positions=[1, 2],
    widths=0.6,
    patch_artist=True,
    showfliers=False,
    boxprops=dict(linewidth=1.0, color="#2C3E50"),
    medianprops=dict(linewidth=2.0, color="#2C3E50"),
    whiskerprops=dict(linewidth=1.0, color="#2C3E50"),
    capprops=dict(linewidth=1.0, color="#2C3E50"),
    zorder=1,
)
bp["boxes"][0].set_facecolor("#70AB5B") # Native
bp["boxes"][1].set_facecolor("#bebcb7") # IVT
bp["boxes"][0].set_alpha(0.8)
bp["boxes"][1].set_alpha(0.8)

np.random.seed(42)
ax.scatter(np.random.normal(1, 0.05, len(native_rates)), native_rates, color="#2C3E50", s=6, zorder=2, alpha=0.7)
ax.scatter(np.random.normal(2, 0.05, len(ivt_rates)), ivt_rates, color="#2C3E50", s=6, zorder=2, alpha=0.7)

ax.set_xticks([1, 2])
ax.set_xticklabels(["Native", "IVT"], fontweight="bold")
ax.set_ylabel("Substitution Error Rate (%)", fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(OUTDIR / "panel_d_substitution_error.png", dpi=300)
plt.savefig(OUTDIR / "panel_d_substitution_error.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.close()
