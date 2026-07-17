#!/usr/bin/env python3
"""
Substitution spectrum for strictly 13 native samples (no IVT), one bar per substitution type.

Bars are coloured by reference base and show the median rate across strictly
13 native samples, with IQR error bars.

Outputs to a distinct new file to avoid overwriting.
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ────────────────────────────────────────────────────────────────────
GITHUB_ROOT = Path(__file__).resolve().parents[2]
NATIVE_DIR = GITHUB_ROOT / "inputs"
OUTDIR = GITHUB_ROOT / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

EXCLUDE_NATIVE = {"HRP_A_013", "HRP_A_014"}

# ── Plotting setup (Human RNome / Nature guidelines) ──────────────────────────
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.labelsize"] = 8
plt.rcParams["axes.titlesize"] = 9
plt.rcParams["xtick.labelsize"] = 7
plt.rcParams["ytick.labelsize"] = 7
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

# ── Load substitution spectra ────────────────────────────────────────────────
def parse_sample_name(path: Path) -> str:
    """Extract HRP_A_xxx part from filename."""
    name = path.stem.replace("_substitution_spectrum", "")
    parts = name.split("_")
    return "_".join(parts[:3])

all_files = sorted(NATIVE_DIR.glob("*_substitution_spectrum.tsv"))
native_files = []
for f in all_files:
    sample = parse_sample_name(f)
    if any(ex in sample for ex in EXCLUDE_NATIVE):
        continue
    # Strictly exclude IVT samples!
    if "IVT" in f.name:
        continue
    native_files.append(f)

print(f"Native spectrum files explicitly filtered (expected 13): {len(native_files)}")

records = []
for f in native_files:
    sample = parse_sample_name(f)
    df = pd.read_csv(f, sep="\t")
    df = df[["substitution", "rate"]].copy()
    df["sample"] = sample
    df["rate_pct"] = df["rate"] * 100
    records.append(df)

native_df = pd.concat(records, ignore_index=True)

# ── Aggregate per substitution type ──────────────────────────────────────────
sub_order = [f"{r}>{a}" for r in ["A", "C", "G", "T"] for a in ["A", "C", "G", "T"] if r != a]
sub_display = [s.replace("T>", "U>").replace(">T", ">U") for s in sub_order]

stats = []
for sub in sub_order:
    vals = native_df[native_df["substitution"] == sub]["rate_pct"].dropna()
    stats.append({
        "substitution": sub,
        "median": vals.median(),
        "q1": vals.quantile(0.25),
        "q3": vals.quantile(0.75),
        "n": len(vals),
    })
stats = pd.DataFrame(stats)

print("\nMedian substitution rates (%) across 13 native samples:")
print(stats[["substitution", "median"]].to_string(index=False))

# ── Plot ─────────────────────────────────────────────────────────────────────
mm = 1 / 25.4
fig, ax = plt.subplots(figsize=(110 * mm, 70 * mm))

x = np.arange(len(sub_order))
medians = stats["median"].values
errs = np.array([medians - stats["q1"].values,
                 stats["q3"].values - medians])

ref_colours = {"A": "#C0392B", "C": "#2C3E8A", "G": "#27AE60", "T": "#E67E22"}
bar_colours = [ref_colours[s[0]] for s in sub_order]

ax.bar(
    x,
    medians,
    yerr=errs,
    color=bar_colours,
    alpha=0.9,
    width=0.72,
    capsize=2.5,
    edgecolor="white",
    linewidth=0.3,
    error_kw={"linewidth": 0.9, "ecolor": "#333333", "capthick": 0.9},
    zorder=3,
)

np.random.seed(42)
for i, sub in enumerate(sub_order):
    y_vals = native_df[native_df["substitution"] == sub]["rate_pct"].values
    jitter = np.random.uniform(-0.15, 0.15, size=len(y_vals))
    ax.scatter(x[i] + jitter, y_vals, color="#666666", edgecolor="black", linewidth=0.3, s=5, alpha=0.7, zorder=5)

# Reference-base group separators and labels
for sep in [3, 6, 9]:
    ax.axvline(sep - 0.5, color="#cccccc", lw=0.8, zorder=1)

for ref, start_pos in [("A", 0), ("C", 3), ("G", 6), ("T", 9)]:
    disp = "U" if ref == "T" else ref
    ax.text(start_pos + 1, 1.02, disp,
            ha="center", fontsize=10, fontweight="bold",
            color=ref_colours[ref],
            transform=ax.get_xaxis_transform())

ax.set_xticks(x)
ax.set_xticklabels(sub_display, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Substitution rate (%)", fontweight="bold", fontsize=8)
ax.set_xlim(-0.6, len(sub_order) - 0.4)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.8)
ax.spines["bottom"].set_linewidth(0.8)
ax.grid(True, axis="y", alpha=0.25, lw=0.5, zorder=0)

plt.tight_layout()
out_prefix = OUTDIR / "panel_f_substitution_spectrum"
plt.savefig(f"{out_prefix}.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{out_prefix}.png", format="png", bbox_inches="tight", dpi=300)
plt.close()

print(f"\nSaved: {out_prefix}.pdf")
print(f"Saved: {out_prefix}.png")
