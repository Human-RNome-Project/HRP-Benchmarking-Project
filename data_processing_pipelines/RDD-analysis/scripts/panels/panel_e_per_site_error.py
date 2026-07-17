#!/usr/bin/env python3
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore", category=FutureWarning)

GITHUB_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = GITHUB_ROOT / "inputs"
OUT_DIR = GITHUB_ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EXCLUDE_NATIVE = {"HRP_A_013", "HRP_A_014"}

MODIFICATION_COLORS = {
    "Am": "#D44F3E", "i6A": "#C0392B", "m1A": "#A52020", "m6A": "#721817",
    "m6,6A": "#F5C4B8", "m6t6A": "#F0A898", "ms2i6A": "#8B1A1A", "mxA": "#E8907A",
    "t6A": "#E06050", "ac4C": "#6AAED6", "Cm": "#0D3B6E", "f5C": "#B8D9F0",
    "f5Cm": "#D4E4F0", "hm5C": "#2A7FC4", "hm5Cm": "#559CC0", "m3C": "#B8D9F0",
    "m5C": "#001427", "m5Cm": "#3A8FD4", "mxC": "#1E6EB5", "Gm": "#74B354",
    "I": "#2D6E1E", "m1G": "#3A7A28", "m1I": "#BAE09E", "m2,2,7G": "#A8D48A",
    "m2,2G": "#8DC46A", "m2G": "#5C9840", "m7G": "#D4EEC4", "mxG": "#4A8532",
    "Q": "#ACC49B", "acp3D": "#FDE9B8", "acp3U": "#A86A00", "cm5U": "#F6C830",
    "D": "#F9D47E", "m5U": "#F5BE45", "mchm5U": "#E8950A", "mcm5U": "#FBD96A",
    "mcmo5U": "#FCE08A", "mxU": "#F5BE45", "ncm5U": "#FCE08A", "s2U": "#FAD55A",
    "U*": "#E0A858", "Um": "#C47A02", "Y": "#F0A202",
}

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 8
plt.rcParams["axes.labelsize"] = 9
plt.rcParams["axes.titlesize"] = 10
plt.rcParams["xtick.labelsize"] = 9
plt.rcParams["ytick.labelsize"] = 8
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0

def parse_sample_name(path: Path) -> str:
    name = path.stem.replace("_modkit_stratified_rates", "")
    parts = name.split("_")
    return "_".join(parts[:3])

files = sorted(DATA_DIR.glob("*_modkit_stratified_rates.tsv"))
files = [f for f in files if not any(ex in parse_sample_name(f) for ex in EXCLUDE_NATIVE)]

records = []
for f in files:
    sample = parse_sample_name(f)
    df = pd.read_csv(f, sep="\t")
    df["sample_id"] = sample
    df["total_mismatch_rate_pct"] = df["total_mismatch_rate"] * 100
    records.append(df)

df = pd.concat(records, ignore_index=True)
df = df.rename(columns={"label": "category"})

unmod_vals = df[df["category"] == "unmodified_sites"]["total_mismatch_rate_pct"].values

stats_records = []
for cat, grp in df.groupby("category"):
    vals = grp["total_mismatch_rate_pct"].values
    med = np.median(vals)
    q1 = np.percentile(vals, 25)
    q3 = np.percentile(vals, 75)
    
    if cat == "unmodified_sites":
        pval = 1.0
    else:
        # Mann-Whitney U test compared to unmodified
        if len(vals) > 0 and len(unmod_vals) > 0:
            stat, pval = mannwhitneyu(vals, unmod_vals, alternative='two-sided')
        else:
            pval = 1.0
            
    stats_records.append({
        "category": cat,
        "median_rate": med,
        "q1": q1,
        "q3": q3,
        "pval": pval
    })

stats = pd.DataFrame(stats_records)

unmod = stats[stats["category"] == "unmodified_sites"].copy()
mods = stats[stats["category"].str.startswith("mod_")].copy()
mods["mod_name"] = mods["category"].str.replace("mod_", "", regex=False)
mods = mods.sort_values("median_rate", ascending=True)

unmod_median = unmod["median_rate"].iloc[0]
unmod_q1 = unmod["q1"].iloc[0]
unmod_q3 = unmod["q3"].iloc[0]

mods["fold_change"] = mods["median_rate"] / unmod_median

plot_df = pd.concat([
    pd.DataFrame({
        "category": ["unmodified_sites"],
        "mod_name": ["Unmodified"],
        "median_rate": [unmod_median],
        "q1": [unmod_q1],
        "q3": [unmod_q3],
        "pval": [1.0],
        "fold_change": [1.0],
    }),
    mods[["category", "mod_name", "median_rate", "q1", "q3", "pval", "fold_change"]],
], ignore_index=True)

def get_sig_stars(p):
    if pd.isna(p): return "ns"
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"

mm = 1 / 25.4
fig, (ax_top, ax_bottom) = plt.subplots(
    2, 1,
    sharex=True,
    figsize=(120 * mm, 100 * mm),
    gridspec_kw={"height_ratios": [3.5, 1]},
)

x = np.arange(len(plot_df))
medians = plot_df["median_rate"].values
errs = np.array([medians - plot_df["q1"].values, plot_df["q3"].values - medians])

colors = ["#7F8C8D"]
for mod in plot_df["mod_name"].iloc[1:]:
    colors.append(MODIFICATION_COLORS.get(mod, "#2C3E50"))

unmod_mask = plot_df["mod_name"] == "Unmodified"

np.random.seed(42)

for ax in (ax_top, ax_bottom):
    ax.bar(x, medians, color=colors, edgecolor="black", linewidth=0.8, zorder=3)
    
    for i, cat in enumerate(plot_df["category"]):
        y_vals = df[df["category"] == cat]["total_mismatch_rate_pct"].values
        jitter = np.random.uniform(-0.15, 0.15, size=len(y_vals))
        ax.scatter(x[i] + jitter, y_vals, color="#666666", edgecolor="black", linewidth=0.3, s=5, alpha=0.7, zorder=5)

ax_bottom.errorbar(
    x[unmod_mask], medians[unmod_mask], yerr=errs[:, unmod_mask],
    fmt="none", ecolor="#333333", capsize=3, capthick=1.0, elinewidth=1.0, zorder=4
)
ax_top.errorbar(
    x[~unmod_mask], medians[~unmod_mask], yerr=errs[:, ~unmod_mask],
    fmt="none", ecolor="#333333", capsize=3, capthick=1.0, elinewidth=1.0, zorder=4
)

for i, (fold, pval) in enumerate(zip(plot_df["fold_change"], plot_df["pval"])):
    if i == 0: continue
    height = medians[i]
    err = errs[1, i]
    stars = get_sig_stars(pval)
    text_str = f"{fold:.0f}x\n{stars}"
    ax_top.text(
        x[i], height + err + 0.1, text_str,
        ha="center", va="bottom", fontsize=7, fontweight="bold", color="black"
    )

ax_top.set_ylim(0.7, 4.8)
ax_top.spines["top"].set_visible(False)
ax_top.spines["right"].set_visible(False)
ax_top.spines["bottom"].set_visible(False)
ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
ax_top.set_yticks([0.7, 1, 2, 3, 4])
ax_top.set_yticklabels(["0.7", "1", "2", "3", "4"])
ax_top.grid(True, axis="y", alpha=0.25, lw=0.5, zorder=0)

ax_bottom.set_ylim(0, 0.1)
ax_bottom.spines["top"].set_visible(False)
ax_bottom.spines["right"].set_visible(False)
ax_bottom.set_yticks([0, 0.05, 0.1])
ax_bottom.grid(True, axis="y", alpha=0.25, lw=0.5, zorder=0)

ax_bottom.set_xticks(x)
ax_bottom.set_xticklabels(plot_df["mod_name"], rotation=45, ha="right")
ax_bottom.set_xlabel("Modification type", fontweight="bold")

d = 0.02
ax_top.plot((-d, d), (-d, d), transform=ax_top.transAxes, color="black", linewidth=1.0, clip_on=False)
ax_bottom.plot((-d, d), (1 - d, 1 + d), transform=ax_bottom.transAxes, color="black", linewidth=1.0, clip_on=False)

fig.text(0.02, 0.5, "Per-site mismatch error rate (%)", va="center", ha="center", rotation="vertical", fontweight="bold", fontsize=9)

plt.tight_layout(rect=[0.06, 0.0, 1.0, 1.0])
out_prefix = OUT_DIR / "panel_e_per_site_error"
plt.savefig(f"{out_prefix}.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{out_prefix}.png", format="png", bbox_inches="tight", dpi=300)
plt.close()

print(f"Saved: {out_prefix}.png")
