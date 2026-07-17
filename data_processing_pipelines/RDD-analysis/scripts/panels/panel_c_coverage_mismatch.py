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
COV_DIR = GITHUB_ROOT / "inputs"
OUTDIR = GITHUB_ROOT / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

NATIVE_SAMPLES = [
    "HRP_A_001_native_polyA_RNA_001", "HRP_A_002_native_polyA_RNA_001",
    "HRP_A_003_native_polyA_RNA_001", "HRP_A_005_native_polyA_RNA_001",
    "HRP_A_006_native_polyA_RNA_004", "HRP_A_007_native_polyA_RNA_001",
    "HRP_A_009_native_polyA_RNA_001", "HRP_A_009_native_polyA_RNA_002",
    "HRP_A_010_native_polyA_RNA_001", "HRP_A_010_native_polyA_RNA_002",
    "HRP_A_011_native_polyA_RNA_001", "HRP_A_011_native_polyA_RNA_002",
    "HRP_A_011_native_polyA_RNA_003",
]
IVT_SAMPLES = [
    "HRP_A_001_IVT_polyA_RNA_002", "HRP_A_006_IVT_polyA_RNA_005"
]

BIN_ORDER = ["1–5", "5–10", "10–20", "20–30", "30–50", "50–75", "75–100",
             "100–150", "150–200", "200–300", "300–500", "500–1k", "1k+"]
BIN_MID = {
    "1–5": 3, "5–10": 7.5, "10–20": 15, "20–30": 25, "30–50": 40,
    "50–75": 62.5, "75–100": 87.5, "100–150": 125, "150–200": 175,
    "200–300": 250, "300–500": 400, "500–1k": 750, "1k+": 1500,
}

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.labelsize"] = 8
plt.rcParams["axes.titlesize"] = 8
plt.rcParams["xtick.labelsize"] = 7
plt.rcParams["ytick.labelsize"] = 7
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

def aggregate_samples(samples):
    sample_curves = {}
    for sample in samples:
        files = sorted(COV_DIR.glob(f"{sample}_chr*_coverage_binned_error_rates.tsv"))
        if not files: continue
        
        chrom_records = []
        for f in files:
            df = pd.read_csv(f, sep="\t")
            chrom_records.append(df)
            
        df_all = pd.concat(chrom_records, ignore_index=True)
        
        agg_records = []
        for bin_label in BIN_ORDER:
            sub = df_all[df_all["cov_bin"] == bin_label]
            if sub.empty: continue
            total_bases = sub["total_bases"].sum()
            if total_bases == 0: continue
            weighted_rate = (sub["sub_weighted_rate"] * sub["total_bases"]).sum() / total_bases
            agg_records.append({
                "cov_bin": bin_label,
                "weighted_rate": weighted_rate,
            })
        sample_curves[sample] = pd.DataFrame(agg_records)
    return sample_curves

native_curves = aggregate_samples(NATIVE_SAMPLES)
ivt_curves = aggregate_samples(IVT_SAMPLES)

def get_plot_data(curves, samples):
    median_x, median_y, err_low, err_high = [], [], [], []
    rates_by_bin = {}
    for i, bin_label in enumerate(BIN_ORDER):
        rates = []
        for sample in samples:
            if sample not in curves: continue
            row = curves[sample][curves[sample]["cov_bin"] == bin_label]
            if not row.empty:
                rates.append(row["weighted_rate"].values[0] * 100)
        if rates:
            median_x.append(i)
            med = np.median(rates)
            median_y.append(med)
            q25, q75 = np.percentile(rates, [25, 75])
            err_low.append(med - q25)
            err_high.append(q75 - med)
            rates_by_bin[bin_label] = med
    return median_x, median_y, err_low, err_high, rates_by_bin

nx, ny, n_low, n_high, n_rates = get_plot_data(native_curves, NATIVE_SAMPLES)
ix, iy, i_low, i_high, i_rates = get_plot_data(ivt_curves, IVT_SAMPLES)

n_min_bin = min(n_rates, key=n_rates.get)
n_max_bin = max(n_rates, key=n_rates.get)
i_min_bin = min(i_rates, key=i_rates.get)
i_max_bin = max(i_rates, key=i_rates.get)

print(f"Native Lowest Error Bin: {n_min_bin} ({n_rates[n_min_bin]:.2f}%)")
print(f"Native Highest Error Bin: {n_max_bin} ({n_rates[n_max_bin]:.2f}%)")
print(f"IVT Lowest Error Bin: {i_min_bin} ({i_rates[i_min_bin]:.2f}%)")
print(f"IVT Highest Error Bin: {i_max_bin} ({i_rates[i_max_bin]:.2f}%)")

mm = 1 / 25.4
fig, ax = plt.subplots(figsize=(100 * mm, 85 * mm))

native_color = "#70AB5B"
ivt_color = "#bebcb7"

ax.errorbar(nx, ny, yerr=[n_low, n_high], marker="o", markersize=5, linewidth=2.0, color=native_color, ecolor=native_color, elinewidth=1.2, capsize=3, capthick=1.0, alpha=0.95, label="Native", zorder=4)
ax.errorbar(ix, iy, yerr=[i_low, i_high], marker="o", markersize=5, linewidth=2.0, color=ivt_color, ecolor=ivt_color, elinewidth=1.2, capsize=3, capthick=1.0, alpha=0.95, label="IVT", zorder=4)

ax.set_xlabel("Coverage depth (reads)", fontweight="bold", fontsize=11)
ax.set_ylabel("Mismatch error rate (%)", fontweight="bold", fontsize=11)

ax.set_xticks(range(len(BIN_ORDER)))
ax.set_xticklabels(BIN_ORDER, rotation=45, ha="right", fontsize=6)

ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.2f}"))
ax.tick_params(axis="y", labelsize=7)
ax.set_ylim(bottom=0)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.8)
ax.spines["bottom"].set_linewidth(0.8)

ax.legend(loc="lower center", fontsize=6, frameon=False, ncol=2)

plt.tight_layout()
plt.savefig(OUTDIR / "panel_c_coverage_mismatch.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.savefig(OUTDIR / "panel_c_coverage_mismatch.png", format="png", bbox_inches="tight", dpi=300)
plt.close()

print(f"Saved: {OUTDIR / 'panel_c_coverage_mismatch.png'}")
print(f"Saved: {OUTDIR / 'panel_c_coverage_mismatch.pdf'}")
