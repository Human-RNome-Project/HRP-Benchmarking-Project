#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib import font_manager

# Set Arial globally
mpl.rcParams['font.family'] = 'Arial'
fonts = [f.name for f in font_manager.fontManager.ttflist]

mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['svg.fonttype'] = "none"


CAP_COLORS = {
    "mγGTP": "#5E4FA2",
    "m7GpppAm":        "#7B6FD0",
    "m2,2,7GpppAm":    "#A55CC5",
    "AppG":            "#C17BCF",
    "AppC":            "#D7A2DD",
    "AppppA":          "#E6C5EB",
    "NAD":             "#8A6BBE",
    "UDPGlc":          "#B38DD9",
    "UDPGlcNAc":       "#CDB4E8",
}

#load data file
df = pd.read_excel("HRP_C_009_cap_structures.xlsx", skiprows=10, engine="openpyxl")


# drop all Unnamed columns
df = df.loc[:, ~df.columns.astype(str).str.contains("^Unnamed")]
df.columns = df.columns.str.strip().str.replace(":", "", regex=False)
keep_cols = ["Modification", "Amount fmol / 1k rN"]
df = df[keep_cols]
df["Amount fmol / 1k rN"] = df["Amount fmol / 1k rN"]*6.022*(10**8)
df["Modification"] = df["Modification"].replace("gammamethyl GTP","mγGTP")
df["Modification"] = df["Modification"].replace("m227GpppAm","m2,2,7GpppAm")
print(df)

df = df.sort_values("Amount fmol / 1k rN", ascending=False)

break_point = 0.00015

width_mm = 90
height_mm = 50

colors = [CAP_COLORS.get(mod, "#808080") for mod in df["Modification"]]

# lower panel = 1/3 of visible plot area -> height ratio upper:lower = 2:1
fig, (ax_top, ax_bot) = plt.subplots(
    nrows=2,
    ncols=1,
    sharex=True,
    figsize=(width_mm / 25.4, height_mm / 25.4),
    gridspec_kw={"height_ratios": [2, 1], "hspace": 0.05},
)

x = np.arange(len(df))

for ax in (ax_top, ax_bot):
    bars = ax.bar(
        x,
        df["Amount fmol / 1k rN"],
        color=colors,
        width=0.75,
        zorder=2,
    )

# value labels — place in whichever panel the label height actually falls in
for i, val in enumerate(df["Amount fmol / 1k rN"]):
    label_y = val * 1.1
    target_ax = ax_top if label_y >= break_point else ax_bot
    target_ax.text(
        x[i],
        label_y,
        f"{val:.1e}",
        ha="center",
        va="bottom",
        fontsize=5,
        zorder=3,
    )

# y-limits: split at break_point
top_max = df["Amount fmol / 1k rN"].max()
ax_top.set_ylim(break_point, top_max * 1.2)
ax_bot.set_ylim(0, break_point)

# full border on both panels
for ax in (ax_top, ax_bot):
    for spine in ax.spines.values():
        spine.set_visible(True)

# hide only the spines between the two panels (where the break is)
ax_top.spines["bottom"].set_visible(False)
ax_bot.spines["top"].set_visible(False)
ax_top.tick_params(labeltop=False, bottom=False)
ax_bot.xaxis.tick_bottom()

# diagonal break marks
d = 0.5  # size of diagonal lines in points
kwargs = dict(
    marker=[(-1, -d), (1, d)],
    markersize=6,
    linestyle="none",
    color="k",
    mec="k",
    mew=1,
    clip_on=False,
)
ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **kwargs)
ax_bot.plot([0, 1], [1, 1], transform=ax_bot.transAxes, **kwargs)

# Nature-style aesthetics
ax_top.tick_params(direction="out", length=3, width=0.8)
ax_bot.tick_params(direction="out", length=3, width=0.8)
ax_top.tick_params(axis="y", labelsize=5.2)
ax_bot.tick_params(axis="y", labelsize=5.2)

ax_bot.set_xticks(x)
ax_bot.set_xticklabels(df["Modification"], rotation=45, ha="right", fontsize=5.2)
ax_bot.set_xlabel("CAP structures", fontsize=6)

# shared y-label centered across both panels
fig.text(
    0.0, 0.5, "Abundance per 1,000 nt",
    va="center", rotation="vertical", fontsize=6,
)

 
plt.tight_layout()
plt.savefig("cap_modifications_barplot.svg", bbox_inches="tight")
plt.show()


