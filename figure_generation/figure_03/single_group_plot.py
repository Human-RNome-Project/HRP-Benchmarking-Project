import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib import font_manager
from translation_dicts import MOD_COLORS

# Set Arial globally
mpl.rcParams['font.family'] = 'Arial'
fonts = [f.name for f in font_manager.fontManager.ttflist]

mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['svg.fonttype'] = "none"


def plot_5S_and_5_8S(df: pd.DataFrame,mean_df:pd.DataFrame, break_point:float, out_name:str, ratio:tuple, type:str):

    mean_col = df.columns[df.columns.str.contains("mean", case=False, na=False)][0]
    df = df.rename(columns={mean_col: "mean"})
    mean_col = df.columns[df.columns.str.contains("std", case=False, na=False)][0]
    df = df.rename(columns={mean_col: "std"})
    canonical = ["A", "C", "G", "U"]
    df = df[~df["name"].isin(canonical)]
    valid_names = df.loc[df["mean"].notna(), "name"]
    df = df[df["name"].isin(valid_names)]
    
   
    mean_df= mean_df[["name","mean","std"]]
    df = df.merge(mean_df[['name', 'mean', 'std']], on='name', how='left', suffixes=('', '_new'))
    df[['mean', 'std']] = df[['mean_new', 'std_new']]
    df = df.drop(columns=['mean_new', 'std_new'])
    df = df.dropna()
    
    width_mm = ratio[0]
    height_mm = ratio[1]

    order = [m for m in MOD_COLORS.keys() if m in df["name"].values]
    reps = [c for c in df.columns if c not in ["name", "std", "mean"]]

    x = np.arange(len(order))
    rng = np.random.default_rng(0)  # fixed seed so jitter is reproducible

    # lower panel = 1/3 of visible plot area -> height ratio upper:lower = 2:1
    fig, (ax_top, ax_bot) = plt.subplots(
        nrows=2,
        ncols=1,
        sharex=True,
        figsize=(width_mm / 25.4, height_mm / 25.4),
        gridspec_kw={"height_ratios": [2, 1], "hspace": 0.05},
    )

    x = np.arange(len(order))
    rng = np.random.default_rng(0)

    # generate jitter once, so every row uses the same relative offsets
    n_reps = len(reps)
    jitter = rng.uniform(-0.25, 0.25, size=n_reps)

    for ax in (ax_top, ax_bot):
        for i, name in enumerate(order):
            row = df[df["name"] == name].iloc[0]
            color = MOD_COLORS[name]

            # bar for the mean, low alpha
            ax.bar(
                x[i],
                row["mean"],
                color=color,
                alpha=0.4,
                width=0.9,
                zorder=1,
            )

            # std as errorbar, drawn on top of the bar
            ax.errorbar(
                x[i],
                row["mean"],
                yerr=row["std"],
                fmt="none",
                ecolor="black",
                elinewidth=1,
                capsize=3,
                zorder=3,
            )

            # individual replicate values as points, jittered along x
            vals = row[reps].values.astype(float)
            mask = ~np.isnan(vals)
            ax.scatter(
                x[i] + jitter[mask],
                vals[mask],
                color=color,
                alpha=1.0,
                s=12,
                edgecolor="black",
                linewidth=0.3,
                zorder=2,
            )

    # y-limits: split at break_point
    top_max = max(
        df["mean"].max() + df["std"].max(),
        df[reps].max(numeric_only=True).max(),
    )
    ax_top.set_ylim(break_point, top_max * 1.05)
    ax_bot.set_ylim(0, break_point)

    # hide the spines between the two panels
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

    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(order, rotation=45, ha="right", fontsize=5.5)
    ax_top.tick_params(axis="y", labelsize=5.5)
    ax_bot.tick_params(axis="y", labelsize=5.5)
    ax_bot.set_xlabel(f"Modifications {type}", fontsize=6)

    # shared y-label centered across both panels
    fig.text(
        0.0, 0.5, "Abundance per 1,000 nt",
        va="center", rotation="vertical", fontsize=6,
    )

    plt.tight_layout()
    plt.savefig(out_name, bbox_inches="tight")
    plt.show()