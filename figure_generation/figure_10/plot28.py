import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import StrMethodFormatter



def plot_28(overlap, non_overlap):

    mod_colors= {
    # A — m6A family (deep crimson → light blush)
    "Am":      "#D44F3E",
    "I":       "#98343E",
    "i6A":      "#C0392B",
    "m1A":     "#A52020",
    "m6A":     "#721817",
    "m6,6A":   "#F5C4B8",
    "m6t6A":    "#F0A898",
    "ms2i6A":   "#8B1A1A",
    "mxA":     "#E8907A",
    "t6A":      "#E06050", 

    # C — m5C family (near-black navy → pale sky)
    "ac4C":    "#6AAED6",
    "Cm":      "#0D3B6E",
    "f5C":      "#B8D9F0",
    "f5Cm":     "#D4E4F0",
    "hm5C":     "#2A7FC4",
    "hm5Cm":    "#559CC0",
    "m3C":     "#B8D9F0",
    "m5C":     "#001427",
    "m5Cm":     "#3A8FD4",
    "mxC":     "#1E6EB5",

    # G — Inosine family (dark forest → pale mint)
    "Gm":      "#74B354",
    "m1G":     "#3A7A28",
    "m1I":     "#BAE09E",
    "m2,2,7G": "#A8D48A",
    "m2,2G":    "#8DC46A",
    "m2G":     "#5C9840",
    "m7G":     "#D4EEC4",
    "mxG":     "#4A8532",
    "Q": "#ACC49B",

    # U — Psi family (deep amber → pale gold)
    "acp3D":    "#FDE9B8",
    "acp3U":    "#A86A00",
    "cm5U":     "#F6C830",
    "D": "#F9D47E",
    "m5U":          "#FBDD80",
    "mchm5U":   "#E8950A",
    "mcm5U":    "#FAD55A",
    "mcmo5U":   "#FBD96A",
    "mxU":     "#F5BE45",
    "ncm5U":    "#FCE08A",
    "s2U":          "#E0A858",
    "U*":           "#FBE8AE",
    "Um":      "#C47A02",
    "Y":       "#F0A202"}


    # =========================
    # FIGURE LAYOUT (BOTTOM FRAGMENTS, SMALL BAND)
    # =========================

    width_mm = 183
    height_mm = 35

    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
   

    # =========================
    # 1. MODIFICATION TRACK (LOLLIPOPS FIXED)
    # =========================

    #merge overlap and non overlap to order after freq for vizualisation
    overlap_plot = overlap
    overlap_plot["is_overlap"] = 1

    non_overlap_plot = non_overlap
    non_overlap_plot["is_overlap"] = 0

    plot_df = pd.concat([overlap_plot, non_overlap_plot], ignore_index=True)

    print(len(plot_df))
    # order fpr plotting
    plot_df = plot_df.sort_values("frequency", ascending=False).reset_index(drop=True)

    chrom = plot_df["chromStart"].values
    freq = plot_df["frequency"].values
    is_overlap = plot_df["is_overlap"].values
    #multi_samples = plot_df["multi_s"].values

    # initial rank
    order = np.arange(len(plot_df), dtype=float)

    # -----------------------------
    # LOCAL CORRECTION RULE
    # -----------------------------
    for i in range(len(plot_df)):

        if is_overlap[i] == 1:

            for j in range(len(plot_df)):

                if is_overlap[j] == 0:

                    # condition 1: spatial proximity
                    if abs(chrom[i] - chrom[j]) <= 100:

                        # condition 2: frequency proximity
                        if abs(freq[i] - freq[j]) <= 5:

                            # enforce overlap ABOVE non-overlap
                            if order[i] < order[j]:
                                order[i] = order[j] + 0.1

    plot_df["order"] = order

    # -----------------------------
    # FINAL ORDER
    # -----------------------------
    plot_df = plot_df.sort_values(
        ["order", "frequency"],
        ascending=[True, False]
    ).reset_index(drop=True)


    
    for i, (_, row) in enumerate(plot_df.iterrows()):

        if row["is_overlap"] == 1:
            # OVERLAP 
            color = mod_colors[row["name"]]

            stem_z = i * 2
            head_z = stem_z + 1

            ax.vlines(
                row["chromStart"],
                ymin=0,
                ymax=row["frequency"],
                color=color,
                linewidth=1.5,
                alpha=0.9,
                zorder=stem_z
            )

            ax.scatter(
                row["chromStart"],
                row["frequency"],
                #marker="s",
                c=color,
                edgecolors="none",
                s=30,
                zorder=head_z
            )

            ax.scatter(
                row["chromStart"],
                row["frequency"],
                marker="x",
                c="black",
                s=10,
                linewidths=1,
                zorder=head_z + 0.1
            )
        

        elif row["is_overlap"] ==0:
            # NON-OVERLAP 

            color = mod_colors[row["name"]]

            stem_z = i*2
            head_z = stem_z + 1

            ax.vlines(
                row["chromStart"],
                ymin=0,
                ymax=row["frequency"],
                color=color,
                linewidth=1.5,
                alpha=0.9,
                zorder=stem_z
            )

            ax.scatter(
                row["chromStart"],
                row["frequency"],
                c=color,
                s=20,
                zorder=head_z
            )

    ax.set_ylabel("Modification level (%)", fontsize=6)
    ax.grid(axis="y", alpha=0.6)
    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    #ax_mod.tick_params(axis="x", labelbottom=False)

    ax.set_ylim(0, 105)
    ax.tick_params(axis='both', labelsize=5.5)
    ax.set_xlabel("28S Position", fontsize=6)



    # =========================
    # 3. SHARED X LIMITS
    # =========================
    xmin = min(
        overlap["chromStart"].min(),
        non_overlap["chromStart"].min()
    )

    xmax = max(
        overlap["chromStart"].max(),
        non_overlap["chromStart"].max()
    )

    ax.set_xlim(xmin - 10, xmax + 50)
    
    ax.set_xticks(np.arange(0, xmax + 1, 250))

    # =========================
    # LEGEND (COMPOSITE CORRECT)
    # =========================

    # X marker only (no square around it)
    overlap_cross = Line2D([], [], marker='x', color='black', markeredgewidth=1.8, linestyle='None', markersize=5)

 

    present_mods = sorted(
        set(overlap["name"]).union(non_overlap["name"]),
        key=lambda m: (m[-1], m)
    )

    mods_in_legend = [
        mod_name
        for mod_name in mod_colors
        if mod_name in present_mods
    ]

    mod_handles = [
        Line2D(
            [0], [0],
            marker='o',
            color='none',
            label=mod_name,
            markerfacecolor=mod_colors[mod_name],
            markeredgecolor='none',
            markersize=5,
        )
        for mod_name in mods_in_legend
    ]

    legend_elements = [overlap_cross] + mod_handles

    legend = ax.legend(
        legend_elements,
        ["Overlap"] + mods_in_legend,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True,
        fontsize=5.5,
    )

    legend.get_frame().set_edgecolor("#808080")
    legend.get_frame().set_linewidth(1.0)
    legend.get_frame().set_alpha(1)

    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    for ax in [ax]:
        ax.spines["bottom"].set_color("#808080")
        ax.spines["left"].set_color("#808080")

        ax.tick_params(
            axis="both",
            color="#808080",      # tick marks
            labelcolor="black"    # tick labels (explicitly keep black)
        )

    plt.tight_layout()
    plt.savefig("rRNA_28S_SRS.svg", bbox_inches="tight")
    plt.show()