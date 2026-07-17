#%%
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np

from matplotlib import font_manager
import matplotlib as mpl
from matplotlib.patches import Patch
from matplotlib.ticker import StrMethodFormatter

from figure_6_helper_functions import extract_fragments, reference_overlap, MOD_COLORS

# Set Arial globally
mpl.rcParams['font.family'] = 'Arial'
fonts = [f.name for f in font_manager.fontManager.ttflist]

mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['svg.fonttype'] = "none"
print("Arial" in fonts)





def plot_18S(overlap_multi_s, overlap_single, non_overlap_multi_s, non_overlap_single, frag_dict):
    # =========================
    # FIGURE LAYOUT (BOTTOM FRAGMENTS, SMALL BAND)
    # =========================

    width_mm = 90
    height_mm = 45

 
    fig_width_in = 12
    fig_height_in = 6

    # factor by which the figure shrinks when placed at width_mm
    scale_factor = (width_mm / 25.4) / fig_width_in   

    target_final_pt = 5.1                              
    base_fontsize = target_final_pt / scale_factor     

    label_fontsize  = base_fontsize          
    tick_fontsize   = base_fontsize * 0.92   
    legend_fontsize = base_fontsize * 0.92   

    fig = plt.figure(figsize=(fig_width_in, fig_height_in))
    gs = gridspec.GridSpec(2, 1, height_ratios=[10, 1], hspace=0.05)

    ax_mod = fig.add_subplot(gs[0])
    ax_frag = fig.add_subplot(gs[1], sharex=ax_mod)

    # =========================
    # 1. MODIFICATION TRACK (LOLLIPOPS FIXED)
    # =========================

    #merge overlap and non overlap to order after freq for vizualisation
    overlap_multi_s_plot = overlap_multi_s
    overlap_multi_s_plot["is_overlap"] = 1
    overlap_multi_s_plot["multi_s"] = 1

    overlap_single_plot = overlap_single
    overlap_single_plot["is_overlap"] = 1
    overlap_single_plot["multi_s"] = 0

    non_overlap_multi_s_plot = non_overlap_multi_s
    non_overlap_multi_s_plot["is_overlap"] = 0
    non_overlap_multi_s_plot["multi_s"] = 1

    non_overlap_single_plot = non_overlap_single
    non_overlap_single_plot["is_overlap"] = 0
    non_overlap_single_plot["multi_s"] = 0

    plot_df = pd.concat([overlap_multi_s_plot, overlap_single_plot,non_overlap_multi_s_plot, non_overlap_single_plot], ignore_index=True)

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

        if row["is_overlap"] == 1 and row["multi_s"] == 1:
            # OVERLAP 
            color = MOD_COLORS[row["name"]]

            stem_z = i * 2
            head_z = stem_z + 1

            ax_mod.vlines(
                row["chromStart"],
                ymin=0,
                ymax=row["frequency"],
                color=color,
                linewidth=4.0,
                alpha=0.9,
                zorder=stem_z
            )

            ax_mod.scatter(
                row["chromStart"],
                row["frequency"],
                marker="s",
                c=color,
                edgecolors="black",
                linewidths=1.5,
                s=110,
                zorder=head_z
            )

            ax_mod.scatter(
                row["chromStart"],
                row["frequency"],
                marker="x",
                c="black",
                s=80,
                linewidths=1.5,
                zorder=head_z + 0.1
            )
        
        elif row["is_overlap"] == 1 and row["multi_s"] == 0:
            # OVERLAP 
            color = MOD_COLORS[row["name"]]

            stem_z = i * 2
            head_z = stem_z + 1

            ax_mod.vlines(
                row["chromStart"],
                ymin=0,
                ymax=row["frequency"],
                color=color,
                linewidth=4.0,
                alpha=0.9,
                zorder=stem_z
            )

            ax_mod.scatter(
                row["chromStart"],
                row["frequency"],
                #marker="s",
                c=color,
                edgecolors="none",
                s=110,
                zorder=head_z
            )

            ax_mod.scatter(
                row["chromStart"],
                row["frequency"],
                marker="x",
                c="black",
                s=65,
                linewidths=1.5,
                zorder=head_z + 0.1
            )
        
        
        elif row["is_overlap"] ==0 and row["multi_s"] == 1:
            # NON-OVERLAP 

            color = MOD_COLORS[row["name"]]

            stem_z = i*2
            head_z = stem_z + 1

            ax_mod.vlines(
                row["chromStart"],
                ymin=0,
                ymax=row["frequency"],
                color=color,
                linewidth=4.0,
                alpha=0.9,
                zorder=stem_z
            )

            ax_mod.scatter(
                row["chromStart"],
                row["frequency"],
                marker = "s",
                c=color,
                edgecolors="black",
                linewidths=1.5,
                s=110,
                zorder=head_z
            )

        elif row["is_overlap"] ==0 and row["multi_s"] == 0:
            # NON-OVERLAP 

            color = MOD_COLORS[row["name"]]

            stem_z = i*2
            head_z = stem_z + 1

            ax_mod.vlines(
                row["chromStart"],
                ymin=0,
                ymax=row["frequency"],
                color=color,
                linewidth=4.0,
                alpha=0.9,
                zorder=stem_z
            )

            ax_mod.scatter(
                row["chromStart"],
                row["frequency"],
                c=color,
                s=110,
                zorder=head_z
            )

    ax_mod.set_ylabel("Modification level (%)", fontsize=label_fontsize)
    ax_mod.grid(axis="y", alpha=0.6)
    ax_mod.set_axisbelow(True)

    ax_mod.spines["top"].set_visible(False)
    ax_mod.spines["right"].set_visible(False)

    ax_mod.tick_params(axis="x", labelbottom=False)

    ax_mod.set_ylim(0, 105)

    # =========================
    # 2. FRAGMENTS TRACK (BOTTOM, COMPACT)
    # =========================
    for (start, end) in frag_dict.keys():
        ax_frag.broken_barh(
            #fragmetns 1 indexed; bed file coordinates 0 indexed; end inclusive
            [(start-1, (end - start) +1)],
            (0.2, 0.6),   
            facecolors="dimgrey",
            alpha=0.99
        )

    ax_frag.set_ylabel("Fragments", rotation=0, fontsize=label_fontsize)
    ax_frag.yaxis.set_label_coords(-0.065, 0.25)
    ax_frag.set_yticks([])

    ax_frag.spines["top"].set_visible(False)
    ax_frag.spines["right"].set_visible(False)
    ax_frag.spines["left"].set_visible(False)

    ax_frag.set_xlabel("18S Position", fontsize=label_fontsize)
    ax_frag.spines["bottom"].set_color("#808080")
    

    # =========================
    # 3. SHARED X LIMITS
    # =========================
    xmin = min(
        overlap_multi_s["chromStart"].min(),
        overlap_single["chromStart"].min(),
        non_overlap_multi_s["chromStart"].min(),
        non_overlap_single["chromStart"].min(),
        min(start for start, _ in frag_dict.keys())
    )

    xmax = max(
        overlap_multi_s["chromStart"].max(),
        overlap_single["chromStart"].max(),
        non_overlap_multi_s["chromStart"].max(),
        non_overlap_single["chromStart"].max(),
        max(end for _, end in frag_dict.keys())
    )

    ax_mod.set_xlim(xmin - 10, xmax + 10)
    
    ax_mod.set_xticks(np.arange(0, xmax + 1, 250))

    ax_mod.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    # =========================
    # LEGEND (COMPOSITE CORRECT)
    # =========================

    # X marker only (no square around it)
    overlap_cross = Line2D([], [], marker='x', color='black', markeredgewidth=1.8, linestyle='None', markersize=7, label="Overlap")

    # Square with black edge (multi-sample indicator)
    overlap_square_edge = Line2D([], [], marker='s', color='none', label="Multiple-samples", markerfacecolor='lightgrey', markeredgecolor='black', markeredgewidth=1.5, linestyle='None', markersize=7)

    no_overlap = Line2D(
        [0], [0],
        marker='o',
        color='none',
        markerfacecolor='lightgrey',
        markeredgecolor='none',
        markersize=8,
        label="Single-sample"
    )

    present_mods = sorted(
        set(overlap_multi_s["name"]).union(overlap_single["name"]).union(non_overlap_multi_s["name"]).union(non_overlap_single["name"]),
        key=lambda m: (m[-1], m)
    )

    mod_handles = [
        Line2D(
            [0], [0],
            marker='o',
            color='none',
            label=mod,
            markerfacecolor=MOD_COLORS[mod],
            markeredgecolor='none',
            markersize=8
        )
        for mod in present_mods
    ]

    blank = Patch(facecolor="none", edgecolor="none", label="")
    ac4C = Patch(facecolor=MOD_COLORS["ac4C"], label="ac4C")
    m227G = Patch(facecolor=MOD_COLORS["m2,2,7G"], label="m2,2,7G")
    mxU = Patch(facecolor=MOD_COLORS["mxU"], label="mxU")
    mxA = Patch(facecolor=MOD_COLORS["mxA"], label="mxA")
    mxC = Patch(facecolor=MOD_COLORS["mxC"], label="mxC")
    mxG = Patch(facecolor=MOD_COLORS["mxG"], label="mxG")

    grid = [[overlap_cross,overlap_square_edge, no_overlap, blank],
            [mxA,   ac4C ,  m227G,  mxU, ],
            [blank, mxC,     mxG,  blank]]
    
    handles = []
    nrows = 3
    ncols = len(grid[0])

    for col in range(ncols):
        for row in range(nrows):
            item = grid[row][col]
            if item is not None:
                handles.append(item)

    legend = ax_mod.legend(
        handles=handles,
        ncol=ncols,
        columnspacing=0.5,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        frameon=True,
        prop={"size": legend_fontsize},
    )   

    for ax in [ax_mod, ax_frag]:
        ax.spines["bottom"].set_color("#808080")
        ax.spines["left"].set_color("#808080")

        ax.tick_params(
            axis="both",
            color="#808080",      # tick marks
            labelcolor="black",   # tick labels (explicitly keep black)
            labelsize=tick_fontsize
        )

    plt.tight_layout()
    plt.savefig("rRNA_18S_modifications.svg", bbox_inches="tight")
    plt.show()
        


def main():
    
    taoka_ref = "/home/alex/Documents/PHD/RNOME/Data-Jamboree/rRNA_session/Reference_Files/H.sapiens_ref_bedRmod_All.bed" #insert your reference 

    file_names = ["HRP_C_001_rRNA_001","HRP_C_002_rRNA_001","HRP_C_005_rRNA_004+HRP_C_006_rRNA_004","HRP_C_014_rRNA_001"]
    base = "/home/alex/Documents/PHD/RNOME/massspec/final_data/new_sunmission_18_06/raw_rRNA/" #select your data storage folder
    end = ".bed"

    rRNA_detailed_consensus = "/home/alex/Documents/PHD/RNOME/massspec/final_data/new_sunmission_18_06/rRNA_consensus.bed"
    
    frag_dict_18S, frag_dict_28S = extract_fragments(file_names,base,end)

    # Find overlaps
    overlap_multi_s_18S, overlap_single_18S, non_overlap_multi_s_18S, non_overlap_single_18S = reference_overlap(rRNA_detailed_consensus,taoka_ref, "18S")
    

    plot_18S(overlap_multi_s_18S, overlap_single_18S, non_overlap_multi_s_18S, non_overlap_single_18S, frag_dict_18S)
 


main()