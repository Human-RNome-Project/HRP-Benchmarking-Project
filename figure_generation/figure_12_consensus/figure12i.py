#%%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import re
from matplotlib.colors import to_rgb
import matplotlib.patches as patches
from matplotlib import font_manager
import matplotlib as mpl
from matplotlib.patches import Rectangle

# Set Arial globally
mpl.rcParams['font.family'] = 'Arial'
fonts = [f.name for f in font_manager.fontManager.ttflist]

mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['svg.fonttype'] = "none"
print("Arial" in fonts)



def build_sprinzl_mapping(df:pd.DataFrame, id_col="tRNA_name"):
    """
    Returns:
        dict[tRNA_name] -> dict[sequential_index] -> sprinzl_coordinate
    """

    df = df.drop(columns=["tRNAs"])

    # columns that represent Sprinzl positions (everything except metadata)
    sprinzl_cols = [c for c in df.columns if c not in [id_col, "tRNAs"]]
    result = {}

    for _, row in df.iterrows():
        tRNA = row[id_col]

        mapping = {}
        idx = 0  

        for col in sprinzl_cols:
            val = row[col]

            if pd.isna(val):
                continue

            mapping[idx] = col  
            idx += 1

        result[tRNA] = mapping

    return result



def transform_tRNA_name(df: pd.DataFrame):
    # tRNA
    name_dict = {'hs_tRNAAla_CGC': 'hs_tRNAAla_CGC',
    'hs_tRNAAla_TGC': 'hs_tRNAAla_TGC',
    'hs_tRNAAla_AGC': 'hs_tRNAAla_AGC',
    'hs_tRNAArg_ACG': 'hs_tRNAArg_ACG',
    'hs_tRNAArg_TCG': 'hs_tRNAArg_TCG',
    'hs_tRNAArg_CCG2': 'hs_tRNAArg_CCG2',
    'hs_tRNAArg_TCT2': 'hs_tRNAArg_TCT2',
    'hs_tRNAArg_CCT': 'hs_tRNAArg_CCT',
    'hs_tRNAArg_CCG1': 'hs_tRNAArg_CCG1',
    'hs_tRNAArg_TCT1': 'hs_tRNAArg_TCT1',
    'hs_tRNAAsn_GTT': 'hs_tRNAAsn_GTT',
    'hs_tRNAAsp_GTC': 'hs_tRNAAsp_GTC',
    'hs_tRNACys_GCA': 'hs_tRNACys_GCA',
    'hs_tRNAGln_CTG_TTG': 'hs_tRNAGln_CTG_TTG',
    'hs_tRNAGlu_CTC': 'hs_tRNAGlu_CTC',
    'hs_tRNAGlu_TTC': 'hs_tRNAGlu_TTC',
    'hs_tRNAGly_GCC': 'hs_tRNAGly_GCC',
    'hs_tRNAGly_TCC': 'hs_tRNAGly_TCC',
    'hs_tRNAGly_CCC': 'hs_tRNAGly_CCC',
    'hs_tRNAHis_GTG': 'hs_tRNAHis_GTG',
    'hs_tRNAIle_AAT': 'hs_tRNAIle_AAT',
    'hs_tRNAIle_TAT': 'hs_tRNAIle_TAT',
    'hs_tRNAIle_GAT': 'hs_tRNAIle_GAT',
    'hs_tRNALeu_AAG': 'hs_tRNALeu_AAG',
    'hs_tRNALeu_CAG': 'hs_tRNALeu_CAG',
    'hs_tRNALeu_CAA': 'hs_tRNALeu_CAA',
    'hs_tRNALeu_TAA': 'hs_tRNALeu_TAA',
    'hs_tRNALys_CTT': 'hs_tRNALys_CTT',
    'hs_tRNALys_TTT': 'hs_tRNALys_TTT',
    'hs_tRNAMet_CAT': 'hs_tRNAMet_CAT',
    'hs_tRNAPhe_GAA': 'hs_tRNAPhe_GAA',
    'hs_tRNAPro_AGG_CGG_TGG': 'hs_tRNAPro_AGG_CGG_TGG',
    'hs_tRNASeC_TCA': 'hs_tRNASeC_TCA',
    'hs_tRNASer_AGA': 'hs_tRNASer_AGA',
    'hs_tRNASer_CGA': 'hs_tRNASer_CGA',
    'hs_tRNASer_GCT': 'hs_tRNASer_GCT',
    'hs_tRNAThr_AGT': 'hs_tRNAThr_AGT',
    'hs_tRNAThr_CGT1': 'hs_tRNAThr_CGT1',
    'hs_tRNAThr_TGT': 'hs_tRNAThr_TGT',
    'hs_tRNAThr_CGT2': 'hs_tRNAThr_CGT2',
    'hs_tRNATrp_CCA': 'hs_tRNATrp_CCA',
    'hs_tRNATyr_GTA2': 'hs_tRNATyr_GTA2',
    'hs_tRNATyr_ATA': 'hs_tRNATyr_ATA',
    'hs_tRNATyr_GTA1': 'hs_tRNATyr_GTA1',
    'hs_tRNAVal_AAC_CAC': 'hs_tRNAVal_AAC_CAC',
    'hs_tRNAVal_TAC': 'hs_tRNAVal_TAC',
    'hs_tRNAArg_CCG': 'hs_tRNAArg_CCG1',
    'hs_tRNAArg_TCT': 'hs_tRNAArg_TCT1',
    'hs-tRNAArg-TCG': 'hs_tRNAArg_TCG',
    'hs_tRNATyr_GTA' : 'hs_tRNATyr_GTA1',
    'hs_tRNAiMet_CAT': 'hs_tRNAiMet_CAT',
    'hs_tRNALeu_TAG': 'hs_tRNALeu_TAG', #Not in hg38 combined ref
    'hs_tRNASer_TGA': 'hs_tRNASer_TGA', #Not in hg38 combined ref
    'hs_tRNAThr_CGT': 'hs_tRNAThr_CGT1',
    'hs_mttRNAAla_TGC': 'hs_mttRNAAla_TGC',
    'hs_mttRNAArg_TCG': 'hs_mttRNAArg_TCG',
    'hs_mttRNAAsn_GTT': 'hs_mttRNAAsn_GTT',
    'hs_mttRNAAsp_GTC': 'hs_mttRNAAsp_GTC',
    'hs_mttRNACys_GCA': 'hs_mttRNACys_GCA',
    'hs_mttRNAGln_TTG': 'hs_mttRNAGln_TTG',
    'hs_mttRNAGlu_TTC': 'hs_mttRNAGlu_TTC',
    'hs_mttRNAGly_TCC': 'hs_mttRNAGly_TCC',
    'hs_mttRNAHis_GTG': 'hs_mttRNAHis_GTG',
    'hs_mttRNAIle_GAT': 'hs_mttRNAIle_GAT',
    'hs_mttRNALeu_TAA': 'hs_mttRNALeu_TAA',
    'hs_mttRNALeu_TAG': 'hs_mttRNALeu_TAG',
    'hs_mttRNALys_TTT': 'hs_mttRNALys_TTT',
    'hs_mttRNAMet_CAT': 'hs_mttRNAMet_CAT',
    'hs_mttRNAPhe_GAA': 'hs_mttRNAPhe_GAA',
    'hs_mttRNAPro_TGG': 'hs_mttRNAPro_TGG',
    'hs_mttRNASer_GCT': 'hs_mttRNASer_GCT',
    'hs_mttRNASer_TGA': 'hs_mttRNASer_TGA',
    'hs_mttRNAThr_TGT': 'hs_mttRNAThr_TGT',
    'hs_mttRNATrp_TCA': 'hs_mttRNATrp_TCA',
    'hs_mttRNATyr_GTA': 'hs_mttRNATyr_GTA',
    'hs_mttRNAVal_TAC': 'hs_mttRNAVal_TAC',
    'hs_tRNAVal_AAC_CAC_G34=I_introduced': 'hs_tRNAVal_AAC_CAC',
    'hs_tRNASer_AGA_G34=I_introduced': 'hs_tRNASer_AGA',
    'hs_tRNAPro_AGG_CGG_TGG_G34=I_introduced': 'hs_tRNAPro_AGG_CGG_TGG',
    'hs_tRNALeu_AAG_G34=I_introduced': 'hs_tRNALeu_AAG',
    'hs_tRNAArg_ACG_G34=I_introduced': 'hs_tRNAArg_ACG',
    'hs_tRNAAla_AGC_G34=I_introduced': 'hs_tRNAAla_AGC',
    'hs_tRNAVal_AAC_CAC G34=I_introduced': 'hs_tRNAVal_AAC_CAC',
    'hs_tRNASer_AGA G34=I_introduced': 'hs_tRNASer_AGA',
    'hs_tRNAPro_AGG_CGG_TGG G34=I_introduced': 'hs_tRNAPro_AGG_CGG_TGG',
    'hs_tRNALeu_AAG G34=I_introduced': 'hs_tRNALeu_AAG',
    'hs_tRNAArg_ACG G34=I_introduced': 'hs_tRNAArg_ACG',
    'hs_tRNAAla_AGC G34=I_introduced': 'hs_tRNAAla_AGC',

    'hs_tRNAArg_ACG G34=I':'hs_tRNAArg_ACG',
    'hs_tRNAAla_AGC G34=I':'hs_tRNAAla_AGC',
    'hs_tRNASer_AGA G34=I introduced':'hs_tRNASer_AGA',

    'tRNA-Pro-AGG-1-1':'hs_tRNAPro_AGG_CGG_TGG',
    'tRNA-VAl-AAC-1-1':'hs_tRNAVal_AAC_CAC',
    }
    try:
        df["tRNA_name"] = df["tRNA_name"].map(name_dict).fillna(df["tRNA_name"])
    except:
        df["chrom"] = df["chrom"].map(name_dict).fillna(df["chrom"])

    return df



def create_data_matrix(positions:list,df_cons:pd.DataFrame,sprinzl_dict: dict,tRNAs:list,sprinzl_df:pd.DataFrame)->pd.DataFrame:
    # map genomic position -> Sprinzl coordinate
    
    df_cons["sprinzl_pos"] = df_cons.apply(
        lambda r: sprinzl_dict.get(r["chrom"], {}).get(r["chromStart"]),
        axis=1
    )

    #empty matrix

    data_matrix = pd.DataFrame(
        [[(0) for _ in range(len(positions))] for _ in range(len(tRNAs))],
        index=tRNAs,
        columns=positions
    )

    data_matrix = data_matrix.astype(object)
    #fill mods
    for _, row in df_cons.iterrows():
        
        mod = row["name"]
        data_matrix.loc[row["chrom"],str(row["sprinzl_pos"])] = mod



    # fill no sequence
    for _, row in sprinzl_df.iterrows():
        trna = row["tRNA_name"]

        for col in sprinzl_df.columns:
            if col == "tRNA_name":
                continue

            if pd.isna(row[col]):
                data_matrix.loc[trna, str(col)] = "ns"

    return data_matrix

def cell_to_color(cell_value, mod_colors):
    if cell_value == 0 or cell_value == "ns":  
        return "#FFFFFF"  
    else: 
        return mod_colors[cell_value] 

def plot_heatmap(data_matrix: pd.DataFrame,present_mods:list):
    print(present_mods)
    data_matrix.index = data_matrix.index.str.removeprefix('hs_tRNA')
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


    colour_matrix = np.empty((data_matrix.shape[0], data_matrix.shape[1], 3), dtype=float)

    for i in range(data_matrix.shape[0]):
        for j in range(data_matrix.shape[1]):
            color = cell_to_color(data_matrix.iloc[i, j],mod_colors)
            rgb_color = to_rgb(color)
            colour_matrix[i, j] = rgb_color

    width_mm = 183
    height_mm = 100

    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))

    ax.imshow(colour_matrix, aspect='auto', interpolation='none', resample=False,zorder=1)

    no_seq_mask = np.vectorize(lambda x: x == "ns")(data_matrix.values)

    for i in range(no_seq_mask.shape[0]):
        for j in range(no_seq_mask.shape[1]):
            if no_seq_mask[i, j]:
                ax.plot(
                    j, i,
                    marker="o",
                    markersize=2,
                    color="black",
                    linestyle="None",
                    markeredgewidth=0,
                    antialiased=True
                )


    #compartment header
    compartments = {
        "AS_5": (1, 7),
        "D_arm": (10, 25),
        "AC": (27, 43),
        "Variable": (44, 48),
        "TYC_arm": (49, 65),
        "AS_3": (66, 72),
        "CCA": (74, 76),
    }

    compartments_names = {
        "AS_5": "AS",
        "D_arm": "D-stem-loop",
        "AC": "AC-stem-loop",
        "Variable": "Variable region",
        "TYC_arm": "TYC-stem-loop",
        "AS_3": "AS",
        "CCA": "CCA",
    }

    # map label -> x index
    pos_to_x = {str(pos): i for i, pos in enumerate(data_matrix.columns)}

    def get_x(pos):
        """Return x-index if exists, else None."""
        return pos_to_x.get(str(pos), None)
    
    for name, (start, end) in compartments.items():
        x_start = get_x(start)
        x_end = get_x(end)

        if x_start is None or x_end is None:
            continue

        # vertical boundary lines
        ax.axvline(x_start - 0.5, color="grey", linewidth=0.8, alpha=0.8)
        ax.axvline(x_end + 0.5, color="grey", linewidth=0.8, alpha=0.8)

    header_ax = ax.inset_axes([0, 1.02, 1, 0.12])  # x, y, w, h (relative to heatmap)
    header_ax.set_xlim(ax.get_xlim())
    header_ax.set_ylim(0, 1)
    header_ax.axis("off")
    header_ax.set_clip_on(False)
    for spine in header_ax.spines.values():
        spine.set_visible(False)

    for name, (start, end) in compartments.items():
        x_start = get_x(start)
        x_end = get_x(end)

        if x_start is None or x_end is None:
            continue

        width = (x_end - x_start)

        # grey header block
        header_ax.add_patch(
            Rectangle(
                (x_start - 0.5, 0.05),
                width+1,
                0.3,
                facecolor="lightgrey",
                edgecolor="darkgrey",
                alpha=0.35,
                zorder = 5
            )
        )

        # label
        header_ax.text(
            x_start + width / 2,
            0.17,
            compartments_names[name],
            ha="center",
            va="center",
            fontsize=5.5,
            color="black",
            zorder = 6
        )

        ax.axvline(
            x_start - 0.5,
            color="lightgrey",
            linewidth=0.9,
            zorder=0
        )

        ax.axvline(
            x_end + 0.5,
            color= "lightgrey",
            linewidth=0.9,
            zorder=0
        )

    # Build legend
    legend_elements = [
        Line2D(
        [0], [0],
        marker=r"$\bullet$", #"x",
        linestyle="None",
        color="black",
        #fontsize=5,
        markersize=2,
        markeredgewidth=2,
        label="No seq.")
    ]

    # add one entry per modification color
    for mod_name, color in mod_colors.items():
        if mod_name in present_mods:
            legend_elements.append(Patch(facecolor=color, label=mod_name))

    ax.legend(
        handles=legend_elements,
        #title="Coverage / modification",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
        prop={"size": 5.1}
    )


    plt.xticks(
        range(len(data_matrix.columns)),
        data_matrix.columns,
        rotation=90,
        fontsize=5.0
    )

    plt.yticks(
        range(len(data_matrix.index)),
        data_matrix.index,
        fontsize=5.1,
    )

    ax.tick_params(axis='y', pad=1)


    plt.xlabel("Sprinzl coordinate", fontsize=6)
    ax.set_ylabel("tRNAs", fontsize=6, rotation=0)
    ax.yaxis.set_label_coords(-0.04, 1.022)
    ax.yaxis.label.set_verticalalignment('bottom')


    plt.tight_layout()
    plt.savefig("tRNA_Method_Consensus_heatmap.svg")
    plt.show()


############# Plotting #######################
consensus = "/home/alex/Documents/PHD/RNOME/Paper/consensus_Kandarp/tiered_tRNA.bed"
df_cons = pd.read_csv(consensus, sep="\t", comment='#',header=None)
df_cons.columns =["chrom","chromStart","chromEnd","name","tier","strand"]
df_cons = transform_tRNA_name(df_cons)

tRNAs = df_cons["chrom"].unique().tolist()

# Load sprinzl Excel file
df_sprinzl = pd.read_excel("tRNA_sprinzl.xlsx", engine="openpyxl")
df_sprinzl = df_sprinzl[df_sprinzl["tRNA_name"].notna()]
df_sprinzl = transform_tRNA_name(df_sprinzl)

df_sprinzl = df_sprinzl[df_sprinzl["tRNA_name"].isin(tRNAs)]
conversion_dict = build_sprinzl_mapping(df_sprinzl)

position_order = [
    "-1", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
    "16", "17", "17a", "18", "19", "20", "20a", "20b", "21", "22", "23", "24", "25", "26",
    "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
    "41", "42", "43", "44", "45",
    "e11", "e12", "e13", "e14", "e15", "e16", "e17",
    "e1", "e2", "e3", "e4", "e5",
    "e27", "e26", "e25", "e24", "e23", "e22", "e21",
    "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
    "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73",
    "74", "75", "76"
]

data_matrix = create_data_matrix(position_order, df_cons, conversion_dict, tRNAs,df_sprinzl)
plot_heatmap(data_matrix,df_cons["name"].unique().tolist())

