#%%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
import matplotlib.patches as patches
from matplotlib import font_manager
import matplotlib as mpl
from matplotlib.patches import Rectangle
from figure_6_helper_functions import MOD_COLORS, align_modification_names

# Set Arial globally
mpl.rcParams['font.family'] = 'Arial'
fonts = [f.name for f in font_manager.fontManager.ttflist]

mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['svg.fonttype'] = "none"
print("Arial" in fonts)



def remove_base_ends(df:pd.DataFrame):
    """Remove trailing base annotation from the 'chrom' column.

    Strips patterns like '_A', '_T', '_U', '_C', '_G' from chromosome names.

    Parameters:
    df (pd.DataFrame): DataFrame with a 'chrom' column.

    Returns:
    pd.DataFrame: DataFrame with cleaned 'chrom' values.
    """
    df["chrom"] = df['chrom'].str.replace(r'_[ATUCG]$', '', regex=True)
    return df

def load_bedrmod_tRNA(path: str ="") -> pd.DataFrame:
    """
    Loads a bedrmod file into a pandas DataFrame.
    Adds header to the loaded DataFrame.
    Parameters:
    path (str): The path to the bedrmod file. Default is an empty string, which will load the default file.
    Returns:
    pd.DataFrame: A DataFrame containing the data from the bedrmod file.
    dict: A dictionary containing the comment lines from the bedrmod file as key-value pairs.
     The key is the part before the = and the value is the part after the =.
     If a line does not contain an =, it will be skipped and a message will be printed.
     The dictionary will be returned along with the DataFrame.
    """
    #Read bedrmod file into a DataFrame, skipping comment lines and adding header
    df = pd.read_csv(path, sep="\t", header=None, comment='#')
    try:
        if len(df.columns) == 11:
            df.columns = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","coverage","frequency"]
        elif len(df.columns) == 18:
            df.columns = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","coverage","frequency", "count_modified", "count_canonical", "count_other_mod", "count_delete", "count_fail", "count_diff", "count_nocall"]
        elif len(df.columns) == 14:
            df.columns = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","coverage","frequency", "multi_mapping", "frag_start", "frag_end"]
        elif len(df.columns) == 16:
            df.columns =["chrom","chromStart","chromEnd","name","strand","score","thickStart","thickEnd","itemRgb","coverage","frequency","n_dataset","score_std","coverage_std","frequency_std","method"]
        elif len(df.columns) == 15:
            df.columns =["chrom","chromStart","chromEnd","name","strand","score","thickStart","thickEnd","itemRgb","coverage","frequency","frequency_std","score_std","overlap_perc","overlap"]
        
    except ValueError as e:
        print(f"Error: {e}. The number of columns in the bedrmod file does not match the expected format.")
    df = align_modification_names(df)
    
    #Extract comments from bedrmod file and store them in a list
    with open(path, 'r') as f:
        lines = f.readlines()
    comment_lines = [l.strip() for l in lines if l.startswith('#')]

    #Dictionairy with comment lines as key value pairs, 
    #where the key is the part before the = and the value is the part after the =
    comment_dict = {}
    for line in comment_lines:
        if line.startswith('#'):
            key_value = line[1:].split('=', 1)
            if len(key_value) == 2:
                key, value = key_value
                try:
                    comment_dict[key.strip()] = value.strip()
                except Exception as e:
                    print(f"Skipped line: {line}")
    return df, comment_dict


def filter_bedfile_massspec_tRNA(df: pd.DataFrame, chrom:str = None, chromStart:int = None, chromEnd:int = None, score_threshold:float=None,strand:str=None,coverage_threshold:int=None,frequency_threshold:float=None,multi_mapping:int=None,frag_start:int=None,frag_end:int=None) -> pd.DataFrame:
    """
    Filter a bedrmod DataFrame based on specified thresholds.
    Parameters:
    df (pd.DataFrame): The DataFrame to filter.
    chrom (str): The chromosome to filter by.
    chromStart (int): The start position of the chromosome to filter by.
    chromEnd (int): The end position of the chromosome to filter by.
    score_threshold (float): The minimum score threshold.
    strand (str): The strand to filter by.
    coverage_threshold (int): The minimum coverage threshold.
    frequency_threshold (float): The minimum frequency threshold.
    Returns:
    pd.DataFrame: A filtered DataFrame.
    """
    if chrom is not None:
        df = df[df["chrom"] == chrom]
    if chromStart is not None:
        df = df[df["chromStart"] >= chromStart]
    if chromEnd is not None:
        df = df[df["chromEnd"] <= chromEnd]
    if score_threshold is not None:
        df = df[df["score"] <= score_threshold]
    if strand is not None:
        df = df[df["strand"] == strand]
    if coverage_threshold is not None:
        df = df[df["coverage"] >= coverage_threshold]
    if frequency_threshold is not None:
        df = df[df["frequency"] >= frequency_threshold]
    if multi_mapping is not None:
        df = df[df["multi_mapping"] <= multi_mapping]
    if frag_start is not None:
        df = df[df["frag_start"] >= frag_start]
    if frag_end is not None:
        df = df[df["frag_end"] >= frag_end]

    return df


def get_fragment_coverage(file_names,base,end,tRNAs):
    fragment_dict = {}

    for tRNA in tRNAs:
        #[start,end,m_mapping]
        fragment_dict[tRNA] = []
    

    for file_name in file_names:
        df,comment_dict = load_bedrmod_tRNA(base + file_name + end)
        df = remove_base_ends(df)
        df = filter_bedfile_massspec_tRNA(df, score_threshold=0.05,frequency_threshold=0.1)
        df = transform_tRNA_name(df)
        
        for _, row in df.iterrows():
            tRNA = row["chrom"]
            frag_start = row["frag_start"]
            frag_end = row["frag_end"]
            multi_mappers = row["multi_mapping"]

            frag = [frag_start,frag_end,multi_mappers]
           
            # create fragment if not existing
            if tRNA in fragment_dict.keys():
                if frag not in fragment_dict[tRNA]:
                    fragment_dict[tRNA].append(frag)
            else:
                print(f"Warning: tRNA {tRNA} from file {file_name} not in tRNA list. Skipping fragment.")

    return fragment_dict

def build_sprinzl_mapping(df:pd.DataFrame, id_col="tRNA_name"):
    """
    Returns:
        dict[tRNA_name] -> dict[sequential_index] -> sprinzl_coordinate
    """
    #print(df.columns)

    df = df.drop(columns=["tRNAs"])

    # columns that represent Sprinzl positions (everything except metadata)
    sprinzl_cols = [c for c in df.columns if c not in [id_col, "tRNAs"]]
    #print("Sprinzl columns:", sprinzl_cols)
    result = {}

    for _, row in df.iterrows():
        tRNA = row[id_col]

        mapping = {}
        idx = 0  # sequential index over observed (non-NaN) bases

        for col in sprinzl_cols:
            val = row[col]

            # skip missing positions
            if pd.isna(val):
                continue

            mapping[idx] = col  # sprinzl coordinate label
            idx += 1

        result[tRNA] = mapping

    return result

def remove_base_ends(df:pd.DataFrame):
    """Remove trailing base annotation from the 'chrom' column.

    Strips patterns like '_A', '_T', '_U', '_C', '_G' from chromosome names.

    Parameters:
    df (pd.DataFrame): DataFrame with a 'chrom' column.

    Returns:
    pd.DataFrame: DataFrame with cleaned 'chrom' values.
    """
    df["chrom"] = df['chrom'].str.replace(r'_[ATUCG]$', '', regex=True)
    return df

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



def create_data_matrix(positions:list,df_cons:pd.DataFrame, fragments:dict,sprinzl_dict: dict,tRNAs:list,sprinzl_df:pd.DataFrame)->pd.DataFrame:
    # map genomic position -> Sprinzl coordinate
    df_cons["sprinzl_pos"] = df_cons.apply(
        lambda r: sprinzl_dict.get(r["chrom"], {}).get(r["chromStart"]),
        axis=1
    )

    #empty matrix

    data_matrix = pd.DataFrame(
        [[(0,0) for _ in range(len(positions))] for _ in range(len(tRNAs))],
        index=tRNAs,
        columns=positions
    )

    # fill matrix
    # (0,0) = not covered,(0,x)=no sequence, (1,1) = multi-mapping covered, (1,2) = unique mapping covered, (X,1)= multi-mapping modified, (X,2)= unique-mapping modified
    
    #fill mods
    for _, row in df_cons.iterrows():
 
        mapping = 2 if row["multi_mapping"] == 1 else 1
        mod = row["name"]
        value = (mod,mapping)
        data_matrix.loc[row["chrom"],str(row["sprinzl_pos"])] = value
        
    #fill coverage
    for tRNA, frags in fragments.items():
        for frag in frags:
            frag_start, frag_end, multi_mapping = frag

            # find sprinzl positions that fall within the fragment
            for pos in range(frag_start, frag_end + 1):
                #print(tRNA,pos)
                #CCA adapters are not in sprinzl
                try:
                    sprinzl_pos = sprinzl_dict[tRNA][pos]
                except:
                    continue

                if sprinzl_pos is not None:
                    current_value = data_matrix.loc[tRNA, str(sprinzl_pos)]
                    #print(current_value)
                    if current_value == (0,0):  # not covered
                        new_mapping = 1 if multi_mapping > 1 else 2
                        #print(new_mapping)
                        #print(tRNA, sprinzl_pos)
                        data_matrix.loc[tRNA, str(sprinzl_pos)] = (1,new_mapping) 
                else:
                    print(f"Warning: No Sprinzl mapping for {tRNA} position {pos}. Skipping coverage update.")
    #fill no sequence

    for _, row in sprinzl_df.iterrows():
        trna = row["tRNA_name"]

        for col in sprinzl_df.columns:
            if col == "tRNA_name":
                continue

            if pd.isna(row[col]):
                #print(trna, col)
                #print(data_matrix[trna, str(col)])
                data_matrix.loc[trna, str(col)] = (0,"x")

    return data_matrix

def cell_to_color(cell_value, mod_colors):
    if cell_value == (0, 0) or cell_value == (0, "x"):
        return "#FFFFFF"  # white for not covered
    elif cell_value[0] == 1:  # covered but not modified
        return "#D3D3D3"
    else:  # modified
        return mod_colors[cell_value[0]]  # use the modification color


def plot_heatmap_header(data_matrix: pd.DataFrame,present_mods:list):
    print("Present modifications:", present_mods)
    print(len(present_mods))
    data_matrix.index = data_matrix.index.str.removeprefix('hs_tRNA')
   

    # (0,0) = not covered,(0,x)=no sequence ,(1,1) = multi-mapping covered, (1,2) = unique mapping covered, (X,1)= multi-mapping modified, (X,2)= unique-mapping modified
    colour_matrix = np.empty((data_matrix.shape[0], data_matrix.shape[1], 3), dtype=float)

    for i in range(data_matrix.shape[0]):
        for j in range(data_matrix.shape[1]):
            color = cell_to_color(data_matrix.iloc[i, j],MOD_COLORS)
            rgb_color = to_rgb(color)
            colour_matrix[i, j] = rgb_color

    width_mm = 183
    height_mm = 100

    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
    
    ax.imshow(colour_matrix, aspect='auto', interpolation='none', resample=False)

    #multimapper mask
    multimapper_mask = np.vectorize(lambda x: x[1] == 1)(data_matrix.values)

    for i in range(multimapper_mask.shape[0]):
        for j in range(multimapper_mask.shape[1]):
            if multimapper_mask[i, j]:
                ax.add_patch(
                    patches.Rectangle(
                        (j - 0.5, i - 0.5),
                        1, 1,
                        fill=False,
                        hatch="/////",
                        edgecolor="white",
                        linewidth=0
                    )
                )
    
    #no seq mask
    no_seq_mask = np.vectorize(lambda x: x == (0,"x"))(data_matrix.values)
 
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

    #legend
    legend_elements = [
        Patch(facecolor="#D3D3D3", label="Covered"),
        Patch(facecolor="#FFFFFF", edgecolor="black", hatch="/////", label="Multi-mapping"),
        Line2D([0], [0],marker=r"$\bullet$",linestyle="None",color="black",markersize=2,markeredgewidth=2,label="No seq."),

    ]

    # add one entry per modification color
    for mod_name, color in MOD_COLORS.items():
        if mod_name in present_mods:
            legend_elements.append(Patch(facecolor=color, label=mod_name))

    ax.legend(
        handles=legend_elements,
        #title="Coverage / modification",
        bbox_to_anchor=(1.005, 1),
        loc="upper left",
        frameon=True,
        borderaxespad=0.1 ,
        prop={"size": 5.0}
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


    plt.tight_layout(pad=0.05)
    plt.savefig("tRNA_MS_heatmap.svg")
    plt.show()


def switch_names(df:pd.DataFrame):
    df["name"] = df["name"].replace({"mC?":"mxC","mG?":"mxG","mU?":"mxU","mA?":"mxA"})
    return df

############# Plotting #######################

raw_file_names = ["HRP_C_001_tRNA_003","HRP_C_008_tRNA_004","HRP_C_015_tRNA"]
raw_base = "/home/alex/Documents/PHD/RNOME/massspec/final_data/new_sunmission_18_06/raw_tRNA/" # path to data folder
end = ".bed"

consensus = "/home/alex/Documents/PHD/RNOME/Paper/data/tRNA_with_multimap.bed" #detailed tRNA cosnensus map including multimapping column

df_cons = pd.read_csv(consensus, sep="\t", comment='#',header=None)
df_cons.columns =["chrom","chromStart","chromEnd","name","strand","score","thickStart","thickEnd","itemRgb","coverage","frequency","frequency_std","score_std","multi_mapping","overlap_perc","overlap"]
df_cons = transform_tRNA_name(df_cons)

df_cons = switch_names(df_cons)

tRNAs = df_cons["chrom"].unique().tolist()

# Load sprinzl Excel file
df_sprinzl = pd.read_excel("tRNA_sprinzl.xlsx", engine="openpyxl")
df_sprinzl = df_sprinzl[df_sprinzl["tRNA_name"].notna()]
df_sprinzl = transform_tRNA_name(df_sprinzl)
df_sprinzl = df_sprinzl[df_sprinzl["tRNA_name"].isin(tRNAs)]

conversion_dict = build_sprinzl_mapping(df_sprinzl)
coverage_dict = get_fragment_coverage(raw_file_names,raw_base,end,tRNAs)

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

data_matrix = create_data_matrix(position_order, df_cons, coverage_dict, conversion_dict, tRNAs,df_sprinzl)
plot_heatmap_header(data_matrix,df_cons["name"].unique().tolist())

