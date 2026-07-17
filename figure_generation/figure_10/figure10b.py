#%%
import pandas as pd
from matplotlib.lines import Line2D

from matplotlib.legend_handler import HandlerBase
from matplotlib import font_manager
import matplotlib as mpl

from plot28 import plot_28
from plot_5_8 import plot_5_8
from plot_18 import plot_18


class HandlerSquareCross(HandlerBase):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height,
                       fontsize, trans):

        square, cross = orig_handle

        x = xdescent + width / 2
        y = ydescent + height / 2

        sq = Line2D(
            [x], [y],
            marker='s',
            markersize=9,
            markerfacecolor='lightgrey',
            markeredgecolor='black',
            linestyle='None',
            transform=trans
        )

        cr = Line2D(
            [x], [y],
            marker='x',
            markersize=9,
            color='black',
            markeredgewidth=1.8,
            linestyle='None',
            transform=trans
        )

        return [sq, cr]

def mature_ribosome_chromname_renaming(code:str)->str:
    chromname_dict = {
    "hs_rRNA_18S":"18S",
    "hs_rRNA_28S":"28S",
    "hs_rRNA_5.8S":"5.8S",
    "hs_rRNA_5S":"5S",
    "18S":"18S",
    "28S":"28S",
    "5.8S":"5.8S",
    "5S":"5S"}
    try:
        return chromname_dict[code]
    except KeyError:
        return code

def align_chromnames(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aligns the chromosome names in the 'chrom' column of a DataFrame to a standard format using a mapping function.
    Parameters:
    df (pd.DataFrame): The DataFrame containing the 'chrom' column with chromosome codes.
    Returns:
    pd.DataFrame: A DataFrame with aligned chromosome names in the 'chrom' column.
    """
    df["chrom"] = df["chrom"].apply(mature_ribosome_chromname_renaming)
    return df

def modbasecode_to_mod(code:str)->str:
    base_modification_type_dict = {
    "a": "m6A",
    "69426": "Am",
    "17596": "I",
    "Ino":"I",
    "17802": "Y",
    "Y":"Y",
    "A":"A",
    "G":"G",
    "C":"C",
    "U":"U",
    "T":"U",
    "Tm":"Um",
    "m3C":"m3C",
    "psU":"Y",
    "Um":"Um",
    "Gm":"Gm",
    "Am":"Am",
    "Cm":"Cm",
    "69426": "Am",
    "19227": "Um",
    "19229": "Gm",
    "19228": "Cm",
    "m5C":"m5C",
    "20607": "m5C",
    "m6A":"m6A",
    "x6A":"m6A",
    "21891": "m6A",
    "xp3Cm":"xp3Cm",
    "m2xp7G":"m2xp7G",
    "m6G":"m6G",
    "xp4U":"xp4U",
    "ac7G":"ac7G",
    "89487":"m3U",
    "m3U":"m3U",
    "m5U":"m5U",
    "m4C":"m4C",
    "28284":"m6,6A",
    "0":"mchm5U",
    "xp6G":"xp6G",
    "xp7G":"xp7G",
    "m7A":"m7A",
    "m1acp3Y":"m1acp3Y",
    "m1ap3U":"m1acp3Y",
    "70989":"ac4C",
    "ac4C":"ac4C",
    "20794": "m7G",
    "m7G":"m7G",
    "m66A":"m6,6A",
    "Ym":"Ym",
    "16020":"m1A",
    "m1A":"m1A",
    "m1A ":"m1A",
    "m": "m5C",
    "143283":"m2,2,7G",
    "17562":"C",
    "16335":"A",
    "16704":"U",
    "16750":"G",
    "99990":"mxA",
    "99993":"mxU",
    "99991":"mxC",
    "99992":"mxG",
    "99995":"ncm5s2U",
    "184012":"m5Cm",
    "99996":"f5Cm",
    "99998":"mchm5U",
    "99997":"hm5Cm",
    "62875":"ms2i6A",
    "19289":"m2,2G",
    "62005":"ncm5U",
    "133071":"m6t6A",
    "75654":"cm5U",
    "234279":"f5C",
    "23774":"D",
    "20598":"mcm5U",
    "60193":"Q",
    "99994":"ncm5Um",
    "71693":"io6A",
    "62879":"ms2io6A",
    "143283":"m2,2,7G",
    "21440":"t6A",
    "27241":"mcmo5U",
    "71588":"acp3D",
    "62881":"i6A",
    "191041":"hm5C",
    "hm5C":"hm5C",
    "Q":"Q",
    "m2,2G":"m2,2G",
    "D":"D",
    "t6A":"t6A",
    "mchm5U":"mchm5U",
    "m5Cm":"m5Cm",
    "ms2i6A":"ms2i6A",
    "mcmo5U":"mcmo5U",
    "acp3D":"acp3D",
    "m6t6A":"m6t6A",
    "cm5U":"cm5U",
    "hm5Cm":"hm5Cm",
    "f5Cm":"f5Cm",
    "f5C":"f5C",
    "mcm5U":"mcm5U",
    "i6A":"i6A",
    "ncm5U":"ncm5U",
    "mA?":"mxA",
    "mU?":"mxU",
    "mC?":"mxC",
    "mG?":"mxG",
    "mxA":"mxA",
    "mxU":"mxU",
    "mxC":"mxC",
    "mxG":"mxG",
    "I":"I",
    "m2,2,7G":"m2,2,7G",
    "m6,6A":"m6,6A",
    "m62A":"m6,6A",
    'm1I':'m1I',
    "acp3U":"acp3U",
    "m1G":"m1G",
    "s2U":"s2U",
    "m2G":"m2G"
    }
    return base_modification_type_dict[str(code)]

def align_modification_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aligns the modification names in the 'name' column of a DataFrame to a standard format using a mapping function.
    Parameters:
    df (pd.DataFrame): The DataFrame containing the 'name' column with modification codes.
    Returns:
    pd.DataFrame: A DataFrame with aligned modification names in the 'name' column.
    """
    df["name"] = df["name"].apply(modbasecode_to_mod)
    return df

def load_bedrmod(path: str ="") -> pd.DataFrame:
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
        elif len(df.columns) == 13:
            df.columns = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","coverage","frequency","single_letter_code","mod_id"]
        elif len(df.columns) == 18:
            df.columns = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","coverage","frequency", "count_modified", "count_canonical", "count_other_mod", "count_delete", "count_fail", "count_diff", "count_nocall"]
        elif len(df.columns) == 14:
            df.columns = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","coverage","frequency", "unique_mapping", "frag_start", "frag_end"]
        elif len(df.columns) == 16:
            df.columns =["chrom","chromStart","chromEnd","name","strand","score","thickStart","thickEnd","itemRgb","coverage","frequency","n_dataset","score_std","coverage_std","frequency_std","method"]
        elif len(df.columns) == 15:
            df.columns =["chrom","chromStart","chromEnd","name","strand","score","thickStart","thickEnd","itemRgb","coverage","frequency","frequency_std","score_std","overlap_perc","overlap"]
        
    except ValueError as e:
        print(f"Error: {e}. The number of columns in the bedrmod file does not match the expected format.")
    df = align_modification_names(df)
    df = align_chromnames(df)
    
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

def filter_bedfile_massspec_rRNA(df: pd.DataFrame, chrom:str = None, chromStart:int = None, chromEnd:int = None, score_threshold:float=None,strand:str=None,coverage_threshold:int=None,frequency_threshold:float=None,unique_mapping:int=None,frag_start:int=None,frag_end:int=None) -> pd.DataFrame:
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
    if unique_mapping is not None:
        df = df[df["unique_mapping"] == unique_mapping]
    if frag_start is not None:
        df = df[df["frag_start"] >= frag_start]
    if frag_end is not None:
        df = df[df["frag_end"] >= frag_end]

    return df



def assign_tracks(fragment_dict):

    fragments = sorted(
        fragment_dict.values(),
        key=lambda x: x["start"]
    )

    track_ends = []

    for fragment in fragments:

        placed = False

        for track_idx, last_end in enumerate(track_ends):

            # no overlap
            if fragment["start"] > last_end:

                fragment["track"] = track_idx
                track_ends[track_idx] = fragment["end"]

                placed = True
                break

        # new track
        if not placed:

            fragment["track"] = len(track_ends)
            track_ends.append(fragment["end"])

def reference_overlap(consensus_file,reference_file,biotype):
    df_cons,comment_dict = load_bedrmod(consensus_file)
    df_cons["chromStart"] = df_cons["chromStart"].astype(int)

    df_ref,comment_dict_ref = load_bedrmod(reference_file)
    df_ref["chromStart"] = df_ref["chromStart"].astype(int)

    df_cons = df_cons[df_cons["chrom"] == biotype]
    #print(len(df_cons))
    df_ref = df_ref[df_ref["chrom"] == biotype]
    #print(df_ref)

    overlap_df = df_cons.merge(
        df_ref[["chromStart", "name"]],
        left_on=["chromStart", "name"],
        right_on=["chromStart", "name"],
        how="inner"
    )

    overlap_df = overlap_df.rename(columns={"name_x": "name"})

    non_overlap_df = df_cons.merge(
        df_ref[["chromStart", "name"]],
        left_on=["chromStart", "name"],
        right_on=["chromStart", "name"],
        how="left", indicator=True
    ).query('_merge == "left_only"').drop(columns=["_merge"])

    non_overlap_df = non_overlap_df.rename(columns={"name_x": "name"})



    return overlap_df, non_overlap_df

def reference_overlap_merge(consensus_files,reference_file,biotype):
    key_cols = ["chrom", "name", "chromStart"]
    mean_cols = ["frequency", "coverage"]

    
    dfs = []
    for file in consensus_files:
        df,comment_dict = load_bedrmod(file)
        dfs.append(df)
    
    combined = pd.concat(dfs, ignore_index=True)

    means = (
        combined
        .groupby(key_cols, as_index=False)[mean_cols]
        .mean()
    )

    # Keep first row (df1 has priority)
    df_cons = (
        combined
        .drop_duplicates(subset=key_cols, keep="first")
        .drop(columns=mean_cols)
        .merge(means, on=key_cols, how="left")
    )

    df_cons["chromStart"] = df_cons["chromStart"].astype(int)

    df_ref,comment_dict_ref = load_bedrmod(reference_file)
    df_ref["chromStart"] = df_ref["chromStart"].astype(int)

    df_cons = df_cons[df_cons["chrom"] == biotype]
    df_ref = df_ref[df_ref["chrom"] == biotype]
    #print(df_ref)

    overlap_df = df_cons.merge(
        df_ref[["chromStart", "name"]],
        left_on=["chromStart", "name"],
        right_on=["chromStart", "name"],
        how="inner"
    )

    overlap_df = overlap_df.rename(columns={"name_x": "name"})

    non_overlap_df = df_cons.merge(
        df_ref[["chromStart", "name"]],
        left_on=["chromStart", "name"],
        right_on=["chromStart", "name"],
        how="left", indicator=True
    ).query('_merge == "left_only"').drop(columns=["_merge"])

    non_overlap_df = non_overlap_df.rename(columns={"name_x": "name"})

    return overlap_df, non_overlap_df



def main():
    taoka_ref = "/home/alex/Documents/PHD/RNOME/Data-Jamboree/rRNA_session/Reference_Files/H.sapiens_ref_bedRmod_All.bed"
    SRS_consensus = "/home/alex/Documents/PHD/RNOME/Paper/Illumina/Illumina_combined_polyARNA_tRNA_rRNA.bed"
   
    overlap_5_8S, non_overlap_5_8S = reference_overlap(SRS_consensus,taoka_ref, "5.8S")
    overlap_18S, non_overlap_18S = reference_overlap(SRS_consensus,taoka_ref, "18S")
    overlap_28S, non_overlap_28S = reference_overlap(SRS_consensus,taoka_ref, "28S")

    plot_5_8(overlap_5_8S,non_overlap_5_8S)
    plot_18(overlap_18S, non_overlap_18S)
    plot_28(overlap_28S, non_overlap_28S)



main()