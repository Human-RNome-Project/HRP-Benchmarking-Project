#%%
import pandas as pd
import matplotlib as mpl
from matplotlib import font_manager
from general_plot import plot_broken_axis 
from translation_dicts import RENAME

# Set Arial globally
mpl.rcParams['font.family'] = 'Arial'
fonts = [f.name for f in font_manager.fontManager.ttflist]

mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['svg.fonttype'] = "none"
print("Arial" in fonts)


def loader_HRP_C_003(path:str):
    df = pd.read_excel(path, skiprows=10, engine="openpyxl")
    # drop all Unnamed columns
    df = df.loc[:, ~df.columns.astype(str).str.contains("^Unnamed")]
    df.columns = df.columns.str.strip().str.replace(":", "", regex=False)
    df = df.rename(columns={"Modification":"name"})
    df["name"] = df["name"].replace(RENAME)
    df = df[["name","replica 1", "replica 2", "replica 3"]]
    
    return df

def loader_HRP_C_021_022(path:str, type:str):
    def load_table(skiprows, nrows):
        df = pd.read_excel(path,skiprows=skiprows, nrows=nrows, usecols="A:W").rename(columns={"Unnamed: 0": "sample"})
        df.set_index("sample", inplace=True)

        df = df.T
        df.index.name = "name"

        df.columns = [f"021_022_rep_{i+1}" for i in range(df.shape[1])]
        df = df.reset_index()
        df["name"] = df["name"].replace(RENAME)
        return df

    df_total = load_table(skiprows=4, nrows=4)
    df_28S   = load_table(skiprows=12, nrows=3)
    df_18S   = load_table(skiprows=19, nrows=3)

    if type=="RNA":
        return df_total
    
    elif type == "28S":
        return df_28S
    
    elif type =="18S":
        return df_18S


def load_mean_data(path:str):
    raw = pd.read_excel(path, header=None)
    raw = raw.dropna(axis=1, how="all")   # drop the empty leading column
    raw.columns = range(raw.shape[1])

    dfs = {}
    current_rna = None
    header_row = None
    rows = []

    for i, row in raw.iterrows():
        non_null = row.dropna()

        # Section header row (e.g. "5S rRNA") — exactly one cell filled
        if len(non_null) == 1 and str(non_null.iloc[0]).strip() != "Modification":
            if current_rna is not None and rows:
                dfs[current_rna] = pd.DataFrame(rows, columns=header_row)
            current_rna = str(non_null.iloc[0]).strip()
            rows, header_row = [], None
            continue

        # Column header row ("Modification", "Mod / 1k rN", "Error", "CV")
        if len(non_null) and non_null.iloc[0] == "Modification":
            header_row = [str(c).strip() for c in row.tolist()]
            continue

        # Data row
        if current_rna is not None and header_row is not None and pd.notna(row.iloc[0]):
            rows.append(row.tolist())

    if current_rna is not None and rows:
        dfs[current_rna] = pd.DataFrame(rows, columns=header_row)

    for name, df in dfs.items():
        df["Modification"] = df["Modification"].astype(str).str.strip()
        df =  df.rename(columns={"Modification":"name","Mod / 1k rN":"mean","Error":"std"}) 
        df["name"] = df["name"].replace(RENAME)
        dfs[name] = df.reset_index(drop=True)

    return dfs


def merge_replicates(dfs:list):
    merged = pd.concat(dfs, axis=0).groupby("name").first()
    merged = merged.fillna(0)

    return merged

def loader_HRP_C_004_005_006_007_008_010(path:str, type:str):
    df = pd.read_excel(path, engine="openpyxl",skiprows=2, header=None)
    columns = [
        "name",

        # 5.8S
        "5.8S_1", "5.8S_2", "5.8S_3", "mean_5.8S", "std_5.8S",

        # 5S
        "5S_1", "5S_2", "5S_3", "mean_5S", "std_5S",

        # 28S
        "28S_1", "28S_2", "mean_28S", "std_28S",

        # 18S
        "18S_1", "18S_2", "18S_3", "mean_18S", "std_18S",

        # tRNA
        "tRNA_1", "tRNA_2", "tRNA_3", "mean_tRNA", "std_tRNA",

        # total RNA
        "RNA_1", "RNA_2", "RNA_3", "mean_RNA", "stab_RNA"]
    df.columns = columns
    cols = df.columns[df.columns.str.contains(type)]
    df = df[["name"] + list(cols)]
    df["name"] = df["name"].replace(RENAME)

    #print(df)
    return df
    

##import files
#All files were generated via the normalization excell sheets provided in the overview; data load can therefore be different since, manual  outputs were generated

HRP_C_004_005_006_007_008_010_details = "HRP-C-004_005_006_007_008_0010_intermediate_results.xlsx"
mean_data = "HRP-C-mean_std_combined.xlsx"
HRP_C_021_022_details = "HRP-C-021_022_intermediate_results.xlsx"
HRP_C_003_df = loader_HRP_C_003("HRP-C-003_intermediate_results.xlsx")

cut_point= 9  
height_ratio = (1,2)


#load mean data
dfs = load_mean_data(mean_data)

#5S
mean_df = dfs["5S rRNA"]
df = loader_HRP_C_004_005_006_007_008_010(HRP_C_004_005_006_007_008_010_details,"5S")
ratio = (183/3,50)
plot_broken_axis(df,mean_df,cut_point,"MS_RQ2_5S.svg",ratio,"5S",height_ratio)

#5.8S
mean_df = dfs["5.8S rRNA"]
df = loader_HRP_C_004_005_006_007_008_010(HRP_C_004_005_006_007_008_010_details,"5.8S")
ratio = (183/3,50)
plot_broken_axis(df,mean_df,cut_point,"MS_RQ2_5.8S.svg",ratio,"5.8S",height_ratio)

#18S
mean_df = dfs["18S rRNA"]
detailed_df = loader_HRP_C_004_005_006_007_008_010(HRP_C_004_005_006_007_008_010_details,"18S")
detailed_df = detailed_df[["name","18S_1", "18S_2", "18S_3"]]
HRP_C_021_022_df = loader_HRP_C_021_022(HRP_C_021_022_details,"18S")
merged_df = merge_replicates([detailed_df,HRP_C_021_022_df])
#exception for m66A, error in measurement in HRPC_022 only data from HRP_C_005 was taken; rest manually deleted
ratio = (183/3,50)

plot_broken_axis(merged_df,mean_df,cut_point,"MS_RQ2_18S.svg",ratio,"18S",height_ratio)

#28S
mean_df = dfs["28S rRNA"]
detailed_df = loader_HRP_C_004_005_006_007_008_010(HRP_C_004_005_006_007_008_010_details,"28S")
detailed_df = detailed_df[["name","28S_1", "28S_2"]]
HRP_C_021_022_df = loader_HRP_C_021_022(HRP_C_021_022_details,"28S")
merged_df = merge_replicates([detailed_df,HRP_C_021_022_df])
ratio = (90,50)

plot_broken_axis(merged_df,mean_df,cut_point,"MS_RQ2_28S.svg",ratio,"28S",height_ratio)

#tRNA
mean_df = dfs["tRNA"]
detailed_df = loader_HRP_C_004_005_006_007_008_010(HRP_C_004_005_006_007_008_010_details,"tRNA")
detailed_df = detailed_df[["name","tRNA_1", "tRNA_2","tRNA_3"]]
merged_df = merge_replicates([detailed_df,HRP_C_003_df])
ratio = (90,50)

plot_broken_axis(merged_df,mean_df,cut_point,"MS_RQ2_tRNA.svg",ratio,"tRNA",height_ratio)

#total
mean_df = dfs["total RNA"]
detailed_df = loader_HRP_C_004_005_006_007_008_010(HRP_C_004_005_006_007_008_010_details,"RNA")
detailed_df = detailed_df[["name","RNA_1", "RNA_2","RNA_3"]]
HRP_C_021_022_df = loader_HRP_C_021_022(HRP_C_021_022_details,"RNA")

merged_df = merge_replicates([detailed_df,HRP_C_021_022_df])
ratio = (90,50)

plot_broken_axis(merged_df,mean_df,cut_point,"MS_RQ2_totalRNA.svg",ratio,"totalRNA",height_ratio)
