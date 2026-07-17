#%%
import pandas as pd

from helper_functions import load_bedrmod, filter_bedfile_massspec, filter_consensus
from merging_functions import merge_positions_in_sample, create_consensus
from renaming_functions import reshape_to_final_bed, remove_base_ends

def run_consensus_creation(samples:list,out_file:str,min_samples: int = 1, q_score:float =0.05, freq:float = 10, min_overlap: int = 100, unique_mapping:bool = True, tRNA:bool=False) -> pd.DataFrame:
    """
    Full pipeline: load bedrmod samples, filter and merge per-sample positions, build consensus and write output.

    Parameters:
    samples (list): List of input bedrmod file paths.
    out_file (str): Path to output consensus bedrmod file.
    min_samples (int): Minimum samples to require a site in consensus.
    q_score (float): Score threshold used for filtering.
    freq (float): Frequency threshold used for filtering.
    min_overlap (int): Minimum overlap percentage for final filtering.
    unique_mapping (bool): If True, require unique mapping (unique_mapping==1) during filtering.

    Returns:
    None
    """
    dfs = {}
    for file_name in samples:
        df,comment_dict = load_bedrmod(file_name)
        if unique_mapping:
            df = filter_bedfile_massspec(df, score_threshold=q_score,unique_mapping=1,frequency_threshold=freq)
        else:
            df = filter_bedfile_massspec(df, score_threshold=q_score,frequency_threshold=freq)

        if tRNA:
            df = remove_base_ends(df)
        df = merge_positions_in_sample(df)
        df = df.drop(columns=["frag_start","frag_end"])
        # coverage is fragment  intensity. This is not comparable to other fragments or equvialent to read coverage. So we set it to 0. 
        df["coverage"] = 0
        #drop Inosine
        df = df[df["name"] != "I"]
     
        dfs[file_name] = df

    consensus, sample_count_list = create_consensus(dfs, min_samples=min_samples)
    filtered_consensus = filter_consensus(consensus, min_overlap_perc=min_overlap)
    #filtered_consensus = switch_methylation_naming(filtered_consensus)

    #save infromation rich file
    name = out_file.rsplit(".", 1)[0]
    extension = out_file.rsplit(".", 1)[-1]
    with open(name+"_statistics."+extension, "w") as f:
        for key, value in comment_dict.items():
            f.write(f"#{key}={value}\n")

        f.write("#" + "\t".join(filtered_consensus.columns) + "\n")

        # write BED data
        filtered_consensus.to_csv(f, sep="\t",header=False,index=False,lineterminator="\n")
    
    #save final file
    final_consensus = reshape_to_final_bed(filtered_consensus)

    with open(out_file, "w") as f:
        for key, value in comment_dict.items():
            f.write(f"#{key}={value}\n")

        f.write("#" + "\t".join(final_consensus.columns) + "\n")

        # write BED data
        final_consensus.to_csv(f, sep="\t",header=False,index=False,lineterminator="\n")



####################################################################################################
#running example 
#rRNA
samples_names = ["HRP_C_001_rRNA_001","HRP_C_002_rRNA_001",
              "HRP_C_005_rRNA_004+HRP_C_006_rRNA_004","HRP_C_014_rRNA_001"]

base = "" #fill in your datas strorage folder
end = ".bed"

samples=[]
for n in samples_names:
    samples.append(base+n+end)

out_file = "test_consensus.bed"
min_samples = 1
q_score = 0.05
frequency = 10

min_overlap = 100
unique_mapping = True

run_consensus_creation(samples,out_file,min_samples,q_score,frequency,min_overlap,unique_mapping)

#tRNA
samples_names =["HRP_C_001_tRNA_003","HRP_C_008_tRNA_004","HRP_C_015_tRNA"]

base = "" #fill in your datas strorage folder
end =  ".bed"
samples=[]
for n in samples_names:
    samples.append(base+n+end)

out_file = "test_consensus.bed"
min_samples = 1
q_score = 0.05
frequency = 10

min_overlap = 100
unique_mapping = False
tRNA = True

run_consensus_creation(samples,out_file,min_samples,q_score,frequency,min_overlap,unique_mapping,tRNA)
