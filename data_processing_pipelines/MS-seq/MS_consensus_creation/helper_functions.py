import pandas as pd

from renaming_functions import align_chromnames,align_modification_names

def filter_bedfile_massspec(df: pd.DataFrame, chrom:str = None, chromStart:int = None, chromEnd:int = None, score_threshold:float=None,strand:str=None,coverage_threshold:int=None,frequency_threshold:float=None,unique_mapping:int=None,frag_start:int=None,frag_end:int=None) -> pd.DataFrame:
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
    #Read bed file into a DataFrame, skipping comment lines and adding header
    df = pd.read_csv(path, sep="\t", header=None, comment='#')
    try:
        if len(df.columns) == 14:
            df.columns = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","coverage","frequency", "unique_mapping", "frag_start", "frag_end"]
    
    except ValueError as e:
        print(f"Error: {e}. The number of columns in the bedrmod file does not match the expected format.")
    df = align_modification_names(df)
    #df = align_chromnames(df)
    
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

def filter_consensus(consensus:pd.DataFrame, min_overlap_perc:float)->pd.DataFrame:
    """
    Filter consensus DataFrame by minimum overlap percentage.

    Parameters:
    consensus (pd.DataFrame): Consensus table containing an 'overlap_perc' column.
    min_overlap_perc (float): Minimum overlap percentage to keep a row.

    Returns:
    pd.DataFrame: Filtered consensus DataFrame.
    """
    return consensus[consensus["overlap_perc"] >= min_overlap_perc]