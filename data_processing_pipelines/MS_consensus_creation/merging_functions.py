import pandas as pd

def merge_positions_in_sample(df:pd.DataFrame)->pd.DataFrame:
    """
    Merge multiple entries at the same genomic position within a single sample.

    Groups by chrom, chromStart, chromEnd, name and strand and aggregates fields
    using sensible defaults (min for score, first for coordinates/ids, mean for
    coverage/frequency).

    Parameters:
    df (pd.DataFrame): Sample DataFrame possibly containing duplicate position entries.

    Returns:
    pd.DataFrame: Aggregated DataFrame with unique positions per sample.
    """
    merged_df = df.groupby(["chrom", "chromStart", "chromEnd", "name", "strand"]).agg({
        "score": "min",
        "thickStart": "first",
        "thickEnd": "first",
        "itemRgb": "first",
        "coverage": "mean",
        "frequency": "mean",
        "unique_mapping": "max",
        "frag_start": "first",
        "frag_end": "first"
    }).reset_index()

    return merged_df


def create_consensus(dfs: dict[str:pd.DataFrame], min_samples: int = 1) -> pd.DataFrame:
    """
    Create a consensus DataFrame across multiple sample DataFrames.

    Combines per-sample tables, counts how many samples have a given position,
    computes mean frequency/coverage and basic statistics, and computes overlap
    metrics for downstream filtering.

    Parameters:
    dfs (dict[str, pd.DataFrame]): Mapping sample name -> sample DataFrame.
    min_samples (int): Minimum number of samples that must contain a modification at a position.

    Returns:
    tuple(pd.DataFrame, list): (consensus DataFrame, list of counts of samples per position)
    """
    #print(dfs)

    # add sample name
    tagged =  [df.assign(_sample=name) for name, df in dfs.items()]
    #print(tagged)
    combined = pd.concat(tagged, ignore_index=True)

    # Which samples have *any* entry at each position?
    samples_per_pos = (
        combined.groupby(["chrom","strand","chromStart"])["_sample"]
        .nunique()
        .rename("_samples_with_pos")
    )

    # Keep only modification sites
    filtered = combined[~combined["name"].isin(["A", "C", "G", "U"])].copy()

    # For each (position,strand, name): mean freq + count how many samples carry it
    grouped = (
        filtered.groupby(["chrom", "chromStart", "chromEnd", "name", "strand"])
        .agg(
        score=("score", "min"),
        thickStart=("thickStart", "first"),
        thickEnd=("thickEnd", "first"),
        itemRgb=("itemRgb", "first"),
        coverage=("coverage", "mean"),
        frequency=("frequency", "mean"),
        frequency_std=("frequency", lambda x: x.std(ddof=0)),
        scores_std=("score", lambda x: x.std(ddof=0)),
        unique_mapping=("unique_mapping", "max"),
        _sample_count=("_sample", "nunique"))
        .reset_index()
    )

    grouped = grouped.join(samples_per_pos, on=["chrom","strand", "chromStart"])
    grouped["overlap_perc"] = (grouped["_sample_count"] / grouped["_samples_with_pos"]) *100
    grouped["overlap"] = grouped["_sample_count"].astype(str) + "/" + grouped["_samples_with_pos"].astype(str)


    #filter for positions that are present in at least 2 samples
    grouped = grouped[grouped["_sample_count"] >= min_samples]

    #count how many samples are called 
    sample_count = grouped["_samples_with_pos"].value_counts().sort_index()
    sample_count_list = sample_count.tolist()

    return grouped.drop(columns=["_sample_count", "_samples_with_pos"]),sample_count_list