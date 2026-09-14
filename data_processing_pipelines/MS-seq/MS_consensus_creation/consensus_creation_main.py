#%%
from pathlib import Path
import pandas as pd
import argparse

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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a consensus bedrmod file from multiple sample bedrmod files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input-folder",
        required=True,
        help="Folder containing OpenMS fragment-sample files; all .bed and .bedrmod files in it are used as samples.",
    )

    parser.add_argument("--out-file", required=True, help="Path to output consensus bedrmod file.")
    parser.add_argument("--min-samples", type=int, default=1, help="Minimum samples covering a modification site required for accepting it in consensus.")
    parser.add_argument("--q-score", type=float, default=0.05, help="Score/q-value threshold used for filtering.")
    parser.add_argument("--freq", type=float, default=10, help="Frequency threshold used for filtering.")
    parser.add_argument("--min-overlap", type=int, default=100, help="Minimum overlap percentage for final filtering. (Percentage of samples that must have a modification at a specific site for it to be included in the consensus.)")

    mapping_group = parser.add_mutually_exclusive_group()
    mapping_group.add_argument(
        "--unique-mapping",
        dest="unique_mapping",
        action="store_true",
        help="Require unique mapping during filtering (default).",
    )
    mapping_group.add_argument(
        "--no-unique-mapping",
        dest="unique_mapping",
        action="store_false",
        help="Do not require unique mapping during filtering.",
    )
    parser.set_defaults(unique_mapping=True)

    parser.add_argument(
        "--tRNA",
        action="store_true",
        default=False,
        help="Treat samples as tRNA data (applies remove_base_ends).",
    )

    return parser.parse_args()

def find_sample_files(input_folder: str) -> list:
    """Return sorted paths of all .bed and .bedrmod files directly inside input_folder."""
    folder = Path(input_folder)
    if not folder.is_dir():
        raise NotADirectoryError(f"Input folder not found: {input_folder}")

    files = sorted(
        str(p) for p in folder.iterdir()
        if p.is_file() and p.suffix.lower() in (".bed", ".bedrmod")
    )

    if not files:
        raise FileNotFoundError(f"No .bed or .bedrmod files found in: {input_folder}")

    return files

def main():
    args = parse_args()

    samples = find_sample_files(args.input_folder)
    print(f"Found {len(samples)} sample file(s) in {args.input_folder}:")
    for s in samples:
        print(f"  - {s}")

    run_consensus_creation(
        samples=samples,
        out_file=args.out_file,
        min_samples=args.min_samples,
        q_score=args.q_score,
        freq=args.freq,
        min_overlap=args.min_overlap,
        unique_mapping=args.unique_mapping,
        tRNA=args.tRNA,
    )


if __name__ == "__main__":
    main()
