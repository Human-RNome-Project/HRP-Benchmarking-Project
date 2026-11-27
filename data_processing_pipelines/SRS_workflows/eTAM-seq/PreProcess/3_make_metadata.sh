#!/bin/bash
# usage: ./3_make_metadata.sh <annot_dir> <data_dir> [sample_sheet]
#
# Puts the <sample>.annot.txt tables written by `run` where /HRP-Benchmarking-Project/data_processing_pipelines/SRS_workflows/eTAM-seq/DownstreamProcess/pipeline.R
# reads them. Group and replicate come from the sample sheet (samples.tsv by default), not
# from the file names:
#   group mRNA, replicate N -> <data_dir>/ETAM_RNA_mrna/metadata/metadata_form_mRNA_<N>.tsv
#   group IVT,  replicate N -> <data_dir>/ETAM_RNA_IVT/metadata/metadata_form_IVT_<N>.tsv
# The content is unchanged (same columns). Samples of any other group are skipped.
annot_dir=${1:-}
data_dir=${2:-}
sheet=${3:-${SAMPLE_SHEET:-$(dirname "$0")/samples.tsv}}
if [[ -z $annot_dir || -z $data_dir ]]; then
  echo "usage: $0 <annot_dir> <data_dir> [sample_sheet]" >&2
  exit 1
fi
[[ -r $sheet ]] || { echo "Sample sheet not found: $sheet" >&2; exit 1; }

mkdir -p "$data_dir/ETAM_RNA_mrna/metadata" "$data_dir/ETAM_RNA_IVT/metadata" || exit 1
n=0
missing=0
while IFS=$'\t' read -r s_name s_group s_rep _; do
  [[ -z $s_name || $s_name == \#* || $s_name == sample ]] && continue
  case $s_group in
    mRNA) out="$data_dir/ETAM_RNA_mrna/metadata/metadata_form_mRNA_${s_rep}.tsv" ;;
    IVT)  out="$data_dir/ETAM_RNA_IVT/metadata/metadata_form_IVT_${s_rep}.tsv" ;;
    *) echo "skipped (group '$s_group' not used by pipeline.R): $s_name"; continue ;;
  esac
  f="$annot_dir/$s_name.annot.txt"
  if [[ ! -s $f ]]; then
    echo "Missing: $f" >&2
    missing=$((missing + 1))
    continue
  fi
  cp "$f" "$out" || exit 1
  echo "$s_name ($s_group, replicate $s_rep) -> $out"
  n=$((n + 1))
done < "$sheet"
[[ $missing -eq 0 ]] || { echo "$missing table(s) of $sheet missing in $annot_dir" >&2; exit 1; }
[[ $n -gt 0 ]] || { echo "No mRNA/IVT sample in $sheet" >&2; exit 1; }
