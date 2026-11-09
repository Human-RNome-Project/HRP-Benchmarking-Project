#!/bin/bash
# usage: ./3_make_metadata.sh <annot_dir> <data_dir>
#
# Puts the <sample>.annot.txt tables written by `run` where ../pipeline.R reads them:
#   Rnome-mRNA-<N>.annot.txt -> <data_dir>/ETAM_RNA_mrna/metadata/metadata_form_mRNA_<N>.tsv
#   Rnome-IVT-<N>.annot.txt  -> <data_dir>/ETAM_RNA_IVT/metadata/metadata_form_IVT_<N>.tsv
# The content is unchanged (same columns). Other samples (Rnome-FTO*) are not used by
# pipeline.R and are skipped.
annot_dir=${1:-}
data_dir=${2:-}
if [[ -z $annot_dir || -z $data_dir ]]; then
  echo "usage: $0 <annot_dir> <data_dir>" >&2
  exit 1
fi

mkdir -p "$data_dir/ETAM_RNA_mrna/metadata" "$data_dir/ETAM_RNA_IVT/metadata"
n=0
for f in "$annot_dir"/Rnome-*.annot.txt; do
  [[ -e $f ]] || continue
  s=$(basename "$f" .annot.txt)          # e.g. Rnome-mRNA-1
  case $s in
    Rnome-mRNA-*) out="$data_dir/ETAM_RNA_mrna/metadata/metadata_form_mRNA_${s#Rnome-mRNA-}.tsv" ;;
    Rnome-IVT-*)  out="$data_dir/ETAM_RNA_IVT/metadata/metadata_form_IVT_${s#Rnome-IVT-}.tsv" ;;
    *) echo "skipped (not used by pipeline.R): $s"; continue ;;
  esac
  cp "$f" "$out" && echo "$s -> $out" && n=$((n + 1))
done
[[ $n -gt 0 ]] || { echo "No Rnome-mRNA-*/Rnome-IVT-* .annot.txt in $annot_dir" >&2; exit 1; }
