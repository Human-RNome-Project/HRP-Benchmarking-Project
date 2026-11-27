#!/bin/bash
# usage: ./match_hub_names.sh <original_fq_dir> <hub_fq_dir>
#
# Finds which released file (HRP_B_002_<N>_L<l>_R<r>.fastq[.gz]) is which original file
# (Rnome-<group>-<N>_L<l>_R<r>_001.fastq.gz) by comparing the first read name of each file.
# An Illumina read name holds the instrument, run, flowcell, lane, tile and x:y position, so
# two files sharing their first read name hold the same reads. Prints one line per hub prefix
# with the original sample it matches, ready to paste into the fastq_prefix column of samples.tsv.
orig_dir=${1:-}
hub_dir=${2:-}
if [[ -z $orig_dir || -z $hub_dir ]]; then
  echo "usage: $0 <original_fq_dir> <hub_fq_dir>" >&2
  exit 1
fi

first_name() { # first read name of a FASTQ file, without the comment after the space
  case $1 in
    *.gz) gzip -dc "$1" 2>/dev/null | head -1 ;;
    *)    head -1 "$1" ;;
  esac | cut -d' ' -f1
}
# <prefix>_L<l>_R<r>[_001].fastq[.gz] -> "<prefix> L<l>_R<r>"
split_name() {
  basename "$1" | sed -E 's/_(L[0-9]+)_(R[12])(_001)?\.(fastq|fq)(\.gz)?$/ \1_\2/'
}

declare -A ORIG
while read -r f; do
  ORIG[$(first_name "$f")]=$f
done < <(find "$orig_dir" -maxdepth 2 -regextype posix-extended -regex '.*_L[0-9]+_R[12](_001)?\.(fastq|fq)(\.gz)?')
[[ ${#ORIG[@]} -gt 0 ]] || { echo "No FASTQ in $orig_dir" >&2; exit 1; }

declare -A MATCH
status=0
while read -r f; do
  read -r hub_prefix hub_part < <(split_name "$f")
  o=${ORIG[$(first_name "$f")]:-}
  if [[ -z $o ]]; then
    echo "no match: $f" >&2; status=1; continue
  fi
  read -r orig_prefix orig_part < <(split_name "$o")
  if [[ $hub_part != "$orig_part" ]]; then
    echo "lane/read differ: $f <-> $o" >&2; status=1; continue
  fi
  if [[ -n ${MATCH[$hub_prefix]:-} && ${MATCH[$hub_prefix]} != "$orig_prefix" ]]; then
    echo "inconsistent: $hub_prefix matches ${MATCH[$hub_prefix]} and $orig_prefix" >&2; status=1
  fi
  MATCH[$hub_prefix]=$orig_prefix
done < <(find "$hub_dir" -maxdepth 2 -regextype posix-extended -regex '.*_L[0-9]+_R[12](_001)?\.(fastq|fq)(\.gz)?' | sort)

printf 'hub_prefix\toriginal_sample\n'
for k in $(printf '%s\n' "${!MATCH[@]}" | sort); do
  printf '%s\t%s\n' "$k" "${MATCH[$k]}"
done
exit $status
