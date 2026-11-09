#!/bin/bash

set -euo pipefail
export LC_ALL=C          # sort and join must use the same byte-wise ordering

usage() {
  echo "Usage: $0 BAYES_FILE KL_FILE NATIVE_FILE OUTPUT_FILE" >&2
  exit 1
}

[ $# -eq 4 ] || usage

BAY=$1
KL=$2
NATIVE=$3
OUT=$4

# input files must exist and not be empty
for f in "$BAY" "$KL" "$NATIVE"; do
  [ -s "$f" ] || { echo "Error: '$f' does not exist or is empty" >&2; exit 1; }
done

# output must not overwrite an input
for f in "$BAY" "$KL" "$NATIVE"; do
  if [ -e "$OUT" ] && [ "$OUT" -ef "$f" ]; then
    echo "Error: output file '$OUT' is the same as input '$f'" >&2
    exit 1
  fi
done

TAB=$'\t'
TMP=$(mktemp -d "${SCRATCH:-.}/join_tmp_XXXXXX")   # temp space on scratch, not $HOME
trap 'rm -rf "$TMP"' EXIT                          # clean up even if the script fails
SORTOPTS=(-t "$TAB" -S "${SORT_MEM:-2G}" -T "$TMP")

# lookup files -> "chrom|start|end|name|strand <TAB> score", sorted by key
# (modkit codes in the name column are translated via the #modification_names header)
mkkey='BEGIN{FS=OFS="\t"}
  /^#modification_names=/ {
    s = $1; sub(/^#modification_names=/, "", s)
    n = split(s, a, ",")
    for (i = 1; i <= n; i++) { split(a[i], p, ":"); map[p[1]] = p[2] }
    next
  }
  /^#/ { next }
  { nm = ($4 in map) ? map[$4] : $4
    print $1"|"$2"|"$3"|"nm"|"$6, $5 }'
awk "$mkkey" "$BAY" | sort "${SORTOPTS[@]}" -k1,1 > "$TMP/bay"
awk "$mkkey" "$KL"  | sort "${SORTOPTS[@]}" -k1,1 > "$TMP/kl"

# native file -> "key <TAB> original line number <TAB> full line", sorted by key
awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1"|"$2"|"$3"|"$4"|"$6, NR, $0}' "$NATIVE" \
  | sort "${SORTOPTS[@]}" -k1,1 > "$TMP/native"

# two left joins, restore the original order, drop the helper columns
join -t "$TAB" -a 1 -e . -o auto "$TMP/native" "$TMP/bay" \
  | join -t "$TAB" -a 1 -e . -o auto - "$TMP/kl" \
  | sort "${SORTOPTS[@]}" -k2,2n \
  | cut -f3- > "$TMP/body"

# header: copy all # lines from the native file, fixing the column-name line
awk -v q="'" '
  /^#\[/ { sub(q "score" q ",", q "score_native" q ",")
           sub(/\][ \t\r]*$/, ", " q "score_bayesian" q ", " q "score_kl" q "]") }
  /^#/   { print }
' "$NATIVE" > "$OUT"
cat "$TMP/body" >> "$OUT"

echo "Done: wrote $OUT" >&2
