#!/usr/bin/env bash
# Dispatches to one of the three HRP MS-seq pipelines based on the first
# argument, so a single image can process one input file per container run
# (e.g. one HPC array-job task per mzML/sample folder).
set -euo pipefail

PIPELINE_ROOT=/opt/hrp-ms-seq

usage() {
    cat <<'EOF'
Usage: docker run [docker-opts] <image> <mode> [args...]

Modes:
  rrna        Run the rRNA NASE wrapper (NASE_rRNA/run_rRNA_analysis.py)
  trna        Run the tRNA NASE wrapper (NASE_tRNA/run_tRNA_analysis.py)
  consensus   Run the consensus creation pipeline (MS_consensus_creation/consensus_creation_main.py)

Any other first argument is executed as-is (e.g. "bash" for an interactive shell).

Run "docker run <image> <mode> help" (or with no further args) to print that
mode's own --help text.

Examples:
  docker run --rm -v "$PWD":/data image rrna /data/sample.mzML /data/ref.fasta \
      --precursor-tolerance 10 --product-tolerance 20 --output-dir /data/results

  docker run --rm -v "$PWD":/data image trna /data/sample.mzML /data/ref.fasta \
      --precursor-tolerance 10 --product-tolerance 20 --output-dir /data/results

  docker run --rm -v "$PWD":/data image consensus \
      --input-folder /data/samples --out-file /data/consensus.bed
EOF
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ "$1" = "help" ]; then
    usage
    exit 0
fi

mode="$1"
shift

# Forward bare "help"/no-args as --help so the script's own argparse help prints,
# instead of "help" being consumed as a positional argument.
if [ "$#" -eq 0 ] || [ "$1" = "help" ]; then
    set -- --help
fi

case "$mode" in
    rrna|rRNA|NASE_rRNA)
        exec python3 "${PIPELINE_ROOT}/NASE_rRNA/run_rRNA_analysis.py" "$@"
        ;;
    trna|tRNA|NASE_tRNA)
        exec python3 "${PIPELINE_ROOT}/NASE_tRNA/run_tRNA_analysis.py" "$@"
        ;;
    consensus|MS_consensus_creation)
        exec python3 "${PIPELINE_ROOT}/MS_consensus_creation/consensus_creation_main.py" "$@"
        ;;
    *)
        exec "$mode" "$@"
        ;;
esac
