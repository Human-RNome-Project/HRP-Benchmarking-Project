#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# UMBS-seq (UBS-seq) core pipeline for the HRP m5C RNome samples
#   m5C_UMBS_seq_rep1 = YSL-5 , m5C_UMBS_seq_rep2 = YSL-6
#
# Reproduces 02_core_pipeline.ipynb as a standalone script.
# Steps: trim -> map -> mark_duplicates -> dedup -> split_refs -> dedup_filter
#        -> conv_unconv3n -> conv_unconv3n_filter -> qc_conversion
#
# The output of `conv_unconv3n_filter` is the input of m5C_calling.py; the
# output of `conv_unconv3n` (unfiltered) is used to estimate the background
# non-conversion rate.
# =============================================================================

# --- Configuration -----------------------------------------------------------
SAMPLES=("YSL-5" "YSL-6")
DATA_DIR="${DATA_DIR:-data}"
HISAT3N_INDEX="../reference/hisat3n/refs_gen"
REF_FASTA="../reference/refs_gen.fa"
SPIKEINS=("lambda" "pUC19")
THREADS=$(nproc)
START_STEP="trim"
STOP_STEP=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

STEPS=("trim" "map" "mark_duplicates" "dedup" "split_refs" "dedup_filter" \
       "conv_unconv3n" "conv_unconv3n_filter" "qc_conversion")

usage() {
    cat <<USAGE
Usage: $0 [options]

  --start-from STEP   resume from a step (default: trim)
  --stop-after STEP   stop after a step (default: run to the end)
                      steps: ${STEPS[*]}
  --sample SAMPLE     run a single sample (e.g. YSL-5)
  --threads N         number of threads (default: all cores)
  --index PATH        hisat-3n index base path
  --ref PATH          reference FASTA used for the index
  --data-dir PATH     directory holding {SAMPLE}_R1.fq.gz / _R2.fq.gz
  -h, --help
USAGE
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --start-from) START_STEP="$2"; shift 2 ;;
        --stop-after) STOP_STEP="$2"; shift 2 ;;
        --sample)     SAMPLES=("$2"); shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --index)      HISAT3N_INDEX="$2"; shift 2 ;;
        --ref)        REF_FASTA="$2"; shift 2 ;;
        --data-dir)   DATA_DIR="$2"; shift 2 ;;
        -h|--help)    usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

valid_step() { local s; for s in "${STEPS[@]}"; do [[ "$s" == "$1" ]] && return 0; done; return 1; }
valid_step "${START_STEP}" || { echo "Error: invalid --start-from '${START_STEP}'"; exit 1; }
[[ -z "${STOP_STEP}" ]] || valid_step "${STOP_STEP}" || { echo "Error: invalid --stop-after '${STOP_STEP}'"; exit 1; }

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# --- One-time index build (uncomment if the index does not exist yet) --------
# hisat-3n-build --base-change C,T -p ${THREADS} ${REF_FASTA} ${HISAT3N_INDEX}

# =============================================================================

step_trim() {
    log "=== trim ==="
    mkdir -p trim
    for S in "${SAMPLES[@]}"; do
        cutadapt -j 0 -n 2 \
            -a 'AGATCGGAAGAGCACACGTCT;e=0.15;o=4;anywhere;' \
            -A 'AGATCGGAAGAGCGTCGTGT;e=0.15;o=4;anywhere;' \
            --max-n=0 -q 15 --nextseq-trim=15 -m 20 \
            -u 5 -U 5 --rename='{id}_{r1.cut_prefix}{r2.cut_prefix}' \
            --too-short-output=trim/${S}_R1.fa_short \
            --too-short-paired-output=trim/${S}_R2.fa_short \
            -o trim/${S}_R1.fq.gz -p trim/${S}_R2.fq.gz \
            ${DATA_DIR}/${S}_R1.fq.gz ${DATA_DIR}/${S}_R2.fq.gz \
            > trim/${S}.report
        log "trim done: ${S}"
    done
}

step_map() {
    log "=== map ==="
    mkdir -p map
    for S in "${SAMPLES[@]}"; do
        hisat-3n --index ${HISAT3N_INDEX} -p ${THREADS} \
            --base-change C,T --mp 8,2 --no-spliced-alignment \
            --summary-file map/${S}.summary --new-summary \
            -1 trim/${S}_R1.fq.gz -2 trim/${S}_R2.fq.gz \
            --un-conc-gz map/${S}_R%.fq.gz \
            -S map/${S}.sam
        samtools view -@ ${THREADS} -F4 -b map/${S}.sam \
          | samtools sort -@ ${THREADS} --write-index -O BAM -o map/${S}.bam -
        rm map/${S}.sam
        log "map done: ${S}"
    done
}

step_mark_duplicates() {
    log "=== mark_duplicates ==="
    mkdir -p mark_duplicates
    export _JAVA_OPTIONS="-Xmx4g"
    for S in "${SAMPLES[@]}"; do
        gatk AddOrReplaceReadGroups \
            -I map/${S}.bam -O mark_duplicates/${S}.rg.bam \
            -LB 1 -PL ILLUMINA -PU 1 -SM ${S}
        gatk MarkDuplicates \
            -I mark_duplicates/${S}.rg.bam -O mark_duplicates/${S}.bam \
            -M mark_duplicates/${S}.tmp \
            --ASSUME_SORT_ORDER coordinate --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
        gatk CollectDuplicateMetrics \
            -I mark_duplicates/${S}.bam -M mark_duplicates/${S}.metrics
        grep -v -E '^#|^$' mark_duplicates/${S}.tmp > mark_duplicates/${S}.tsv
        samtools index mark_duplicates/${S}.bam
        rm mark_duplicates/${S}.rg.bam mark_duplicates/${S}.tmp
        log "mark_duplicates done: ${S}"
    done
}

step_dedup() {
    log "=== dedup ==="
    mkdir -p dedup
    for S in "${SAMPLES[@]}"; do
        samtools view -h -F 1024 mark_duplicates/${S}.bam -o dedup/${S}.bam
        samtools index dedup/${S}.bam
        log "dedup done: ${S}"
    done
}

step_split_refs() {
    log "=== split_refs ==="
    mkdir -p split
    local spike_re
    spike_re="$(IFS='|'; echo "${SPIKEINS[*]}")"
    for S in "${SAMPLES[@]}"; do
        local human_contigs
        human_contigs=$(samtools idxstats dedup/${S}.bam | cut -f1 \
                        | grep -v -E "^(${spike_re}|\*)$" | tr '\n' ' ')
        samtools view -@ ${THREADS} -b dedup/${S}.bam ${human_contigs} \
            -o split/${S}.human.bam
        samtools index split/${S}.human.bam
        for SP in "${SPIKEINS[@]}"; do
            if samtools idxstats dedup/${S}.bam | cut -f1 | grep -qx "${SP}"; then
                samtools view -@ ${THREADS} -b dedup/${S}.bam "${SP}" \
                    -o split/${S}.${SP}.bam
                samtools index split/${S}.${SP}.bam
            fi
        done
        log "split_refs done: ${S}"
    done
}

step_dedup_filter() {
    log "=== dedup_filter (read-level quality filter) ==="
    mkdir -p dedup_filter
    for S in "${SAMPLES[@]}"; do
        samtools view -@ ${THREADS} \
            -e '[XM] * 20 <= (qlen - sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]' \
            -O BAM -o dedup_filter/${S}.bam split/${S}.human.bam
        samtools index dedup_filter/${S}.bam
        log "dedup_filter done: ${S}"
    done
}

# per-position table; $1 = input BAM dir, $2 = output dir, $3 = "motif"|"nomotif"
conv_table() {
    local S="$1" IN="$2" OUT="$3" WITH_MOTIF="$4"
    mkdir -p "${OUT}"
    if [[ "${WITH_MOTIF}" == "motif" ]]; then
        samtools view -e "rlen<100000" -h "${IN}" \
          | hisat-3n-table -p ${THREADS} --unique-only --alignments - \
                --ref ${REF_FASTA} --output-name /dev/stdout --base-change C,T \
          | python "${SCRIPT_DIR}/annotate_motif.py" --ref ${REF_FASTA} --sample "${S}" \
          | bgzip -@ ${THREADS} -c > "${OUT}/${S}.tsv.gz"
    else
        samtools view -e "rlen<100000" -h "${IN}" \
          | hisat-3n-table -p ${THREADS} --unique-only --alignments - \
                --ref ${REF_FASTA} --output-name /dev/stdout --base-change C,T \
          | python "${SCRIPT_DIR}/annotate_motif.py" --sample "${S}" --no-motif \
          | bgzip -@ ${THREADS} -c > "${OUT}/${S}.tsv.gz"
    fi
}

step_conv_unconv3n() {
    log "=== conv_unconv3n (unfiltered; used for the background rate) ==="
    for S in "${SAMPLES[@]}"; do
        conv_table "${S}" "split/${S}.human.bam" "conv_unconv3n" "nomotif"
        log "conv_unconv3n done: ${S}"
    done
}

step_conv_unconv3n_filter() {
    log "=== conv_unconv3n_filter (input of the m5C calling) ==="
    for S in "${SAMPLES[@]}"; do
        conv_table "${S}" "dedup_filter/${S}.bam" "conv_unconv3n_filter" "motif"
        log "conv_unconv3n_filter done: ${S}"
    done
}

step_qc_conversion() {
    log "=== qc_conversion (spike-in conversion rate) ==="
    mkdir -p qc
    for S in "${SAMPLES[@]}"; do
        for SP in "${SPIKEINS[@]}"; do
            [[ -f "split/${S}.${SP}.bam" ]] || continue
            samtools view -h "split/${S}.${SP}.bam" \
              | hisat-3n-table -p ${THREADS} --unique-only --alignments - \
                    --ref ${REF_FASTA} --output-name /dev/stdout --base-change C,T \
              | awk -F'\t' 'NR>1 {conv+=$5; unconv+=$7}
                   END {if (conv+unconv > 0)
                          printf "%s\t%s\tconversion_rate\t%.6f\tdepth\t%d\n",
                                 "'"${S}"'", "'"${SP}"'", conv/(conv+unconv), conv+unconv}' \
              | tee -a qc/conversion_rate.tsv
        done
    done
}

# =============================================================================
log "samples: ${SAMPLES[*]} | threads: ${THREADS}"
log "index: ${HISAT3N_INDEX} | ref: ${REF_FASTA}"

running=false
for step in "${STEPS[@]}"; do
    [[ "$step" == "$START_STEP" ]] && running=true
    if [[ "$running" == "true" ]]; then
        step_${step}
        [[ -n "$STOP_STEP" && "$step" == "$STOP_STEP" ]] && { log "stopping after '${STOP_STEP}'"; break; }
    fi
done
log "pipeline complete."
