#!/usr/bin/env bash
set -euo pipefail
echo "Generating input JSON for LRS_tRNA pipeline"
POD5_DIR="/global/cfs/cdirs/m5243/benchmarking/ont/NOVOA_E/ONT_tRNA_native_vs_gentegra/raw/bc_1"
SAMPLE_ID="tRNA_example_1_fast"
REF_TRANSCRIPTOME="/global/cfs/cdirs/m5243/references/Novoa_V2/hg38.rtRNA.with_oligos.fa"
CPUS=12
BARCODED=false
BARCODING_MODEL="b04_RNA004"
BASECALLING_TYPE="fast"
USE_GPU="cuda:0,1,2,3"
GPUS=4

POD5_JSON=$(printf '%s\n' "$POD5_DIR"/*.pod5 | jq -R . | jq -s .)

OUTPUT_JSON="tRNA_inputs_${SAMPLE_ID}.json"

cat > "${OUTPUT_JSON}" <<EOF
{
  "ont_tRNA.pod5_files": ${POD5_JSON},
  "ont_tRNA.sample_id": "${SAMPLE_ID}",
  "ont_tRNA.ref_transcriptome": "${REF_TRANSCRIPTOME}",
  "ont_tRNA.cpus": ${CPUS},
  "ont_tRNA.barcoded": ${BARCODED},
  "ont_tRNA.barcoding_model": "${BARCODING_MODEL}",
  "ont_tRNA.basecalling_type": "${BASECALLING_TYPE}",
  "ont_tRNA.use_gpu": "${USE_GPU}",
  "ont_tRNA.gpus": "${GPUS}",
  "ont_tRNA.modeldir": "/global/cfs/cdirs/m5243/analysis/ont/ont_mRNA_pilot/pipeline/modeldirectory"
}
EOF
