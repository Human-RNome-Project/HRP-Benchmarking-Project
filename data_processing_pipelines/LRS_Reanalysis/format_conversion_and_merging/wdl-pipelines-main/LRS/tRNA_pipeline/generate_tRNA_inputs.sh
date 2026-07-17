#!/usr/bin/env bash
set -euo pipefail
echo "Generating input JSON for ont_tRNA workflow..."
POD5_DIR="/global/cfs/cdirs/m5243/analysis/ont/all_runs/raw/NOVOA_E/ONT_tRNA_native_vs_IVT/bc_9"
SAMPLE_ID="HRP_A_021_3_IVT_tRNA_003_SUP"
REF_TX="/global/cfs/cdirs/m5243/references/hg38.rtRNA.with_oligos.fa"
REF_GENOME="/global/cfs/cdirs/m5243/references/GRCh38.primary_assembly.genome.fa"
CPUS=12
BARCODED=false
BARCODING_MODEL="b04_RNA004"
BASECALLING_TYPE="sup"
USE_GPU="cuda:0,1,2,3"
GPUS=4

POD5_JSON=$(printf '%s\n' "$POD5_DIR"/*.pod5 | jq -R . | jq -s .)

OUTPUT_JSON="tRNA/inputs_${SAMPLE_ID}.json"

cat > "${OUTPUT_JSON}" <<EOF
{
  "ont_tRNA.pod5_files": ${POD5_JSON},
  "ont_tRNA.sample_id": "${SAMPLE_ID}",
  "ont_tRNA.ref_genome": "${REF_GENOME}",
  "ont_tRNA.ref_transcriptome": "${REF_TX}",
  "ont_tRNA.cpus": ${CPUS},
  "ont_tRNA.barcoded": ${BARCODED},
  "ont_tRNA.barcoding_model": "${BARCODING_MODEL}",
  "ont_tRNA.basecalling_type": "${BASECALLING_TYPE}",
  "ont_tRNA.use_gpu": "${USE_GPU}",
  "ont_tRNA.gpus": "${GPUS}",
  "ont_tRNA.modeldir": "/global/cfs/cdirs/m5243/analysis/ont/ont_mRNA_pilot/pipeline/modeldirectory"
}
EOF
