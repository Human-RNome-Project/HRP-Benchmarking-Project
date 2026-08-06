#!/usr/bin/env bash
set -euo pipefail
echo "Generating input JSON for ont_rRNA workflow..."
#POD5_DIR="/global/cfs/cdirs/m5243/analysis/ont/ont_rRNA/raw/GERBER_S_pod5/ONT_Gerber_IVT_rRNA_251218_r1_v1/media/johannes/Tower_2/RNOME/IVT_rRNA/"
#SAMPLE_ID="IVT_rRNA_SUP_Gerber_IVT_rRNA_251218_r1_v1"

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <POD5_DIR> <SAMPLE_ID>"
    exit 1
fi

POD5_DIR=$1
SAMPLE_ID=$2
REF_GENOME="/global/cfs/cdirs/m5243/references/GRCh38.primary_assembly.genome.fa"
REF_TX="/global/cfs/cdirs/m5243/references/hs_rRNAs_NR_046235.fa"
CPUS=12
BARCODED=false
BARCODING_MODEL="b04_RNA004"
BASECALLING_TYPE="sup"
USE_GPU="cuda:0,1,2,3"
GPUS=4

POD5_JSON=$(printf '%s\n' "$POD5_DIR"/*.pod5 | jq -R . | jq -s .)

OUTPUT_JSON="rRNA/inputs_${SAMPLE_ID}.json"

cat > "${OUTPUT_JSON}" <<EOF
{
  "ont_rRNA.pod5_files": ${POD5_JSON},
  "ont_rRNA.sample_id": "${SAMPLE_ID}",
  "ont_rRNA.ref_genome": "${REF_GENOME}",
  "ont_rRNA.ref_transcriptome": "${REF_TX}",
  "ont_rRNA.cpus": ${CPUS},
  "ont_rRNA.barcoded": ${BARCODED},
  "ont_rRNA.barcoding_model": "${BARCODING_MODEL}",
  "ont_rRNA.basecalling_type": "${BASECALLING_TYPE}",
  "ont_rRNA.use_gpu": "${USE_GPU}",
  "ont_rRNA.gpus": "${GPUS}",
  "ont_rRNA.modeldir": "/global/cfs/cdirs/m5243/analysis/ont/ont_mRNA_pilot/pipeline/modeldirectory"
}
EOF
