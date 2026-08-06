#!/usr/bin/env bash
set -euo pipefail
echo "Generating input JSON for ont_mRNA_pilot workflow..."
POD5_DIR="/global/cfs/cdirs/m5243/analysis/ont/all_runs/raw/RAGOUSSIS_J/ONT_GM12878_directRNA_12022026_03"
SAMPLE_ID="polyA_SUP_RAGOUSSIS_J_ONT_GM12878_directRNA_12022026_03"
REF_GENOME="/global/cfs/cdirs/m5243/references/GRCh38.primary_assembly.genome.fa"
REF_TX="/global/cfs/cdirs/m5243/references/longest_gencode.v49.transcripts.fa"
CPUS=12
BARCODED=false
BARCODING_MODEL="b04_RNA004"
BASECALLING_TYPE="sup"
USE_GPU="cuda:0,1,2,3"
GPUS=4

POD5_JSON=$(printf '%s\n' "$POD5_DIR"/*.pod5 | jq -R . | jq -s .)

OUTPUT_JSON="polyA_RNA/inputs_${SAMPLE_ID}.json"

cat > "${OUTPUT_JSON}" <<EOF
{
  "ont_mRNA_pilot.pod5_files": ${POD5_JSON},
  "ont_mRNA_pilot.sample_id": "${SAMPLE_ID}",
  "ont_mRNA_pilot.ref_genome": "${REF_GENOME}",
  "ont_mRNA_pilot.ref_transcriptome": "${REF_TX}",
  "ont_mRNA_pilot.cpus": ${CPUS},
  "ont_mRNA_pilot.barcoded": ${BARCODED},
  "ont_mRNA_pilot.barcoding_model": "${BARCODING_MODEL}",
  "ont_mRNA_pilot.basecalling_type": "${BASECALLING_TYPE}",
  "ont_mRNA_pilot.use_gpu": "${USE_GPU}",
  "ont_mRNA_pilot.gpus": "${GPUS}",
  "ont_mRNA_pilot.modeldir": "/global/cfs/cdirs/m5243/analysis/ont/ont_mRNA_pilot/pipeline/modeldirectory"
}
EOF
