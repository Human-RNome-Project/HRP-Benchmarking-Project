#!/usr/bin/env bash
#SBATCH --job-name=bidseq_HRP_pU
#SBATCH --account=pi-taopan
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=36:00:00
#SBATCH --output=logs/bidseq_%j.out
#SBATCH --error=logs/bidseq_%j.err

set -euo pipefail

WORKDIR="$(realpath "$(dirname "$0")")"
SIF="${WORKDIR}/bidseq.sif"          # apptainer pull bidseq.sif docker://y9ch/bidseq
CONF="${WORKDIR}/data.yaml"

cd "${WORKDIR}"
mkdir -p logs

module load apptainer 2>/dev/null || true

# Bind every real path that appears in data.yaml (references, FASTQ, workdir).
BINDS="/project,/scratch"

echo "[$(date '+%F %T')] image:"
apptainer inspect "${SIF}" || true

apptainer run -B "${BINDS}" "${SIF}" -c "${CONF}" -j "${SLURM_CPUS_PER_TASK}"

echo "[$(date '+%F %T')] done"
ls -la workspace/filter_sites/
