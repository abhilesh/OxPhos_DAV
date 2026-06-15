#!/bin/bash
#SBATCH --job-name=pagel_retry
#SBATCH --array=3,4,6,7,8,9,10,11,12,13,16,19,23,26,30
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=/work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs/pagel_retry_%A_%a.out
#SBATCH --error=/work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs/pagel_retry_%A_%a.err

set -euo pipefail

module load miniconda/2024-02-20
eval "$(conda shell.bash hook)"
conda activate OXPHOS_DAV

CHUNK=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
HPC_HOME="/home/ad2347/OxPhos_DAV"
HPC_WORK="/work/WorkDevelopmentC/ad2347/OxPhos_DAV"
PAGEL_JOBS="$HPC_WORK/data/phylo/pagel_jobs"
OUT_DIR="$HPC_WORK/results/phylo/pagel_results"
mkdir -p "$OUT_DIR" "$HPC_WORK/logs"

cd "$PAGEL_JOBS"

Rscript "$HPC_HOME/src/phylo/pagel_discrete.R" \
    "chunks/chunk_${CHUNK}.tsv" \
    "$OUT_DIR/results_${CHUNK}.tsv"

echo "Done: chunk $CHUNK"
