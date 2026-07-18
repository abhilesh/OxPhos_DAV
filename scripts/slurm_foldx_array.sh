#!/bin/bash
#SBATCH --job-name=foldx_ddg
#SBATCH --array=0-59
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs/foldx_%A_%a.out
#SBATCH --error=/work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs/foldx_%A_%a.err
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=abhilesh7@gmail.com
#
# SLURM array for FoldX BuildModel/AnalyseComplex ΔΔG (01_foldx_ddg.py).
# 60 tasks, round-robin over the input pair list; each pair costs ~8 min
# (intraprotein: 3 BuildModel runs) to ~16 min (interprotein: +AnalyseComplex).
# For the full 694-pair prioritized_pairs.csv with 629 pairs remaining after
# the local 64-pair run, 60 chunks × ~11 pairs/chunk ≈ 2 hr per task at the
# mixed intra/inter average — --time=04:00:00 above leaves headroom.
#
# Writes: results/mutagenesis/foldx_ddg_chunk{NNNN}.csv (one per array task)
#
# ── SETUP (run locally first) ────────────────────────────────────────────────
#   rsync -av results/mutagenesis/prioritized_pairs.csv \
#       lambda:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/mutagenesis/
#   rsync -av results/mutagenesis/foldx_ddg.csv \
#       lambda:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/mutagenesis/  # for --skip-existing
#   rsync -av data/structures/ \
#       lambda:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/data/structures/
#   rsync -av tools/foldx/ lambda:/home/ad2347/OxPhos_DAV/tools/foldx/
#   rsync -av src/    lambda:/home/ad2347/OxPhos_DAV/src/
#   rsync -av setup.py lambda:/home/ad2347/OxPhos_DAV/
#   rsync -av scripts/slurm_foldx_array.sh lambda:/home/ad2347/OxPhos_DAV/scripts/
#
#   # FoldX binary must be executable on the remote host (statically-linked
#   # Linux x86-64 build — no rebuild needed, but re-chmod after rsync):
#   ssh lambda 'chmod +x /home/ad2347/OxPhos_DAV/tools/foldx/foldx_*'
#
# ── SUBMIT ───────────────────────────────────────────────────────────────────
#   mkdir -p /work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs
#   mkdir -p /work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/mutagenesis
#   sbatch /home/ad2347/OxPhos_DAV/scripts/slurm_foldx_array.sh
#
# ── COPY CHUNKS BACK ─────────────────────────────────────────────────────────
#   rsync -av \
#       lambda:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/mutagenesis/foldx_ddg_chunk*.csv \
#       results/mutagenesis/
#
# ── MERGE LOCALLY ────────────────────────────────────────────────────────────
#   docker run --rm -v $(pwd):/app oxphos_dav_analysis conda run -n oxphos_dav \
#       python src/mutagenesis/04_merge_foldx_chunks.py
#   docker run --rm -v $(pwd):/app oxphos_dav_analysis conda run -n oxphos_dav \
#       python src/mutagenesis/03_compile_targets.py   # refresh final_targets.csv
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

HPC_HOME="/home/ad2347/OxPhos_DAV"
HPC_WORK="/work/WorkDevelopmentC/ad2347/OxPhos_DAV"
N_CHUNKS=60

# Input pair list for this array run. Swap to tier1_pairs.csv (82 pairs, only
# 24 of which currently have PDB coordinates) or a custom subset as needed —
# 01_foldx_ddg.py requires the input CSV to carry a populated pdb_id column.
INPUT_CSV="${HPC_WORK}/results/mutagenesis/prioritized_pairs.csv"

module load miniconda/2024-02-20
eval "$(conda shell.bash hook)"
conda activate OXPHOS_DAV

mkdir -p "${HPC_WORK}/logs" "${HPC_WORK}/results/mutagenesis"

echo "Starting FoldX chunk ${SLURM_ARRAY_TASK_ID} of ${N_CHUNKS}"
echo "Host: $(hostname)  Date: $(date)"
echo "Input: ${INPUT_CSV}"

python "${HPC_HOME}/src/mutagenesis/01_foldx_ddg.py" \
    --input        "${INPUT_CSV}" \
    --chunk-idx    "${SLURM_ARRAY_TASK_ID}" \
    --n-chunks     "${N_CHUNKS}" \
    --work-dir     "${HPC_WORK}" \
    --skip-existing

echo "Finished FoldX chunk ${SLURM_ARRAY_TASK_ID}  $(date)"
