#!/bin/bash
#SBATCH --job-name=plmdca_oxphos
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs/plmdca_%j.out
#SBATCH --error=/work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs/plmdca_%j.err
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=abhilesh7@gmail.com
#
# Single job for plmDCA (pseudolikelihood DCA) via 02_mi_analysis.py.
# Runs PlmDCA for all 693 compensatory pairs:
#   - 594 intraprotein pairs  → DI matrix cached per gene (~68 unique genes)
#   - 99 interprotein pairs   → concatenated MSA; M_eff check; mt-nuc flagged
#   - Top-50 by dca_di        → pyvolve WAG tree-null (n_sim=1000) [unless --skip-tree-null]
# Expected runtime: ~4–6 h with 8 CPUs (pydca uses OpenBLAS threads).
#
# ── SETUP (run locally first) ────────────────────────────────────────────────
#   rsync -av results/mutagenesis/prioritized_pairs.csv \
#       hpc:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/mutagenesis/
#   rsync -av results/structural/compensatory_partners.csv \
#       hpc:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/structural/
#   rsync -av results/structural/all_tested_pairs.csv \
#       hpc:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/structural/
#   rsync -av data/alignments/ \
#       hpc:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/data/alignments/
#   rsync -av data/phylo/iqtree_jobs/ \
#       hpc:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/data/phylo/iqtree_jobs/
#   rsync -av src/    hpc:/home/ad2347/OxPhos_DAV/src/
#   rsync -av setup.py hpc:/home/ad2347/OxPhos_DAV/
#   rsync -av scripts/slurm_plmdca.sh hpc:/home/ad2347/OxPhos_DAV/scripts/
#
# ── SUBMIT ───────────────────────────────────────────────────────────────────
#   mkdir -p /work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs
#   mkdir -p /work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/mutagenesis
#   sbatch /home/ad2347/OxPhos_DAV/scripts/slurm_plmdca.sh
#
#   To skip the pyvolve top-50 tree null (faster, ~2 h):
#   sbatch --export=ALL,SKIP_TREE_NULL=1 /home/ad2347/OxPhos_DAV/scripts/slurm_plmdca.sh
#
# ── COPY RESULTS BACK ────────────────────────────────────────────────────────
#   rsync -av \
#       hpc:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/mutagenesis/mi_scores.csv \
#       results/mutagenesis/
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

HPC_HOME="/home/ad2347/OxPhos_DAV"
HPC_WORK="/work/WorkDevelopmentC/ad2347/OxPhos_DAV"

module load miniconda/2024-02-20
eval "$(conda shell.bash hook)"
conda activate OXPHOS_DAV

# Let pydca / OpenBLAS use all allocated CPUs
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "${HPC_WORK}/logs" "${HPC_WORK}/results/mutagenesis"

echo "Starting plmDCA run  Host: $(hostname)  Date: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-8}   Mem: $(free -h | awk '/^Mem/{print $2}')"

EXTRA_FLAGS=""
if [[ "${SKIP_TREE_NULL:-0}" == "1" ]]; then
    EXTRA_FLAGS="--skip-tree-null"
    echo "Tree-null skipped (SKIP_TREE_NULL=1)"
fi

python "${HPC_HOME}/src/mutagenesis/02_mi_analysis.py" ${EXTRA_FLAGS}

echo "Finished plmDCA  $(date)"
