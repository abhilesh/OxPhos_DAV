#!/bin/bash
#SBATCH --job-name=merge_perm
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
#SBATCH --output=/work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs/merge_perm_%j.out
#SBATCH --error=/work/WorkDevelopmentC/ad2347/OxPhos_DAV/logs/merge_perm_%j.err
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=abhilesh7@gmail.com
#
# Post-array merge: collect all conditional_permissiveness_chunk*.csv files
# into a single conditional_permissiveness.csv and print the pre-registered
# decision summary.
#
# Submit AFTER the pyvolve array completes (use --dependency):
#   ARRAY_JOB=$(sbatch --parsable /home/ad2347/OxPhos_DAV/scripts/slurm_pyvolve_perm.sh)
#   sbatch --dependency=afterok:${ARRAY_JOB} /home/ad2347/OxPhos_DAV/scripts/slurm_merge_perm.sh
#
# Or submit manually after confirming all array tasks succeeded:
#   sbatch /home/ad2347/OxPhos_DAV/scripts/slurm_merge_perm.sh
#
# ── COPY RESULTS BACK (after merge job completes) ────────────────────────────
#   rsync -av \
#       hpc:/work/WorkDevelopmentC/ad2347/OxPhos_DAV/results/phylo/conditional_permissiveness.csv \
#       results/phylo/
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

HPC_HOME="/home/ad2347/OxPhos_DAV"
HPC_WORK="/work/WorkDevelopmentC/ad2347/OxPhos_DAV"

module load miniconda/2024-02-20
eval "$(conda shell.bash hook)"
conda activate OXPHOS_DAV

echo "Merging pyvolve permutation chunks  Host: $(hostname)  Date: $(date)"

python "${HPC_HOME}/src/phylo/05_merge_perm_chunks.py" \
    --work-dir "${HPC_WORK}"

echo "Merge complete  $(date)"
