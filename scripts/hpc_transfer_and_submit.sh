#!/bin/bash
# scripts/hpc_transfer_and_submit.sh
#
# Complete HPC workflow for OxPhos DAV analysis.
# Run each section sequentially from the project root on your local machine.
#
# HPC layout:
#   Scripts/code: /home/ad2347/OxPhos_DAV        (HPC_HOME)
#   Data/results: /work/WorkDevelopmentC/ad2347/OxPhos_DAV  (HPC_WORK)
# HPC alias:      hpc   (set in ~/.ssh/config, or replace with user@hostname)
#
# ─────────────────────────────────────────────────────────────────────────────
# OVERVIEW OF HPC-BOUND STAGES
# ─────────────────────────────────────────────────────────────────────────────
#
#  Stage A: Pagel's discrete test
#    Local:  python src/phylo/prepare_pagel_hpc.py
#    HPC:    sbatch /home/ad2347/OxPhos_DAV/scripts/slurm_pagel_array.sh
#    Local:  python src/phylo/merge_pagel_results.py
#
#  Stage B: DCA for all DAVs (per-gene plmDCA)
#    Local:  conda run -n oxphos_dav python src/mutagenesis/03a_prep_dca_jobs.py
#    HPC:    sbatch /home/ad2347/OxPhos_DAV/scripts/slurm_dca_array.sh
#    Local:  conda run -n oxphos_dav python src/mutagenesis/03_dca_all_davs.py
#              --aggregate results/mutagenesis/dca_gene_results
#
#  Stage C: Pyvolve conditional permissiveness
#    Local:  conda run -n oxphos_dav python src/phylo/02_phylogenetic_timing.py
#    HPC:    sbatch /home/ad2347/OxPhos_DAV/scripts/slurm_pyvolve_perm.sh
#    Local:  conda run -n oxphos_dav python src/phylo/05_merge_perm_chunks.py
#
# ─────────────────────────────────────────────────────────────────────────────

HPC_HOME="/home/ad2347/OxPhos_DAV"
HPC_WORK="/work/WorkDevelopmentC/ad2347/OxPhos_DAV"

# =============================================================================
# STAGE A — PAGEL'S DISCRETE TEST
# =============================================================================

# ── A1. Local prep ────────────────────────────────────────────────────────────
#
# Regenerate Pagel job manifests from the current classified parquet.
# (Old manifests in data/phylo/pagel_jobs/ were built from legacy JSON source.)
#
#   docker run --rm -v $(pwd):/app oxphos_dav_analysis \
#     conda run -n oxphos_dav python src/phylo/prepare_pagel_hpc.py

# ── A2. Transfer to HPC ───────────────────────────────────────────────────────
#
#   rsync -av data/phylo/pagel_jobs/       hpc:$HPC_WORK/data/phylo/pagel_jobs/
#   rsync -av data/phylo/ancestral_states/ hpc:$HPC_WORK/data/phylo/ancestral_states/
#   rsync -av data/reference/              hpc:$HPC_WORK/data/reference/
#   rsync -av src/                         hpc:$HPC_HOME/src/
#   rsync -av setup.py                     hpc:$HPC_HOME/
#   rsync -av scripts/slurm_pagel_array.sh hpc:$HPC_HOME/scripts/

# ── A3. Submit on HPC ─────────────────────────────────────────────────────────
#
#   ssh hpc "mkdir -p $HPC_WORK/logs $HPC_WORK/results/phylo/pagel_results && \
#     sbatch $HPC_HOME/scripts/slurm_pagel_array.sh"

# ── A4. Copy results back ─────────────────────────────────────────────────────
#
#   rsync -av hpc:$HPC_WORK/results/phylo/pagel_results/ \
#       data/phylo/pagel_results/

# ── A5. Merge locally ─────────────────────────────────────────────────────────
#
#   docker run --rm -v $(pwd):/app oxphos_dav_analysis \
#     conda run -n oxphos_dav python src/phylo/merge_pagel_results.py


# =============================================================================
# STAGE B — DCA FOR ALL DAVs
# =============================================================================

# ── B1. Local prep ────────────────────────────────────────────────────────────
#
# Writes sanitized per-gene FASTA files and generates slurm_dca_array.sh.
#
#   docker run --rm -v $(pwd):/app oxphos_dav_analysis \
#     conda run -n oxphos_dav python src/mutagenesis/03a_prep_dca_jobs.py

# ── B2. Transfer to HPC ───────────────────────────────────────────────────────
#
#   rsync -av data/dca_jobs/                     hpc:$HPC_WORK/data/dca_jobs/
#   rsync -av results/structural/dar_contacts_cbcb8A.csv \
#       hpc:$HPC_WORK/results/structural/
#   rsync -av src/                               hpc:$HPC_HOME/src/
#   rsync -av setup.py                           hpc:$HPC_HOME/
#   rsync -av scripts/slurm_dca_array.sh         hpc:$HPC_HOME/scripts/

# ── B3. Submit on HPC ─────────────────────────────────────────────────────────
#
#   ssh hpc "mkdir -p $HPC_WORK/logs $HPC_WORK/results/mutagenesis/dca_gene_results && \
#     sbatch $HPC_HOME/scripts/slurm_dca_array.sh"
#
# Monitor progress:
#   ssh hpc "squeue -u ad2347"
#   ssh hpc "ls $HPC_WORK/results/mutagenesis/dca_gene_results/ | wc -l"

# ── B4. Copy results back ─────────────────────────────────────────────────────
#
#   rsync -av hpc:$HPC_WORK/results/mutagenesis/dca_gene_results/ \
#       results/mutagenesis/dca_gene_results/

# ── B5. Aggregate locally ─────────────────────────────────────────────────────
#
#   docker run --rm -v $(pwd):/app oxphos_dav_analysis \
#     conda run -n oxphos_dav python src/mutagenesis/03_dca_all_davs.py \
#       --aggregate results/mutagenesis/dca_gene_results
#
# Outputs:
#   results/mutagenesis/dca_all_davs.csv
#   results/mutagenesis/dca_cdav_vs_udav_comparison.csv


# =============================================================================
# STAGE C — PYVOLVE CONDITIONAL PERMISSIVENESS
# =============================================================================

# ── C1. Local prep ────────────────────────────────────────────────────────────
#
# Runs phylogenetic timing analysis; produces timing_annotations.csv required
# by the pyvolve permutation stage.
#
#   docker run --rm -v $(pwd):/app oxphos_dav_analysis \
#     conda run -n oxphos_dav python src/phylo/02_phylogenetic_timing.py

# ── C2. Transfer to HPC ───────────────────────────────────────────────────────
#
#   rsync -av results/phylo/timing_annotations.csv  hpc:$HPC_WORK/results/phylo/
#   rsync -av data/phylo/ancestral_states/          hpc:$HPC_WORK/data/phylo/ancestral_states/
#   rsync -av data/reference/                        hpc:$HPC_WORK/data/reference/
#   rsync -av src/                                   hpc:$HPC_HOME/src/
#   rsync -av setup.py                               hpc:$HPC_HOME/
#   rsync -av scripts/slurm_pyvolve_perm.sh          hpc:$HPC_HOME/scripts/

# ── C3. Submit on HPC ─────────────────────────────────────────────────────────
#
#   ssh hpc "mkdir -p $HPC_WORK/logs $HPC_WORK/results/phylo && \
#     sbatch $HPC_HOME/scripts/slurm_pyvolve_perm.sh"

# ── C4. Copy results back ─────────────────────────────────────────────────────
#
#   rsync -av hpc:$HPC_WORK/results/phylo/conditional_permissiveness_chunk*.csv \
#       results/phylo/

# ── C5. Merge locally ─────────────────────────────────────────────────────────
#
#   docker run --rm -v $(pwd):/app oxphos_dav_analysis \
#     conda run -n oxphos_dav python src/phylo/05_merge_perm_chunks.py


# =============================================================================
# QUICK REFERENCE — FILE CHECKLIST
# =============================================================================
#
# Before each stage, verify these local files exist:
#
# Stage A (Pagel):
#   data/derived/classified/variants_master_classified.parquet  ← from classify stage
#   data/phylo/ancestral_states/{gene}.state                    ← from IQTree stage
#   data/reference/vert_mammal_tree.nwk                         ← downloaded
#
# Stage B (DCA):
#   results/structural/dar_contacts_cbcb8A.csv                  ← from structural stage
#   data/alignments/toga_hg38_aa/{gene}_aa_alignment.fasta      ← from align stage
#   data/alignments/mtdna_aa/{gene}_aa_alignment.fasta          ← from align stage
#
# Stage C (Pyvolve):
#   results/phylo/timing_annotations.csv                         ← from C1 local prep
#   data/phylo/ancestral_states/{gene}.state                    ← from IQTree stage
#
# =============================================================================
