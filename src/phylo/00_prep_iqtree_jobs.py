"""
src/phylo/00_prep_iqtree_jobs.py

Prepares per-gene IQTree job directories for ancestral state reconstruction
on the HPC, and writes a SLURM array script.

Gene scope: ALL genes that have at least one cDAV variant in the classified
parquet (mtDNA + nucDNA). This is broader than the structural contact analysis —
it enables downstream intra-protein compensation scanning across all positions
in each protein.

Species scope: the 405-species cross-genome overlap (species with both nucDNA
TOGA and mtDNA data), represented in the pruned VertLife tree. Using a common
species set across all genes ensures:
  - Consistent phylogenetic background for Pagel tests
  - Branch co-occurrence is comparable across different gene pairs
  - Intra-protein compensation can be assessed in a single unified tree

Inputs (must exist before running):
  data/derived/classified/variants_master_classified.parquet
  data/phylo/species_tree/mammaltree_crossgenome.nwk  -- pruned to 405-species overlap
  data/phylo/cross_genome_species.txt                 -- species list
  data/alignments/toga_hg38_aa/                       -- nucDNA AA alignments
  data/alignments/mtdna_aa/                           -- mtDNA AA alignments

Run src/phylo/00a_download_vertlife_tree.py first to produce the pruned tree.

Outputs:
  data/phylo/iqtree_jobs/{gene}/
      {gene}.fasta         -- AA alignment filtered to cross-genome species
      {gene}_tree.nwk      -- cross-genome pruned tree (symlink or copy)
  scripts/slurm_iqtree_array.sh

Run from project root inside the Docker container:
    python src/phylo/00_prep_iqtree_jobs.py
"""

import sys
from pathlib import Path

import pandas as pd

from Bio import Phylo, SeqIO

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT               = Path(__file__).resolve().parents[2]
DATA_DIR           = ROOT / "data"
SCRIPTS_DIR        = ROOT / "scripts"

CLASSIFIED_PARQUET = ROOT / "data" / "derived" / "classified" / "variants_master_classified.parquet"
CG_TREE            = DATA_DIR / "phylo" / "species_tree" / "mammaltree_crossgenome.nwk"
CROSS_SPP    = DATA_DIR / "phylo" / "cross_genome_species.txt"
TOGA_AA_DIR  = DATA_DIR / "alignments" / "toga_hg38_aa"
MT_AA_DIR    = DATA_DIR / "alignments" / "mtdna_aa"
JOBS_DIR     = DATA_DIR / "phylo" / "iqtree_jobs"
SLURM_SCRIPT = SCRIPTS_DIR / "slurm_iqtree_array.sh"

_TOGA_TO_CANONICAL = {"COXFA4": "NDUFA4"}


# ── Gene discovery ─────────────────────────────────────────────────────────────

def genes_with_cdavs() -> dict[str, str]:
    """
    Returns {gene_symbol: genome} for all genes that have at least one cDAV.
    genome is "mtDNA" or "nucDNA".
    """
    df = pd.read_parquet(CLASSIFIED_PARQUET)
    cdavs = df[
        (df["classification_status"] == "classified") &
        (df["is_cdav_amino_acid"] == True)
    ].to_dict(orient="records")
    genes: dict[str, str] = {}
    for var in cdavs:
        gene = (
            var.get("interpreted_gene")
            or var.get("classification_gene")
            or str(var.get("locus", "")).split("/")[0]
        )
        if gene:
            genes[gene] = var.get("genome", "nucDNA")
    return genes


# ── Alignment helpers ─────────────────────────────────────────────────────────

def find_alignment(gene: str) -> tuple[Path | None, str]:
    """
    Returns (fasta_path, genome) for a gene, checking TOGA then mtDNA directories.
    Handles known gene aliases.
    """
    for aln_dir, genome in ((TOGA_AA_DIR, "nucDNA"), (MT_AA_DIR, "mtDNA")):
        for name in (gene, *[k for k, v in _TOGA_TO_CANONICAL.items() if v == gene]):
            fasta = aln_dir / f"{name}_aa_alignment.fasta"
            if fasta.exists():
                return fasta, genome
    return None, ""


def write_filtered_fasta(fasta_in: Path, keep_spp: set[str], out_path: Path) -> int:
    """Write alignment FASTA keeping only cross-genome species + Homo_sapiens reference.

    IQ-TREE rejects:
    - Sequence names with special characters (|, spaces, colons, parentheses) —
      headers are sanitized to just the species name (part before first '|').
    - Duplicate sequence names — TOGA alignments may have multiple transcripts per
      species. When duplicates exist, keep the sequence with the fewest gap characters.
    """
    # Collect best sequence per species (fewest gaps wins)
    best: dict[str, object] = {}
    for rec in SeqIO.parse(fasta_in, "fasta"):
        sp = rec.id.split("|")[0]
        if sp != "Homo_sapiens" and sp not in keep_spp:
            continue
        rec.id          = sp
        rec.name        = sp
        rec.description = ""
        if sp not in best:
            best[sp] = rec
        else:
            # Prefer the sequence with fewer gap/ambiguous characters
            prev_gaps = str(best[sp].seq).count("-") + str(best[sp].seq).count("X")
            this_gaps = str(rec.seq).count("-") + str(rec.seq).count("X")
            if this_gaps < prev_gaps:
                best[sp] = rec

    with open(out_path, "w") as f:
        for rec in best.values():
            SeqIO.write(rec, f, "fasta")
    return len(best)


# ── SLURM script ───────────────────────────────────────────────────────────────

def write_slurm_script(job_dirs: list[Path], out_path: Path):
    gene_list = " ".join(d.name for d in sorted(job_dirs))
    n = len(job_dirs)
    script = f"""\
#!/bin/bash
#SBATCH --job-name=iqtree_oxphos
#SBATCH --array=1-{n}
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=logs/iqtree_%A_%a.out
#SBATCH --error=logs/iqtree_%A_%a.err

# Auto-generated by src/phylo/00_prep_iqtree_jobs.py — do not edit manually.
#
# Transfer inputs:
#   rsync -av data/phylo/iqtree_jobs/ hpc:~/oxphos/iqtree_jobs/
# Submit:
#   sbatch scripts/slurm_iqtree_array.sh
# Copy results back:
#   rsync -av hpc:~/oxphos/ancestral_states/ data/phylo/ancestral_states/

set -euo pipefail

# Activate IQTree via conda (adjust env name as needed on your HPC)
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate iqtree 2>/dev/null || true

GENES=({gene_list})
GENE="${{GENES[${{SLURM_ARRAY_TASK_ID}}-1]}}"

JOB_DIR="$HOME/oxphos/iqtree_jobs/$GENE"
OUT_DIR="$HOME/oxphos/ancestral_states/$GENE"
mkdir -p "$OUT_DIR" logs

cd "$JOB_DIR" || {{ echo "Cannot cd to $JOB_DIR"; exit 1; }}

iqtree \\
    -s "${{GENE}}.fasta" \\
    -te "${{GENE}}_tree.nwk" \\
    -m TEST \\
    --ancestral \\
    --prefix "$OUT_DIR/$GENE" \\
    -T 4 \\
    --quiet

echo "Done: $GENE"
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(script)
    out_path.chmod(0o755)
    print(f"SLURM script ({n} jobs) → {out_path}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("Preparing IQTree HPC jobs (scope: all genes with cDAVs)\n")

    # ── Check prerequisites ────────────────────────────────────────────────────
    if not CLASSIFIED_PARQUET.exists():
        print(
            f"ERROR: Classified parquet not found at {CLASSIFIED_PARQUET}\n"
            "Run src/classify/00_classify_DAV.py first."
        )
        sys.exit(1)

    if not CG_TREE.exists():
        print(
            f"ERROR: Cross-genome pruned tree not found at {CG_TREE}\n"
            "Run src/phylo/00a_download_vertlife_tree.py first."
        )
        sys.exit(1)

    if not CROSS_SPP.exists():
        print(f"ERROR: {CROSS_SPP} not found — run 00a_download_vertlife_tree.py first.")
        sys.exit(1)

    # ── Load cross-genome species (TOGA names, TaxID-matched) ─────────────────
    cross_spp = set(CROSS_SPP.read_text().splitlines())
    cross_spp.discard("")   # remove blank lines
    print(f"Cross-genome species (TaxID-matched): {len(cross_spp)}")

    # Load pruned tree tip names (may be subset of cross_spp if some absent from VertLife)
    cg_tree = Phylo.read(str(CG_TREE), "newick")
    tree_tips = {c.name for c in cg_tree.get_terminals()}
    print(f"Tips in cross-genome tree: {len(tree_tips)}\n")

    # ── Discover genes with cDAVs ─────────────────────────────────────────────
    gene_genome = genes_with_cdavs()
    print(f"Genes with cDAVs: {len(gene_genome)}")
    mt_genes  = [g for g, gn in gene_genome.items() if gn == "mtDNA"]
    nuc_genes = [g for g, gn in gene_genome.items() if gn == "nucDNA"]
    print(f"  mtDNA : {len(mt_genes)}")
    print(f"  nucDNA: {len(nuc_genes)}\n")

    # ── Build job directories ─────────────────────────────────────────────────
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    job_dirs: list[Path] = []
    skipped: list[str]   = []

    for gene in sorted(gene_genome):
        fasta_in, genome = find_alignment(gene)
        if fasta_in is None:
            print(f"  [MISS]  {gene}: no alignment found")
            skipped.append(gene)
            continue

        # Filter to cross-genome species present in alignment
        spp_in_aln = {
            rec.id.split("|")[0]
            for rec in SeqIO.parse(fasta_in, "fasta")
            if not rec.id.startswith("Homo_sapiens")
        }
        keep = spp_in_aln & tree_tips   # only species in pruned tree
        n_keep = len(keep)

        if n_keep < 20:
            print(f"  [SKIP]  {gene}: only {n_keep} cross-genome species in alignment — too few")
            skipped.append(gene)
            continue

        # Always include Homo_sapiens as reference in both tree and FASTA
        keep_with_human = keep | {"Homo_sapiens"}

        # Write job directory
        job_dir  = JOBS_DIR / gene
        job_dir.mkdir(exist_ok=True)

        fasta_out = job_dir / f"{gene}.fasta"
        tree_out  = job_dir / f"{gene}_tree.nwk"

        # Filtered FASTA (cross-genome species + Homo_sapiens reference)
        n_written = write_filtered_fasta(fasta_in, keep, fasta_out)

        # Gene-specific pruned tree: prune cross-genome tree to species in this alignment
        # keep_with_human ensures Homo_sapiens is retained in the tree
        gene_tree = Phylo.read(str(CG_TREE), "newick")
        all_gene_tips = {c.name for c in gene_tree.get_terminals()}
        for tip in all_gene_tips - keep_with_human:
            gene_tree.prune(tip)
        Phylo.write(gene_tree, str(tree_out), "newick")

        n_tips = len([c for c in gene_tree.get_terminals()])
        print(f"  [OK]    {gene} ({genome}): {n_written} seqs, {n_tips} tree tips")
        job_dirs.append(job_dir)

    print(f"\n{'='*55}")
    print(f"Job directories created : {len(job_dirs)}")
    if skipped:
        print(f"Skipped ({len(skipped)})           : {', '.join(skipped)}")

    if job_dirs:
        write_slurm_script(job_dirs, SLURM_SCRIPT)
        print(f"\nNext steps:")
        print(f"  1. rsync -av data/phylo/iqtree_jobs/ hpc:~/oxphos/iqtree_jobs/")
        print(f"  2. sbatch scripts/slurm_iqtree_array.sh")
        print(f"  3. rsync -av hpc:~/oxphos/ancestral_states/ data/phylo/ancestral_states/")
        print(f"  4. python src/phylo/01_parse_ancestral_states.py")


if __name__ == "__main__":
    main()