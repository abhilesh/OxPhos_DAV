#!/usr/bin/env python3
"""
src/mutagenesis/03b_dca_gene_worker.py

HPC per-gene plmDCA worker. Called once per SLURM array task.
Runs full plmDCA on a sanitized alignment and writes results as .npz.

Usage (auto-invoked by scripts/slurm_dca_array.sh):
  python src/mutagenesis/03b_dca_gene_worker.py GENE ALN_FASTA OUT_NPZ

Arguments:
  GENE      — gene symbol (for logging only)
  ALN_FASTA — path to sanitized alignment (stop codons already replaced by 03a)
  OUT_NPZ   — output path for numpy .npz with arrays col_i, col_j, di_score
"""

import sys
import time
from pathlib import Path

import numpy as np
from Bio import SeqIO


def _compute_meff(seqs: list[str], seqid: float = 0.8) -> float:
    """Vectorized Meff; O(n²) but fast for ≤1000 sequences."""
    n = len(seqs)
    if n == 0:
        return 0.0
    gap = ord("-")
    mat = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(n, -1)
    sim_count = np.ones(n, dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            nongap = (mat[i] != gap) & (mat[j] != gap)
            n_ng = int(nongap.sum())
            if n_ng > 0 and int(((mat[i] == mat[j]) & nongap).sum()) / n_ng > seqid:
                sim_count[i] += 1.0
                sim_count[j] += 1.0
    return float((1.0 / sim_count).sum())


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} GENE ALN_FASTA OUT_NPZ")

    gene, aln_fasta, out_npz = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])

    if not aln_fasta.exists():
        sys.exit(f"ERROR: {aln_fasta} not found")

    recs = list(SeqIO.parse(str(aln_fasta), "fasta"))
    if len(recs) < 10:
        sys.exit(f"ERROR: {gene} has only {len(recs)} sequences — skipping")

    seqs = [str(r.seq).upper() for r in recs]
    n, L = len(seqs), len(seqs[0])
    print(f"Gene: {gene}  seqs: {n}  length: {L}", flush=True)

    meff = _compute_meff(seqs)
    print(f"Meff: {meff:.1f}", flush=True)

    # Import here so import errors surface clearly in SLURM logs
    from pydca.plmdca.plmdca import PlmDCA

    t0 = time.time()
    plm = PlmDCA(
        biomolecule="protein",
        msa_file=str(aln_fasta),
        seqid=0.8,
        lambda_h=0.01,
        lambda_J=0.05,
    )
    raw = plm.compute_sorted_DI()
    elapsed = time.time() - t0
    print(f"plmDCA done in {elapsed:.1f}s  pairs: {len(raw)}", flush=True)

    # raw is list of ((col_i, col_j), score)
    col_is   = np.empty(len(raw), dtype=np.int32)
    col_js   = np.empty(len(raw), dtype=np.int32)
    di_scores = np.empty(len(raw), dtype=np.float64)
    for k, entry in enumerate(raw):
        ci, cj = int(entry[0][0]), int(entry[0][1])
        if ci > cj:
            ci, cj = cj, ci
        col_is[k]    = ci
        col_js[k]    = cj
        di_scores[k] = float(entry[1])

    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        str(out_npz),
        col_i=col_is,
        col_j=col_js,
        di_score=di_scores,
        meff=np.array([meff]),
        n_seqs=np.array([n]),
        aln_len=np.array([L]),
    )
    print(f"Saved → {out_npz}", flush=True)


if __name__ == "__main__":
    main()
