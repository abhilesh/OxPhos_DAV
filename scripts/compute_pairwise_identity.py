"""
Compute pairwise sequence identity distributions for all OXPHOS alignment files.

For each gene, computes:
- Pairwise identity distribution (sampled for large genes to stay fast)
- Distribution statistics (mean, median, P5, P25, P75, P95)
- Meff at multiple seqid thresholds (0.6, 0.7, 0.8, 0.9, 0.95)

Sampling: for genes with > MAX_SAMPLE sequences, draw MAX_SAMPLE sequences
uniformly at random (seed=42) for the pairwise distribution. Meff is computed
on the full alignment.

Outputs:
  results/mutagenesis/pairwise_identity_summary.csv

Usage (inside Docker):
    python scripts/compute_pairwise_identity.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
ALN_NUC = ROOT / "data" / "alignments" / "toga_hg38_aa"
ALN_MT  = ROOT / "data" / "alignments" / "mtdna_aa"
OUT_DIR = ROOT / "results" / "mutagenesis"

SEQID_THRESHOLDS = [0.60, 0.70, 0.80, 0.90, 0.95]
GAP_BYTES = frozenset(ord(c) for c in "-X*")
MAX_SAMPLE = 1000   # max seqs for pairwise dist estimation
MEFF_BLOCK  = 100   # block size for Meff computation


def load_mat(fasta_path: Path) -> np.ndarray:
    """Load FASTA → (n_seqs, L) uint8 matrix."""
    seqs = [np.frombuffer(str(r.seq).encode(), dtype=np.uint8)
            for r in SeqIO.parse(str(fasta_path), "fasta")]
    return np.vstack(seqs)


def non_gap_mask(mat: np.ndarray) -> np.ndarray:
    """Return (n, L) bool mask: True = non-gap character."""
    mask = np.ones(mat.shape, dtype=bool)
    for b in GAP_BYTES:
        mask &= (mat != b)
    return mask


def pairwise_identity_sample(mat: np.ndarray, ng: np.ndarray,
                              rng: np.random.Generator) -> np.ndarray:
    """
    Compute all pairwise identities among a (possibly sub-sampled) sequence set.
    Uses numpy broadcasting block-by-block; returns flat array of upper-triangle values.
    """
    n = len(mat)
    if n > MAX_SAMPLE:
        idx = rng.choice(n, MAX_SAMPLE, replace=False)
        mat = mat[idx]
        ng  = ng[idx]
        n   = MAX_SAMPLE

    BLOCK = 256
    idents = []
    for i in range(0, n, BLOCK):
        i_end = min(i + BLOCK, n)
        for j in range(i, n, BLOCK):
            j_end = min(j + BLOCK, n)
            # shapes: (bi, 1, L), (1, bj, L) → (bi, bj, L)
            and_mask = ng[i:i_end, None, :] & ng[j:j_end, None, :].swapaxes(0,1)
            # wait — ng[j:j_end] is (bj, L), need (1, bj, L)
            and_mask = ng[i:i_end][:, None, :] & ng[j:j_end][None, :, :]  # (bi, bj, L)
            match    = (mat[i:i_end][:, None, :] == mat[j:j_end][None, :, :]) & and_mask
            denom = and_mask.sum(axis=2).astype(np.float32)
            numer = match.sum(axis=2).astype(np.float32)
            with np.errstate(invalid="ignore"):
                pid = np.where(denom > 0, numer / denom, 0.0)  # (bi, bj)
            # Extract upper triangle relative to global indices
            bi_range = np.arange(i, i_end)
            bj_range = np.arange(j, j_end)
            for li, gi in enumerate(bi_range):
                start = max(0, gi + 1 - j)  # first bj index where gj > gi
                if start < len(bj_range):
                    idents.append(pid[li, start:])
    return np.concatenate(idents).astype(np.float32) if idents else np.array([], dtype=np.float32)


def compute_meff(mat: np.ndarray, ng: np.ndarray, threshold: float) -> float:
    """
    Meff at given seqid threshold (AND criterion, full alignment).
    weights[i] = number of sequences with pid(i,j) >= threshold (including self).
    Meff = sum(1/weights).
    """
    n = len(mat)
    weights = np.zeros(n, dtype=np.float64)
    for i in range(0, n, MEFF_BLOCK):
        i_end = min(i + MEFF_BLOCK, n)
        and_mask = ng[i:i_end][:, None, :] & ng[None, :, :]   # (bi, n, L)
        match    = (mat[i:i_end][:, None, :] == mat[None, :, :]) & and_mask
        denom = and_mask.sum(axis=2).astype(np.float64)
        numer = match.sum(axis=2).astype(np.float64)
        with np.errstate(invalid="ignore"):
            pid = np.where(denom > 0, numer / denom, 0.0)     # (bi, n)
        weights[i:i_end] = (pid >= threshold).sum(axis=1)
    weights = np.maximum(weights, 1.0)
    return float((1.0 / weights).sum())


def process_gene(fasta_path: Path, genome: str, rng: np.random.Generator) -> dict:
    gene = fasta_path.stem.replace("_aa_alignment", "")
    mat  = load_mat(fasta_path)
    n_seqs, L = mat.shape
    print(f"  {gene:20s} ({genome}) n={n_seqs:5d}  L={L:4d}", end="  ", flush=True)

    if n_seqs < 2:
        return {"gene": gene, "genome": genome, "n_seqs": n_seqs, "aln_length": L,
                "mean_pid": np.nan, "median_pid": np.nan,
                "p5_pid": np.nan, "p25_pid": np.nan, "p75_pid": np.nan, "p95_pid": np.nan,
                "pid_sampled": False,
                **{f"meff_{int(t*100)}": np.nan for t in SEQID_THRESHOLDS}}

    ng = non_gap_mask(mat)
    pids = pairwise_identity_sample(mat, ng, rng)

    meff_vals = {}
    for t in SEQID_THRESHOLDS:
        meff_vals[f"meff_{int(t*100)}"] = compute_meff(mat, ng, t)

    row = {
        "gene": gene,
        "genome": genome,
        "n_seqs": n_seqs,
        "aln_length": L,
        "mean_pid":   float(pids.mean()),
        "median_pid": float(np.median(pids)),
        "p5_pid":     float(np.percentile(pids, 5)),
        "p25_pid":    float(np.percentile(pids, 25)),
        "p75_pid":    float(np.percentile(pids, 75)),
        "p95_pid":    float(np.percentile(pids, 95)),
        "pid_sampled": n_seqs > MAX_SAMPLE,
        **meff_vals,
    }
    print(f"median_pid={row['median_pid']:.3f}  "
          f"meff_60={row['meff_60']:.1f}  meff_80={row['meff_80']:.1f}  meff_90={row['meff_90']:.1f}")
    return row


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    rows = []

    print("=== nucDNA alignments ===")
    for fasta in sorted(ALN_NUC.glob("*_aa_alignment.fasta")):
        rows.append(process_gene(fasta, "nucDNA", rng))

    print("\n=== mtDNA alignments ===")
    for fasta in sorted(ALN_MT.glob("*_aa_alignment.fasta")):
        rows.append(process_gene(fasta, "mtDNA", rng))

    df = pd.DataFrame(rows)
    out_csv = OUT_DIR / "pairwise_identity_summary.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {len(df)} rows → {out_csv.relative_to(ROOT)}")

    for genome, grp in df.groupby("genome"):
        print(f"\n── {genome} (n={len(grp)} genes) ──")
        print(f"  median_pid : {grp['median_pid'].mean():.3f} ± {grp['median_pid'].std():.3f}"
              f"  [{grp['median_pid'].min():.3f}, {grp['median_pid'].max():.3f}]")
        for t in SEQID_THRESHOLDS:
            col = f"meff_{int(t*100)}"
            print(f"  {col}     : {grp[col].mean():.1f} ± {grp[col].std():.1f}"
                  f"  [{grp[col].min():.1f}, {grp[col].max():.1f}]")

    # Genes where meff_80 < 10
    low = df[df["meff_80"] < 10].sort_values("meff_80")
    print(f"\nGenes with Meff_80 < 10 (DCA unreliable): {len(low)}")
    if len(low):
        cols = ["gene", "genome", "n_seqs", "median_pid", "meff_60", "meff_80", "meff_90"]
        print(low[cols].to_string(index=False))


if __name__ == "__main__":
    main()
