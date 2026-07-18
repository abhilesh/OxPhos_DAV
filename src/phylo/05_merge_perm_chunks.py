#!/usr/bin/env python3
"""
src/phylo/05_merge_perm_chunks.py

Merge per-chunk conditional permissiveness results produced by a SLURM
job array into a single results/phylo/conditional_permissiveness.csv.

Prints the pre-registered decision summary after merging.

Usage (run after all SLURM tasks complete):
    python src/phylo/05_merge_perm_chunks.py
"""
import csv
import sys
from pathlib import Path

import numpy as np

ROOT    = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "phylo"

CHUNK_PATTERN = "conditional_permissiveness_chunk*.csv"
FINAL_OUT     = OUT_DIR / "conditional_permissiveness.csv"


def summarise(rows: list[dict]) -> None:
    from scipy.stats import binomtest

    valid = [r for r in rows if r.get("perm_p") not in ("", None)]
    if not valid:
        print("No permutation results with perm_p found.")
        return

    perm_ps = np.array([float(r["perm_p"]) for r in valid])
    obs_ors = [float(r["observed_or"]) for r in valid
               if r.get("observed_or") not in ("", None)]

    n_total = len(valid)
    n_sig   = int((perm_ps < 0.05).sum())
    frac_sig = n_sig / n_total

    print(f"\nPermutation summary ({n_total} cDAV contact pairs tested):")
    print(f"  Significant (perm_p < 0.05):  {n_sig}/{n_total} = {100*frac_sig:.1f}%")
    print(f"  Significant (perm_p < 0.01):  {int((perm_ps<0.01).sum())}/{n_total} = {100*np.mean(perm_ps<0.01):.1f}%")
    print(f"  Median observed OR:            {np.median(obs_ors):.2f}")
    print()
    # Compare observed sig rate to 5% null (what permutation test calibrates against)
    bt = binomtest(n_sig, n_total, 0.05, alternative="greater")
    print(f"  Binomial test vs 5% null: {n_sig}/{n_total} sig → p = {bt.pvalue:.2e}")
    print(f"  (perm_p is already calibrated against neutral evolution; the 5% null")
    print(f"   represents the background rate for contacts without co-evolutionary signal)")

    for asr in ("high", "low", "root"):
        sub = [r for r in valid if r.get("asr_confidence") == asr]
        if not sub:
            continue
        sub_ps  = np.array([float(r["perm_p"]) for r in sub])
        sub_ors = [float(r["observed_or"]) for r in sub if r.get("observed_or")]
        ns = int((sub_ps < 0.05).sum())
        print(f"\n  ASR confidence = {asr} (n={len(sub)}):")
        print(f"    perm_p < 0.05:  {ns}/{len(sub)} = {100*ns/len(sub):.1f}%")
        print(f"    Median OR = {np.median(sub_ors):.2f}")


def main() -> None:
    chunk_files = sorted(OUT_DIR.glob(CHUNK_PATTERN))
    if not chunk_files:
        sys.exit(f"ERROR: No chunk files found matching {OUT_DIR}/{CHUNK_PATTERN}\n"
                 "  Run the SLURM array first.")

    print(f"Merging {len(chunk_files)} chunk files...")

    all_rows: list[dict] = []
    fieldnames: list[str] = []

    for cf in chunk_files:
        with open(cf) as f:
            reader = csv.DictReader(f)
            if not fieldnames:
                fieldnames = reader.fieldnames or []
            all_rows.extend(reader)

    print(f"  Total rows: {len(all_rows)}")

    with open(FINAL_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_rows)
    print(f"  Saved: {FINAL_OUT}")

    summarise(all_rows)

    # Clean up chunk files
    for cf in chunk_files:
        cf.unlink()
    print(f"\nChunk files deleted.")


if __name__ == "__main__":
    main()
