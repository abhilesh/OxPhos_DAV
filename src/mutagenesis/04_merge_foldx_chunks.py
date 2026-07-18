#!/usr/bin/env python3
"""
src/mutagenesis/04_merge_foldx_chunks.py

Merge per-chunk FoldX ΔΔG results produced by a SLURM job array
(01_foldx_ddg.py --n-chunks N --chunk-idx i) into a single
results/mutagenesis/foldx_ddg.csv, matching the schema/behaviour of a
non-chunked run.

If results/mutagenesis/foldx_ddg.csv already exists (e.g. from an earlier
local 50-pair run), chunk rows are appended/merged on top of it — existing
rows for the same physical pair (dar_gene, dar_aa_coord, dar_alt_aa,
contact_gene, contact_refseq_pos, contact_alt_aa) are overwritten by the
new chunk result (the array run is assumed authoritative for any pair it
touched), everything else is kept as-is.

Usage (run locally after rsync-ing chunk files back from HPC):
    python src/mutagenesis/04_merge_foldx_chunks.py [--input <pairs_csv used for the array>]
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT    = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "mutagenesis"

CHUNK_PATTERN = "foldx_ddg_chunk*.csv"
FINAL_OUT     = OUT_DIR / "foldx_ddg.csv"

PHYS_KEY = [
    "dar_gene", "dar_aa_coord", "dar_alt_aa",
    "contact_gene", "contact_refseq_pos", "contact_alt_aa",
]


def _print_summary(df: pd.DataFrame) -> None:
    print("\n── FoldX status breakdown (merged) ──────────────────────────")
    print(df["foldx_status"].value_counts(dropna=False).to_string())
    if "ddg_rescue_stab" in df.columns and "ddg_stab_dar" in df.columns:
        n_rescue_10 = (df["ddg_rescue_stab"] > 1.0).sum()
        n_dar_10    = (df["ddg_stab_dar"] > 1.0).sum()
        n_both_10   = ((df["ddg_rescue_stab"] > 1.0) & (df["ddg_stab_dar"] > 1.0)).sum()
        print(f"\n  ΔΔG_DAR > 1.0 kcal/mol:                {n_dar_10}")
        print(f"  ΔΔG_rescue > 1.0 kcal/mol:             {n_rescue_10}")
        print(f"  Both > 1.0 (foldx_tier='strong'):      {n_both_10}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--chunk-dir", type=Path, default=OUT_DIR,
        help="Directory containing foldx_ddg_chunk*.csv (default: results/mutagenesis/)",
    )
    parser.add_argument(
        "--keep-chunks", action="store_true",
        help="Do not delete chunk files after a successful merge",
    )
    args = parser.parse_args()

    chunk_files = sorted(args.chunk_dir.glob(CHUNK_PATTERN))
    if not chunk_files:
        sys.exit(f"ERROR: No chunk files found matching {args.chunk_dir}/{CHUNK_PATTERN}\n"
                  "  Run the SLURM array first (scripts/slurm_foldx_array.sh) and rsync results back.")

    print(f"Merging {len(chunk_files)} chunk files...")
    chunk_dfs = [pd.read_csv(cf, dtype={"ann_id": str}) for cf in chunk_files]
    new_df = pd.concat(chunk_dfs, ignore_index=True)
    print(f"  New rows from chunks: {len(new_df)}")

    # Drop the bookkeeping column used only to keep foldx_work/pairNNNN_* unique
    # across array tasks — not part of the standard foldx_ddg.csv schema.
    new_df = new_df.drop(columns=["orig_idx"], errors="ignore")

    missing_key = [c for c in PHYS_KEY if c not in new_df.columns]
    if missing_key:
        sys.exit(f"ERROR: chunk files missing join-key columns: {missing_key}")

    if FINAL_OUT.exists():
        old_df = pd.read_csv(FINAL_OUT, dtype={"ann_id": str})
        print(f"  Existing foldx_ddg.csv: {len(old_df)} rows (will be updated, not replaced)")
        for col in PHYS_KEY:
            old_df[col] = old_df[col].astype(str)
            new_df[col] = new_df[col].astype(str)
        old_keys = old_df[PHYS_KEY].agg("|".join, axis=1)
        new_keys = new_df[PHYS_KEY].agg("|".join, axis=1)
        old_df = old_df[~old_keys.isin(set(new_keys))]
        merged = pd.concat([old_df, new_df], ignore_index=True)
    else:
        merged = new_df

    merged.to_csv(FINAL_OUT, index=False)
    print(f"\n  Total rows in merged foldx_ddg.csv: {len(merged)}")
    print(f"  Saved: {FINAL_OUT}")

    _print_summary(merged)

    if not args.keep_chunks:
        for cf in chunk_files:
            cf.unlink()
        print(f"\nChunk files deleted ({len(chunk_files)} files).")
    else:
        print(f"\nChunk files kept (--keep-chunks): {len(chunk_files)} files.")

    print("\nNext step: re-run src/mutagenesis/03_compile_targets.py to fold these")
    print("FoldX results into results/mutagenesis/final_targets.csv.")


if __name__ == "__main__":
    main()
