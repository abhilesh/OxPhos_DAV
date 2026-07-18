"""
src/phylo/merge_pagel_results.py

Merges HPC Pagel job results back into all_tested_pairs.csv and
recomputes compensatory_partners.csv with Pagel FDR column filled.

Run after:
  1. sbatch slurm_pagel_array.sh  (on HPC)
  2. rsync -av hpc:~/oxphos/pagel_results/ data/phylo/pagel_results/

Inputs:
  results/structural/all_tested_pairs.csv   -- existing output from partners script
  data/phylo/pagel_results/results_*.tsv    -- one file per SLURM array task
  data/phylo/pagel_jobs/all_records_index.json  -- identity_key → rec_idx (unused;
                                                   we join on identity_key directly)

Outputs:
  results/structural/all_tested_pairs.csv   -- updated with pagel_p and pagel_fdr
  results/structural/compensatory_partners.csv -- refiltered view

Run from project root inside the Docker container:
    python src/phylo/merge_pagel_results.py
"""

import csv
import sys
from pathlib import Path

import numpy as np

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "results" / "structural"
ALL_PAIRS   = RESULTS_DIR / "all_tested_pairs.csv"
PARTNERS    = RESULTS_DIR / "compensatory_partners.csv"
PAGEL_DIR   = ROOT / "data" / "phylo" / "pagel_results"

PAGEL_FDR_THRESH  = 0.10
BRANCH_FDR_THRESH = 0.10


# ── BH FDR ────────────────────────────────────────────────────────────────────

def bh_fdr(p_values: list[float]) -> list[float]:
    n = len(p_values)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: p_values[i])
    fdrs  = [0.0] * n
    min_so_far = 1.0
    for rank, i in enumerate(reversed(order), 1):
        raw = p_values[i] * n / (n - rank + 1)
        min_so_far = min(raw, min_so_far)
        fdrs[i] = min(min_so_far, 1.0)
    return fdrs


# ── Significance filter (mirrors partners script _get_sig) ────────────────────

def is_significant(row: dict) -> bool:
    if row.get("low_power", "").lower() in ("true", "1"):
        return False

    pagel_fdr  = _safe_float(row.get("pagel_fdr"))
    branch_fdr = _safe_float(row.get("branch_cooccur_fdr"))
    fisher_fdr = _safe_float(row.get("fisher_fdr"))

    has_pagel  = pagel_fdr  is not None
    has_branch = branch_fdr is not None

    if has_pagel or has_branch:
        return (has_pagel  and pagel_fdr  <= PAGEL_FDR_THRESH) or \
               (has_branch and branch_fdr <= BRANCH_FDR_THRESH)
    return fisher_fdr is not None and fisher_fdr <= PAGEL_FDR_THRESH


def _safe_float(val) -> float | None:
    if val is None or val in ("", "None", "NA", "nan"):
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    if not ALL_PAIRS.exists():
        print(f"ERROR: {ALL_PAIRS} not found — run 01_find_compensating_partners.py first.")
        sys.exit(1)

    if not PAGEL_DIR.exists():
        print(f"ERROR: {PAGEL_DIR} not found — copy HPC results there first.")
        sys.exit(1)

    result_files = sorted(PAGEL_DIR.glob("results_*.tsv"))
    if not result_files:
        print(f"ERROR: No results_*.tsv files in {PAGEL_DIR}.")
        sys.exit(1)

    print(f"Loading {len(result_files)} Pagel result files ...")

    # ── Load Pagel p-values: identity_key → p_value ───────────────────────────
    pagel_pvals: dict[str, float | None] = {}
    n_read = 0
    for rf in result_files:
        with open(rf) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) != 3 or parts[1] != "pagel_p":
                    continue
                identity_key = parts[0]
                val_str      = parts[2]
                p = None if val_str == "NA" else _safe_float(val_str)
                pagel_pvals[identity_key] = p
                n_read += 1

    n_non_na = sum(1 for p in pagel_pvals.values() if p is not None)
    print(f"  Pairs read: {n_read}  ({n_non_na} with finite p-value)")

    # ── Load all_tested_pairs.csv ─────────────────────────────────────────────
    print(f"Loading {ALL_PAIRS.name} ...")
    with open(ALL_PAIRS) as f:
        reader    = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows       = list(reader)
    print(f"  {len(rows)} pairs.")

    # ── Build identity keys and inject pagel_p ────────────────────────────────
    matched = 0
    for row in rows:
        key = (f"{row['ann_id']}|{row['contact_gene']}|"
               f"{row['contact_refseq_pos']}|{row['contact_alt_aa']}")
        if key in pagel_pvals:
            row["pagel_p"] = pagel_pvals[key]
            matched += 1

    print(f"  Matched: {matched} / {len(rows)}")
    not_matched = len(rows) - matched
    if not_matched:
        print(f"  Not matched (low_power or Fisher >= 0.20): {not_matched}")

    # ── Compute BH FDR globally across all pairs with finite pagel_p ─────────
    valid_idx = [i for i, r in enumerate(rows) if _safe_float(r.get("pagel_p")) is not None]
    if valid_idx:
        ps   = [_safe_float(rows[i]["pagel_p"]) for i in valid_idx]
        fdrs = bh_fdr(ps)
        for idx, fdr in zip(valid_idx, fdrs):
            rows[idx]["pagel_fdr"] = f"{fdr:.3e}"
        n_sig = sum(1 for fdr in fdrs if fdr <= PAGEL_FDR_THRESH)
        print(f"  Pagel FDR computed for {len(valid_idx)} pairs; "
              f"{n_sig} significant at FDR ≤ {PAGEL_FDR_THRESH}")
    else:
        print("  No finite pagel_p values — pagel_fdr remains empty.")

    # ── Format floats for writing ─────────────────────────────────────────────
    _fmt = lambda v: f"{v:.3e}" if isinstance(v, float) else ("" if v is None else str(v))

    def fmt_row(row: dict) -> dict:
        row = dict(row)
        for fld in ("fisher_p", "fisher_fdr", "pagel_p", "pagel_fdr",
                    "branch_cooccur_p", "branch_cooccur_fdr"):
            val = row.get(fld)
            if isinstance(val, float):
                row[fld] = f"{val:.3e}"
        return row

    # ── Rewrite all_tested_pairs.csv ─────────────────────────────────────────
    with open(ALL_PAIRS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(fmt_row(row))
    print(f"\nUpdated → {ALL_PAIRS}")

    # ── Rewrite compensatory_partners.csv ─────────────────────────────────────
    sig_rows = [r for r in rows if is_significant(r)]
    partner_fields = [f for f in fieldnames
                      if f not in ("n_dar_only_branches", "n_contact_only_branches")]
    with open(PARTNERS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=partner_fields, extrasaction="ignore")
        w.writeheader()
        for row in sig_rows:
            w.writerow(fmt_row(row))

    # ── Summary ───────────────────────────────────────────────────────────────
    pagel_sig  = sum(1 for r in sig_rows if _safe_float(r.get("pagel_fdr"))  is not None
                     and _safe_float(r.get("pagel_fdr")) <= PAGEL_FDR_THRESH)
    branch_sig = sum(1 for r in sig_rows if _safe_float(r.get("branch_cooccur_fdr")) is not None
                     and _safe_float(r.get("branch_cooccur_fdr")) <= BRANCH_FDR_THRESH)
    print(f"Updated → {PARTNERS}")
    print(f"\nSignificant pairs: {len(sig_rows)}")
    print(f"  Pagel FDR ≤ {PAGEL_FDR_THRESH}:        {pagel_sig}")
    print(f"  Branch co-occur FDR ≤ {BRANCH_FDR_THRESH}: {branch_sig}")


if __name__ == "__main__":
    main()
