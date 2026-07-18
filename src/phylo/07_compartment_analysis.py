"""
Genomic Compartment Analysis — Phase 1 of downstream analytical pipeline.

Characterises the cDAV/uDAV compensation landscape separately for:
  - intra-mtDNA (mt-mt)
  - intra-nucDNA (nuc-nuc)
  - mito-nuclear (mt-nuc)

Uses only locally available results; no new HPC runs required.

Inputs:
  results/structural/dar_contacts_cbcb8A.csv     — all variant-contact pairs; contact_type, is_cdav_amino_acid
  results/structural/compensatory_partners.csv   — pagel_fdr, branch_cooccur_fdr
  results/mutagenesis/dca_all_davs.csv           — mi_apc, dca_di_percentile, is_cdav_amino_acid
  results/phylo/conditional_permissiveness.csv   — perm_p, observed_or

Outputs:
  results/phylo/compartment_analysis.csv         — per-compartment × test summary (one row per stratum)
  results/phylo/mito_nuclear_detail.csv          — one row per mito-nuclear pair with all test results

Usage (inside Docker):
  python src/phylo/07_compartment_analysis.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu, fisher_exact

ROOT = Path(__file__).resolve().parents[2]
CONTACTS  = ROOT / "results" / "structural" / "dar_contacts_cbcb8A.csv"
COMP_PART = ROOT / "results" / "structural" / "compensatory_partners.csv"
DCA       = ROOT / "results" / "mutagenesis" / "dca_all_davs.csv"
PYVOLVE   = ROOT / "results" / "phylo" / "conditional_permissiveness.csv"
OUT_DIR   = ROOT / "results" / "phylo"

JOIN = ["dar_gene", "dar_aa_coord", "contact_gene", "contact_refseq_pos"]
PAGEL_ALPHA  = 0.10
BRANCH_ALPHA = 0.10
PYVOLVE_ALPHA = 0.05
APC_PERCENTILE_THRESHOLD = 75.0   # top quartile within gene


def load_universe() -> pd.DataFrame:
    """
    Build the universe of unique position-pairs from dar_contacts_cbcb8A.csv.
    Contacts file uses dar_locus; rename to dar_gene for consistent join keys.
    One row per (dar_gene, dar_aa_coord, contact_gene, contact_refseq_pos, contact_type).
    is_cdav = True if ANY source row has is_cdav_amino_acid == True.
    """
    raw = pd.read_csv(CONTACTS)
    raw = raw.rename(columns={"dar_locus": "dar_gene"})

    universe = (
        raw.groupby(JOIN + ["contact_type", "dar_genome"], as_index=False)
        .agg(
            is_cdav=("is_cdav_amino_acid", "any"),
            n_contact_classes=("contact_class", "nunique"),
            contact_classes=("contact_class", lambda s: ",".join(sorted(s.unique()))),
        )
    )
    return universe


def add_pagel(universe: pd.DataFrame) -> pd.DataFrame:
    """Join best Pagel/branch FDR per pair from compensatory_partners.csv."""
    cp = pd.read_csv(COMP_PART)
    best = (
        cp.groupby(JOIN, as_index=False)
        .agg(pagel_fdr=("pagel_fdr", "min"), branch_fdr=("branch_cooccur_fdr", "min"))
    )
    return universe.merge(best, on=JOIN, how="left")


def add_pyvolve(universe: pd.DataFrame) -> pd.DataFrame:
    """Join best (min perm_p) pyvolve result per pair."""
    pv = pd.read_csv(PYVOLVE)
    best = (
        pv.sort_values("perm_p")
        .groupby(JOIN, as_index=False)
        .first()
        .rename(columns={"observed_or": "pyvolve_or", "perm_p": "pyvolve_p"})
    )
    return universe.merge(
        best[JOIN + ["pyvolve_or", "pyvolve_p", "perm_method"]], on=JOIN, how="left"
    )


def add_apc_mi(universe: pd.DataFrame) -> pd.DataFrame:
    """
    Join within-gene APC-MI percentile from dca_all_davs.csv.
    Percentile is computed here across all pairs (cDAV + uDAV) per gene.
    """
    dca = pd.read_csv(DCA)
    # Compute within-gene percentile over all pairs (not just cDAV)
    dca["mi_apc_pct"] = dca.groupby("dar_gene")["mi_apc"].rank(pct=True) * 100
    dca_sub = (
        dca.groupby(JOIN, as_index=False)
        .agg(mi_apc_pct=("mi_apc_pct", "first"), mi_apc=("mi_apc", "first"),
             dca_di_pct=("dca_di_percentile", "first"), dca_meff=("dca_meff", "first"))
    )
    return universe.merge(dca_sub, on=JOIN, how="left")


def sig_flags(df: pd.DataFrame) -> pd.DataFrame:
    """Add boolean significance columns for each test."""
    df = df.copy()
    df["sig_pagel"]   = (df["pagel_fdr"] <= PAGEL_ALPHA) | (df["branch_fdr"] <= BRANCH_ALPHA)
    df["sig_pyvolve"] = df["pyvolve_p"] < PYVOLVE_ALPHA
    df["sig_apc_mi"]  = df["mi_apc_pct"] >= APC_PERCENTILE_THRESHOLD
    df["sig_any"]     = df["sig_pagel"] | df["sig_pyvolve"] | df["sig_apc_mi"]
    df["n_sig_tests"] = df["sig_pagel"].astype(int) + df["sig_pyvolve"].astype(int) + df["sig_apc_mi"].astype(int)
    return df


def mw_test(a: pd.Series, b: pd.Series) -> tuple[float, float]:
    """Mann-Whitney U, one-sided (a > b). Returns (U, p). NaN if insufficient data."""
    a = a.dropna()
    b = b.dropna()
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan
    u, p = mannwhitneyu(a, b, alternative="greater")
    return float(u), float(p)


def compartment_summary(df: pd.DataFrame, label: str) -> list[dict]:
    """Compute summary statistics for a single compartment slice."""
    rows = []
    cdav = df[df["is_cdav"]]
    udav = df[~df["is_cdav"]]

    base = {
        "compartment": label,
        "n_total_pairs": len(df),
        "n_cdav_pairs": len(cdav),
        "n_udav_pairs": len(udav),
        "cdav_fraction": len(cdav) / len(df) if len(df) > 0 else np.nan,
    }

    # Compensation rates (cDAV only)
    for test_col, test_label in [
        ("sig_pagel",   "pagel_or_branch"),
        ("sig_pyvolve", "pyvolve"),
        ("sig_apc_mi",  "apc_mi_top_quartile"),
        ("sig_any",     "any_test"),
    ]:
        n_sig = cdav[test_col].sum() if len(cdav) > 0 else 0
        base[f"cdav_sig_{test_label}_n"] = int(n_sig)
        base[f"cdav_sig_{test_label}_rate"] = n_sig / len(cdav) if len(cdav) > 0 else np.nan

    # APC-MI: cDAV vs uDAV (within-gene percentile)
    u_apc, p_apc = mw_test(cdav["mi_apc_pct"], udav["mi_apc_pct"])
    base["apc_mi_pct_median_cdav"] = cdav["mi_apc_pct"].median()
    base["apc_mi_pct_median_udav"] = udav["mi_apc_pct"].median()
    base["apc_mi_mw_p_cdav_gt_udav"] = p_apc

    # DCA: cDAV vs uDAV
    u_dca, p_dca = mw_test(cdav["dca_di_pct"], udav["dca_di_pct"])
    base["dca_di_pct_median_cdav"] = cdav["dca_di_pct"].median()
    base["dca_di_pct_median_udav"] = udav["dca_di_pct"].median()
    base["dca_di_mw_p_cdav_gt_udav"] = p_dca

    # Pyvolve OR distribution (cDAV only)
    base["pyvolve_or_median_cdav"] = cdav["pyvolve_or"].median()
    base["pyvolve_or_p95_cdav"] = cdav["pyvolve_or"].quantile(0.95) if len(cdav) > 0 else np.nan
    base["pyvolve_tested_cdav_n"] = int(cdav["pyvolve_p"].notna().sum())

    rows.append(base)
    return rows


def mito_nuclear_detail(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return one row per mito-nuclear pair with all test results and direction annotation.
    Direction: 'mt_dar_nuc_contact' = mtDNA DAR contacts a nuclear subunit
               'nuc_dar_mt_contact' = nuclear DAR contacts an mtDNA subunit
    """
    mt_nuc = df[df["contact_type"] == "mt-nuc"].copy()
    mt_nuc["direction"] = np.where(
        mt_nuc["dar_genome"] == "mtDNA",
        "mt_dar_nuc_contact",
        "nuc_dar_mt_contact",
    )
    return mt_nuc


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading and joining data...")
    universe = load_universe()
    print(f"  Universe: {len(universe)} unique pairs")

    universe = add_pagel(universe)
    universe = add_pyvolve(universe)
    universe = add_apc_mi(universe)
    universe = sig_flags(universe)

    print(f"  Pairs with pyvolve: {universe['pyvolve_p'].notna().sum()}")
    print(f"  Pairs with APC-MI:  {universe['mi_apc_pct'].notna().sum()}")
    print(f"  Pairs with Pagel:   {universe['pagel_fdr'].notna().sum()}")

    # ── Compartment summaries ──────────────────────────────────────────────────
    summary_rows = []
    compartments = {
        "all":     universe,
        "mt-mt":   universe[universe["contact_type"] == "mt-mt"],
        "nuc-nuc": universe[universe["contact_type"] == "nuc-nuc"],
        "mt-nuc":  universe[universe["contact_type"] == "mt-nuc"],
    }
    for label, subset in compartments.items():
        if len(subset) == 0:
            continue
        summary_rows.extend(compartment_summary(subset, label))

    # Mito-nuclear direction sub-strata
    mt_nuc = universe[universe["contact_type"] == "mt-nuc"].copy()
    mt_nuc["direction"] = np.where(
        mt_nuc["dar_genome"] == "mtDNA", "mt_dar_nuc_contact", "nuc_dar_mt_contact"
    )
    for direction, subset in mt_nuc.groupby("direction"):
        summary_rows.extend(compartment_summary(subset, f"mt-nuc:{direction}"))

    summary = pd.DataFrame(summary_rows)
    summary_path = OUT_DIR / "compartment_analysis.csv"
    summary.to_csv(summary_path, index=False)
    print(f"\nWrote {len(summary)} rows → {summary_path.relative_to(ROOT)}")

    # ── Mito-nuclear detail ────────────────────────────────────────────────────
    detail = mito_nuclear_detail(universe)
    detail_path = OUT_DIR / "mito_nuclear_detail.csv"
    detail.to_csv(detail_path, index=False)
    print(f"Wrote {len(detail)} rows → {detail_path.relative_to(ROOT)}")

    # ── Print headline results ─────────────────────────────────────────────────
    print("\n═══ Compartment Analysis ═══")
    print(f"{'Compartment':<25} {'n_cdav':>8} {'n_udav':>8} "
          f"{'cdav_frac':>10} {'sig_any%':>9} "
          f"{'APC-MI cDAV':>12} {'APC-MI uDAV':>12} {'APC-MI p':>10}")
    print("-" * 100)
    for _, r in summary.iterrows():
        print(
            f"{r['compartment']:<25} {int(r['n_cdav_pairs']):>8} {int(r['n_udav_pairs']):>8} "
            f"{r['cdav_fraction']:>10.3f} "
            f"{100*r['cdav_sig_any_test_rate']:>8.1f}% "
            f"{r['apc_mi_pct_median_cdav']:>12.1f} "
            f"{r['apc_mi_pct_median_udav']:>12.1f} "
            f"{r['apc_mi_mw_p_cdav_gt_udav']:>10.2e}"
        )

    # ── Mito-nuclear top hits ──────────────────────────────────────────────────
    mt_sig = detail[detail["sig_any"]].sort_values("n_sig_tests", ascending=False)
    print(f"\nMito-nuclear pairs significant in ≥1 test: {len(mt_sig)}")
    if len(mt_sig) > 0:
        cols = JOIN + ["contact_type", "direction", "is_cdav",
                       "pagel_fdr", "branch_fdr", "pyvolve_p", "pyvolve_or",
                       "mi_apc_pct", "n_sig_tests"]
        print(mt_sig[cols].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
