"""
src/phylo/10_nuclear_dca_stratified.py

Step 3.3 — High-Meff nuclear DCA stratified analysis.

Asks: at what Meff threshold does nuclear DCA-DI signal emerge?
Stratifies nuclear DAR-contact pairs by per-gene median Meff into bins
(≤2, 2-5, 5-10, 10-20, >20) and runs Mann-Whitney U for cDAV vs uDAV
DCA-DI percentile within each bin.

Also adds branch-change density (n_branches_with_changes from ASR audit)
as an independent evolutionary rate proxy and tests if it mediates the
Meff-OR relationship.

Output: results/mutagenesis/dca_nuclear_meff_stratified.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

DCA_FILE   = ROOT / "results" / "mutagenesis" / "dca_all_davs.csv"
AUDIT_FILE = ROOT / "data" / "phylo" / "asr_coordinate_harmonization_audit.tsv"
OUT_FILE   = ROOT / "results" / "mutagenesis" / "dca_nuclear_meff_stratified.csv"

BINS = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 10000)]
BIN_LABELS = ["≤2", "2–5", "5–10", "10–20", ">20"]


def main() -> None:
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    dca = pd.read_csv(DCA_FILE)
    audit = pd.read_csv(AUDIT_FILE, sep="\t")

    # Nuclear pairs only
    nuc = dca[dca["dar_genome"] == "nucDNA"].copy()
    print(f"Nuclear DAR-contact pairs: {len(nuc)}")

    # Per-gene median Meff and branch-change count
    gene_meff = (nuc.groupby("dar_gene")["dca_meff"]
                 .median()
                 .reset_index()
                 .rename(columns={"dca_meff": "gene_median_meff"}))
    audit_sub = audit[["gene", "n_branches_with_changes"]].rename(
        columns={"gene": "dar_gene"})
    gene_meff = gene_meff.merge(audit_sub, on="dar_gene", how="left")

    nuc = nuc.merge(gene_meff, on="dar_gene", how="left")

    print(f"\nGene median Meff distribution:")
    print(gene_meff["gene_median_meff"].describe().round(2))
    print(f"\nGenes per Meff bin:")
    for (lo, hi), label in zip(BINS, BIN_LABELS):
        ng = gene_meff[(gene_meff["gene_median_meff"] > lo) &
                       (gene_meff["gene_median_meff"] <= hi)]
        print(f"  {label}: {len(ng)} genes")

    print("\n=== Mann-Whitney (cDAV vs uDAV DCA-DI percentile) by Meff bin ===")
    rows = []
    for (lo, hi), label in zip(BINS, BIN_LABELS):
        grp = nuc[(nuc["gene_median_meff"] > lo) & (nuc["gene_median_meff"] <= hi)]
        cdav = grp[grp["is_cdav_amino_acid"]]["dca_di_percentile"].dropna()
        udav = grp[~grp["is_cdav_amino_acid"]]["dca_di_percentile"].dropna()
        if len(cdav) < 5 or len(udav) < 5:
            p_mw, u_stat = np.nan, np.nan
            direction = "—"
        else:
            u_stat, p_mw = mannwhitneyu(cdav, udav, alternative="greater")
            direction = "cDAV>uDAV" if p_mw < 0.05 else "ns"

        n_genes = int(nuc[(nuc["gene_median_meff"] > lo) &
                          (nuc["gene_median_meff"] <= hi)]["dar_gene"].nunique())
        med_meff = grp["gene_median_meff"].median() if len(grp) > 0 else np.nan
        rows.append({
            "meff_bin": label, "n_genes": n_genes,
            "n_cdav": len(cdav), "n_udav": len(udav),
            "med_meff": med_meff,
            "median_di_pct_cdav": cdav.median() if len(cdav) > 0 else np.nan,
            "median_di_pct_udav": udav.median() if len(udav) > 0 else np.nan,
            "mw_u": u_stat, "mw_p": p_mw, "direction": direction,
        })
        print(f"  Meff {label:<8}  genes={n_genes:>3}  "
              f"n_cDAV={len(cdav):>5}  n_uDAV={len(udav):>5}  "
              f"med_DI_cDAV={cdav.median():.1f}  med_DI_uDAV={udav.median():.1f}  "
              f"MW_p={p_mw:.3e}  [{direction}]")

    df = pd.DataFrame(rows)

    # Spearman between gene median Meff and per-gene MW p-value
    print("\n=== Spearman r (gene median Meff vs per-gene MW p for cDAV>uDAV) ===")
    per_gene_rows = []
    for gene, grp in nuc.groupby("dar_gene"):
        cdav = grp[grp["is_cdav_amino_acid"]]["dca_di_percentile"].dropna()
        udav = grp[~grp["is_cdav_amino_acid"]]["dca_di_percentile"].dropna()
        if len(cdav) < 3 or len(udav) < 3:
            continue
        _, p = mannwhitneyu(cdav, udav, alternative="greater")
        med_meff = grp["gene_median_meff"].iloc[0]
        nbranch  = grp["n_branches_with_changes"].iloc[0]
        per_gene_rows.append({"dar_gene": gene, "gene_median_meff": med_meff,
                               "n_branches": nbranch, "mw_p": p,
                               "n_cdav": len(cdav), "n_udav": len(udav)})
    pg = pd.DataFrame(per_gene_rows)
    if len(pg) > 5:
        log_p = -np.log10(pg["mw_p"].clip(lower=1e-10))
        r_meff, p_meff = spearmanr(pg["gene_median_meff"].fillna(0), log_p)
        print(f"  r(Meff, -log10(p)) = {r_meff:.3f}, p = {p_meff:.3e}")
        r_br, p_br = spearmanr(pg["n_branches"].fillna(0), log_p)
        print(f"  r(n_branches, -log10(p)) = {r_br:.3f}, p = {p_br:.3e}")
        print(f"\n  Top nuclear genes by Meff with cDAV>uDAV signal:")
        top = pg.sort_values("gene_median_meff", ascending=False).head(10)
        for _, row in top.iterrows():
            star = "*" if row["mw_p"] < 0.05 else ""
            print(f"    {row['dar_gene']:<12} Meff={row['gene_median_meff']:.1f}  "
                  f"branches={row['n_branches']:.0f}  "
                  f"cDAV_med={nuc[nuc['dar_gene']==row['dar_gene']][nuc['is_cdav_amino_acid']]['dca_di_percentile'].median():.1f}  "
                  f"p={row['mw_p']:.3e}{star}")

    df.to_csv(OUT_FILE, index=False)
    print(f"\nSaved {len(df)} rows → {OUT_FILE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
