"""
src/phylo/08_pagel_complex_stratification.py

Step 3.4 — Pagel complex stratification and contact-type enrichment.

Assigns each of the 44 Pagel-significant pairs to an OXPHOS complex, tests
whether specific structural contact types (H-bond, electrostatic, hydrophobic,
vdW) are enriched among significant vs non-significant pairs, and checks whether
SDH-excluded nuclear signal spans multiple complexes.

Output: results/phylo/pagel_by_complex.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

COMP_PART = ROOT / "results" / "structural" / "compensatory_partners.csv"
OUT_FILE  = ROOT / "results" / "phylo" / "pagel_by_complex.csv"

MTDNA = {"MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
         "MT-CYB","MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8"}
SDH   = {"SDHA","SDHB","SDHC","SDHD","SDHAF2"}

GENE_TO_COMPLEX = {
    # Complex I (NADH:ubiquinone oxidoreductase) — nucDNA
    "NDUFV1":"CI","NDUFV2":"CI","NDUFV3":"CI",
    "NDUFS1":"CI","NDUFS2":"CI","NDUFS3":"CI","NDUFS4":"CI","NDUFS5":"CI",
    "NDUFS6":"CI","NDUFS7":"CI","NDUFS8":"CI",
    "NDUFA1":"CI","NDUFA2":"CI","NDUFA3":"CI","NDUFA4":"CI","NDUFA5":"CI",
    "NDUFA6":"CI","NDUFA7":"CI","NDUFA8":"CI","NDUFA9":"CI","NDUFA10":"CI",
    "NDUFA11":"CI","NDUFA12":"CI","NDUFA13":"CI",
    "NDUFB1":"CI","NDUFB2":"CI","NDUFB3":"CI","NDUFB4":"CI","NDUFB5":"CI",
    "NDUFB6":"CI","NDUFB7":"CI","NDUFB8":"CI","NDUFB9":"CI","NDUFB10":"CI",
    "NDUFB11":"CI",
    "NDUFAB1":"CI",
    # Complex I — mtDNA
    "MT-ND1":"CI","MT-ND2":"CI","MT-ND3":"CI","MT-ND4":"CI","MT-ND4L":"CI",
    "MT-ND5":"CI","MT-ND6":"CI",
    # Complex II (succinate dehydrogenase) — all nucDNA
    "SDHA":"CII","SDHB":"CII","SDHC":"CII","SDHD":"CII","SDHAF2":"CII",
    # Complex III (cytochrome bc1) — mixed
    "MT-CYB":"CIII",
    "CYC1":"CIII","CYCS":"CIII",
    "UQCRC1":"CIII","UQCRC2":"CIII","UQCRB":"CIII","UQCRH":"CIII",
    "UQCR10":"CIII","UQCR11":"CIII","UQCRFS1":"CIII","UQCRQ":"CIII",
    "UQCRSF1":"CIII",
    # Complex IV (cytochrome c oxidase) — mixed
    "MT-CO1":"CIV","MT-CO2":"CIV","MT-CO3":"CIV",
    "COX4I1":"CIV","COX4I2":"CIV","COX5A":"CIV","COX5B":"CIV",
    "COX6A1":"CIV","COX6A2":"CIV","COX6B1":"CIV","COX6B2":"CIV",
    "COX6C":"CIV","COX7A1":"CIV","COX7A2":"CIV","COX7B":"CIV",
    "COX8A":"CIV","COX8C":"CIV",
    "COX4I1":"CIV","COXFA4":"CIV","COXFA4L2":"CIV","COXFA4L3":"CIV",
    # Complex V (ATP synthase) — mixed
    "MT-ATP6":"CV","MT-ATP8":"CV",
    "ATP5F1A":"CV","ATP5F1B":"CV","ATP5F1C":"CV","ATP5F1D":"CV","ATP5F1E":"CV",
    "ATP5IF1":"CV","ATP5MC1":"CV","ATP5MC2":"CV","ATP5MC3":"CV",
    "ATP5MF":"CV","ATP5MG":"CV","ATP5MJ":"CV","ATP5MK":"CV",
    "ATP5PB":"CV","ATP5PD":"CV","ATP5PF":"CV","ATP5PO":"CV",
}


def assign_complex(gene: str) -> str:
    return GENE_TO_COMPLEX.get(gene, "unknown")


def main() -> None:
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    cp = pd.read_csv(COMP_PART)
    cp["dar_compartment"] = cp["dar_gene"].apply(
        lambda g: "mtDNA" if g in MTDNA else "nucDNA")
    cp["is_sdh"] = cp["dar_gene"].isin(SDH)
    cp["complex"] = cp["dar_gene"].apply(assign_complex)
    cp["pagel_sig"] = cp["pagel_fdr"] <= 0.10

    sig  = cp[cp["pagel_sig"]]
    nsig = cp[~cp["pagel_sig"]]

    # ── 1. Significant pairs by complex ────────────────────────────────────────
    print("=== Pagel significant pairs by complex ===")
    by_complex = (
        sig.groupby(["complex","dar_compartment","is_sdh"])
        .agg(n_sig=("dar_gene","count"),
             genes=("dar_gene", lambda x: ", ".join(sorted(set(x)))))
        .reset_index()
    )
    print(by_complex.to_string(index=False))

    print("\n=== non-SDH nuclear pairs by complex ===")
    non_sdh_nuc = sig[(sig.dar_compartment=="nucDNA") & ~sig.is_sdh]
    print(non_sdh_nuc.groupby("complex").agg(
        n=("dar_gene","count"),
        genes=("dar_gene", lambda x: ", ".join(sorted(set(x)))),
        pagel_p_range=("pagel_p", lambda x: f"{x.min():.2e}–{x.max():.2e}")
    ).to_string())

    # ── 2. Contact-class enrichment (H-bond, electrostatic, hydrophobic, vdW) ──
    print("\n=== Contact-class distribution: sig vs non-sig Pagel pairs ===")
    ct_sig  = sig["contact_class"].value_counts()
    ct_nsig = nsig["contact_class"].value_counts()
    all_types = sorted(set(list(ct_sig.index) + list(ct_nsig.index)))
    ct_df = pd.DataFrame({
        "contact_class": all_types,
        "n_sig":   [ct_sig.get(t, 0) for t in all_types],
        "n_nonsig":[ct_nsig.get(t, 0) for t in all_types],
    })
    ct_df["pct_sig"]   = ct_df["n_sig"]   / ct_df["n_sig"].sum() * 100
    ct_df["pct_nonsig"]= ct_df["n_nonsig"]/ ct_df["n_nonsig"].sum() * 100
    print(ct_df.to_string(index=False))

    # Chi-squared over the contingency table
    table = ct_df[["n_sig","n_nonsig"]].values
    if table.min() > 0:
        chi2, p, dof, _ = chi2_contingency(table)
        print(f"\nChi-squared ({len(all_types)} contact classes, sig vs non-sig): χ²={chi2:.2f}, df={dof}, p={p:.3e}")
    else:
        chi2, p, dof, _ = chi2_contingency(table)
        print(f"\nChi-squared ({len(all_types)} contact classes, sig vs non-sig): χ²={chi2:.2f}, df={dof}, p={p:.3e}")
        print("  (some cells are zero; Fisher test per class is more reliable)")

    # Per-class Fisher
    total_sig   = ct_df["n_sig"].sum()
    total_nonsig= ct_df["n_nonsig"].sum()
    print("\nPer-class Fisher exact (enriched in sig?):")
    for _, row in ct_df.iterrows():
        a, b = int(row["n_sig"]), int(row["n_nonsig"])
        c, d = total_sig - a, total_nonsig - b
        _, p_fe = fisher_exact([[a, b],[c, d]], alternative="greater")
        print(f"  {row['contact_class']:<15} sig={a}/{total_sig} ({100*a/total_sig:.1f}%)  "
              f"nonsig={b}/{total_nonsig} ({100*b/total_nonsig:.1f}%)  fisher_p={p_fe:.3e}")

    # ── 3. Save output ─────────────────────────────────────────────────────────
    out = sig[["dar_gene","dar_aa_coord","dar_ref_aa","dar_alt_aa",
               "contact_gene","contact_refseq_pos","contact_human_aa","contact_alt_aa",
               "contact_class","contact_type","dar_compartment","is_sdh","complex",
               "pagel_p","pagel_fdr","fisher_p","fisher_fdr"]].copy()
    out.to_csv(OUT_FILE, index=False)
    print(f"\nSaved {len(out)} rows → {OUT_FILE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
