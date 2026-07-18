"""
src/phylo/11_tier1_extraction.py

Step 3.5 — Tier-1 compensatory pair extraction for mutagenesis prioritization.

Extracts pairs satisfying: Pagel-sig (FDR ≤ 0.10) + favorable temporal ordering
(dominant_refined_timing ∈ {contact_first, permissive_background}).

Pyvolve (perm_p < 0.05) is added as a bonus column where available, but NOT used
as a filter because compensatory_partners.csv (Pagel pairs) and
conditional_permissiveness.csv (Pyvolve pairs) are largely disjoint pair sets.
Both datasets are included as supplemental tiers.

Tiers:
  Tier-1A: Pagel-sig + permissive timing (canonical)
  Tier-1B: Pagel-sig only (timing pending / insufficient)
  Tier-1C: Mito-nuclear Pyvolve-sig + permissive timing (cross-interface)

Output: results/phylo/tier1_pairs.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

COMP_PART  = ROOT / "results" / "structural" / "compensatory_partners.csv"
PERM_FILE  = ROOT / "results" / "phylo" / "conditional_permissiveness.csv"
TIMING_V2  = ROOT / "results" / "phylo" / "timing_annotations_v2.csv"
TIMING_MN  = ROOT / "results" / "phylo" / "timing_mitonuclear_extended.csv"
OUT_FILE   = ROOT / "results" / "phylo" / "tier1_pairs.csv"

PAGEL_ALPHA = 0.10
PERM_ALPHA  = 0.05
PERMISSIVE_CLASSES = {"contact_first", "permissive_background"}


def main() -> None:
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    comp   = pd.read_csv(COMP_PART)
    perm   = pd.read_csv(PERM_FILE)
    timing = pd.read_csv(TIMING_V2)
    timing_mn = pd.read_csv(TIMING_MN)

    join4 = ["dar_gene","dar_aa_coord","contact_gene","contact_refseq_pos"]

    # ── Tier-1A/1B: Pagel-sig pairs with timing ────────────────────────────────
    pagel_sig = comp[comp["pagel_fdr"] <= PAGEL_ALPHA].copy()
    print(f"Pagel-significant pairs: {len(pagel_sig)}")

    # Merge timing into Pagel-sig pairs
    timing_sub = timing[join4 + ["dominant_refined_timing","timing_confidence",
                                  "n_dar_gain_branches","n_contact_first",
                                  "n_permissive_background","n_contact_after",
                                  "phys_blosum_median","phys_miyata_median"]].drop_duplicates()
    pagel_timed = pagel_sig.merge(timing_sub, on=join4, how="left")
    print(f"  With timing data: {pagel_timed['dominant_refined_timing'].notna().sum()}")

    # Optionally add perm_p
    perm_sub = perm[join4 + ["perm_p","observed_or"]].drop_duplicates(subset=join4)
    pagel_timed = pagel_timed.merge(perm_sub, on=join4, how="left")
    print(f"  With Pyvolve data: {pagel_timed['perm_p'].notna().sum()}")

    # Tier-1A: Pagel + permissive timing
    tier1a = pagel_timed[
        pagel_timed["dominant_refined_timing"].isin(PERMISSIVE_CLASSES)
    ].copy()
    tier1a["tier"] = "1A"

    # Tier-1B: Pagel-sig without permissive timing (missing timing or non-permissive)
    tier1b = pagel_timed[
        ~pagel_timed["dominant_refined_timing"].isin(PERMISSIVE_CLASSES)
    ].copy()
    tier1b["tier"] = "1B"

    print(f"\nTier-1A (Pagel + permissive timing): {len(tier1a)}")
    print(f"Tier-1B (Pagel, timing absent/non-permissive): {len(tier1b)}")

    # ── Tier-1C: Mito-nuclear Pyvolve-sig + permissive timing ──────────────────
    # timing_mn comes from 06b_temporal_ordering_mitonuc.py
    mn_sig = timing_mn[
        (timing_mn["perm_p"] < PERM_ALPHA) &
        timing_mn["dominant_refined_timing"].isin(PERMISSIVE_CLASSES)
    ].copy()
    mn_sig["tier"] = "1C"
    mn_sig["pagel_fdr"] = np.nan  # not Pagel-tested
    print(f"Tier-1C (mito-nuclear Pyvolve + permissive timing): {len(mn_sig)}")

    # ── Combine and format output ───────────────────────────────────────────────
    output_cols = ["tier","dar_gene","dar_aa_coord","dar_alt_aa",
                   "contact_gene","contact_refseq_pos","contact_alt_aa",
                   "contact_class","contact_type",
                   "pagel_fdr","fisher_fdr","perm_p","observed_or",
                   "dominant_refined_timing","timing_confidence",
                   "n_dar_gain_branches","n_contact_first",
                   "n_permissive_background","n_contact_after",
                   "phys_blosum_median","phys_miyata_median"]

    def safe_select(df: pd.DataFrame, cols: list) -> pd.DataFrame:
        return df[[c for c in cols if c in df.columns]]

    tier1a_out = safe_select(tier1a, output_cols)
    tier1b_out = safe_select(tier1b, output_cols)
    tier1c_out = safe_select(mn_sig, output_cols)

    combined = pd.concat([tier1a_out, tier1c_out, tier1b_out], ignore_index=True)
    combined.to_csv(OUT_FILE, index=False)
    print(f"\nSaved {len(combined)} total tier-1 pairs → {OUT_FILE.relative_to(ROOT)}")

    # ── Print summary ───────────────────────────────────────────────────────────
    print("\n=== Tier-1A (mutagenesis priority) ===")
    if len(tier1a) > 0:
        sdh_mask = tier1a["dar_gene"].isin({"SDHA","SDHB","SDHC","SDHD","SDHAF2"})
        print(f"  Total: {len(tier1a)} ({(~sdh_mask).sum()} non-SDH)")
        disp = tier1a.sort_values("pagel_fdr")[
            ["dar_gene","dar_aa_coord","dar_alt_aa","contact_gene","contact_refseq_pos",
             "pagel_fdr","dominant_refined_timing","timing_confidence",
             "n_contact_first","n_permissive_background","phys_blosum_median"]]
        print(disp.to_string(index=False))

    print("\n=== Tier-1C (mito-nuclear cross-interface) ===")
    if len(mn_sig) > 0:
        disp = mn_sig.sort_values("perm_p")[
            ["dar_gene","dar_aa_coord","dar_alt_aa","contact_gene","contact_refseq_pos",
             "perm_p","observed_or","dominant_refined_timing","timing_confidence"]]
        print(disp.to_string(index=False))

    print("\n=== Tier-1B (Pagel-sig, timing absent or non-permissive) ===")
    if len(tier1b) > 0:
        disp = tier1b.sort_values("pagel_fdr")[
            ["dar_gene","dar_aa_coord","dar_alt_aa","contact_gene","contact_refseq_pos",
             "pagel_fdr","dominant_refined_timing","timing_confidence"]]
        print(disp.to_string(index=False))


if __name__ == "__main__":
    main()
