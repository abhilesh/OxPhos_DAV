"""
src/phylo/06b_temporal_ordering_mitonuc.py

Step 3.1 — Extended temporal ordering for all Pyvolve-significant mito-nuclear pairs.

Runs the same per-branch timing analysis as 06_temporal_ordering.py but on
the 45 mito-nuclear DAR-contact pairs where perm_p < 0.05 (from
conditional_permissiveness.csv), stratified by direction:
  mt-dar:  mtDNA DAR gene → nuclear contact gene
  nuc-dar: nuclear DAR gene → mtDNA contact gene

Primary question: do mito-nuclear pairs show structural pre-adaptation
(contact_first + permissive_background) vs rescue (contact_after), and
does this differ by direction?

Output: results/phylo/timing_mitonuclear_extended.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from phylo._timing_core import analyze_pair, TIMING_CATS, REFINED_CATS

ASR_MAPS   = ROOT / "data" / "phylo" / "ancestral_state_maps.json"
PERM_FILE  = ROOT / "results" / "phylo" / "conditional_permissiveness.csv"
CONTACTS   = ROOT / "results" / "structural" / "dar_contacts_cbcb8A.csv"
OUT_FILE   = ROOT / "results" / "phylo" / "timing_mitonuclear_extended.csv"

MTDNA = {"MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
         "MT-CYB","MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8"}

WITHIN_GENOME_BASELINE = 0.632  # permissive rate from 44 Pagel pairs (63.2%)


def summarise_direction(df: pd.DataFrame, direction: str) -> None:
    sub = df[df["direction"] == direction].copy()
    testable = sub[~sub["is_ancestral_cdav"] & (sub["n_dar_gain_branches"] > 0)]
    print(f"\n── Direction: {direction} (pairs={len(sub)}, testable={len(testable)}) ──")
    if len(testable) == 0:
        print("  No testable pairs.")
        return

    n_cf  = int(testable["n_contact_first"].sum())
    n_co  = int(testable["n_co_occurring"].sum())
    n_ca  = int(testable["n_contact_after"].sum())
    n_ncc = int(testable["n_no_contact_change"].sum())
    n_pb  = int(testable.get("n_permissive_background", pd.Series([0])).sum()
                if "n_permissive_background" in testable.columns else 0)
    n_coa = int(testable.get("n_co_adaptation", pd.Series([0])).sum()
                if "n_co_adaptation" in testable.columns else 0)
    n_cp  = int(testable.get("n_constitutively_permissive", pd.Series([0])).sum()
                if "n_constitutively_permissive" in testable.columns else 0)

    total_events = n_cf + n_co + n_ca + n_ncc
    permissive = n_cf + n_pb + n_cp
    print(f"  Total branch events: {total_events}")

    print(f"  Refined counts:")
    remaining_ncc = n_ncc - n_cp
    for cat, n in [("contact_first", n_cf), ("permissive_background", n_pb),
                   ("co_adaptation", n_coa), ("constitutively_permissive", n_cp),
                   ("contact_after", n_ca), ("no_contact_change", remaining_ncc)]:
        pct = 100 * n / total_events if total_events > 0 else 0
        print(f"    {cat:<30} {n:>5}  ({pct:.1f}%)")

    print(f"\n  Permissive signal (cf+pb+cp): {permissive}/{total_events} "
          f"({100*permissive/total_events:.1f}%)")
    print(f"  Rescue (contact_after):       {n_ca}/{total_events} "
          f"({100*n_ca/total_events:.1f}%)")

    if n_cf + n_ca > 0:
        bt = binomtest(n_cf, n_cf + n_ca, 0.5, alternative="greater")
        print(f"  Directional test (contact_first > contact_after): "
              f"n_first={n_cf} n_after={n_ca}  p={bt.pvalue:.3e}")

    # Rate control vs within-genome baseline
    bt2 = binomtest(permissive, total_events, WITHIN_GENOME_BASELINE, alternative="two-sided")
    print(f"  vs within-genome baseline ({100*WITHIN_GENOME_BASELINE:.1f}%): "
          f"p={bt2.pvalue:.3e}")

    # Dominant refined timing per pair
    dom = testable["dominant_refined_timing"].value_counts()
    print(f"  Dominant refined timing per pair:")
    for k, v in dom.items():
        print(f"    {k:<30} {v:>4}  ({100*v/len(testable):.1f}%)")

    # Top pairs
    top = testable.sort_values("n_dar_gain_branches", ascending=False).head(5)
    print(f"\n  Top pairs by n_gain_branches:")
    for _, r in top.iterrows():
        print(f"    {r['dar_gene']}:{r['dar_aa_coord']}{r['dar_alt_aa']}"
              f"→{r['contact_gene']}:{r['contact_refseq_pos']}"
              f"  branches={r['n_dar_gain_branches']}"
              f"  dom={r['dominant_refined_timing']}")


def main() -> None:
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    perm = pd.read_csv(PERM_FILE)
    contacts = pd.read_csv(CONTACTS, low_memory=False)

    sig_mn = perm[(perm["perm_p"] < 0.05) & (perm["contact_type"] == "mt-nuc")].copy()
    print(f"  Pyvolve-significant mt-nuc pairs: {len(sig_mn)}")

    # Get reference AAs from structural contacts
    contacts2 = contacts.rename(columns={
        "dar_structure_gene": "dar_gene",
        "contact_aa": "contact_human_aa"
    })
    merge_keys = ["dar_gene","dar_aa_coord","dar_alt_aa","contact_gene","contact_refseq_pos"]
    ref_info = (contacts2[merge_keys + ["dar_ref_aa","contact_human_aa"]]
                .drop_duplicates(subset=merge_keys))

    sig_mn = sig_mn.merge(ref_info, on=merge_keys, how="left")
    n_matched = sig_mn["dar_ref_aa"].notna().sum()
    print(f"  Reference AA lookup: {n_matched}/{len(sig_mn)} matched")

    # Assign direction
    sig_mn["direction"] = sig_mn["dar_gene"].apply(
        lambda g: "mt-dar" if g in MTDNA else "nuc-dar")
    print(f"  mt-dar:  {(sig_mn['direction']=='mt-dar').sum()}")
    print(f"  nuc-dar: {(sig_mn['direction']=='nuc-dar').sum()}")

    print("\nLoading ASR maps...")
    with open(ASR_MAPS) as f:
        asr_data = json.load(f)
    print(f"  Genes in ASR maps: {len(asr_data)}")

    print("\nRunning temporal ordering on mito-nuclear pairs...")
    node_states_cache: dict = {}
    results = []
    for i, (_, row) in enumerate(sig_mn.iterrows()):
        if i % 10 == 0:
            print(f"  {i+1}/{len(sig_mn)}  {row['dar_gene']} "
                  f"{row['dar_aa_coord']}→{row['contact_gene']} "
                  f"{row['contact_refseq_pos']}")
        res = analyze_pair(row, asr_data, node_states_cache)
        res["direction"] = row["direction"]
        res["perm_p"]    = row["perm_p"]
        res["observed_or"] = row.get("observed_or", np.nan)
        results.append(res)

    df = pd.DataFrame(results)
    df.to_csv(OUT_FILE, index=False)
    print(f"\nWrote {len(df)} rows → {OUT_FILE.relative_to(ROOT)}")

    # ── Summary ────────────────────────────────────────────────────────────────
    print("\n═══ Mito-Nuclear Temporal Ordering Summary ═══")
    print(f"Total pairs: {len(df)}")
    print(f"Ancestral cDAV: {df['is_ancestral_cdav'].sum()}")
    print(f"No gain branches: {(df['timing_confidence']=='no_gains_found').sum()}")
    print(f"Missing ASR gene: {(df['timing_confidence']=='missing_gene').sum()}")
    testable = df[~df["is_ancestral_cdav"] & (df["n_dar_gain_branches"] > 0)]
    print(f"Testable (≥1 gain branch): {len(testable)}")

    for direction in ["mt-dar", "nuc-dar"]:
        summarise_direction(df, direction)

    # Combined
    print("\n═══ Combined mito-nuclear ═══")
    n_cf  = int(testable["n_contact_first"].sum())
    n_pb  = int(testable.get("n_permissive_background", pd.Series([0])).sum()
                if "n_permissive_background" in testable.columns else 0)
    n_cp  = int(testable.get("n_constitutively_permissive", pd.Series([0])).sum()
                if "n_constitutively_permissive" in testable.columns else 0)
    n_ca  = int(testable["n_contact_after"].sum())
    n_co  = int(testable["n_co_occurring"].sum())
    n_ncc = int(testable["n_no_contact_change"].sum())
    total = n_cf + n_co + n_ca + n_ncc
    permissive = n_cf + n_pb + n_cp
    if total > 0:
        print(f"  Permissive signal: {permissive}/{total} ({100*permissive/total:.1f}%)")
        print(f"  Rescue: {n_ca}/{total} ({100*n_ca/total:.1f}%)")
        bt = binomtest(permissive, total, WITHIN_GENOME_BASELINE, alternative="two-sided")
        print(f"  vs within-genome baseline ({100*WITHIN_GENOME_BASELINE:.1f}%): p={bt.pvalue:.3e}")


if __name__ == "__main__":
    main()
