"""
Temporal ordering analysis for cDAV compensation — Phase 2b.

For each Pagel/branch-significant DAR-contact pair, determines whether the
compensatory contact amino acid arose BEFORE the DAV (pre-adaptation), on the
SAME branch (co-occurring), or AFTER (rescue), using IQTree MAP ancestral states
pre-processed in ancestral_state_maps.json.

Coordinate contract: ASR state keys are human ungapped protein positions. IQTree
alignment-site IDs are harmonized in src/phylo/01_parse_ancestral_states.py.

Inputs:
  results/structural/compensatory_partners.csv   — Pagel/branch test results
  data/phylo/ancestral_state_maps.json           — pre-processed ASR per gene

Outputs:
  results/phylo/timing_annotations_v2.csv        — per-pair timing classification

Usage (inside Docker):
  python src/phylo/06_temporal_ordering.py
"""

import json
import numpy as np
import pandas as pd
from collections import deque
from pathlib import Path
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[2]
ASR_MAPS  = ROOT / "data" / "phylo" / "ancestral_state_maps.json"
COMP_PART = ROOT / "results" / "structural" / "compensatory_partners.csv"
OUT_DIR   = ROOT / "results" / "phylo"

PAGEL_ALPHA  = 0.10
BRANCH_ALPHA = 0.10
TIMING_CATS  = ["contact_first", "co_occurring", "contact_after", "no_contact_change"]


# ── Tree utilities ─────────────────────────────────────────────────────────────

def build_parent_map(node_to_children: dict, root: str) -> dict:
    """BFS from root, return {child: parent}."""
    parent_map: dict = {}
    q: deque = deque([root])
    while q:
        node = q.popleft()
        for child in node_to_children.get(node, []):
            parent_map[child] = node
            q.append(child)
    return parent_map


def compute_node_states(gene_asr: dict, site_str: str) -> dict:
    """
    Return {node: aa_state} for all nodes at a human protein position.
    Leaf nodes use observed states at the harmonized human alignment column;
    internal nodes use
    tree traversal from root_states + branch substitutions.
    None indicates a gap/ambiguous root state.
    """
    root = gene_asr["root_node"]
    leaf_nodes = set(gene_asr.get("leaf_nodes", []))
    leaf_states = gene_asr.get("leaf_states", {})

    states: dict = {}
    states[root] = gene_asr["root_states"].get(site_str)  # may be None if gapped

    q: deque = deque([root])
    while q:
        parent = q.popleft()
        parent_state = states[parent]
        for child in gene_asr["node_to_children"].get(parent, []):
            if child in leaf_nodes:
                # Use observed alignment state (most reliable for leaves)
                states[child] = leaf_states.get(child, {}).get(site_str, parent_state)
            else:
                bkey = f"{parent}|{child}"
                subs = gene_asr["branches"].get(bkey, {})
                if site_str in subs:
                    _, child_aa = subs[site_str]
                    states[child] = child_aa
                else:
                    states[child] = parent_state
            q.append(child)
    return states


def find_dav_gain_branches(gene_asr: dict, site_str: str, alt_aa: str,
                            parent_map: dict, dar_node_states: dict) -> list[str]:
    """
    Return list of branch keys "parent|child" where ANY amino acid → alt_aa transition
    occurs (i.e., any gain of the pathogenic amino acid regardless of starting state).

    This is intentionally broader than ref→alt: non-human cDAV species often acquired
    the pathogenic amino acid from the ancestral state (which is frequently different
    from the human reference amino acid), so requiring parent==ref_aa would miss the
    vast majority of gain events.

    Primary source: branches dict (internal transitions).
    Secondary: leaf nodes with alt_aa whose inferred parent has a different state.
    """
    gains: set = set()
    for bkey, subs in gene_asr["branches"].items():
        if site_str in subs:
            p_aa, c_aa = subs[site_str]
            if c_aa == alt_aa and p_aa != alt_aa:
                gains.add(bkey)

    for leaf in gene_asr.get("leaf_nodes", []):
        if gene_asr["leaf_states"].get(leaf, {}).get(site_str) == alt_aa:
            parent = parent_map.get(leaf)
            if parent and dar_node_states.get(parent) != alt_aa:
                bkey = f"{parent}|{leaf}"
                if bkey not in gains:
                    gains.add(bkey)

    return list(gains)


def get_clade_states(node_to_children: dict, node_states: dict,
                     start_node: str) -> set:
    """Return all AA states found in the clade rooted at start_node (DFS)."""
    vals: set = set()
    stack = [start_node]
    seen: set = set()
    while stack:
        n = stack.pop()
        if n in seen:
            continue
        seen.add(n)
        v = node_states.get(n)
        if v is not None:
            vals.add(v)
        stack.extend(node_to_children.get(n, []))
    return vals


def classify_branch_timing(contact_node_states: dict, node_to_children: dict,
                            parent_node: str, child_node: str,
                            contact_alt: str) -> str:
    """
    Classify timing of contact_alt appearance relative to the DAV gain on a single branch.
      contact_first  — parent already has contact_alt (pre-adaptation)
      co_occurring   — contact_alt appears on the same branch as the DAV
      contact_after  — contact_alt appears in a descendant branch (rescue)
      no_contact_change — contact_alt never appears in this clade
    """
    parent_contact = contact_node_states.get(parent_node)
    child_contact  = contact_node_states.get(child_node)

    if parent_contact == contact_alt:
        return "contact_first"
    if child_contact == contact_alt:
        return "co_occurring"

    clade_states = get_clade_states(node_to_children, contact_node_states, child_node)
    if contact_alt in clade_states:
        return "contact_after"
    return "no_contact_change"


# ── Per-pair analysis ──────────────────────────────────────────────────────────

def analyze_pair(row: pd.Series, asr_data: dict,
                 node_states_cache: dict) -> dict:
    dar_gene     = row["dar_gene"]
    dar_col      = str(int(row["dar_aa_coord"]))
    dar_ref      = row["dar_ref_aa"]
    dar_alt      = row["dar_alt_aa"]
    contact_gene = row["contact_gene"]
    contact_col  = str(int(row["contact_refseq_pos"]))
    contact_alt  = row["contact_alt_aa"]

    base = {
        "ann_id":            row.get("ann_id"),
        "dar_gene":          dar_gene,
        "dar_aa_coord":      row["dar_aa_coord"],
        "dar_ref_aa":        dar_ref,
        "dar_alt_aa":        dar_alt,
        "contact_gene":      contact_gene,
        "contact_refseq_pos": row["contact_refseq_pos"],
        "contact_human_aa":  row["contact_human_aa"],
        "contact_alt_aa":    contact_alt,
        "is_ancestral_cdav": False,
        "n_dar_gain_branches": 0,
        "n_contact_first":   0,
        "n_co_occurring":    0,
        "n_contact_after":   0,
        "n_no_contact_change": 0,
        "dominant_timing":   np.nan,
        "timing_confidence": "insufficient",
        "dar_origin_age_mya": np.nan,
    }

    if dar_gene not in asr_data or contact_gene not in asr_data:
        base["timing_confidence"] = "missing_gene"
        return base

    dar_asr     = asr_data[dar_gene]
    contact_asr = asr_data[contact_gene]

    # Cache: node states per gene per site
    g = node_states_cache.setdefault(dar_gene, {})
    if dar_col not in g:
        g[dar_col] = compute_node_states(dar_asr, dar_col)
    dar_states = g[dar_col]

    cg = node_states_cache.setdefault(contact_gene, {})
    if contact_col not in cg:
        cg[contact_col] = compute_node_states(contact_asr, contact_col)
    contact_states = cg[contact_col]

    # Build parent map (cached on the asr_data dict)
    if "_parent_map" not in dar_asr:
        dar_asr["_parent_map"] = build_parent_map(
            dar_asr["node_to_children"], dar_asr["root_node"])
    parent_map = dar_asr["_parent_map"]

    # Ancestral cDAV: dar_alt already at root (cannot be timed)
    if dar_states.get(dar_asr["root_node"]) == dar_alt:
        base["is_ancestral_cdav"] = True
        base["timing_confidence"] = "ancestral"
        return base

    # Find DAV gain branches (any gain of dar_alt, not just from dar_ref)
    gains = find_dav_gain_branches(
        dar_asr, dar_col, dar_alt, parent_map, dar_states)
    base["n_dar_gain_branches"] = len(gains)

    if not gains:
        base["timing_confidence"] = "no_gains_found"
        return base

    # Classify each gain branch
    for bkey in gains:
        parts = bkey.split("|", 1)
        if len(parts) != 2:
            continue
        parent_node, child_node = parts
        timing = classify_branch_timing(
            contact_states, contact_asr["node_to_children"],
            parent_node, child_node, contact_alt)
        base[f"n_{timing}"] += 1

    # Dominant timing and confidence level
    counts = {cat: base[f"n_{cat}"] for cat in TIMING_CATS}
    total  = sum(counts.values())
    if total > 0:
        dominant = max(counts, key=counts.get)
        base["dominant_timing"] = dominant
        majority_frac = counts[dominant] / total
        lc_branches = set(dar_asr.get("branches_lc", {}).keys())
        low_conf_gain = any(bkey in lc_branches for bkey in gains)
        base["timing_confidence"] = (
            "high" if majority_frac >= 0.70 and not low_conf_gain else "low")

    return base


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    cp  = pd.read_csv(COMP_PART)
    sig = cp[(cp["pagel_fdr"] <= PAGEL_ALPHA) | (cp["branch_cooccur_fdr"] <= BRANCH_ALPHA)].copy()
    print(f"  Pagel/branch significant pairs: {len(sig)}")

    print("Loading ASR maps...")
    with open(ASR_MAPS) as f:
        asr_data = json.load(f)
    print(f"  Genes in ASR maps: {len(asr_data)}")

    # Run analysis
    print("Running temporal ordering...")
    node_states_cache: dict = {}
    results = []
    for i, (_, row) in enumerate(sig.iterrows()):
        if i % 50 == 0:
            print(f"  {i+1}/{len(sig)}  {row['dar_gene']} {row['dar_aa_coord']}"
                  f"→{row['contact_gene']} {row['contact_refseq_pos']}")
        results.append(analyze_pair(row, asr_data, node_states_cache))

    df = pd.DataFrame(results)
    out_path = OUT_DIR / "timing_annotations_v2.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {len(df)} rows → {out_path.relative_to(ROOT)}")

    # ── Summary stats ──────────────────────────────────────────────────────────
    print("\n═══ Temporal Ordering Summary ═══")
    print(f"Total pairs analyzed:          {len(df)}")
    print(f"Ancestral cDAV (excluded):     {df['is_ancestral_cdav'].sum()}")
    print(f"No gain branches found:        {(df['timing_confidence'] == 'no_gains_found').sum()}")
    print(f"Missing gene in ASR:           {(df['timing_confidence'] == 'missing_gene').sum()}")

    testable = df[~df["is_ancestral_cdav"] & (df["n_dar_gain_branches"] > 0)].copy()
    print(f"Testable pairs (≥1 gain):      {len(testable)}")

    if len(testable) == 0:
        return

    print(f"\nBranch-event totals (all testable pairs):")
    for cat in TIMING_CATS:
        n = int(testable[f"n_{cat}"].sum())
        print(f"  {cat:<25} {n:>6}")

    print(f"\nDominant timing per pair:")
    dom = testable["dominant_timing"].value_counts()
    for k, v in dom.items():
        print(f"  {k:<25} {v:>5}  ({100*v/len(testable):.1f}%)")

    high_conf = testable[testable["timing_confidence"] == "high"].copy()
    print(f"\nHigh-confidence pairs (majority ≥70%, no lc):  {len(high_conf)}")
    if len(high_conf) > 0:
        n_first = int(high_conf["n_contact_first"].sum())
        n_after = int(high_conf["n_contact_after"].sum())
        n_cooc  = int(high_conf["n_co_occurring"].sum())
        n_total = n_first + n_after + n_cooc
        print(f"  contact_first events:  {n_first}")
        print(f"  contact_after events:  {n_after}")
        print(f"  co_occurring events:   {n_cooc}")
        if n_first + n_after > 0:
            btest = binomtest(n_first, n_first + n_after, 0.5, alternative="greater")
            print(f"  Binomial test (contact_first > contact_after): "
                  f"p = {btest.pvalue:.3e}")

    # Top hits (contact_first dominant, most gain branches)
    top = (testable[testable["dominant_timing"] == "contact_first"]
           .sort_values("n_dar_gain_branches", ascending=False)
           .head(20))
    print(f"\nTop contact_first pairs (n={len(top)}):")
    cols = ["dar_gene", "dar_aa_coord", "dar_alt_aa",
            "contact_gene", "contact_refseq_pos", "contact_alt_aa",
            "n_dar_gain_branches", "n_contact_first", "n_contact_after",
            "timing_confidence"]
    print(top[cols].to_string(index=False))


if __name__ == "__main__":
    main()
