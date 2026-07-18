"""
Temporal ordering analysis for cDAV compensation — Phase 2b.

For each Pagel/branch-significant DAR-contact pair, determines whether the
compensatory contact amino acid arose BEFORE the DAV (pre-adaptation), on the
SAME branch (co-occurring), or AFTER (rescue), using IQTree MAP ancestral states
pre-processed in ancestral_state_maps.json.

Physicochemical refinement: co_occurring and no_contact_change events are
further split using BLOSUM62 score between the contact's ancestral amino acid
and the Pagel-identified compensatory alt amino acid:

  co_occurring + BLOSUM(contact_parent → contact_alt) ≥ 1
      → permissive_background  (ancestral contact already physicochemically
                                similar to alt; structural environment was
                                accommodating before the formal change)
  co_occurring + BLOSUM < 1
      → co_adaptation          (genuine secondary site mutation arose with DAV)

  no_contact_change + BLOSUM(current_contact → contact_alt) ≥ 1
      → constitutively_permissive  (contact never needed to change; current
                                    AA already similar to compensatory state)
  no_contact_change + BLOSUM < 1
      → no_contact_change          (contact genuinely incompatible with alt)

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
import sys
import numpy as np
import pandas as pd
from collections import deque
from pathlib import Path
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from utils.variant_record import blosum62 as _blosum62_fn, miyata_distance as _miyata_fn

ASR_MAPS  = ROOT / "data" / "phylo" / "ancestral_state_maps.json"
COMP_PART = ROOT / "results" / "structural" / "compensatory_partners.csv"
OUT_DIR   = ROOT / "results" / "phylo"

PAGEL_ALPHA  = 0.10
BRANCH_ALPHA = 0.10

# Physicochemical similarity threshold for refined timing classification.
# BLOSUM62 ≥ 1 = positive (conservative) substitution: the two AAs are
# evolutionarily substitutable with no fitness cost. Anything below (including
# 0) is considered physicochemically distinct for this analysis.
BLOSUM_PERMISSIVE_THRESHOLD = 1
MIYATA_PERMISSIVE_THRESHOLD = 1.0   # stored as supporting metric

TIMING_CATS = ["contact_first", "co_occurring", "contact_after", "no_contact_change"]
REFINED_CATS = [
    "contact_first",
    "permissive_background",   # was co_occurring; ancestral contact ~ alt
    "co_adaptation",           # was co_occurring; ancestral contact != alt
    "constitutively_permissive",  # was no_contact_change; current contact ~ alt
    "contact_after",
    "no_contact_change",       # no change AND contact != alt
]


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
    Leaf nodes use observed states; internal nodes use tree traversal from
    root_states + branch substitutions. None indicates a gap/ambiguous root.
    """
    root = gene_asr["root_node"]
    leaf_nodes = set(gene_asr.get("leaf_nodes", []))
    leaf_states = gene_asr.get("leaf_states", {})

    states: dict = {}
    states[root] = gene_asr["root_states"].get(site_str)

    q: deque = deque([root])
    while q:
        parent = q.popleft()
        parent_state = states[parent]
        for child in gene_asr["node_to_children"].get(parent, []):
            if child in leaf_nodes:
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
    Return list of branch keys "parent|child" where ANY → alt_aa transition
    occurs (i.e., any gain of the pathogenic amino acid regardless of starting
    state). Broader than ref→alt: non-human cDAV species often acquired the
    pathogenic amino acid from a non-human ancestral state.
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


# ── Physicochemical helpers ────────────────────────────────────────────────────

def phys_metrics(aa1: str | None, aa2: str | None) -> dict:
    """
    Return BLOSUM62 score and Miyata distance between two amino acids.
    Returns NaN values if either AA is None, a gap, or not in the matrix.
    """
    if not aa1 or not aa2 or aa1 == "-" or aa2 == "-":
        return {"blosum": float("nan"), "miyata": float("nan")}
    if aa1 == aa2:
        return {"blosum": _blosum62_fn(aa1, aa2) or 0, "miyata": 0.0}
    b = _blosum62_fn(aa1, aa2)
    m = _miyata_fn(aa1, aa2)
    return {
        "blosum": float(b) if b is not None else float("nan"),
        "miyata": float(m) if m is not None else float("nan"),
    }


def refine_timing(timing: str, blosum: float, threshold: int = BLOSUM_PERMISSIVE_THRESHOLD) -> str:
    """
    Apply physicochemical refinement to a raw timing class.
    Only co_occurring and no_contact_change are split; others are passed through.
    """
    if np.isnan(blosum):
        return timing   # cannot refine without metrics
    if timing == "co_occurring":
        return "permissive_background" if blosum >= threshold else "co_adaptation"
    if timing == "no_contact_change":
        return "constitutively_permissive" if blosum >= threshold else "no_contact_change"
    return timing


# ── Per-branch timing classification ──────────────────────────────────────────

def classify_branch_timing(
    contact_node_states: dict,
    node_to_children: dict,
    parent_node: str,
    child_node: str,
    contact_alt: str,
) -> tuple[str, str | None]:
    """
    Classify timing of contact_alt appearance relative to a single DAV gain
    branch (parent_node → child_node).

    Returns (timing_class, contact_parent_aa):
      contact_first  — parent already has contact_alt
      co_occurring   — contact_alt appears on this same branch
      contact_after  — contact_alt appears in a descendant branch (rescue)
      no_contact_change — contact_alt never appears in this clade

    contact_parent_aa is the AA at the contact position in parent_node; used
    downstream to compute physicochemical distance to contact_alt.
    """
    parent_contact = contact_node_states.get(parent_node)
    child_contact  = contact_node_states.get(child_node)

    if parent_contact == contact_alt:
        return "contact_first", parent_contact
    if child_contact == contact_alt:
        return "co_occurring", parent_contact
    clade_states = get_clade_states(node_to_children, contact_node_states, child_node)
    if contact_alt in clade_states:
        return "contact_after", parent_contact
    return "no_contact_change", parent_contact


# ── Per-pair analysis ──────────────────────────────────────────────────────────

def analyze_pair(row: pd.Series, asr_data: dict, node_states_cache: dict) -> dict:
    dar_gene     = row["dar_gene"]
    dar_col      = str(int(row["dar_aa_coord"]))
    dar_ref      = row["dar_ref_aa"]
    dar_alt      = row["dar_alt_aa"]
    contact_gene = row["contact_gene"]
    contact_col  = str(int(row["contact_refseq_pos"]))
    contact_alt  = row["contact_alt_aa"]

    base = {
        "ann_id":                  row.get("ann_id"),
        "dar_gene":                dar_gene,
        "dar_aa_coord":            row["dar_aa_coord"],
        "dar_ref_aa":              dar_ref,
        "dar_alt_aa":              dar_alt,
        "contact_gene":            contact_gene,
        "contact_refseq_pos":      row["contact_refseq_pos"],
        "contact_human_aa":        row["contact_human_aa"],
        "contact_alt_aa":          contact_alt,
        "is_ancestral_cdav":       False,
        "n_dar_gain_branches":     0,
        # Raw timing counts
        "n_contact_first":         0,
        "n_co_occurring":          0,
        "n_contact_after":         0,
        "n_no_contact_change":     0,
        # Refined timing counts
        "n_permissive_background": 0,
        "n_co_adaptation":         0,
        "n_constitutively_permissive": 0,
        # Physicochemical metrics (median across gain branches)
        "contact_parent_aa_modal": np.nan,
        "phys_blosum_median":      np.nan,
        "phys_miyata_median":      np.nan,
        # Summary
        "dominant_timing":         np.nan,
        "dominant_refined_timing": np.nan,
        "timing_confidence":       "insufficient",
        "dar_origin_age_mya":      np.nan,
    }

    if dar_gene not in asr_data or contact_gene not in asr_data:
        base["timing_confidence"] = "missing_gene"
        return base

    dar_asr     = asr_data[dar_gene]
    contact_asr = asr_data[contact_gene]

    g = node_states_cache.setdefault(dar_gene, {})
    if dar_col not in g:
        g[dar_col] = compute_node_states(dar_asr, dar_col)
    dar_states = g[dar_col]

    cg = node_states_cache.setdefault(contact_gene, {})
    if contact_col not in cg:
        cg[contact_col] = compute_node_states(contact_asr, contact_col)
    contact_states = cg[contact_col]

    if "_parent_map" not in dar_asr:
        dar_asr["_parent_map"] = build_parent_map(
            dar_asr["node_to_children"], dar_asr["root_node"])
    parent_map = dar_asr["_parent_map"]

    if dar_states.get(dar_asr["root_node"]) == dar_alt:
        base["is_ancestral_cdav"] = True
        base["timing_confidence"] = "ancestral"
        return base

    gains = find_dav_gain_branches(dar_asr, dar_col, dar_alt, parent_map, dar_states)
    base["n_dar_gain_branches"] = len(gains)

    if not gains:
        base["timing_confidence"] = "no_gains_found"
        return base

    blosum_vals: list[float] = []
    miyata_vals: list[float] = []
    parent_aas:  list[str]   = []

    for bkey in gains:
        parts = bkey.split("|", 1)
        if len(parts) != 2:
            continue
        parent_node, child_node = parts
        timing, contact_parent_aa = classify_branch_timing(
            contact_states, contact_asr["node_to_children"],
            parent_node, child_node, contact_alt)

        base[f"n_{timing}"] += 1

        # Physicochemical metrics: compare contact_parent_aa to contact_alt
        pm = phys_metrics(contact_parent_aa, contact_alt)
        blosum_vals.append(pm["blosum"])
        miyata_vals.append(pm["miyata"])
        if contact_parent_aa:
            parent_aas.append(contact_parent_aa)

        # Refined timing using the branch-level physicochemical distance.
        # Only increment the refined counter when the class actually changes;
        # categories that pass through unchanged (contact_first, contact_after,
        # no_contact_change with BLOSUM < threshold) are already counted in the
        # raw n_{timing} increment above.
        refined = refine_timing(timing, pm["blosum"])
        if refined != timing:
            base[f"n_{refined}"] += 1

    # Store median physicochemical metrics across all gain branches
    valid_b = [v for v in blosum_vals if not np.isnan(v)]
    valid_m = [v for v in miyata_vals if not np.isnan(v)]
    if valid_b:
        base["phys_blosum_median"] = float(np.median(valid_b))
    if valid_m:
        base["phys_miyata_median"] = float(np.median(valid_m))
    if parent_aas:
        from collections import Counter
        base["contact_parent_aa_modal"] = Counter(parent_aas).most_common(1)[0][0]

    # Raw dominant timing
    raw_counts = {cat: base[f"n_{cat}"] for cat in TIMING_CATS}
    raw_total  = sum(raw_counts.values())
    if raw_total > 0:
        base["dominant_timing"] = max(raw_counts, key=raw_counts.get)

    # Refined dominant timing
    refined_counts = {cat: base.get(f"n_{cat}", 0) for cat in REFINED_CATS}
    refined_total  = sum(refined_counts.values())
    if refined_total > 0:
        dominant_refined = max(refined_counts, key=refined_counts.get)
        base["dominant_refined_timing"] = dominant_refined
        majority_frac = refined_counts[dominant_refined] / refined_total
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

    print(f"\nRaw branch-event totals (all testable pairs):")
    for cat in TIMING_CATS:
        n = int(testable[f"n_{cat}"].sum())
        print(f"  {cat:<25} {n:>6}")

    print(f"\nRefined branch-event totals (mutually exclusive, sum = total branch events):")
    for cat in REFINED_CATS:
        col = f"n_{cat}"
        if cat == "no_contact_change":
            # Remaining after subtracting constitutively_permissive subset
            n_total = int(testable["n_no_contact_change"].sum()) if "n_no_contact_change" in testable.columns else 0
            n_cp    = int(testable["n_constitutively_permissive"].sum()) if "n_constitutively_permissive" in testable.columns else 0
            n = n_total - n_cp
        else:
            n = int(testable[col].sum()) if col in testable.columns else 0
        print(f"  {cat:<30} {n:>6}")

    print(f"\nPhysicochemical metrics (median across all gain branches per pair):")
    print(f"  BLOSUM62 median (contact_parent→alt):  {testable['phys_blosum_median'].median():.1f}")
    print(f"  Miyata  median (contact_parent→alt):   {testable['phys_miyata_median'].median():.2f}")

    print(f"\nDominant refined timing per pair:")
    dom = testable["dominant_refined_timing"].value_counts()
    for k, v in dom.items():
        print(f"  {k:<30} {v:>5}  ({100*v/len(testable):.1f}%)")

    high_conf = testable[testable["timing_confidence"] == "high"].copy()
    print(f"\nHigh-confidence pairs (majority ≥70%, no lc):  {len(high_conf)}")
    if len(high_conf) > 0:
        # Raw
        n_first = int(high_conf["n_contact_first"].sum())
        n_after = int(high_conf["n_contact_after"].sum())
        n_cooc  = int(high_conf["n_co_occurring"].sum())
        print(f"  [raw] contact_first={n_first}  co_occurring={n_cooc}  contact_after={n_after}")
        if n_first + n_after > 0:
            bt = binomtest(n_first, n_first + n_after, 0.5, alternative="greater")
            print(f"  Binomial test contact_first > contact_after: p = {bt.pvalue:.3e}")

        # Refined
        n_pb  = int(high_conf.get("n_permissive_background", pd.Series([0])).sum())
        n_coa = int(high_conf.get("n_co_adaptation", pd.Series([0])).sum())
        n_cp  = int(high_conf.get("n_constitutively_permissive", pd.Series([0])).sum())
        print(f"\n  [refined] permissive_background={n_pb}  co_adaptation={n_coa}"
              f"  constitutively_permissive={n_cp}")
        print(f"  → Permissive signal (contact_first + permissive_background + constitutively_permissive):"
              f" {n_first + n_pb + n_cp}")
        print(f"  → Secondary mutation signal (co_adaptation): {n_coa}")

    # Top contact_first pairs
    top = (testable[testable["dominant_refined_timing"].isin(
               ["contact_first", "permissive_background", "co_adaptation"])]
           .sort_values(["dominant_refined_timing", "n_dar_gain_branches"], ascending=[True, False])
           .head(20))
    if len(top) > 0:
        print(f"\nTop pairs with directional refined timing (n={len(top)}):")
        cols = ["dar_gene", "dar_aa_coord", "dar_alt_aa",
                "contact_gene", "contact_refseq_pos", "contact_alt_aa",
                "dominant_refined_timing",
                "n_dar_gain_branches", "n_contact_first",
                "n_permissive_background", "n_co_adaptation",
                "phys_blosum_median", "phys_miyata_median",
                "timing_confidence"]
        present = [c for c in cols if c in top.columns]
        print(top[present].to_string(index=False))


if __name__ == "__main__":
    main()
