"""
src/phylo/_timing_core.py

Shared per-branch temporal ordering logic used by both
06_temporal_ordering.py (Pagel pairs) and
06b_temporal_ordering_mitonuc.py (mito-nuclear Pyvolve pairs).
"""
import sys
from collections import Counter, deque
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from utils.variant_record import blosum62 as _blosum62_fn, miyata_distance as _miyata_fn

BLOSUM_PERMISSIVE_THRESHOLD = 1

TIMING_CATS = ["contact_first", "co_occurring", "contact_after", "no_contact_change"]
REFINED_CATS = [
    "contact_first",
    "permissive_background",
    "co_adaptation",
    "constitutively_permissive",
    "contact_after",
    "no_contact_change",
]


def build_parent_map(node_to_children: dict, root: str) -> dict:
    parent_map: dict = {}
    q: deque = deque([root])
    while q:
        node = q.popleft()
        for child in node_to_children.get(node, []):
            parent_map[child] = node
            q.append(child)
    return parent_map


def compute_node_states(gene_asr: dict, site_str: str) -> dict:
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
                            parent_map: dict, dar_node_states: dict) -> list:
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


def get_clade_states(node_to_children: dict, node_states: dict, start_node: str) -> set:
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


def phys_metrics(aa1, aa2) -> dict:
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
    if np.isnan(blosum):
        return timing
    if timing == "co_occurring":
        return "permissive_background" if blosum >= threshold else "co_adaptation"
    if timing == "no_contact_change":
        return "constitutively_permissive" if blosum >= threshold else "no_contact_change"
    return timing


def classify_branch_timing(contact_node_states, node_to_children, parent_node,
                            child_node, contact_alt):
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


def analyze_pair(row: pd.Series, asr_data: dict, node_states_cache: dict) -> dict:
    dar_gene     = row["dar_gene"]
    dar_col      = str(int(row["dar_aa_coord"]))
    dar_ref      = row.get("dar_ref_aa")
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
        "contact_human_aa":        row.get("contact_human_aa"),
        "contact_alt_aa":          contact_alt,
        "is_ancestral_cdav":       False,
        "n_dar_gain_branches":     0,
        "n_contact_first":         0,
        "n_co_occurring":          0,
        "n_contact_after":         0,
        "n_no_contact_change":     0,
        "n_permissive_background": 0,
        "n_co_adaptation":         0,
        "n_constitutively_permissive": 0,
        "contact_parent_aa_modal": np.nan,
        "phys_blosum_median":      np.nan,
        "phys_miyata_median":      np.nan,
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

    blosum_vals: list = []
    miyata_vals: list = []
    parent_aas:  list = []

    for bkey in gains:
        parts = bkey.split("|", 1)
        if len(parts) != 2:
            continue
        parent_node, child_node = parts
        timing, contact_parent_aa = classify_branch_timing(
            contact_states, contact_asr["node_to_children"],
            parent_node, child_node, contact_alt)

        base[f"n_{timing}"] += 1

        pm = phys_metrics(contact_parent_aa, contact_alt)
        blosum_vals.append(pm["blosum"])
        miyata_vals.append(pm["miyata"])
        if contact_parent_aa:
            parent_aas.append(contact_parent_aa)

        refined = refine_timing(timing, pm["blosum"])
        if refined != timing:
            base[f"n_{refined}"] += 1

    valid_b = [v for v in blosum_vals if not np.isnan(v)]
    valid_m = [v for v in miyata_vals if not np.isnan(v)]
    if valid_b:
        base["phys_blosum_median"] = float(np.median(valid_b))
    if valid_m:
        base["phys_miyata_median"] = float(np.median(valid_m))
    if parent_aas:
        base["contact_parent_aa_modal"] = Counter(parent_aas).most_common(1)[0][0]

    raw_counts = {cat: base[f"n_{cat}"] for cat in TIMING_CATS}
    raw_total  = sum(raw_counts.values())
    if raw_total > 0:
        base["dominant_timing"] = max(raw_counts, key=raw_counts.get)

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
