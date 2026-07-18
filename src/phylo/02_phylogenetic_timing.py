"""
src/phylo/02_phylogenetic_timing.py

Annotates all tested compensatory pairs with phylogenetic timing:
  - Is the cDAV ancestral (alt_aa is the mammalian ancestral state) or derived?
  - For derived cDAVs: on how many branches did the alt_aa arise (gains)?
  - For ancestral cDAVs: on how many branches was the alt_aa lost?
  - For each gain/origin branch: did the compensatory contact substitution
    occur before (contact_first), concurrently (co_occurring), or after
    (contact_after)?
  - Age of origin node estimated from MamPhy branch lengths (Mya).

Cross-gene node matching is done by subtended leaf species rather than
internal node IDs (which differ between IQTree runs per gene).

Inputs:
  results/structural/all_tested_pairs.csv   -- all tested pairs
  data/phylo/ancestral_state_maps.json       -- from 01_parse_ancestral_states.py
  data/phylo/species_tree/mammaltree_crossgenome.nwk  -- time-calibrated tree

Output:
  results/phylo/timing_annotations.csv

Run from project root inside the Docker container:
    python src/phylo/02_phylogenetic_timing.py
"""

import csv
import json
import sys
from collections import Counter
from pathlib import Path

from Bio import Phylo

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parents[2]
DATA_DIR    = ROOT / "data"
RESULTS_DIR = ROOT / "results"

PAIRS_CSV  = RESULTS_DIR / "structural" / "all_tested_pairs.csv"
ANC_MAPS   = DATA_DIR / "phylo" / "ancestral_state_maps.json"
MAMM_TREE  = DATA_DIR / "phylo" / "species_tree" / "mammaltree_crossgenome.nwk"
OUT_DIR    = RESULTS_DIR / "phylo"
OUT_CSV    = OUT_DIR / "timing_annotations.csv"

ANCESTRAL_THRESHOLD = 0.5   # leaf-majority fallback when root state is ambiguous
LC_CONFIDENCE_LABEL = "low_confidence"   # label for 0.50–0.80 posterior tier

# ── Physicochemical property tables ───────────────────────────────────────────
# Charge groups
_POS_CHARGED  = {"R", "K", "H"}
_NEG_CHARGED  = {"D", "E"}
_POLAR_UNCH   = {"S", "T", "N", "Q", "Y", "C"}
_NONPOLAR     = {"A", "V", "L", "I", "M", "F", "W", "P", "G"}

# Van der Waals volumes (Å³) — Pontius et al. approximations
_VDW_VOL = {
    "G": 60,  "A": 88,  "V": 140, "L": 166, "I": 166,
    "P": 112, "F": 190, "W": 227, "M": 162, "S": 89,
    "T": 116, "C": 108, "Y": 193, "H": 153, "D": 111,
    "E": 138, "N": 114, "Q": 143, "K": 168, "R": 173,
}


# ── Tree age helpers ───────────────────────────────────────────────────────────

def load_timetree(path: Path):
    """Load MamPhy ultrametric tree. Branch lengths are in Mya."""
    if not path.exists():
        return None
    return Phylo.read(str(path), "newick")


def node_age_mya(tree, species_set: set[str]) -> float | None:
    """
    Estimate the age of the DAR origin node using the MamPhy ultrametric tree.
    Age = total tree depth - distance from root to MRCA of species_set.

    For internal-node origins (>=2 species): returns MRCA age directly.
    For leaf-node origins (1 species): the DAR arose on the terminal branch
    leading to that species; age is reported as 0.0 Mya (present-day lineage).
    Returns None if the species is absent from the timetree or tree is None.
    """
    if tree is None or not species_set:
        return None
    tip_names = {c.name for c in tree.get_terminals()}
    present = [s for s in species_set if s in tip_names]
    if not present:
        return None
    if len(present) == 1:
        # Leaf-origin: DAR arose on the terminal branch of a single extant species
        return 0.0
    try:
        mrca = tree.common_ancestor(*[{"name": s} for s in present])
        root_to_mrca = tree.distance(tree.root, mrca)
        root_to_leaf = tree.distance(tree.root, tree.get_terminals()[0])
        return round(root_to_leaf - root_to_mrca, 2)
    except Exception:
        return None


# ── Subtended species ──────────────────────────────────────────────────────────

def subtended_species(node: str, node_to_children: dict, leaf_nodes: set) -> frozenset[str]:
    """BFS from node; return frozenset of all leaf nodes in the subtree."""
    queue  = [node]
    leaves = set()
    while queue:
        n = queue.pop()
        if n in leaf_nodes:
            leaves.add(n)
        for child in node_to_children.get(n, []):
            queue.append(child)
    return frozenset(leaves)


def build_subtree_map(gene_map: dict) -> dict[str, frozenset[str]]:
    """
    For every node in the gene's topology, pre-compute the frozenset of
    subtended leaf species. Used for cross-gene node matching.
    """
    node_to_children = gene_map.get("node_to_children", {})
    leaf_nodes       = set(gene_map.get("leaf_nodes", []))
    all_nodes        = set(node_to_children.keys()) | leaf_nodes
    return {n: subtended_species(n, node_to_children, leaf_nodes) for n in all_nodes}


def match_node_by_species(
    dar_species_set: frozenset[str],
    contact_subtree_map: dict[str, frozenset[str]],
) -> str | None:
    """
    Find the contact-gene node whose subtended species best matches
    dar_species_set. Best = largest intersection / smallest symmetric difference.
    Returns node label or None if no reasonable match (< 2 shared species).
    """
    best_node  = None
    best_score = -1
    for node, spp in contact_subtree_map.items():
        shared = len(dar_species_set & spp)
        if shared < 2:
            continue
        # Jaccard similarity
        union = len(dar_species_set | spp)
        score = shared / union if union else 0
        if score > best_score:
            best_score = score
            best_node  = node
    return best_node


def build_parent_map(node_to_children: dict[str, list[str]]) -> dict[str, str]:
    """Return child node -> parent node for a parsed IQTree topology."""
    parents = {}
    for parent, children in node_to_children.items():
        for child in children:
            parents[child] = parent
    return parents


def branch_change_for_pos(branches: dict, parent: str, child: str, pos_str: str):
    """Return a branch change entry for a position, if IQTree reported one."""
    return branches.get(f"{parent}|{child}", {}).get(pos_str)


def infer_node_states_for_pos(gene_map: dict, aa_pos: int, branches_key: str = "branches") -> dict[str, str]:
    """
    Reconstruct node states for one amino-acid position by propagating the root
    state through the tree and applying branch changes.

    The ASR map stores root states and branch deltas, not every internal-node
    state. Timing needs the state at the parent and child of the DAR-origin
    branch; using deltas alone misses pre-existing permissive backgrounds.
    """
    pos_str = str(aa_pos)
    root_node = gene_map.get("root_node")
    root_state = gene_map.get("root_states", {}).get(pos_str)
    if not root_node or not root_state:
        return {}

    node_to_children = gene_map.get("node_to_children", {})
    branches = gene_map.get(branches_key, {})
    states = {root_node: root_state}
    queue = [root_node]
    while queue:
        parent = queue.pop(0)
        parent_state = states.get(parent)
        for child in node_to_children.get(parent, []):
            entry = branch_change_for_pos(branches, parent, child, pos_str)
            states[child] = entry[1] if entry else parent_state
            queue.append(child)
    return states


# ── Ancestral cDAV detection ───────────────────────────────────────────────────

def is_ancestral_cdav(gene_map: dict, aa_pos: int, alt_aa: str) -> tuple[bool, str]:
    """
    Returns (is_ancestral, method) where method describes how the call was made.

    Strategy (in priority order):
    1. Root state from IQTree reconstruction (PP >= 0.80): if the root node has a
       high-confidence reconstructed state at this position, use it directly.
       This avoids bias from species-rich clades (e.g. Rodentia) inflating leaf counts.
    2. Leaf majority vote (fallback): if root state is absent or ambiguous (no HC call),
       fall back to the fraction of leaf species carrying alt_aa.

    Returns:
      (True,  "root_state")    -- root node reconstructed state == alt_aa (HC)
      (False, "root_state")    -- root node reconstructed state != alt_aa (HC)
      (True,  "leaf_majority") -- root ambiguous; >=50% leaves carry alt_aa
      (False, "leaf_majority") -- root ambiguous; <50% leaves carry alt_aa
      (False, "no_data")       -- no leaf data at this position
    """
    pos_str   = str(aa_pos)
    root_states = gene_map.get("root_states", {})

    # Strategy 1: direct root state (preferred)
    root_aa = root_states.get(pos_str)
    if root_aa:
        return (root_aa == alt_aa), "root_state"

    # Strategy 2: leaf majority fallback
    leaf_states = gene_map.get("leaf_states", {})
    total = 0
    n_alt = 0
    for states in leaf_states.values():
        s = states.get(pos_str)
        if s:
            total += 1
            if s == alt_aa:
                n_alt += 1
    if total == 0:
        return False, "no_data"
    return (n_alt / total) >= ANCESTRAL_THRESHOLD, "leaf_majority"


def find_dar_gains(gene_map: dict, aa_pos: int, alt_aa: str) -> list[tuple[str, str]]:
    """Branches where alt_aa arose: parent != alt_aa, child == alt_aa."""
    pos_str = str(aa_pos)
    origins = []
    for branch_id, changes in gene_map.get("branches", {}).items():
        entry = changes.get(pos_str)
        if entry and entry[1] == alt_aa and entry[0] != alt_aa:
            parent, child = branch_id.split("|", 1)
            origins.append((parent, child))
    return origins


def find_dar_losses(gene_map: dict, aa_pos: int, alt_aa: str) -> list[tuple[str, str]]:
    """Branches where alt_aa was lost: parent == alt_aa, child != alt_aa."""
    pos_str = str(aa_pos)
    losses = []
    for branch_id, changes in gene_map.get("branches", {}).items():
        entry = changes.get(pos_str)
        if entry and entry[0] == alt_aa and entry[1] != alt_aa:
            parent, child = branch_id.split("|", 1)
            losses.append((parent, child))
    return losses


# ── Physicochemical complementarity ───────────────────────────────────────────

def _charge_group(aa: str) -> str:
    if aa in _POS_CHARGED:  return "positive"
    if aa in _NEG_CHARGED:  return "negative"
    if aa in _POLAR_UNCH:   return "polar"
    if aa in _NONPOLAR:     return "nonpolar"
    return "unknown"


def physicochemical_complementarity(
    dar_ref_aa:    str,
    dar_alt_aa:    str,
    contact_human: str,
    contact_alt:   str,
) -> str:
    """
    Classify the physicochemical relationship between the DAR substitution and
    the contact substitution as a compensatory mechanism type.

    Returns one of:
      'charge_reversal'     -- DAR changes charge sign, contact reverses to compensate
                               (e.g. DAR: K→E introduces negative; contact: D→K restores salt bridge)
      'charge_rescue'       -- DAR neutralises a charge; contact changes to restore electrostatics
      'volume_swap'         -- DAR and contact show reciprocal volume changes (one grows, one shrinks)
                               indicative of packing compensation
      'polarity_swap'       -- one changes from polar to nonpolar (or vice versa), other reciprocates
      'same_direction'      -- both changes shift properties in the same direction (less specific)
      'unclassified'        -- insufficient data or no clear pattern
    """
    if not all([dar_ref_aa, dar_alt_aa, contact_human, contact_alt]):
        return "unclassified"
    if dar_alt_aa == dar_ref_aa or contact_alt == contact_human:
        return "unclassified"

    dar_ref_charge  = _charge_group(dar_ref_aa)
    dar_alt_charge  = _charge_group(dar_alt_aa)
    cont_ref_charge = _charge_group(contact_human)
    cont_alt_charge = _charge_group(contact_alt)

    # Charge reversal: DAR switches +/-, contact switches the opposite way
    if ({dar_ref_charge, dar_alt_charge} == {"positive", "negative"} and
            {cont_ref_charge, cont_alt_charge} == {"positive", "negative"} and
            dar_alt_charge != cont_alt_charge):
        return "charge_reversal"

    # Charge rescue: DAR loses/gains charge, contact compensates
    dar_charge_change  = dar_ref_charge  != dar_alt_charge
    cont_charge_change = cont_ref_charge != cont_alt_charge
    if dar_charge_change and cont_charge_change:
        return "charge_rescue"

    # Volume swap: reciprocal change in side-chain size
    dar_vol_ref  = _VDW_VOL.get(dar_ref_aa,   0)
    dar_vol_alt  = _VDW_VOL.get(dar_alt_aa,   0)
    cont_vol_ref = _VDW_VOL.get(contact_human, 0)
    cont_vol_alt = _VDW_VOL.get(contact_alt,   0)
    dar_vol_delta  = dar_vol_alt  - dar_vol_ref
    cont_vol_delta = cont_vol_alt - cont_vol_ref
    # Reciprocal: one grows, other shrinks, both by > 20 Å³
    if (abs(dar_vol_delta) > 20 and abs(cont_vol_delta) > 20 and
            dar_vol_delta * cont_vol_delta < 0):
        return "volume_swap"

    # Polarity swap
    dar_polar_ref  = dar_ref_aa  in _POLAR_UNCH or dar_ref_aa  in _POS_CHARGED or dar_ref_aa  in _NEG_CHARGED
    dar_polar_alt  = dar_alt_aa  in _POLAR_UNCH or dar_alt_aa  in _POS_CHARGED or dar_alt_aa  in _NEG_CHARGED
    cont_polar_ref = contact_human in _POLAR_UNCH or contact_human in _POS_CHARGED or contact_human in _NEG_CHARGED
    cont_polar_alt = contact_alt  in _POLAR_UNCH or contact_alt  in _POS_CHARGED or contact_alt  in _NEG_CHARGED
    if (dar_polar_ref != dar_polar_alt) and (cont_polar_ref != cont_polar_alt) and (dar_polar_alt != cont_polar_alt):
        return "polarity_swap"

    # Same direction (both get more polar, or both get larger, etc.)
    if (dar_polar_ref != dar_polar_alt) and (cont_polar_ref != cont_polar_alt):
        return "same_direction"

    return "unclassified"


# ── Contact state characterization ────────────────────────────────────────────

def contact_state_in_cdav_spp(
    contact_gene_map: dict,
    contact_pos: int,
    contact_alt_aa: str,
    contact_human_aa: str,
    cdav_species: set[str],
) -> str:
    """
    For species that carry the cDAV alt_aa (or in the ancestral cDAV subtree),
    characterize the contact residue state:
      'conserved_human'  -- contact residue == human AA in all cDAV species
      'conserved_alt'    -- contact alt_aa fixed in all cDAV species
      'polymorphic'      -- mixed; both or other states present
      'unknown'          -- no data
    """
    if not cdav_species:
        return "unknown"
    leaf_states = contact_gene_map.get("leaf_states", {})
    pos_str = str(contact_pos)
    states_in_cdav = []
    for sp in cdav_species:
        s = leaf_states.get(sp, {}).get(pos_str)
        if s:
            states_in_cdav.append(s)
    if not states_in_cdav:
        return "unknown"
    unique = set(states_in_cdav)
    if unique == {contact_human_aa}:
        return "conserved_human"
    if unique == {contact_alt_aa}:
        return "conserved_alt"
    return "polymorphic"


# ── Timing per origin ──────────────────────────────────────────────────────────

def timing_for_origin(
    dar_origin_child: str,
    dar_subtree_map:  dict[str, frozenset[str]],
    contact_gene_map: dict,
    contact_subtree_map: dict[str, frozenset[str]],
    contact_pos: int,
    contact_alt_aa: str,
) -> str:
    """
    Determine whether the contact substitution occurred relative to the cDAV
    origin at dar_origin_child. Matches nodes across genes by subtended species.

    Returns one of: 'contact_first', 'co_occurring', 'contact_after', 'no_contact_change'
    """
    dar_spp = dar_subtree_map.get(dar_origin_child, frozenset())
    pos_str = str(contact_pos)
    c_branches = contact_gene_map.get("branches", {})
    contact_states = infer_node_states_for_pos(contact_gene_map, contact_pos)
    contact_parent_map = build_parent_map(contact_gene_map.get("node_to_children", {}))

    # Find the contact-gene node that corresponds to dar_origin_child
    # by matching subtended leaf species
    matched_contact_node = match_node_by_species(dar_spp, contact_subtree_map)

    # ── Check state on the matched DAR-origin branch ──────────────────────────
    # parent already has contact_alt_aa: permissive background/contact_first
    # parent lacks it and child gains it: co-occurring on the DAR-origin branch
    if matched_contact_node:
        matched_parent = contact_parent_map.get(matched_contact_node)
        parent_state = contact_states.get(matched_parent)
        child_state = contact_states.get(matched_contact_node)
        if parent_state == contact_alt_aa:
            return "contact_first"
        if parent_state and parent_state != contact_alt_aa and child_state == contact_alt_aa:
            return "co_occurring"

    # ── Check contact_after: contact alt_aa arises in the DAR subtree ──────────
    for branch_id, changes in c_branches.items():
        parent, child = branch_id.split("|", 1)
        child_spp = contact_subtree_map.get(child, frozenset())
        if child_spp and child_spp.issubset(dar_spp) and child_spp != dar_spp:
            entry = changes.get(pos_str)
            if entry and entry[0] != contact_alt_aa and entry[1] == contact_alt_aa:
                return "contact_after"

    return "no_contact_change"


def directional_class(is_ancestral_cdav: bool, timing_counts: Counter) -> str:
    """
    Interpret timing counts for the biological mechanism question.

    Directionality is only meaningful for derived cDAV gains. If the disease
    amino acid is reconstructed as ancestral, there is no DAV-gain event whose
    responder/permissive order can be inferred.
    """
    if is_ancestral_cdav:
        return "ancestral_cdav_not_directional"
    first = timing_counts["contact_first"]
    co = timing_counts["co_occurring"]
    after = timing_counts["contact_after"]
    none = timing_counts["no_contact_change"]
    directional = first + co + after
    if directional == 0:
        return "no_detected_partner_change" if none else "ambiguous_no_origin"
    if first and not co and not after:
        return "permissive_background"
    if after and not first and not co:
        return "responding_secondary"
    if co and not first and not after:
        return "co_occurring_unresolved"
    return "mixed_timing"


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    for required in (PAIRS_CSV, ANC_MAPS):
        if not required.exists():
            print(f"ERROR: {required} not found.")
            sys.exit(1)

    timetree = load_timetree(MAMM_TREE)
    if timetree is None:
        print(f"WARN: MamPhy tree not found at {MAMM_TREE} — age estimation skipped.")
    else:
        print(f"Loaded MamPhy tree: {sum(1 for _ in timetree.get_terminals())} tips")

    print("Loading ancestral state maps...")
    with open(ANC_MAPS) as f:
        anc_maps: dict = json.load(f)
    print(f"  {len(anc_maps)} genes loaded.")

    # Pre-compute subtree maps (node → frozenset of subtended species) per gene
    print("Pre-computing subtree maps...")
    subtree_maps: dict[str, dict] = {
        gene: build_subtree_map(gmap) for gene, gmap in anc_maps.items()
    }
    print(f"  Done.\n")

    print("Loading tested pairs...")
    with open(PAIRS_CSV) as f:
        pairs = list(csv.DictReader(f))
    print(f"  {len(pairs)} pairs.\n")

    # Deduplicate by (ann_id, contact_gene, contact_pos, contact_alt_aa)
    seen: set = set()
    unique_pairs: list[dict] = []
    for row in pairs:
        key = (row["ann_id"], row.get("contact_gene"),
               row.get("contact_refseq_pos"), row.get("contact_alt_aa"))
        if key not in seen:
            seen.add(key)
            unique_pairs.append(row)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    out_fields = [
        "ann_id", "dar_gene", "dar_genome", "dar_aa_coord", "dar_ref_aa", "dar_alt_aa",
        "contact_gene", "contact_refseq_pos", "contact_human_aa", "contact_alt_aa",
        "tier",
        "is_ancestral_cdav", "ancestral_method",
        "confidence_tier",
        "n_dar_gain_branches", "n_dar_loss_branches",
        "n_contact_first", "n_co_occurring", "n_contact_after", "n_no_contact_change",
        "dominant_timing",
        "directional_class",
        "contact_state_in_cdav_spp",
        "physicochemical_type",
        "dar_origin_node",
        "dar_origin_age_mya",
    ]

    rows_out:  list[dict] = []
    processed  = 0
    skip_counts = Counter()

    for row in unique_pairs:
        dar_gene     = row["dar_gene"]
        contact_gene = row.get("contact_gene", "")
        if not contact_gene or contact_gene.startswith("Unknown("):
            skip_counts["no_contact_gene"] += 1
            continue

        try:
            dar_pos     = int(row["dar_aa_coord"])
            contact_pos = int(row["contact_refseq_pos"])
        except (ValueError, TypeError):
            skip_counts["bad_coords"] += 1
            continue

        dar_alt_aa     = row.get("dar_alt_aa", "")
        contact_alt_aa = row.get("contact_alt_aa", "")
        contact_human_aa = row.get("contact_human_aa", "")

        if dar_gene not in anc_maps or contact_gene not in anc_maps:
            skip_counts["missing_anc_map"] += 1
            continue

        dar_map          = anc_maps[dar_gene]
        contact_map      = anc_maps[contact_gene]
        dar_stmap        = subtree_maps[dar_gene]
        contact_stmap    = subtree_maps[contact_gene]
        node_to_children = dar_map.get("node_to_children", {})
        leaf_nodes       = set(dar_map.get("leaf_nodes", []))

        # ── Classify ancestral vs derived cDAV ────────────────────────────────
        ancestral, ancestral_method = is_ancestral_cdav(dar_map, dar_pos, dar_alt_aa)

        gains  = find_dar_gains(dar_map, dar_pos, dar_alt_aa)
        losses = find_dar_losses(dar_map, dar_pos, dar_alt_aa)

        # For timing, use gains for derived cDAVs, losses for ancestral cDAVs
        # (for ancestral cDAVs, "contact_first" means the contact changed before
        #  the alt_aa was lost, i.e., a permissive pre-adaptation)
        timing_branches = gains if not ancestral else losses
        confidence_tier = "high_confidence"

        if not timing_branches and not ancestral:
            # Derived cDAV with no HC gain branches — try low-confidence tier
            lc_map = {"branches": dar_map.get("branches_lc", {}),
                      "node_to_children": dar_map.get("node_to_children", {}),
                      "leaf_nodes": dar_map.get("leaf_nodes", [])}
            lc_gains = find_dar_gains(lc_map, dar_pos, dar_alt_aa)
            if lc_gains:
                timing_branches  = lc_gains
                confidence_tier  = LC_CONFIDENCE_LABEL
            else:
                # Truly unresolvable even at low confidence
                skip_counts["no_dar_origin_in_tree"] += 1
                continue

        # ── Timing per origin/loss branch ─────────────────────────────────────
        timing_counts    = Counter()
        first_event_node = timing_branches[0][1] if timing_branches else None
        cdav_species: set[str] = set()

        for parent_node, child_node in timing_branches:
            event_spp = dar_stmap.get(child_node, frozenset())
            cdav_species.update(event_spp)

            timing = timing_for_origin(
                child_node, dar_stmap,
                contact_map, contact_stmap,
                contact_pos, contact_alt_aa,
            )
            timing_counts[timing] += 1

        dominant = timing_counts.most_common(1)[0][0] if timing_counts else "no_contact_change"
        direction = directional_class(ancestral, timing_counts)

        # ── Contact state characterization ────────────────────────────────────
        c_state = contact_state_in_cdav_spp(
            contact_map, contact_pos, contact_alt_aa, contact_human_aa, cdav_species
        )

        # ── Age estimation from MamPhy tree ───────────────────────────────────
        age_mya = None
        if first_event_node and timetree is not None:
            origin_spp = dar_stmap.get(first_event_node, frozenset())
            age_mya = node_age_mya(timetree, origin_spp)

        dar_ref_aa = row.get("dar_ref_aa", "")
        pc_type = physicochemical_complementarity(
            dar_ref_aa, dar_alt_aa, contact_human_aa, contact_alt_aa
        )

        rows_out.append({
            "ann_id":                  row["ann_id"],
            "dar_gene":                dar_gene,
            "dar_genome":              row.get("dar_genome", ""),
            "dar_aa_coord":            dar_pos,
            "dar_ref_aa":              dar_ref_aa,
            "dar_alt_aa":              dar_alt_aa,
            "contact_gene":            contact_gene,
            "contact_refseq_pos":      contact_pos,
            "contact_human_aa":        contact_human_aa,
            "contact_alt_aa":          contact_alt_aa,
            "tier":                    row.get("tier", ""),
            "is_ancestral_cdav":       ancestral,
            "ancestral_method":        ancestral_method,
            "confidence_tier":         confidence_tier,
            "n_dar_gain_branches":     len(gains),
            "n_dar_loss_branches":     len(losses),
            "n_contact_first":         timing_counts["contact_first"],
            "n_co_occurring":          timing_counts["co_occurring"],
            "n_contact_after":         timing_counts["contact_after"],
            "n_no_contact_change":     timing_counts["no_contact_change"],
            "dominant_timing":         dominant,
            "directional_class":        direction,
            "contact_state_in_cdav_spp": c_state,
            "physicochemical_type":    pc_type,
            "dar_origin_node":         first_event_node or "",
            "dar_origin_age_mya":      f"{age_mya:.1f}" if age_mya is not None else "",
        })
        processed += 1

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fields)
        w.writeheader()
        w.writerows(rows_out)

    total_skipped = sum(skip_counts.values())
    print(f"Annotated {processed} pairs  ({total_skipped} skipped)")
    if skip_counts:
        for reason, n in skip_counts.most_common():
            print(f"  skip/{reason}: {n}")
    print(f"Output → {OUT_CSV}\n")

    if rows_out:
        total = len(rows_out)
        ancestral_n = sum(1 for r in rows_out if r["is_ancestral_cdav"])
        derived_n   = total - ancestral_n
        lc_n        = sum(1 for r in rows_out if r["confidence_tier"] == LC_CONFIDENCE_LABEL)
        print(f"cDAV type:")
        print(f"  ancestral (alt_aa is mammalian ancestral state): {ancestral_n} ({100*ancestral_n/total:.1f}%)")
        print(f"  derived   (alt_aa arose in one or more lineages): {derived_n} ({100*derived_n/total:.1f}%)")
        print(f"\nAncestral classification method:")
        for meth in ("root_state", "leaf_majority", "no_data"):
            n = sum(1 for r in rows_out if r["ancestral_method"] == meth)
            print(f"  {meth:<15}: {n:>5}  ({100*n/total:.1f}%)")
        print(f"\nConfidence tier:")
        print(f"  high_confidence: {total - lc_n:>5}  ({100*(total-lc_n)/total:.1f}%)")
        print(f"  low_confidence : {lc_n:>5}  ({100*lc_n/total:.1f}%)")
        print(f"\nDominant timing (across all annotated pairs):")
        for cat in ("contact_first", "co_occurring", "contact_after", "no_contact_change"):
            n = sum(1 for r in rows_out if r["dominant_timing"] == cat)
            print(f"  {cat:<22}: {n:>5}  ({100*n/total:.1f}%)")
        print(f"\nDirectional interpretation:")
        for cat in ("permissive_background", "responding_secondary", "co_occurring_unresolved",
                    "mixed_timing", "no_detected_partner_change", "ancestral_cdav_not_directional",
                    "ambiguous_no_origin"):
            n = sum(1 for r in rows_out if r["directional_class"] == cat)
            print(f"  {cat:<30}: {n:>5}  ({100*n/total:.1f}%)")
        print(f"\nContact state in cDAV species:")
        for cat in ("conserved_human", "conserved_alt", "polymorphic", "unknown"):
            n = sum(1 for r in rows_out if r["contact_state_in_cdav_spp"] == cat)
            print(f"  {cat:<20}: {n:>5}  ({100*n/total:.1f}%)")
        print(f"\nPhysicochemical complementarity type:")
        for cat in ("charge_reversal", "charge_rescue", "volume_swap",
                    "polarity_swap", "same_direction", "unclassified"):
            n = sum(1 for r in rows_out if r["physicochemical_type"] == cat)
            print(f"  {cat:<20}: {n:>5}  ({100*n/total:.1f}%)")


if __name__ == "__main__":
    main()
