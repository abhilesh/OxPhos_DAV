"""
src/phylo/01_parse_ancestral_states.py

Parses IQTree --ancestral output (.state files + .treefile) into a compact
JSON map used by the compensatory partners analysis and timing annotation.

IQTree .state file format (per-gene):
    Node    Site    State   p_A  p_R  ...
    Node1   1       A       0.95 0.01 ...
    ...
    Sp_name 1       A       ...       (leaf nodes = species names)

IQTree .treefile: Newick with internal node labels (Node1, Node2, ...).

Output: data/phylo/ancestral_state_maps.json
{
  "gene_name": {
    "branches": {
      "parent_node|child_node": {
        "aa_pos": [parent_state, child_state],
        ...   (only positions where state changed)
      },
      ...
    },
    "leaf_to_node": {"Mus_musculus": "Mus_musculus", ...},
    "node_to_children": {"Node1": ["Node2", "Mus_musculus"], ...}
  }
}

All active state keys are harmonized to human ungapped protein positions. IQTree
reports internal-node states by alignment-site ID, so those sites are first
mapped through the Homo sapiens sequence in the exact FASTA used for IQTree.

Run from project root inside the Docker container:
    python src/phylo/01_parse_ancestral_states.py
"""

import json
import sys
from pathlib import Path

from Bio import Phylo, SeqIO

PROB_THRESHOLD    = 0.80   # minimum MAP state posterior — high-confidence calls
PROB_THRESHOLD_LC = 0.50   # low-confidence tier (0.50–0.80); flagged in output

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parents[2]
ANCS_DIR    = ROOT / "data" / "phylo" / "ancestral_states"
JOBS_DIR    = ROOT / "data" / "phylo" / "iqtree_jobs"   # for reading leaf FASTA
OUT_FILE    = ROOT / "data" / "phylo" / "ancestral_state_maps.json"
AUDIT_FILE  = ROOT / "data" / "phylo" / "asr_coordinate_harmonization_audit.tsv"


# ── Tree helpers ───────────────────────────────────────────────────────────────

def parse_tree_topology(treefile: Path) -> tuple[dict[str, list[str]], dict[str, str]]:
    """
    Parse IQTree .treefile (Newick with labelled internal nodes).

    Returns:
      node_to_children  {node_label: [child_label, ...]}
      child_to_parent   {child_label: parent_label}
    """
    tree = Phylo.read(str(treefile), "newick")

    node_to_children: dict[str, list[str]] = {}
    child_to_parent:  dict[str, str]       = {}

    for clade in tree.find_clades(order="level"):
        parent_name = clade.name or clade.confidence
        if parent_name is None:
            parent_name = "ROOT"
        parent_name = str(parent_name)

        children = []
        for child in clade.clades:
            child_name = child.name or child.confidence
            if child_name is None:
                child_name = f"anon_{id(child)}"
            child_name = str(child_name)
            children.append(child_name)
            child_to_parent[child_name] = parent_name

        if children:
            node_to_children[parent_name] = children

    return node_to_children, child_to_parent


# ── Coordinate harmonization and leaf states ──────────────────────────────────

def build_human_coordinate_maps(fasta_path: Path) -> dict:
    """
    Build coordinate maps from the Homo sapiens gapped alignment sequence.

    IQTree .state files use 1-based alignment-site IDs. Downstream variant and
    contact tables use 1-based human protein positions. Alignment columns where
    Homo sapiens is a gap/ambiguous residue are excluded from the active ASR map.
    """
    human_seq = None
    n_alignment_sites = 0
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if rec.id.split("|")[0] == "Homo_sapiens":
            human_seq = str(rec.seq).upper()
            n_alignment_sites = len(human_seq)
            break
    if human_seq is None:
        raise ValueError(f"Homo_sapiens sequence not found in {fasta_path}")

    protein_pos_to_alignment_site: dict[str, str] = {}
    alignment_site_to_protein_pos: dict[str, str] = {}
    protein_pos = 0
    dropped = 0
    for alignment_site, aa in enumerate(human_seq, 1):
        if aa in ("-", "X", "*"):
            dropped += 1
            continue
        protein_pos += 1
        protein_pos_to_alignment_site[str(protein_pos)] = str(alignment_site)
        alignment_site_to_protein_pos[str(alignment_site)] = str(protein_pos)

    return {
        "coordinate_system": "human_protein_position",
        "protein_pos_to_alignment_site": protein_pos_to_alignment_site,
        "alignment_site_to_protein_pos": alignment_site_to_protein_pos,
        "n_alignment_sites": n_alignment_sites,
        "n_human_protein_positions": protein_pos,
        "n_human_gap_alignment_sites_dropped": dropped,
    }


def read_leaf_states(
    fasta_path: Path,
    alignment_site_to_protein_pos: dict[str, str],
) -> dict[str, dict[str, str]]:
    """
    Read observed (leaf) AA states from the FASTA alignment used for IQTree.
    IQTree .state files only contain internal nodes; leaf states come from the input FASTA.

    Returns {species: {human_protein_pos: state}}. Every species is sampled at
    the same alignment columns defined by the Homo sapiens coordinate map; species
    gaps/ambiguous residues at those columns are omitted.
    """
    leaf_states: dict[str, dict[str, str]] = {}
    if not fasta_path.exists():
        return leaf_states
    for rec in SeqIO.parse(fasta_path, "fasta"):
        sp  = rec.id.split("|")[0]
        seq = str(rec.seq).upper()
        states: dict[str, str] = {}
        for alignment_site, ch in enumerate(seq, 1):
            protein_pos = alignment_site_to_protein_pos.get(str(alignment_site))
            if protein_pos is None:
                continue
            if ch in ("-", "X", "*"):
                continue
            states[protein_pos] = ch
        if states:
            leaf_states[sp] = states
    return leaf_states


def remap_node_states_to_human_positions(
    node_states: dict[str, dict[int, str]],
    alignment_site_to_protein_pos: dict[str, str],
) -> dict[str, dict[int, str]]:
    """Convert parsed IQTree alignment-site keys to human protein-position keys."""
    remapped: dict[str, dict[int, str]] = {}
    for node, states in node_states.items():
        out: dict[int, str] = {}
        for alignment_site, aa in states.items():
            protein_pos = alignment_site_to_protein_pos.get(str(alignment_site))
            if protein_pos is None:
                continue
            out[int(protein_pos)] = aa
        if out:
            remapped[node] = out
    return remapped


# ── State file parser ──────────────────────────────────────────────────────────

def parse_state_file(
    state_file: Path,
) -> tuple[dict[str, dict[int, str]], dict[str, dict[int, str]]]:
    """
    Parse IQTree .state file.

    Returns:
      node_states_hc  -- high-confidence: MAP prob >= PROB_THRESHOLD (0.80)
                         {node_label: {alignment_site (1-based): MAP_state}}
      node_states_lc  -- low-confidence:  PROB_THRESHOLD_LC <= MAP prob < PROB_THRESHOLD
                         Same structure; used for a separate "low_confidence" tier.

    Positions below PROB_THRESHOLD_LC (0.50) are discarded entirely.
    Root node (Node1 by IQTree convention) states are included in both dicts
    and are used directly for ancestral state classification.
    """
    node_states_hc: dict[str, dict[int, str]] = {}
    node_states_lc: dict[str, dict[int, str]] = {}
    prob_col_indices: list[int] = []

    with open(state_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            # Header: "Node  Site  State  p_A  p_R  p_N  ..."
            if parts[0] == "Node":
                prob_col_indices = list(range(3, len(parts)))
                continue

            if len(parts) < 3:
                continue

            node = parts[0]
            try:
                site = int(parts[1])
            except ValueError:
                continue

            state = parts[2]
            if state == "-" or len(state) != 1:
                continue   # gap or ambiguous

            max_prob = 1.0   # default if no prob columns
            if prob_col_indices:
                try:
                    probs = [float(parts[i]) for i in prob_col_indices if i < len(parts)]
                    max_prob = max(probs) if probs else 1.0
                except (ValueError, IndexError):
                    pass

            if max_prob >= PROB_THRESHOLD:
                node_states_hc.setdefault(node, {})[site] = state
            elif max_prob >= PROB_THRESHOLD_LC:
                node_states_lc.setdefault(node, {})[site] = state
            # else: below 0.50 — discard

    return node_states_hc, node_states_lc


# ── Build branch-change map ────────────────────────────────────────────────────

def build_branch_changes(
    node_states:      dict[str, dict[int, str]],
    child_to_parent:  dict[str, str],
) -> dict[str, dict[int, tuple[str, str]]]:
    """
    For each (parent, child) branch: record positions where state changed.

    Returns:
      {"parent|child": {aa_pos: (parent_state, child_state), ...}}

    Only positions that differ are stored (sparse representation).
    """
    branches: dict[str, dict[int, tuple[str, str]]] = {}

    for child, parent in child_to_parent.items():
        if child not in node_states or parent not in node_states:
            continue
        parent_map = node_states[parent]
        child_map  = node_states[child]
        changes: dict[int, tuple[str, str]] = {}
        for pos, child_state in child_map.items():
            parent_state = parent_map.get(pos)
            if parent_state and parent_state != child_state:
                changes[pos] = (parent_state, child_state)
        if changes:
            branch_id = f"{parent}|{child}"
            branches[branch_id] = {str(k): list(v) for k, v in changes.items()}

    return branches


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    if not ANCS_DIR.exists():
        print(f"ERROR: {ANCS_DIR} not found — copy IQTree outputs from HPC first.")
        sys.exit(1)

    results: dict = {}
    audit_rows: list[dict[str, object]] = []
    n_genes = 0

    for gene_dir in sorted(ANCS_DIR.iterdir()):
        if not gene_dir.is_dir():
            continue
        gene = gene_dir.name

        state_files = list(gene_dir.glob("*.state"))
        tree_files  = list(gene_dir.glob("*.treefile"))

        if not state_files or not tree_files:
            print(f"  [SKIP]  {gene}: missing .state or .treefile")
            continue

        state_file = state_files[0]
        tree_file  = tree_files[0]

        try:
            node_to_children, child_to_parent = parse_tree_topology(tree_file)
            node_states_hc_aln, node_states_lc_aln = parse_state_file(state_file)
            fasta_path  = JOBS_DIR / gene / f"{gene}.fasta"
            coord_maps = build_human_coordinate_maps(fasta_path)
            node_states_hc = remap_node_states_to_human_positions(
                node_states_hc_aln,
                coord_maps["alignment_site_to_protein_pos"],
            )
            node_states_lc = remap_node_states_to_human_positions(
                node_states_lc_aln,
                coord_maps["alignment_site_to_protein_pos"],
            )

            # Leaf states from FASTA — IQTree .state only covers internal nodes;
            # observed species states are read here and merged into the HC node
            # state dict so that terminal branches (parent_internal→leaf) are
            # included in build_branch_changes.  Observed data has certainty 1.0,
            # so leaf states are always treated as high-confidence.
            leaf_states = read_leaf_states(
                fasta_path,
                coord_maps["alignment_site_to_protein_pos"],
            )
            for sp, states_str in leaf_states.items():
                # leaf_states uses str keys; node_states_hc uses int keys
                node_states_hc[sp] = {int(pos): aa for pos, aa in states_str.items()}

            branches    = build_branch_changes(node_states_hc, child_to_parent)
            branches_lc = build_branch_changes(node_states_lc, child_to_parent)

            # Leaf nodes: tips in tree topology with no children
            leaf_nodes = set(child_to_parent.keys()) - set(node_to_children.keys())

            # Root node: parent that is never a child
            all_children = {c for cs in node_to_children.values() for c in cs}
            root_nodes   = set(node_to_children.keys()) - all_children
            root_node    = next(iter(root_nodes)) if root_nodes else None

            # Root node states (high-confidence only) — used for ancestral cDAV detection
            root_states: dict[str, str] = {}
            if root_node and root_node in node_states_hc:
                root_states = {str(k): v for k, v in node_states_hc[root_node].items()}

            results[gene] = {
                "coordinate_system":   coord_maps["coordinate_system"],
                "protein_pos_to_alignment_site": coord_maps["protein_pos_to_alignment_site"],
                "alignment_site_to_protein_pos": coord_maps["alignment_site_to_protein_pos"],
                "n_alignment_sites":   coord_maps["n_alignment_sites"],
                "n_human_protein_positions": coord_maps["n_human_protein_positions"],
                "n_human_gap_alignment_sites_dropped": coord_maps["n_human_gap_alignment_sites_dropped"],
                "branches":            branches,
                "branches_lc":         branches_lc,   # low-confidence tier
                "node_to_children":    node_to_children,
                "leaf_nodes":          sorted(leaf_nodes),
                "leaf_states":         leaf_states,
                "root_node":           root_node,
                "root_states":         root_states,   # {str(human protein pos): state} at PP>=0.80
                "n_branches_with_changes": len(branches),
                "n_sites":             coord_maps["n_human_protein_positions"],
            }
            audit_rows.append({
                "gene": gene,
                "coordinate_system": coord_maps["coordinate_system"],
                "n_alignment_sites": coord_maps["n_alignment_sites"],
                "n_human_protein_positions": coord_maps["n_human_protein_positions"],
                "n_human_gap_alignment_sites_dropped": coord_maps["n_human_gap_alignment_sites_dropped"],
                "n_hc_internal_nodes": len(node_states_hc_aln),
                "n_lc_internal_nodes": len(node_states_lc_aln),
                "n_leaf_species": len(leaf_states),
                "n_branches_with_changes": len(branches),
                "n_branches_lc_with_changes": len(branches_lc),
            })
            n_genes += 1
            print(
                f"  [OK]    {gene}: {len(node_states_hc)} nodes (HC), "
                f"{len(branches)} branches with AA changes"
                + (f", root={root_node}" if root_node else "")
            )

        except Exception as e:
            print(f"  [ERROR] {gene}: {e}")

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_FILE, "w") as f:
        json.dump(results, f)   # no indent — file can be large

    if audit_rows:
        import csv
        with open(AUDIT_FILE, "w", newline="") as f:
            fieldnames = list(audit_rows[0])
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(audit_rows)

    print(f"\nParsed {n_genes} genes → {OUT_FILE}")
    print(f"Coordinate harmonization audit → {AUDIT_FILE}")


if __name__ == "__main__":
    main()
