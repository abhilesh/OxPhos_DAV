#!/usr/bin/env python3
"""
src/structural/01_find_compensating_partners.py

For every non-discarded cDAV (Tier A, B, C), test whether structurally
contacting residues show co-evolutionary enrichment in species carrying
the disease amino acid.

Statistical tests
─────────────────
Three tests are run per (DAR × contact position × alt AA) tuple, in order
of phylogenetic validity:

  1. Fisher's exact (species as units) — retained as comparison column.
     INVALID as a primary test because cDAV-carrying species are phylogenetically
     clustered; included for backward compatibility only.

  2. Pagel's discrete (R phytools::fitPagel) — LRT for correlated binary
     character evolution on the pruned VertLife tree. Valid phylogenetic test.
     Requires data/phylo/species_tree/mammaltree.nwk and R + phytools.
     Skipped gracefully if tree is absent.

  3. Branch co-occurrence (IQTree ancestral states) — counts branches where
     both the DAR alt AA and the contact alt AA arose. Fisher's exact over
     branches is valid because branches are approximately independent events.
     Requires data/phylo/ancestral_state_maps.json.
     Skipped gracefully if ancestral maps are absent.

Outputs
───────
  results/structural/all_tested_pairs.csv   — ALL tested pairs, all test columns
  results/structural/compensatory_partners.csv — derived significant view
  results/structural/concordance_summary.csv   — per-cDAV enrichment counts

Significance threshold for compensatory_partners.csv:
  (pagel_fdr ≤ 0.10 OR branch_cooccur_fdr ≤ 0.10) AND low_power == False
  Falls back to fisher_fdr ≤ 0.10 if neither phylogenetic test has run yet.

Run from project root:
  python src/structural/01_find_compensating_partners.py
"""

import copy
import csv
import json
import os
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from scipy.stats import fisher_exact

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parents[2]
CLASSIFIED_DIR = ROOT / "data" / "derived" / "classified"
CLASSIFIED_PARQUET = CLASSIFIED_DIR / "variants_master_classified.parquet"
CONTACTS    = ROOT / "results" / "structural" / "dar_contacts_cbcb8A.csv"
TOGA_AA_DIR = ROOT / "data" / "alignments" / "toga_hg38_aa"
MT_AA_DIR   = ROOT / "data" / "alignments" / "mtdna_aa"
OUT_DIR     = ROOT / "results" / "structural"

# Phylogenetic test inputs (optional — graceful degradation if absent)
VERT_TREE   = ROOT / "data" / "phylo" / "species_tree" / "mammaltree_crossgenome.nwk"
ANC_MAPS    = ROOT / "data" / "phylo" / "ancestral_state_maps.json"
CROSS_SPP   = ROOT / "data" / "phylo" / "cross_genome_species.txt"
PAGEL_R     = ROOT / "src" / "phylo" / "pagel_discrete.R"

_MASK            = {"-", "X", "!", "*"}
_INCOMP_CLASSES  = {"hbond", "electrostatic"}
_INCOMP_MIN_SENS = 0.5
_INCOMP_MIN_SPEC = 0.5
_PAGEL_MIN_SPP   = 20    # minimum species for Pagel test

_ALIAS_EQUIV = {"COXFA4": "NDUFA4", "NDUFA4": "COXFA4"}
_TIER_ORDER  = {"Tier A": 0, "Tier B": 1, "Tier C": 2}


# ── Alignment loading ─────────────────────────────────────────────────────────

def load_alignments(genes: set) -> dict[str, dict[str, str]]:
    """Load per-gene AA alignments. Returns {gene: {species_shortname: gapped_seq}}.

    Alignments are stored under both the FASTA filename key and its alias (if any)
    so that structural-layer names (e.g. NDUFA4 from contacts CSV) and parquet-layer
    names (e.g. COXFA4 from interpreted_gene) both resolve to the same sequences.
    """
    alns: dict[str, dict[str, str]] = {}
    for aln_dir in (TOGA_AA_DIR, MT_AA_DIR):
        if not aln_dir.exists():
            continue
        for fasta in aln_dir.glob("*_aa_alignment.fasta"):
            gene = fasta.name.replace("_aa_alignment.fasta", "")
            if gene not in genes and _ALIAS_EQUIV.get(gene) not in genes:
                continue
            seqs = {}
            for rec in SeqIO.parse(fasta, "fasta"):
                sp = rec.id.split("|")[0]
                seqs[sp] = str(rec.seq).upper()
            alns[gene] = seqs
            if gene in _ALIAS_EQUIV:
                alns[_ALIAS_EQUIV[gene]] = seqs
    return alns


def build_pos_to_col(ref_seq: str) -> dict[int, int]:
    """1-based ungapped biological position → 0-based alignment column."""
    pos_map: dict[int, int] = {}
    pos = 0
    for col, ch in enumerate(ref_seq):
        if ch not in _MASK:
            pos += 1
            pos_map[pos] = col
    return pos_map


def get_ref_seq(aln: dict[str, str]) -> tuple[str, str]:
    for sp, seq in aln.items():
        if sp == "Homo_sapiens":
            return sp, seq
    raise ValueError("Human reference not found in alignment.")


# ── BH FDR ────────────────────────────────────────────────────────────────────

def bh_fdr(p_values: list[float]) -> list[float]:
    n = len(p_values)
    if n == 0:
        return []
    order = np.argsort(p_values)
    adj   = np.array(p_values, dtype=float)
    for rank, idx in enumerate(order, 1):
        adj[idx] = min(1.0, p_values[idx] * n / rank)
    for i in range(n - 2, -1, -1):
        adj[order[i]] = min(adj[order[i]], adj[order[i + 1]])
    return adj.tolist()


# ── Cross-genome species overlap ──────────────────────────────────────────────

def load_cross_genome_species() -> set[str]:
    if not CROSS_SPP.exists():
        return set()
    return set(CROSS_SPP.read_text().splitlines())


# ── VertLife tree helpers ─────────────────────────────────────────────────────

_master_tree = None   # loaded once

def get_master_tree():
    global _master_tree
    if _master_tree is None and VERT_TREE.exists():
        _master_tree = Phylo.read(str(VERT_TREE), "newick")
    return _master_tree


def prune_tree_to_species(species: set[str]):
    """Return a pruned copy of the VertLife tree for the given species set."""
    tree = get_master_tree()
    if tree is None:
        return None
    all_tips = {c.name for c in tree.get_terminals()}
    keep = species & all_tips
    if len(keep) < _PAGEL_MIN_SPP:
        return None, keep
    # Work on a fresh read each time (Bio.Phylo prune is in-place)
    t = Phylo.read(str(VERT_TREE), "newick")
    remove = {c.name for c in t.get_terminals()} - keep
    for name in remove:
        t.prune(name)
    remaining = {c.name for c in t.get_terminals()}
    return t, remaining


# ── Pagel's discrete test ─────────────────────────────────────────────────────

def pagel_discrete_test(
    dar_spp: set[str],
    all_readable_spp: set[str],
    contact_alt_spp: set[str],
    pair_type: str,
    cross_genome_spp: set[str],
    tmp_dir: str,
) -> dict:
    """
    Run Pagel's discrete LRT via R subprocess.
    Returns dict with keys: pagel_p, n_species_in_test, low_power
    """
    null_result = {"pagel_p": None, "n_species_in_test": 0, "low_power": True}

    if not VERT_TREE.exists() or not PAGEL_R.exists():
        return null_result

    # Resolve species set for this pair
    test_spp = all_readable_spp
    if pair_type == "mt-nuc" and cross_genome_spp:
        test_spp = test_spp & cross_genome_spp

    n = len(test_spp)
    if n < _PAGEL_MIN_SPP:
        return {"pagel_p": None, "n_species_in_test": n, "low_power": True}

    pruned_tree, pruned_spp = prune_tree_to_species(test_spp)
    if pruned_tree is None or len(pruned_spp) < _PAGEL_MIN_SPP:
        return {"pagel_p": None, "n_species_in_test": len(pruned_spp or []), "low_power": True}

    # Write temp files
    tree_file   = os.path.join(tmp_dir, "tree.nwk")
    trait1_file = os.path.join(tmp_dir, "trait1.txt")
    trait2_file = os.path.join(tmp_dir, "trait2.txt")

    Phylo.write(pruned_tree, tree_file, "newick")

    with open(trait1_file, "w") as f:
        for sp in pruned_spp:
            f.write(f"{sp}\t{1 if sp in dar_spp else 0}\n")
    with open(trait2_file, "w") as f:
        for sp in pruned_spp:
            f.write(f"{sp}\t{1 if sp in contact_alt_spp else 0}\n")

    try:
        result = subprocess.run(
            ["Rscript", str(PAGEL_R), tree_file, trait1_file, trait2_file],
            capture_output=True, text=True, timeout=120,
        )
        for line in result.stdout.splitlines():
            if line.startswith("pagel_p\t"):
                val = line.split("\t", 1)[1].strip()
                p = None if val == "NA" else float(val)
                return {
                    "pagel_p": p,
                    "n_species_in_test": len(pruned_spp),
                    "low_power": len(pruned_spp) < _PAGEL_MIN_SPP,
                }
    except Exception:
        pass

    return {"pagel_p": None, "n_species_in_test": len(pruned_spp), "low_power": False}


# ── Branch co-occurrence test ─────────────────────────────────────────────────

def _subtended_leaves(node: str, node_to_children: dict, leaf_nodes: set) -> frozenset:
    queue = [node]
    leaves = set()
    while queue:
        n = queue.pop()
        if n in leaf_nodes:
            leaves.add(n)
        for child in node_to_children.get(n, []):
            queue.append(child)
    return frozenset(leaves)


def branch_cooccurrence_test(
    anc_maps: dict,
    dar_gene: str,
    dar_pos: int,
    dar_alt_aa: str,
    contact_gene: str,
    contact_pos: int,
    contact_alt_aa: str,
) -> dict:
    """
    Count branches where DAR alt AA and/or contact alt AA arose.
    Returns Fisher p-value over branches + raw counts.
    """
    null = {
        "branch_cooccur_p": None,
        "n_cooccur_branches": 0,
        "n_dar_only_branches": 0,
        "n_contact_only_branches": 0,
    }

    if dar_gene not in anc_maps or contact_gene not in anc_maps:
        return null

    dar_map     = anc_maps[dar_gene]
    contact_map = anc_maps[contact_gene]

    # Build leaf subtree fingerprint per branch (for alignment across genes)
    def origin_branches(gene_map: dict, pos: int, alt_aa: str) -> set[frozenset]:
        """Return set of leaf-subtree frozensets for branches where alt_aa arose."""
        origins = set()
        n2c  = gene_map.get("node_to_children", {})
        leaves = set(gene_map.get("leaf_nodes", []))
        for branch_id, changes in gene_map.get("branches", {}).items():
            entry = changes.get(str(pos))
            if entry and entry[1] == alt_aa and entry[0] != alt_aa:
                _, child = branch_id.split("|", 1)
                subtree = _subtended_leaves(child, n2c, leaves)
                if subtree:
                    origins.add(subtree)
        return origins

    dar_origins     = origin_branches(dar_map,     dar_pos,     dar_alt_aa)
    contact_origins = origin_branches(contact_map, contact_pos, contact_alt_aa)

    if not dar_origins:
        return null

    # Match branches by subtree overlap (same clade = same branch)
    cooccur = contact_only = 0
    dar_only = 0

    for d_sub in dar_origins:
        matched = any(
            len(d_sub & c_sub) / max(len(d_sub), len(c_sub)) > 0.8
            for c_sub in contact_origins
        )
        if matched:
            cooccur += 1
        else:
            dar_only += 1

    for c_sub in contact_origins:
        matched = any(
            len(d_sub & c_sub) / max(len(d_sub), len(c_sub)) > 0.8
            for d_sub in dar_origins
        )
        if not matched:
            contact_only += 1

    # Fisher: cooccur vs dar_only in DAR branches
    # (is co-change more than expected given total branch counts?)
    total_dar_branches = len(dar_origins)
    a = cooccur
    c = dar_only
    total_contact_branches = len(contact_origins)
    b = contact_only
    d = max(0, total_dar_branches - a)   # DAR branches with no co-change

    if a == 0:
        p_val = 1.0
    else:
        try:
            _, p_val = fisher_exact([[a, b], [c, d]], alternative="greater")
        except Exception:
            p_val = 1.0

    return {
        "branch_cooccur_p": p_val,
        "n_cooccur_branches": cooccur,
        "n_dar_only_branches": dar_only,
        "n_contact_only_branches": contact_only,
    }


def _species_set(raw) -> set[str]:
    """
    Coerce a `lineages_with_disease_allele`-style field to a species set.

    The Parquet-loaded value is normally a native list, but if it were ever
    a JSON-encoded string instead, `set(raw)` would silently iterate its
    characters rather than raise — corrupting every downstream Fisher/Pagel/
    branch-co-occurrence test. Guard against both representations explicitly.
    """
    if isinstance(raw, str):
        raw = json.loads(raw) if raw else []
    return set(raw or [])


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # ── Load cDAVs ─────────────────────────────────────────────────────────────
    cdavs: dict[str, dict] = {}
    if not CLASSIFIED_PARQUET.exists():
        raise FileNotFoundError(f"Missing classified master table: {CLASSIFIED_PARQUET}")
    for var in pd.read_parquet(CLASSIFIED_PARQUET).to_dict(orient="records"):
        if var.get("classification_status") != "classified":
            continue
        if not var.get("is_cdav_amino_acid"):
            continue
        cdavs[var["variant_id"]] = var

    by_tier = Counter(v["tier"] for v in cdavs.values())
    print(f"Non-discarded cDAVs: {len(cdavs)}")
    for t in sorted(by_tier, key=lambda x: _TIER_ORDER.get(x, 99)):
        print(f"  {t}: {by_tier[t]}")

    if not cdavs:
        print("No cDAVs found.")
        return

    # ── Load contacts ─────────────────────────────────────────────────────────
    contacts: list[dict] = []
    with open(CONTACTS) as f:
        for row in csv.DictReader(f):
            contact_key = row.get("variant_id") or row.get("ann_id")
            if contact_key in cdavs:
                row["variant_id"] = contact_key
                contacts.append(row)
    print(f"Contacts for these DARs: {len(contacts)}")

    needed_genes = {
        g for row in contacts
        for g in (row["dar_locus"], row["contact_gene"])
        if not g.startswith("Unknown(")
    }
    alns = load_alignments(needed_genes)
    print(f"Alignments loaded for {len(alns)} genes.")

    # ── Optional: load ancestral maps and VertLife tree ───────────────────────
    anc_maps: dict = {}
    if ANC_MAPS.exists():
        with open(ANC_MAPS) as f:
            anc_maps = json.load(f)
        print(f"Ancestral state maps loaded for {len(anc_maps)} genes.")
    else:
        print("Ancestral state maps not found — branch co-occurrence test will be skipped.")

    cross_genome_spp = load_cross_genome_species()
    if cross_genome_spp:
        print(f"Cross-genome species overlap: {len(cross_genome_spp)}")

    vert_available = VERT_TREE.exists() and PAGEL_R.exists()
    if vert_available:
        print(f"VertLife tree found — Pagel's discrete test will run (batch mode).")
    else:
        print("VertLife tree or pagel_discrete.R not found — Pagel test will be skipped.")

    # ── Per-DAR tests ─────────────────────────────────────────────────────────
    by_dar: dict[str, list] = defaultdict(list)
    for row in contacts:
        by_dar[row["variant_id"]].append(row)

    all_records: list[dict] = []
    summary_rows: list[dict] = []
    # Pagel batch: list of (pair_idx_in_all_records, dar_spp, test_spp, contact_alt_spp)
    pagel_pending: list[tuple[int, set, set, set]] = []

    master_tree = get_master_tree() if vert_available else None

    # ── Pass 1: Fisher + branch co-occurrence (fast), collect Pagel inputs ────
    for variant_id, dar_contacts in by_dar.items():
        var          = cdavs[variant_id]
        dar_gene     = var.get("interpreted_gene") or var.get("classification_gene") or var.get("locus", "")
        tier         = var["tier"]
        pair_type    = dar_contacts[0].get("contact_type", "")
        cdav_spp     = _species_set(var.get("lineages_with_disease_allele", []))
        dar_aa_coord = int(dar_contacts[0]["dar_aa_coord"])

        if dar_gene not in alns:
            continue

        dar_aln = alns[dar_gene]
        try:
            _, dar_ref_gapped = get_ref_seq(dar_aln)
        except ValueError:
            continue
        dar_col_map = build_pos_to_col(dar_ref_gapped)
        if dar_aa_coord not in dar_col_map:
            continue
        dar_col = dar_col_map[dar_aa_coord]

        readable_spp = {
            sp for sp, seq in dar_aln.items()
            if sp != "Homo_sapiens"
            and len(seq) > dar_col
            and seq[dar_col] not in _MASK
        }
        cdav_in_aln = cdav_spp & readable_spp
        bg_spp      = readable_spp - cdav_spp

        if not cdav_in_aln:
            continue

        seen_contacts: set = set()
        raw_tests: list[dict] = []

        for row in dar_contacts:
            contact_gene = row["contact_gene"]
            dar_structure_gene = row.get("dar_structure_gene", row.get("dar_locus", ""))
            contact_structure_gene = row.get("contact_structure_gene", contact_gene)
            if contact_gene.startswith("Unknown(") or contact_gene not in alns:
                continue
            try:
                contact_pos = int(row["contact_refseq_pos"])
            except (ValueError, TypeError):
                continue

            key = (contact_gene, contact_pos)
            if key in seen_contacts:
                continue
            seen_contacts.add(key)

            c_aln = alns[contact_gene]
            try:
                _, c_ref_gapped = get_ref_seq(c_aln)
            except ValueError:
                continue
            c_col_map = build_pos_to_col(c_ref_gapped)
            if contact_pos not in c_col_map:
                continue
            c_col = c_col_map[contact_pos]

            human_contact_aa = row["contact_aa"]
            contact_class    = row["contact_class"]
            contact_type     = row["contact_type"]

            cdav_contact_aas = [
                c_aln[sp][c_col]
                for sp in cdav_in_aln
                if sp in c_aln and len(c_aln[sp]) > c_col
                and c_aln[sp][c_col] not in _MASK
            ]
            bg_contact_aas = [
                c_aln[sp][c_col]
                for sp in bg_spp
                if sp in c_aln and len(c_aln[sp]) > c_col
                and c_aln[sp][c_col] not in _MASK
            ]

            if not cdav_contact_aas:
                continue

            for alt_aa in {aa for aa in cdav_contact_aas if aa != human_contact_aa}:
                a = cdav_contact_aas.count(alt_aa)
                c = len(cdav_contact_aas) - a
                b = bg_contact_aas.count(alt_aa)
                d = len(bg_contact_aas) - b

                if a == 0:
                    continue

                _, fisher_p = fisher_exact([[a, b], [c, d]], alternative="greater")
                sensitivity = a / (a + c) if (a + c) > 0 else 0.0
                specificity = a / (a + b) if (a + b) > 0 else 0.0

                incompatible = (
                    contact_class in _INCOMP_CLASSES
                    and sensitivity >= _INCOMP_MIN_SENS
                    and specificity >= _INCOMP_MIN_SPEC
                )

                # Branch co-occurrence test (fast — no subprocess)
                branch_result = branch_cooccurrence_test(
                    anc_maps, dar_gene, dar_aa_coord, var.get("alt_aa", ""),
                    contact_gene, contact_pos, alt_aa,
                ) if anc_maps else {
                    "branch_cooccur_p": None,
                    "n_cooccur_branches": 0,
                    "n_dar_only_branches": 0,
                    "n_contact_only_branches": 0,
                }

                # Pagel species set (computed now; R call deferred to batch)
                contact_alt_spp = {
                    sp for sp in readable_spp
                    if sp in c_aln and len(c_aln[sp]) > c_col
                    and c_aln[sp][c_col] == alt_aa
                }
                test_spp = readable_spp
                if contact_type == "mt-nuc" and cross_genome_spp:
                    test_spp = test_spp & cross_genome_spp
                n_test = len(test_spp)
                low_power = n_test < _PAGEL_MIN_SPP

                rec_idx = len(all_records) + len(raw_tests)
                # Only run Pagel for pairs with some Fisher signal (p < 0.20).
                # Pairs with no signal at the species level will also lack
                # phylogenetic signal; running Pagel on all 8087 pairs takes
                # ~45 min due to fitPagel convergence overhead.
                if vert_available and not low_power and fisher_p < 0.20:
                    pagel_pending.append((rec_idx, cdav_in_aln, test_spp, contact_alt_spp))

                raw_tests.append({
                    "variant_id":              variant_id,
                    "ann_id":                  var.get("ann_id", ""),
                    "tier":                    tier,
                    "dar_gene":                dar_gene,
                    "dar_structure_gene":      dar_structure_gene,
                    "dar_genome":              var.get("genome", ""),
                    "dar_aa_coord":            dar_aa_coord,
                    "dar_ref_aa":              var.get("ref_aa", ""),
                    "dar_alt_aa":              var.get("alt_aa", ""),
                    "contact_gene":            contact_gene,
                    "contact_structure_gene":  contact_structure_gene,
                    "contact_refseq_pos":      contact_pos,
                    "contact_human_aa":        human_contact_aa,
                    "contact_alt_aa":          alt_aa,
                    "contact_class":           contact_class,
                    "contact_type":            contact_type,
                    "likely_incompatible":     incompatible,
                    "n_cdav_spp":              len(cdav_in_aln),
                    "n_bg_spp":                len(bg_spp),
                    "n_cdav_with_alt":         a,
                    "n_bg_with_alt":           b,
                    "sensitivity":             sensitivity,
                    "specificity":             specificity,
                    "fisher_p":                fisher_p,
                    "fisher_fdr":              None,
                    "n_species_in_test":       n_test,
                    "low_power":               low_power,
                    "pagel_p":                 None,
                    "pagel_fdr":               None,
                    "branch_cooccur_p":        branch_result["branch_cooccur_p"],
                    "branch_cooccur_fdr":      None,
                    "n_cooccur_branches":      branch_result["n_cooccur_branches"],
                    "n_dar_only_branches":     branch_result["n_dar_only_branches"],
                    "n_contact_only_branches": branch_result["n_contact_only_branches"],
                })

        if not raw_tests:
            continue

        # fisher_fdr / branch_cooccur_fdr are left unset here and corrected
        # globally below, once every DAR has been processed — see the
        # "Global FDR correction" block. Per-DAR correction would give each
        # DAR its own independent FDR scope while _get_sig() checks the
        # result against a single fixed threshold, which is not a valid
        # global false-discovery-rate guarantee.
        all_records.extend(raw_tests)

    # ── Global FDR correction (Fisher + branch co-occurrence) ──────────────────
    # Corrected once across every tested pair, matching the Pagel FDR pass
    # below (Pass 2) rather than per-DAR. Per-DAR correction previously let
    # _get_sig()'s fixed 0.10 threshold silently mean different things
    # depending on how many contacts a given DAR happened to have.
    valid_fisher = [(i, r["fisher_p"]) for i, r in enumerate(all_records)
                    if r["fisher_p"] is not None]
    if valid_fisher:
        idxs, ps = zip(*valid_fisher)
        fdrs = bh_fdr(list(ps))
        for idx, fdr in zip(idxs, fdrs):
            all_records[idx]["fisher_fdr"] = fdr
        print(f"Fisher FDR computed globally for {len(valid_fisher)} pairs.")

    valid_branch = [(i, r["branch_cooccur_p"]) for i, r in enumerate(all_records)
                    if r["branch_cooccur_p"] is not None]
    if valid_branch:
        idxs, ps = zip(*valid_branch)
        fdrs = bh_fdr(list(ps))
        for idx, fdr in zip(idxs, fdrs):
            all_records[idx]["branch_cooccur_fdr"] = fdr
        print(f"Branch co-occurrence FDR computed globally for {len(valid_branch)} pairs.")

    # ── Pass 2: Pagel batch (one R process for all pairs) ─────────────────────
    if vert_available and pagel_pending:
        print(f"\nRunning Pagel batch test for {len(pagel_pending)} pairs (single R session)...")
        with tempfile.TemporaryDirectory() as tmp_dir:
            manifest_lines: list[str] = []
            all_tips = {c.name for c in master_tree.get_terminals()}

            # Cache pruned trees by their species frozenset (read tree once per unique species set)
            _pruned_tree_cache: dict[frozenset, tuple] = {}

            def get_pruned(test_spp_set: frozenset):
                """Returns (pruned_newick_path, pruned_spp_set) or (None, None)."""
                if test_spp_set in _pruned_tree_cache:
                    return _pruned_tree_cache[test_spp_set]
                keep = test_spp_set & all_tips
                if len(keep) < _PAGEL_MIN_SPP:
                    _pruned_tree_cache[test_spp_set] = (None, None)
                    return None, None
                t = copy.deepcopy(master_tree)
                remove = {c.name for c in t.get_terminals()} - keep
                for name in remove:
                    t.prune(name)
                pruned_spp = frozenset(c.name for c in t.get_terminals())
                if len(pruned_spp) < _PAGEL_MIN_SPP:
                    _pruned_tree_cache[test_spp_set] = (None, None)
                    return None, None
                # Write once, reuse path for all pairs with same species set
                tree_f = os.path.join(tmp_dir, f"tree_{len(_pruned_tree_cache)}.nwk")
                Phylo.write(t, tree_f, "newick")
                result = (tree_f, pruned_spp)
                _pruned_tree_cache[test_spp_set] = result
                return result

            for idx, (rec_idx, dar_spp, test_spp, contact_alt_spp) in enumerate(pagel_pending):
                tree_f, pruned_spp = get_pruned(frozenset(test_spp))
                if tree_f is None:
                    continue

                trait1_f = os.path.join(tmp_dir, f"t1_{idx}.txt")
                trait2_f = os.path.join(tmp_dir, f"t2_{idx}.txt")

                with open(trait1_f, "w") as f:
                    for sp in pruned_spp:
                        f.write(f"{sp}\t{1 if sp in dar_spp else 0}\n")
                with open(trait2_f, "w") as f:
                    for sp in pruned_spp:
                        f.write(f"{sp}\t{1 if sp in contact_alt_spp else 0}\n")

                manifest_lines.append(f"{rec_idx}\t{tree_f}\t{trait1_f}\t{trait2_f}")
            print(f"  Unique tree topologies: {len(_pruned_tree_cache)}")

            if manifest_lines:
                manifest_f = os.path.join(tmp_dir, "manifest.tsv")
                with open(manifest_f, "w") as f:
                    f.write("\n".join(manifest_lines) + "\n")

                # Stream results line-by-line so partial output is saved if
                # the wall-time budget is exhausted before all pairs complete.
                # The R script flushes stdout after each pair, so every completed
                # pair is recoverable even on an early termination.
                _PAGEL_WALL_BUDGET = 7200   # 2 h; use HPC chunk mode for larger runs
                pagel_lines: list[str] = []
                deadline = time.monotonic() + _PAGEL_WALL_BUDGET
                try:
                    proc = subprocess.Popen(
                        ["Rscript", str(PAGEL_R), manifest_f],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                    )
                    for line in proc.stdout:
                        pagel_lines.append(line)
                        if time.monotonic() > deadline:
                            proc.terminate()
                            print(f"  Pagel wall-time budget ({_PAGEL_WALL_BUDGET}s) reached "
                                  f"— {len(pagel_lines)} pairs completed, remainder skipped.")
                            break
                    proc.stdout.close()
                    proc.wait()
                    stderr_out = proc.stderr.read()
                    print(f"  Pagel batch complete ({len(pagel_lines)} pairs processed).")
                    if stderr_out:
                        print(f"  R stderr (first 500 chars): {stderr_out[:500]}")
                except Exception as e:
                    print(f"  Pagel batch failed: {e}")

                for line in pagel_lines:
                    parts = line.strip().split("\t")
                    if len(parts) == 3 and parts[1] == "pagel_p":
                        try:
                            rec_idx = int(parts[0])
                            val_str = parts[2]
                            p = None if val_str == "NA" else float(val_str)
                            all_records[rec_idx]["pagel_p"] = p
                        except (ValueError, IndexError):
                            pass

        # Pagel FDR across all records
        valid_pagel = [(i, r["pagel_p"]) for i, r in enumerate(all_records)
                       if r["pagel_p"] is not None]
        if valid_pagel:
            idxs, ps = zip(*valid_pagel)
            fdrs = bh_fdr(list(ps))
            for idx, fdr in zip(idxs, fdrs):
                all_records[idx]["pagel_fdr"] = fdr
            print(f"  Pagel FDR computed for {len(valid_pagel)} pairs.")

    # ── Per-DAR concordance summary (uses fully annotated all_records) ─────────
    by_ann_id = defaultdict(list)
    for rec in all_records:
        by_ann_id[rec["variant_id"]].append(rec)

    for variant_id, recs in by_ann_id.items():
        var      = cdavs[variant_id]
        dar_gene = var.get("interpreted_gene") or var.get("classification_gene") or var.get("locus", "")
        dar_structure_gene = recs[0].get("dar_structure_gene", "")
        tier     = var["tier"]

        dar_aln = alns.get(dar_gene)
        if dar_aln is None:
            continue
        try:
            _, dar_ref_gapped = get_ref_seq(dar_aln)
        except ValueError:
            continue
        dar_col_map  = build_pos_to_col(dar_ref_gapped)
        dar_aa_coord = int(recs[0]["dar_aa_coord"])
        dar_col      = dar_col_map.get(dar_aa_coord)

        cdav_spp    = _species_set(var.get("lineages_with_disease_allele", []))
        readable_spp = set()
        if dar_col is not None:
            readable_spp = {
                sp for sp, seq in dar_aln.items()
                if sp != "Homo_sapiens" and len(seq) > dar_col
                and seq[dar_col] not in _MASK
            }

        cdav_in_aln_summary = cdav_spp & readable_spp
        bg_spp_summary      = readable_spp - cdav_spp
        sig = _get_sig(recs)
        summary_rows.append({
            "variant_id":        variant_id,
            "ann_id":            var.get("ann_id", ""),
            "tier":              tier,
            "dar_gene":          dar_gene,
            "dar_structure_gene": dar_structure_gene,
            "dar_genome":        var.get("genome", ""),
            "dar_aa_coord":      dar_aa_coord,
            "dar_ref_aa":        var.get("ref_aa", ""),
            "dar_alt_aa":        var.get("alt_aa", ""),
            "n_cdav_spp":        len(cdav_in_aln_summary),
            "n_bg_spp":          len(bg_spp_summary),
            "n_contacts_tested": len(recs),
            "n_sig_contacts":    len(sig),
            "n_incompatible":    sum(1 for t in sig if t["likely_incompatible"]),
            "intra_genomic_sig": sum(1 for t in sig if t["contact_type"] != "mt-nuc"),
            "inter_genomic_sig": sum(1 for t in sig if t["contact_type"] == "mt-nuc"),
            "sig_contact_genes": ",".join(sorted({t["contact_gene"] for t in sig})),
        })

    # ── Write all_tested_pairs.csv ────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    all_fields = [
        "variant_id", "ann_id", "tier", "dar_gene", "dar_structure_gene", "dar_genome",
        "dar_aa_coord", "dar_ref_aa", "dar_alt_aa",
        "contact_gene", "contact_structure_gene", "contact_refseq_pos", "contact_human_aa", "contact_alt_aa",
        "contact_class", "contact_type", "likely_incompatible",
        "n_cdav_spp", "n_bg_spp", "n_cdav_with_alt", "n_bg_with_alt",
        "sensitivity", "specificity",
        "fisher_p", "fisher_fdr",
        "n_species_in_test", "low_power",
        "pagel_p", "pagel_fdr",
        "branch_cooccur_p", "branch_cooccur_fdr",
        "n_cooccur_branches", "n_dar_only_branches", "n_contact_only_branches",
    ]

    _fmt_float = lambda v, d=3: f"{v:.{d}e}" if isinstance(v, float) else ("" if v is None else str(v))

    with open(OUT_DIR / "all_tested_pairs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=all_fields)
        w.writeheader()
        for row in sorted(all_records, key=lambda x: (_TIER_ORDER.get(x["tier"], 99), x["variant_id"])):
            row = dict(row)
            for fld in ("fisher_p", "fisher_fdr", "pagel_p", "pagel_fdr",
                        "branch_cooccur_p", "branch_cooccur_fdr"):
                row[fld] = _fmt_float(row[fld])
            for fld in ("sensitivity", "specificity"):
                row[fld] = f"{row[fld]:.3f}" if isinstance(row[fld], float) else row[fld]
            w.writerow(row)

    # ── Write compensatory_partners.csv (derived significant view) ────────────
    sig_candidates = _get_sig(all_records)
    pair_fields = [f for f in all_fields if f not in ("n_dar_only_branches", "n_contact_only_branches")]
    with open(OUT_DIR / "compensatory_partners.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=pair_fields, extrasaction="ignore")
        w.writeheader()
        for row in sorted(sig_candidates, key=lambda x: (_TIER_ORDER.get(x["tier"], 99), x["variant_id"])):
            row = dict(row)
            for fld in ("fisher_p", "fisher_fdr", "pagel_p", "pagel_fdr",
                        "branch_cooccur_p", "branch_cooccur_fdr"):
                row[fld] = _fmt_float(row[fld])
            for fld in ("sensitivity", "specificity"):
                row[fld] = f"{row[fld]:.3f}" if isinstance(row[fld], float) else row[fld]
            w.writerow(row)

    # ── Write concordance_summary.csv ─────────────────────────────────────────
    summary_fields = [
        "variant_id", "ann_id", "tier", "dar_gene", "dar_structure_gene", "dar_genome",
        "dar_aa_coord", "dar_ref_aa", "dar_alt_aa",
        "n_cdav_spp", "n_bg_spp",
        "n_contacts_tested", "n_sig_contacts",
        "n_incompatible", "intra_genomic_sig", "inter_genomic_sig",
        "sig_contact_genes",
    ]
    with open(OUT_DIR / "concordance_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary_fields)
        w.writeheader()
        w.writerows(sorted(summary_rows, key=lambda r: (_TIER_ORDER.get(r["tier"], 99), -r["n_sig_contacts"])))

    # ── Summary ───────────────────────────────────────────────────────────────
    sig_by_tier    = Counter(t["tier"] for t in sig_candidates)
    intra_by_tier  = Counter(t["tier"] for t in sig_candidates if t["contact_type"] != "mt-nuc")
    inter_by_tier  = Counter(t["tier"] for t in sig_candidates if t["contact_type"] == "mt-nuc")
    incomp_by_tier = Counter(t["tier"] for t in sig_candidates if t["likely_incompatible"] is True)

    has_pagel  = any(r["pagel_fdr"]         is not None for r in all_records)
    has_branch = any(r["branch_cooccur_fdr"] is not None for r in all_records)
    sig_method = "Fisher (fallback)" if not (has_pagel or has_branch) else \
                 "Pagel + branch co-occurrence"

    print(f"\n{'='*65}")
    print(f"COMPENSATORY PARTNER ANALYSIS  [{sig_method}]")
    print(f"{'='*65}")
    print(f"  All tested pairs  : {len(all_records)}")
    print(f"  Significant pairs : {len(sig_candidates)}")
    print(f"{'Tier':<10} {'cDAVs':>7} {'sig pairs':>10} {'intra':>7} {'inter':>7} {'incompatible':>13}")
    print("-" * 60)
    for t in sorted(by_tier, key=lambda x: _TIER_ORDER.get(x, 99)):
        n_cdavs = by_tier.get(t, 0)
        print(f"  {t:<8} {n_cdavs:>7} {sig_by_tier[t]:>10} "
              f"{intra_by_tier[t]:>7} {inter_by_tier[t]:>7} {incomp_by_tier[t]:>13}")
    print("-" * 60)
    print(f"  {'Total':<8} {len(cdavs):>7} {len(sig_candidates):>10} "
          f"{sum(intra_by_tier.values()):>7} {sum(inter_by_tier.values()):>7} "
          f"{sum(incomp_by_tier.values()):>13}")
    print(f"\nOutputs:")
    print(f"  all_tested_pairs.csv      ({len(all_records)} rows)")
    print(f"  compensatory_partners.csv ({len(sig_candidates)} rows)")
    print(f"  concordance_summary.csv   ({len(summary_rows)} rows)")


def _get_sig(records: list[dict]) -> list[dict]:
    """
    Return significant records.
    Priority: pagel_fdr OR branch_cooccur_fdr (phylogenetically valid).
    Fallback: fisher_fdr if neither phylogenetic test has run.
    """
    has_pagel  = any(r.get("pagel_fdr")         is not None for r in records)
    has_branch = any(r.get("branch_cooccur_fdr") is not None for r in records)

    def is_sig(r):
        if r.get("low_power"):
            return False
        if has_pagel or has_branch:
            p_ok = (has_pagel  and r.get("pagel_fdr")         is not None and r["pagel_fdr"]         <= 0.10)
            b_ok = (has_branch and r.get("branch_cooccur_fdr") is not None and r["branch_cooccur_fdr"] <= 0.10)
            return p_ok or b_ok
        # fallback
        fdr = r.get("fisher_fdr")
        return fdr is not None and fdr <= 0.10

    return [r for r in records if is_sig(r)]


if __name__ == "__main__":
    main()
