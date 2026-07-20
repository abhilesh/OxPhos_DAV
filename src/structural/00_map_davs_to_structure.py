#!/usr/bin/env python3
"""
src/structural/00_map_davs_to_structure.py

Map classified DAV rows to OXPHOS complex cryo-EM structures and extract
Cβ–Cβ 8 Å contact shells across a manifest-backed structure panel.

Contact definition: Cβ–Cβ ≤ 8 Å (Cα for Gly).
Contact class: hbond / electrostatic / hydrophobic / vdw.

  hbond:         any N/O heavy-atom pair ≤ 3.5 Å
  electrostatic: oppositely charged pair (K/R/H vs D/E), any heavy ≤ 5 Å
  hydrophobic:   both residues hydrophobic, sidechain C–C ≤ 5 Å
  vdw:           all other Cβ–Cβ contacts

Chain→gene assignment uses LOCAL alignment (best coverage of the shorter
sequence).  Position mapping uses GLOBAL alignment (preserves numbering
outside the high-scoring core).  After mapping, ref_aa is verified against
the PDB residue; a ±10 AA sliding window corrects for isoform offsets.

Input:
  data/derived/classified/variants_master_classified.parquet
  data/reference/structure_model_manifest.tsv
  data/structures/{PDB_ID}.cif

Output:
  results/structural/dar_structure_map.csv      — one row per (variant_id, pdb_id)
  results/structural/dar_contacts_cbcb8A.csv    — all Cβ–Cβ 8 Å contacts
  results/structural/dar_mito_nuc_contacts.csv  — cross-genome contacts only
  results/structural/structure_model_summary.csv — per-variant structure-panel summary
  results/structural/structure_transcript_reconciliation_audit.csv — transcript-space audit

Run from project root:
  python src/structural/00_map_davs_to_structure.py
"""

import csv
import json
import re
import time
import urllib.request
import urllib.error
import numpy as np
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import MMCIFParser, NeighborSearch, is_aa
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]
STRUC_DIR = ROOT / "data" / "structures"
CLASSIFIED_DIR = ROOT / "data" / "derived" / "classified"
CLASSIFIED_PARQUET = CLASSIFIED_DIR / "variants_master_classified.parquet"
STRUCTURE_MANIFEST = ROOT / "data" / "reference" / "structure_model_manifest.tsv"
if not STRUCTURE_MANIFEST.exists():
    STRUCTURE_MANIFEST = ROOT / "data" / "derived" / "reference" / "structure_model_manifest.tsv"
TRANSCRIPT_MAPS_JSON = ROOT / "data" / "derived" / "curated" / "transcript_position_maps.json"
ANCHOR_EXCEPTION_REGISTRY = ROOT / "data" / "reference" / "structural_anchor_exception_registry.tsv"
CHAIN_GENE_OVERRIDE_REGISTRY = ROOT / "data" / "reference" / "structure_chain_gene_overrides.tsv"
STRUCTURE_PANEL_ELIGIBILITY_REGISTRY = ROOT / "data" / "reference" / "structure_panel_eligibility.tsv"
TOGA_AA_DIR = ROOT / "data" / "alignments" / "toga_hg38_aa"
MT_AA_DIR   = ROOT / "data" / "alignments" / "mtdna_aa"
OUT_DIR = ROOT / "results" / "structural"

# ── Isoform proxies for structural mapping ─────────────────────────────────────
# Tissue-specific or isoform-2 subunits absent from the preferred PDB are mapped
# to their isoform-1 structural equivalent. Position numbering is shared (these
# isoforms differ mainly in N-terminal targeting sequences, not core structure).
# Status field will be tagged "proxy=<gene>" so the substitution is traceable.
_ISOFORM_PROXY = {
    "COX4I2":  "COX4I1",   # lung isoform; proxy chain present in 9I7U/9I6F
    "COX6A2":  "COX6A1",   # heart/muscle isoform; proxy chain present in 9I7U/9I6F
    "COX7A1":  "COX7A2",   # heart/muscle isoform; proxy chain present in 9I7U/9I6F
    "ATP5MC2": "ATP5MC1",  # c-subunit isoform 2; 8H9S has ATP5MC1
    "ATP5MC3": "ATP5MC1",  # c-subunit isoform 3; 8H9S has ATP5MC1
}

# ── Chain-assignment fallback threshold ────────────────────────────────────────
# When neither the RCSB API nor a manual chain-gene override resolves a chain,
# assign_chain_to_gene() falls back to local pairwise-alignment scoring. A
# chain is accepted only if its best-scoring gene clears this fraction of a
# perfect self-alignment score for the shorter sequence. 0.30 is a
# conservative, commonly-cited "safe zone" for confident local-alignment
# homology detection, but it has NOT been empirically calibrated against this
# project's own chain-assignment outcomes (flagged in the 2026-05-12
# structural audit, bug P8) — there is no labeled correct/incorrect dataset
# to calibrate against yet. Every fallback assignment is logged with its
# score ratio so a future calibration pass has real data to work from.
_CHAIN_ASSIGNMENT_IDENTITY_THRESHOLD = 0.30

# ── TOGA filename → canonical HGNC symbol ─────────────────────────────────────
# TOGA alignment files are stored under the old gene name; map them to the
# current approved symbol used everywhere else in the pipeline.
_TOGA_TO_CANONICAL = {
    # HGNC/classification uses COXFA4; PDB/UniProt resources still commonly use
    # the previous symbol NDUFA4. Keep both identities explicit in outputs.
    "COXFA4": "NDUFA4",
}
_STRUCTURE_TO_ANALYSIS = {v: k for k, v in _TOGA_TO_CANONICAL.items()}

# UniProt accessions/mnemonics observed in the active structure panel. This
# local resolver makes chain assignment reproducible even when the RCSB API is
# unavailable or omits gene names for short/repeated ATP synthase chains.
_STRUCTURE_REF_TO_GENE = {
    "P24539": "ATP5PB", "AT5F1_HUMAN": "ATP5PB",
    "P18859": "ATP5PF", "ATP5J_HUMAN": "ATP5PF",
    "O75947": "ATP5PD", "ATP5H_HUMAN": "ATP5PD",
    "P25705": "ATP5F1A", "ATPA_HUMAN": "ATP5F1A",
    "P06576": "ATP5F1B", "ATPB_HUMAN": "ATP5F1B",
    "P36542": "ATP5F1C", "ATPG_HUMAN": "ATP5F1C",
    "Q9UII2": "ATP5IF1", "ATIF1_HUMAN": "ATP5IF1",
    "P48047": "ATP5PO", "ATPO_HUMAN": "ATP5PO",
    "P05496": "ATP5MC1", "AT5G1_HUMAN": "ATP5MC1",
    "P30049": "ATP5F1D", "ATPD_HUMAN": "ATP5F1D",
    "P56381": "ATP5F1E", "ATP5E_HUMAN": "ATP5F1E",
    "P00846": "MT-ATP6", "ATP6_HUMAN": "MT-ATP6",
    "P56378": "ATP5MJ", "ATP68_HUMAN": "ATP5MJ",
    "P03928": "MT-ATP8", "ATP8_HUMAN": "MT-ATP8",
    "P56134": "ATP5MF", "ATPK_HUMAN": "ATP5MF",
    "O75964": "ATP5MG", "ATP5L_HUMAN": "ATP5MG",
    "P56385": "ATP5ME", "ATP5I_HUMAN": "ATP5ME",
    "P00395": "MT-CO1", "COX1_HUMAN": "MT-CO1",
    "P00403": "MT-CO2", "COX2_HUMAN": "MT-CO2",
    "P00414": "MT-CO3", "COX3_HUMAN": "MT-CO3",
    "P13073": "COX4I1", "COX41_HUMAN": "COX4I1",
    "P20674": "COX5A", "COX5A_HUMAN": "COX5A",
    "P10606": "COX5B", "COX5B_HUMAN": "COX5B",
    "P12074": "COX6A1", "CX6A1_HUMAN": "COX6A1",
    "P14854": "COX6B1", "CX6B1_HUMAN": "COX6B1",
    "P09669": "COX6C", "COX6C_HUMAN": "COX6C",
    "P14406": "COX7A2", "CX7A2_HUMAN": "COX7A2",
    "P24311": "COX7B", "COX7B_HUMAN": "COX7B",
    "P15954": "COX7C", "COX7C_HUMAN": "COX7C",
    "P10176": "COX8A", "COX8A_HUMAN": "COX8A",
    "O00483": "NDUFA4", "NDUA4_HUMAN": "NDUFA4",
}

# ── Gene → complex ─────────────────────────────────────────────────────────────
GENE_COMPLEX = {}
for g in [
    "MT-ND1",
    "MT-ND2",
    "MT-ND3",
    "MT-ND4",
    "MT-ND4L",
    "MT-ND5",
    "MT-ND6",
    "NDUFA1",
    "NDUFA2",
    "NDUFA3",
    "NDUFA5",
    "NDUFA6",
    "NDUFA7",
    "NDUFA8",
    "NDUFA9",
    "NDUFA10",
    "NDUFA11",
    "NDUFA12",
    "NDUFA13",
    "NDUFAB1",
    "NDUFB1",
    "NDUFB2",
    "NDUFB3",
    "NDUFB4",
    "NDUFB5",
    "NDUFB6",
    "NDUFB7",
    "NDUFB8",
    "NDUFB9",
    "NDUFB10",
    "NDUFB11",
    "NDUFC1",
    "NDUFC2",
    "NDUFS1",
    "NDUFS2",
    "NDUFS3",
    "NDUFS4",
    "NDUFS5",
    "NDUFS6",
    "NDUFS7",
    "NDUFS8",
    "NDUFV1",
    "NDUFV2",
    "NDUFV3",
]:
    GENE_COMPLEX[g] = "CI"
for g in ["SDHA", "SDHB", "SDHC", "SDHD"]:
    GENE_COMPLEX[g] = "CII"
for g in [
    "MT-CYB",
    "UQCRB",
    "UQCRC1",
    "UQCRC2",
    "UQCRFS1",
    "UQCRH",
    "UQCRQ",
    "UQCR10",
    "UQCR11",
    "CYC1",
]:
    GENE_COMPLEX[g] = "CIII"
for g in [
    "MT-CO1",
    "MT-CO2",
    "MT-CO3",
    "COX4I1",
    "COX4I2",
    "COX5A",
    "COX5B",
    "COX6A1",
    "COX6A2",
    "COX6B1",
    "COX6C",
    "COX7A1",
    "COX7A2",
    "COX7B",
    "COX7C",
    "COX8A",
    "NDUFA4",
]:
    GENE_COMPLEX[g] = "CIV"
for g in [
    "MT-ATP6",
    "MT-ATP8",
    "ATP5F1A",
    "ATP5F1B",
    "ATP5F1C",
    "ATP5F1D",
    "ATP5F1E",
    "ATP5MC1",
    "ATP5MC2",
    "ATP5MC3",
    "ATP5PB",
    "ATP5PD",
    "ATP5IF1",
    "ATP5ME",
    "ATP5MF",
    "ATP5MG",
    "ATP5MJ",
    "ATP5PF",
    "ATP5PO",
]:
    GENE_COMPLEX[g] = "CV"

MT_GENES = {g for g in GENE_COMPLEX if g.startswith("MT-")}

# ── Contact classification constants ──────────────────────────────────────────
_BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
_CHARGED_POS = {"ARG", "LYS", "HIS"}
_CHARGED_NEG = {"ASP", "GLU"}
_HYDROPHOBIC = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "CYS"}


# ── Structure helpers ──────────────────────────────────────────────────────────


def get_cb_atom(residue):
    """Return Cβ (or Cα for Gly), or None."""
    resname = residue.get_resname().strip().upper()
    try:
        return residue["CA"] if resname == "GLY" else residue["CB"]
    except KeyError:
        return None


def get_chain_residues(chain) -> list:
    return [r for r in chain if is_aa(r, standard=True)]


def get_chain_sequence(res_list: list) -> str:
    return "".join(
        protein_letters_3to1.get(r.get_resname().strip().capitalize(), "X")
        for r in res_list
    )


# ── RCSB API helpers ──────────────────────────────────────────────────────────

_RCSB_BASE = "https://data.rcsb.org/rest/v1/core"
_RCSB_CACHE: dict[str, dict] = {}  # pdb_id → {chain_id: gene_symbol}


def _rcsb_get(url: str, retries: int = 3) -> dict | None:
    """Fetch JSON from RCSB REST API with simple retry."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt < retries - 1:
                time.sleep(1.5)
        except Exception:
            if attempt < retries - 1:
                time.sleep(1.5)
    return None


def fetch_chain_gene_map(pdb_id: str, known_genes: set) -> dict[str, str]:
    """
    Query RCSB for {chain_id: gene_symbol} for a given PDB entry.

    Steps:
      1. Fetch entry → polymer entity IDs.
      2. Per entity: get gene name(s) from rcsb_entity_source_organism and
         chain IDs from entity_poly.pdbx_strand_id.
      3. Resolve the gene name against known_genes (exact → prefix → skip).

    Returns an empty dict on network failure (caller falls back to alignment).
    """
    pdb_upper = pdb_id.upper()
    if pdb_upper in _RCSB_CACHE:
        return _RCSB_CACHE[pdb_upper]

    result: dict[str, str] = {}

    # Step 1: entry → entity IDs
    entry = _rcsb_get(f"{_RCSB_BASE}/entry/{pdb_upper}")
    if not entry:
        _RCSB_CACHE[pdb_upper] = result
        return result

    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "polymer_entity_ids", []
    )

    # Step 2: per entity → gene names + chain IDs
    for eid in entity_ids:
        entity = _rcsb_get(f"{_RCSB_BASE}/polymer_entity/{pdb_upper}/{eid}")
        if not entity:
            continue

        # Chain IDs for this entity (comma-separated string in pdbx_strand_id)
        strand_str = entity.get("entity_poly", {}).get("pdbx_strand_id", "") or ""
        chain_ids = [c.strip() for c in strand_str.split(",") if c.strip()]

        # Gene names from source organism annotation
        gene_names: list[str] = []
        for org in entity.get("rcsb_entity_source_organism", []):
            for gn in org.get("rcsb_gene_name", []):
                val = gn.get("value", "")
                if val:
                    gene_names.append(val.upper())

        # Also check the entity description as a fallback name source
        desc = entity.get("rcsb_polymer_entity", {}).get("pdbx_description", "") or ""

        # Resolve to a known gene symbol
        resolved_gene = None
        for gname in gene_names:
            if gname in known_genes:
                resolved_gene = gname
                break
        if resolved_gene is None:
            # Try prefix match (e.g. RCSB returns "MT-CO1" vs our "MT-CO1")
            for gname in gene_names:
                for kg in known_genes:
                    if kg == gname or kg.upper() == gname:
                        resolved_gene = kg
                        break
                if resolved_gene:
                    break
        if resolved_gene is None:
            # Last resort: scan description for a known gene token
            desc_upper = desc.upper()
            for kg in sorted(known_genes, key=len, reverse=True):
                if re.search(rf"\b{re.escape(kg)}\b", desc_upper):
                    resolved_gene = kg
                    break

        if resolved_gene:
            for cid in chain_ids:
                result[cid] = resolved_gene

    _RCSB_CACHE[pdb_upper] = result
    return result


def _as_list(value) -> list:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def load_chain_gene_overrides() -> dict[tuple[str, str], str]:
    if not CHAIN_GENE_OVERRIDE_REGISTRY.exists():
        return {}
    overrides: dict[tuple[str, str], str] = {}
    with open(CHAIN_GENE_OVERRIDE_REGISTRY, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pdb_id = (row.get("pdb_id") or "").strip().upper()
            chain_id = (row.get("chain_id") or "").strip()
            structure_gene = (row.get("structure_gene") or "").strip()
            if pdb_id and chain_id and structure_gene:
                overrides[(pdb_id, chain_id)] = structure_gene
    return overrides


def load_structure_panel_eligibility() -> dict[str, str]:
    """
    Load the structure-panel eligibility registry: per-gene ground truth on
    whether structural absence reflects a real panel gap (`structurally_unrepresented`),
    no classified/eligible variants ever reaching this stage (`no_eligible_variants`),
    or normal direct/proxy coverage. Built from the historical mapping run in
    data/reference/structure_panel_eligibility.tsv (see tools/ for regeneration).
    Missing file degrades gracefully to "no prior evidence" for every gene.
    """
    if not STRUCTURE_PANEL_ELIGIBILITY_REGISTRY.exists():
        return {}
    eligibility: dict[str, str] = {}
    with open(STRUCTURE_PANEL_ELIGIBILITY_REGISTRY, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = (row.get("gene") or "").strip()
            status = (row.get("eligibility_status") or "").strip()
            if gene and status:
                eligibility[gene] = status
    return eligibility


def resolve_structure_ref_gene(accession: str, db_code: str, description: str, known_genes: set) -> str | None:
    for key in (accession, db_code):
        gene = _STRUCTURE_REF_TO_GENE.get(str(key).strip().upper())
        if gene:
            return gene if gene in known_genes else None
    desc_upper = str(description or "").upper()
    for gene in sorted(known_genes, key=len, reverse=True):
        if re.search(rf"\b{re.escape(gene)}\b", desc_upper):
            return gene
    return None


def load_mmcif_chain_gene_map(cif_path: Path, known_genes: set) -> dict[str, str]:
    """
    Build Bio.PDB chain-id → gene mapping from local mmCIF metadata.

    Bio.PDB exposes author chain IDs by default. The most reliable local route is
    auth_asym_id → label_asym_id → entity_id → UniProt/description → gene.
    """
    d = MMCIF2Dict(str(cif_path))

    descriptions = {
        entity_id: desc
        for entity_id, desc in zip(
            _as_list(d.get("_entity.id")),
            _as_list(d.get("_entity.pdbx_description")),
        )
    }

    entity_to_gene: dict[str, str] = {}
    for ref_id, entity_id, accession, db_code in zip(
        _as_list(d.get("_struct_ref.id")),
        _as_list(d.get("_struct_ref.entity_id")),
        _as_list(d.get("_struct_ref.pdbx_db_accession")),
        _as_list(d.get("_struct_ref.db_code")),
    ):
        gene = resolve_structure_ref_gene(
            accession,
            db_code,
            descriptions.get(entity_id, ""),
            known_genes,
        )
        if gene:
            entity_to_gene[entity_id] = gene

    label_to_entity = {
        label_asym: entity_id
        for label_asym, entity_id in zip(
            _as_list(d.get("_struct_asym.id")),
            _as_list(d.get("_struct_asym.entity_id")),
        )
    }

    auth_to_label: dict[str, str] = {}
    for label_asym, auth_asym in zip(
        _as_list(d.get("_atom_site.label_asym_id")),
        _as_list(d.get("_atom_site.auth_asym_id")),
    ):
        auth_to_label.setdefault(auth_asym, label_asym)

    chain_map: dict[str, str] = {}
    for auth_asym, label_asym in auth_to_label.items():
        gene = entity_to_gene.get(label_to_entity.get(label_asym, ""))
        if gene:
            chain_map[auth_asym] = gene
            chain_map[label_asym] = gene
    return chain_map


# ── Alignment helpers ──────────────────────────────────────────────────────────


def assign_chain_to_gene(ref_seq: str, chain_seq: str) -> float:
    """Local alignment score — fallback chain→gene identification."""
    aln = PairwiseAligner()
    aln.mode = "local"
    aln.match_score = 2
    aln.mismatch_score = -1
    aln.open_gap_score = -3
    aln.extend_gap_score = -0.5
    try:
        return next(iter(aln.align(ref_seq, chain_seq))).score
    except StopIteration:
        return 0.0


def map_refseq_to_chain(ref_seq: str, chain_seq: str) -> dict:
    """
    Global alignment → {refseq_1based_pos: chain_0based_idx}.
    Global mode maps positions outside the high-scoring core, preventing
    silent off-by-one errors when the RefSeq and PDB sequence differ slightly.
    """
    aln = PairwiseAligner()
    aln.mode = "global"
    aln.match_score = 2
    aln.mismatch_score = -1
    aln.open_gap_score = -5
    aln.extend_gap_score = -0.5
    try:
        best = next(iter(aln.align(ref_seq, chain_seq)))
    except StopIteration:
        return {}

    pos_map = {}
    ref_pos = chain_pos = 0
    for r_aa, c_aa in zip(best[0], best[1]):
        if r_aa != "-":
            ref_pos += 1
        if c_aa != "-":
            chain_pos += 1
        if r_aa != "-" and c_aa != "-":
            pos_map[ref_pos] = chain_pos - 1
    return pos_map


def find_anchor(
    pos_map: dict,
    aa_coord: int,
    ref_aa: str,
    ref_seq: str,
    window: int = 10,
) -> int | None:
    """
    Return the RefSeq position whose mapped PDB residue matches ref_aa.

    1. Direct lookup at aa_coord — fastest path.
    2. Slide ±window AA to correct for isoform N-terminal offsets.
    Returns None if the residue cannot be confidently anchored.
    """
    if aa_coord in pos_map and 1 <= aa_coord <= len(ref_seq):
        if ref_seq[aa_coord - 1] == ref_aa:
            return aa_coord

    for delta in range(1, window + 1):
        for candidate in (aa_coord - delta, aa_coord + delta):
            if 1 <= candidate <= len(ref_seq) and candidate in pos_map:
                if ref_seq[candidate - 1] == ref_aa:
                    return candidate
    return None


def nearest_mapped_positions(pos_map: dict, aa_coord: int, limit: int = 2) -> list[int]:
    if not pos_map:
        return []
    ranked = sorted(pos_map.keys(), key=lambda pos: (abs(pos - aa_coord), pos))
    return ranked[:limit]


def classify_anchor_failure(
    pos_map: dict,
    aa_coord: int,
    ref_aa: str,
    ref_seq: str,
    default_anchor: int | None,
    extended_anchor: int | None,
    default_window: int,
    extended_window: int,
) -> tuple[str, str]:
    if default_anchor is not None:
        return "anchored_default_window", "default_window"
    if extended_anchor is not None:
        delta = extended_anchor - aa_coord
        return f"large_offset_candidate_{delta:+d}", "extended_window"

    if aa_coord not in pos_map:
        nearest = nearest_mapped_positions(pos_map, aa_coord, limit=2)
        if nearest:
            nearest_delta = min(abs(pos - aa_coord) for pos in nearest)
            if nearest_delta <= 3:
                return "unresolved_structure_segment_candidate", "nearest_mapped_gap"
            if nearest_delta <= extended_window:
                return "possible_large_offset_or_gap", "nearest_mapped_gap"
        return "unmapped_reference_position", "position_absent_from_pos_map"

    expected = ref_seq[aa_coord - 1] if 1 <= aa_coord <= len(ref_seq) else None
    if expected != ref_aa:
        return "reference_sequence_conflict", "refseq_residue_mismatch"

    return "anchor_window_exhausted", "window_exhausted"


# ── Variant field helpers ──────────────────────────────────────────────────────


def parse_aa_coord(aa_change: str) -> int | None:
    """Extract integer AA position from 'S45F' → 45."""
    m = re.search(r"[A-Za-z]+(\d+)[A-Za-z]+", aa_change or "")
    return int(m.group(1)) if m else None


def parse_variant_aa_coord(record: dict) -> int | None:
    coord = record.get("frame_specific_codon_index")
    if coord not in (None, ""):
        try:
            return int(coord)
        except (TypeError, ValueError):
            pass
    for field in ("frame_specific_hgvs_p", "hgvs_p", "aa_change"):
        aa_coord = parse_aa_coord(str(record.get(field, "") or ""))
        if aa_coord is not None:
            return aa_coord
    return None


def clean_value(value):
    return None if pd.isna(value) else value


def first_non_null(*values):
    for value in values:
        value = clean_value(value)
        if value not in (None, ""):
            return value
    return None


def load_structure_manifest() -> list[dict]:
    if not STRUCTURE_MANIFEST.exists():
        raise FileNotFoundError(f"Missing structure model manifest: {STRUCTURE_MANIFEST}")
    rows: list[dict] = []
    with open(STRUCTURE_MANIFEST, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if str(row.get("active", "")).strip().lower() != "true":
                continue
            row["priority"] = int(row.get("priority") or 999)
            row["resolution_angstrom"] = (
                float(row["resolution_angstrom"]) if row.get("resolution_angstrom") else None
            )
            rows.append(row)
    rows.sort(key=lambda item: (item["complex_id"], item["priority"], item["pdb_id"]))
    return rows


def load_transcript_maps() -> dict[str, dict]:
    if not TRANSCRIPT_MAPS_JSON.exists():
        return {}
    with open(TRANSCRIPT_MAPS_JSON, encoding="utf-8") as handle:
        raw = json.load(handle)
    tx_maps: dict[str, dict] = {}
    for gene, entry in raw.items():
        if entry.get("type") not in {"mapped", "identity"}:
            continue
        pos_map = {}
        for key, value in (entry.get("map", {}) or {}).items():
            if value in (None, ""):
                continue
            try:
                pos_map[int(key)] = int(value)
            except (TypeError, ValueError):
                continue
        tx_maps[gene] = {
            "map": pos_map,
            "type": entry.get("type"),
            "identity_fraction": entry.get("identity_fraction"),
            "coverage_fraction": entry.get("coverage_fraction"),
            "nm": entry.get("nm"),
            "enst": entry.get("enst"),
        }
    return tx_maps


def load_anchor_exception_registry() -> dict[str, dict]:
    if not ANCHOR_EXCEPTION_REGISTRY.exists():
        return {}
    registry = {}
    with open(ANCHOR_EXCEPTION_REGISTRY, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = (row.get("gene") or "").strip()
            if not gene:
                continue
            registry[gene] = {
                "allow_extended_anchor": str(row.get("allow_extended_anchor", "")).strip().lower() == "true",
                "max_offset": int(row.get("max_offset") or 10),
                "rationale": row.get("rationale", ""),
                "evidence_status": row.get("evidence_status", ""),
                "notes": row.get("notes", ""),
            }
    return registry


def anchor_policy_for_gene(gene: str, registry: dict[str, dict]) -> dict:
    policy = registry.get(gene, {})
    allow_extended = bool(policy.get("allow_extended_anchor"))
    max_offset = int(policy.get("max_offset") or 10)
    if max_offset < 10:
        max_offset = 10
    return {
        "allow_extended_anchor": allow_extended,
        "max_offset": max_offset,
        "policy_label": "registry_enabled" if allow_extended else "default_only",
    }


def remap_aa_coord_to_structure_space(
    gene: str,
    genome: str,
    aa_coord: int | None,
    tx_maps: dict[str, dict],
) -> tuple[int | None, str | None, str | None]:
    if aa_coord is None:
        return None, None, None
    if genome != "nucDNA":
        return aa_coord, "native", None
    tx_entry = tx_maps.get(gene)
    if not tx_entry:
        return aa_coord, "anchor_fallback", "no_transcript_map"
    tx_map = tx_entry.get("map") or {}
    mapped = tx_map.get(int(aa_coord))
    if mapped is None:
        return None, "transcript_map", "position_not_in_enst"
    if mapped == aa_coord:
        return mapped, "transcript_identity", None
    return mapped, "transcript_map", None


def build_transcript_reconciliation(
    record: dict,
    gene: str,
    genome: str,
    aa_coord_raw: int | None,
    aa_coord_mapped: int | None,
    aa_coord_method: str | None,
    aa_coord_status: str | None,
    tx_maps: dict[str, dict],
) -> dict:
    if genome != "nucDNA":
        return {
            "clinvar_transcript_id": clean_value(record.get("transcript_id")),
            "preferred_nm": None,
            "toga_enst": None,
            "transcript_map_type": None,
            "transcript_map_identity": None,
            "transcript_map_coverage": None,
            "transcript_reconciliation_status": "mt_native",
            "transcript_coord_delta": None,
        }

    tx_entry = tx_maps.get(gene) or {}
    if aa_coord_status == "position_not_in_enst":
        reconciliation_status = "position_not_in_enst"
    elif aa_coord_method == "transcript_identity":
        reconciliation_status = "transcript_identity"
    elif aa_coord_method == "transcript_map":
        reconciliation_status = "transcript_remapped"
    elif aa_coord_method == "anchor_fallback":
        reconciliation_status = "no_transcript_map"
    else:
        reconciliation_status = "native"

    delta = None
    if aa_coord_raw is not None and aa_coord_mapped is not None:
        delta = int(aa_coord_mapped) - int(aa_coord_raw)

    return {
        "clinvar_transcript_id": clean_value(record.get("transcript_id")),
        "preferred_nm": clean_value(tx_entry.get("nm")),
        "toga_enst": clean_value(tx_entry.get("enst")),
        "transcript_map_type": clean_value(tx_entry.get("type")),
        "transcript_map_identity": clean_value(tx_entry.get("identity_fraction")),
        "transcript_map_coverage": clean_value(tx_entry.get("coverage_fraction")),
        "transcript_reconciliation_status": reconciliation_status,
        "transcript_coord_delta": delta,
    }


def status_category(status: str, proxy_gene: str = "", gene_eligibility: str = "") -> str:
    if status == "structural_ineligible":
        return "structural_ineligible"
    if status == "secondary_not_attempted_no_primary_support":
        return "secondary_not_attempted"
    if status == "missing_structure_file":
        return "missing_model_file"
    if status == "position_not_in_enst":
        return "sequence_model_mismatch"
    if status == "no_chain_in_pdb":
        if proxy_gene:
            return "isoform_proxy_gap"
        # Split the previously-conflated "chain_assignment_or_model_gap" bucket
        # (bug P9, 2026-05-12 audit) using the structure-panel eligibility
        # registry: a gene with no historical successful mapping anywhere is a
        # known panel gap, not an unexpected per-variant failure.
        if gene_eligibility == "structurally_unrepresented":
            return "gene_unrepresented_in_panel"
        return "chain_assignment_or_model_gap"
    if status in {
        "unresolved_structure_segment_candidate",
        "possible_large_offset_or_gap",
        "unmapped_reference_position",
        "anchor_window_exhausted",
        "reference_sequence_conflict",
    }:
        return "residue_anchoring_failure"
    if status.startswith("large_offset_candidate_"):
        return "mature_offset_candidate"
    if status.startswith("extended_offset_rescue_"):
        return "mapped_with_extended_offset_rescue"
    if status == "position_not_anchored":
        return "residue_anchoring_failure"
    if status.startswith("isoform_offset"):
        return "mapped_with_isoform_offset"
    if status == "ok":
        return "mapped_direct"
    return "other"


# ── Contact classification ─────────────────────────────────────────────────────


def classify_contact(res_a, res_b) -> str:
    """
    Priority: hbond > electrostatic > hydrophobic > vdw.
    """
    name_a = res_a.get_resname().strip().upper()
    name_b = res_b.get_resname().strip().upper()
    atoms_a = list(res_a.get_atoms())
    atoms_b = list(res_b.get_atoms())

    # 1. H-bond: any N/O pair ≤ 3.5 Å
    polar_a = [a for a in atoms_a if a.element in ("N", "O")]
    polar_b = [a for a in atoms_b if a.element in ("N", "O")]
    for pa in polar_a:
        for pb in polar_b:
            if np.linalg.norm(pa.coord - pb.coord) <= 3.5:
                return "hbond"

    # 2. Electrostatic: oppositely charged, any heavy ≤ 5 Å
    if (name_a in _CHARGED_POS and name_b in _CHARGED_NEG) or (
        name_a in _CHARGED_NEG and name_b in _CHARGED_POS
    ):
        for aa in atoms_a:
            for ab in atoms_b:
                if np.linalg.norm(aa.coord - ab.coord) <= 5.0:
                    return "electrostatic"

    # 3. Hydrophobic: both nonpolar, sidechain C–C ≤ 5 Å
    if name_a in _HYDROPHOBIC and name_b in _HYDROPHOBIC:
        sc_a = [
            a for a in atoms_a if a.element == "C" and a.name not in _BACKBONE_ATOMS
        ]
        sc_b = [
            a for a in atoms_b if a.element == "C" and a.name not in _BACKBONE_ATOMS
        ]
        for ca in sc_a:
            for cb in sc_b:
                if np.linalg.norm(ca.coord - cb.coord) <= 5.0:
                    return "hydrophobic"

    return "vdw"


# ── Main ───────────────────────────────────────────────────────────────────────


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cif_parser = MMCIFParser(QUIET=True)
    structure_manifest = load_structure_manifest()
    transcript_maps = load_transcript_maps()
    anchor_registry = load_anchor_exception_registry()
    chain_gene_overrides = load_chain_gene_overrides()
    panel_eligibility = load_structure_panel_eligibility()
    manifest_by_complex: dict[str, list[dict]] = {}
    for row in structure_manifest:
        manifest_by_complex.setdefault(row["complex_id"], []).append(row)

    # Load reference protein sequences from per-gene alignment FASTAs.
    # The human reference (Homo_sapiens|9606) is extracted from each alignment
    # and stored gap-free as the canonical RefSeq for position mapping.
    known_genes = set(GENE_COMPLEX.keys())
    ref_seqs = {}

    for aln_dir, genome_label in [(TOGA_AA_DIR, "nucDNA"), (MT_AA_DIR, "mtDNA")]:
        if not aln_dir.exists():
            print(f"  WARNING: alignment dir not found: {aln_dir}")
            continue
        for fasta in sorted(aln_dir.glob("*_aa_alignment.fasta")):
            # Gene name is the file stem before _aa_alignment;
            # remap TOGA old names to current canonical HGNC symbols.
            gene = fasta.name.replace("_aa_alignment.fasta", "")
            gene = _TOGA_TO_CANONICAL.get(gene, gene)
            if gene not in known_genes:
                continue
            # Extract gap-free human reference sequence
            for rec in SeqIO.parse(fasta, "fasta"):
                if rec.id.startswith("Homo_sapiens|9606"):
                    ref_seqs[gene] = str(rec.seq).replace("-", "")
                    break

    print(f"Loaded {len(ref_seqs)} reference protein sequences.")
    print(f"Loaded {len(transcript_maps)} transcript position maps.")
    print(f"Loaded {len(anchor_registry)} anchor-exception registry entries.")
    print(f"Loaded {len(chain_gene_overrides)} chain-gene override entries.")

    if not CLASSIFIED_PARQUET.exists():
        raise FileNotFoundError(f"Missing classified master table: {CLASSIFIED_PARQUET}")
    classified_df = pd.read_parquet(CLASSIFIED_PARQUET)
    classified_records = classified_df.to_dict(orient="records")

    structural_rows = []
    for r in classified_records:
        interpreted_gene = r.get("interpreted_gene") or r.get("classification_gene") or r.get("locus", "")
        structure_gene = _TOGA_TO_CANONICAL.get(interpreted_gene, interpreted_gene)
        complex_id = GENE_COMPLEX.get(structure_gene)
        aa_coord_raw = parse_variant_aa_coord(r)
        aa_coord, aa_coord_method, aa_coord_status = remap_aa_coord_to_structure_space(
            interpreted_gene,
            r.get("genome", ""),
            aa_coord_raw,
            transcript_maps,
        )
        transcript_info = build_transcript_reconciliation(
            r,
            interpreted_gene,
            r.get("genome", ""),
            aa_coord_raw,
            aa_coord,
            aa_coord_method,
            aa_coord_status,
            transcript_maps,
        )
        anchor_policy = anchor_policy_for_gene(structure_gene, anchor_registry)
        structural_rows.append({
            "variant_id": r.get("variant_id", ""),
            "ann_id": r.get("ann_id", ""),
            "interpreted_gene": interpreted_gene,
            "structure_gene": structure_gene,
            "complex_id": complex_id,
            "genome": r.get("genome", ""),
            "tier": r.get("tier", ""),
            "nc_change": r.get("nc_change", ""),
            "aa_change": r.get("aa_change", ""),
            "ref_aa": first_non_null(r.get("frame_specific_ref_aa"), r.get("ref_aa")) or "",
            "alt_aa": first_non_null(r.get("frame_specific_alt_aa"), r.get("alt_aa")) or "",
            "is_cdav_amino_acid": bool(r.get("is_cdav_amino_acid")),
            "is_cdav_nucleotide": bool(r.get("is_cdav_nucleotide")),
            "n_species_with_disease_allele": r.get("n_species_with_disease_allele", 0),
            "classification_status": r.get("classification_status"),
            "classification_subset": clean_value(r.get("classification_subset")),
            "exception_applied": bool(r.get("exception_applied")),
            "classification_warning_reason": clean_value(r.get("classification_warning_reason")),
            "aa_coord_raw": aa_coord_raw,
            "aa_coord": aa_coord,
            "aa_coord_method": aa_coord_method,
            "aa_coord_status": aa_coord_status,
            "classification_coordinate_method": clean_value(r.get("classification_coordinate_method")),
            "classification_coordinate_status": clean_value(r.get("classification_coordinate_status")),
            **transcript_info,
            "anchor_policy_label": anchor_policy["policy_label"],
            "anchor_allow_extended": anchor_policy["allow_extended_anchor"],
            "anchor_max_offset": anchor_policy["max_offset"],
            "anchor_method": None,
            "anchor_failure_detail": None,
            "structure_mapping_eligible": bool(
                r.get("classification_status") == "classified"
                and structure_gene in ref_seqs
                and complex_id in manifest_by_complex
                and aa_coord is not None
            ),
        })
    print(
        f"Loaded {len(structural_rows)} classified rows "
        f"({sum(1 for d in structural_rows if d['genome']=='mtDNA')} mtDNA, "
        f"{sum(1 for d in structural_rows if d['genome']=='nucDNA')} nucDNA)."
    )

    # ── Per-PDB structure cache ────────────────────────────────────────────────
    # Built once per PDB, reused for all variants in that complex.
    struct_cache = {}

    map_rows = []  # one row per DAR — position mapping result + status
    contact_rows = []  # one row per (DAR, contact residue)
    for row in structural_rows:
        gene = row["interpreted_gene"]
        structure_gene = row["structure_gene"]
        complex_id = row["complex_id"]
        if not row["structure_mapping_eligible"]:
            map_rows.append({
                **row,
                "pdb_id": "",
                "model_role": "",
                "state_label": "",
                "chain_id": "",
                "struct_resnum": "",
                "status": "structural_ineligible",
                "proxy_gene": "",
                "status_category": status_category("structural_ineligible"),
                "anchor_method": None,
                "anchor_failure_detail": None,
            })
            continue

        model_rows = manifest_by_complex.get(complex_id, [])
        primary_models = [m for m in model_rows if m.get("role") == "primary"]
        secondary_models = [m for m in model_rows if m.get("role") != "primary"]
        ordered_models = primary_models + secondary_models
        primary_mapped_any = False

        for model_meta in ordered_models:
            pdb_id = model_meta["pdb_id"]
            if model_meta.get("role") != "primary" and not primary_mapped_any:
                map_rows.append({
                    **row,
                    "pdb_id": pdb_id,
                    "model_role": model_meta["role"],
                    "state_label": model_meta["state_label"],
                    "chain_id": "",
                    "struct_resnum": "",
                    "status": "secondary_not_attempted_no_primary_support",
                    "proxy_gene": "",
                    "status_category": status_category("secondary_not_attempted_no_primary_support"),
                    "anchor_method": None,
                    "anchor_failure_detail": "primary_support_required",
                })
                continue

            cif_path = STRUC_DIR / f"{pdb_id}.cif"
            if not cif_path.exists():
                map_rows.append({
                    **row,
                    "pdb_id": pdb_id,
                    "model_role": model_meta["role"],
                    "state_label": model_meta["state_label"],
                        "chain_id": "",
                        "struct_resnum": "",
                        "status": "missing_structure_file",
                        "proxy_gene": "",
                        "status_category": status_category("missing_structure_file"),
                        "anchor_method": None,
                        "anchor_failure_detail": None,
                    })
                continue

            if pdb_id not in struct_cache:
                print(f"\nIndexing {pdb_id} ({complex_id})...")
                model = cif_parser.get_structure(pdb_id, str(cif_path))[0]

                chain_to_gene = {}
                gene_to_chains = {}
                chain_residues = {}
                chain_seqs = {}
                cb_atoms = []

                for chain in model.get_chains():
                    res_list = get_chain_residues(chain)
                    if not res_list:
                        continue
                    chain_residues[chain.id] = res_list
                    chain_seqs[chain.id] = get_chain_sequence(res_list)
                    cb_atoms.extend(a for r in res_list for a in [get_cb_atom(r)] if a is not None)

                known_structure_genes = set(ref_seqs.keys())
                mmcif_map = load_mmcif_chain_gene_map(cif_path, known_structure_genes)
                if mmcif_map:
                    print(f"  Using local mmCIF chain-gene map ({len(mmcif_map)} chain IDs annotated).")
                    for chain_id, mapped_gene in mmcif_map.items():
                        if chain_id in chain_seqs and mapped_gene in ref_seqs:
                            chain_to_gene[chain_id] = mapped_gene
                            gene_to_chains.setdefault(mapped_gene, []).append(chain_id)

                for (override_pdb, override_chain), mapped_gene in chain_gene_overrides.items():
                    if override_pdb == pdb_id.upper() and override_chain in chain_seqs and mapped_gene in ref_seqs:
                        chain_to_gene[override_chain] = mapped_gene
                        gene_to_chains.setdefault(mapped_gene, []).append(override_chain)

                api_map = fetch_chain_gene_map(pdb_id, known_structure_genes)
                if api_map:
                    print(f"  Using RCSB API gene-chain map ({len(api_map)} chains annotated).")
                    for chain_id, mapped_gene in api_map.items():
                        if chain_id in chain_seqs and mapped_gene in ref_seqs and chain_id not in chain_to_gene:
                            chain_to_gene[chain_id] = mapped_gene
                            gene_to_chains.setdefault(mapped_gene, []).append(chain_id)
                elif not mmcif_map:
                    print("  RCSB API unavailable — falling back to local alignment.")

                for chain_id, chain_seq in chain_seqs.items():
                    if chain_id in chain_to_gene:
                        continue
                    best_gene, best_score, best_ratio = None, 0.0, 0.0
                    for r_gene, r_seq in ref_seqs.items():
                        if GENE_COMPLEX.get(r_gene) != complex_id:
                            continue
                        score = assign_chain_to_gene(r_seq, chain_seq)
                        max_possible = min(len(r_seq), len(chain_seq)) * 2
                        threshold = max_possible * _CHAIN_ASSIGNMENT_IDENTITY_THRESHOLD
                        if score > best_score and score > threshold:
                            best_score = score
                            best_gene = r_gene
                            best_ratio = score / max_possible if max_possible else 0.0
                    if best_gene:
                        chain_to_gene[chain_id] = best_gene
                        gene_to_chains.setdefault(best_gene, []).append(chain_id)
                        print(f"    Fallback local-alignment chain assignment: "
                              f"{pdb_id} chain {chain_id} -> {best_gene} "
                              f"(score_ratio={best_ratio:.2f}, "
                              f"threshold={_CHAIN_ASSIGNMENT_IDENTITY_THRESHOLD:.2f})")

                pos_maps = {}
                for mapped_gene, chain_ids in gene_to_chains.items():
                    r_seq = ref_seqs.get(mapped_gene, "")
                    for chain_id in chain_ids:
                        chain_seq = chain_seqs.get(chain_id, "")
                        if r_seq and chain_seq:
                            pos_maps[(mapped_gene, chain_id)] = map_refseq_to_chain(r_seq, chain_seq)

                struct_cache[pdb_id] = {
                    "chain_to_gene": chain_to_gene,
                    "gene_to_chains": gene_to_chains,
                    "pos_maps": pos_maps,
                    "chain_residues": chain_residues,
                    "ns": NeighborSearch(cb_atoms),
                }

            cache = struct_cache[pdb_id]
            aa_coord = row["aa_coord"]
            ref_aa = row["ref_aa"]

            proxy_gene = None
            if structure_gene not in cache["gene_to_chains"]:
                proxy_gene = _ISOFORM_PROXY.get(structure_gene)
                if not (proxy_gene and proxy_gene in cache["gene_to_chains"]):
                    map_rows.append({
                        **row,
                        "pdb_id": pdb_id,
                        "model_role": model_meta["role"],
                        "state_label": model_meta["state_label"],
                        "chain_id": "",
                        "struct_resnum": "",
                        "status": "no_chain_in_pdb",
                        "proxy_gene": proxy_gene or "",
                        "status_category": status_category(
                            "no_chain_in_pdb", proxy_gene or "",
                            panel_eligibility.get(row.get("interpreted_gene", ""), ""),
                        ),
                        "anchor_method": None,
                        "anchor_failure_detail": None,
                    })
                    continue

            lookup_gene = proxy_gene if proxy_gene else structure_gene
            anchor_policy = anchor_policy_for_gene(lookup_gene, anchor_registry)
            candidate = None
            failure_evidence = None
            for chain_id in cache["gene_to_chains"].get(lookup_gene, []):
                pos_map = cache["pos_maps"].get((lookup_gene, chain_id), {})
                ref_seq = ref_seqs.get(lookup_gene, "")
                default_anchor = find_anchor(pos_map, aa_coord, ref_aa, ref_seq, window=10)
                extended_anchor = None
                if (
                    default_anchor is None
                    and anchor_policy["allow_extended_anchor"]
                    and anchor_policy["max_offset"] > 10
                ):
                    extended_anchor = find_anchor(
                        pos_map, aa_coord, ref_aa, ref_seq, window=anchor_policy["max_offset"]
                    )
                diagnostic_anchor = None if default_anchor is not None else find_anchor(
                    pos_map, aa_coord, ref_aa, ref_seq, window=80
                )
                anchor_status, anchor_method = classify_anchor_failure(
                    pos_map,
                    aa_coord,
                    ref_aa,
                    ref_seq,
                    default_anchor,
                    diagnostic_anchor,
                    10,
                    anchor_policy["max_offset"],
                )
                true_coord = default_anchor if default_anchor is not None else extended_anchor
                if true_coord is None:
                    if failure_evidence is None or (
                        anchor_status.startswith("large_offset_candidate_")
                        and not str(failure_evidence[0]).startswith("large_offset_candidate_")
                    ):
                        failure_evidence = (anchor_status, anchor_method, chain_id)
                    continue
                struct_idx = pos_map[true_coord]
                res_list = cache["chain_residues"][chain_id]
                if struct_idx >= len(res_list):
                    continue
                target_res = res_list[struct_idx]
                struct_aa = protein_letters_3to1.get(
                    target_res.get_resname().strip().capitalize(), "?"
                )
                if struct_aa != ref_aa:
                    continue
                offset = abs(true_coord - aa_coord)
                candidate = (
                    offset,
                    chain_id,
                    true_coord,
                    struct_idx,
                    target_res,
                    anchor_method,
                ) if candidate is None or offset < candidate[0] else candidate

            if candidate is None:
                failed_status = "position_not_in_enst" if row.get("aa_coord_status") == "position_not_in_enst" else "position_not_anchored"
                anchor_method = None
                anchor_failure_detail = None
                if row.get("aa_coord_status") != "position_not_in_enst" and failure_evidence is not None:
                    failed_status, anchor_failure_detail, _ = failure_evidence
                elif row.get("aa_coord_status") == "position_not_in_enst":
                    anchor_failure_detail = "transcript_position_missing"
                map_rows.append({
                    **row,
                    "pdb_id": pdb_id,
                    "model_role": model_meta["role"],
                    "state_label": model_meta["state_label"],
                    "chain_id": "",
                    "struct_resnum": "",
                    "status": failed_status,
                    "proxy_gene": proxy_gene or "",
                    "status_category": status_category(failed_status, proxy_gene or ""),
                    "anchor_method": anchor_method,
                    "anchor_failure_detail": anchor_failure_detail,
                })
                continue

            offset, chain_id, true_coord, struct_idx, target_res, anchor_method = candidate
            struct_resnum = target_res.get_id()[1]
            coord_delta = true_coord - aa_coord
            if coord_delta != 0:
                status = (
                    f"extended_offset_rescue_{coord_delta:+d}"
                    if anchor_method == "extended_window"
                    else f"isoform_offset_{coord_delta:+d}"
                )
            else:
                status = "ok"
            map_rows.append({
                **row,
                "pdb_id": pdb_id,
                "model_role": model_meta["role"],
                "state_label": model_meta["state_label"],
                "chain_id": chain_id,
                "struct_resnum": struct_resnum,
                "status": status,
                "proxy_gene": proxy_gene or "",
                "status_category": status_category(status, proxy_gene or ""),
                "anchor_method": anchor_method,
                "anchor_failure_detail": None,
            })
            if model_meta.get("role") == "primary":
                primary_mapped_any = True

            target_cb = get_cb_atom(target_res)
            if target_cb is None:
                continue

            for n_res in cache["ns"].search(target_cb.coord, 8.0, "R"):
                if n_res is target_res or not is_aa(n_res, standard=True):
                    continue

                n_chain_id = n_res.get_parent().id
                n_structure_gene = cache["chain_to_gene"].get(n_chain_id, f"Unknown({n_chain_id})")
                n_gene = _STRUCTURE_TO_ANALYSIS.get(n_structure_gene, n_structure_gene)
                contact_is_mt = n_structure_gene in MT_GENES
                dar_is_mt = gene in MT_GENES
                contact_type = (
                    "mt-mt"
                    if dar_is_mt and contact_is_mt
                    else "nuc-nuc" if not dar_is_mt and not contact_is_mt else "mt-nuc"
                )

                n_res_list = cache["chain_residues"].get(n_chain_id, [])
                try:
                    n_chain_idx = n_res_list.index(n_res)
                except ValueError:
                    n_chain_idx = -1
                n_pos_map = cache["pos_maps"].get((n_structure_gene, n_chain_id), {})
                rev_map = {v: k for k, v in n_pos_map.items()}
                contact_refseq_pos = rev_map.get(n_chain_idx, "")

                contact_rows.append({
                    "variant_id": row["variant_id"],
                    "ann_id": row["ann_id"],
                    "dar_locus": gene,
                    "dar_structure_gene": structure_gene,
                    "dar_genome": row["genome"],
                    "tier": row["tier"],
                    "classification_subset": row["classification_subset"],
                    "exception_applied": row["exception_applied"],
                    "is_cdav_amino_acid": row["is_cdav_amino_acid"],
                    "is_cdav_nucleotide": row["is_cdav_nucleotide"],
                    "n_species_with_disease_allele": row["n_species_with_disease_allele"],
                    "dar_aa_coord": aa_coord,
                    "dar_ref_aa": ref_aa,
                    "dar_alt_aa": row["alt_aa"],
                    "pdb_id": pdb_id,
                    "model_role": model_meta["role"],
                    "state_label": model_meta["state_label"],
                    "dar_chain": chain_id,
                    "dar_struct_res": struct_resnum,
                    "contact_chain": n_chain_id,
                    "contact_gene": n_gene,
                    "contact_structure_gene": n_structure_gene,
                    "contact_resnum": n_res.get_id()[1],
                    "contact_refseq_pos": contact_refseq_pos,
                    "contact_aa": protein_letters_3to1.get(
                        n_res.get_resname().strip().capitalize(), "?"
                    ),
                    "contact_type": contact_type,
                    "contact_class": classify_contact(target_res, n_res),
                    "contact_source": "canonical" if model_meta["role"] == "primary" else "validation",
                })

    # ── Write outputs ──────────────────────────────────────────────────────────
    map_fields = [
        "variant_id",
        "ann_id",
        "interpreted_gene",
        "structure_gene",
        "complex_id",
        "genome",
        "tier",
        "classification_status",
        "classification_subset",
        "exception_applied",
        "classification_warning_reason",
        "classification_coordinate_method",
        "classification_coordinate_status",
        "clinvar_transcript_id",
        "preferred_nm",
        "toga_enst",
        "transcript_map_type",
        "transcript_map_identity",
        "transcript_map_coverage",
        "transcript_reconciliation_status",
        "transcript_coord_delta",
        "anchor_policy_label",
        "anchor_allow_extended",
        "anchor_max_offset",
        "structure_mapping_eligible",
        "aa_coord_raw",
        "is_cdav_amino_acid",
        "is_cdav_nucleotide",
        "aa_coord_method",
        "aa_coord_status",
        "aa_coord",
        "ref_aa",
        "alt_aa",
        "nc_change",
        "aa_change",
        "pdb_id",
        "model_role",
        "state_label",
        "chain_id",
        "struct_resnum",
        "proxy_gene",
        "anchor_method",
        "anchor_failure_detail",
        "status_category",
        "status",
    ]
    with open(OUT_DIR / "dar_structure_map.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=map_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(map_rows)
    print(
        f"\nPosition map   : {OUT_DIR/'dar_structure_map.csv'}  ({len(map_rows)} rows)"
    )

    contact_fields = [
        "variant_id",
        "ann_id",
        "dar_locus",
        "dar_structure_gene",
        "dar_genome",
        "tier",
        "classification_subset",
        "exception_applied",
        "is_cdav_amino_acid",
        "is_cdav_nucleotide",
        "n_species_with_disease_allele",
        "dar_aa_coord",
        "dar_ref_aa",
        "dar_alt_aa",
        "pdb_id",
        "model_role",
        "contact_source",
        "state_label",
        "dar_chain",
        "dar_struct_res",
        "contact_chain",
        "contact_gene",
        "contact_structure_gene",
        "contact_resnum",
        "contact_refseq_pos",
        "contact_aa",
        "contact_type",
        "contact_class",
    ]

    # ── Split canonical vs validation contacts ────────────────────────────────
    # Canonical: contacts from primary-role models (one per complex, highest
    #   quality).  Validation: contacts from validation-role models, deduplicated
    #   to one row per (variant_id, contact_gene, contact_resnum) — the first
    #   occurrence wins (models are already in ascending priority order).
    canonical_contacts = [r for r in contact_rows if r["contact_source"] == "canonical"]
    validation_contacts_raw = [r for r in contact_rows if r["contact_source"] == "validation"]
    seen_val: set[tuple] = set()
    validation_contacts: list[dict] = []
    for r in validation_contacts_raw:
        key = (r["variant_id"], r["contact_gene"], r["contact_resnum"])
        if key not in seen_val:
            seen_val.add(key)
            validation_contacts.append(r)

    with open(OUT_DIR / "dar_contacts_canonical.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=contact_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(canonical_contacts)
    print(
        f"Canonical contacts : {OUT_DIR/'dar_contacts_canonical.csv'}  "
        f"({len(canonical_contacts)} rows)"
    )

    with open(OUT_DIR / "dar_contacts_validation.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=contact_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(validation_contacts)
    print(
        f"Validation contacts: {OUT_DIR/'dar_contacts_validation.csv'}  "
        f"({len(validation_contacts)} rows, deduplicated)"
    )

    # dar_contacts_cbcb8A.csv — full union (canonical + deduplicated validation)
    # for backward-compatibility with downstream scripts
    all_contacts = canonical_contacts + validation_contacts
    with open(OUT_DIR / "dar_contacts_cbcb8A.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=contact_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(all_contacts)
    print(
        f"All contacts (union): {OUT_DIR/'dar_contacts_cbcb8A.csv'}  "
        f"({len(all_contacts)} rows)"
    )

    mito_nuc = [r for r in canonical_contacts if r["contact_type"] == "mt-nuc"]
    with open(OUT_DIR / "dar_mito_nuc_contacts.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=contact_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(mito_nuc)
    print(
        f"Mito↔nuc (canonical): {OUT_DIR/'dar_mito_nuc_contacts.csv'}  ({len(mito_nuc)} rows)"
    )

    summary_rows = []
    by_variant: dict[str, list] = {}
    for row in map_rows:
        by_variant.setdefault(row["variant_id"], []).append(row)
    for variant_id, rows in by_variant.items():
        attempted = [r for r in rows if r.get("pdb_id")]
        mapped = [
            r for r in rows
            if (
                str(r.get("status", "")).startswith("ok")
                or str(r.get("status", "")).startswith("isoform_offset")
                or str(r.get("status", "")).startswith("extended_offset_rescue")
            )
        ]
        primary_attempted = [r for r in attempted if r.get("model_role") == "primary"]
        validation_attempted = [r for r in attempted if r.get("model_role") == "validation"]
        reference_attempted = [r for r in attempted if r.get("model_role") == "reference"]
        primary_mapped = [r for r in mapped if r.get("model_role") == "primary"]
        validation_mapped = [r for r in mapped if r.get("model_role") == "validation"]
        reference_mapped = [r for r in mapped if r.get("model_role") == "reference"]
        mapping_support_class = (
            "primary_and_validation" if primary_mapped and validation_mapped else
            "primary_only" if primary_mapped else
            "validation_only" if validation_mapped else
            "unmapped"
        )
        summary_rows.append({
            "variant_id": variant_id,
            "interpreted_gene": rows[0].get("interpreted_gene"),
            "structure_gene": rows[0].get("structure_gene"),
            "complex_id": rows[0].get("complex_id"),
            "genome": rows[0].get("genome"),
            "classification_subset": rows[0].get("classification_subset"),
            "is_cdav_amino_acid": rows[0].get("is_cdav_amino_acid"),
            "is_cdav_nucleotide": rows[0].get("is_cdav_nucleotide"),
            "clinvar_transcript_id": rows[0].get("clinvar_transcript_id"),
            "preferred_nm": rows[0].get("preferred_nm"),
            "toga_enst": rows[0].get("toga_enst"),
            "transcript_reconciliation_status": rows[0].get("transcript_reconciliation_status"),
            "transcript_coord_delta": rows[0].get("transcript_coord_delta"),
            "anchor_policy_label": rows[0].get("anchor_policy_label"),
            "anchor_allow_extended": rows[0].get("anchor_allow_extended"),
            "anchor_max_offset": rows[0].get("anchor_max_offset"),
            "aa_coord_method": rows[0].get("aa_coord_method"),
            "aa_coord_status": rows[0].get("aa_coord_status"),
            "anchor_methods": "|".join(sorted({str(r["anchor_method"]) for r in rows if r.get("anchor_method")})),
            "anchor_failure_details": "|".join(sorted({str(r["anchor_failure_detail"]) for r in rows if r.get("anchor_failure_detail")})),
            "n_models_attempted": len(attempted),
            "n_models_mapped": len(mapped),
            "n_primary_models_attempted": len(primary_attempted),
            "n_validation_models_attempted": len(validation_attempted),
            "n_reference_models_attempted": len(reference_attempted),
            "n_primary_models_mapped": len(primary_mapped),
            "n_validation_models_mapped": len(validation_mapped),
            "n_reference_models_mapped": len(reference_mapped),
            "mapping_support_class": mapping_support_class,
            "extended_offset_rescue_used": any(
                str(r.get("status", "")).startswith("extended_offset_rescue")
                for r in mapped
            ),
            "primary_pdb_ids_mapped": ",".join(sorted({r["pdb_id"] for r in primary_mapped if r.get("pdb_id")})),
            "validation_pdb_ids_mapped": ",".join(sorted({r["pdb_id"] for r in validation_mapped if r.get("pdb_id")})),
            "reference_pdb_ids_mapped": ",".join(sorted({r["pdb_id"] for r in reference_mapped if r.get("pdb_id")})),
            "mapped_pdb_ids": ",".join(sorted({r["pdb_id"] for r in mapped if r.get("pdb_id")})),
            "status_categories": "|".join(sorted({r["status_category"] for r in rows if r.get("status_category")})),
            "mapping_statuses": "|".join(sorted({r["status"] for r in rows if r.get("status")})),
            "proxy_mapping_used": any(bool(r.get("proxy_gene")) for r in mapped),
        })
    summary_fields = [
        "variant_id", "interpreted_gene", "structure_gene", "complex_id", "genome",
        "classification_subset", "is_cdav_amino_acid", "is_cdav_nucleotide",
        "clinvar_transcript_id", "preferred_nm", "toga_enst",
        "transcript_reconciliation_status", "transcript_coord_delta",
        "anchor_policy_label", "anchor_allow_extended", "anchor_max_offset",
        "aa_coord_method", "aa_coord_status", "anchor_methods", "anchor_failure_details",
        "n_models_attempted", "n_models_mapped",
        "n_primary_models_attempted", "n_validation_models_attempted", "n_reference_models_attempted",
        "n_primary_models_mapped", "n_validation_models_mapped", "n_reference_models_mapped",
        "mapping_support_class", "extended_offset_rescue_used",
        "primary_pdb_ids_mapped", "validation_pdb_ids_mapped", "reference_pdb_ids_mapped", "mapped_pdb_ids",
        "status_categories", "mapping_statuses", "proxy_mapping_used",
    ]
    with open(OUT_DIR / "structure_model_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary_fields)
        w.writeheader()
        w.writerows(summary_rows)

    transcript_audit_rows = []
    transcript_audit_counter = Counter(
        (
            r.get("genome", ""),
            r.get("interpreted_gene", ""),
            r.get("clinvar_transcript_id", ""),
            r.get("preferred_nm", ""),
            r.get("toga_enst", ""),
            r.get("transcript_map_type", ""),
            r.get("transcript_reconciliation_status", ""),
            r.get("aa_coord_method", ""),
            r.get("aa_coord_status", ""),
            r.get("anchor_method", ""),
            r.get("anchor_failure_detail", ""),
            r.get("status_category", ""),
            r.get("status", ""),
        )
        for r in map_rows
        if r.get("genome") == "nucDNA"
    )
    for key, count in sorted(
        transcript_audit_counter.items(),
        key=lambda item: (-item[1], tuple("" if value is None else str(value) for value in item[0])),
    ):
        (
            genome,
            gene,
            clinvar_tx,
            preferred_nm,
            toga_enst,
            transcript_map_type,
            transcript_reconciliation_status,
            aa_coord_method,
            aa_coord_status,
            anchor_method,
            anchor_failure_detail,
            status_category_value,
            status,
        ) = key
        transcript_audit_rows.append({
            "genome": genome,
            "interpreted_gene": gene,
            "clinvar_transcript_id": clinvar_tx,
            "preferred_nm": preferred_nm,
            "toga_enst": toga_enst,
            "transcript_map_type": transcript_map_type,
            "transcript_reconciliation_status": transcript_reconciliation_status,
            "aa_coord_method": aa_coord_method,
            "aa_coord_status": aa_coord_status,
            "anchor_method": anchor_method,
            "anchor_failure_detail": anchor_failure_detail,
            "status_category": status_category_value,
            "status": status,
            "n_rows": count,
        })
    transcript_audit_fields = [
        "genome",
        "interpreted_gene",
        "clinvar_transcript_id",
        "preferred_nm",
        "toga_enst",
        "transcript_map_type",
        "transcript_reconciliation_status",
        "aa_coord_method",
        "aa_coord_status",
        "anchor_method",
        "anchor_failure_detail",
        "status_category",
        "status",
        "n_rows",
    ]
    with open(OUT_DIR / "structure_transcript_reconciliation_audit.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=transcript_audit_fields)
        w.writeheader()
        w.writerows(transcript_audit_rows)

    audit_rows = []
    failure_rows = [
        r for r in map_rows
        if r.get("status_category") not in {
            "mapped_direct",
            "mapped_with_isoform_offset",
            "mapped_with_extended_offset_rescue",
            "secondary_not_attempted",
            "structural_ineligible",
        }
    ]
    audit_counter = Counter(
        (
            r.get("status_category", ""),
            r.get("status", ""),
            r.get("genome", ""),
            r.get("complex_id", ""),
            r.get("pdb_id", ""),
            r.get("interpreted_gene", ""),
            r.get("aa_coord_method", ""),
            r.get("aa_coord_status", ""),
            r.get("anchor_method", ""),
            r.get("anchor_failure_detail", ""),
            r.get("proxy_gene", ""),
        )
        for r in failure_rows
    )
    for key, count in sorted(
        audit_counter.items(),
        key=lambda item: (-item[1], tuple("" if value is None else str(value) for value in item[0])),
    ):
        status_cat, status, genome, complex_id, pdb_id, gene, aa_method, aa_status, anchor_method, anchor_failure_detail, proxy_gene = key
        audit_rows.append({
            "status_category": status_cat,
            "status": status,
            "genome": genome,
            "complex_id": complex_id,
            "pdb_id": pdb_id,
            "interpreted_gene": gene,
            "aa_coord_method": aa_method,
            "aa_coord_status": aa_status,
            "anchor_method": anchor_method,
            "anchor_failure_detail": anchor_failure_detail,
            "proxy_gene": proxy_gene,
            "n_rows": count,
        })
    audit_fields = [
        "status_category",
        "status",
        "genome",
        "complex_id",
        "pdb_id",
        "interpreted_gene",
        "aa_coord_method",
        "aa_coord_status",
        "anchor_method",
        "anchor_failure_detail",
        "proxy_gene",
        "n_rows",
    ]
    with open(OUT_DIR / "structure_mapping_failure_audit.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=audit_fields)
        w.writeheader()
        w.writerows(audit_rows)

    # ── Summary ──────────��─────────────────────────────────────────────────────
    ok_maps = [
        r
        for r in map_rows
        if (
            r["status"] == "ok"
            or r["status"].startswith("ok(")
            or r["status"].startswith("isoform_offset")
            or r["status"].startswith("extended_offset_rescue")
        )
    ]
    isoform_c = sum(1 for r in ok_maps if r["status"].startswith("isoform"))
    status_ctr = Counter(r["status"] for r in map_rows)

    _MAPPED_CATS = {
        "mapped_direct",
        "mapped_with_isoform_offset",
        "mapped_with_extended_offset_rescue",
    }
    n_skipped_policy = sum(
        1 for d in structural_rows
        if d.get("classification_status") == "skipped_by_policy"
    )
    n_eligible = sum(1 for d in structural_rows if d.get("structure_mapping_eligible"))
    mapped_variant_ids = {
        r["variant_id"]
        for r in map_rows
        if r.get("status_category") in _MAPPED_CATS
    }
    n_mapped = len(mapped_variant_ids)
    n_eligible_failed = n_eligible - n_mapped
    rate_pct = (n_mapped / n_eligible * 100) if n_eligible else 0.0

    print(f"\n{'='*50}")
    print("STRUCTURAL MAPPING SUMMARY")
    print(f"{'='*50}")
    print(f"Classified rows loaded      : {len(structural_rows)}")
    print(f"  skipped_by_policy         : {n_skipped_policy}  (benign/VUS, not attempted)")
    print(f"  mapping-eligible          : {n_eligible}")
    print(f"  mapped (unique variants)  : {n_mapped}"
          f"  (incl. {isoform_c} with isoform offset correction)")
    print(f"  eligible but failed       : {n_eligible_failed}")
    print(f"  True mapping rate         : {rate_pct:.1f}%  (mapped / eligible)")
    print(f"Variant summaries           : {len(summary_rows)}")
    print(f"\nStatus breakdown (model rows):")
    for status, n in sorted(status_ctr.items(), key=lambda x: -x[1]):
        print(f"  {status:<40}: {n}")

    cc = Counter(r["contact_class"] for r in all_contacts)
    print(f"\nContact class breakdown:")
    for k, v in sorted(cc.items(), key=lambda x: -x[1]):
        print(f"  {k:<15}: {v:>7,}")

    ct = Counter(r["contact_type"] for r in all_contacts)
    print(f"\n{'Contact type':<10}  {'contacts':>10}  {'unique DARs':>12}")
    for ctype in ("mt-mt", "nuc-nuc", "mt-nuc"):
        n_dars = len({r["variant_id"] for r in all_contacts if r["contact_type"] == ctype})
        print(f"  {ctype:<10}  {ct[ctype]:>10,}  {n_dars:>12,}")

    iface_dars = {r["variant_id"] for r in mito_nuc}
    cdav_iface = len({
        r["variant_id"]
        for r in map_rows
        if r["variant_id"] in iface_dars
        and (
            r["status"] == "ok"
            or r.get("status", "").startswith("ok(")
            or r.get("status", "").startswith("isoform")
            or r.get("status", "").startswith("extended_offset_rescue")
        )
        and r.get("is_cdav_amino_acid")
    })
    print(
        f"\nDARs at mito↔nuc interface : {len(iface_dars):,} total, "
        f"{cdav_iface} are AA-level cDAVs"
    )


if __name__ == "__main__":
    main()
