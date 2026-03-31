#!/usr/bin/env python3
"""
src/structural/00_map_davs_to_structure.py

Map DAR/c-DAR residues to OXPHOS complex cryo-EM structures and extract
Cβ–Cβ 8 Å contact shells.

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
  data/annotations/curated/cdar_classifications_mtDNA.json
  data/annotations/curated/cdar_classifications_nucDNA.json
  data/structures/{PDB_ID}.cif
  data/Homo_sapiens_OXPHOS_seqs/Homo_sapiens_OXPHOS_mtDNA.fasta
  data/Homo_sapiens_OXPHOS_seqs/Homo_sapiens_OXPHOS_nucDNA.fasta

Output:
  results/structural/dar_structure_map.csv      — position → chain mapping, status per DAR
  results/structural/dar_contacts_cbcb8A.csv    — all Cβ–Cβ 8 Å contacts
  results/structural/dar_mito_nuc_contacts.csv  — cross-genome contacts only

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

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import MMCIFParser, NeighborSearch, is_aa

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT = Path("/app")
STRUC_DIR = ROOT / "data" / "structures"
CURATED_DIR = ROOT / "data" / "annotations" / "curated"
MT_JSON = CURATED_DIR / "cdar_classifications_mtDNA.json"
NUC_JSON = CURATED_DIR / "cdar_classifications_nucDNA.json"
TOGA_AA_DIR = ROOT / "data" / "alignments" / "toga_hg38_aa"
MT_AA_DIR   = ROOT / "data" / "alignments" / "mtdna_aa"
OUT_DIR = ROOT / "results" / "structural"

# ── Isoform proxies for structural mapping ─────────────────────────────────────
# Tissue-specific or isoform-2 subunits absent from the preferred PDB are mapped
# to their isoform-1 structural equivalent. Position numbering is shared (these
# isoforms differ mainly in N-terminal targeting sequences, not core structure).
# Status field will be tagged "proxy=<gene>" so the substitution is traceable.
_ISOFORM_PROXY = {
    "COX4I2":  "COX4I1",   # lung isoform; 9I6F has COX4I1
    "COX6A2":  "COX6A1",   # heart/muscle isoform; 9I6F has COX6A1
    "COX7A1":  "COX7A2",   # heart/muscle isoform; 9I6F has COX7A2
    "ATP5MC2": "ATP5MC1",  # c-subunit isoform 2; 8H9S has ATP5MC1
    "ATP5MC3": "ATP5MC1",  # c-subunit isoform 3; 8H9S has ATP5MC1
}

# ── TOGA filename → canonical HGNC symbol ─────────────────────────────────────
# TOGA alignment files are stored under the old gene name; map them to the
# current approved symbol used everywhere else in the pipeline.
_TOGA_TO_CANONICAL = {
    "COXFA4": "NDUFA4",   # renamed: NDUFA4 is the ClinVar/HGNC canonical name
}

# ── Complex → preferred PDB ────────────────────────────────────────────────────
PREFERRED_PDB = {
    "CI": "9I4I",
    "CII": "8GS8",
    "CIII": "9HZL",
    "CIV": "9I6F",
    "CV": "8H9S",
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


# ── Variant field helpers ──────────────────────────────────────────────────────


def parse_aa_coord(aa_change: str) -> int | None:
    """Extract integer AA position from 'S45F' → 45."""
    m = re.search(r"[A-Za-z]+(\d+)[A-Za-z]+", aa_change or "")
    return int(m.group(1)) if m else None


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

    # Load DARs from c-DAR classification JSONs (carry all fields downstream)
    all_dars = []
    for json_path, genome in [(MT_JSON, "mtDNA"), (NUC_JSON, "nucDNA")]:
        if not json_path.exists():
            print(f"  WARNING: {json_path.name} not found — skipping {genome}.")
            continue
        with open(json_path) as f:
            records = json.load(f)
        for r in records:
            aa_coord = parse_aa_coord(r.get("aa_change", ""))
            if aa_coord is None:
                continue
            all_dars.append(
                {
                    "ann_id": r.get("ann_id", ""),
                    "locus": r.get("locus", ""),
                    "genome": genome,
                    "tier": r.get("tier", ""),
                    "nc_change": r.get("nc_change", ""),
                    "aa_change": r.get("aa_change", ""),
                    "ref_aa": r.get("ref_aa", ""),
                    "alt_aa": r.get("alt_aa", ""),
                    "cdar_aa": r.get("cdar_aa", False),
                    "cdar_nt": r.get("cdar_nt", False),
                    "compensating_species_count": r.get(
                        "compensating_species_count", 0
                    ),
                    "aa_coord": aa_coord,
                }
            )
    print(
        f"Loaded {len(all_dars)} DARs "
        f"({sum(1 for d in all_dars if d['genome']=='mtDNA')} mtDNA, "
        f"{sum(1 for d in all_dars if d['genome']=='nucDNA')} nucDNA)."
    )

    # ── Per-PDB structure cache ────────────────────────────────────────────────
    # Built once per PDB, reused for all variants in that complex.
    struct_cache = {}

    map_rows = []  # one row per DAR — position mapping result + status
    contact_rows = []  # one row per (DAR, contact residue)

    for dar in all_dars:
        gene = dar["locus"].split("/")[0]
        complex_id = GENE_COMPLEX.get(gene)
        if not complex_id:
            continue

        pdb_id = PREFERRED_PDB[complex_id]
        cif_path = STRUC_DIR / f"{pdb_id}.cif"
        if not cif_path.exists():
            continue

        # ── Build structure cache (once per PDB) ──────────────────────────────
        if pdb_id not in struct_cache:
            print(f"\nIndexing {pdb_id} ({complex_id})...")
            model = cif_parser.get_structure(pdb_id, str(cif_path))[0]

            chain_to_gene = {}
            gene_to_chain = {}
            chain_residues = {}
            chain_seqs = {}
            cb_atoms = []

            for chain in model.get_chains():
                res_list = get_chain_residues(chain)
                if not res_list:
                    continue
                chain_residues[chain.id] = res_list
                chain_seqs[chain.id] = get_chain_sequence(res_list)
                cb_atoms.extend(
                    a for r in res_list for a in [get_cb_atom(r)] if a is not None
                )

            # RCSB API → chain→gene assignment (fallback: local alignment) ──────
            api_map = fetch_chain_gene_map(pdb_id, set(ref_seqs.keys()))
            if api_map:
                print(
                    f"  Using RCSB API gene-chain map ({len(api_map)} chains annotated)."
                )
                for chain_id, gene in api_map.items():
                    if chain_id in chain_seqs and gene in ref_seqs:
                        chain_to_gene[chain_id] = gene
                        gene_to_chain[gene] = chain_id
                        print(f"  Chain {chain_id:2} → {gene:10}  (RCSB API)")
            else:
                # Fallback: local alignment (slower, used when offline or PDB not in RCSB)
                print(f"  RCSB API unavailable — falling back to local alignment.")
                for chain_id, chain_seq in chain_seqs.items():
                    best_gene, best_score = None, 0.0
                    for r_gene, r_seq in ref_seqs.items():
                        if GENE_COMPLEX.get(r_gene) != complex_id:
                            continue
                        score = assign_chain_to_gene(r_seq, chain_seq)
                        threshold = min(len(r_seq), len(chain_seq)) * 2 * 0.30
                        if score > best_score and score > threshold:
                            best_score, best_gene = score, r_gene
                    if best_gene:
                        chain_to_gene[chain_id] = best_gene
                        gene_to_chain[best_gene] = chain_id
                        print(
                            f"  Chain {chain_id:2} → {best_gene:10}  (alignment score {best_score:.0f})"
                        )

            # GLOBAL alignment → per-gene position maps ────────────────────────
            pos_maps = {}
            for r_gene, chain_id in gene_to_chain.items():
                r_seq = ref_seqs.get(r_gene, "")
                chain_seq = chain_seqs.get(chain_id, "")
                if r_seq and chain_seq:
                    pos_maps[r_gene] = map_refseq_to_chain(r_seq, chain_seq)

            struct_cache[pdb_id] = {
                "chain_to_gene": chain_to_gene,
                "gene_to_chain": gene_to_chain,
                "pos_maps": pos_maps,
                "chain_residues": chain_residues,
                "ns": NeighborSearch(cb_atoms),
            }

        cache = struct_cache[pdb_id]
        aa_coord = dar["aa_coord"]
        ref_aa = dar["ref_aa"]

        # Try isoform proxy if the gene has no chain in this PDB
        proxy_gene = None
        if gene not in cache["gene_to_chain"]:
            proxy_gene = _ISOFORM_PROXY.get(gene)
            if not (proxy_gene and proxy_gene in cache["gene_to_chain"]):
                map_rows.append(
                    {
                        **dar,
                        "pdb_id": pdb_id,
                        "chain_id": "",
                        "struct_resnum": "",
                        "status": "no_chain_in_pdb",
                    }
                )
                continue

        lookup_gene = proxy_gene if proxy_gene else gene
        chain_id = cache["gene_to_chain"][lookup_gene]
        pos_map = cache["pos_maps"].get(lookup_gene, {})
        ref_seq = ref_seqs.get(lookup_gene, "")

        # ── Anchor: verify ref_aa at aa_coord, slide ±10 for isoform offsets ──
        true_coord = find_anchor(pos_map, aa_coord, ref_aa, ref_seq, window=10)
        if true_coord is None:
            map_rows.append(
                {
                    **dar,
                    "pdb_id": pdb_id,
                    "chain_id": chain_id,
                    "struct_resnum": "",
                    "status": "position_not_anchored",
                }
            )
            continue

        struct_idx = pos_map[true_coord]
        res_list = cache["chain_residues"][chain_id]
        if struct_idx >= len(res_list):
            map_rows.append(
                {
                    **dar,
                    "pdb_id": pdb_id,
                    "chain_id": chain_id,
                    "struct_resnum": "",
                    "status": "index_out_of_range",
                }
            )
            continue

        target_res = res_list[struct_idx]
        struct_resnum = target_res.get_id()[1]
        struct_aa = protein_letters_3to1.get(
            target_res.get_resname().strip().capitalize(), "?"
        )

        # Final residue identity check (should pass after anchoring)
        if struct_aa != ref_aa:
            map_rows.append(
                {
                    **dar,
                    "pdb_id": pdb_id,
                    "chain_id": chain_id,
                    "struct_resnum": struct_resnum,
                    "status": f"aa_mismatch(struct={struct_aa},ref={ref_aa})",
                }
            )
            continue

        offset = true_coord - aa_coord
        status = f"isoform_offset_{offset:+d}" if offset else "ok"
        if proxy_gene:
            status += f"(proxy={proxy_gene})"
        map_rows.append(
            {
                **dar,
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "struct_resnum": struct_resnum,
                "status": status,
            }
        )

        # ── Cβ–Cβ 8 Å neighbor search ─────────────────────────────────────────
        target_cb = get_cb_atom(target_res)
        if target_cb is None:
            continue

        for n_res in cache["ns"].search(target_cb.coord, 8.0, "R"):
            if n_res is target_res or not is_aa(n_res, standard=True):
                continue

            n_chain_id = n_res.get_parent().id
            n_gene = cache["chain_to_gene"].get(n_chain_id, f"Unknown({n_chain_id})")
            contact_is_mt = n_gene in MT_GENES
            dar_is_mt = gene in MT_GENES
            contact_type = (
                "mt-mt"
                if dar_is_mt and contact_is_mt
                else "nuc-nuc" if not dar_is_mt and not contact_is_mt else "mt-nuc"
            )

            # Reverse-map PDB chain index → RefSeq position for the contact residue
            n_res_list = cache["chain_residues"].get(n_chain_id, [])
            try:
                n_chain_idx = n_res_list.index(n_res)
            except ValueError:
                n_chain_idx = -1
            n_pos_map = cache["pos_maps"].get(n_gene, {})
            _rev_map = {v: k for k, v in n_pos_map.items()}
            contact_refseq_pos = _rev_map.get(n_chain_idx, "")

            contact_rows.append(
                {
                    "ann_id": dar["ann_id"],
                    "dar_locus": gene,
                    "dar_genome": dar["genome"],
                    "tier": dar["tier"],
                    "cdar_aa": dar["cdar_aa"],
                    "cdar_nt": dar["cdar_nt"],
                    "compensating_species_count": dar["compensating_species_count"],
                    "dar_aa_coord": aa_coord,
                    "dar_ref_aa": ref_aa,
                    "dar_alt_aa": dar["alt_aa"],
                    "pdb_id": pdb_id,
                    "dar_chain": chain_id,
                    "dar_struct_res": struct_resnum,
                    "contact_chain": n_chain_id,
                    "contact_gene": n_gene,
                    "contact_resnum": n_res.get_id()[1],
                    "contact_refseq_pos": contact_refseq_pos,
                    "contact_aa": protein_letters_3to1.get(
                        n_res.get_resname().strip().capitalize(), "?"
                    ),
                    "contact_type": contact_type,
                    "contact_class": classify_contact(target_res, n_res),
                }
            )

    # ── Write outputs ──────────────────────────────────────────────────────────
    map_fields = [
        "ann_id",
        "locus",
        "genome",
        "tier",
        "cdar_aa",
        "cdar_nt",
        "aa_coord",
        "ref_aa",
        "alt_aa",
        "nc_change",
        "aa_change",
        "pdb_id",
        "chain_id",
        "struct_resnum",
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
        "ann_id",
        "dar_locus",
        "dar_genome",
        "tier",
        "cdar_aa",
        "cdar_nt",
        "compensating_species_count",
        "dar_aa_coord",
        "dar_ref_aa",
        "dar_alt_aa",
        "pdb_id",
        "dar_chain",
        "dar_struct_res",
        "contact_chain",
        "contact_gene",
        "contact_resnum",
        "contact_refseq_pos",
        "contact_aa",
        "contact_type",
        "contact_class",
    ]
    with open(OUT_DIR / "dar_contacts_cbcb8A.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=contact_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(contact_rows)
    print(
        f"All contacts   : {OUT_DIR/'dar_contacts_cbcb8A.csv'}  ({len(contact_rows)} rows)"
    )

    mito_nuc = [r for r in contact_rows if r["contact_type"] == "mt-nuc"]
    with open(OUT_DIR / "dar_mito_nuc_contacts.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=contact_fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(mito_nuc)
    print(
        f"Mito↔nuc       : {OUT_DIR/'dar_mito_nuc_contacts.csv'}  ({len(mito_nuc)} rows)"
    )

    # ── Summary ──────────��─────────────────────────────────────────────────────
    ok_maps = [
        r
        for r in map_rows
        if r["status"] == "ok" or r["status"].startswith("ok(") or r["status"].startswith("isoform_offset")
    ]
    isoform_c = sum(1 for r in ok_maps if r["status"].startswith("isoform"))
    status_ctr = Counter(r["status"] for r in map_rows)

    print(f"\n{'='*50}")
    print("STRUCTURAL MAPPING SUMMARY")
    print(f"{'='*50}")
    print(f"DARs processed          : {len(all_dars)}")
    print(
        f"Mapped successfully     : {len(ok_maps)}"
        f"  (incl. {isoform_c} with isoform offset correction)"
    )
    print(f"\nStatus breakdown:")
    for status, n in sorted(status_ctr.items(), key=lambda x: -x[1]):
        print(f"  {status:<40}: {n}")

    cc = Counter(r["contact_class"] for r in contact_rows)
    print(f"\nContact class breakdown:")
    for k, v in sorted(cc.items(), key=lambda x: -x[1]):
        print(f"  {k:<15}: {v:>7,}")

    ct = Counter(r["contact_type"] for r in contact_rows)
    print(f"\n{'Contact type':<10}  {'contacts':>10}  {'unique DARs':>12}")
    for ctype in ("mt-mt", "nuc-nuc", "mt-nuc"):
        n_dars = len({r["ann_id"] for r in contact_rows if r["contact_type"] == ctype})
        print(f"  {ctype:<10}  {ct[ctype]:>10,}  {n_dars:>12,}")

    iface_dars = {r["ann_id"] for r in mito_nuc}
    cdar_iface = sum(
        1
        for r in map_rows
        if r["ann_id"] in iface_dars
        and (r["status"] == "ok" or r.get("status", "").startswith("ok(") or r.get("status", "").startswith("isoform"))
        and r.get("cdar_aa")
    )
    print(
        f"\nDARs at mito↔nuc interface : {len(iface_dars):,} total, "
        f"{cdar_iface} are AA-level c-DARs"
    )


if __name__ == "__main__":
    main()
