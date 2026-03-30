"""
NT-Level Codon Diagnostic — mtDNA and nucDNA
=============================================
For each genome, samples up to N_SAMPLE variants that have an AA-level c-DAR
and performs the Triangle of Truth: inject alt_nt into the WT codon extracted
from the alignment, translate, and verify it matches the expected alt_aa.

Prints a per-variant trace so coordinate or allele mismatches are immediately
visible.

Run from project root:
    python tests/test_diagnose_codons.py
"""

import csv
import json
import re
import sys
from pathlib import Path

from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from utils.alignment_parser import AlignmentParser

# ==== Configuration ====
CURATED_DIR   = ROOT / "data" / "annotations" / "curated"
MT_JSON       = CURATED_DIR / "cdar_classifications_mtDNA.json"
NUC_JSON      = CURATED_DIR / "cdar_classifications_nucDNA.json"

TOGA_AA_DIR   = ROOT / "data" / "alignments" / "toga_hg38_aa"
TOGA_NT_DIR   = ROOT / "data" / "alignments" / "toga_hg38_codon"
MT_AA_DIR     = ROOT / "data" / "alignments" / "mtdna_aa"
MT_NT_DIR     = ROOT / "data" / "alignments" / "mtdna_codon"
MT_COORD_FILE = ROOT / "data" / "reference" / "mtdna_gene_coordinates.tsv"

N_SAMPLE = 20  # variants to diagnose per genome


# ==== Helpers ====

def aa_dir(genome: str) -> Path:
    return MT_AA_DIR if genome == "mtDNA" else TOGA_AA_DIR


def nt_dir(genome: str) -> Path:
    return MT_NT_DIR if genome == "mtDNA" else TOGA_NT_DIR


_MT_ALIAS = {"MT-COX1": "MT-CO1", "MT-COX2": "MT-CO2", "MT-COX3": "MT-CO3", "MT-CYTB": "MT-CYB"}


def _load_mt_coords() -> dict:
    """Returns {gene: (start, end, strand)} from mtdna_gene_coordinates.tsv."""
    coords = {}
    if not MT_COORD_FILE.exists():
        return coords
    with open(MT_COORD_FILE) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            start, end = int(row["start"]), int(row["end"])
            entry = (min(start, end), max(start, end), row["strand"])
            coords[row["gene"]] = entry
            if row["gene"] in _MT_ALIAS:
                coords[_MT_ALIAS[row["gene"]]] = entry
    return coords


MT_COORDS = _load_mt_coords()


def _genomic_to_cds(genomic_pos: int, locus: str) -> int | None:
    """Converts a rCRS genomic position to a 1-indexed CDS position."""
    if locus not in MT_COORDS:
        return None
    start, end, strand = MT_COORDS[locus]
    if strand == "+":
        return genomic_pos - start + 1
    else:
        return end - genomic_pos + 1


def parse_variant_coordinates(variant: dict) -> tuple:
    """Returns (aa_pos, wt_aa, mut_aa, nt_pos).

    For mtDNA: derives the exact CDS position from the genomic position stored
    in the variant, which correctly handles second- and third-base mutations.
    For nucDNA: uses the c. CDS coordinate from nc_change.
    Fallback: (aa_pos * 3) - 2  — always the first base of the codon.
    """
    aa_str = variant.get("aa_change", "")
    nc_str = variant.get("nc_change", "")

    aa_match = re.search(r"([a-zA-Z]+)(\d+)([a-zA-Z]+)", aa_str)
    if not aa_match:
        return None, None, None, None

    wt_aa  = aa_match.group(1)
    aa_pos = int(aa_match.group(2))
    mut_aa = aa_match.group(3)

    # nucDNA: exact CDS position from c. notation
    if variant.get("genome") == "nucDNA":
        nc_match = re.search(r"c\.(\d+)", nc_str)
        if nc_match:
            return aa_pos, wt_aa, mut_aa, int(nc_match.group(1))

    # mtDNA: derive exact CDS position from genomic coordinate
    genomic_pos = variant.get("genomic_pos")
    locus       = variant.get("locus", "").split("/")[0]
    if genomic_pos and locus:
        cds_pos = _genomic_to_cds(genomic_pos, locus)
        if cds_pos:
            return aa_pos, wt_aa, mut_aa, cds_pos

    # Final fallback — first base of the codon only
    return aa_pos, wt_aa, mut_aa, (aa_pos * 3) - 2


def diagnose_genome(json_path: Path, genome: str):
    if not json_path.exists():
        print(f"  {json_path.name} not found — run 00_classify_DAV.py first.\n")
        return

    with open(json_path) as f:
        variants = json.load(f)

    candidates = [v for v in variants if v.get("cdar_aa")]
    print(f"Loaded {len(variants)} variants, {len(candidates)} with AA c-DARs.")

    loaded = {}
    tested = passed = failed = skipped = 0

    for var in candidates:
        if tested >= N_SAMPLE:
            break

        locus = var["locus"].split("/")[0]
        aa_pos, wt_aa, mut_aa, nt_pos = parse_variant_coordinates(var)
        if not aa_pos:
            skipped += 1
            continue

        aa_path = aa_dir(genome) / f"{locus}_aa_alignment.fasta"
        nt_path = nt_dir(genome) / f"{locus}_codon_alignment.fasta"

        if locus not in loaded:
            if aa_path.exists() and nt_path.exists():
                loaded[locus] = AlignmentParser(aa_path, nt_path, genome)
            else:
                loaded[locus] = None

        parser = loaded[locus]
        if not parser:
            skipped += 1
            continue

        # Use find_sequence_anchor to correct for alignment start offsets (e.g. MACSE)
        true_aa_pos = parser.find_sequence_anchor(aa_pos, wt_aa)
        if true_aa_pos is None:
            skipped += 1
            continue

        # Shift nt_pos by the same amount as the AA anchor correction
        aa_shift    = true_aa_pos - aa_pos
        nt_pos_corr = nt_pos + aa_shift * 3

        if nt_pos_corr not in parser.nt_map:
            skipped += 1
            continue

        # Extract WT codon from alignment using anchor-corrected position
        codon_start_bio = nt_pos_corr - ((nt_pos_corr - 1) % 3)
        col_start       = parser.nt_map[codon_start_bio]
        wt_codon        = parser.nt_alignment[parser.ref_header][col_start: col_start + 3]

        alt_nt       = var.get("alt_nt", "MISSING")
        pos_in_codon = (nt_pos_corr - 1) % 3
        mut_codon    = wt_codon[:pos_in_codon] + alt_nt + wt_codon[pos_in_codon + 1:]

        _table = 2 if genome == "mtDNA" else 1
        try:
            translated_wt  = str(Seq(wt_codon).translate(table=_table))
            translated_mut = str(Seq(mut_codon).translate(table=_table))
        except Exception:
            translated_wt = translated_mut = "ERR"

        alt_aa   = var.get("alt_aa") or mut_aa
        is_pass  = translated_mut == alt_aa
        wt_pass  = translated_wt  == wt_aa

        if is_pass:
            passed += 1
        else:
            failed += 1

        tested += 1
        status  = "PASS" if is_pass else "FAIL"
        wt_flag = "" if wt_pass else "  [WT MISMATCH]"

        print(f"{'='*55}")
        print(f"[{status}] {locus}  {var.get('nc_change','')}  ({var.get('aa_change','')})")
        print(f"  alt_nt       : '{alt_nt}'")
        print(f"  WT  codon    : {wt_codon} → {translated_wt}  (expected {wt_aa}){wt_flag}")
        print(f"  Mut codon    : {mut_codon} → {translated_mut}  (expected {alt_aa})")
        if not is_pass:
            print(f"  *** Triangle of Truth FAILED ***")

    print(f"\n{'─'*55}")
    print(f"{genome} — tested={tested}  pass={passed}  fail={failed}  skipped={skipped}")


def main():
    for genome, path in [("mtDNA", MT_JSON), ("nucDNA", NUC_JSON)]:
        print(f"\n{'#'*55}")
        print(f"  {genome} Codon Diagnostic  (up to {N_SAMPLE} c-DARs)")
        print(f"{'#'*55}\n")
        diagnose_genome(path, genome)


if __name__ == "__main__":
    main()
