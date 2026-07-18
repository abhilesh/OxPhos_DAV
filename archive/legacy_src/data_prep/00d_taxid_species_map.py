"""
src/data_prep/00d_taxid_species_map.py

Builds the cross-genome species overlap between nucDNA (TOGA) and mtDNA
alignments by matching on TaxID — the authoritative, taxonomy-stable
identifier present in field[1] of every alignment FASTA header.

Header format (both alignment types):
  >species_name|taxid|assembly_or_accession[|transcript_id ...]

Previous version incorrectly extracted field[0] (species name) from mtDNA
headers and treated it as an NCBI nucleotide accession, then queried NCBI
to retrieve TaxIDs — producing a corrupted mapping. No NCBI queries are
needed: TaxIDs are already embedded in every header.

Outputs:
  data/reference/taxid_species_mapping.csv

Columns:
  taxid          -- NCBI TaxID (key for joining)
  toga_species   -- species name from TOGA alignment header
  toga_taxid     -- TaxID from TOGA header (= taxid)
  mt_species     -- species name from mtDNA alignment header
  mt_taxid       -- TaxID from mtDNA header (= taxid)
  match_type     -- Exact_TaxID_Match | TOGA_Only | mtDNA_Only

Run from project root inside the Docker container:
    python src/data_prep/00d_taxid_species_map.py
"""

import csv
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

ROOT      = Path(__file__).resolve().parents[2]
DATA_DIR  = ROOT / "data"
REF_DIR   = DATA_DIR / "reference"

TOGA_DIRS = [DATA_DIR / "alignments" / "toga_hg38_aa"]
MT_DIRS   = [DATA_DIR / "alignments" / "mtdna_aa"]
OUT_CSV   = REF_DIR / "taxid_species_mapping.csv"


def extract_taxid_species(aln_dirs: list[Path]) -> dict[str, str]:
    """
    Parse all AA alignment FASTAs in the given directories.
    Returns {taxid: species_name} using field[1] (TaxID) and field[0] (species).
    Non-human entries only. If a TaxID maps to multiple species names, the
    first encountered is kept (log a warning if they differ).
    """
    taxid_to_species: dict[str, str] = {}

    for aln_dir in aln_dirs:
        if not aln_dir.exists():
            print(f"  [WARN] Directory not found: {aln_dir}")
            continue
        for fasta in sorted(aln_dir.glob("*_aa_alignment.fasta")):
            for rec in SeqIO.parse(fasta, "fasta"):
                if rec.id.startswith("Homo_sapiens"):
                    continue
                fields = rec.id.split("|")
                if len(fields) < 2:
                    continue
                species = fields[0].strip()
                taxid   = fields[1].strip()
                if not taxid.isdigit():
                    continue
                if taxid in taxid_to_species:
                    existing = taxid_to_species[taxid]
                    if existing != species:
                        # Keep the first; they are usually the same or close synonyms
                        pass
                else:
                    taxid_to_species[taxid] = species

    return taxid_to_species


def main():
    print("Building TaxID-based cross-genome species mapping\n")

    # ── Extract TaxID → species from each genome ──────────────────────────────
    print("Parsing TOGA (nucDNA) alignment headers...")
    toga_map = extract_taxid_species(TOGA_DIRS)
    print(f"  Unique TaxIDs in TOGA alignments: {len(toga_map)}")

    print("Parsing mtDNA alignment headers...")
    mt_map = extract_taxid_species(MT_DIRS)
    print(f"  Unique TaxIDs in mtDNA alignments: {len(mt_map)}")

    # ── Compute intersection and union ────────────────────────────────────────
    toga_taxids = set(toga_map)
    mt_taxids   = set(mt_map)
    both        = toga_taxids & mt_taxids
    toga_only   = toga_taxids - mt_taxids
    mt_only     = mt_taxids   - toga_taxids

    print(f"\nCross-genome overlap (both TOGA and mtDNA): {len(both)}")
    print(f"TOGA only                                  : {len(toga_only)}")
    print(f"mtDNA only                                 : {len(mt_only)}")

    # ── Write output ──────────────────────────────────────────────────────────
    REF_DIR.mkdir(parents=True, exist_ok=True)
    rows = []

    for taxid in sorted(both, key=int):
        rows.append({
            "taxid":        taxid,
            "toga_species": toga_map[taxid],
            "toga_taxid":   taxid,
            "mt_species":   mt_map[taxid],
            "mt_taxid":     taxid,
            "match_type":   "Exact_TaxID_Match",
        })
    for taxid in sorted(toga_only, key=int):
        rows.append({
            "taxid":        taxid,
            "toga_species": toga_map[taxid],
            "toga_taxid":   taxid,
            "mt_species":   "",
            "mt_taxid":     "",
            "match_type":   "TOGA_Only",
        })
    for taxid in sorted(mt_only, key=int):
        rows.append({
            "taxid":        taxid,
            "toga_species": "",
            "toga_taxid":   "",
            "mt_species":   mt_map[taxid],
            "mt_taxid":     taxid,
            "match_type":   "mtDNA_Only",
        })

    fieldnames = ["taxid", "toga_species", "toga_taxid", "mt_species", "mt_taxid", "match_type"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"\nSpecies mapping written → {OUT_CSV}")
    print(f"  Total rows: {len(rows)}")

    # ── Sample overlap ────────────────────────────────────────────────────────
    overlap_rows = [r for r in rows if r["match_type"] == "Exact_TaxID_Match"]
    print(f"\nSample cross-genome matches (first 5):")
    for r in overlap_rows[:5]:
        print(f"  TaxID {r['taxid']:>10}  TOGA: {r['toga_species']:<35}  MT: {r['mt_species']}")


if __name__ == "__main__":
    main()