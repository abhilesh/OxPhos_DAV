#!/usr/bin/env python3
"""
Build canonical TaxID and species overlap reference products.

This step keeps accession/header normalization separate from comparative
species-overlap bookkeeping.
"""

from pathlib import Path
import csv
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from utils.download_paths import DERIVED_REFERENCE_DIR, LEGACY_REFERENCE_DIR, ensure_layout, sync_compat_file

DATA_DIR = LEGACY_REFERENCE_DIR.parents[0]
TOGA_DIRS = [DATA_DIR / "alignments" / "toga_hg38_aa"]
MT_DIRS = [DATA_DIR / "alignments" / "mtdna_aa"]


def extract_taxid_species(aln_dirs: list[Path]) -> dict[str, str]:
    taxid_to_species: dict[str, str] = {}
    for aln_dir in aln_dirs:
        if not aln_dir.exists():
            continue
        for fasta in sorted(aln_dir.glob("*_aa_alignment.fasta")):
            with open(fasta, encoding="utf-8") as handle:
                for line in handle:
                    if not line.startswith(">"):
                        continue
                    rec_id = line[1:].strip().split()[0]
                    if rec_id.startswith("Homo_sapiens"):
                        continue
                    fields = rec_id.split("|")
                    if len(fields) < 2:
                        continue
                    species = fields[0].strip()
                    taxid = fields[1].strip()
                    if taxid.isdigit() and taxid not in taxid_to_species:
                        taxid_to_species[taxid] = species
    return taxid_to_species


def main() -> None:
    ensure_layout()
    toga_map = extract_taxid_species(TOGA_DIRS)
    mt_map = extract_taxid_species(MT_DIRS)
    both = set(toga_map) & set(mt_map)
    toga_only = set(toga_map) - set(mt_map)
    mt_only = set(mt_map) - set(toga_map)

    out = DERIVED_REFERENCE_DIR / "taxid_species_mapping.csv"
    rows = []
    for taxid in sorted(both, key=int):
        rows.append({
            "taxid": taxid,
            "toga_species": toga_map[taxid],
            "toga_taxid": taxid,
            "mt_species": mt_map[taxid],
            "mt_taxid": taxid,
            "match_type": "Exact_TaxID_Match",
        })
    for taxid in sorted(toga_only, key=int):
        rows.append({
            "taxid": taxid,
            "toga_species": toga_map[taxid],
            "toga_taxid": taxid,
            "mt_species": "",
            "mt_taxid": "",
            "match_type": "TOGA_Only",
        })
    for taxid in sorted(mt_only, key=int):
        rows.append({
            "taxid": taxid,
            "toga_species": "",
            "toga_taxid": "",
            "mt_species": mt_map[taxid],
            "mt_taxid": taxid,
            "match_type": "mtDNA_Only",
        })
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["taxid", "toga_species", "toga_taxid", "mt_species", "mt_taxid", "match_type"],
        )
        writer.writeheader()
        writer.writerows(rows)
    sync_compat_file(out, LEGACY_REFERENCE_DIR / "taxid_species_mapping.csv")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
