#!/usr/bin/env python3
"""
Record deterministic alignment-sanitation metadata without recomputation.

This step preserves existing expensive alignment products and only writes
sanitation and normalization metadata required downstream.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

DATA_DIR = ROOT / "data"
DERIVED_CURATED_DIR = DATA_DIR / "derived" / "curated"
DERIVED_REFERENCE_DIR = DATA_DIR / "derived" / "reference"
TOGA_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_aa"
TOGA_CODON_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"
MT_AA_DIR = DATA_DIR / "alignments" / "mtdna_aa"
MT_CODON_DIR = DATA_DIR / "alignments" / "mtdna_codon"
ALIGNMENT_SANITATION_PARQUET = DERIVED_CURATED_DIR / "alignment_sanitation_manifest.parquet"
MT_ACCESSION_NORMALIZATION_PARQUET = DERIVED_CURATED_DIR / "mt_accession_header_normalization.parquet"


def ensure_layout() -> None:
    DERIVED_CURATED_DIR.mkdir(parents=True, exist_ok=True)


def repo_rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def classify_header(header: str) -> tuple[str, str, str, str]:
    content = header[1:].strip()
    left, _, gene = content.partition(" | ")
    fields = left.split("|")
    species = fields[0] if fields else ""
    taxid = fields[1] if len(fields) > 1 else ""
    accession = fields[2] if len(fields) > 2 else ""
    if len(fields) >= 4:
        schema = "species_taxid_assembly_transcript"
    elif accession.startswith("NC_"):
        schema = "species_taxid_accession"
    else:
        schema = "unclassified"
    return schema, species, taxid, accession


def scan_alignment_dir(directory: Path, alignment_type: str) -> tuple[list[dict], list[dict]]:
    manifest_rows: list[dict] = []
    accession_rows: list[dict] = []
    for fasta in sorted(directory.glob("*.fasta")):
        with open(fasta, encoding="utf-8") as handle:
            n_headers = 0
            unresolved = 0
            schemas: set[str] = set()
            while True:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline().strip()
                if not header.startswith(">"):
                    continue
                n_headers += 1
                schema, species, taxid, accession = classify_header(header)
                schemas.add(schema)
                if "UNKNOWN" in header or not taxid:
                    unresolved += 1
                if accession.startswith("NC_"):
                    accession_rows.append({
                        "alignment_file": fasta.name,
                        "alignment_type": alignment_type,
                        "accession": accession,
                        "species": species,
                        "taxid": taxid,
                        "header_schema": schema,
                    })
            manifest_rows.append({
                "alignment_file": fasta.name,
                "alignment_type": alignment_type,
                "input_path": repo_rel(fasta),
                "header_schema": ",".join(sorted(schemas)),
                "n_headers": n_headers,
                "n_unresolved_headers": unresolved,
                "sanitation_rule_version": "phase2_manifest_only_v1",
                "dropped_records": 0,
                "status": "scanned",
            })
    return manifest_rows, accession_rows


def main() -> None:
    import pandas as pd

    ensure_layout()
    manifest_rows: list[dict] = []
    accession_rows: list[dict] = []
    for directory, label in [
        (TOGA_AA_DIR, "toga_aa"),
        (TOGA_CODON_DIR, "toga_codon"),
        (MT_AA_DIR, "mtdna_aa"),
        (MT_CODON_DIR, "mtdna_codon"),
    ]:
        rows, accessions = scan_alignment_dir(directory, label)
        manifest_rows.extend(rows)
        accession_rows.extend(accessions)

    pd.DataFrame(manifest_rows).to_parquet(ALIGNMENT_SANITATION_PARQUET, index=False)
    pd.DataFrame(accession_rows).to_parquet(MT_ACCESSION_NORMALIZATION_PARQUET, index=False)
    print(f"Wrote alignment sanitation manifest: {ALIGNMENT_SANITATION_PARQUET}")


if __name__ == "__main__":
    main()
