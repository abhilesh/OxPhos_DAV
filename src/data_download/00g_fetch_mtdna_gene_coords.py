#!/usr/bin/env python3
"""
Build mitochondrial gene coordinate references.

This reference is used for gene assignment and overlap-aware mtDNA curation.
"""

from pathlib import Path
import csv
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from Bio import Entrez, SeqIO

from utils.download_paths import DERIVED_REFERENCE_DIR, LEGACY_REFERENCE_DIR, ensure_layout, sync_compat_file

Entrez.email = "pipeline@analysis.local"


def fetch_mtdna_coords() -> dict[str, tuple[int, int, str]]:
    handle = Entrez.efetch(db="nucleotide", id="NC_012920.1", rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    gene_map: dict[str, tuple[int, int, str]] = {}
    for feature in record.features:
        if feature.type != "CDS":
            continue
        gene_name = feature.qualifiers.get("gene", [None])[0]
        if not gene_name:
            continue
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand = "+" if feature.location.strand == 1 else "-"
        standard_name = f"MT-{gene_name}" if not gene_name.startswith("MT-") else gene_name
        gene_map[standard_name] = (start, end, strand)
    return gene_map


def main() -> None:
    ensure_layout()
    out = DERIVED_REFERENCE_DIR / "mtdna_gene_coordinates.tsv"
    coords = fetch_mtdna_coords()
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "start", "end", "strand"])
        for gene in sorted(coords):
            start, end, strand = coords[gene]
            writer.writerow([gene, start, end, strand])
    sync_compat_file(out, LEGACY_REFERENCE_DIR / "mtdna_gene_coordinates.tsv")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
