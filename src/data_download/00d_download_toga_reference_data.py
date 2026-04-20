#!/usr/bin/env python3
"""
Inspect TOGA reference support for nuclear OXPHOS genes.

This step inventories local TOGA alignment support together with MANE and TOGA
overview references needed by downstream comparative analyses.
"""

from pathlib import Path
import csv
import gzip
import sys
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from utils.gene_reference import GeneReference
from utils.download_paths import LEGACY_REFERENCE_DIR, RAW_REFERENCE_DIR, ensure_layout, latest_existing

DATA_DIR = LEGACY_REFERENCE_DIR.parents[0]
CODON_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"


def fetch_assembly_map() -> dict:
    overview = latest_existing([(RAW_REFERENCE_DIR, "TOGA_overview_table_hg38_*.tsv")])
    if overview is None:
        return {}
    mapping = {}
    with open(overview, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            assembly = row.get("Assembly name", "").strip()
            species = row.get("Species", "").strip().replace(" ", "_")
            if assembly and species:
                mapping[assembly] = species
    return mapping


def fetch_mane_mapping() -> dict:
    mane = latest_existing([(RAW_REFERENCE_DIR, "MANE_GRCh38_v1.5_*.txt.gz"), (RAW_REFERENCE_DIR, "MANE_GRCh38_v1.5.txt.gz")])
    if mane is None:
        return {}
    mapping = {}
    with gzip.open(mane, "rt", encoding="utf-8") as handle:
        for line in handle:
            cols = line.lstrip("#").rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            nm_versioned = cols[5].strip()
            enst_versioned = cols[7].strip()
            if nm_versioned.startswith("NM_") and enst_versioned.startswith("ENST"):
                mapping[nm_versioned.split(".")[0]] = enst_versioned.split(".")[0]
    return mapping


def local_toga_alignment_index() -> dict[str, list[str]]:
    index: dict[str, list[str]] = defaultdict(list)
    if not CODON_DIR.exists():
        return {}
    for fasta in sorted(CODON_DIR.glob("*_codon_alignment.fasta")):
        gene = fasta.name.removesuffix("_codon_alignment.fasta")
        index[gene].append(fasta.name)
    return dict(index)


def mane_gene_symbols(mane_map: dict[str, str]) -> set[str]:
    mane = latest_existing([(RAW_REFERENCE_DIR, "MANE_GRCh38_v1.5_*.txt.gz"), (RAW_REFERENCE_DIR, "MANE_GRCh38_v1.5.txt.gz")])
    if mane is None:
        return set()
    symbols: set[str] = set()
    with gzip.open(mane, "rt", encoding="utf-8") as handle:
        for line in handle:
            cols = line.lstrip("#").rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            symbol = cols[3].strip()
            nm = cols[5].strip().split(".")[0]
            if symbol and nm and nm in mane_map:
                symbols.add(symbol)
    return symbols


def main() -> None:
    ensure_layout()
    hgnc_path = latest_existing([(RAW_REFERENCE_DIR, "Canonical_OXPHOS_Subunits_HGNC_*.csv")])
    if hgnc_path is None:
        raise FileNotFoundError("HGNC gene list not found in data/raw/reference")
    hgnc_ref = GeneReference(hgnc_path)
    assembly_map = fetch_assembly_map()
    mane_map = fetch_mane_mapping()
    toga_index = local_toga_alignment_index()

    target_genes = sorted({
        data["symbol"] for data in hgnc_ref.lookup.values()
        if data.get("symbol") and not data["symbol"].startswith("MT-")
    })
    mane_gene_count = len(set(target_genes) & set(mane_gene_symbols(mane_map)))
    indexed_target_genes = [gene for gene in target_genes if gene in toga_index]

    print(f"TOGA genes indexed: {len(toga_index)}")
    print(f"Target genes: {len(target_genes)}")
    print(f"Target genes present in TOGA index: {len(indexed_target_genes)}")
    print(f"Target genes with MANE transcript support: {mane_gene_count}")
    print(f"Assembly map entries: {len(assembly_map)}")
    print("Existing TOGA alignments are preserved in place; this step inventories reference support only.")


if __name__ == "__main__":
    main()
