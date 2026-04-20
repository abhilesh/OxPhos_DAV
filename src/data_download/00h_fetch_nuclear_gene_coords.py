#!/usr/bin/env python3
"""
Build nuclear gene coordinate references for OXPHOS targets.

These coordinates support deterministic transcript and genomic mapping layers.
"""

from pathlib import Path
import csv
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from Bio import Entrez

from utils.gene_reference import GeneReference
from utils.download_paths import (
    DERIVED_REFERENCE_DIR,
    LEGACY_REFERENCE_DIR,
    RAW_REFERENCE_DIR,
    ensure_layout,
    latest_existing,
    sync_compat_file,
)

Entrez.email = "pipeline@analysis.local"


def fetch_gene_metadata(symbol: str) -> dict | None:
    try:
        search_term = f"({symbol}[Gene Name]) AND (Homo sapiens[Organism])"
        search_handle = Entrez.esearch(db="gene", term=search_term)
        search_results = Entrez.read(search_handle)
        if not search_results["IdList"]:
            return None
        gene_id = search_results["IdList"][0]
        summary_handle = Entrez.esummary(db="gene", id=gene_id)
        summary = Entrez.read(summary_handle)
        genomic_info = summary["DocumentSummarySet"]["DocumentSummary"][0].get("GenomicInfo", [])
        for info in genomic_info:
            if info.get("ChrLoc"):
                start = int(info["ChrStart"]) + 1
                end = int(info["ChrStop"]) + 1
                return {
                    "gene": symbol,
                    "chr": info["ChrLoc"],
                    "start": start,
                    "end": end,
                    "strand": "-" if start > end else "+",
                }
    except Exception:
        return None
    return None


def main() -> None:
    ensure_layout()
    hgnc_path = latest_existing([(RAW_REFERENCE_DIR, "Canonical_OXPHOS_Subunits_HGNC_*.csv")])
    if hgnc_path is None:
        raise FileNotFoundError("HGNC gene list not found in data/raw/reference")
    hgnc_ref = GeneReference(hgnc_path)
    genes = sorted({
        data["symbol"] for data in hgnc_ref.lookup.values()
        if data.get("symbol") and not data["symbol"].startswith("MT-")
    })

    rows = []
    for gene in genes:
        meta = fetch_gene_metadata(gene)
        if meta:
            rows.append(meta)
        time.sleep(0.34)

    out = DERIVED_REFERENCE_DIR / "nucdna_gene_coordinates.tsv"
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "chr", "start", "end", "strand"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    sync_compat_file(out, LEGACY_REFERENCE_DIR / "nucdna_gene_coordinates.tsv")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
