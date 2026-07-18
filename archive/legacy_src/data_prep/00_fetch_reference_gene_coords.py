import csv
import time
from pathlib import Path
from datetime import date
from Bio import Entrez

from utils.parsers import GeneReference
from utils.utils import get_latest

# ==== Configuration ====
# Get today's date for metadata
today = date.today().isoformat()

Entrez.email = "abhilesh7@gmail.com"

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
REF_DIR = DATA_DIR / "reference"
REF_DIR.mkdir(parents=True, exist_ok=True)

# Path to your specific HGNC file snippet
HGNC_FILE = get_latest(REF_DIR, "Canonical_OXPHOS_Subunits_HGNC*.csv")
OUT_TSV = REF_DIR / f"reference_gene_coordinates{today}.tsv"


def fetch_gene_metadata(symbol):
    """Fetches GRCh38 coordinates and strand for a nuclear gene symbol."""
    try:
        # Search for the gene in the 'gene' database restricted to Human
        search_term = f"({symbol}[Gene Name]) AND (Homo sapiens[Organism])"
        search_handle = Entrez.esearch(db="gene", term=search_term)
        search_results = Entrez.read(search_handle)

        if not search_results["IdList"]:
            return None

        gene_id = search_results["IdList"][0]
        summary_handle = Entrez.esummary(db="gene", id=gene_id)
        summary = Entrez.read(summary_handle)

        # Extract genomic location from the GRCh38 annotation
        genomic_info = summary["DocumentSummarySet"]["DocumentSummary"][0].get(
            "GenomicInfo", []
        )

        for info in genomic_info:
            # Ensure we are pulling coordinates for the primary GRCh38 assembly
            if info.get("ChrLoc"):
                start = int(info["ChrStart"]) + 1  # Convert 0-indexed to 1-indexed
                end = int(info["ChrStop"]) + 1

                return {
                    "gene": symbol,
                    "chr": info["ChrLoc"],
                    "start": start,
                    "end": end,
                    # If start > end, NCBI denotes the minus strand
                    "strand": "-" if start > end else "+",
                }
    except Exception as e:
        print(f"    Error querying {symbol}: {e}")
    return None


def main():
    print(f"Reading HGNC list from {HGNC_FILE.name}...")

    gene_list = []
    with open(HGNC_FILE, "r", encoding="utf-8-sig") as f:
        # Fixed: Explicitly using tab delimiter to match your file structure
        reader = csv.DictReader(f, delimiter="\t")

        if "Approved symbol" not in reader.fieldnames:
            print(
                f"Error: Could not find 'Approved symbol'. Found: {reader.fieldnames}"
            )
            return

        gene_list = [
            row["Approved symbol"].strip() for row in reader if row["Approved symbol"]
        ]

    print(f"Found {len(gene_list)} target genes. Fetching GRCh38 coordinates...")

    results = []
    for i, gene in enumerate(gene_list):
        print(f"  [{i+1}/{len(gene_list)}] Querying {gene}...")
        meta = fetch_gene_metadata(gene)
        if meta:
            results.append(meta)
        else:
            print(f"    Warning: No GRCh38 metadata found for {gene}")

        time.sleep(0.34)  # Respect NCBI API rate limits (3 requests per second)

    # Save to TSV for use in 01_curate_variants.py
    with open(OUT_TSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["gene", "chr", "start", "end", "strand"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSuccessfully saved {len(results)} gene coordinates to {OUT_TSV.name}")


if __name__ == "__main__":
    main()
