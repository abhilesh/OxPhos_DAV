import ssl
import time
import urllib.parse
import urllib.request
import re
from pathlib import Path
import xml.etree.ElementTree as ET
from collections import defaultdict

# Import utility classes
from utils.parsers import GeneReference
from utils.utils import get_latest

# ==== Configuration ====
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
RAW_MT_DIR = DATA_DIR / "alignments" / "mtdna_raw_cds"
RAW_MT_DIR.mkdir(parents=True, exist_ok=True)

# NCBI Entrez Configuration
ENTREZ_EMAIL = "abhilesh7@gmail.com"  # Replace with your email
BATCH_SIZE = 500
RATE_LIMIT = (
    0.5  # Seconds between requests (NCBI requires <= 3 requests/sec without API key)
)

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


def search_ncbi_mammals() -> tuple:
    """Uses Entrez ESearch to find all complete mammalian RefSeq mitochondrial genomes."""
    print("Querying NCBI for mammalian RefSeq mitochondrial genomes...")
    query = '"Mammalia"[Organism] AND "RefSeq"[Keyword] AND "complete genome"[All Fields] AND "mitochondrion"[Filter]'
    encoded_query = urllib.parse.quote(query)

    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={encoded_query}&usehistory=y&email={ENTREZ_EMAIL}"

    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, context=ctx) as response:
        xml_data = response.read()

    root = ET.fromstring(xml_data)
    count = int(root.find("Count").text)
    webenv = root.find("WebEnv").text
    query_key = root.find("QueryKey").text

    print(f"Found {count} mammalian mitochondrial genomes.")
    return count, webenv, query_key


def fetch_and_route_cds(
    count: int, webenv: str, query_key: str, hgnc_ref: GeneReference
):
    """Fetches CDS records in batches and routes them to canonical gene files."""

    # Pre-open file handles for the 13 canonical mtDNA genes
    mt_genes = [
        data["primary_symbol"]
        for sym, data in hgnc_ref.lookup.items()
        if data.get("primary_symbol", "").startswith("MT-")
    ]
    mt_genes = sorted(list(set(mt_genes)))

    file_handles = {}
    for gene in mt_genes:
        file_path = RAW_MT_DIR / f"{gene}_raw_cds.fasta"
        # Overwrite if exists to ensure clean run
        file_handles[gene] = open(file_path, "w", encoding="utf-8")

    total_records = 0
    genes_captured = defaultdict(int)

    print(f"\nDownloading CDS sequences in batches of {BATCH_SIZE}...")

    try:
        for start in range(0, count, BATCH_SIZE):
            print(f"  Fetching records {start} to {min(start + BATCH_SIZE, count)}...")

            url = (
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
                f"db=nuccore&query_key={query_key}&WebEnv={webenv}"
                f"&retstart={start}&retmax={BATCH_SIZE}&rettype=fasta_cds_na&retmode=text"
                f"&email={ENTREZ_EMAIL}"
            )

            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, context=ctx) as response:
                fasta_data = response.read().decode("utf-8")

            # Parse the Multi-FASTA response
            current_header = ""
            current_seq = []

            for line in fasta_data.splitlines():
                line = line.strip()
                if line.startswith(">"):
                    # Process previous record before starting new one
                    if current_header and current_seq:
                        process_record(
                            current_header,
                            "".join(current_seq),
                            hgnc_ref,
                            file_handles,
                            genes_captured,
                        )
                        total_records += 1

                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)

            # Process the very last record in the batch
            if current_header and current_seq:
                process_record(
                    current_header,
                    "".join(current_seq),
                    hgnc_ref,
                    file_handles,
                    genes_captured,
                )
                total_records += 1

            time.sleep(RATE_LIMIT)  # Respect NCBI servers

    finally:
        # Guarantee file handles are closed even if script crashes
        for fh in file_handles.values():
            fh.close()

    return total_records, genes_captured


def process_record(
    header: str, sequence: str, hgnc_ref: GeneReference, file_handles: dict, stats: dict
):
    """Extracts the gene name from NCBI header, standardizes it, and writes to disk."""

    # NCBI fasta_cds_na headers look like:
    # >lcl|NC_001665.2_cds_NP_008222.1_1 [gene=ND1] [db_xref=GeneID:807892] ...
    match = re.search(r"\[gene=([^\]]+)\]", header)
    if not match:
        return  # Skip if no gene annotation

    ncbi_gene_name = match.group(1).upper()

    # Use the HGNC dictionary to resolve aliases (e.g. "ND1" -> "MT-ND1")
    gene_data = hgnc_ref.get_gene_data(ncbi_gene_name)
    if not gene_data:
        return

    canonical_symbol = gene_data.get("primary_symbol")

    if canonical_symbol in file_handles:
        # Simplify header for downstream aligners (e.g. ">NC_001665.2_Pan_troglodytes")
        # We extract the main accession for the ID
        accession = header.split("_cds_")[0].replace(">lcl|", "")
        clean_header = f">{accession} | {canonical_symbol}"

        file_handles[canonical_symbol].write(f"{clean_header}\n{sequence}\n")
        stats[canonical_symbol] += 1


def main():
    print("Initializing Mitochondrial CDS Download...\n")

    try:
        hgnc_file = get_latest(DATA_DIR, "Canonical_OXPHOS_Subunits_HGNC*.csv")
        hgnc_ref = GeneReference(hgnc_file)
    except Exception as e:
        print(f"Error loading HGNC reference: {e}")
        return

    count, webenv, query_key = search_ncbi_mammals()

    if count == 0:
        print("Search failed to find any genomes.")
        return

    total_cds, stats = fetch_and_route_cds(count, webenv, query_key, hgnc_ref)

    print(f"\n{'='*50}")
    print("mtDNA DOWNLOAD SUMMARY")
    print(f"{'='*50}")
    print(f"Total valid CDS records processed: {total_cds}")

    for gene, captured_count in sorted(stats.items()):
        print(f"  {gene:<10}: {captured_count} mammalian sequences")

    print(f"\nSaved raw unaligned FASTA files to: {RAW_MT_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
