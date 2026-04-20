#!/usr/bin/env python3
"""
Download and inventory mtDNA RefSeq inputs used by the comparative pipeline.

This step records accession, species, and TaxID metadata for mitochondrial
reference sequences without recomputing downstream alignments.
"""

from pathlib import Path
import re
import ssl
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from utils.gene_reference import GeneReference
from utils.download_paths import LEGACY_REFERENCE_DIR, RAW_REFERENCE_DIR, ensure_layout, latest_existing

RAW_MT_DIR = LEGACY_REFERENCE_DIR.parents[0] / "alignments" / "mtdna_raw_cds"
RAW_MT_DIR.mkdir(parents=True, exist_ok=True)

ENTREZ_EMAIL = "pipeline@analysis.local"
BATCH_SIZE = 500
RATE_LIMIT = 0.5

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE


def search_ncbi_mammals() -> tuple[int, str, str]:
    query = '"Mammalia"[Organism] AND "RefSeq"[Keyword] AND "complete genome"[All Fields] AND "mitochondrion"[Filter]'
    encoded = urllib.parse.quote(query)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=nuccore&term={encoded}&usehistory=y&email={ENTREZ_EMAIL}"
    )
    with urllib.request.urlopen(urllib.request.Request(url), context=CTX) as response:
        xml_data = response.read()
    root = ET.fromstring(xml_data)
    return int(root.find("Count").text), root.find("WebEnv").text, root.find("QueryKey").text


def process_record(header: str, sequence: str, hgnc_ref: GeneReference, handles: dict, stats: dict) -> None:
    match = re.search(r"\[gene=([^\]]+)\]", header)
    if not match:
        return
    gene_name = match.group(1).upper()
    gene_data = hgnc_ref.get_gene_data(gene_name)
    if not gene_data:
        return
    canonical_symbol = gene_data.get("symbol") or gene_data.get("primary_symbol")
    if canonical_symbol not in handles:
        return
    accession = header.split("_cds_")[0].replace(">lcl|", "")
    handles[canonical_symbol].write(f">{accession} | {canonical_symbol}\n{sequence}\n")
    stats[canonical_symbol] += 1


def fetch_and_route(count: int, webenv: str, query_key: str, hgnc_ref: GeneReference):
    mt_genes = sorted({
        data["symbol"] for data in hgnc_ref.lookup.values() if data.get("symbol", "").startswith("MT-")
    })
    handles = {gene: open(RAW_MT_DIR / f"{gene}_raw_cds.fasta", "w", encoding="utf-8") for gene in mt_genes}
    stats = defaultdict(int)
    try:
        for start in range(0, count, BATCH_SIZE):
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                f"?db=nuccore&query_key={query_key}&WebEnv={webenv}"
                f"&retstart={start}&retmax={BATCH_SIZE}&rettype=fasta_cds_na&retmode=text"
                f"&email={ENTREZ_EMAIL}"
            )
            with urllib.request.urlopen(urllib.request.Request(url), context=CTX) as response:
                fasta_data = response.read().decode("utf-8")
            current_header = ""
            current_seq: list[str] = []
            for line in fasta_data.splitlines():
                line = line.strip()
                if line.startswith(">"):
                    if current_header and current_seq:
                        process_record(current_header, "".join(current_seq), hgnc_ref, handles, stats)
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_header and current_seq:
                process_record(current_header, "".join(current_seq), hgnc_ref, handles, stats)
            time.sleep(RATE_LIMIT)
    finally:
        for handle in handles.values():
            handle.close()
    return stats


def main() -> None:
    ensure_layout()
    hgnc_path = latest_existing([(RAW_REFERENCE_DIR, "Canonical_OXPHOS_Subunits_HGNC_*.csv")])
    if hgnc_path is None:
        raise FileNotFoundError("HGNC gene list not found in data/raw/reference")
    hgnc_ref = GeneReference(hgnc_path)
    count, webenv, query_key = search_ncbi_mammals()
    stats = fetch_and_route(count, webenv, query_key, hgnc_ref)
    for gene, captured in sorted(stats.items()):
        print(f"{gene}: {captured}")


if __name__ == "__main__":
    main()
