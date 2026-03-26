#!/usr/bin/env python3
"""
Download all disease associated variant data

Primary sources:
- MITOMAP: https://www.mitomap.org/MITOMAP/MutationsCoding
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/docs/data_file_download/

"""

import csv
import json
import io
import re
import ssl
import sys
import time
import zipfile
import requests
import urllib.request
import pandas as pd
from datetime import date
from pathlib import Path
from utils.utils import get_cached_path, make_read_only

# Get today's date for metadata
today = date.today().isoformat()

# Paths
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "annotations" / "raw"

RAW_DIR.mkdir(parents=True, exist_ok=True)

# Files
GENE_LIST_TARGET = DATA_DIR / f"Canonical_OXPHOS_Subunits_HGNC_{today}.csv"
MITOMAP_FILE_TARGET = RAW_DIR / f"MITOMAP_CodingVariants_{today}.tsv"
CLINVAR_FILE_TARGET = RAW_DIR / f"ClinVar_VariantSummary_{today}.txt.gz"
MITIMPACT_FILE_TARGET = RAW_DIR / f"MitImpact_db_{today}.txt.zip"
PHYLOTREE_FILE_TARGET = RAW_DIR / f"PhyloTree_build_17_{today}.zip"
GNOMAD_FILE_TARGET = RAW_DIR / f"gnomAD_{today}.vcf.bgz"
MYVARIANT_FILE_TARGET = RAW_DIR / f"MyVariant_dbNSFP_gnomAD_{today}.json"


GENE_LIST = get_cached_path(GENE_LIST_TARGET, today)
MITOMAP_FILE = get_cached_path(MITOMAP_FILE_TARGET, today)
CLINVAR_FILE = get_cached_path(CLINVAR_FILE_TARGET, today)
MITIMPACT_FILE = get_cached_path(MITIMPACT_FILE_TARGET, today)
PHYLOTREE_FILE = get_cached_path(PHYLOTREE_FILE_TARGET, today)
MYVARIANT_FILE = get_cached_path(MYVARIANT_FILE_TARGET, today)


# URLs
MITOMAP_URL = "https://www.mitomap.org/cgi-bin/disease.cgi"
CLINVAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)
MITIMPACT_URL = "https://mitimpact.css-mendel.it/cdn/MitImpact_db_3.1.3.txt.zip"
PHYLOTREE_URL = "https://www.phylotree.org/builds/mtDNA_tree_Build_17%20-%20rCRS-oriented%20version.zip"
MYVARIANT_URL = "https://myvariant.info/v1/query"


# SSL context (macOS sometimes needs relaxed verification for NCBI/Senckenberg)
_CTX = ssl.create_default_context()
_CTX.check_hostname = False
_CTX.verify_mode = ssl.CERT_NONE


def download(url, dest, description=""):
    """Stream-download url to destination"""

    print(f" Downloading {description or url} ...", end=" ", flush=True)
    try:
        with urllib.request.urlopen(url, context=_CTX) as resp, open(dest, "wb") as f:
            chunk = 65536
            total = 0
            while True:
                buf = resp.read(chunk)
                if not buf:
                    break
                f.write(buf)
                total += len(buf)
        size_mb = total / 1_048_576
        print(f"Done ({size_mb:.2f} MB).")
        return True, total
    except Exception as e:
        dest.unlink(missing_ok=True)
        print(f"Failed ({e}).")
        return False, str(e)


# --- 1. Download MITOMAP coding variants ---

if not MITOMAP_FILE.exists():
    print("Downloading MITOMAP variant summary...")
    headers = {
        "User-Agent": "Bioinformatics-Research-Pipeline/1.0 (Contact: research@example.com)"
    }

    try:
        # Using requests for MITOMAP as it handles CGI streams better
        response = requests.get(MITOMAP_URL, headers=headers, timeout=60)
        response.raise_for_status()

        with open(MITOMAP_FILE_TARGET, "wb") as f:
            f.write(response.content)
        print(f"Successfully saved MITOMAP data to {MITOMAP_FILE_TARGET.name}")
        make_read_only(MITOMAP_FILE_TARGET)
        MITOMAP_FILE = MITOMAP_FILE_TARGET

    except requests.exceptions.RequestException as e:
        print(f"FAILED MITOMAP download: {e}")
else:
    print(f"Using cached MITOMAP: {MITOMAP_FILE.name}")

# --- 2. Download ClinVar variant summary ---

if not CLINVAR_FILE.exists():
    print("Downloading ClinVar variant summary...")
    try:
        success, detail = download(
            CLINVAR_URL, CLINVAR_FILE_TARGET, description="ClinVar Variant Summary"
        )
        if success:
            make_read_only(CLINVAR_FILE_TARGET)
            CLINVAR_FILE = CLINVAR_FILE_TARGET
    except Exception as e:
        print(f"FAILED: {e}")
        sys.exit(1)
else:
    print(f"Using cached: {CLINVAR_FILE.name}")

# --- 3. MitImpact3D APOGEE2 scores ---

if not MITIMPACT_FILE.exists():
    print("Downloading MitImpact3D...")
    try:
        success, detail = download(
            MITIMPACT_URL, MITIMPACT_FILE_TARGET, description="MitImpact3D Database"
        )
        if success:
            make_read_only(MITIMPACT_FILE_TARGET)
            MITIMPACT_FILE = MITIMPACT_FILE_TARGET
    except Exception as e:
        print(f"FAILED: {e}")
        sys.exit(1)

# --- 4. PhyloTree download ---
if not PHYLOTREE_FILE.exists():
    print("Downloading PhyloTree build 17...")
    try:
        success, detail = download(
            PHYLOTREE_URL, PHYLOTREE_FILE_TARGET, description="PhyloTree Build 17"
        )
        if success:
            make_read_only(PHYLOTREE_FILE_TARGET)
    except Exception as e:
        print(f"FAILED: {e}")
        sys.exit(1)


def fetch_myvariant_annotations(url):
    if not GENE_LIST.exists():
        print(
            f"\n[WARNING] GENE_LIST ({GENE_LIST.name}) not found! Skipping MyVariant fetch."
        )
        return False

    print("\nFetching dbNSFP & nucDNA gnomAD data via MyVariant API...")
    genes = []
    with open(GENE_LIST, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("Approved symbol", "").strip()
            if sym:
                genes.append(sym)

    if not genes:
        print("No genes found. Cannot query MyVariant.")
        return False

    all_results = []

    # Query one gene at a time to safely handle pagination
    for gene in genes:
        skip = 0
        limit = 1000

        while True:
            params = {
                "q": gene,
                "scopes": "symbol",
                "fields": "dbnsfp,gnomad_genome,gnomad_exome,vcf",
                "species": "human",
                "size": limit,
                "from": skip,
            }

            try:
                resp = requests.post(url, data=params, timeout=120)
                resp.raise_for_status()
                data = resp.json()

                # MyVariant returns a list of hits for POST requests
                if isinstance(data, list):
                    hits = data
                else:
                    hits = data.get("hits", [])

                if not hits:
                    break

                all_results.extend(hits)

                # If we retrieved fewer hits than the limit, we have reached the end for this gene
                if len(hits) < limit:
                    break

                skip += limit
                time.sleep(0.2)

            except Exception as e:
                print(f"      [Error] Failed on {gene} at skip {skip}: {e}")
                break

        time.sleep(0.5)

    try:
        with open(MYVARIANT_FILE_TARGET, "w") as f:
            json.dump(all_results, f, indent=2)
        make_read_only(MYVARIANT_FILE_TARGET)
        print(
            f"  Success! Saved {len(all_results)} variant records to {MYVARIANT_FILE_TARGET.name}"
        )
        return True
    except Exception as e:
        print(f"  FAILED to save MyVariant data: {e}")
        return False


if not MYVARIANT_FILE.exists():
    success = fetch_myvariant_annotations(MYVARIANT_URL)
    if success:
        MYVARIANT_FILE = MYVARIANT_FILE_TARGET
else:
    print(f"\nUsing cached MyVariant API annotations: {MYVARIANT_FILE.name}")
