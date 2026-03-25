#!/usr/bin/env python3
"""
Download all disease associated variant data

Primary sources:
- MITOMAP: https://www.mitomap.org/MITOMAP/MutationsCoding
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/docs/data_file_download/

"""

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


# Get today's date for metadata
today = date.today().isoformat()

# Paths
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "annotations" / "raw"

RAW_DIR.mkdir(parents=True, exist_ok=True)


def get_cached_path(target_path: Path, date_str: str) -> Path:
    """
    Checks if a previously downloaded version of the file exists by replacing
    the date string with a wildcard '*' and finding the most recent match.
    """
    pattern = target_path.name.replace(date_str, "*")
    matches = list(target_path.parent.glob(pattern))
    if matches:
        return max(matches, key=lambda p: p.stat().st_mtime)
    return target_path


# Files
GENE_LIST_TARGET = DATA_DIR / f"Canonical_OXPHOS_Subunits_HGNC_{today}.csv"
MITOMAP_FILE_TARGET = DATA_DIR / f"MITOMAP_CodingVariants_{today}.tsv"
CLINVAR_FILE_TARGET = DATA_DIR / f"ClinVar_VariantSummary_{today}.txt.gz"

GENE_LIST = get_cached_path(GENE_LIST_TARGET, today)
MITOMAP_FILE = get_cached_path(MITOMAP_FILE_TARGET, today)
CLINVAR_FILE = get_cached_path(CLINVAR_FILE_TARGET, today)

# URLs
MITOMAP_URL = "https://www.mitomap.org/cgi-bin/disease.cgi"
CLINVAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)

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
            CLINVAR_URL, CLINVAR_FILE, description="ClinVar Variant Summary"
        )
    except Exception as e:
        print(f"FAILED: {e}")
        sys.exit(1)
else:
    print(f"Using cached: {CLINVAR_FILE.name}")
