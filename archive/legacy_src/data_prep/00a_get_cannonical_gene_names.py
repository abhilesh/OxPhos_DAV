#!/usr/bin/env python3
"""
Download cannonical gene names for OXPHOS genes from HGNC
Source: https://www.genenames.org/data/genegroup/#!/group/639
"""

import urllib.request
import urllib.parse
import json
from pathlib import Path
from datetime import date

# Get today's date for metadata
today = date.today().isoformat()

# Paths
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
REF_DIR = DATA_DIR / "reference"

REF_DIR.mkdir(parents=True, exist_ok=True)
GENE_LIST = REF_DIR / f"Canonical_OXPHOS_Subunits_HGNC_{today}.csv"


def fetch_hgnc_group_genes(group_id, output_path):
    # Use the website's explicit download endpoint (returns TSV)
    # type=branch includes the group and any subgroups
    url = f"https://www.genenames.org/cgi-bin/genegroup/download?id={group_id}&type=branch"

    req = urllib.request.Request(url)

    try:
        with urllib.request.urlopen(req) as response:
            content = response.read().decode("utf-8")

            if not content.strip():
                print("Warning: Received empty response.")
                return False

            # Ensure the output directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Write the complete TSV content to the file
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(content)

            return True

    except Exception as e:
        print(f"Error fetching data: {e}")
        return False


if __name__ == "__main__":
    group_id = 639
    print(f"Downloading TSV for HGNC group ID {group_id} to {GENE_LIST}...")

    success = fetch_hgnc_group_genes(group_id, GENE_LIST)

    if success:
        print("Download complete.")
    else:
        print("Download failed.")
