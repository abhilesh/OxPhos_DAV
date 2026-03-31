import urllib.request
import urllib.error
from pathlib import Path

# ── Config ──────────────────────────────────────────────────────────────────
# Resolve data/structures relative to this script's location
ROOT = Path(__file__).resolve().parent.parent.parent
STRUC_DIR = ROOT / "data" / "structures"

STRUC_DIR.mkdir(parents=True, exist_ok=True)

RCSB_CIF_URL = "https://files.rcsb.org/download/{}.cif"

# 2024-2026 Human OXPHOS structures [cite: 4]
NEW_PDBS = {
    "CI": "9I4I",
    "CII": "8GS8",
    "CIII": "9HZL",
    "CIV": "9I6F",
    "CV": "8H9S",
}


def download_cif(pdb_id, out_path):
    """Download mmCIF file from RCSB. Returns True on success."""
    url = RCSB_CIF_URL.format(pdb_id.upper())
    try:
        print(f"  Downloading {pdb_id.upper()} from {url}...")
        urllib.request.urlretrieve(url, out_path)
        return True
    except urllib.error.HTTPError as e:
        print(f"  HTTP Error {e.code} for {pdb_id}: {e.reason}")
    except Exception as e:
        print(f"  Unexpected error downloading {pdb_id}: {e}")
    return False


def main():
    # Ensure target directory exists
    STRUC_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Target directory: {STRUC_DIR}")
    print(f"Starting download of {len(NEW_PDBS)} structures...\n")

    results = {"success": [], "failed": []}

    for complex_name, pdb_id in NEW_PDBS.items():
        filename = f"{pdb_id.upper()}.cif"
        out_path = STRUC_DIR / filename

        print(f"Processing {complex_name} ({pdb_id.upper()}):")

        if out_path.exists():
            print(f"  [Skipped] {filename} already exists.")
            results["success"].append(pdb_id)
            continue

        if download_cif(pdb_id, out_path):
            print(f"  [Success] Saved to {out_path}")
            results["success"].append(pdb_id)
        else:
            results["failed"].append(pdb_id)

    # ── Summary ──────────────────────────────────────────────────────────────
    print(f"\n{'='*40}")
    print(f"DOWNLOAD SUMMARY")
    print(f"  Total structures: {len(NEW_PDBS)}")
    print(f"  Successful:       {len(results['success'])}")
    print(f"  Failed:           {len(results['failed'])}")

    if results["failed"]:
        print(f"  Check connectivity for: {', '.join(results['failed'])}")
    print(f"{'='*40}\n")


if __name__ == "__main__":
    main()
