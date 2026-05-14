#!/usr/bin/env python3
"""
Download OXPHOS structure mmCIF files listed in the structure manifest.

Primary, validation, and reference structures are all downloaded into
`data/structures/` so the structural mapping stage can evaluate contact
robustness across models while retaining the curated structural provenance map.

Run from project root:
  python src/data_download/00j_download_structures.py
"""

import csv
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STRUCTURE_DIR = ROOT / "data" / "structures"
MANIFEST = ROOT / "data" / "reference" / "structure_model_manifest.tsv"
if not MANIFEST.exists():
    MANIFEST = ROOT / "data" / "derived" / "reference" / "structure_model_manifest.tsv"

RCSB_CIF_URL = "https://files.rcsb.org/download/{}.cif"


def load_manifest() -> list[dict]:
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Missing structure manifest: {MANIFEST}")
    rows = []
    with open(MANIFEST, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if str(row.get("active", "")).strip().lower() != "true":
                continue
            rows.append(row)
    rows.sort(key=lambda item: (item["complex_id"], int(item.get("priority") or 999), item["pdb_id"]))
    return rows


def download_cif(pdb_id: str, out_path: Path) -> tuple[bool, str]:
    url = RCSB_CIF_URL.format(pdb_id.upper())
    try:
        print(f"  Downloading {pdb_id.upper()} from {url} ...")
        urllib.request.urlretrieve(url, out_path)
        return True, "downloaded"
    except urllib.error.HTTPError as exc:
        return False, f"http_{exc.code}"
    except Exception as exc:
        return False, f"error:{exc}"


def main() -> None:
    STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)
    manifest_rows = load_manifest()

    print(f"Target directory : {STRUCTURE_DIR}")
    print(f"Structure panel  : {len(manifest_rows)} active manifest rows")

    success = []
    failed = []

    for row in manifest_rows:
        complex_id = row["complex_id"]
        pdb_id = row["pdb_id"].upper()
        role = row["role"]
        filename = f"{pdb_id}.cif"
        out_path = STRUCTURE_DIR / filename

        print(f"\nProcessing {complex_id} | {pdb_id} | {role}")
        if out_path.exists():
            print(f"  [Skipped] {filename} already exists")
            success.append((pdb_id, "already_present"))
            continue

        ok, detail = download_cif(pdb_id, out_path)
        if ok:
            print(f"  [Success] Saved to {out_path}")
            success.append((pdb_id, detail))
        else:
            print(f"  [Failed] {pdb_id} ({detail})")
            failed.append((pdb_id, detail))

    print(f"\n{'=' * 48}")
    print("STRUCTURE DOWNLOAD SUMMARY")
    print(f"  Active manifest rows : {len(manifest_rows)}")
    print(f"  Successful           : {len(success)}")
    print(f"  Failed               : {len(failed)}")
    if failed:
        for pdb_id, detail in failed:
            print(f"    - {pdb_id}: {detail}")
    print(f"{'=' * 48}")


if __name__ == "__main__":
    main()
