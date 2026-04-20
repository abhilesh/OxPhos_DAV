#!/usr/bin/env python3
"""
Fetch and cache amino-acid property reference tables.

These tables support downstream substitution-property annotation.
"""

from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

try:
    from aaindex import aaindex1, aaindex2
except Exception:
    aaindex1 = aaindex2 = None

from utils.download_paths import DERIVED_REFERENCE_DIR, LEGACY_REFERENCE_DIR, TODAY, ensure_layout, latest_existing, sync_compat_file

OUT_JSON = DERIVED_REFERENCE_DIR / f"aaindex_properties_{TODAY}.json"

TARGETS = {
    "hydrophobicity_kd": "KYTJ820101",
    "volume_fauvereaux": "FAUJ880103",
    "blosum62": "HENS920102",
    "miyata_distance": "MIYT790101",
}


def extract_data(record):
    for key in ["values", "matrix", "index", "mutations"]:
        if key in record:
            return record[key]
    for _, value in record.items():
        if isinstance(value, dict) and "A" in value:
            return value
    raise KeyError("Could not locate AAIndex payload.")


def main() -> None:
    ensure_layout()
    if aaindex1 is None or aaindex2 is None:
        existing = latest_existing(
            [
                (DERIVED_REFERENCE_DIR, "aaindex_properties_*.json"),
                (LEGACY_REFERENCE_DIR, "aaindex_properties_*.json"),
            ]
        )
        if existing is None:
            raise ModuleNotFoundError(
                "aaindex package is unavailable and no existing aaindex_properties_*.json file was found."
            )
        sync_compat_file(existing, DERIVED_REFERENCE_DIR / existing.name)
        sync_compat_file(existing, LEGACY_REFERENCE_DIR / existing.name)
        print(f"Reused existing AAIndex cache: {existing}")
        return

    payload = {
        "metadata": {
            "source": "https://www.genome.jp/aaindex/",
            "download_date": TODAY,
            "accessions": TARGETS,
        },
        "indices": {
            "hydrophobicity_kd": extract_data(aaindex1[TARGETS["hydrophobicity_kd"]]),
            "volume": extract_data(aaindex1[TARGETS["volume_fauvereaux"]]),
        },
        "matrices": {
            "blosum62": extract_data(aaindex2[TARGETS["blosum62"]]),
            "miyata_distance": extract_data(aaindex2[TARGETS["miyata_distance"]]),
        },
    }
    with open(OUT_JSON, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    sync_compat_file(OUT_JSON, LEGACY_REFERENCE_DIR / OUT_JSON.name)
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
