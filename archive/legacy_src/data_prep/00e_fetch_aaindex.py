import json
from datetime import date
from pathlib import Path
from aaindex import aaindex1, aaindex2

# ==== Configuration ====
today = date.today().isoformat()

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
REF_DIR = DATA_DIR / "reference"
REF_DIR.mkdir(parents=True, exist_ok=True)

OUT_JSON = REF_DIR / f"aaindex_properties_{today}.json"

# Target AAindex Accessions
TARGETS = {
    "hydrophobicity_kd": "KYTJ820101",  # Kyte-Doolittle hydrophobicity
    "volume_fauvereaux": "FAUJ880103",  # van der Waals volume (Fauvereaux 1988)
    "blosum62": "HENS920102",  # BLOSUM62 substitution matrix
    "miyata_distance": "MIYT790101",  # Miyata distance matrix
}


def extract_data(record):
    """Safely extracts the numerical data/matrix from an AAIndex record
    regardless of the key name used by the underlying library."""

    # Common keys the library might use for the data payload
    for key in ["values", "matrix", "index", "mutations"]:
        if key in record:
            return record[key]

    # Fallback: Find the first sub-dictionary that contains amino acid keys (e.g., 'A' for Alanine)
    for key, val in record.items():
        if isinstance(val, dict) and "A" in val:
            return val

    # If all else fails, print what keys are actually available to help debug
    raise KeyError(
        f"Could not locate data matrix. Available keys: {list(record.keys())}"
    )


def fetch_and_build_cache():
    """Fetches specific scales from AAindex and caches them as JSON."""

    cache_data = {
        "metadata": {
            "source": "https://www.genome.jp/aaindex/",
            "download_date": today,
            "accessions": TARGETS,
        },
        "indices": {},
        "matrices": {},
    }

    print("Extracting AAindex1 properties (single-value scales)...")
    cache_data["indices"]["hydrophobicity_kd"] = extract_data(
        aaindex1[TARGETS["hydrophobicity_kd"]]
    )
    cache_data["indices"]["volume"] = extract_data(
        aaindex1[TARGETS["volume_fauvereaux"]]
    )

    print("Extracting AAindex2 properties (matrices)...")
    cache_data["matrices"]["blosum62"] = extract_data(aaindex2[TARGETS["blosum62"]])
    cache_data["matrices"]["miyata_distance"] = extract_data(
        aaindex2[TARGETS["miyata_distance"]]
    )

    print(f"Saving data to {OUT_JSON}...")
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(cache_data, f, indent=4)

    return True


def main():
    print(f"Targeting reference directory: {REF_DIR}")
    success = fetch_and_build_cache()
    if success:
        print(f"\nSuccessfully saved AAindex properties to {OUT_JSON.name}")


if __name__ == "__main__":
    main()
