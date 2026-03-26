import pandas as pd
import re
import io
from pathlib import Path

# Standard rCRS coordinates for MT-OXPHOS genes
MT_GENES = [
    ("MT-ND1", 3307, 4262),
    ("MT-ND2", 4470, 5511),
    ("MT-CO1", 5904, 7445),
    ("MT-CO2", 7586, 8269),
    ("MT-ATP8", 8366, 8572),
    ("MT-ATP6", 8527, 9207),
    ("MT-CO3", 9207, 9990),
    ("MT-ND3", 10059, 10404),
    ("MT-ND4L", 10470, 10766),
    ("MT-ND4", 10760, 12137),
    ("MT-ND5", 12337, 14148),
    ("MT-ND6", 14149, 14673),
    ("MT-CYB", 14747, 15887),
]


def get_locus(pos: int) -> str:
    """Map standard rCRS coordinate to mitochondrial gene."""
    genes = [name for name, start, end in MT_GENES if start <= pos <= end]
    return "/".join(genes) if genes else "Non-OXPHOS"


def categorize_aachange(aa_str: str) -> str:
    """Classifies raw MITOMAP aachange strings into standard categories."""
    aa = str(aa_str).strip()
    aa_lower = aa.lower()

    if not aa or aa_lower == "nan":
        return "Missing"
    if "noncoding" in aa_lower:
        return "Noncoding"
    if "syn" in aa_lower:
        return "Synonymous"
    if "frameshift" in aa_lower or "fs" in aa_lower:
        return "Frameshift"
    if "*" in aa or "ter" in aa_lower or "stop" in aa_lower:
        return "Nonsense"

    # Check for standard protein change format (e.g., V113A)
    match = re.match(r"^([a-zA-Z]+)(\d+)([a-zA-Z]+)$", aa)
    if match:
        if match.group(1).upper() == match.group(3).upper():
            return "Synonymous"
        return "Missense"

    return "Other"


def parse_and_generate_stats(file_source) -> pd.DataFrame:
    """Parses MITOMAP CGI output and prints summary statistics."""

    # Added encoding="windows-1252" to handle 0xa0 bytes (non-breaking spaces)
    df = pd.read_csv(
        file_source, sep="\t", on_bad_lines="skip", encoding="windows-1252"
    )
    df.columns = [c.strip().lower() for c in df.columns]

    # 1. Safely parse coordinates and assign Locus
    df["pos_clean"] = pd.to_numeric(df["pos"], errors="coerce").fillna(-1).astype(int)
    df["locus"] = df["pos_clean"].apply(get_locus)

    # 2. Categorize amino acid changes
    df["variant_type"] = df["aachange"].apply(categorize_aachange)

    # --- Print Statistics ---
    print("=" * 50)
    print("MITOMAP RAW DATA STATISTICS")
    print("=" * 50)
    print(f"Total Rows Parsed: {len(df)}\n")

    print("--- By Gene (Locus) ---")
    print(df["locus"].value_counts().to_string())
    print()

    print("--- By Variant Type ---")
    print(df["variant_type"].value_counts().to_string())
    print()

    print("--- By Clinical Status ---")
    status_counts = df["status"].fillna("Missing").value_counts()
    print(status_counts.to_string())
    print()

    return df


if __name__ == "__main__":
    # Define the path to the raw MITOMAP CGI file you downloaded earlier
    # Update this path if your file is named differently or in a different folder
    target_file = Path("data/annotations/raw/MITOMAP_CodingVariants_2026-03-26.tsv")

    if target_file.exists():
        print(f"Processing {target_file.name}...")
        df = parse_and_generate_stats(target_file)

        # 1. Export Non-OXPHOS variants
        non_oxphos = df[df["locus"] == "Non-OXPHOS"]
        non_oxphos.to_csv("inspect_non_oxphos.csv", index=False)
        print(
            f"Exported {len(non_oxphos)} Non-OXPHOS variants to inspect_non_oxphos.csv"
        )

        # 2. Export Missing amino acid variants
        missing_aa = df[df["variant_type"] == "Missing"]
        missing_aa.to_csv("inspect_missing_aa.csv", index=False)
        print(
            f"Exported {len(missing_aa)} Missing AA variants to inspect_missing_aa.csv"
        )

        # 3. Export Synonymous variants
        synonymous = df[df["variant_type"] == "Synonymous"]
        synonymous.to_csv("inspect_synonymous.csv", index=False)
        print(
            f"Exported {len(synonymous)} Synonymous variants to inspect_synonymous.csv"
        )
