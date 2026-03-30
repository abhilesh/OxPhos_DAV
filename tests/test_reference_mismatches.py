import json
import re
from pathlib import Path
from collections import Counter
from utils.alignment_parser import AlignmentParser

# ==== Configuration ====
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
CURATED_DIR = DATA_DIR / "annotations" / "curated"
MT_JSON = CURATED_DIR / "cdar_classifications_mtDNA.json"

MT_AA_DIR = DATA_DIR / "alignments" / "mtdna_aa"
MT_NT_DIR = DATA_DIR / "alignments" / "mtdna_codon"

# Canonical rCRS boundaries
MT_GENE_COORDS = {
    "MT-ND1": (3307, 4262, "+"),
    "MT-ND2": (4470, 5511, "+"),
    "MT-CO1": (5904, 7445, "+"),
    "MT-CO2": (7586, 8269, "+"),
    "MT-ATP8": (8366, 8572, "+"),
    "MT-ATP6": (8527, 9207, "+"),
    "MT-CO3": (9207, 9990, "+"),
    "MT-ND3": (10059, 10404, "+"),
    "MT-ND4L": (10470, 10766, "+"),
    "MT-ND4": (10760, 12137, "+"),
    "MT-ND5": (12337, 14148, "+"),
    "MT-ND6": (14149, 14673, "-"),
    "MT-CYB": (14747, 15887, "+"),
}


def resolve_nt_pos(genomic_pos, locus):
    if locus not in MT_GENE_COORDS:
        return None
    start, end, strand = MT_GENE_COORDS[locus]
    return (genomic_pos - start + 1) if strand == "+" else (end - genomic_pos + 1)


def main():
    if not MT_JSON.exists():
        print("Run 00_classify_DAV.py first to generate the classification JSON.")
        return

    with open(MT_JSON, "r") as f:
        variants = json.load(f)

    # Focus only on the 137 mismatches
    mismatches = [v for v in variants if not v.get("ref_allele_match", True)]
    print(f"Auditing {len(mismatches)} mtDNA reference mismatches...\n")

    loaded_parsers = {}
    error_types = Counter()

    for var in mismatches:
        locus = var["locus"].split("/")[0]
        g_pos = var["genomic_pos"]
        expected = var["ref_nt"]

        if locus not in loaded_parsers:
            aa_p = MT_AA_DIR / f"{locus}_aa_alignment.fasta"
            nt_p = MT_NT_DIR / f"{locus}_codon_alignment.fasta"
            if not aa_p.exists():
                continue
            loaded_parsers[locus] = AlignmentParser(aa_p, nt_p, "mtDNA")

        parser = loaded_parsers[locus]
        nt_pos = resolve_nt_pos(g_pos, locus)

        if not nt_pos or nt_pos not in parser.nt_map:
            error_types["Position out of alignment bounds"] += 1
            continue

        col = parser.nt_map[nt_pos]
        actual = parser.nt_alignment[parser.ref_header][col]

        # Check for systematic shifts
        prev_base = parser.nt_alignment[parser.ref_header][col - 1] if col > 0 else None
        next_base = (
            parser.nt_alignment[parser.ref_header][col + 1]
            if col < len(parser.nt_alignment[parser.ref_header]) - 1
            else None
        )

        if actual == expected:
            error_types["False Positive Mismatch (Logic Error)"] += 1
        elif prev_base == expected:
            error_types["Systematic Shift: Off-by-one (-1)"] += 1
        elif next_base == expected:
            error_types["Systematic Shift: Off-by-one (+1)"] += 1
        else:
            error_types["Total Sequence/Coordinate Disconnect"] += 1
            if error_types["Total Sequence/Coordinate Disconnect"] <= 5:
                print(f"Sample Mismatch in {locus}:")
                print(f"  Genomic {g_pos} -> CDS {nt_pos}")
                print(f"  Expected {expected}, found {actual}")
                print(
                    f"  Surrounding sequence: {parser.nt_alignment[parser.ref_header][col-2:col+3]}"
                )

    print("\n" + "=" * 40)
    print("MISMATCH DIAGNOSTIC SUMMARY")
    print("=" * 40)
    for err, count in error_types.items():
        print(f"{err:<40}: {count}")


if __name__ == "__main__":
    main()
