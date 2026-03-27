import json
import re
from pathlib import Path
from Bio.Seq import Seq
import sys

# ==== Configuration ====
ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))  # Ensure we can import src modules

from src.utils.alignment_parser import AlignmentParser

DATA_DIR = ROOT / "data"
# Use the output from your previous run so we can target known AA c-DARs
IN_JSON = DATA_DIR / "annotations" / "curated" / "cdar_classifications_nucDNA.json"

TOGA_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_aa"
TOGA_NT_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"


def parse_variant_coordinates(variant: dict) -> tuple:
    aa_str = variant.get("aa_change", "")
    nc_str = variant.get("nc_change", "")

    aa_match = re.search(r"([a-zA-Z]+)(\d+)([a-zA-Z]+)", aa_str)
    if not aa_match:
        return None, None, None, None

    wt_aa = aa_match.group(1)
    aa_pos = int(aa_match.group(2))
    mut_aa = aa_match.group(3)

    nt_pos = None
    if variant["genome"] == "nucDNA":
        nc_match = re.search(r"c\.(\d+)", nc_str)
        if nc_match:
            nt_pos = int(nc_match.group(1))

    if not nt_pos:
        nt_pos = (aa_pos * 3) - 2

    return aa_pos, wt_aa, mut_aa, nt_pos


def main():
    print("Initializing NT-Level Diagnostic Engine...\n")

    with open(IN_JSON, "r") as f:
        variants = json.load(f)

    tested = 0
    loaded_alignments = {}

    for var in variants:
        if tested >= 20:
            break

        # Only test nucDNA variants that we KNOW have AA compensation
        if var.get("genome") != "nucDNA" or not var.get("cdar_aa"):
            continue

        locus = var["locus"].split("/")[0]
        aa_pos, wt_aa, mut_aa, nt_pos = parse_variant_coordinates(var)

        if not aa_pos or not nt_pos:
            continue

        if locus not in loaded_alignments:
            aa_path = TOGA_AA_DIR / f"{locus}_aa_alignment.fasta"
            nt_path = TOGA_NT_DIR / f"{locus}_codon_alignment.fasta"
            if not aa_path.exists() or not nt_path.exists():
                continue
            loaded_alignments[locus] = AlignmentParser(aa_path, nt_path, "nucDNA")

        parser = loaded_alignments[locus]

        # --- THE DIAGNOSTIC EXTRACTION ---
        if nt_pos not in parser.nt_map:
            continue

        col_idx = parser.nt_map[nt_pos]
        codon_start_bio = nt_pos - ((nt_pos - 1) % 3)
        col_start = parser.nt_map[codon_start_bio]

        # 1. Extract WT Codon
        wt_codon = parser.nt_alignment[parser.ref_header][col_start : col_start + 3]

        # 2. Get the ALT allele directly from your JSON
        alt_allele = var.get("alt_nt", "MISSING")

        # 3. Fabricate Mutant Codon
        pos_in_codon = (nt_pos - 1) % 3
        mut_codon = wt_codon[:pos_in_codon] + alt_allele + wt_codon[pos_in_codon + 1 :]

        # 4. Perform the Triangle of Truth Translations
        try:
            translated_wt = str(Seq(wt_codon).translate())
            translated_mut = str(Seq(mut_codon).translate())
        except Exception:
            translated_wt = "ERR"
            translated_mut = "ERR"

        print(f"{'='*50}")
        print(f"VARIANT: {locus} {var.get('nc_change')} ({var.get('aa_change')})")
        print(f"JSON 'alt_nt' value : '{alt_allele}'")
        print(f"-" * 50)
        print(
            f"Alignment WT Codon  : {wt_codon} -> Translates to: {translated_wt} (Expected: {wt_aa})"
        )

        # The ultimate check:
        translation_match = "PASS" if translated_mut == mut_aa else "FAIL"
        print(
            f"Fabricated Mut Codon: {mut_codon} -> Translates to: {translated_mut} (Expected: {mut_aa})"
        )
        print(f"Triangle of Truth   : {translation_match}")

        tested += 1


if __name__ == "__main__":
    main()
