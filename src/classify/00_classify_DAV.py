import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Seq import Seq
from utils.alignment_parser import AlignmentParser

# ==== Configuration ====
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
CURATED_DIR = DATA_DIR / "annotations" / "curated"

# Input: separate curated files
MT_CURATED = CURATED_DIR / "mtDNA_annotations_curated.json"
NUC_CURATED = CURATED_DIR / "nucDNA_annotations_curated.json"

# Alignment Directories
TOGA_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_aa"
TOGA_NT_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"
MT_AA_DIR = DATA_DIR / "alignments" / "mtdna_aa"
MT_NT_DIR = DATA_DIR / "alignments" / "mtdna_codon"

# Output: split files only
OUT_JSON_MT = CURATED_DIR / "cdar_classifications_mtDNA.json"
OUT_JSON_NUC = CURATED_DIR / "cdar_classifications_nucDNA.json"


def parse_variant_coordinates(variant: dict) -> tuple:
    """Extracts biological AA and exact NT coordinates from standard variant strings."""
    aa_str = variant.get("aa_change", "")
    nc_str = variant.get("nc_change", "")

    # Extract AA Position and Mutant AA (e.g., "L109F" -> pos: 109, mut: "F")
    aa_match = re.search(r"[a-zA-Z]+(\d+)([a-zA-Z]+)", aa_str)
    if not aa_match:
        return None, None, None

    aa_pos = int(aa_match.group(1))
    mut_aa = aa_match.group(2)

    # Extract Exact NT Position
    nt_pos = None
    if variant["genome"] == "nucDNA":
        # nucDNA variants use c. coordinates (e.g., c.327G>C), which perfectly match our 1-indexed CDS alignment
        nc_match = re.search(r"c\.(\d+)", nc_str)
        if nc_match:
            nt_pos = int(nc_match.group(1))

    # Fallback (This will keep the script from crashing if a variant is missing a c. string,
    # but we will need to handle mtDNA absolute mapping separately when MACSE finishes)
    if not nt_pos:
        nt_pos = (aa_pos * 3) - 2

    return aa_pos, mut_aa, nt_pos


def get_alignment_paths(locus: str, genome: str) -> tuple:
    """Routes to the correct FASTA files based on the genome."""
    if genome == "nucDNA":
        return (
            TOGA_AA_DIR / f"{locus}_aa_alignment.fasta",
            TOGA_NT_DIR / f"{locus}_codon_alignment.fasta",
        )
    else:
        return (
            MT_AA_DIR / f"{locus}_aa_alignment.fasta",
            MT_NT_DIR / f"{locus}_codon_alignment.fasta",
        )


def check_ref_allele(parser: AlignmentParser, nt_pos: int, expected_ref: str) -> bool:
    """Verifies the variant ref allele matches the human reference in the alignment."""
    if nt_pos not in parser.nt_map:
        return False
    col_idx = parser.nt_map[nt_pos]
    ref_allele_in_alignment = parser.nt_alignment[parser.ref_header][col_idx]
    return ref_allele_in_alignment == expected_ref.upper()


def process_variants(variants: list, genome: str, loaded_alignments: dict) -> tuple:
    """
    Classifies a list of curated variants as cDAR/uDAR.
    Returns (enriched_list, cdar_stats, ref_mismatch_count).
    cdar_stats: Counter keyed by tier -> {"Total", "aa_cDAR", "nt_cDAR"}.
    """
    enriched = []
    ref_mismatch = 0
    # Per-tier counters: tier -> {Total, aa_cDAR, nt_cDAR}
    tier_stats = defaultdict(lambda: {"Total": 0, "aa_cDAR": 0, "nt_cDAR": 0})
    global_stats = {"Total": 0, "aa_cDAR": 0, "nt_cDAR": 0}

    total = len(variants)
    for idx, var in enumerate(variants):
        if idx % 500 == 0:
            print(f"  [{genome}] Processed {idx} / {total}...")

        tier = var.get("tier", "Discarded")

        # Skip discards and synonymous variants
        if "Discarded" in tier or var["is_synonymous"]:
            continue

        aa_pos, mut_aa, nt_pos = parse_variant_coordinates(var)
        if not aa_pos:
            continue

        locus = var["locus"].split("/")[0]  # Handle overlapping mtDNA genes

        # Load AlignmentParser lazily (cached across variants in the same locus)
        if locus not in loaded_alignments:
            aa_path, nt_path = get_alignment_paths(locus, genome)
            if not aa_path.exists() or not nt_path.exists():
                loaded_alignments[locus] = None
            else:
                loaded_alignments[locus] = AlignmentParser(aa_path, nt_path, genome)

        parser = loaded_alignments[locus]
        if not parser:
            continue

        # === Ref Allele Verification ===
        is_ref_match = check_ref_allele(parser, nt_pos, var["ref_nt"])
        if not is_ref_match:
            ref_mismatch += 1
            var["ref_allele_match"] = False
            var["mismatch_reason"] = "Genomic reference mismatch"
        else:
            var["ref_allele_match"] = True
            var["mismatch_reason"] = None

        # Calculate mutant codon and test compensation
        mut_codon = parser.extract_mutant_codon(nt_pos, var["alt_nt"])
        if not mut_codon:
            continue

        # === Triangle of Truth Validation (Transcript Check) ===
        if genome == "nucDNA":
            try:
                translated_aa = str(Seq(mut_codon).translate())
            except Exception:
                translated_aa = "ERR"

            expected_aa = var.get("alt_aa", mut_aa)

            if translated_aa != expected_aa:
                var["ref_allele_match"] = False
                var["mismatch_reason"] = (
                    f"Transcript mismatch (TOGA uses {parser.transcript_id})"
                )

                # If it passed the pure base check but failed translation, log the mismatch
                if is_ref_match:
                    ref_mismatch += 1

                # Assign default negative values to preserve JSON schema
                var["cdar_aa"] = False
                var["cdar_nt"] = False
                var["cdar_aa_species"] = []
                var["cdar_nt_species"] = []
                var["compensating_species_count"] = 0

                enriched.append(var)
                global_stats["Total"] += 1
                tier_stats[tier]["Total"] += 1
                continue

        comp_results = parser.check_compensation(aa_pos, mut_aa, nt_pos, mut_codon)

        # Update variant object
        var["cdar_aa"] = comp_results["aa_cdar"]
        var["cdar_nt"] = comp_results["nt_cdar"]
        var["cdar_aa_species"] = comp_results["aa_species"]
        var["cdar_nt_species"] = comp_results["nt_species"]
        var["compensating_species_count"] = len(comp_results["aa_species"])

        enriched.append(var)

        # Update stats
        global_stats["Total"] += 1
        tier_stats[tier]["Total"] += 1
        if comp_results["aa_cdar"]:
            global_stats["aa_cDAR"] += 1
            tier_stats[tier]["aa_cDAR"] += 1
        if comp_results["nt_cdar"]:
            global_stats["nt_cDAR"] += 1
            tier_stats[tier]["nt_cDAR"] += 1

    return enriched, global_stats, dict(tier_stats), ref_mismatch


def print_genome_summary(
    genome: str, enriched: list, global_stats: dict, tier_stats: dict, ref_mismatch: int
):
    print(f"\n{genome}:")
    print(f"  Analyzed Variants  : {global_stats['Total']}")
    print(f"  Ref Allele Mismatches: {ref_mismatch}")
    if global_stats["Total"] > 0:
        aa_rate = (global_stats["aa_cDAR"] / global_stats["Total"]) * 100
        nt_rate = (global_stats["nt_cDAR"] / global_stats["Total"]) * 100
        print(f"  AA-Level c-DARs    : {global_stats['aa_cDAR']} ({aa_rate:.1f}%)")
        print(f"  NT-Level c-DARs    : {global_stats['nt_cDAR']} ({nt_rate:.1f}%)")

    print(f"\n  {'Tier':<12} {'Total':>7} {'aa_cDAR':>9} {'nt_cDAR':>9}")
    print(f"  {'-'*40}")
    for tier in sorted(tier_stats.keys()):
        s = tier_stats[tier]
        aa_pct = f"({(s['aa_cDAR']/s['Total'])*100:.0f}%)" if s["Total"] else ""
        nt_pct = f"({(s['nt_cDAR']/s['Total'])*100:.0f}%)" if s["Total"] else ""
        print(
            f"  {tier:<12} {s['Total']:>7} {s['aa_cDAR']:>5} {aa_pct:<5} {s['nt_cDAR']:>5} {nt_pct}"
        )
    print("-" * 50)


def main():
    print("Initializing c-DAR Identification Engine...\n")

    with open(MT_CURATED, "r") as f:
        mt_variants = json.load(f)
    with open(NUC_CURATED, "r") as f:
        nuc_variants = json.load(f)

    print(
        f"Loaded {len(mt_variants)} mtDNA variants, {len(nuc_variants)} nucDNA variants."
    )

    # Alignment cache is shared so overlapping loci are only loaded once
    loaded_alignments = {}

    print("\nProcessing mtDNA variants...")
    mt_enriched, mt_global, mt_tier, mt_mismatch = process_variants(
        mt_variants, "mtDNA", loaded_alignments
    )

    print("\nProcessing nucDNA variants...")
    nuc_enriched, nuc_global, nuc_tier, nuc_mismatch = process_variants(
        nuc_variants, "nucDNA", loaded_alignments
    )

    # Save outputs
    with open(OUT_JSON_MT, "w", encoding="utf-8") as f:
        json.dump(mt_enriched, f, indent=2)
    with open(OUT_JSON_NUC, "w", encoding="utf-8") as f:
        json.dump(nuc_enriched, f, indent=2)

    print(f"\n{'='*50}")
    print("c-DAR DISCOVERY SUMMARY")
    print(f"{'='*50}")
    print_genome_summary("mtDNA", mt_enriched, mt_global, mt_tier, mt_mismatch)
    print_genome_summary("nucDNA", nuc_enriched, nuc_global, nuc_tier, nuc_mismatch)


if __name__ == "__main__":
    main()
