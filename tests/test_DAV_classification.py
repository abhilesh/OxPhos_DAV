"""
4-Layer c-DAR/u-DAR Classification Verification Diagnostic
============================================================
Verifies the output of 00_classify_DAV.py at every logical level
for BOTH mtDNA and nucDNA datasets:

  Layer 1 — Triangle of Truth        : codon injection + translation self-check
  Layer 2 — Species Claim Validation : direct alignment re-scan for each c-DAR call
  Layer 3 — NT ⊆ AA Consistency      : cdar_nt=True must imply cdar_aa=True
  Layer 4 — Ref Allele Mismatch Audit: informational mismatch rate report

Layers 1 and 2 use random sampling over multiple trials for robust coverage.
Layers 3 and 4 scan the full dataset.

Run from project root:
    python tests/test_DAV_classification.py
"""

import csv
import json
import random
import re
import sys
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT / "src"))

from utils.alignment_parser import AlignmentParser

# ==== Configuration ====
CURATED_DIR  = ROOT / "data" / "annotations" / "curated"
MT_JSON      = CURATED_DIR / "cdar_classifications_mtDNA.json"
NUC_JSON     = CURATED_DIR / "cdar_classifications_nucDNA.json"

TOGA_AA_DIR  = ROOT / "data" / "alignments" / "toga_hg38_aa"
TOGA_NT_DIR  = ROOT / "data" / "alignments" / "toga_hg38_codon"
MT_AA_DIR    = ROOT / "data" / "alignments" / "mtdna_aa"
MT_NT_DIR    = ROOT / "data" / "alignments" / "mtdna_codon"
MT_COORD_FILE = ROOT / "data" / "reference" / "mtdna_gene_coordinates.tsv"


_MT_ALIAS = {"MT-COX1": "MT-CO1", "MT-COX2": "MT-CO2", "MT-COX3": "MT-CO3", "MT-CYTB": "MT-CYB"}


def _load_mt_coords() -> dict:
    coords = {}
    if not MT_COORD_FILE.exists():
        return coords
    with open(MT_COORD_FILE) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            start, end = int(row["start"]), int(row["end"])
            entry = (min(start, end), max(start, end), row["strand"])
            coords[row["gene"]] = entry
            if row["gene"] in _MT_ALIAS:
                coords[_MT_ALIAS[row["gene"]]] = entry
    return coords


_MT_COORDS = _load_mt_coords()


def _genomic_to_cds(genomic_pos: int, locus: str) -> int | None:
    if locus not in _MT_COORDS:
        return None
    start, end, strand = _MT_COORDS[locus]
    return genomic_pos - start + 1 if strand == "+" else end - genomic_pos + 1

# Sampling: SAMPLE_SIZE variants drawn per trial, N_TRIALS independent draws
SAMPLE_SIZE = 75
N_TRIALS    = 5

# Known genuine NT⊄AA violations — Layer 3 fails only if this count INCREASES
LAYER3_KNOWN_VIOLATIONS = 0


# ==== Helpers ====

def parse_variant_coordinates(variant: dict) -> tuple:
    """Mirrors parse_variant_coordinates in 00_classify_DAV.py exactly.

    Returns (aa_pos, wt_aa, mut_aa, nt_pos).
    """
    aa_str = variant.get("aa_change", "")
    nc_str = variant.get("nc_change", "")

    aa_match = re.search(r"([a-zA-Z]+)(\d+)([a-zA-Z]+)", aa_str)
    if not aa_match:
        return None, None, None, None

    wt_aa  = aa_match.group(1)
    aa_pos = int(aa_match.group(2))
    mut_aa = aa_match.group(3)

    if variant.get("genome") == "nucDNA":
        nc_match = re.search(r"c\.(\d+)", nc_str)
        if nc_match:
            return aa_pos, wt_aa, mut_aa, int(nc_match.group(1))
    else:
        genomic_pos = variant.get("genomic_pos")
        locus = variant.get("locus", "").split("/")[0]
        if genomic_pos and locus:
            cds_pos = _genomic_to_cds(genomic_pos, locus)
            if cds_pos:
                return aa_pos, wt_aa, mut_aa, cds_pos

    return aa_pos, wt_aa, mut_aa, (aa_pos * 3) - 2


def aa_dir(genome: str) -> Path:
    return MT_AA_DIR if genome == "mtDNA" else TOGA_AA_DIR


def nt_dir(genome: str) -> Path:
    return MT_NT_DIR if genome == "mtDNA" else TOGA_NT_DIR


def load_parser(locus: str, genome: str, cache: dict):
    key = (genome, locus)
    if key not in cache:
        aa_path = aa_dir(genome) / f"{locus}_aa_alignment.fasta"
        nt_path = nt_dir(genome) / f"{locus}_codon_alignment.fasta"
        cache[key] = (
            AlignmentParser(aa_path, nt_path, genome)
            if aa_path.exists() and nt_path.exists()
            else None
        )
    return cache[key]


def _print_trial_bar(trial_results: list, n_eligible: int):
    for i, (passed, total) in enumerate(trial_results):
        bar = "#" * passed + "-" * (total - passed)
        print(f"    Trial {i+1:2d}: {passed:>4}/{total:<4}  [{bar}]")
    total_p = sum(p for p, _ in trial_results)
    total_t = sum(t for _, t in trial_results)
    rate    = (total_p / total_t * 100) if total_t else 0
    print(f"\n  Consolidated : {total_p}/{total_t} passed ({rate:.1f}%)  |  "
          f"Eligible pool: {n_eligible}")


# ==== Layer 1 — Triangle of Truth ====

def _build_layer1_eligible(variants: list, genome: str, cache: dict) -> list:
    eligible = []
    for var in variants:
        if not var.get("cdar_aa"):
            continue
        locus = var["locus"].split("/")[0]
        aa_pos, wt_aa, mut_aa, nt_pos = parse_variant_coordinates(var)
        if not aa_pos:
            continue
        parser = load_parser(locus, genome, cache)
        if not parser or nt_pos not in parser.nt_map:
            continue
        codon_start_bio = nt_pos - ((nt_pos - 1) % 3)
        if codon_start_bio not in parser.nt_map:
            continue
        eligible.append((var, mut_aa, nt_pos))
    return eligible


def _run_layer1_trial(sample: list, genome: str, cache: dict) -> tuple:
    passed = 0
    failures = {}
    for var, mut_aa, nt_pos in sample:
        locus  = var["locus"].split("/")[0]
        parser = load_parser(locus, genome, cache)

        codon_start_bio = nt_pos - ((nt_pos - 1) % 3)
        col_start       = parser.nt_map[codon_start_bio]
        wt_codon        = parser.nt_alignment[parser.ref_header][col_start: col_start + 3]
        alt_nt          = var.get("alt_nt", "")
        pos_in_codon    = (nt_pos - 1) % 3
        mut_codon       = wt_codon[:pos_in_codon] + alt_nt + wt_codon[pos_in_codon + 1:]

        _table = 2 if genome == "mtDNA" else 1
        try:
            translated = str(Seq(mut_codon).translate(table=_table))
        except Exception:
            translated = "ERR"

        alt_aa = var.get("alt_aa") or mut_aa
        if translated == alt_aa:
            passed += 1
        else:
            key = var.get("ann_id", "") + "|" + var.get("nc_change", "")
            failures[key] = (var, wt_codon, mut_codon, translated, alt_aa)

    return passed, len(sample), failures


def layer1_triangle_of_truth(variants: list, genome: str, cache: dict) -> bool:
    print(f"\n{'='*60}")
    print(f"LAYER 1: Triangle of Truth  [{genome}]")
    print(f"{'='*60}")
    print(f"Inject alt_nt into WT codon → translate → must equal alt_aa.")
    print(f"{N_TRIALS} trials × {SAMPLE_SIZE} random c-DARs per trial.\n")

    eligible = _build_layer1_eligible(variants, genome, cache)
    if not eligible:
        print("  No eligible c-DARs found — skipping.")
        return True

    trial_results = []
    all_failures  = {}

    for _ in range(N_TRIALS):
        sample = random.sample(eligible, min(SAMPLE_SIZE, len(eligible)))
        passed, total, failures = _run_layer1_trial(sample, genome, cache)
        trial_results.append((passed, total))
        all_failures.update(failures)

    _print_trial_bar(trial_results, len(eligible))

    if all_failures:
        print(f"\n  Unique failing variants across all trials ({len(all_failures)}):")
        for key, (var, wt, mut, trans, expected) in list(all_failures.items())[:15]:
            locus = var["locus"].split("/")[0]
            print(
                f"    {locus:10s} {var.get('nc_change','?'):15s} ({var.get('aa_change','?'):8s})"
                f"  WT:{wt} → mut:{mut} → {trans}  (expected {expected})"
            )

    result = len(all_failures) == 0
    print(f"\nLayer 1 result: {'PASS' if result else 'FAIL'}")
    return result


# ==== Layer 2 — Species Claim Validation ====

def _build_layer2a_eligible(variants: list, genome: str, cache: dict) -> list:
    eligible = []
    for var in variants:
        if not var.get("cdar_aa") or not var.get("cdar_aa_species"):
            continue
        locus = var["locus"].split("/")[0]
        aa_pos, wt_aa, mut_aa, _ = parse_variant_coordinates(var)
        if not aa_pos:
            continue
        parser = load_parser(locus, genome, cache)
        if not parser or aa_pos not in parser.aa_map:
            continue
        eligible.append((var, aa_pos, mut_aa))
    return eligible


def _build_layer2b_eligible(variants: list, genome: str, cache: dict) -> list:
    eligible = []
    for var in variants:
        if var.get("cdar_aa") or not var.get("ref_allele_match"):
            continue
        locus = var["locus"].split("/")[0]
        aa_pos, wt_aa, mut_aa, _ = parse_variant_coordinates(var)
        if not aa_pos:
            continue
        parser = load_parser(locus, genome, cache)
        if not parser or aa_pos not in parser.aa_map:
            continue
        eligible.append((var, aa_pos, mut_aa))
    return eligible


def _run_layer2a_trial(sample: list, genome: str, cache: dict) -> tuple:
    passed = 0
    failures = {}
    for var, aa_pos, mut_aa in sample:
        locus  = var["locus"].split("/")[0]
        parser = load_parser(locus, genome, cache)
        aa_col = parser.aa_map[aa_pos]
        variant_ok = True

        for species in var["cdar_aa_species"]:
            sp_seq = parser.aa_alignment.get(species)
            if sp_seq is None:
                continue
            actual_aa = sp_seq[aa_col]
            if actual_aa != mut_aa:
                variant_ok = False
                key = var.get("ann_id", "") + "|" + species
                failures[key] = (var, species, actual_aa, mut_aa, aa_col)

        if variant_ok:
            passed += 1

    return passed, len(sample), failures


def _run_layer2b_trial(sample: list, genome: str, cache: dict) -> tuple:
    passed = 0
    failures = {}
    MASKED = {"X", "-", "!", "*"}

    for var, aa_pos, mut_aa in sample:
        locus  = var["locus"].split("/")[0]
        parser = load_parser(locus, genome, cache)
        aa_col = parser.aa_map[aa_pos]

        imposters = [
            sp for sp, seq in parser.aa_alignment.items()
            if sp != parser.ref_header
            and seq[aa_col] not in MASKED
            and seq[aa_col] == mut_aa
        ]
        if imposters:
            key = var.get("ann_id", "") + "|" + var.get("nc_change", "")
            failures[key] = (var, imposters)
        else:
            passed += 1

    return passed, len(sample), failures


def layer2_species_claim_validation(variants: list, genome: str, cache: dict) -> bool:
    print(f"\n{'='*60}")
    print(f"LAYER 2: Species Claim Validation  [{genome}]")
    print(f"{'='*60}")
    print(f"{N_TRIALS} trials × {SAMPLE_SIZE} random variants per part.\n")

    # --- Part A: c-DAR species re-scan ---
    print(f"Part A — Verify cdar_aa_species claims (c-DARs):")
    eligible_a = _build_layer2a_eligible(variants, genome, cache)

    a_trial_results = []
    a_all_failures  = {}
    for _ in range(N_TRIALS):
        sample  = random.sample(eligible_a, min(SAMPLE_SIZE, len(eligible_a)))
        passed, total, failures = _run_layer2a_trial(sample, genome, cache)
        a_trial_results.append((passed, total))
        a_all_failures.update(failures)

    _print_trial_bar(a_trial_results, len(eligible_a))

    if a_all_failures:
        print(f"\n  Failing species claims ({len(a_all_failures)}):")
        for key, (var, sp, actual, expected, col) in list(a_all_failures.items())[:10]:
            print(f"    {var['locus']:10s} {var.get('aa_change','?'):8s}  "
                  f"species '{sp}': has '{actual}' at col {col}, expected '{expected}'")

    # --- Part B: u-DAR negative scan ---
    print(f"\nPart B — Negative scan (u-DARs, ref_allele_match=True):")
    eligible_b = _build_layer2b_eligible(variants, genome, cache)

    b_trial_results = []
    b_all_failures  = {}
    for _ in range(N_TRIALS):
        sample  = random.sample(eligible_b, min(SAMPLE_SIZE, len(eligible_b)))
        passed, total, failures = _run_layer2b_trial(sample, genome, cache)
        b_trial_results.append((passed, total))
        b_all_failures.update(failures)

    _print_trial_bar(b_trial_results, len(eligible_b))

    if b_all_failures:
        print(f"\n  u-DARs with unexpected compensating species ({len(b_all_failures)}):")
        for key, (var, imposters) in list(b_all_failures.items())[:10]:
            print(f"    {var['locus']:10s} {var.get('aa_change','?'):8s}  "
                  f"found in: {imposters[:3]}{'...' if len(imposters) > 3 else ''}")

    result = not a_all_failures and not b_all_failures
    print(f"\nLayer 2 result: {'PASS' if result else 'FAIL'}")
    return result


# ==== Layer 3 — NT ⊆ AA Logical Consistency ====

def layer3_nt_subset_aa_consistency(variants: list, genome: str) -> bool:
    print(f"\n{'='*60}")
    print(f"LAYER 3: NT ⊆ AA Logical Consistency  [{genome}]")
    print(f"{'='*60}")
    print("Rule: cdar_nt=True must imply cdar_aa=True for every variant.\n")

    expected_violations = []
    genuine_violations  = []

    for var in variants:
        if var.get("cdar_nt") and not var.get("cdar_aa"):
            bucket = expected_violations if not var.get("ref_allele_match") else genuine_violations
            bucket.append(var)

    print(f"  Expected violations (ref_allele_match=False, spurious codon hit): {len(expected_violations)}")
    print(f"  Genuine violations  (ref_allele_match=True,  real logic error)  : {len(genuine_violations)}")

    if genuine_violations:
        print(f"\n  Genuine violation detail (first 10):")
        for v in genuine_violations[:10]:
            print(
                f"    {v['locus']:10s} {v.get('nc_change','?'):15s} ({v.get('aa_change','?'):8s})"
                f"  nt_species={v.get('cdar_nt_species', [])[:3]}"
            )

    if len(genuine_violations) <= LAYER3_KNOWN_VIOLATIONS:
        print(
            f"\n  [{len(genuine_violations)}/{LAYER3_KNOWN_VIOLATIONS} known violations] "
            f"Layer 3 PASS — no regression detected."
        )
        return True
    else:
        diff = len(genuine_violations) - LAYER3_KNOWN_VIOLATIONS
        print(f"\n  REGRESSION: {diff} new genuine violations above baseline of {LAYER3_KNOWN_VIOLATIONS}.")
        return False


# ==== Layer 4 — Ref Allele Mismatch Audit ====

def layer4_ref_allele_mismatch_audit(variants: list, genome: str) -> bool:
    print(f"\n{'='*60}")
    print(f"LAYER 4: Ref Allele Mismatch Audit  [{genome}]  [informational]")
    print(f"{'='*60}")

    total      = len(variants)
    mismatches = [v for v in variants if not v.get("ref_allele_match")]
    n_mismatch = len(mismatches)
    pct        = (n_mismatch / total * 100) if total else 0

    print(f"  Total variants       : {total}")
    print(f"  Ref allele mismatches: {n_mismatch} ({pct:.1f}%)")

    if pct > 5:
        print(
            f"\n  WARNING: Mismatch rate {pct:.1f}% exceeds 5% threshold.\n"
            f"  Check ref_nt/alt_nt allele orientation vs CDS alignment strand."
        )

    print(f"\n  First 5 mismatches:")
    print(f"  {'Locus':10s} {'nc_change':20s} {'ref_nt':6s} {'alt_nt':6s} {'tier'}")
    print(f"  {'-'*57}")
    for v in mismatches[:5]:
        print(
            f"  {v['locus']:10s} {v.get('nc_change','?'):20s} "
            f"{v.get('ref_nt','?'):6s} {v.get('alt_nt','?'):6s} {v.get('tier','?')}"
        )

    locus_counts = Counter(v["locus"] for v in mismatches)
    print(f"\n  Top 5 loci by mismatch count:")
    for locus, count in locus_counts.most_common(5):
        print(f"    {locus:12s}: {count:5d} ({count/total*100:.1f}% of all variants)")

    note = "NOTE" if pct > 5 else "PASS"
    print(f"\nLayer 4 result: {note}")
    return True


# ==== Main ====

def run_genome(json_path: Path, genome: str, cache: dict) -> list:
    """Run all 4 layers for one genome. Returns list of (name, passed) tuples."""
    if not json_path.exists():
        print(f"\nERROR: {json_path.name} not found — run 00_classify_DAV.py first.")
        return [(f"Layer {i}: {genome}", False) for i in range(1, 5)]

    with open(json_path) as f:
        variants = json.load(f)
    print(f"\nLoaded {len(variants)} {genome} variants from {json_path.name}")

    results = []
    results.append((f"Layer 1: Triangle of Truth           [{genome}]",
                     layer1_triangle_of_truth(variants, genome, cache)))
    results.append((f"Layer 2: Species Claim Validation    [{genome}]",
                     layer2_species_claim_validation(variants, genome, cache)))
    results.append((f"Layer 3: NT ⊆ AA Consistency         [{genome}]",
                     layer3_nt_subset_aa_consistency(variants, genome)))
    results.append((f"Layer 4: Ref Allele Mismatch Audit   [{genome}]",
                     layer4_ref_allele_mismatch_audit(variants, genome)))
    return results


def main():
    print("=" * 60)
    print("DAV Classification Verification — 4-Layer Diagnostic")
    print(f"Sampling: {N_TRIALS} trials × {SAMPLE_SIZE} variants")
    print("=" * 60)

    cache = {}  # shared across both genomes — loci loaded once

    all_results = []
    all_results += run_genome(MT_JSON,  "mtDNA",  cache)
    all_results += run_genome(NUC_JSON, "nucDNA", cache)

    print(f"\n{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    passed = sum(1 for _, ok in all_results if ok)
    for name, ok in all_results:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}]  {name}")
    print(f"\n{passed} / {len(all_results)} layers passed.")


if __name__ == "__main__":
    main()
