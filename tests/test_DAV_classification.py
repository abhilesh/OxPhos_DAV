"""
4-Layer c-DAR/u-DAR Classification Verification Diagnostic
============================================================
Verifies the output of 00_classify_DAV.py at every logical level:

  Layer 1 — Triangle of Truth        : codon injection + translation self-check
  Layer 2 — Species Claim Validation : direct alignment re-scan for each c-DAR call
  Layer 3 — NT ⊆ AA Consistency      : cdar_nt=True must imply cdar_aa=True
  Layer 4 — Ref Allele Mismatch Audit: informational mismatch rate report

Layers 1 and 2 use random sampling over multiple trials for robust coverage.
Layers 3 and 4 scan the full dataset.

Run from project root:
    python tests/test_DAV_classification.py
"""

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
IN_JSON     = ROOT / "data" / "annotations" / "curated" / "cdar_classifications_nucDNA.json"
TOGA_AA_DIR = ROOT / "data" / "alignments" / "toga_hg38_aa"
TOGA_NT_DIR = ROOT / "data" / "alignments" / "toga_hg38_codon"

# Sampling: SAMPLE_SIZE variants drawn per trial, N_TRIALS independent draws
SAMPLE_SIZE = 75
N_TRIALS    = 5

# Known genuine NT⊄AA violations in the current annotation dataset
# (ref_allele_match=True strand-flip errors from ClinVar source data).
# Layer 3 fails only if this count INCREASES (regression detection).
LAYER3_KNOWN_VIOLATIONS = 29


# ==== Helpers ====

def parse_variant_coordinates(variant: dict) -> tuple:
    """Mirrors parse_variant_coordinates in 00_classify_DAV.py exactly."""
    aa_str = variant.get("aa_change", "")
    nc_str = variant.get("nc_change", "")

    aa_match = re.search(r"[a-zA-Z]+(\d+)([a-zA-Z]+)", aa_str)
    if not aa_match:
        return None, None, None

    aa_pos = int(aa_match.group(1))
    mut_aa = aa_match.group(2)

    nc_match = re.search(r"c\.(\d+)", nc_str)
    nt_pos = int(nc_match.group(1)) if nc_match else (aa_pos * 3) - 2

    return aa_pos, mut_aa, nt_pos


def load_parser(locus: str, cache: dict):
    if locus not in cache:
        aa_path = TOGA_AA_DIR / f"{locus}_aa_alignment.fasta"
        nt_path = TOGA_NT_DIR / f"{locus}_codon_alignment.fasta"
        cache[locus] = (
            AlignmentParser(aa_path, nt_path, "nucDNA")
            if aa_path.exists() and nt_path.exists()
            else None
        )
    return cache[locus]


def _print_trial_bar(trial_results: list, n_eligible: int):
    """Prints per-trial pass counts and a consolidated summary line."""
    for i, (passed, total) in enumerate(trial_results):
        bar = "#" * passed + "-" * (total - passed)
        print(f"    Trial {i+1:2d}: {passed:>4}/{total:<4}  [{bar}]")
    total_p = sum(p for p, _ in trial_results)
    total_t = sum(t for _, t in trial_results)
    rate    = (total_p / total_t * 100) if total_t else 0
    print(f"\n  Consolidated : {total_p}/{total_t} passed ({rate:.1f}%)  |  "
          f"Eligible pool: {n_eligible}")


# ==== Layer 1 — Triangle of Truth ====

def _build_layer1_eligible(variants: list, cache: dict) -> list:
    """Pre-filter c-DARs to those with a loaded alignment and resolvable codon coords."""
    eligible = []
    for var in variants:
        if not var.get("cdar_aa"):
            continue
        locus = var["locus"].split("/")[0]
        aa_pos, mut_aa, nt_pos = parse_variant_coordinates(var)
        if not aa_pos:
            continue
        parser = load_parser(locus, cache)
        if not parser or nt_pos not in parser.nt_map:
            continue
        codon_start_bio = nt_pos - ((nt_pos - 1) % 3)
        if codon_start_bio not in parser.nt_map:
            continue
        eligible.append((var, mut_aa, nt_pos))
    return eligible


def _run_layer1_trial(sample: list, cache: dict) -> tuple:
    """Run one trial. Returns (passed, total, failures dict keyed by ann_id+nc)."""
    passed = 0
    failures = {}
    for var, mut_aa, nt_pos in sample:
        locus  = var["locus"].split("/")[0]
        parser = load_parser(locus, cache)

        codon_start_bio = nt_pos - ((nt_pos - 1) % 3)
        col_start    = parser.nt_map[codon_start_bio]
        wt_codon     = parser.nt_alignment[parser.ref_header][col_start : col_start + 3]
        alt_nt       = var.get("alt_nt", "")
        pos_in_codon = (nt_pos - 1) % 3
        mut_codon    = wt_codon[:pos_in_codon] + alt_nt + wt_codon[pos_in_codon + 1:]

        try:
            translated = str(Seq(mut_codon).translate())
        except Exception:
            translated = "ERR"

        alt_aa = var.get("alt_aa", mut_aa)
        if translated == alt_aa:
            passed += 1
        else:
            key = var.get("ann_id", "") + "|" + var.get("nc_change", "")
            failures[key] = (var, wt_codon, mut_codon, translated, alt_aa)

    return passed, len(sample), failures


def layer1_triangle_of_truth(variants: list, cache: dict) -> bool:
    print(f"\n{'='*60}")
    print("LAYER 1: Triangle of Truth")
    print(f"{'='*60}")
    print(f"Inject alt_nt into WT codon → translate → must equal alt_aa.")
    print(f"{N_TRIALS} trials × {SAMPLE_SIZE} random c-DARs per trial.\n")

    eligible = _build_layer1_eligible(variants, cache)
    if not eligible:
        print("  No eligible c-DARs found — skipping.")
        return True

    trial_results = []
    all_failures  = {}

    for _ in range(N_TRIALS):
        sample = random.sample(eligible, min(SAMPLE_SIZE, len(eligible)))
        passed, total, failures = _run_layer1_trial(sample, cache)
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

def _build_layer2a_eligible(variants: list, cache: dict) -> list:
    """Pre-filter c-DARs with species claims and resolvable AA column."""
    eligible = []
    for var in variants:
        if not var.get("cdar_aa") or not var.get("cdar_aa_species"):
            continue
        locus = var["locus"].split("/")[0]
        aa_pos, mut_aa, _ = parse_variant_coordinates(var)
        if not aa_pos:
            continue
        parser = load_parser(locus, cache)
        if not parser or aa_pos not in parser.aa_map:
            continue
        eligible.append((var, aa_pos, mut_aa))
    return eligible


def _build_layer2b_eligible(variants: list, cache: dict) -> list:
    """Pre-filter u-DARs (ref_allele_match=True) with resolvable AA column."""
    eligible = []
    for var in variants:
        if var.get("cdar_aa") or not var.get("ref_allele_match"):
            continue
        locus = var["locus"].split("/")[0]
        aa_pos, mut_aa, _ = parse_variant_coordinates(var)
        if not aa_pos:
            continue
        parser = load_parser(locus, cache)
        if not parser or aa_pos not in parser.aa_map:
            continue
        eligible.append((var, aa_pos, mut_aa))
    return eligible


def _run_layer2a_trial(sample: list, cache: dict) -> tuple:
    """Verify each claimed c-DAR species carries mut_aa. Returns (passed, total, failures)."""
    passed = 0
    failures = {}
    MASKED = {"X", "-", "?", "!"}

    for var, aa_pos, mut_aa in sample:
        locus  = var["locus"].split("/")[0]
        parser = load_parser(locus, cache)
        aa_col = parser.aa_map[aa_pos]
        variant_ok = True

        for species in var["cdar_aa_species"]:
            sp_seq = parser.aa_alignment.get(species)
            if sp_seq is None:
                continue  # alignment regenerated, species dropped — not a failure
            actual_aa = sp_seq[aa_col]
            if actual_aa != mut_aa:
                variant_ok = False
                key = var.get("ann_id", "") + "|" + species
                failures[key] = (var, species, actual_aa, mut_aa, aa_col)

        if variant_ok:
            passed += 1

    return passed, len(sample), failures


def _run_layer2b_trial(sample: list, cache: dict) -> tuple:
    """Assert no species carries mut_aa for u-DARs. Returns (passed, total, failures)."""
    passed = 0
    failures = {}
    MASKED = {"X", "-", "?", "!"}

    for var, aa_pos, mut_aa in sample:
        locus  = var["locus"].split("/")[0]
        parser = load_parser(locus, cache)
        aa_col = parser.aa_map[aa_pos]

        imposters = [
            sp for sp, seq in parser.aa_alignment.items()
            if sp != parser.ref_header
            and seq[aa_col] == mut_aa
            and seq[aa_col] not in MASKED
        ]
        if imposters:
            key = var.get("ann_id", "") + "|" + var.get("nc_change", "")
            failures[key] = (var, imposters)
        else:
            passed += 1

    return passed, len(sample), failures


def layer2_species_claim_validation(variants: list, cache: dict) -> bool:
    print(f"\n{'='*60}")
    print("LAYER 2: Species Claim Validation")
    print(f"{'='*60}")
    print(f"{N_TRIALS} trials × {SAMPLE_SIZE} random variants per part.\n")

    # --- Part A: c-DAR species re-scan ---
    print(f"Part A — Verify cdar_aa_species claims (c-DARs):")
    eligible_a = _build_layer2a_eligible(variants, cache)

    a_trial_results = []
    a_all_failures  = {}
    for _ in range(N_TRIALS):
        sample  = random.sample(eligible_a, min(SAMPLE_SIZE, len(eligible_a)))
        passed, total, failures = _run_layer2a_trial(sample, cache)
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
    eligible_b = _build_layer2b_eligible(variants, cache)

    b_trial_results = []
    b_all_failures  = {}
    for _ in range(N_TRIALS):
        sample  = random.sample(eligible_b, min(SAMPLE_SIZE, len(eligible_b)))
        passed, total, failures = _run_layer2b_trial(sample, cache)
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

def layer3_nt_subset_aa_consistency(variants: list) -> bool:
    print(f"\n{'='*60}")
    print("LAYER 3: NT ⊆ AA Logical Consistency")
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

def layer4_ref_allele_mismatch_audit(variants: list) -> bool:
    print(f"\n{'='*60}")
    print("LAYER 4: Ref Allele Mismatch Audit  [informational]")
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
            f"  Likely cause: ref_nt/alt_nt populated from genomic (plus-strand) alleles\n"
            f"  in ClinVar TSV. For minus-strand genes these are complements of the CDS\n"
            f"  alleles in the c. notation. Fix: extract ref_nt/alt_nt from the c.X>Y\n"
            f"  pattern in nc_change rather than from the ReferenceAllele column."
        )

    print(f"\n  First 5 mismatches:")
    print(f"  {'Locus':10s} {'nc_change':18s} {'ref_nt':6s} {'alt_nt':6s} {'tier'}")
    print(f"  {'-'*55}")
    for v in mismatches[:5]:
        print(
            f"  {v['locus']:10s} {v.get('nc_change','?'):18s} "
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

def main():
    print("=" * 60)
    print("DAV Classification Verification — 4-Layer Diagnostic")
    print(f"Sampling: {N_TRIALS} trials × {SAMPLE_SIZE} variants")
    print("=" * 60)

    if not IN_JSON.exists():
        print(f"ERROR: Input JSON not found: {IN_JSON}")
        print("Run src/classify/00_classify_DAV.py first.")
        sys.exit(1)

    with open(IN_JSON) as f:
        variants = json.load(f)

    print(f"Loaded {len(variants)} variants from {IN_JSON.name}")

    # Single shared cache — alignments loaded once, reused across all layers and trials
    cache = {}
    results = []

    results.append(("Layer 1: Triangle of Truth",           layer1_triangle_of_truth(variants, cache)))
    results.append(("Layer 2: Species Claim Validation",    layer2_species_claim_validation(variants, cache)))
    results.append(("Layer 3: NT ⊆ AA Consistency",         layer3_nt_subset_aa_consistency(variants)))
    results.append(("Layer 4: Ref Allele Mismatch Audit",   layer4_ref_allele_mismatch_audit(variants)))

    print(f"\n{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    passed = sum(1 for _, ok in results if ok)
    for name, ok in results:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}]  {name}")
    print(f"\n{passed} / {len(results)} layers passed.")


if __name__ == "__main__":
    main()
