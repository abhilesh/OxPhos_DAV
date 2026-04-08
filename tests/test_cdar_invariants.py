"""
Data-driven invariant tests for cDAV classification output.

Validates the biological and logical constraints that must hold over the full
cdav_classifications_*.json datasets produced by 00_classify_DAV.py:

  INV-1  NT ⊆ AA (subset):       is_cdav_nucleotide=True  →  is_cdav_amino_acid=True
  INV-2  Count field consistency: n_species_with_disease_allele == len(lineages_with_disease_allele)
  INV-3  Global count:            total AA cDAVs >= total NT cDAVs  (per genome, per tier)
  INV-4  No empty species lists on positive AA flag:
             is_cdav_amino_acid=True  →  len(lineages_with_disease_allele) > 0
  INV-5  Discarded/synonymous variants must not appear in cDAV output

Failures are reported per-variant so every broken record is visible at once.

Run from project root:
    pytest tests/test_cdar_invariants.py -v
"""

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CURATED = ROOT / "data" / "annotations" / "curated"
MT_JSON  = CURATED / "cdav_classifications_mtDNA.json"
NUC_JSON = CURATED / "cdav_classifications_nucDNA.json"


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _load(path: Path) -> list:
    if not path.exists():
        pytest.skip(f"Classification output not found: {path.name} — run 00_classify_DAV.py first")
    with open(path) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def mt_variants():
    return _load(MT_JSON)


@pytest.fixture(scope="module")
def nuc_variants():
    return _load(NUC_JSON)


# ── Helper ────────────────────────────────────────────────────────────────────

def _label(v: dict) -> str:
    return f"{v.get('locus','?')} | {v.get('nc_change','?')} | {v.get('aa_change','?')}"


# ── INV-1: NT ⊆ AA (subset rule) ─────────────────────────────────────────────

class TestNtSubsetAa:

    def test_nt_implies_aa_mtdna(self, mt_variants):
        """Every mtDNA variant with is_cdav_nucleotide=True must also have is_cdav_amino_acid=True."""
        violations = [
            _label(v) for v in mt_variants
            if v.get("is_cdav_nucleotide") and not v.get("is_cdav_amino_acid")
        ]
        assert not violations, (
            f"INV-1 violated for {len(violations)} mtDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_nt_implies_aa_nucdna(self, nuc_variants):
        """Every nucDNA variant with is_cdav_nucleotide=True must also have is_cdav_amino_acid=True."""
        violations = [
            _label(v) for v in nuc_variants
            if v.get("is_cdav_nucleotide") and not v.get("is_cdav_amino_acid")
        ]
        assert not violations, (
            f"INV-1 violated for {len(violations)} nucDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )


# ── INV-2: Count field consistency ───────────────────────────────────────────

class TestSpeciesCountField:

    def _check(self, variants: list):
        violations = []
        for v in variants:
            stored  = v.get("n_species_with_disease_allele", -1)
            derived = len(v.get("lineages_with_disease_allele") or [])
            if stored != derived:
                violations.append(
                    f"{_label(v)} — stored={stored}, len(lineages)={derived}"
                )
        return violations

    def test_count_field_consistent_mtdna(self, mt_variants):
        violations = self._check(mt_variants)
        assert not violations, (
            f"INV-2 violated for {len(violations)} mtDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_count_field_consistent_nucdna(self, nuc_variants):
        violations = self._check(nuc_variants)
        assert not violations, (
            f"INV-2 violated for {len(violations)} nucDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )


# ── INV-3: Global and per-tier count monotonicity ─────────────────────────────

class TestGlobalCountMonotonicity:

    def _tier_counts(self, variants: list) -> dict:
        """Returns {tier: {"aa": int, "nt": int}} plus a "TOTAL" key."""
        counts = defaultdict(lambda: {"aa": 0, "nt": 0})
        for v in variants:
            tier = v.get("tier", "Unknown")
            if v.get("is_cdav_amino_acid"):
                counts[tier]["aa"] += 1
                counts["TOTAL"]["aa"] += 1
            if v.get("is_cdav_nucleotide"):
                counts[tier]["nt"] += 1
                counts["TOTAL"]["nt"] += 1
        return dict(counts)

    def test_global_aa_gte_nt_mtdna(self, mt_variants):
        c = self._tier_counts(mt_variants)
        total = c.get("TOTAL", {"aa": 0, "nt": 0})
        assert total["aa"] >= total["nt"], (
            f"Global mtDNA: aa_cDAVs={total['aa']} < nt_cDAVs={total['nt']}"
        )

    def test_global_aa_gte_nt_nucdna(self, nuc_variants):
        c = self._tier_counts(nuc_variants)
        total = c.get("TOTAL", {"aa": 0, "nt": 0})
        assert total["aa"] >= total["nt"], (
            f"Global nucDNA: aa_cDAVs={total['aa']} < nt_cDAVs={total['nt']}"
        )

    def test_per_tier_aa_gte_nt_mtdna(self, mt_variants):
        counts = self._tier_counts(mt_variants)
        violations = [
            f"{tier}: aa={c['aa']} < nt={c['nt']}"
            for tier, c in counts.items()
            if tier != "TOTAL" and c["nt"] > c["aa"]
        ]
        assert not violations, (
            "INV-3 per-tier violations (mtDNA):\n"
            + "\n".join(f"  {x}" for x in violations)
        )

    def test_per_tier_aa_gte_nt_nucdna(self, nuc_variants):
        counts = self._tier_counts(nuc_variants)
        violations = [
            f"{tier}: aa={c['aa']} < nt={c['nt']}"
            for tier, c in counts.items()
            if tier != "TOTAL" and c["nt"] > c["aa"]
        ]
        assert not violations, (
            "INV-3 per-tier violations (nucDNA):\n"
            + "\n".join(f"  {x}" for x in violations)
        )


# ── INV-4: Non-empty species lists on positive AA flag ───────────────────────

class TestNonEmptySpeciesOnPositiveFlag:

    def test_aa_flag_has_species_mtdna(self, mt_variants):
        violations = [
            _label(v) for v in mt_variants
            if v.get("is_cdav_amino_acid") and not v.get("lineages_with_disease_allele")
        ]
        assert not violations, (
            f"INV-4: {len(violations)} mtDNA variants flagged is_cdav_amino_acid=True "
            f"but have empty lineages_with_disease_allele:\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_aa_flag_has_species_nucdna(self, nuc_variants):
        violations = [
            _label(v) for v in nuc_variants
            if v.get("is_cdav_amino_acid") and not v.get("lineages_with_disease_allele")
        ]
        assert not violations, (
            f"INV-4: {len(violations)} nucDNA variants flagged is_cdav_amino_acid=True "
            f"but have empty lineages_with_disease_allele:\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_nt_flag_implies_nonempty_lineages_mtdna(self, mt_variants):
        """is_cdav_nucleotide=True implies is_cdav_amino_acid=True (INV-1), so lineages must be non-empty."""
        violations = [
            _label(v) for v in mt_variants
            if v.get("is_cdav_nucleotide") and not v.get("lineages_with_disease_allele")
        ]
        assert not violations, (
            f"INV-4 (nt): {len(violations)} mtDNA variants flagged is_cdav_nucleotide=True "
            f"but have empty lineages_with_disease_allele:\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_nt_flag_implies_nonempty_lineages_nucdna(self, nuc_variants):
        violations = [
            _label(v) for v in nuc_variants
            if v.get("is_cdav_nucleotide") and not v.get("lineages_with_disease_allele")
        ]
        assert not violations, (
            f"INV-4 (nt): {len(violations)} nucDNA variants flagged is_cdav_nucleotide=True "
            f"but have empty lineages_with_disease_allele:\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )


# ── INV-5: Discarded / synonymous variants absent from output ─────────────────

class TestDiscardedAndSynonymousAbsent:

    def test_no_discarded_in_mtdna_output(self, mt_variants):
        discarded = [_label(v) for v in mt_variants if v.get("tier") == "Discarded"]
        assert not discarded, (
            f"INV-5: {len(discarded)} Discarded mtDNA variants appear in cDAV output:\n"
            + "\n".join(f"  {x}" for x in discarded[:10])
        )

    def test_no_discarded_in_nucdna_output(self, nuc_variants):
        discarded = [_label(v) for v in nuc_variants if v.get("tier") == "Discarded"]
        assert not discarded, (
            f"INV-5: {len(discarded)} Discarded nucDNA variants appear in cDAV output:\n"
            + "\n".join(f"  {x}" for x in discarded[:10])
        )

    def test_no_synonymous_in_mtdna_output(self, mt_variants):
        synonymous = [_label(v) for v in mt_variants if v.get("is_synonymous")]
        assert not synonymous, (
            f"INV-5: {len(synonymous)} synonymous mtDNA variants appear in cDAV output:\n"
            + "\n".join(f"  {x}" for x in synonymous[:10])
        )

    def test_no_synonymous_in_nucdna_output(self, nuc_variants):
        synonymous = [_label(v) for v in nuc_variants if v.get("is_synonymous")]
        assert not synonymous, (
            f"INV-5: {len(synonymous)} synonymous nucDNA variants appear in cDAV output:\n"
            + "\n".join(f"  {x}" for x in synonymous[:10])
        )


# ── Summary report (printed even on pass) ────────────────────────────────────

def test_summary_report(mt_variants, nuc_variants):
    """Non-failing summary printed alongside results for quick orientation."""
    for label, variants in [("mtDNA", mt_variants), ("nucDNA", nuc_variants)]:
        total   = len(variants)
        n_aa    = sum(1 for v in variants if v.get("is_cdav_amino_acid"))
        n_nt    = sum(1 for v in variants if v.get("is_cdav_nucleotide"))
        by_tier = defaultdict(lambda: {"total": 0, "aa": 0, "nt": 0})
        for v in variants:
            t = v.get("tier", "?")
            by_tier[t]["total"] += 1
            if v.get("is_cdav_amino_acid"): by_tier[t]["aa"] += 1
            if v.get("is_cdav_nucleotide"): by_tier[t]["nt"] += 1

        print(f"\n{label} ({total} variants in cDAV output):")
        print(f"  AA cDAVs  : {n_aa}  ({100*n_aa/total:.1f}%)")
        print(f"  NT cDAVs  : {n_nt}  ({100*n_nt/total:.1f}%)")
        print(f"  NT/AA ratio: {n_nt/n_aa*100:.1f}%" if n_aa else "  NT/AA ratio: N/A")
        print(f"  {'Tier':<12} {'Total':>7} {'aa':>7} {'nt':>7} {'nt/aa %':>9}")
        for tier in sorted(by_tier):
            c = by_tier[tier]
            ratio = f"{c['nt']/c['aa']*100:.0f}%" if c["aa"] else "—"
            print(f"  {tier:<12} {c['total']:>7} {c['aa']:>7} {c['nt']:>7} {ratio:>9}")
    assert True
