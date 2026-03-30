"""
Data-driven invariant tests for c-DAR classification output.

Validates the biological and logical constraints that must hold over the full
cdar_classifications_*.json datasets produced by 00_classify_DAV.py:

  INV-1  NT ⊆ AA (subset):       cdar_nt=True  →  cdar_aa=True
  INV-2  Species subset:          cdar_nt_species ⊆ cdar_aa_species
  INV-3  Count monotonicity:      len(cdar_aa_species) >= len(cdar_nt_species)
  INV-4  Count field consistency: compensating_species_count == len(cdar_aa_species)
  INV-5  Global count:            total aa_cDARs >= total nt_cDARs  (per genome, per tier)
  INV-6  No empty species lists on positive flags:
             cdar_aa=True  →  len(cdar_aa_species) > 0
             cdar_nt=True  →  len(cdar_nt_species) > 0
  INV-7  Discarded/synonymous variants must not appear in c-DAR output

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
MT_JSON  = CURATED / "cdar_classifications_mtDNA.json"
NUC_JSON = CURATED / "cdar_classifications_nucDNA.json"


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


@pytest.fixture(scope="module", params=["mtDNA", "nucDNA"])
def all_variants(request, mt_variants, nuc_variants):
    return mt_variants if request.param == "mtDNA" else nuc_variants


# ── Helper ────────────────────────────────────────────────────────────────────

def _label(v: dict) -> str:
    return f"{v.get('locus','?')} | {v.get('nc_change','?')} | {v.get('aa_change','?')}"


# ── INV-1: NT ⊆ AA (subset rule) ─────────────────────────────────────────────

class TestNtSubsetAa:

    def test_nt_implies_aa_mtdna(self, mt_variants):
        """Every mtDNA variant with cdar_nt=True must also have cdar_aa=True."""
        violations = [
            _label(v) for v in mt_variants
            if v.get("cdar_nt") and not v.get("cdar_aa")
        ]
        assert not violations, (
            f"INV-1 violated for {len(violations)} mtDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_nt_implies_aa_nucdna(self, nuc_variants):
        """Every nucDNA variant with cdar_nt=True must also have cdar_aa=True."""
        violations = [
            _label(v) for v in nuc_variants
            if v.get("cdar_nt") and not v.get("cdar_aa")
        ]
        assert not violations, (
            f"INV-1 violated for {len(violations)} nucDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )


# ── INV-2: Species subset ─────────────────────────────────────────────────────

class TestSpeciesSubset:

    def _check(self, variants: list, genome: str):
        violations = []
        for v in variants:
            aa_sp = set(v.get("cdar_aa_species") or [])
            nt_sp = set(v.get("cdar_nt_species") or [])
            extra = nt_sp - aa_sp
            if extra:
                violations.append(
                    f"{_label(v)} — nt_species not in aa_species: {sorted(extra)[:3]}"
                )
        return violations

    def test_nt_species_subset_of_aa_species_mtdna(self, mt_variants):
        violations = self._check(mt_variants, "mtDNA")
        assert not violations, (
            f"INV-2 violated for {len(violations)} mtDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_nt_species_subset_of_aa_species_nucdna(self, nuc_variants):
        violations = self._check(nuc_variants, "nucDNA")
        assert not violations, (
            f"INV-2 violated for {len(violations)} nucDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )


# ── INV-3: Count monotonicity ─────────────────────────────────────────────────

class TestCountMonotonicity:

    def _check(self, variants: list):
        violations = []
        for v in variants:
            n_aa = len(v.get("cdar_aa_species") or [])
            n_nt = len(v.get("cdar_nt_species") or [])
            if n_nt > n_aa:
                violations.append(
                    f"{_label(v)} — aa_species={n_aa} < nt_species={n_nt}"
                )
        return violations

    def test_aa_count_gte_nt_count_mtdna(self, mt_variants):
        violations = self._check(mt_variants)
        assert not violations, (
            f"INV-3 violated for {len(violations)} mtDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_aa_count_gte_nt_count_nucdna(self, nuc_variants):
        violations = self._check(nuc_variants)
        assert not violations, (
            f"INV-3 violated for {len(violations)} nucDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )


# ── INV-4: compensating_species_count field consistency ───────────────────────

class TestSpeciesCountField:

    def _check(self, variants: list):
        violations = []
        for v in variants:
            stored  = v.get("compensating_species_count", -1)
            derived = len(v.get("cdar_aa_species") or [])
            if stored != derived:
                violations.append(
                    f"{_label(v)} — stored={stored}, len(aa_species)={derived}"
                )
        return violations

    def test_count_field_consistent_mtdna(self, mt_variants):
        violations = self._check(mt_variants)
        assert not violations, (
            f"INV-4 violated for {len(violations)} mtDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_count_field_consistent_nucdna(self, nuc_variants):
        violations = self._check(nuc_variants)
        assert not violations, (
            f"INV-4 violated for {len(violations)} nucDNA variant(s):\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )


# ── INV-5: Global and per-tier count monotonicity ─────────────────────────────

class TestGlobalCountMonotonicity:

    def _tier_counts(self, variants: list) -> dict:
        """Returns {tier: {"aa": int, "nt": int}} plus a "TOTAL" key."""
        counts = defaultdict(lambda: {"aa": 0, "nt": 0})
        for v in variants:
            tier = v.get("tier", "Unknown")
            if v.get("cdar_aa"):
                counts[tier]["aa"] += 1
                counts["TOTAL"]["aa"] += 1
            if v.get("cdar_nt"):
                counts[tier]["nt"] += 1
                counts["TOTAL"]["nt"] += 1
        return dict(counts)

    def test_global_aa_gte_nt_mtdna(self, mt_variants):
        c = self._tier_counts(mt_variants)
        total = c.get("TOTAL", {"aa": 0, "nt": 0})
        assert total["aa"] >= total["nt"], (
            f"Global mtDNA: aa_cDARs={total['aa']} < nt_cDARs={total['nt']}"
        )

    def test_global_aa_gte_nt_nucdna(self, nuc_variants):
        c = self._tier_counts(nuc_variants)
        total = c.get("TOTAL", {"aa": 0, "nt": 0})
        assert total["aa"] >= total["nt"], (
            f"Global nucDNA: aa_cDARs={total['aa']} < nt_cDARs={total['nt']}"
        )

    def test_per_tier_aa_gte_nt_mtdna(self, mt_variants):
        counts = self._tier_counts(mt_variants)
        violations = [
            f"{tier}: aa={c['aa']} < nt={c['nt']}"
            for tier, c in counts.items()
            if tier != "TOTAL" and c["nt"] > c["aa"]
        ]
        assert not violations, (
            "INV-5 per-tier violations (mtDNA):\n"
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
            "INV-5 per-tier violations (nucDNA):\n"
            + "\n".join(f"  {x}" for x in violations)
        )


# ── INV-6: Non-empty species lists on positive flags ─────────────────────────

class TestNonEmptySpeciesOnPositiveFlag:

    def _check_aa(self, variants: list):
        return [
            _label(v) for v in variants
            if v.get("cdar_aa") and not v.get("cdar_aa_species")
        ]

    def _check_nt(self, variants: list):
        return [
            _label(v) for v in variants
            if v.get("cdar_nt") and not v.get("cdar_nt_species")
        ]

    def test_aa_flag_has_species_mtdna(self, mt_variants):
        violations = self._check_aa(mt_variants)
        assert not violations, (
            f"INV-6 (aa): {len(violations)} mtDNA variants flagged cdar_aa=True "
            f"but have empty aa_species list:\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_aa_flag_has_species_nucdna(self, nuc_variants):
        violations = self._check_aa(nuc_variants)
        assert not violations, (
            f"INV-6 (aa): {len(violations)} nucDNA variants flagged cdar_aa=True "
            f"but have empty aa_species list:\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_nt_flag_has_species_mtdna(self, mt_variants):
        violations = self._check_nt(mt_variants)
        assert not violations, (
            f"INV-6 (nt): {len(violations)} mtDNA variants flagged cdar_nt=True "
            f"but have empty nt_species list:\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )

    def test_nt_flag_has_species_nucdna(self, nuc_variants):
        violations = self._check_nt(nuc_variants)
        assert not violations, (
            f"INV-6 (nt): {len(violations)} nucDNA variants flagged cdar_nt=True "
            f"but have empty nt_species list:\n"
            + "\n".join(f"  {x}" for x in violations[:20])
        )


# ── INV-7: Discarded / synonymous variants absent from output ─────────────────

class TestDiscardedAndSynonymousAbsent:

    def test_no_discarded_in_mtdna_output(self, mt_variants):
        discarded = [_label(v) for v in mt_variants if v.get("tier") == "Discarded"]
        assert not discarded, (
            f"INV-7: {len(discarded)} Discarded mtDNA variants appear in c-DAR output "
            f"(should have been filtered before classification):\n"
            + "\n".join(f"  {x}" for x in discarded[:10])
        )

    def test_no_discarded_in_nucdna_output(self, nuc_variants):
        discarded = [_label(v) for v in nuc_variants if v.get("tier") == "Discarded"]
        assert not discarded, (
            f"INV-7: {len(discarded)} Discarded nucDNA variants appear in c-DAR output:\n"
            + "\n".join(f"  {x}" for x in discarded[:10])
        )

    def test_no_synonymous_in_mtdna_output(self, mt_variants):
        synonymous = [_label(v) for v in mt_variants if v.get("is_synonymous")]
        assert not synonymous, (
            f"INV-7: {len(synonymous)} synonymous mtDNA variants appear in c-DAR output:\n"
            + "\n".join(f"  {x}" for x in synonymous[:10])
        )

    def test_no_synonymous_in_nucdna_output(self, nuc_variants):
        synonymous = [_label(v) for v in nuc_variants if v.get("is_synonymous")]
        assert not synonymous, (
            f"INV-7: {len(synonymous)} synonymous nucDNA variants appear in c-DAR output:\n"
            + "\n".join(f"  {x}" for x in synonymous[:10])
        )


# ── Summary report (printed even on pass) ────────────────────────────────────

def test_summary_report(mt_variants, nuc_variants):
    """Non-failing summary printed alongside results for quick orientation."""
    for label, variants in [("mtDNA", mt_variants), ("nucDNA", nuc_variants)]:
        total   = len(variants)
        n_aa    = sum(1 for v in variants if v.get("cdar_aa"))
        n_nt    = sum(1 for v in variants if v.get("cdar_nt"))
        by_tier = defaultdict(lambda: {"total": 0, "aa": 0, "nt": 0})
        for v in variants:
            t = v.get("tier", "?")
            by_tier[t]["total"] += 1
            if v.get("cdar_aa"): by_tier[t]["aa"] += 1
            if v.get("cdar_nt"): by_tier[t]["nt"] += 1

        print(f"\n{label} ({total} variants in c-DAR output):")
        print(f"  AA c-DARs : {n_aa}  ({100*n_aa/total:.1f}%)")
        print(f"  NT c-DARs : {n_nt}  ({100*n_nt/total:.1f}%)")
        print(f"  NT/AA ratio: {n_nt/n_aa*100:.1f}%" if n_aa else "  NT/AA ratio: N/A")
        print(f"  {'Tier':<12} {'Total':>7} {'aa':>7} {'nt':>7} {'nt/aa %':>9}")
        for tier in sorted(by_tier):
            c = by_tier[tier]
            ratio = f"{c['nt']/c['aa']*100:.0f}%" if c["aa"] else "—"
            print(f"  {tier:<12} {c['total']:>7} {c['aa']:>7} {c['nt']:>7} {ratio:>9}")
    # Always passes — purely informational
    assert True
