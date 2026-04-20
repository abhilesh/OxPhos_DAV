"""
Invariant tests for the current classified master table.

These tests validate the active Parquet-based classification contract rather
than the old compatibility JSON outputs.

Run from project root:
    pytest tests/test_cdav_invariants.py -v
"""

from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
CLASSIFIED_DIR = ROOT / "data" / "derived" / "classified"
CLASSIFIED_MASTER = CLASSIFIED_DIR / "variants_master_classified.parquet"
CLASSIFIED_ALL = CLASSIFIED_DIR / "classified_all.parquet"
CLASSIFIED_CLEAN = CLASSIFIED_DIR / "classified_clean.parquet"
CLASSIFIED_WARNING = CLASSIFIED_DIR / "classified_warning.parquet"
QC_SUMMARY = CLASSIFIED_DIR / "classification_qc_summary.json"


def _load_parquet(path: Path) -> pd.DataFrame:
    if not path.exists():
        pytest.skip(f"Required classified output missing: {path.name}")
    return pd.read_parquet(path)


@pytest.fixture(scope="module")
def classified_master() -> pd.DataFrame:
    return _load_parquet(CLASSIFIED_MASTER)


@pytest.fixture(scope="module")
def classified_all() -> pd.DataFrame:
    return _load_parquet(CLASSIFIED_ALL)


@pytest.fixture(scope="module")
def classified_clean() -> pd.DataFrame:
    return _load_parquet(CLASSIFIED_CLEAN)


@pytest.fixture(scope="module")
def classified_warning() -> pd.DataFrame:
    return _load_parquet(CLASSIFIED_WARNING)


def _label(row: pd.Series) -> str:
    return f"{row.get('variant_id', '?')} | {row.get('classification_gene', '?')} | {row.get('aa_change', '?')}"


def _safe_len(value) -> int:
    if value is None:
        return 0
    try:
        return len(value)
    except TypeError:
        return 0


class TestClassifiedInvariants:

    def test_nt_implies_aa_for_classified_rows(self, classified_all):
        violations = classified_all[
            (classified_all["is_cdav_nucleotide"] == True)
            & (classified_all["is_cdav_amino_acid"] != True)
        ]
        assert violations.empty, (
            "NT-level cDAV must imply AA-level cDAV for classified rows:\n"
            + "\n".join(_label(row) for _, row in violations.head(20).iterrows())
        )

    def test_species_count_field_matches_lineage_list(self, classified_all):
        allele_lineages = classified_all["lineages_with_disease_allele"].apply(
            _safe_len
        )
        codon_lineages = classified_all["lineages_with_disease_codon"].apply(
            _safe_len
        )
        bad_aa = classified_all[
            classified_all["n_species_with_disease_allele"].fillna(-1) != allele_lineages
        ]
        bad_nt = classified_all[
            classified_all["n_species_with_disease_codon"].fillna(-1) != codon_lineages
        ]
        assert bad_aa.empty and bad_nt.empty, (
            "Stored species counts do not match lineage lists.\n"
            f"AA mismatches: {len(bad_aa)}; NT mismatches: {len(bad_nt)}"
        )

    def test_global_aa_cdav_count_gte_nt_cdav_count(self, classified_all):
        for genome in ("mtDNA", "nucDNA"):
            subset = classified_all[classified_all["genome"] == genome]
            aa_count = int(subset["is_cdav_amino_acid"].fillna(False).sum())
            nt_count = int(subset["is_cdav_nucleotide"].fillna(False).sum())
            assert aa_count >= nt_count, (
                f"{genome}: aa_cDAVs={aa_count} < nt_cDAVs={nt_count}"
            )

    def test_positive_cdav_flags_require_nonempty_lineages(self, classified_all):
        aa_violations = classified_all[
            (classified_all["is_cdav_amino_acid"] == True)
            & (classified_all["lineages_with_disease_allele"].apply(_safe_len) == 0)
        ]
        nt_violations = classified_all[
            (classified_all["is_cdav_nucleotide"] == True)
            & (classified_all["lineages_with_disease_codon"].apply(_safe_len) == 0)
        ]
        assert aa_violations.empty and nt_violations.empty, (
            "Positive cDAV flags must have non-empty lineage lists.\n"
            f"AA violations: {len(aa_violations)}; NT violations: {len(nt_violations)}"
        )

    def test_no_unresolved_row_counted_as_udav(self, classified_master):
        violations = classified_master[
            (classified_master["classification_status"] == "unresolved")
            & (
                (classified_master["is_udav_amino_acid"] == True)
                | (classified_master["is_udav_nucleotide"] == True)
            )
        ]
        assert violations.empty, (
            "Unresolved rows must not be counted as uDAV.\n"
            + "\n".join(_label(row) for _, row in violations.head(20).iterrows())
        )

    def test_classified_rows_have_valid_basis(self, classified_all):
        valid_basis = {"nt_and_aa", "nt_only", "aa_only", "no_disease_allele_detected"}
        violations = classified_all[~classified_all["classification_basis"].isin(valid_basis)]
        assert violations.empty, (
            "Classified rows have invalid classification_basis values.\n"
            + "\n".join(_label(row) for _, row in violations.head(20).iterrows())
        )

    def test_no_mitochondrial_clinvar_row_is_eligible_for_nuclear_branch(self, classified_master):
        violations = classified_master[
            (classified_master["source_db"] == "ClinVar")
            & (classified_master["encoded_by"] == "mitochondrial")
            & (classified_master["eligible_core_comparative_pipeline"] == True)
        ]
        assert violations.empty, (
            "Mitochondrial ClinVar rows must remain ineligible for the nuclear branch.\n"
            + "\n".join(_label(row) for _, row in violations.head(20).iterrows())
        )


class TestClassifiedSubsetOutputs:

    def test_subset_outputs_partition_classified_rows(self, classified_all, classified_clean, classified_warning):
        assert len(classified_all) == len(classified_clean) + len(classified_warning)

    def test_clean_subset_has_no_warning_reason(self, classified_clean):
        violations = classified_clean[
            classified_clean["classification_warning_reason"].notna()
        ]
        assert violations.empty, "classified_clean must not contain warning reasons"

    def test_warning_subset_has_warning_reason(self, classified_warning):
        violations = classified_warning[
            classified_warning["classification_warning_reason"].isna()
        ]
        assert violations.empty, "classified_warning rows must carry a warning reason"

    def test_warning_subset_matches_expected_flags(self, classified_warning):
        violations = classified_warning[
            (classified_warning["exception_applied"] != True)
            & (classified_warning["mismatch_code"] != "REF_ALLELE_MISMATCH")
        ]
        assert violations.empty, (
            "classified_warning rows must be driven by exception metadata or REF_ALLELE_MISMATCH"
        )


def test_summary_report(classified_master, classified_all, classified_clean, classified_warning):
    mt = classified_master[classified_master["genome"] == "mtDNA"]
    nuc = classified_master[classified_master["genome"] == "nucDNA"]
    print("\nCurrent classified master summary:")
    print(f"  rows_total: {len(classified_master)}")
    print(f"  mtDNA rows: {len(mt)}")
    print(f"  nucDNA rows: {len(nuc)}")
    print(f"  classified_all: {len(classified_all)}")
    print(f"  classified_clean: {len(classified_clean)}")
    print(f"  classified_warning: {len(classified_warning)}")
    assert QC_SUMMARY.exists(), "classification_qc_summary.json should exist after classification"
