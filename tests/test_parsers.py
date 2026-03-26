import pytest
import pandas as pd
from pathlib import Path
from utils.hgnc_parser import GeneReference
from utils.mitomap_parser import MitomapParser
from utils.clinvar_parser import ClinvarParser


@pytest.fixture
def mock_hgnc_file(tmp_path):
    """Creates a temporary HGNC TSV file for testing."""
    hgnc_data = (
        "HGNC ID\tApproved symbol\tGroup name\tAlias symbols\tPrevious symbols\n"
        "HGNC:10680\tSDHA\tMitochondrial complex II\tFP\tSDH2\n"
        "HGNC:12585\tUQCRC1\tMitochondrial complex III\t\tQCR1\n"
    )
    file_path = tmp_path / "mock_hgnc.tsv"
    file_path.write_text(hgnc_data, encoding="utf-8")
    return file_path


@pytest.fixture
def mock_mitomap_file(tmp_path):
    """Creates a temporary MITOMAP TSV file with controlled edge cases."""
    mitomap_data = (
        "id\tpos\tref\talt\taachange\tdisease\tstatus\n"
        "1\t72\tT\tC\tnoncoding\tNone\tReported\n"  # Should be dropped: Non-OXPHOS & noncoding
        "2\t3308\tG\tA\tA1T\tLHON\tCfrm [P]\n"  # Should be kept: Valid missense in MT-ND1
        "3\t3309\tC\tT\tA1A\tNone\tReported [B]\n"  # Should be kept: Valid synonymous in MT-ND1
        "4\t4470\tA\tG\tframeshift\tNone\tReported\n"  # Should be dropped: Frameshift
        "5\t99999\tA\tG\tM1V\tNone\tReported\n"  # Should be dropped: Invalid coordinate
    )
    file_path = tmp_path / "mock_mitomap.tsv"
    file_path.write_text(mitomap_data, encoding="windows-1252")
    return file_path


def test_gene_reference_logic(mock_hgnc_file):
    """Verifies HGNC parsing and alias resolution."""
    ref = GeneReference(mock_hgnc_file)

    # Test canonical names
    assert ref.is_target("SDHA") is True
    assert ref.is_target("UQCRC1") is True

    # Test aliases and previous symbols
    assert ref.is_target("FP") is True
    assert ref.is_target("QCR1") is True

    # Test negative case
    assert ref.is_target("BRCA1") is False

    # Test metadata retrieval
    data = ref.get_gene_data("SDHA")
    assert data["group"] == "Mitochondrial complex II"


def test_mitomap_parser_logic(mock_hgnc_file, mock_mitomap_file):
    """Verifies structural filtering of MITOMAP variants."""
    hgnc_ref = GeneReference(mock_hgnc_file)
    parser = MitomapParser(hgnc_reference=hgnc_ref)

    results = parser.parse(mock_mitomap_file)

    # Out of 5 mock rows, only 2 should survive the strict structural filters
    assert len(results) == 2

    missense_var = next(v for v in results if v["aa_change"] == "A1T")
    assert missense_var["is_synonymous"] is False
    assert missense_var["locus"] == "MT-ND1"

    synonymous_var = next(v for v in results if v["aa_change"] == "A1A")
    assert synonymous_var["is_synonymous"] is True


def test_clinvar_protein_parser():
    """Unit tests the internal regex of the ClinvarParser without needing a file."""
    parser = ClinvarParser()

    # Valid Missense
    aa, is_syn = parser._parse_protein("NM_005169.3(ATP5F1A):c.327G>C (p.Leu109Phe)")
    assert aa == "L109F"
    assert is_syn is False

    # Valid Synonymous
    aa, is_syn = parser._parse_protein("p.Arg207Arg")
    assert aa == "R207R"
    assert is_syn is True

    # Stop codon (Should be rejected)
    aa, is_syn = parser._parse_protein("p.Arg207*")
    assert aa is None

    # Invalid string
    aa, is_syn = parser._parse_protein("c.327G>C")
    assert aa is None
