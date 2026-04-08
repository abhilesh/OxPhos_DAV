"""
tests/test_parsers_v2.py

Unit tests for src/utils/parsers_v2.py

Covers:
  - MitomapParserV2: version extraction, plasmy classification, PubMed counting,
    MT-ND6 strand flip, Non-OXPHOS / stop-codon / non-standard AA-change filtering
  - ClinvarParserV2: HGVS-c extraction, HGVS-p extraction, protein missense /
    synonymous / stop-codon parsing, star ratings
  - MyVariantParserV2: metric extraction from nested dbNSFP / gnomAD dicts,
    missing-data returns None, popmax calculation
  - VariantRecord helpers: populate_substitution_properties, populate_gene_context

Run from project root:
    pytest tests/test_parsers_v2.py -v
"""

import json
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from utils.parsers_v2 import (
    MitomapParserV2,
    ClinvarParserV2,
    MyVariantParserV2,
    _first,
    _float,
    _int,
)
from utils.variant_record import VariantRecord


# ── Helpers ──────────────────────────────────────────────────────────────────


class TestHelpers:
    def test_first_scalar(self):
        assert _first(3.14) == 3.14

    def test_first_list(self):
        assert _first([7, 8, 9]) == 7

    def test_first_empty_list(self):
        assert _first([]) is None

    def test_first_none_returns_default(self):
        assert _first(None, default=0) == 0

    def test_float_from_scalar(self):
        assert _float("0.75") == pytest.approx(0.75)

    def test_float_from_list(self):
        assert _float([0.5, 0.9]) == pytest.approx(0.5)

    def test_float_none_returns_default(self):
        assert _float(None) is None

    def test_float_invalid_string(self):
        assert _float("n/a") is None

    def test_int_from_scalar(self):
        assert _int("42") == 42

    def test_int_none_returns_default(self):
        assert _int(None) is None


# ── MitomapParserV2 ───────────────────────────────────────────────────────────


class TestMitomapParserV2:

    def test_extract_version_from_filename(self):
        p = Path("MITOMAP_CodingVariants_2026-01-15.tsv")
        assert MitomapParserV2._extract_version(p) == "2026-01-15"

    def test_extract_version_missing(self):
        p = Path("MITOMAP_CodingVariants.tsv")
        assert MitomapParserV2._extract_version(p) == ""

    @pytest.mark.parametrize("homo,hetero,expected", [
        ("yes", "",    "homo"),
        ("",    "yes", "hetero"),
        ("yes", "yes", "homo/hetero"),
        ("no",  "no",  ""),
        ("nr",  "nan", ""),
        ("",    "",    ""),
    ])
    def test_plasmy(self, homo, hetero, expected):
        assert MitomapParserV2._plasmy(homo, hetero) == expected

    @pytest.mark.parametrize("pubmed_field,expected", [
        ("12345",           1),
        ("12345,67890",     2),
        ("12345, 67890",    2),
        ("",                0),
        ("nan",             0),
        ("12345,,67890",    2),  # empty token after split ignored
    ])
    def test_pubmed_count(self, pubmed_field, expected):
        assert MitomapParserV2._pubmed_count(pubmed_field) == expected

    @pytest.fixture
    def minimal_hgnc(self, tmp_path):
        """GeneReference with just enough MT genes for parse() tests."""
        from utils.parsers import GeneReference
        content = (
            "HGNC ID\tApproved symbol\tApproved name\tStatus\tLocus type\t"
            "Previous symbols\tAlias symbols\tChromosome\tNCBI Gene ID\t"
            "Ensembl gene ID\tVega gene ID\tGroup ID\tGroup name\n"
            "HGNC:1\tMT-ND1\tND1\tApproved\tgene with protein product\t\t\tMT\t4535\tENSG0001\t\t1\tCI\n"
            "HGNC:2\tMT-ND6\tND6\tApproved\tgene with protein product\t\t\tMT\t4541\tENSG0002\t\t1\tCI\n"
        )
        f = tmp_path / "hgnc.tsv"
        f.write_text(content)

        coord_content = (
            "gene\tchr\tstart\tend\tstrand\n"
            "MT-ND1\tMT\t3307\t4262\t+\n"
            "MT-ND6\tMT\t14673\t14149\t-\n"
        )
        c = tmp_path / "coords.tsv"
        c.write_text(coord_content)

        ref = GeneReference(f)
        ref.load_coordinates(c)
        return ref

    @pytest.fixture
    def mitomap_tsv(self, tmp_path):
        """Minimal MITOMAP TSV with a few representative rows."""
        rows = [
            "id\tpos\tref\talt\taachange\tdisease\tstatus\thomoplasmy\theteroplasmy\tpubmed_ids",
            # Valid ND1 variant
            "1\t3697\tG\tA\tA52T\tLHON\tCfrm\tyes\t\t12345",
            # ND6 (minus-strand) variant — ref/alt must be complemented
            "2\t14484\tT\tC\tM64V\tLHON\tCfrm\t\tyes\t11111,22222",
            # Non-OXPHOS position — must be skipped
            "3\t5000\tA\tG\tA10T\tSNP\tRep\t\t\t",
            # Non-coding aachange — must be skipped
            "4\t3697\tG\tA\tnoncoding\ttest\tRep\t\t\t",
            # Stop codon in aachange — must be skipped
            "5\t3697\tG\tT\tA52*\ttest\tRep\t\t\t",
        ]
        f = tmp_path / "MITOMAP_CodingVariants_2026-01-01.tsv"
        f.write_text("\n".join(rows), encoding="windows-1252")
        return f

    def test_parse_valid_nd1_variant(self, minimal_hgnc, mitomap_tsv):
        parser = MitomapParserV2(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        nd1 = [r for r in results if r["locus"] == "MT-ND1"]
        assert len(nd1) >= 1
        r = nd1[0]
        assert r["locus"] == "MT-ND1"
        assert r["aa_change"] == "A52T"
        assert r["ref_aa"] == "A"
        assert r["alt_aa"] == "T"
        assert r["is_synonymous"] is False
        assert r["source_db"] == "MITOMAP"
        assert r["source_db_version"] == "2026-01-01"
        assert r["mitomap_plasmy"] == "homo"
        assert r["mitomap_pubmed_count"] == 1

    def test_parse_nd6_strand_flip(self, minimal_hgnc, mitomap_tsv):
        """MT-ND6 is on the minus strand; T>C must become A>G after complement."""
        parser = MitomapParserV2(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        nd6 = [r for r in results if r["locus"] == "MT-ND6"]
        assert len(nd6) == 1
        # genomic alleles unchanged
        assert nd6[0]["genomic_ref"] == "T"
        assert nd6[0]["genomic_alt"] == "C"
        # CDS alleles complemented
        assert nd6[0]["ref"] == "A"
        assert nd6[0]["alt"] == "G"

    def test_non_oxphos_position_skipped(self, minimal_hgnc, mitomap_tsv):
        parser = MitomapParserV2(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        positions = [r["rCRS_pos"] for r in results]
        assert 5000 not in positions

    def test_noncoding_aachange_skipped(self, minimal_hgnc, mitomap_tsv):
        parser = MitomapParserV2(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        aachanges = [r["aa_change"] for r in results]
        assert "noncoding" not in aachanges

    def test_stop_codon_aachange_skipped(self, minimal_hgnc, mitomap_tsv):
        parser = MitomapParserV2(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        aachanges = [r["aa_change"] for r in results]
        assert not any("*" in a for a in aachanges)

    def test_hgvs_c_format(self, minimal_hgnc, mitomap_tsv):
        parser = MitomapParserV2(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        for r in results:
            assert r["hgvs_c"].startswith("m."), f"bad hgvs_c: {r['hgvs_c']}"


# ── ClinvarParserV2 ───────────────────────────────────────────────────────────


class TestClinvarParserV2:

    @pytest.fixture
    def parser(self):
        return ClinvarParserV2()

    # _extract_hgvs_c
    @pytest.mark.parametrize("name,hgvs_c,transcript,cds_ref,cds_alt", [
        (
            "NM_001369.4(NDUFS1):c.1057C>T (p.Arg353Cys)",
            "NM_001369.4:c.1057C>T", "NM_001369.4", "C", "T",
        ),
        (
            "NM_005169.3(ATP5F1A):c.327G>C (p.Leu109Phe)",
            "NM_005169.3:c.327G>C", "NM_005169.3", "G", "C",
        ),
        (
            "no_transcript:c.100A>G",
            "c.100A>G", "", "A", "G",
        ),
        (
            "NM_000059.4(BRCA2):p.Lys207Lys",
            "", "NM_000059.4", "", "",
        ),
    ])
    def test_extract_hgvs_c(self, parser, name, hgvs_c, transcript, cds_ref, cds_alt):
        got_hgvs_c, got_tx, got_ref, got_alt = parser._extract_hgvs_c(name)
        assert got_hgvs_c == hgvs_c
        assert got_tx == transcript
        assert got_ref == cds_ref
        assert got_alt == cds_alt

    # _extract_hgvs_p
    @pytest.mark.parametrize("name,expected", [
        ("NM_001369.4(NDUFS1):c.1057C>T (p.Arg353Cys)", "p.Arg353Cys"),
        ("NM_005169.3(ATP5F1A):c.327G>C (p.Leu109Phe)", "p.Leu109Phe"),
        ("NM_000059.4(BRCA2):c.621G>A",                 ""),
    ])
    def test_extract_hgvs_p(self, parser, name, expected):
        assert parser._extract_hgvs_p(name) == expected

    # _parse_protein
    @pytest.mark.parametrize("name,expected_aa,expected_syn", [
        ("NM_001369.4(NDUFS1):c.1057C>T (p.Arg353Cys)", "R353C", False),
        ("NM_005169.3(ATP5F1A):c.327G>C (p.Leu109Phe)", "L109F", False),
        ("NM_000059.4(BRCA2):c.621G>A (p.Lys207Lys)",   "K207K", True),
    ])
    def test_protein_missense_and_synonymous(self, parser, name, expected_aa, expected_syn):
        aa, is_syn = parser._parse_protein(name)
        assert aa == expected_aa
        assert is_syn is expected_syn

    def test_protein_stop_codon_rejected(self, parser):
        aa, _ = parser._parse_protein("NM_000059.4:c.100C>T (p.Arg34*)")
        assert aa is None

    def test_protein_no_p_dot_rejected(self, parser):
        aa, _ = parser._parse_protein("NM_000059.4:c.100C>T")
        assert aa is None

    def test_protein_unknown_aa_code(self, parser):
        """Three-letter codes not in the table map to X, causing rejection."""
        aa, _ = parser._parse_protein("NM_000059.4:c.100C>T (p.Xaa34Ala)")
        assert aa is None or "X" in (aa or "")

    # _stars
    @pytest.mark.parametrize("review,expected", [
        ("practice guideline",                                   4),
        ("reviewed by expert panel",                             3),
        ("criteria provided, multiple submitters, no conflicts", 2),
        ("criteria provided, single submitter",                  1),
        ("criteria provided, conflicting classifications",       1),
        ("-",                                                    0),
        ("no assertion criteria provided",                       0),
        ("no classification provided",                           0),
    ])
    def test_stars(self, parser, review, expected):
        assert parser._stars(review) == expected


# ── MyVariantParserV2 ─────────────────────────────────────────────────────────


@pytest.fixture
def myvariant_json(tmp_path):
    """Minimal pre-downloaded MyVariant JSON with one complete record."""
    records = [
        {
            "_id": "chr1:g.100A>G",
            "dbnsfp": {
                "revel": {"score": 0.82},
                "alphamissense": {"score": 0.91, "pred": "lp"},
                "esm1b": {"score": -4.5},
                "mpc": {"score": 1.7},
                "phylop": {"100way_vertebrate": {"score": 5.3}},
                "gerp++": {"rs": 4.1},
            },
            "gnomad_genome": {
                "af": {
                    "af": 0.0012,
                    "af_afr": 0.0020,
                    "af_amr": 0.0005,
                    "af_asj": 0.0010,
                    "af_eas": 0.0001,
                    "af_fin": 0.0030,
                    "af_nfe": 0.0015,
                    "af_oth": 0.0008,
                },
                "ac": {"ac": 30},
                "an": {"an": 25000},
                "hom": {"hom": 2},
            },
        },
        # Record missing from MyVariant
        {"_id": "chr2:g.200C>T", "notfound": True},
    ]
    f = tmp_path / "MyVariant_dbNSFP_gnomAD_2026-01-01.json"
    f.write_text(json.dumps(records))
    return f


class TestMyVariantParserV2:

    @pytest.fixture
    def mv(self, myvariant_json):
        return MyVariantParserV2(myvariant_json)

    def test_revel_extracted(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["revel_score"] == pytest.approx(0.82)

    def test_alphamissense_score_extracted(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["alphamissense_score"] == pytest.approx(0.91)

    def test_alphamissense_class_lp(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["alphamissense_class"] == "likely_pathogenic"

    @pytest.mark.parametrize("pred,expected_class", [
        ("lp", "likely_pathogenic"),
        ("p",  "likely_pathogenic"),
        ("lb", "likely_benign"),
        ("b",  "likely_benign"),
        ("am", "ambiguous"),
    ])
    def test_alphamissense_class_mapping(self, myvariant_json, tmp_path, pred, expected_class):
        records = [{"_id": "chr9:g.1A>G", "dbnsfp": {"alphamissense": {"score": 0.5, "pred": pred}}}]
        f = tmp_path / f"MyVariant_{pred}.json"
        f.write_text(json.dumps(records))
        mv = MyVariantParserV2(f)
        m = mv.get_all_metrics("9", 1, "A", "G")
        assert m["alphamissense_class"] == expected_class

    def test_esm1b_extracted(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["esm1b_score"] == pytest.approx(-4.5)

    def test_mpc_extracted(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["mpc_score"] == pytest.approx(1.7)

    def test_phylop_extracted(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["phylop_100vert"] == pytest.approx(5.3)

    def test_gerp_extracted(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["gerp_rs"] == pytest.approx(4.1)

    def test_gnomad_af_global(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["gnomad_af_global"] == pytest.approx(0.0012)

    def test_gnomad_ac_an(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["gnomad_ac"] == 30
        assert m["gnomad_an"] == 25000

    def test_gnomad_nhomalt(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["gnomad_nhomalt"] == 2

    def test_gnomad_popmax_is_fin(self, mv):
        """af_fin (0.003) is the highest population AF and must be the popmax."""
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["gnomad_af_popmax"] == pytest.approx(0.003)
        assert m["gnomad_popmax_pop"] == "fin"

    def test_missing_variant_returns_none(self, mv):
        """Variants not in the JSON (notfound=True) must return None for all metrics.
        gnomAD keys are absent (not None) when there is no gnomAD block — callers
        use .get() which returns None, matching the VariantRecord defaults."""
        m = mv.get_all_metrics("2", 200, "C", "T")
        assert m["revel_score"] is None
        assert m.get("gnomad_af_global") is None
        assert m["alphamissense_score"] is None

    def test_variant_not_in_lookup_returns_none(self, mv):
        """Unknown variant key must return None for all dbNSFP metrics."""
        m = mv.get_all_metrics("99", 9999, "A", "T")
        assert m["revel_score"] is None
        assert m.get("gnomad_af_global") is None

    def test_no_removed_fields_returned(self, mv):
        """Removed fields (cadd_phred, sift_score, polyphen2, phastcons) must not appear."""
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert "cadd_phred" not in m
        assert "sift_score" not in m
        assert "polyphen2_hdiv_score" not in m
        assert "phastcons_100vert" not in m


# ── VariantRecord helpers ─────────────────────────────────────────────────────


def _make_record(**kwargs) -> VariantRecord:
    defaults = dict(
        ann_id="test",
        source_db="ClinVar",
        source_db_version="2026-01-01",
        source_record_id="12345",
        genome="nucDNA",
        reference_assembly="GRCh38",
        locus="NDUFS1",
        nc_change="c.100C>T",
        aa_change="R34C",
        ref_nt="C",
        alt_nt="T",
        ref_aa="R",
        alt_aa="C",
        genomic_pos=100,
        is_synonymous=False,
    )
    defaults.update(kwargs)
    return VariantRecord(**defaults)


class TestVariantRecordSubstitutionProperties:

    def test_blosum62_populated(self):
        r = _make_record(ref_aa="R", alt_aa="C")
        r.populate_substitution_properties()
        assert r.blosum62 is not None
        assert isinstance(r.blosum62, int)

    def test_miyata_populated(self):
        r = _make_record(ref_aa="R", alt_aa="C")
        r.populate_substitution_properties()
        assert r.miyata_distance is not None
        assert r.miyata_distance >= 0

    def test_synonymous_substitution_blosum_is_positive(self):
        """Same amino acid (synonymous) should have a positive BLOSUM62 self-score."""
        r = _make_record(ref_aa="A", alt_aa="A")
        r.populate_substitution_properties()
        assert r.blosum62 is not None
        assert r.blosum62 > 0

    def test_charge_change_neutral_to_positive(self):
        """Ala (neutral) → Arg (positive) must give charge_change = +1."""
        r = _make_record(ref_aa="A", alt_aa="R")
        r.populate_substitution_properties()
        assert r.charge_change == 1

    def test_charge_change_negative_to_positive(self):
        """Asp (−1) → Arg (+1) = charge_change of +2."""
        r = _make_record(ref_aa="D", alt_aa="R")
        r.populate_substitution_properties()
        assert r.charge_change == 2

    def test_proline_flag(self):
        r = _make_record(ref_aa="A", alt_aa="P")
        r.populate_substitution_properties()
        assert r.is_proline_involved is True

    def test_glycine_flag(self):
        r = _make_record(ref_aa="G", alt_aa="A")
        r.populate_substitution_properties()
        assert r.is_glycine_involved is True

    def test_cysteine_flag(self):
        r = _make_record(ref_aa="C", alt_aa="S")
        r.populate_substitution_properties()
        assert r.is_cysteine_involved is True

    def test_empty_aa_skips_computation(self):
        r = _make_record(ref_aa="", alt_aa="")
        r.populate_substitution_properties()
        assert r.blosum62 is None
        assert r.miyata_distance is None

    def test_removed_fields_absent(self):
        """Fields removed from VariantRecord must not exist as attributes."""
        r = _make_record()
        for removed in ("cadd_phred", "sift_score", "polyphen2_hdiv_score",
                        "phastcons_100vert", "gene_loeuf", "gene_pli",
                        "gene_mis_z", "mitomap_gb_freq", "haplogroup_label"):
            assert not hasattr(r, removed), f"Removed field '{removed}' still present"


class TestVariantRecordGeneContext:

    @pytest.mark.parametrize("locus,expected_complex", [
        ("MT-ND1",   "CI"),
        ("NDUFS1",   "CI"),
        ("SDHA",     "CII"),
        ("MT-CYB",   "CIII"),
        ("MT-CO1",   "CIV"),
        ("MT-ATP6",  "CV"),
        ("NDUFA4",   "CIV"),   # NDUFA4 is annotated as CIV
    ])
    def test_complex_assignment(self, locus, expected_complex):
        r = _make_record(locus=locus)
        r.populate_gene_context()
        assert r.complex_id == expected_complex, f"{locus}: expected {expected_complex}, got {r.complex_id}"

    @pytest.mark.parametrize("locus,expected_encoded", [
        ("MT-ND1",  "mitochondrial"),
        ("MT-CYB",  "mitochondrial"),
        ("NDUFS1",  "nuclear"),
        ("SDHA",    "nuclear"),
    ])
    def test_encoded_by(self, locus, expected_encoded):
        r = _make_record(locus=locus)
        r.populate_gene_context()
        assert r.encoded_by == expected_encoded

    def test_is_sdh_true(self):
        r = _make_record(locus="SDHA")
        r.populate_gene_context()
        assert r.is_sdh is True

    def test_is_sdh_false(self):
        r = _make_record(locus="MT-ND1")
        r.populate_gene_context()
        assert r.is_sdh is False

    def test_subunit_role_core(self):
        r = _make_record(locus="NDUFS1")
        r.populate_gene_context()
        assert r.subunit_role == "core"

    def test_subunit_role_accessory(self):
        r = _make_record(locus="NDUFA1")
        r.populate_gene_context()
        assert r.subunit_role == "accessory"

    def test_locus_slash_split(self):
        """Loci like 'GENE/ALIAS' should resolve via the first token."""
        r = _make_record(locus="MT-ND1/MTND1")
        r.populate_gene_context()
        assert r.complex_id == "CI"
