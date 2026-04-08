"""
tests/test_parsers.py

Unit tests for src/utils/parsers.py and src/utils/alignment_parser.py

Covers:
  - Helper functions: _first, _float, _int
  - GeneReference: canonical lookup, previous-symbol aliasing, coord inversion, MT position
  - MitomapParser: version extraction, plasmy classification, PubMed counting,
    MT-ND6 strand flip, Non-OXPHOS / stop-codon filtering
  - ClinvarParser: HGVS-c/p extraction, protein missense / synonymous / stop-codon
    parsing, star ratings
  - MyVariantParser: metric extraction from nested dbNSFP / gnomAD dicts,
    missing-data returns None, popmax calculation
  - VariantRecord helpers: populate_substitution_properties, populate_gene_context
  - AlignmentParser: mask character filtering (X, !, *, -), mutant codon extraction,
    c-DAR detection counts

Run from project root:
    pytest tests/test_parsers.py -v
"""

import json
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from utils.parsers import (
    GeneReference,
    MitomapParser,
    ClinvarParser,
    MyVariantParser,
    _first,
    _float,
    _int,
)
from utils.alignment_parser import AlignmentParser
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


# ── GeneReference ─────────────────────────────────────────────────────────────


@pytest.fixture
def hgnc_file(tmp_path):
    """Minimal HGNC TSV covering canonical names, previous symbols, and minus-strand coords."""
    content = (
        "HGNC ID\tApproved symbol\tApproved name\tStatus\tLocus type\t"
        "Previous symbols\tAlias symbols\tChromosome\tNCBI Gene ID\t"
        "Ensembl gene ID\tVega gene ID\tGroup ID\tGroup name\n"
        "HGNC:1234\tCOXFA4\tCOX factor A4\tApproved\tgene with protein product\t"
        "NDUFA4\tMLRQ\t7\t100\tENSG0001\t\t1\tMitochondrial complex I\n"
        "HGNC:5678\tSDHB\tsuccinate dehydrogenase B\tApproved\tgene with protein product\t"
        "\t\t1\t200\tENSG0002\t\t2\tMitochondrial complex II\n"
        "HGNC:9999\tMT-ND6\tNADH dehydrogenase 6\tApproved\tgene with protein product\t"
        "\t\tMT\t300\tENSG0003\t\t3\tMitochondrial complex I\n"
    )
    f = tmp_path / "hgnc.tsv"
    f.write_text(content, encoding="utf-8")
    return f


@pytest.fixture
def coord_file(tmp_path):
    """Coordinate file with one plus-strand gene, one minus-strand gene (start > end),
    and the mtDNA MT-ND6 gene."""
    content = (
        "gene\tchr\tstart\tend\tstrand\n"
        "SDHB\t1\t17054032\t17018722\t-\n"   # inverted (minus strand)
        "COXFA4\t7\t10931943\t10940153\t+\n"  # normal (plus strand)
        "MT-ND6\tMT\t14673\t14149\t-\n"        # inverted mtDNA minus strand
    )
    f = tmp_path / "coords.tsv"
    f.write_text(content, encoding="utf-8")
    return f


class TestGeneReference:

    def test_canonical_symbol_in_lookup(self, hgnc_file):
        ref = GeneReference(hgnc_file)
        assert "COXFA4" in ref.lookup
        assert "SDHB" in ref.lookup

    def test_previous_symbol_aliased(self, hgnc_file):
        """NDUFA4 is a previous symbol of COXFA4 and must be accepted."""
        ref = GeneReference(hgnc_file)
        assert "NDUFA4" in ref.lookup

    def test_previous_symbol_resolves_to_canonical(self, hgnc_file):
        """Looking up NDUFA4 and COXFA4 must return the same entry dict."""
        ref = GeneReference(hgnc_file)
        assert ref.lookup["NDUFA4"] is ref.lookup["COXFA4"]

    def test_symbol_key_is_canonical(self, hgnc_file):
        """The 'symbol' key inside the entry must always be the approved name."""
        ref = GeneReference(hgnc_file)
        assert ref.lookup["NDUFA4"]["symbol"] == "COXFA4"
        assert ref.lookup["COXFA4"]["symbol"] == "COXFA4"

    def test_unknown_gene_not_in_lookup(self, hgnc_file):
        ref = GeneReference(hgnc_file)
        assert "BRCA1" not in ref.lookup

    def test_coord_inversion_plus_strand(self, hgnc_file, coord_file):
        """Plus-strand gene: start must be <= end after load_coordinates."""
        ref = GeneReference(hgnc_file)
        ref.load_coordinates(coord_file)
        d = ref.lookup["COXFA4"]
        assert d["start"] <= d["end"]
        assert d["strand"] == "+"

    def test_coord_inversion_minus_strand(self, hgnc_file, coord_file):
        """Minus-strand gene stored as start > end must be normalised to start <= end."""
        ref = GeneReference(hgnc_file)
        ref.load_coordinates(coord_file)
        d = ref.lookup["MT-ND6"]
        assert d["start"] <= d["end"], (
            f"MT-ND6 start={d['start']} > end={d['end']} — inversion not corrected"
        )
        assert d["strand"] == "-"

    def test_mt_nd6_position_lookup(self, hgnc_file, coord_file):
        """Any position within MT-ND6 (14149–14673) must resolve to 'MT-ND6'."""
        ref = GeneReference(hgnc_file)
        ref.load_coordinates(coord_file)
        assert ref.get_mt_locus_by_position(14163) == "MT-ND6"
        assert ref.get_mt_locus_by_position(14149) == "MT-ND6"
        assert ref.get_mt_locus_by_position(14673) == "MT-ND6"

    def test_non_mt_position_returns_non_oxphos(self, hgnc_file, coord_file):
        ref = GeneReference(hgnc_file)
        ref.load_coordinates(coord_file)
        assert ref.get_mt_locus_by_position(1) == "Non-OXPHOS"


# ── MitomapParser ─────────────────────────────────────────────────────────────


class TestMitomapParser:

    def test_extract_version_from_filename(self):
        p = Path("MITOMAP_CodingVariants_2026-01-15.tsv")
        assert MitomapParser._extract_version(p) == "2026-01-15"

    def test_extract_version_missing(self):
        p = Path("MITOMAP_CodingVariants.tsv")
        assert MitomapParser._extract_version(p) == ""

    @pytest.mark.parametrize("homo,hetero,expected", [
        ("yes", "",    "homo"),
        ("",    "yes", "hetero"),
        ("yes", "yes", "homo/hetero"),
        ("no",  "no",  ""),
        ("nr",  "nan", ""),
        ("",    "",    ""),
    ])
    def test_plasmy(self, homo, hetero, expected):
        assert MitomapParser._plasmy(homo, hetero) == expected

    @pytest.mark.parametrize("pubmed_field,expected", [
        ("12345",           1),
        ("12345,67890",     2),
        ("12345, 67890",    2),
        ("",                0),
        ("nan",             0),
        ("12345,,67890",    2),
    ])
    def test_pubmed_count(self, pubmed_field, expected):
        assert MitomapParser._pubmed_count(pubmed_field) == expected

    @pytest.fixture
    def minimal_hgnc(self, tmp_path):
        """GeneReference with just enough MT genes for parse() tests."""
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
            "1\t3697\tG\tA\tA52T\tLHON\tCfrm\tyes\t\t12345",
            "2\t14484\tT\tC\tM64V\tLHON\tCfrm\t\tyes\t11111,22222",
            "3\t5000\tA\tG\tA10T\tSNP\tRep\t\t\t",
            "4\t3697\tG\tA\tnoncoding\ttest\tRep\t\t\t",
            "5\t3697\tG\tT\tA52*\ttest\tRep\t\t\t",
        ]
        f = tmp_path / "MITOMAP_CodingVariants_2026-01-01.tsv"
        f.write_text("\n".join(rows), encoding="windows-1252")
        return f

    def test_parse_valid_nd1_variant(self, minimal_hgnc, mitomap_tsv):
        parser = MitomapParser(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        nd1 = [r for r in results if r["locus"] == "MT-ND1"]
        assert len(nd1) >= 1
        r = nd1[0]
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
        parser = MitomapParser(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        nd6 = [r for r in results if r["locus"] == "MT-ND6"]
        assert len(nd6) == 1
        assert nd6[0]["genomic_ref"] == "T"
        assert nd6[0]["genomic_alt"] == "C"
        assert nd6[0]["ref"] == "A"
        assert nd6[0]["alt"] == "G"

    def test_non_oxphos_position_skipped(self, minimal_hgnc, mitomap_tsv):
        parser = MitomapParser(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        assert 5000 not in [r["rCRS_pos"] for r in results]

    def test_noncoding_aachange_skipped(self, minimal_hgnc, mitomap_tsv):
        parser = MitomapParser(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        assert "noncoding" not in [r["aa_change"] for r in results]

    def test_stop_codon_aachange_skipped(self, minimal_hgnc, mitomap_tsv):
        parser = MitomapParser(minimal_hgnc)
        results = parser.parse(mitomap_tsv)
        assert not any("*" in a for a in [r["aa_change"] for r in results])


# ── ClinvarParser ─────────────────────────────────────────────────────────────


class TestClinvarParser:

    @pytest.fixture
    def parser(self):
        return ClinvarParser()

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

    @pytest.mark.parametrize("name,expected", [
        ("NM_001369.4(NDUFS1):c.1057C>T (p.Arg353Cys)", "p.Arg353Cys"),
        ("NM_005169.3(ATP5F1A):c.327G>C (p.Leu109Phe)", "p.Leu109Phe"),
        ("NM_000059.4(BRCA2):c.621G>A",                 ""),
    ])
    def test_extract_hgvs_p(self, parser, name, expected):
        assert parser._extract_hgvs_p(name) == expected

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

    @pytest.mark.parametrize("review,expected", [
        ("practice guideline",                                   4),
        ("reviewed by expert panel",                             3),
        ("criteria provided, multiple submitters, no conflicts", 2),
        ("criteria provided, single submitter",                  1),
        ("criteria provided, conflicting classifications",       1),
        ("-",                                                    0),
        ("no assertion criteria provided",                       0),
        ("no classification provided",                           0),
        ("no classification for the single variant",             0),
        ("no classifications from unflagged records",            0),
    ])
    def test_stars(self, parser, review, expected):
        assert parser._stars(review) == expected

    def test_locus_canonicalised_via_previous_symbol(self, hgnc_file):
        """ClinvarParser must remap old gene symbols to current approved names."""
        ref = GeneReference(hgnc_file)
        parser = ClinvarParser(hgnc_reference=ref)
        gene = "NDUFA4"
        entry = ref.lookup.get(gene, {})
        canonical = entry.get("symbol", gene)
        assert canonical == "COXFA4"


# ── MyVariantParser ───────────────────────────────────────────────────────────


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
        {"_id": "chr2:g.200C>T", "notfound": True},
    ]
    f = tmp_path / "MyVariant_dbNSFP_gnomAD_2026-01-01.json"
    f.write_text(json.dumps(records))
    return f


class TestMyVariantParser:

    @pytest.fixture
    def mv(self, myvariant_json):
        return MyVariantParser(myvariant_json)

    def test_revel_extracted(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["revel_score"] == pytest.approx(0.82)

    def test_alphamissense_score_extracted(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["alphamissense_score"] == pytest.approx(0.91)

    def test_alphamissense_class_lp(self, mv):
        m = mv.get_all_metrics("1", 100, "A", "G")
        assert m["alphamissense_class"] == "likely_pathogenic"

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
        m = mv.get_all_metrics("2", 200, "C", "T")
        assert m["revel_score"] is None
        assert m.get("gnomad_af_global") is None

    def test_variant_not_in_lookup_returns_none(self, mv):
        m = mv.get_all_metrics("99", 9999, "A", "T")
        assert m["revel_score"] is None

    def test_no_removed_fields_returned(self, mv):
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
        r = _make_record(ref_aa="A", alt_aa="A")
        r.populate_substitution_properties()
        assert r.blosum62 is not None
        assert r.blosum62 > 0

    def test_charge_change_neutral_to_positive(self):
        r = _make_record(ref_aa="A", alt_aa="R")
        r.populate_substitution_properties()
        assert r.charge_change == 1

    def test_charge_change_negative_to_positive(self):
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
        ("NDUFA4",   "CIV"),
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
        r = _make_record(locus="MT-ND1/MTND1")
        r.populate_gene_context()
        assert r.complex_id == "CI"


# ── AlignmentParser ───────────────────────────────────────────────────────────


@pytest.fixture
def aa_fasta(tmp_path):
    """Minimal 3-residue AA alignment for testing mask character filtering.

    Human:     M F A
    Species_1: M F A  (same as ref — never a c-DAR)
    Species_2: M Y A  (Y ≠ F — c-DAR if mut_aa='Y')
    Species_3: M X A  (X — uncertain, must be filtered)
    Species_4: M ! A  (! — frameshift marker, must be filtered)
    Species_5: M * A  (* — stop codon, must be filtered)
    """
    content = (
        ">Homo_sapiens|9606|hg38|ENST001\nMFA\n"
        ">Species_1|100|asm1|ENST001\nMFA\n"
        ">Species_2|200|asm2|ENST001\nMYA\n"
        ">Species_3|300|asm3|ENST001\nMXA\n"
        ">Species_4|400|asm4|ENST001\nM!A\n"
        ">Species_5|500|asm5|ENST001\nM*A\n"
    )
    f = tmp_path / "test_aa.fasta"
    f.write_text(content)
    return f


@pytest.fixture
def nt_fasta(tmp_path):
    """9-nt codon alignment matching the aa_fasta above.

    Human:     ATG TTC GCC  (M F A)
    Species_1: ATG TTC GCC  (M F A — same as ref)
    Species_2: ATG TAC GCC  (M Y A — TAC = Tyr)
    Species_3: ATG XXX GCC  (X-masked)
    Species_4: ATG !!! GCC  (!-masked)
    Species_5: ATG *** GCC  (*-masked)
    """
    content = (
        ">Homo_sapiens|9606|hg38|ENST001\nATGTTCGCC\n"
        ">Species_1|100|asm1|ENST001\nATGTTCGCC\n"
        ">Species_2|200|asm2|ENST001\nATGTACGCC\n"
        ">Species_3|300|asm3|ENST001\nATGXXXGCC\n"
        ">Species_4|400|asm4|ENST001\nATG!!!GCC\n"
        ">Species_5|500|asm5|ENST001\nATG***GCC\n"
    )
    f = tmp_path / "test_nt.fasta"
    f.write_text(content)
    return f


class TestAlignmentParser:

    def test_reference_identified(self, aa_fasta, nt_fasta):
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        assert parser.ref_header.startswith("Homo_sapiens")

    def test_coordinate_map_length(self, aa_fasta, nt_fasta):
        """3-residue protein → aa_map should have exactly 3 entries."""
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        assert len(parser.aa_map) == 3

    def test_mutant_codon_extraction(self, aa_fasta, nt_fasta):
        """Inject alt_nt 'A' at CDS position 5 (second base of codon 2, TTC→TAC = Tyr)."""
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        codon = parser.extract_mutant_codon(5, "A")
        assert codon == "TAC"

    def test_mask_X_filtered_from_aa_cdar(self, aa_fasta, nt_fasta):
        """Species_3 has X at position 2 and must not be counted as a c-DAR."""
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        result = parser.check_compensation(2, "F", "Y", 5, "A")
        assert not any("Species_3" in s for s in result["aa_species"]), (
            "Species with X in AA must be filtered from c-DAR species list"
        )

    def test_mask_frameshift_filtered_from_aa_cdar(self, aa_fasta, nt_fasta):
        """Species_4 has ! at position 2 (frameshift) and must not be counted."""
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        result = parser.check_compensation(2, "F", "Y", 5, "A")
        assert not any("Species_4" in s for s in result["aa_species"]), (
            "Frameshifted species (!) must be filtered from c-DAR species list"
        )

    def test_mask_stop_filtered_from_aa_cdar(self, aa_fasta, nt_fasta):
        """Species_5 has * at position 2 (stop codon) and must not be counted."""
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        result = parser.check_compensation(2, "F", "Y", 5, "A")
        assert not any("Species_5" in s for s in result["aa_species"]), (
            "Stop-codon species (*) must be filtered from c-DAR species list"
        )

    def test_valid_cdar_species_counted(self, aa_fasta, nt_fasta):
        """Species_2 has Y at position 2 with codon TAC — must be an AA and NT c-DAR."""
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        result = parser.check_compensation(2, "F", "Y", 5, "A")
        assert result["aa_cdar"] is True
        assert any("Species_2" in s for s in result["aa_species"])
        assert result["nt_cdar"] is True
        assert any("Species_2" in s for s in result["nt_species"])

    def test_check_compensation_returns_mut_codon(self, aa_fasta, nt_fasta):
        """check_compensation must return the mutant codon it built internally."""
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        result = parser.check_compensation(2, "F", "Y", 5, "A")
        assert result["mut_codon"] == "TAC"

    def test_ref_species_not_counted(self, aa_fasta, nt_fasta):
        """Species_1 has the same AA as human ref — not a c-DAR."""
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        result = parser.check_compensation(2, "F", "F", 5, "A")
        assert not any("Homo_sapiens" in s for s in result["aa_species"])

    def test_nt_cdar_is_subset_of_aa_cdar(self, aa_fasta, nt_fasta):
        """NT c-DAR species must be a subset of AA c-DAR species."""
        parser = AlignmentParser(aa_fasta, nt_fasta, "nucDNA")
        result = parser.check_compensation(2, "F", "Y", 5, "A")
        assert set(result["nt_species"]).issubset(set(result["aa_species"])), (
            "NT c-DAR species must be a strict subset of AA c-DAR species"
        )
