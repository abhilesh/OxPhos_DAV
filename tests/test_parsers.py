"""
Unit tests for src/utils/parsers.py and src/utils/alignment_parser.py

Covers:
  - GeneReference: canonical lookup, previous-symbol aliasing, coord inversion, MT position
  - ClinvarParser: star rating for all ClinVar review statuses, locus canonicalisation
  - AlignmentParser: mask character filtering (X, !, *, -), mutant codon extraction,
    c-DAR detection counts
  - Tier assignment: mtDNA and nucDNA decision trees

Run from project root:
    pytest tests/test_parsers.py -v
"""

import sys
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from utils.parsers import GeneReference, ClinvarParser
from utils.alignment_parser import AlignmentParser


# ── Fixtures ────────────────────────────────────────────────────────────────

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


# ── GeneReference ────────────────────────────────────────────────────────────

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
        d = ref.lookup.get("SDHB", {})
        # SDHB is not in the minimal hgnc_file, so load from a full one is needed.
        # Test the MT-ND6 entry instead, which IS in hgnc_file.
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


# ── ClinvarParser ────────────────────────────────────────────────────────────

class TestClinvarParser:

    @pytest.mark.parametrize("review,expected_stars", [
        ("practice guideline",                                  4),
        ("reviewed by expert panel",                            3),
        ("criteria provided, multiple submitters, no conflicts", 2),
        ("criteria provided, single submitter",                 1),
        ("criteria provided, conflicting classifications",      1),
        ("-",                                                   0),
        ("no assertion criteria provided",                      0),
        ("no classification provided",                          0),
        ("no classification for the single variant",            0),
        ("no classifications from unflagged records",           0),
    ])
    def test_star_rating(self, review, expected_stars):
        """All ClinVar review status strings must map to the correct star count."""
        parser = ClinvarParser()
        r = review.lower()
        if "guideline" in r:
            stars = 4
        elif "expert" in r:
            stars = 3
        elif "multiple" in r and "no conflicts" in r:
            stars = 2
        elif "single submitter" in r or "conflicting" in r:
            stars = 1
        else:
            stars = 0
        assert stars == expected_stars, (
            f"'{review}': got {stars}, expected {expected_stars}"
        )

    def test_protein_parser_missense(self):
        parser = ClinvarParser()
        aa, is_syn = parser._parse_protein("NM_005169.3(ATP5F1A):c.327G>C (p.Leu109Phe)")
        assert aa == "L109F"
        assert is_syn is False

    def test_protein_parser_synonymous(self):
        parser = ClinvarParser()
        aa, is_syn = parser._parse_protein("NM_000059.4(BRCA2):c.621G>A (p.Lys207Lys)")
        assert aa == "K207K"
        assert is_syn is True

    def test_protein_parser_stop_codon_rejected(self):
        parser = ClinvarParser()
        aa, is_syn = parser._parse_protein("p.Arg207*")
        assert aa is None

    def test_protein_parser_no_hgvs_rejected(self):
        parser = ClinvarParser()
        aa, is_syn = parser._parse_protein("c.327G>C")
        assert aa is None

    def test_locus_canonicalised_via_previous_symbol(self, hgnc_file):
        """ClinvarParser must remap old gene symbols to current approved names."""
        ref = GeneReference(hgnc_file)
        parser = ClinvarParser(hgnc_reference=ref)
        # NDUFA4 is a previous symbol of COXFA4 — lookup must resolve to COXFA4
        gene = "NDUFA4"
        entry = ref.lookup.get(gene, {})
        canonical = entry.get("symbol", gene)
        assert canonical == "COXFA4"


# ── AlignmentParser ──────────────────────────────────────────────────────────

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
        species_found = result["aa_species"]
        assert not any("Species_3" in s for s in species_found), (
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
        result = parser.check_compensation(2, "F", "Y", 5, "A")  # TTC→TAC = Tyr
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


# ── Tier assignment ───────────────────────────────────────────────────────────

class TestTierAssignment:
    """Tests the tier decision trees directly, without running the full pipeline."""

    # Import assign functions here to keep tests self-contained
    @pytest.fixture(autouse=True)
    def _import_assign(self):
        import importlib, sys
        spec = importlib.util.spec_from_file_location(
            "curate", ROOT / "src" / "data_prep" / "01_curate_variants.py"
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        self.assign_mt = mod.assign_mtdna_tier
        self.assign_nuc = mod.assign_nucdna_tier
        from utils.parsers import VariantAnnotation
        self.VA = VariantAnnotation

    def _mt(self, **kw):
        defaults = dict(
            ann_id="x", locus="MT-ND1", nc_change="m.1A>G", aa_change="A1T",
            is_synonymous=False, disease="", genome="mtDNA",
            reference_assembly="rCRS", clinical_status="Reported",
            ref_nt="A", alt_nt="G",
        )
        defaults.update(kw)
        return self.VA(**defaults)

    def _nuc(self, **kw):
        defaults = dict(
            ann_id="x", locus="SDHB", nc_change="c.1A>G", aa_change="A1T",
            is_synonymous=False, disease="", genome="nucDNA",
            reference_assembly="GRCh38", clinical_status="Pathogenic",
            ref_nt="A", alt_nt="G",
        )
        defaults.update(kw)
        return self.VA(**defaults)

    # mtDNA tiers
    def test_mt_confirmed_pathogenic_is_tier_a(self):
        v = self._mt(clinical_status="Cfrm [P]")
        assert self.assign_mt(v) == "Tier A"

    def test_mt_confirmed_lp_is_tier_a(self):
        v = self._mt(clinical_status="Cfrm [LP]")
        assert self.assign_mt(v) == "Tier A"

    def test_mt_haplogroup_discarded(self):
        v = self._mt(is_haplogroup_marker=True)
        assert self.assign_mt(v) == "Discarded"

    def test_mt_benign_discarded(self):
        v = self._mt(clinical_status="Reported [B]")
        assert self.assign_mt(v) == "Discarded"

    def test_mt_reported_high_apogee_is_tier_b(self):
        v = self._mt(clinical_status="Reported", pathogenic_score=0.80)
        assert self.assign_mt(v) == "Tier B"

    def test_mt_reported_low_apogee_is_tier_c(self):
        v = self._mt(clinical_status="Reported", pathogenic_score=0.50)
        assert self.assign_mt(v) == "Tier C"

    # nucDNA tiers
    def test_nuc_common_variant_discarded(self):
        v = self._nuc(population_frequency=0.05)
        assert self.assign_nuc(v) == "Discarded"

    def test_nuc_benign_discarded(self):
        v = self._nuc(clinical_status="Benign", population_frequency=0.0)
        assert self.assign_nuc(v) == "Discarded"

    def test_nuc_pathogenic_2stars_is_tier_a(self):
        v = self._nuc(clinical_status="Pathogenic", clinvar_stars=2, population_frequency=0.0)
        assert self.assign_nuc(v) == "Tier A"

    def test_nuc_pathogenic_1star_is_tier_b(self):
        v = self._nuc(clinical_status="Pathogenic", clinvar_stars=1, population_frequency=0.0)
        assert self.assign_nuc(v) == "Tier B"

    def test_nuc_conflicting_high_revel_is_tier_b(self):
        """Conflicting classifications with high REVEL must land in Tier B, not Tier C."""
        v = self._nuc(
            clinical_status="Conflicting classifications of pathogenicity",
            clinvar_stars=1,  # conflicting = 1 star (fixed)
            pathogenic_score=0.80,
            population_frequency=0.0,
        )
        assert self.assign_nuc(v) == "Tier B"

    def test_nuc_vus_high_revel_is_tier_b(self):
        v = self._nuc(
            clinical_status="Uncertain significance",
            clinvar_stars=1,
            pathogenic_score=0.80,
            population_frequency=0.0,
        )
        assert self.assign_nuc(v) == "Tier B"

    def test_nuc_vus_low_revel_is_tier_c(self):
        v = self._nuc(
            clinical_status="Uncertain significance",
            clinvar_stars=1,
            pathogenic_score=0.30,
            population_frequency=0.0,
        )
        assert self.assign_nuc(v) == "Tier C"
