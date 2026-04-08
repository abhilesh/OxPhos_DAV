"""
src/utils/variant_record.py

Rich evidence schema for OXPHOS disease-associated variants.

Design principles:
  - Store everything you might want to condition on downstream.
  - Structural and cross-species fields are None at curation time; they are
    populated by later pipeline steps (structural mapping, classify scripts).
  - All fields have explicit types; Optional[T] = None means "not yet computed."
  - to_dict() / from_dict() round-trip cleanly through JSON.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional


# Load physicochemical properties from local AAindex cache
_ROOT = Path(__file__).resolve().parents[2]
_REF_DIR = _ROOT / "data" / "reference"

_aaindex_files = sorted(_REF_DIR.glob("aaindex_properties_*.json"))
if not _aaindex_files:
    raise FileNotFoundError(f"Could not locate aaindex_properties_*.json in {_REF_DIR}")

with open(_aaindex_files[-1], "r", encoding="utf-8") as _f:
    _AA_DATA = json.load(_f)

_KD_HYDRO = _AA_DATA["indices"]["hydrophobicity_kd"]
_VOLUME = _AA_DATA["indices"]["volume"]
_BLOSUM62 = _AA_DATA["matrices"]["blosum62"]
_MIYATA = _AA_DATA["matrices"]["miyata_distance"]


# Formal charge of residues (pH ~7) remains hard-coded as it is a universal constant
_CHARGE: dict[str, int] = {
    "R": 1,
    "K": 1,
    "H": 0,  # positive (H treated as neutral at pH 7)
    "D": -1,
    "E": -1,  # negative
}


def blosum62(ref: str, alt: str) -> Optional[int]:
    try:
        return int(_BLOSUM62[ref][alt])
    except KeyError:
        return None


def miyata_distance(ref: str, alt: str) -> Optional[float]:
    try:
        return float(_MIYATA[ref][alt])
    except KeyError:
        return None


def charge_change(ref: str, alt: str) -> int:
    return _CHARGE.get(alt, 0) - _CHARGE.get(ref, 0)


def volume_change(ref: str, alt: str) -> Optional[float]:
    try:
        return round(float(_VOLUME[alt]) - float(_VOLUME[ref]), 1)
    except KeyError:
        return None


def hydrophobicity_change(ref: str, alt: str) -> Optional[float]:
    try:
        return round(float(_KD_HYDRO[alt]) - float(_KD_HYDRO[ref]), 2)
    except KeyError:
        return None


# Gene context constants
_COMPLEX_MAP: dict[str, str] = {}
for _g in [
    "MT-ND1",
    "MT-ND2",
    "MT-ND3",
    "MT-ND4",
    "MT-ND4L",
    "MT-ND5",
    "MT-ND6",
    "NDUFA1",
    "NDUFA2",
    "NDUFA3",
    "NDUFA5",
    "NDUFA6",
    "NDUFA7",
    "NDUFA8",
    "NDUFA9",
    "NDUFA10",
    "NDUFA11",
    "NDUFA12",
    "NDUFA13",
    "NDUFAB1",
    "NDUFB1",
    "NDUFB2",
    "NDUFB3",
    "NDUFB4",
    "NDUFB5",
    "NDUFB6",
    "NDUFB7",
    "NDUFB8",
    "NDUFB9",
    "NDUFB10",
    "NDUFB11",
    "NDUFC1",
    "NDUFC2",
    "NDUFS1",
    "NDUFS2",
    "NDUFS3",
    "NDUFS4",
    "NDUFS5",
    "NDUFS6",
    "NDUFS7",
    "NDUFS8",
    "NDUFV1",
    "NDUFV2",
    "NDUFV3",
]:
    _COMPLEX_MAP[_g] = "CI"
for _g in ["SDHA", "SDHB", "SDHC", "SDHD"]:
    _COMPLEX_MAP[_g] = "CII"
for _g in [
    "MT-CYB",
    "UQCRB",
    "UQCRC1",
    "UQCRC2",
    "UQCRFS1",
    "UQCRH",
    "UQCRQ",
    "UQCR10",
    "UQCR11",
    "CYC1",
]:
    _COMPLEX_MAP[_g] = "CIII"
for _g in [
    "MT-CO1",
    "MT-CO2",
    "MT-CO3",
    "COX4I1",
    "COX4I2",
    "COX5A",
    "COX5B",
    "COX6A1",
    "COX6A2",
    "COX6B1",
    "COX6C",
    "COX7A1",
    "COX7A2",
    "COX7B",
    "COX7C",
    "COX8A",
    "NDUFA4",
]:
    _COMPLEX_MAP[_g] = "CIV"
for _g in [
    "MT-ATP6",
    "MT-ATP8",
    "ATP5F1A",
    "ATP5F1B",
    "ATP5F1C",
    "ATP5F1D",
    "ATP5F1E",
    "ATP5MC1",
    "ATP5MC2",
    "ATP5MC3",
    "ATP5PB",
    "ATP5PD",
    "ATP5PF",
    "ATP5PO",
]:
    _COMPLEX_MAP[_g] = "CV"

# Subunit roles: core catalytic vs. accessory vs. assembly_factor
# Core = directly in the catalytic/electron-transfer core; accessory = structural/regulatory
# Source: MitoCarta3 annotations + literature
_SUBUNIT_ROLE: dict[str, str] = {
    # CI core (Zickermann et al., Fiedorczuk et al.)
    "MT-ND1": "core",
    "MT-ND2": "core",
    "MT-ND3": "core",
    "MT-ND4": "core",
    "MT-ND4L": "core",
    "MT-ND5": "core",
    "MT-ND6": "core",
    "NDUFS1": "core",
    "NDUFS2": "core",
    "NDUFS3": "core",
    "NDUFS7": "core",
    "NDUFS8": "core",
    "NDUFV1": "core",
    "NDUFV2": "core",
    # CI accessory
    "NDUFA1": "accessory",
    "NDUFA2": "accessory",
    "NDUFA3": "accessory",
    "NDUFA5": "accessory",
    "NDUFA6": "accessory",
    "NDUFA7": "accessory",
    "NDUFA8": "accessory",
    "NDUFA9": "accessory",
    "NDUFA10": "accessory",
    "NDUFA11": "accessory",
    "NDUFA12": "accessory",
    "NDUFA13": "accessory",
    "NDUFAB1": "accessory",
    "NDUFB1": "accessory",
    "NDUFB2": "accessory",
    "NDUFB3": "accessory",
    "NDUFB4": "accessory",
    "NDUFB5": "accessory",
    "NDUFB6": "accessory",
    "NDUFB7": "accessory",
    "NDUFB8": "accessory",
    "NDUFB9": "accessory",
    "NDUFB10": "accessory",
    "NDUFB11": "accessory",
    "NDUFC1": "accessory",
    "NDUFC2": "accessory",
    "NDUFS4": "accessory",
    "NDUFS5": "accessory",
    "NDUFS6": "accessory",
    "NDUFV3": "accessory",
    # CII all four are core catalytic
    "SDHA": "core",
    "SDHB": "core",
    "SDHC": "core",
    "SDHD": "core",
    # CIII core
    "MT-CYB": "core",
    "UQCRFS1": "core",
    "CYC1": "core",
    "UQCRC1": "core",
    "UQCRC2": "core",
    # CIII accessory
    "UQCRB": "accessory",
    "UQCRH": "accessory",
    "UQCRQ": "accessory",
    "UQCR10": "accessory",
    "UQCR11": "accessory",
    # CIV core
    "MT-CO1": "core",
    "MT-CO2": "core",
    "MT-CO3": "core",
    # CIV accessory
    "COX4I1": "accessory",
    "COX4I2": "accessory",
    "COX5A": "accessory",
    "COX5B": "accessory",
    "COX6A1": "accessory",
    "COX6A2": "accessory",
    "COX6B1": "accessory",
    "COX6C": "accessory",
    "COX7A1": "accessory",
    "COX7A2": "accessory",
    "COX7B": "accessory",
    "COX7C": "accessory",
    "COX8A": "accessory",
    "NDUFA4": "accessory",
    # CV core
    "MT-ATP6": "core",
    "ATP5F1A": "core",
    "ATP5F1B": "core",
    "ATP5F1C": "core",
    "ATP5F1D": "core",
    # CV accessory
    "MT-ATP8": "accessory",
    "ATP5F1E": "accessory",
    "ATP5MC1": "accessory",
    "ATP5MC2": "accessory",
    "ATP5MC3": "accessory",
    "ATP5PB": "accessory",
    "ATP5PD": "accessory",
    "ATP5PF": "accessory",
    "ATP5PO": "accessory",
}

_SDH_GENES = {"SDHA", "SDHB", "SDHC", "SDHD"}


def gene_complex(gene: str) -> Optional[str]:
    return _COMPLEX_MAP.get(gene)


def subunit_role(gene: str) -> str:
    return _SUBUNIT_ROLE.get(gene, "unknown")


def encoded_by(gene: str) -> str:
    return "mitochondrial" if gene.startswith("MT-") else "nuclear"


# Main dataclass


@dataclass
class VariantRecord:
    """
    Full evidence schema for one disease-associated variant.

    Fields are grouped by population stage:
      - Identity / provenance       : always populated at curation
      - Gene context                : always populated at curation
      - Clinical evidence           : populated at curation
      - Computational predictors    : populated at curation (from MyVariant/MitImpact)
      - Population genetics         : populated at curation
      - Substitution properties     : computed at curation from ref_aa / alt_aa
      - Structural context          : None -> populated by structural mapping step
      - Cross-species               : None -> populated by classification step
      - Tier assignment             : assigned at curation
    """

    # Identity / provenance
    ann_id: str  # unique key (nt_change for mtDNA, VariationID for nucDNA)
    source_db: str  # "MITOMAP" | "ClinVar"
    source_db_version: str  # file date, e.g. "2026-03-26"
    source_record_id: str  # MITOMAP id or ClinVar VariationID
    genome: str  # "mtDNA" | "nucDNA"
    reference_assembly: str  # "rCRS" | "GRCh38"

    # Variant coordinates
    locus: str
    nc_change: str  # HGVS-like coding change string
    aa_change: str  # e.g. "S45F"
    ref_nt: str
    alt_nt: str
    ref_aa: str
    alt_aa: str
    genomic_pos: int
    is_synonymous: bool

    # HGVS canonical strings (CDS-relative; from Name field for ClinVar)
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    transcript_id: Optional[str] = None  # transcript used for annotation

    # Gene context
    complex_id: Optional[str] = None  # CI / CII / CIII / CIV / CV
    subunit_role: str = "unknown"  # core / accessory
    encoded_by: str = "unknown"  # nuclear / mitochondrial
    is_sdh: bool = False

    # Clinical evidence
    disease: str = ""
    disease_terms: list = field(default_factory=list)  # pipe-split, deduplicated
    clinical_status: Optional[str] = None

    # ClinVar-specific
    clinvar_stars: Optional[int] = None
    clinvar_review_status: Optional[str] = None
    clinvar_submitters_n: Optional[int] = None
    clinvar_last_evaluated: Optional[str] = None
    clinvar_conflicting: Optional[bool] = None

    # MITOMAP-specific
    mitomap_plasmy: Optional[str] = None  # "homo" | "hetero" | "homo/hetero"
    mitomap_pubmed_count: Optional[int] = None
    is_haplogroup_marker: Optional[bool] = None
    # Computational predictors - conservation-dependent
    revel_score: Optional[float] = None
    phylop_100vert: Optional[float] = None
    gerp_rs: Optional[float] = None

    # Conservation-independent / structural
    alphamissense_score: Optional[float] = None
    alphamissense_class: Optional[str] = None
    esm1b_score: Optional[float] = None
    mpc_score: Optional[float] = None

    # mtDNA-specific
    apogee2_score: Optional[float] = None
    mitoclass: Optional[str] = None

    # Population genetics
    gnomad_af_global: Optional[float] = None
    gnomad_af_popmax: Optional[float] = None
    gnomad_popmax_pop: Optional[str] = None
    gnomad_ac: Optional[int] = None
    gnomad_an: Optional[int] = None
    gnomad_nhomalt: Optional[int] = None  # homozygote count - key for recessive

    # Substitution properties
    blosum62: Optional[int] = None
    miyata_distance: Optional[float] = None
    charge_change: Optional[int] = None
    volume_change: Optional[float] = None
    hydrophobicity_change: Optional[float] = None
    is_proline_involved: bool = False
    is_glycine_involved: bool = False
    is_cysteine_involved: bool = False

    # Structural context (populated by structural mapping step)
    pdb_id: Optional[str] = None
    pdb_chain: Optional[str] = None
    pdb_resnum: Optional[int] = None
    structure_resolved: Optional[bool] = None
    secondary_structure: Optional[str] = None  # H / E / C
    rsa: Optional[float] = None  # relative solvent accessibility
    is_buried: Optional[bool] = None  # rsa < 0.2
    is_interface: Optional[bool] = None
    is_mito_nuclear_interface: Optional[bool] = None
    distance_to_cofactor: Optional[float] = None
    is_cofactor_binding: Optional[bool] = None

    # Cross-species (populated by classification step)
    n_species_aligned: Optional[int] = None
    n_species_with_disease_allele: Optional[int] = None
    lineages_with_disease_allele: list = field(default_factory=list)
    is_cdav_nucleotide: Optional[bool] = None
    is_cdav_amino_acid: Optional[bool] = None

    # Tier assignment
    tier: str = "Unassigned"

    # Helpers

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "VariantRecord":
        known = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in d.items() if k in known})

    def populate_substitution_properties(self) -> None:
        """Compute substitution properties from ref_aa / alt_aa in-place."""
        if self.ref_aa and self.alt_aa:
            self.blosum62 = blosum62(self.ref_aa, self.alt_aa)
            self.miyata_distance = miyata_distance(self.ref_aa, self.alt_aa)
            self.charge_change = charge_change(self.ref_aa, self.alt_aa)
            self.volume_change = volume_change(self.ref_aa, self.alt_aa)
            self.hydrophobicity_change = hydrophobicity_change(self.ref_aa, self.alt_aa)
            aas = {self.ref_aa, self.alt_aa}
            self.is_proline_involved = "P" in aas
            self.is_glycine_involved = "G" in aas
            self.is_cysteine_involved = "C" in aas

    def populate_gene_context(self) -> None:
        """Fill complex_id, subunit_role, encoded_by, is_sdh from gene name."""
        base_gene = self.locus.split("/")[0]
        self.complex_id = gene_complex(base_gene)
        self.subunit_role = subunit_role(base_gene)
        self.encoded_by = encoded_by(base_gene)
        self.is_sdh = base_gene in _SDH_GENES
