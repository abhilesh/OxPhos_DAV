import csv
import json
from collections import Counter
from pathlib import Path

from utils.parsers import (
    ClinvarParser,
    GeneReference,
    MitimpactParser,
    MitomapParser,
    MyVariantParser,
    PhylotreeParser,
    VariantAnnotation,
)
from utils.utils import get_latest

# ==== Path Resolution ====
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
REF_DIR = DATA_DIR / "reference"
CURATED_DIR = DATA_DIR / "annotations" / "curated"
CURATED_DIR.mkdir(parents=True, exist_ok=True)

# Single Source of Truth for coordinates
UNIFIED_COORD_FILE = REF_DIR / "nucdna_gene_coordinates.tsv"

# ==== Configuration ====
APOGEE_LP = 0.716
REVEL_PATHOGENIC_THRESHOLD = 0.75
GNOMAD_MENDELIAN_CUTOFF = 0.01

# Disease relevance filter for nucDNA variants.
# Variants whose disease annotation contains at least one OXPHOS term are kept.
# Variants with only non-OXPHOS disease terms are discarded regardless of tier.
# Generic/unspecified terms ("not provided" etc.) are treated as neutral (kept).
_OXPHOS_DISEASE_TERMS = {
    "mitochondrial complex",
    "leigh syndrome",
    "melas",
    "oxidative phosphorylation",
    "respiratory chain",
    "mitochondrial disease",
    "mitochondrial myopathy",
    "leber",
    "narp",
}
_NON_OXPHOS_DISEASE_TERMS = {
    "pheochromocytoma",
    "paraganglioma",
    "gastrointestinal stromal",
    "hereditary cancer-predisposing",
    "thyroid cancer",
    "renal cell carcinoma",
    "neuroblastoma",
}
_GENERIC_DISEASE_TERMS = {
    "not provided",
    "not specified",
    "inborn genetic diseases",
    "not classified",
}


def _is_oxphos_disease(disease: str) -> bool:
    """Returns True if the disease string is relevant to OXPHOS.

    Logic:
      - Any OXPHOS keyword → keep (True)
      - Generic/unspecified terms only → keep (True, benefit of the doubt)
      - Only non-OXPHOS cancer/tumour terms → discard (False)
    """
    d = disease.lower()
    if any(t in d for t in _OXPHOS_DISEASE_TERMS):
        return True
    if any(t in d for t in _NON_OXPHOS_DISEASE_TERMS):
        return False
    # Generic or unrecognised disease term — keep
    return True


def assign_mtdna_tier(variant: VariantAnnotation) -> str:
    """Implements a strict 3-tier mitochondrial disease decision tree."""
    status = variant.clinical_status.lower()
    apogee = variant.pathogenic_score

    # 1. Confirmed Pathogenic ALWAYS trumps haplogroup status
    if "cfrm [p]" in status or "cfrm [lp]" in status:
        return "Tier A"

    # 2. Haplogroup / Benign Sieves
    if variant.is_haplogroup_marker:
        return "Discarded"
    if "benign" in status or "[lb]" in status or "[b]" in status:
        return "Discarded"

    # 3. Moderate Evidence / Predicted Pathogenic
    is_high_apogee = apogee is not None and apogee >= APOGEE_LP

    if "cfrm [vus*]" in status:
        return "Tier B"
    if ("reported" in status or "conflicting" in status) and is_high_apogee:
        return "Tier B"

    # 4. Low Evidence / Unconfirmed VUS
    return "Tier C"


def assign_nucdna_tier(variant: VariantAnnotation) -> str:
    """Implements a strict 3-tier nuclear disease decision tree mirroring mtDNA."""
    status = variant.clinical_status.lower()
    stars = variant.clinvar_stars or 0
    revel = variant.pathogenic_score or 0.0
    af = variant.population_frequency or 0.0

    # 0. Disease relevance filter — discard variants not linked to OXPHOS disorders
    if not _is_oxphos_disease(variant.disease or ""):
        return "Discarded"

    # 1. Population Sieve and Benign Filters (DISCARD)
    if af > GNOMAD_MENDELIAN_CUTOFF:
        return "Discarded"

    if "benign" in status and "pathogenic" not in status:
        return "Discarded"

    # Identify variants that are Pathogenic/LP without conflicts
    is_pathogenic = "pathogenic" in status and "conflicting" not in status
    is_high_revel = revel >= REVEL_PATHOGENIC_THRESHOLD

    # 2. Strong Evidence / Confirmed (Tier A)
    if is_pathogenic and stars >= 2:
        return "Tier A"

    # 3. Moderate Evidence / Predicted Pathogenic (Tier B)
    if is_pathogenic and stars == 1:
        return "Tier B"

    if (
        "uncertain" in status
        or "vus" in status
        or "conflicting" in status
        or stars == 0
    ):
        if is_high_revel:
            return "Tier B"

    # 4. Low Evidence / Unconfirmed VUS (Tier C)
    return "Tier C"


def save_outputs(annotations: list, genome: str):
    out_file = CURATED_DIR / f"{genome}_annotations_curated.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump([ann.to_dict() for ann in annotations], f, indent=2)


def main():
    print("Locating latest raw data files...")
    hgnc_file = get_latest(DATA_DIR, "Canonical_OXPHOS_Subunits_HGNC*.csv")
    mitomap_file = get_latest(DATA_DIR, "MITOMAP_CodingVariants*.tsv")
    mitimpact_file = get_latest(DATA_DIR, "MitImpact_db*.zip")
    phylotree_file = get_latest(DATA_DIR, "PhyloTree_build_17*.zip")
    clinvar_file = get_latest(DATA_DIR, "ClinVar_VariantSummary*.txt.gz")
    myvariant_file = get_latest(DATA_DIR, "MyVariant_dbNSFP_gnomAD*.json")

    print("Initializing Unified Gene Reference...")
    hgnc_ref = GeneReference(hgnc_file)
    hgnc_ref.load_coordinates(UNIFIED_COORD_FILE)

    print("Loading specialized parsers...")
    mitimpact_parser = MitimpactParser(mitimpact_file)
    phylotree_parser = PhylotreeParser(phylotree_file)
    myvariant_parser = MyVariantParser(myvariant_file)

    mitomap_parser = MitomapParser(hgnc_reference=hgnc_ref)
    clinvar_parser = ClinvarParser(hgnc_reference=hgnc_ref)

    mt_stats = Counter()
    nuc_stats = Counter()

    print("\nProcessing mtDNA (MITOMAP) variants...")
    raw_mito_dicts = mitomap_parser.parse(mitomap_file)
    mt_annotations = []

    for d in raw_mito_dicts:
        variant = VariantAnnotation(
            ann_id=d["nt_change"],
            locus=d["locus"],
            nc_change=d["nt_change"],
            aa_change=d["aa_change"],
            is_synonymous=d["is_synonymous"],
            disease=d["disease"],
            genome="mtDNA",
            reference_assembly="rCRS",
            clinical_status=d["clinical_status"],
            ref_nt=d["ref"],
            alt_nt=d["alt"],
            ref_aa=d.get("ref_aa", ""),
            alt_aa=d.get("alt_aa", ""),
            genomic_pos=d["rCRS_pos"],
        )

        apogee, mitoclass = mitimpact_parser.get_metrics(
            variant.genomic_pos,
            d.get("genomic_ref", variant.ref_nt),
            d.get("genomic_alt", variant.alt_nt),
        )
        variant.pathogenic_score = apogee
        variant.mitoclass = mitoclass
        variant.is_haplogroup_marker = phylotree_parser.is_haplogroup(
            variant.genomic_pos,
            d.get("genomic_ref", variant.ref_nt),
            d.get("genomic_alt", variant.alt_nt),
        )

        variant.tier = assign_mtdna_tier(variant)
        mt_stats[variant.tier] += 1
        mt_annotations.append(variant)

    print("Processing nucDNA (ClinVar) variants...")
    raw_nuc_dicts = clinvar_parser.parse(clinvar_file)
    nuc_annotations = []

    for d in raw_nuc_dicts:
        gene_data = hgnc_ref.get_gene_data(d["locus"])

        variant = VariantAnnotation(
            ann_id=d["allele_id"],
            locus=d["locus"],
            nc_change=d["nt_change"],
            aa_change=d["aa_change"],
            is_synonymous=d["is_synonymous"],
            disease=d["disease"],
            genome="nucDNA",
            reference_assembly="GRCh38",
            clinical_status=d["clinical_status"],
            ref_nt=d["ref"],
            alt_nt=d["alt"],
            ref_aa=d.get("ref_aa", ""),
            alt_aa=d.get("alt_aa", ""),
            genomic_pos=d["grch38_pos"],
            clinvar_stars=d["stars"],
            clinvar_review_status=d["review_status"],
        )

        if gene_data and "strand" in gene_data:
            variant.mitoclass = f"Strand: {gene_data['strand']}"

        # No kwargs here, call sequentially to respect the method signature
        revel, af = myvariant_parser.get_metrics(
            d["chromosome"],
            variant.genomic_pos,
            d.get("genomic_ref", variant.ref_nt),
            d.get("genomic_alt", variant.alt_nt),
        )

        variant.pathogenic_score = revel
        variant.population_frequency = af

        variant.tier = assign_nucdna_tier(variant)
        nuc_stats[variant.tier] += 1
        nuc_annotations.append(variant)

    print("\nSaving curated datasets...")
    save_outputs(mt_annotations, "mtDNA")
    save_outputs(nuc_annotations, "nucDNA")

    print(f"\n{'='*50}")
    print("UNIFIED CURATION SUMMARY")
    print(f"{'='*50}")

    print(f"mtDNA Variants Processed: {len(mt_annotations)}")
    for tier, count in mt_stats.most_common():
        print(f"  {tier:<12}: {count}")

    print(f"\nnucDNA Variants Processed: {len(nuc_annotations)}")
    for tier, count in nuc_stats.most_common():
        print(f"  {tier:<12}: {count}")


if __name__ == "__main__":
    main()
