import csv
import json
from collections import Counter
from pathlib import Path

# Import consolidated utility classes
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
CURATED_DIR = DATA_DIR / "annotations" / "curated"
CURATED_DIR.mkdir(parents=True, exist_ok=True)

# ==== Configuration ====
APOGEE_LP = 0.716
REVEL_PATHOGENIC_THRESHOLD = 0.75
GNOMAD_MENDELIAN_CUTOFF = 0.01


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


def save_outputs(annotations: list, prefix: str):
    """Saves a variant list to JSON and CSV formats."""
    out_json = CURATED_DIR / f"{prefix}_annotations_curated.json"
    out_csv = CURATED_DIR / f"{prefix}_annotations_curated.csv"

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump([ann.to_dict() for ann in annotations], f, indent=2)

    if annotations:
        with open(out_csv, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=annotations[0].to_dict().keys())
            writer.writeheader()
            writer.writerows([ann.to_dict() for ann in annotations])


def main():
    # ==== 1. Dynamically Resolve Input Paths ====
    print("Locating latest raw data files...")
    hgnc_file = get_latest(DATA_DIR, "Canonical_OXPHOS_Subunits_HGNC*.csv")
    mitomap_file = get_latest(DATA_DIR, "MITOMAP_CodingVariants*.tsv")
    mitimpact_file = get_latest(DATA_DIR, "MitImpact_db*.zip")
    phylotree_file = get_latest(DATA_DIR, "PhyloTree_build_17*.zip")
    clinvar_file = get_latest(DATA_DIR, "ClinVar_VariantSummary*.txt.gz")
    myvariant_file = get_latest(DATA_DIR, "MyVariant_dbNSFP_gnomAD*.json")

    # ==== 2. Load Parsers & References ====
    print("Loading databases and reference files into memory...")
    hgnc_ref = GeneReference(hgnc_file)
    mitimpact_parser = MitimpactParser(mitimpact_file)
    phylotree_parser = PhylotreeParser(phylotree_file)
    myvariant_parser = MyVariantParser(myvariant_file)

    mitomap_parser = MitomapParser(hgnc_reference=hgnc_ref)
    clinvar_parser = ClinvarParser(hgnc_reference=hgnc_ref)

    all_annotations = []
    mt_stats = Counter()
    nuc_stats = Counter()

    # ==== 3. Process Mitochondrial DNA ====
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
            genome=d["genome"],
            reference_assembly="rCRS",
            clinical_status=d["clinical_status"],
            ref_nt=d["ref"],
            alt_nt=d["alt"],
            genomic_pos=d["rCRS_pos"],
        )

        apogee, mitoclass = mitimpact_parser.get_metrics(
            variant.genomic_pos, variant.ref_nt, variant.alt_nt
        )
        variant.pathogenic_score = apogee
        variant.mitoclass = mitoclass
        variant.is_haplogroup_marker = phylotree_parser.is_haplogroup(
            variant.genomic_pos, variant.ref_nt, variant.alt_nt
        )

        variant.tier = assign_mtdna_tier(variant)
        mt_stats[variant.tier.split(":")[0]] += 1
        mt_annotations.append(variant)
        all_annotations.append(variant)

    # ==== 4. Process Nuclear DNA ====
    print("Processing nucDNA (ClinVar) variants...")
    raw_nuc_dicts = clinvar_parser.parse(clinvar_file)
    nuc_annotations = []

    for d in raw_nuc_dicts:
        variant = VariantAnnotation(
            ann_id=d["allele_id"],
            locus=d["locus"],
            nc_change=d["nt_change"],
            aa_change=d["aa_change"],
            is_synonymous=d["is_synonymous"],
            disease=d["disease"],
            genome=d["genome"],
            reference_assembly="GRCh38",
            clinical_status=d["clinical_status"],
            ref_nt=d["ref"],
            alt_nt=d["alt"],
            genomic_pos=d["grch38_pos"],
            clinvar_stars=d["stars"],
            clinvar_review_status=d["review_status"],
        )

        revel, af = myvariant_parser.get_metrics(
            chrom=d["chromosome"],
            pos=variant.genomic_pos,
            ref=variant.ref_nt,
            alt=variant.alt_nt,
        )
        variant.pathogenic_score = revel
        variant.population_frequency = af

        variant.tier = assign_nucdna_tier(variant)
        nuc_stats[variant.tier] += 1
        nuc_annotations.append(variant)
        all_annotations.append(variant)

    # ==== 5. Export Outputs ====
    print("\nSaving curated datasets...")
    save_outputs(mt_annotations, "mtDNA")
    save_outputs(nuc_annotations, "nucDNA")
    save_outputs(all_annotations, "all")

    # ==== 6. Print Summary ====
    print(f"\n{'='*50}")
    print("UNIFIED CURATION SUMMARY")
    print(f"{'='*50}")

    print(f"mtDNA Variants Processed: {len(mt_annotations)}")
    for tier, count in mt_stats.most_common():
        print(f"  {tier:<20}: {count}")

    print(f"\nnucDNA Variants Processed: {len(nuc_annotations)}")
    for tier, count in nuc_stats.most_common():
        print(f"  {tier:<20}: {count}")

    print(f"\nTotal Output: {len(all_annotations)} variants.")


if __name__ == "__main__":
    main()
