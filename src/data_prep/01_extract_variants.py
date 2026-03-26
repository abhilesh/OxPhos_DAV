import json
from pathlib import Path
from utils.hgnc_parser import GeneReference
from utils.mitomap_parser import MitomapParser
from utils.variant_annotation import VariantAnnotation


def save_annotations_to_json(annotations: list[VariantAnnotation], output_path: Path):
    """Serializes a list of VariantAnnotation objects to a JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dict_list = [ann.to_dict() for ann in annotations]

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(dict_list, f, indent=2)
    print(f"Saved {len(annotations)} mitochondrial annotation objects to {output_path}")


def extract_mtdna_variants():
    # 1. Load HGNC Reference
    hgnc_ref = GeneReference(Path("data/Canonical_OXPHOS_Subunits_HGNC_2026-03-25.csv"))

    # 2. Parse MITOMAP
    print("Parsing MITOMAP dataset...")
    mito_parser = MitomapParser(hgnc_reference=hgnc_ref)

    # Update this path if your raw data file is named differently
    raw_mito_file = Path("data/annotations/raw/MITOMAP_CodingVariants_2026-03-26.tsv")
    raw_mito_dicts = mito_parser.parse(raw_mito_file)

    # 3. Convert to VariantAnnotation objects
    mt_annotations = []
    for d in raw_mito_dicts:
        obj = VariantAnnotation(
            ann_id=d["nt_change"],
            locus=d["locus"],
            nc_change=d["nt_change"],
            aa_change=d["aa_change"],
            is_synonymous=d["is_synonymous"],
            disease=d["disease"],
            genome=d["genome"],
            clinical_status=d["clinical_status"],
            ref_nt=d["ref"],
            alt_nt=d["alt"],
            genomic_pos=d["rCRS_pos"],
        )
        mt_annotations.append(obj)

    # 4. Save intermediate output
    output_target = Path("data/debug/mt_parsed_base.json")
    save_annotations_to_json(mt_annotations, output_target)


if __name__ == "__main__":
    extract_mtdna_variants()
