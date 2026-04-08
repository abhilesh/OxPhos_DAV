#!/usr/bin/env python3
"""
src/data_prep/02_curate_variants.py

Rich-schema data collection for OXPHOS disease-associated variants.
This script assembles all available evidence into VariantRecord objects.
Tier assignment is NOT done here — that is handled by 01_curate_variants.py.

Produces:
  data/annotations/curated/mtDNA_annotations_v2.json
  data/annotations/curated/nucDNA_annotations_v2.json

Each record is a VariantRecord with:
  - Full provenance (source DB, version, record ID)
  - Disaggregated clinical evidence (submitter count, last evaluated, etc.)
  - HGVS c. / p. strings and transcript ID (extracted from Name field)
  - Computational predictors (REVEL, PhyloP, GERP++, AlphaMissense, ESM1b, MPC)
  - gnomAD population data (AF global, popmax, AC, AN, nhomalt)
  - Substitution properties (BLOSUM62, Miyata distance, charge/volume/hydrophobicity)
  - Gene context (complex, subunit role, encoded_by, is_sdh)
  - Reserved structural + cross-species slots (None; populated by later steps)
  - tier: "Unassigned" (populated by classification step)

Run from project root inside the Docker container:
  python src/data_prep/02_curate_variants.py
"""

import json
from pathlib import Path


from utils.parsers import GeneReference, MitimpactParser, PhylotreeParser
from utils.parsers_v2 import ClinvarParserV2, MitomapParserV2, MyVariantParserV2
from utils.utils import get_latest
from utils.variant_record import VariantRecord

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
REF_DIR = DATA_DIR / "reference"
CURATED_DIR = DATA_DIR / "annotations" / "curated"
CURATED_DIR.mkdir(parents=True, exist_ok=True)
UNIFIED_COORD_FILE = REF_DIR / "nucdna_gene_coordinates.tsv"


# ── Saving ─────────────────────────────────────────────────────────────────────


def save_outputs(records: list[VariantRecord], genome: str) -> None:
    out = CURATED_DIR / f"{genome}_annotations_v2.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump([r.to_dict() for r in records], f, indent=2)
    print(f"  Saved {len(records)} records → {out}")


# ── Main ───────────────────────────────────────────────────────────────────────


def main():
    print("Locating raw data files...")
    hgnc_file = get_latest(REF_DIR, "Canonical_OXPHOS_Subunits_HGNC*.csv")
    mitomap_file = get_latest(DATA_DIR, "MITOMAP_CodingVariants*.tsv")
    mitimpact_file = get_latest(DATA_DIR, "MitImpact_db*.zip")
    phylotree_file = get_latest(DATA_DIR, "PhyloTree_build_17*.zip")
    clinvar_file = get_latest(DATA_DIR, "ClinVar_VariantSummary*.txt.gz")
    myvariant_file = get_latest(DATA_DIR, "MyVariant_dbNSFP_gnomAD*.json")

    print("Initializing gene reference...")
    hgnc_ref = GeneReference(hgnc_file)
    hgnc_ref.load_coordinates(UNIFIED_COORD_FILE)

    print("Loading parsers...")
    mitimpact = MitimpactParser(mitimpact_file)
    phylotree = PhylotreeParser(phylotree_file)
    myvariant = MyVariantParserV2(myvariant_file)
    mitomap_p = MitomapParserV2(hgnc_ref)
    clinvar_p = ClinvarParserV2(hgnc_ref)

    # ── mtDNA ──────────────────────────────────────────────────────────────────
    print("\nProcessing mtDNA (MITOMAP) variants...")
    mt_records: list[VariantRecord] = []

    for d in mitomap_p.parse(mitomap_file):
        pos = d["rCRS_pos"]
        apogee, mitoclass = mitimpact.get_metrics(
            pos, d.get("genomic_ref", d["ref"]), d.get("genomic_alt", d["alt"])
        )
        is_haplo = phylotree.is_haplogroup(
            pos, d.get("genomic_ref", d["ref"]), d.get("genomic_alt", d["alt"])
        )

        rec = VariantRecord(
            # Identity
            ann_id=d["nt_change"],
            source_db=d["source_db"],
            source_db_version=d["source_db_version"],
            source_record_id=d["source_record_id"],
            genome="mtDNA",
            reference_assembly="rCRS",
            # Coordinates
            locus=d["locus"],
            nc_change=d["nt_change"],
            aa_change=d["aa_change"],
            ref_nt=d["ref"],
            alt_nt=d["alt"],
            ref_aa=d.get("ref_aa", ""),
            alt_aa=d.get("alt_aa", ""),
            genomic_pos=pos,
            is_synonymous=d["is_synonymous"],
            hgvs_c=d.get("hgvs_c", ""),
            # Clinical
            disease=d["disease"],
            disease_terms=[t.strip() for t in d["disease"].split("|") if t.strip()],
            clinical_status=d["clinical_status"],
            mitomap_plasmy=d.get("mitomap_plasmy", ""),
            mitomap_pubmed_count=d.get("mitomap_pubmed_count", 0),
            is_haplogroup_marker=is_haplo,
            # Predictors
            apogee2_score=apogee,
            mitoclass=mitoclass or "",
        )
        rec.populate_gene_context()
        rec.populate_substitution_properties()
        mt_records.append(rec)

    # ── nucDNA ─────────────────────────────────────────────────────────────────
    print("Processing nucDNA (ClinVar) variants...")
    nuc_records: list[VariantRecord] = []

    for d in clinvar_p.parse(clinvar_file):
        gene_data = hgnc_ref.get_gene_data(d["locus"])
        metrics = myvariant.get_all_metrics(
            d["chromosome"],
            d["grch38_pos"],
            d.get("genomic_ref", d["ref"]),
            d.get("genomic_alt", d["alt"]),
        )

        rec = VariantRecord(
            # Identity
            ann_id=d["allele_id"],
            source_db=d["source_db"],
            source_db_version=d["source_db_version"],
            source_record_id=d["source_record_id"],
            genome="nucDNA",
            reference_assembly="GRCh38",
            # Coordinates
            locus=d["locus"],
            nc_change=d["nt_change"],
            aa_change=d["aa_change"],
            ref_nt=d["ref"],
            alt_nt=d["alt"],
            ref_aa=d.get("ref_aa", ""),
            alt_aa=d.get("alt_aa", ""),
            genomic_pos=d["grch38_pos"],
            is_synonymous=d["is_synonymous"],
            hgvs_c=d.get("hgvs_c", ""),
            hgvs_p=d.get("hgvs_p", ""),
            transcript_id=d.get("transcript_id", ""),
            # Clinical
            disease=d["disease"],
            disease_terms=d.get("disease_terms", []),
            clinical_status=d["clinical_status"],
            clinvar_stars=d["stars"],
            clinvar_review_status=d["review_status"],
            clinvar_submitters_n=d.get("clinvar_submitters_n", 0),
            clinvar_last_evaluated=d.get("clinvar_last_evaluated", ""),
            clinvar_conflicting=d.get("clinvar_conflicting", False),
            # Computational predictors
            revel_score=metrics.get("revel_score"),
            phylop_100vert=metrics.get("phylop_100vert"),
            gerp_rs=metrics.get("gerp_rs"),
            alphamissense_score=metrics.get("alphamissense_score"),
            alphamissense_class=metrics.get("alphamissense_class"),
            esm1b_score=metrics.get("esm1b_score"),
            mpc_score=metrics.get("mpc_score"),
            # Population
            gnomad_af_global=metrics.get("gnomad_af_global"),
            gnomad_af_popmax=metrics.get("gnomad_af_popmax"),
            gnomad_popmax_pop=metrics.get("gnomad_popmax_pop"),
            gnomad_ac=metrics.get("gnomad_ac"),
            gnomad_an=metrics.get("gnomad_an"),
            gnomad_nhomalt=metrics.get("gnomad_nhomalt"),
        )
        rec.populate_gene_context()
        rec.populate_substitution_properties()
        nuc_records.append(rec)

    # ── Save ───────────────────────────────────────────────────────────────────
    print("\nSaving curated datasets...")
    save_outputs(mt_records, "mtDNA")
    save_outputs(nuc_records, "nucDNA")

    # ── Summary ────────────────────────────────────────────────────────────────
    def _cov(records, field):
        return sum(1 for r in records if getattr(r, field) is not None)

    print(f"\n{'='*55}")
    print("DATA COLLECTION SUMMARY (v2 rich schema)")
    print(f"{'='*55}")
    print(f"  mtDNA variants  : {len(mt_records)}")
    print(f"  nucDNA variants : {len(nuc_records)}")
    print(f"\nField coverage (nucDNA, n={len(nuc_records)}):")
    for fld in [
        "hgvs_c",
        "transcript_id",
        "revel_score",
        "alphamissense_score",
        "alphamissense_class",
        "esm1b_score",
        "mpc_score",
        "phylop_100vert",
        "gerp_rs",
        "gnomad_af_global",
        "gnomad_af_popmax",
        "gnomad_popmax_pop",
        "gnomad_ac",
        "gnomad_an",
        "gnomad_nhomalt",
        "miyata_distance",
        "blosum62",
        "complex_id",
    ]:
        n = _cov(nuc_records, fld)
        print(f"  {fld:<30}: {n:>5}  ({100*n/len(nuc_records):.0f}%)")

    print(f"\nField coverage (mtDNA, n={len(mt_records)}):")
    for fld in [
        "apogee2_score",
        "mitoclass",
        "miyata_distance",
        "blosum62",
        "complex_id",
    ]:
        n = _cov(mt_records, fld)
        print(f"  {fld:<30}: {n:>5}  ({100*n/len(mt_records):.0f}%)")


if __name__ == "__main__":
    main()
