from __future__ import annotations

from .download_paths import RAW_ANNOTATIONS_DIR, RAW_REFERENCE_DIR, TODAY


RESOURCE_CONFIG = {
    "hgnc_oxphos_gene_list": {
        "url": "https://www.genenames.org/cgi-bin/genegroup/download?id=639&type=branch",
        "target": RAW_REFERENCE_DIR / f"Canonical_OXPHOS_Subunits_HGNC_{TODAY}.csv",
        "pattern": "Canonical_OXPHOS_Subunits_HGNC_*.csv",
        "file_kind": "tsv",
    },
    "clinvar_variant_summary": {
        "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
        "target": RAW_ANNOTATIONS_DIR / f"ClinVar_VariantSummary_{TODAY}.txt.gz",
        "pattern": "ClinVar_VariantSummary_*.txt.gz",
        "file_kind": "gzip_tsv",
    },
    "mitomap_coding_variants": {
        "url": "https://www.mitomap.org/cgi-bin/disease.cgi",
        "target": RAW_ANNOTATIONS_DIR / f"MITOMAP_CodingVariants_{TODAY}.tsv",
        "pattern": "MITOMAP_CodingVariants_*.tsv",
        "file_kind": "tsv",
    },
    "mitimpact_db": {
        "url": "https://mitimpact.css-mendel.it/cdn/MitImpact_db_3.1.3.txt.zip",
        "target": RAW_ANNOTATIONS_DIR / f"MitImpact_db_{TODAY}.txt.zip",
        "pattern": "MitImpact_db_*.txt.zip",
        "file_kind": "zip",
        "version_hint": "3.1.3",
    },
    "myvariant_dbnsfp_gnomad": {
        "url": "https://myvariant.info/v1/query",
        "target": RAW_ANNOTATIONS_DIR / f"MyVariant_dbNSFP_gnomAD_{TODAY}.json",
        "pattern": "MyVariant_dbNSFP_gnomAD_*.json",
        "file_kind": "json",
    },
    "phylotree_build17": {
        "url": "https://www.phylotree.org/builds/mtDNA_tree_Build_17%20-%20rCRS-oriented%20version.zip",
        "target": RAW_REFERENCE_DIR / f"PhyloTree_build_17_{TODAY}.zip",
        "pattern": "PhyloTree_build_17_*.zip",
        "file_kind": "zip",
        "version_hint": "17",
    },
    "mane_grch38_summary": {
        "url": "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz",
        "target": RAW_REFERENCE_DIR / f"MANE_GRCh38_v1.5_{TODAY}.txt.gz",
        "pattern": "MANE_GRCh38_v1.5_*.txt.gz",
        "fallback_patterns": ["MANE_GRCh38_v1.5.txt.gz"],
        "file_kind": "gzip_tsv",
        "version_hint": "v1.5",
    },
    "toga_overview_hg38": {
        "url": "https://genome.senckenberg.de/download/TOGA/human_hg38_reference/overview.table.tsv",
        "target": RAW_REFERENCE_DIR / f"TOGA_overview_table_hg38_{TODAY}.tsv",
        "pattern": "TOGA_overview_table_hg38_*.tsv",
        "file_kind": "tsv",
    },
}
