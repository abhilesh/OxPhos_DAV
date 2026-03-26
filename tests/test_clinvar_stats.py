import pandas as pd
from collections import Counter
from pathlib import Path


def generate_clinvar_stats(file_path: Path, target_genes: set = None):
    """Parses ClinVar summary in chunks and generates statistics."""

    print("Parsing ClinVar variant_summary.txt.gz...")

    # Only load the columns we actually need for stats to save memory
    usecols = [
        "Assembly",
        "Chromosome",
        "Type",
        "GeneSymbol",
        "ClinicalSignificance",
        "ReviewStatus",
    ]

    total_snvs = 0
    gene_counts = Counter()
    sig_counts = Counter()
    review_counts = Counter()

    # Process in chunks of 250,000 rows
    chunk_iterator = pd.read_csv(
        file_path,
        sep="\t",
        compression="gzip",
        usecols=usecols,
        chunksize=250000,
        low_memory=False,
    )

    for chunk in chunk_iterator:
        # Filter for GRCh38, single nucleotide variants, and exclude mtDNA
        mask = (
            (chunk["Assembly"] == "GRCh38")
            & (chunk["Type"] == "single nucleotide variant")
            & (chunk["Chromosome"] != "MT")
        )

        # Apply target gene filter if a set is provided
        if target_genes:
            mask = mask & (chunk["GeneSymbol"].isin(target_genes))

        filtered = chunk[mask]
        total_snvs += len(filtered)

        # Update cumulative counters
        gene_counts.update(filtered["GeneSymbol"].dropna())
        sig_counts.update(filtered["ClinicalSignificance"].dropna())
        review_counts.update(filtered["ReviewStatus"].dropna())

    # --- Print Statistics ---
    print("=" * 50)
    print("CLINVAR RAW DATA STATISTICS (GRCh38 SNVs Only)")
    print("=" * 50)
    print(f"Total Target SNVs Found: {total_snvs}\n")

    print("--- Top 15 Genes by Variant Count ---")
    for gene, count in gene_counts.most_common(15):
        print(f"{gene:<15}: {count}")
    print()

    print("--- By Clinical Significance ---")
    for sig, count in sig_counts.most_common():
        print(f"{sig:<35}: {count}")
    print()

    print("--- By Review Status (Star Rating Proxy) ---")
    for rev, count in review_counts.most_common():
        print(f"{rev:<55}: {count}")
    print()


# Usage (Assuming you have extracted your HGNC target genes into a Python set):
# target_hgnc_genes = {"SDHA", "NDUFS1", "UQCRC1", ...}
# generate_clinvar_stats("data/annotations/raw/ClinVar_VariantSummary_2024-05-20.txt.gz", target_hgnc_genes)
