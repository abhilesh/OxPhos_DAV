import csv
import time
from pathlib import Path
from Bio import Entrez

# ==== Configuration ====
# REQUIRES YOUR EMAIL FOR NCBI API ACCESS
Entrez.email = "abhilesh7@gmail.com"

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"

# Alignments and Tables
NUC_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"
MT_AA_DIR = DATA_DIR / "alignments" / "mtdna_aa"

REF_DIR = DATA_DIR / "reference"
REF_DIR.mkdir(parents=True, exist_ok=True)
OVERVIEW_FILE = REF_DIR / "TOGA_overview_table_hg38.tsv"
OUT_CSV = REF_DIR / "taxid_species_mapping.csv"


def get_toga_taxids_from_tsv() -> dict:
    """Parses the local TOGA overview TSV to extract exact TaxIDs."""
    toga_to_taxid = {}
    if not OVERVIEW_FILE.exists():
        print(f"Error: {OVERVIEW_FILE.name} not found. Run TOGA download script first.")
        return toga_to_taxid

    with open(OVERVIEW_FILE, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species_raw = row.get("Species", "").strip()
            taxid = row.get("Species Taxonomy ID", "").strip()

            if species_raw and taxid:
                # Normalize spaces to underscores to match the FASTA headers
                species_norm = species_raw.replace(" ", "_")
                toga_to_taxid[species_norm] = taxid

    print(f"Loaded {len(toga_to_taxid)} TOGA TaxIDs directly from the overview TSV.")
    return toga_to_taxid


def extract_mt_accessions() -> set:
    """Extracts raw NCBI accessions from mitochondrial FASTA headers."""
    accessions = set()
    for fasta_path in MT_AA_DIR.glob("*_aa_alignment.fasta"):
        with open(fasta_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    acc = line[1:].strip().split("|")[0].strip()
                    accessions.add(acc)
    return accessions


def get_mt_taxids(accessions: list) -> dict:
    """Uses NCBI eSummary to rapidly extract the TaxID embedded in nucleotide records."""
    acc_to_taxid = {}
    batch_size = 100

    print(f"Fetching TaxIDs for {len(accessions)} mitochondrial accessions via NCBI...")
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        try:
            handle = Entrez.esummary(db="nuccore", id=",".join(batch))
            records = Entrez.read(handle)
            for record in records:
                acc = record.get("AccessionVersion", "")
                raw_taxid = record.get("TaxId", "")

                # Force cast to a standard Python integer first to strip the BioPython wrapper
                if acc and raw_taxid:
                    taxid = str(int(raw_taxid))
                    acc_to_taxid[acc] = taxid
            time.sleep(0.34)  # Be polite to NCBI servers (max 3 req/sec)
            print(
                f"  Processed {min(i+batch_size, len(accessions))} / {len(accessions)}"
            )
        except Exception as e:
            print(f"  Error fetching batch starting at {i}: {e}")

    return acc_to_taxid


def get_taxonomy_names(taxids: set) -> dict:
    """Fetches the official scientific names for a set of NCBI TaxIDs."""
    taxid_to_name = {}
    batch_size = 50  # efetch prefers slightly smaller batches
    taxid_list = list(taxids)

    print(
        f"Fetching Scientific Names for {len(taxid_list)} unique TaxIDs via NCBI eFetch..."
    )
    for i in range(0, len(taxid_list), batch_size):
        batch = taxid_list[i : i + batch_size]
        try:
            # Switch to 'efetch' to get the definitive full taxonomy record
            handle = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
            records = Entrez.read(handle)
            for record in records:
                # eFetch reliably returns standard 'TaxId' and 'ScientificName' keys
                raw_id = record.get("TaxId", "")
                name = record.get("ScientificName", "")

                # Strip the biopython wrapper and assign
                if raw_id and name:
                    t_id = str(int(raw_id))
                    taxid_to_name[t_id] = name.replace(" ", "_")
            time.sleep(0.34)
        except Exception as e:
            print(f"  Error fetching taxonomy batch starting at {i}: {e}")

    return taxid_to_name


def main():
    if Entrez.email == "your.email@example.com":
        print("ERROR: Please update Entrez.email with your actual email address.")
        return

    # 1. Load TOGA TaxIDs locally
    toga_to_taxid = get_toga_taxids_from_tsv()
    if not toga_to_taxid:
        return

    # Reverse the TOGA dict to map TaxID -> TOGA Species Name
    taxid_to_toga = {taxid: sp for sp, taxid in toga_to_taxid.items()}

    # 2. Extract and Fetch Mitochondrial TaxIDs
    mt_accessions = list(extract_mt_accessions())
    if not mt_accessions:
        print("No mitochondrial accessions found. Check directory paths.")
        return

    mt_to_taxid = get_mt_taxids(mt_accessions)

    # 3. Fetch NCBI Scientific names for the identified TaxIDs
    unique_taxids = {tid for tid in mt_to_taxid.values() if tid}
    taxid_to_ncbi_name = get_taxonomy_names(unique_taxids)

    # 4. Join on TaxID
    mapping_results = []
    exact_matches = 0
    unmapped = 0

    print("\nMapping datasets via immutable TaxIDs...")
    for acc in mt_accessions:
        # Check standard accession, fallback to versionless if needed
        mt_taxid = mt_to_taxid.get(acc, mt_to_taxid.get(acc.split(".")[0], None))

        ncbi_species = taxid_to_ncbi_name.get(mt_taxid, "UNKNOWN")

        if not mt_taxid:
            mapping_results.append(
                {
                    "ncbi_accession": acc,
                    "taxid": "NOT_FOUND",
                    "ncbi_species": "NOT_FOUND",
                    "toga_species": "",
                    "match_type": "No_NCBI_Record",
                }
            )
            unmapped += 1
            continue

        if mt_taxid in taxid_to_toga:
            mapping_results.append(
                {
                    "ncbi_accession": acc,
                    "taxid": mt_taxid,
                    "ncbi_species": ncbi_species,
                    "toga_species": taxid_to_toga[mt_taxid],
                    "match_type": "Exact_TaxID_Match",
                }
            )
            exact_matches += 1
        else:
            mapping_results.append(
                {
                    "ncbi_accession": acc,
                    "taxid": mt_taxid,
                    "ncbi_species": ncbi_species,
                    "toga_species": "",
                    "match_type": "No_TOGA_Overlap",
                }
            )
            unmapped += 1

    # Write to CSV
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "ncbi_accession",
                "taxid",
                "ncbi_species",
                "toga_species",
                "match_type",
            ],
        )
        writer.writeheader()
        writer.writerows(mapping_results)

    print(f"\nTaxID Mapping Summary:")
    print(f"  Exact Biological Matches: {exact_matches}")
    print(f"  No TOGA Overlap         : {unmapped}")
    print(f"\nSpecies name mapping file written to: {OUT_CSV}")


if __name__ == "__main__":
    main()
