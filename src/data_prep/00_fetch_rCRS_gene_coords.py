import csv
from pathlib import Path
from Bio import Entrez, SeqIO

# ==== Configuration ====
# NCBI requires an email to track API usage
Entrez.email = "abhilesh7@gmail.com"

# Setup project paths consistent with your existing pipeline
ROOT = Path(__file__).resolve().parents[2]
REF_DIR = ROOT / "data" / "reference"
REF_DIR.mkdir(parents=True, exist_ok=True)

COORD_FILE = REF_DIR / "mtdna_gene_coordinates.tsv"


def fetch_mtdna_coords():
    """Retrieves canonical CDS boundaries from the rCRS (NC_012920.1) record."""
    print("Connecting to NCBI Entrez...")
    print("Fetching rCRS (NC_012920.1) GenBank record...")

    try:
        handle = Entrez.efetch(
            db="nucleotide", id="NC_012920.1", rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()
    except Exception as e:
        print(f"Error connecting to NCBI: {e}")
        return None

    gene_map = {}
    print(f"Parsing CDS features from {record.id}...")

    for feature in record.features:
        if feature.type == "CDS":
            # Extract official gene name (e.g., 'ND1')
            gene_name = feature.qualifiers.get("gene", [None])[0]
            if not gene_name:
                continue

            # Extract 1-indexed boundaries and strand orientation
            # BioPython uses 0-indexed half-open intervals; we convert to biological 1-indexing
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = "+" if feature.location.strand == 1 else "-"

            # Normalize to HGNC/MITOMAP standard (MT- prefix)
            standard_name = (
                f"MT-{gene_name}" if not gene_name.startswith("MT-") else gene_name
            )
            gene_map[standard_name] = (start, end, strand)

    return gene_map


def save_mtdna_coords(gene_map: dict):
    """Writes the coordinate dictionary to a tab-separated file for downstream use."""
    if not gene_map:
        print("No data to save.")
        return

    print(f"Saving coordinates to {COORD_FILE.name}...")
    with open(COORD_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        # Header matches the structure expected by your potential loaders
        writer.writerow(["gene", "start", "end", "strand"])

        for gene in sorted(gene_map.keys()):
            start, end, strand = gene_map[gene]
            writer.writerow([gene, start, end, strand])

    print("Mitochondrial coordinate map successfully generated.")


def main():
    # 1. Fetch live data from NCBI
    mt_coords = fetch_mtdna_coords()

    # 2. Persist to disk
    if mt_coords:
        save_mtdna_coords(mt_coords)

        # Quick verification of problematic genes like MT-ND6
        if "MT-ND6" in mt_coords:
            s, e, st = mt_coords["MT-ND6"]
            print(
                f"Verification: MT-ND6 correctly identified on strand ({st}) at {s}-{e}."
            )


if __name__ == "__main__":
    main()
