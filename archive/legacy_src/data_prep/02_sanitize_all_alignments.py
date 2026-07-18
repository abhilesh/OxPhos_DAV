import csv
from pathlib import Path

# ==== Configuration ====
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"

# Alignment Directories
TOGA_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_aa"
TOGA_NT_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"
MT_AA_DIR = DATA_DIR / "alignments" / "mtdna_aa"
MT_NT_DIR = DATA_DIR / "alignments" / "mtdna_codon"

# Reference Files
REF_DIR = DATA_DIR / "reference"
OVERVIEW_FILE = REF_DIR / "TOGA_overview_table_hg38.tsv"
MAP_CSV = REF_DIR / "taxid_species_mapping.csv"

# Universal Human ID
HUMAN_TAXID = "9606"


def load_toga_taxids() -> dict:
    """Loads Species -> TaxID mapping directly from the TOGA overview."""
    toga_to_taxid = {}
    if not OVERVIEW_FILE.exists():
        print(f"Warning: {OVERVIEW_FILE.name} not found.")
        return toga_to_taxid

    with open(OVERVIEW_FILE, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species = row.get("Species", "").strip().replace(" ", "_")
            taxid = row.get("Species Taxonomy ID", "").strip()
            if species and taxid:
                toga_to_taxid[species] = taxid
    return toga_to_taxid


def load_mt_mapping() -> dict:
    """Loads the generated NCBI mapping for mitochondrial accessions."""
    acc_map = {}
    if not MAP_CSV.exists():
        print(f"Warning: {MAP_CSV.name} not found.")
        return acc_map

    with open(MAP_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Handle variations in the CSV header name from previous steps
            acc = row.get("ncbi_accession", row.get("accession", "")).strip()
            taxid = row.get("taxid", "").strip()
            ncbi_sp = row.get("ncbi_species", "").strip()
            toga_sp = row.get("toga_species", "").strip()
            match_type = row.get("match_type", "").strip()

            # CRITICAL: Use TOGA species name if it's an exact match to ensure tree intersection
            final_species = (
                toga_sp if match_type == "Exact_TaxID_Match" and toga_sp else ncbi_sp
            )

            if acc:
                acc_map[acc] = {"taxid": taxid, "species": final_species}
    return acc_map


def sanitize_toga_directory(directory: Path, toga_to_taxid: dict):
    """Injects TaxIDs into TOGA nucDNA alignments."""
    if not directory.exists():
        return

    print(f"Sanitizing TOGA directory: {directory.name}...")
    for fasta in directory.glob("*.fasta"):
        lines = []
        with open(fasta, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # Current format: >Homo_sapiens|hg38|ENST... | GENE
                    # Or: >Species_Name|Assembly|ENST... | GENE
                    parts = line[1:].strip().split(" | ")
                    header_blocks = parts[0].split("|")

                    species = header_blocks[0]
                    assembly = header_blocks[1] if len(header_blocks) > 1 else ""
                    transcript = header_blocks[2] if len(header_blocks) > 2 else ""
                    gene_info = parts[1] if len(parts) > 1 else ""

                    if "Homo_sapiens" in species:
                        taxid = HUMAN_TAXID
                    else:
                        taxid = toga_to_taxid.get(species, "UNKNOWN_TAXID")

                    # New format: >Species|TaxID|Assembly|Transcript | GENE
                    new_header = f">{species}|{taxid}|{assembly}|{transcript}"
                    if gene_info:
                        new_header += f" | {gene_info}"

                    lines.append(new_header + "\n")
                else:
                    lines.append(line)

        with open(fasta, "w") as f:
            f.writelines(lines)


def sanitize_mt_directory(directory: Path, mt_map: dict):
    """Sanitizes NCBI mtDNA alignments to match TOGA formatting."""
    if not directory.exists():
        return

    print(f"Sanitizing MT directory: {directory.name}...")
    for fasta in directory.glob("*.fasta"):

        # Ensure filename has MT- prefix
        file_name = fasta.name
        base_gene = file_name.split("_")[0].replace("MT-", "")
        new_file_name = f"MT-{base_gene}_{'_'.join(file_name.split('_')[1:])}"
        new_file_path = directory / new_file_name

        lines = []
        with open(fasta, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # Parse: >NC_045205.1 | MT-CO1
                    parts = line[1:].strip().split("|")
                    acc_full = parts[0].strip()
                    acc_base = acc_full.split(".")[0]  # e.g. NC_045205
                    gene_info = (
                        parts[1].strip() if len(parts) > 1 else f"MT-{base_gene}"
                    )

                    # Explicit human catch
                    if "NC_012920" in acc_full:
                        species = "Homo_sapiens"
                        taxid = HUMAN_TAXID
                    else:
                        map_data = mt_map.get(acc_full, mt_map.get(acc_base, {}))
                        species = map_data.get("species", "UNKNOWN_SPECIES")
                        taxid = map_data.get("taxid", "UNKNOWN_TAXID")

                    # New format: >Species|TaxID|Accession | GENE
                    lines.append(f">{species}|{taxid}|{acc_full} | {gene_info}\n")
                else:
                    lines.append(line)

        with open(new_file_path, "w") as f:
            f.writelines(lines)

        if fasta != new_file_path:
            fasta.unlink()


def main():
    print("Loading taxonomic maps...")
    toga_to_taxid = load_toga_taxids()
    mt_map = load_mt_mapping()

    print("\nExecuting Header Sanitization...")
    sanitize_toga_directory(TOGA_AA_DIR, toga_to_taxid)
    sanitize_toga_directory(TOGA_NT_DIR, toga_to_taxid)

    sanitize_mt_directory(MT_AA_DIR, mt_map)
    sanitize_mt_directory(MT_NT_DIR, mt_map)

    print(
        "\nAlignment sanitization complete. All FASTA files now share a unified taxonomic nomenclature."
    )


if __name__ == "__main__":
    main()
