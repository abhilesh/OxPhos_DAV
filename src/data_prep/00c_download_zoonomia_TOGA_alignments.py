import ssl
import gzip
import urllib.request
import re
import io
import csv
from pathlib import Path

# Import utility classes
from utils.parsers import GeneReference
from utils.utils import get_latest

# ==== Configuration ====
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
ALIGN_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"
REF_DIR = DATA_DIR / "reference"

ALIGN_DIR.mkdir(parents=True, exist_ok=True)
REF_DIR.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://genome.senckenberg.de/download/TOGA/human_hg38_reference/MultipleCodonAlignments"
OVERVIEW_URL = "https://genome.senckenberg.de/download/TOGA/human_hg38_reference/overview.table.tsv"
OVERVIEW_FILE = REF_DIR / "TOGA_overview_table_hg38.tsv"

# Disable SSL verification for academic servers
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


def fetch_and_parse_overview() -> dict:
    """Downloads the TOGA overview table and creates an Assembly -> Species map."""
    if not OVERVIEW_FILE.exists():
        print("Downloading TOGA overview table...")
        try:
            req = urllib.request.Request(OVERVIEW_URL)
            with urllib.request.urlopen(req, context=ctx) as response:
                with open(OVERVIEW_FILE, "wb") as f:
                    f.write(response.read())
        except Exception as e:
            print(f"Error downloading overview table: {e}")
            return {}

    assembly_map = {}
    try:
        with open(OVERVIEW_FILE, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                assembly = row.get("Assembly name", "").strip()
                # Replace spaces with underscores for clean FASTA formatting
                species = row.get("Species", "").strip().replace(" ", "_")
                if assembly and species:
                    assembly_map[assembly] = species
    except Exception as e:
        print(f"Error reading overview table: {e}")

    print(
        f"Loaded {len(assembly_map)} assembly-to-species mappings from the overview table."
    )
    return assembly_map


def fetch_toga_index() -> dict:
    """Scrapes the server index to map Gene Symbols to exact TOGA filenames."""
    print("Fetching TOGA directory index from Senckenberg server...")
    try:
        req = urllib.request.Request(f"{BASE_URL}/")
        with urllib.request.urlopen(req, context=ctx) as response:
            html = response.read().decode("utf-8")

        pattern = re.compile(r'href="(ENST\d+\.([^.]+)\.fasta\.gz)"')
        toga_index = {match[1]: match[0] for match in pattern.findall(html)}
        return toga_index
    except Exception as e:
        print(f"Error fetching directory index: {e}")
        return {}


def main():
    print("Initializing Targeted TOGA Download with Species Mapping...\n")

    toga_map = fetch_and_parse_overview()
    if not toga_map:
        return

    try:
        hgnc_file = get_latest(DATA_DIR, "Canonical_OXPHOS_Subunits_HGNC*.csv")
        hgnc_ref = GeneReference(hgnc_file)
    except Exception as e:
        print(f"Error loading HGNC reference: {e}")
        return

    # Extract strictly primary nuclear genes
    target_genes = sorted(
        list(
            {
                data.get("primary_symbol")
                for symbol, data in hgnc_ref.lookup.items()
                if data.get("primary_symbol")
                and not data.get("primary_symbol").startswith("MT-")
            }
        )
    )

    toga_index = fetch_toga_index()
    if not toga_index:
        return

    success_count = 0
    missing_genes = []

    print(f"\nAttempting to download and map {len(target_genes)} alignments...\n")

    for gene in target_genes:
        remote_filename = None

        if gene in toga_index:
            remote_filename = toga_index[gene]
        else:
            gene_data = hgnc_ref.get_gene_data(gene)
            if gene_data:
                primary = gene_data["primary_symbol"]
                aliases = [
                    k
                    for k, v in hgnc_ref.lookup.items()
                    if v["primary_symbol"] == primary
                ]
                for alias in aliases:
                    if alias in toga_index:
                        remote_filename = toga_index[alias]
                        print(
                            f"  [ALIAS] Mapped local '{gene}' to TOGA remote '{alias}'"
                        )
                        break

        if not remote_filename:
            print(f"  [MISS] {gene} : Not found under any known alias.")
            missing_genes.append(gene)
            continue

        # 1. Extract the ENST transcript ID from the filename
        transcript_id = remote_filename.split(".")[0]

        url = f"{BASE_URL}/{remote_filename}"
        out_file = ALIGN_DIR / f"{gene}_codon_alignment.fasta"

        if out_file.exists() and out_file.stat().st_size > 100:
            print(f"  [SKIP] {gene} (Already downloaded & mapped)")
            success_count += 1
            continue

        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, context=ctx) as response:
                with gzip.GzipFile(fileobj=response) as uncompressed:
                    text_stream = io.TextIOWrapper(uncompressed, encoding="utf-8")
                    with open(out_file, "w", encoding="utf-8") as f_out:
                        for line in text_stream:
                            if line.startswith(">"):
                                if line.startswith(">REFERENCE"):
                                    # 2. Inject transcript_id into the human reference header
                                    f_out.write(
                                        f">Homo_sapiens|hg38|{transcript_id} | {gene}\n"
                                    )
                                else:
                                    parts = line.strip().split()
                                    header_base = parts[0]

                                    if header_base.startswith(">vs_"):
                                        assembly_id = header_base[4:]
                                        species = toga_map.get(assembly_id, assembly_id)

                                        # 3. Inject transcript_id into all species headers
                                        new_header = (
                                            f">{species}|{assembly_id}|{transcript_id}"
                                        )

                                        rest_of_header = " ".join(parts[1:])
                                        if rest_of_header:
                                            new_header += f" | {rest_of_header}"

                                        f_out.write(new_header + "\n")
                                    else:
                                        f_out.write(line)
                            else:
                                f_out.write(line)

            print(f"  [DOWN] {gene} <- Mapped and Saved")
            success_count += 1
        except Exception as e:
            print(f"  [FAIL] {gene} : {e}")

    print(f"\n{'='*50}")
    print("TOGA DOWNLOAD SUMMARY")
    print(f"{'='*50}")
    print(
        f"Successfully downloaded and mapped: {success_count} / {len(target_genes)} alignments."
    )


if __name__ == "__main__":
    main()
