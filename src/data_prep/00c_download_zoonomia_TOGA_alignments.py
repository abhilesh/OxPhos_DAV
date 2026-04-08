"""
src/data_prep/00c_download_zoonomia_TOGA_alignments.py

Downloads codon (and, when available, amino-acid) TOGA alignments from the
Senckenberg server for all canonical nuclear OXPHOS subunits.

Transcript selection strategy
──────────────────────────────
The TOGA server publishes one file per Ensembl transcript (ENST), not per gene.
ClinVar annotates variants against RefSeq NM_ transcripts.  To avoid coordinate
mismatches we must ensure the ENST used in the alignment corresponds to the same
CDS as the NM_ used in our ClinVar records.

1.  Download the NCBI MANE Select summary table (one canonical NM_ ↔ ENST pair
    per gene; the pair is guaranteed to share the same CDS coordinates).
2.  Load `nucDNA_annotations.json` (if present) to find the NM_ actually used
    for each gene in our dataset.
3.  Build a gene → preferred ENST map from NM_ via MANE.
4.  When choosing which TOGA file to download for a gene, prefer the ENST that
    matches the ClinVar NM_.  Fall back to the first available ENST if no match.

Run from project root inside the Docker container:
    python src/data_prep/00c_download_zoonomia_TOGA_alignments.py
"""

import gzip
import io
import json
import re
import ssl
import urllib.request
import csv
from collections import defaultdict, Counter
from datetime import date
from pathlib import Path

from utils.parsers import GeneReference
from utils.utils import get_latest

# ── Paths ──────────────────────────────────────────────────────────────────────
today = date.today().isoformat()

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
CODON_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"
REF_DIR = DATA_DIR / "reference"
CURATED_DIR = DATA_DIR / "annotations" / "curated"

CODON_DIR.mkdir(parents=True, exist_ok=True)
REF_DIR.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://genome.senckenberg.de/download/TOGA/human_hg38_reference/MultipleCodonAlignments"
OVERVIEW_URL = "https://genome.senckenberg.de/download/TOGA/human_hg38_reference/overview.table.tsv"
MANE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz"

OVERVIEW_FILE = REF_DIR / f"TOGA_overview_table_hg38_{today}.tsv"
MANE_FILE = REF_DIR / "MANE_GRCh38_v1.5.txt.gz"
NUC_CURATED = CURATED_DIR / "nucDNA_annotations.json"

# Disable SSL verification for academic servers
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


# ── Network helpers ───────────────────────────────────────────────────────────

def _fetch(url: str) -> bytes:
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, context=ctx) as r:
        return r.read()


# ── Overview table (assembly → species) ──────────────────────────────────────

def fetch_assembly_map() -> dict:
    """Returns {assembly_id: species_name}."""
    if not OVERVIEW_FILE.exists():
        print("Downloading TOGA overview table...")
        try:
            OVERVIEW_FILE.write_bytes(_fetch(OVERVIEW_URL))
        except Exception as e:
            print(f"  Error: {e}")
            return {}

    mapping = {}
    with open(OVERVIEW_FILE, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            assembly = row.get("Assembly name", "").strip()
            species = row.get("Species", "").strip().replace(" ", "_")
            if assembly and species:
                mapping[assembly] = species
    print(f"Loaded {len(mapping)} assembly→species mappings.")
    return mapping


# ── MANE Select table (NM_ → ENST) ───────────────────────────────────────────

def fetch_mane_mapping() -> dict:
    """
    Downloads the NCBI MANE Select summary and returns {NM_base: ENST_base}.

    MANE Select guarantees identical CDS coordinates between the RefSeq NM_
    and the Ensembl ENST transcript.  Versioned IDs are stripped (NM_005006.7
    → NM_005006) so they match what comes out of ClinVar.
    """
    if not MANE_FILE.exists():
        print("Downloading NCBI MANE Select summary...")
        try:
            MANE_FILE.write_bytes(_fetch(MANE_URL))
        except Exception as e:
            print(f"  Error downloading MANE: {e}")
            return {}

    mapping = {}
    with gzip.open(MANE_FILE, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            # Columns (v1.4): NCBI_GeneID, Ensembl_Gene, HGNC_ID, status, name,
            #                  RefSeq_nuc, RefSeq_prot, Ensembl_nuc, Ensembl_prot,
            #                  MANE_status
            if len(cols) < 8:
                continue
            nm_versioned = cols[5].strip()   # e.g. NM_005006.7
            enst_versioned = cols[7].strip() # e.g. ENST00000378781.9
            if nm_versioned.startswith("NM_") and enst_versioned.startswith("ENST"):
                nm_base   = nm_versioned.split(".")[0]
                enst_base = enst_versioned.split(".")[0]
                mapping[nm_base] = enst_base

    print(f"Loaded {len(mapping)} NM_ → ENST MANE Select pairs.")
    return mapping


# ── ClinVar NM_ per gene ──────────────────────────────────────────────────────

def get_clinvar_nm_per_gene() -> dict:
    """
    Reads the curated nucDNA annotations and returns the most common NM_
    transcript (unversioned) used for each gene.

    Returns {gene_symbol: NM_base}.
    """
    if not NUC_CURATED.exists():
        print("  nucDNA_annotations.json not found; skipping NM_ pre-selection.")
        return {}

    nm_counts: dict[str, Counter] = defaultdict(Counter)
    with open(NUC_CURATED) as f:
        variants = json.load(f)
    for v in variants:
        gene = v.get("locus", "")
        tx   = v.get("transcript_id", "").split(".")[0]
        if gene and tx.startswith("NM_"):
            nm_counts[gene][tx] += 1

    # Pick the most-used NM_ per gene (nearly always unanimous)
    result = {gene: counts.most_common(1)[0][0]
              for gene, counts in nm_counts.items()}
    print(f"Found NM_ transcripts for {len(result)} genes in ClinVar data.")
    return result


# ── TOGA directory index ──────────────────────────────────────────────────────

def fetch_toga_index() -> dict:
    """
    Scrapes the Senckenberg index and returns {gene_symbol: [filename, ...]}.

    Collects ALL transcript filenames per gene (there may be several ENST per gene).
    """
    print("Fetching TOGA directory index...")
    try:
        html = _fetch(f"{BASE_URL}/").decode("utf-8")
    except Exception as e:
        print(f"  Error: {e}")
        return {}

    # filename format: ENST00000378781.NDUFS1.fasta.gz
    pattern = re.compile(r'href="(ENST\d+\.([^.]+)\.fasta\.gz)"')
    index: dict[str, list] = defaultdict(list)
    for filename, gene_symbol in pattern.findall(html):
        index[gene_symbol].append(filename)

    total_files = sum(len(v) for v in index.values())
    print(f"Found {total_files} alignment files for {len(index)} gene symbols.")
    return dict(index)


# ── Alignment file writer ─────────────────────────────────────────────────────

def download_and_write(
    gene: str,
    remote_filename: str,
    assembly_map: dict,
    out_file: Path,
) -> bool:
    """
    Downloads one TOGA codon alignment, remaps headers, and writes the FASTA.

    Human REFERENCE header → Homo_sapiens|9606|hg38|{ENST} | {gene}
    Species headers         → {species}|{taxon_id}|{assembly}|{ENST}

    Returns True on success.
    """
    transcript_id = remote_filename.split(".")[0]   # ENST part
    url = f"{BASE_URL}/{remote_filename}"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, context=ctx) as response:
            with gzip.GzipFile(fileobj=response) as uncompressed:
                stream = io.TextIOWrapper(uncompressed, encoding="utf-8")
                with open(out_file, "w", encoding="utf-8") as out:
                    for line in stream:
                        if not line.startswith(">"):
                            out.write(line)
                            continue

                        if line.startswith(">REFERENCE"):
                            out.write(
                                f">Homo_sapiens|9606|hg38|{transcript_id} | {gene}\n"
                            )
                        else:
                            parts = line.strip().split()
                            header_base = parts[0]
                            if header_base.startswith(">vs_"):
                                assembly_id = header_base[4:]
                                species = assembly_map.get(assembly_id, assembly_id)
                                new_header = f">{species}|{assembly_id}|{transcript_id}"
                                rest = " ".join(parts[1:])
                                if rest:
                                    new_header += f" | {rest}"
                                out.write(new_header + "\n")
                            else:
                                out.write(line)
        return True
    except Exception as e:
        print(f"  [FAIL] {gene}: {e}")
        return False


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("Initializing Targeted TOGA Download with MANE Transcript Selection\n")

    assembly_map  = fetch_assembly_map()
    mane_mapping  = fetch_mane_mapping()       # NM_base → ENST_base
    clinvar_nm    = get_clinvar_nm_per_gene()  # gene → NM_base
    toga_index    = fetch_toga_index()         # gene_symbol → [filename, ...]

    if not assembly_map or not toga_index:
        print("Cannot proceed without assembly map and TOGA index.")
        return

    # Build gene → preferred ENST from MANE (via ClinVar NM_)
    preferred_enst: dict[str, str] = {}
    for gene, nm_base in clinvar_nm.items():
        enst = mane_mapping.get(nm_base)
        if enst:
            preferred_enst[gene] = enst

    try:
        hgnc_file = get_latest(DATA_DIR, "Canonical_OXPHOS_Subunits_HGNC*.csv")
        hgnc_ref  = GeneReference(hgnc_file)
    except Exception as e:
        print(f"Error loading HGNC reference: {e}")
        return

    target_genes = sorted({
        data["symbol"]
        for symbol, data in hgnc_ref.lookup.items()
        if data.get("symbol") and not data["symbol"].startswith("MT-")
    })

    print(f"\nProcessing {len(target_genes)} genes...\n")

    success, skipped, missing, wrong_tx = 0, 0, [], []

    for gene in target_genes:
        out_file = CODON_DIR / f"{gene}_codon_alignment.fasta"

        # ── Resolve candidate filenames for this gene ────────────────────────
        candidates = toga_index.get(gene, [])

        # Also check aliases
        if not candidates:
            gene_data = hgnc_ref.get_gene_data(gene)
            if gene_data:
                for alias in hgnc_ref.lookup:
                    if hgnc_ref.lookup[alias].get("symbol") == gene:
                        if alias in toga_index:
                            candidates = toga_index[alias]
                            print(f"  [ALIAS] {gene} → TOGA symbol '{alias}'")
                            break

        if not candidates:
            print(f"  [MISS]  {gene}: no alignment on TOGA server")
            missing.append(gene)
            continue

        # ── Prefer the MANE/ClinVar-matched ENST ────────────────────────────
        pref = preferred_enst.get(gene)
        chosen = None
        if pref:
            for fn in candidates:
                if fn.startswith(pref + "."):
                    chosen = fn
                    break

        if chosen is None:
            chosen = candidates[0]
            if pref and len(candidates) > 1:
                wrong_tx.append(
                    f"{gene}: wanted {pref}, "
                    f"available [{', '.join(f.split('.')[0] for f in candidates)}], "
                    f"using {chosen.split('.')[0]}"
                )

        # ── Skip if already downloaded with the correct transcript ───────────
        current_enst = None
        if out_file.exists() and out_file.stat().st_size > 100:
            m = re.search(r"(ENST\d+)", out_file.read_text(encoding="utf-8")[:200])
            current_enst = m.group(1) if m else None

        wanted_enst = chosen.split(".")[0]
        if current_enst == wanted_enst:
            skipped += 1
            continue
        elif current_enst and current_enst != wanted_enst:
            print(f"  [RETX]  {gene}: replacing {current_enst} → {wanted_enst}")

        # ── Download ─────────────────────────────────────────────────────────
        ok = download_and_write(gene, chosen, assembly_map, out_file)
        if ok:
            nm = clinvar_nm.get(gene, "?")
            tx_note = f" (ClinVar: {nm} → MANE: {wanted_enst})" if pref else f" (ENST: {wanted_enst})"
            print(f"  [DOWN]  {gene}{tx_note}")
            success += 1
        else:
            missing.append(gene)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*55}")
    print("TOGA DOWNLOAD SUMMARY")
    print(f"{'='*55}")
    print(f"  Downloaded  : {success}")
    print(f"  Already up-to-date : {skipped}")
    print(f"  Not on server      : {len(missing)}")
    if missing:
        print(f"    {', '.join(missing)}")

    if wrong_tx:
        print(f"\n  Genes where preferred ENST was not available ({len(wrong_tx)}):")
        for msg in wrong_tx:
            print(f"    {msg}")


if __name__ == "__main__":
    main()
