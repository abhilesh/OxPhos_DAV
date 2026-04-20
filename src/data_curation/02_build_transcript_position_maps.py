#!/usr/bin/env python3
"""
Build deterministic NM-to-ENST position maps for nuclear OXPHOS genes.

These maps resolve transcript mismatches between ClinVar NM accessions and TOGA
ENST alignments.
"""

from pathlib import Path
import csv
import gzip
import json
import re
import socket
import sys
import time
from collections import Counter

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

from utils.gene_reference import GeneReference

DATA_DIR = ROOT / "data"
RAW_REFERENCE_DIR = DATA_DIR / "raw" / "reference"
TOGA_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_aa"
DERIVED_CURATED_DIR = DATA_DIR / "derived" / "curated"
LEGACY_REFERENCE_DIR = DATA_DIR / "reference"
TRANSCRIPT_MAPS_PARQUET = DERIVED_CURATED_DIR / "transcript_position_maps.parquet"
TRANSCRIPT_MAPS_JSON = DERIVED_CURATED_DIR / "transcript_position_maps.json"
COMPAT_TRANSCRIPT_MAPS_JSON = LEGACY_REFERENCE_DIR / "transcript_position_maps.json"

Entrez.email = "pipeline@analysis.local"
socket.setdefaulttimeout(30)

TOGA_TRANSCRIPT_EXCEPTIONS = {
    "ATP5MC2": "ENST00000673498",
    "ATP5MF": "ENST00000449683",
    "ATP5PF": "ENST00000400099",
    "COX5A": "ENST00000568783",
    "COXFA4L2": "ENST00000556732",
    "NDUFA10": "ENST00000307300",
    "NDUFA11": "ENST00000418389",
    "NDUFA13": "ENST00000428459",
    "NDUFB1": "ENST00000617122",
    "NDUFS6": "ENST00000469176",
    "NDUFS7": "ENST00000414651",
    "NDUFV2": "ENST00000400033",
    "UQCRB": "ENST00000523920",
}

ALIGNER = PairwiseAligner()
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.mode = "global"
ALIGNER.open_gap_score = -10
ALIGNER.extend_gap_score = -0.5


def ensure_layout() -> None:
    for path in (DERIVED_CURATED_DIR, LEGACY_REFERENCE_DIR):
        path.mkdir(parents=True, exist_ok=True)


def latest(pattern: str) -> Path:
    matches = sorted(RAW_REFERENCE_DIR.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)
    if not matches:
        raise FileNotFoundError(f"Missing {pattern}")
    return matches[0]


def mane_nm_map(mane_path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(mane_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            cols = line.lstrip("#").rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            symbol = cols[3].strip()
            nm = cols[5].strip()
            if symbol and nm.startswith("NM_"):
                mapping[symbol] = nm.split(".")[0]
    return mapping


def dominant_clinvar_nm(clinvar_counts: dict[str, int]) -> str:
    if not clinvar_counts:
        return ""
    return max(
        clinvar_counts.items(),
        key=lambda item: (item[1], item[0]),
    )[0]


def extract_nm_from_name(name: str) -> str:
    match = re.search(r"(NM_[\d.]+)", str(name))
    return match.group(1).split(".")[0] if match else ""


def iter_clinvar_rows(clinvar_path: Path):
    import pandas as pd

    for chunk in pd.read_csv(clinvar_path, sep="\t", compression="gzip", chunksize=50000, low_memory=False):
        for row in chunk.to_dict(orient="records"):
            yield row


def clinvar_nm_counts(clinvar_path: Path) -> dict[str, dict[str, int]]:
    counts: dict[str, Counter] = {}
    for row in iter_clinvar_rows(clinvar_path):
        gene = str(row.get("GeneSymbol", "")).strip()
        tx = extract_nm_from_name(row.get("Name", ""))
        if gene and tx:
            counts.setdefault(gene, Counter())[tx] += 1
    return {gene: dict(counter) for gene, counter in counts.items()}


def enst_protein(gene: str) -> tuple[str, str]:
    fasta = TOGA_AA_DIR / f"{gene}_aa_alignment.fasta"
    if not fasta.exists():
        return "", ""
    with open(fasta, encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            header = f">{record.description}"
            if header.startswith(">Homo_sapiens|9606|") or header.startswith(">Homo_sapiens|9606"):
                match = re.search(r"(ENST\d+)", header)
                enst = TOGA_TRANSCRIPT_EXCEPTIONS.get(gene, match.group(1) if match else "")
                return enst, str(record.seq).replace("-", "").replace("*", "")
    return "", ""


def nm_protein(nm_id: str) -> str:
    try:
        handle = Entrez.efetch(db="nucleotide", id=nm_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
    except Exception:
        return ""
    for feat in record.features:
        if feat.type == "CDS":
            return feat.qualifiers.get("translation", [""])[0]
    return ""


def build_pos_map(nm_seq: str, enst_seq: str) -> dict[int, int | None]:
    best = next(iter(ALIGNER.align(nm_seq, enst_seq)))
    aligned_nm = str(best[0])
    aligned_enst = str(best[1])
    pos_map: dict[int, int | None] = {}
    nm_pos = enst_pos = 0
    for nm_char, enst_char in zip(aligned_nm, aligned_enst):
        if nm_char != "-":
            nm_pos += 1
        if enst_char != "-":
            enst_pos += 1
        if nm_char != "-":
            pos_map[nm_pos] = None if enst_char == "-" else enst_pos
    return pos_map


def main() -> None:
    import pandas as pd

    ensure_layout()
    hgnc_file = latest("Canonical_OXPHOS_Subunits_HGNC_*.csv")
    mane_file = latest("MANE_GRCh38_v1.5_*.txt.gz")
    clinvar_file = DATA_DIR / "raw" / "annotations"
    clinvar_file = sorted(clinvar_file.glob("ClinVar_VariantSummary_*.txt.gz"), key=lambda p: p.stat().st_mtime, reverse=True)[0]

    hgnc_ref = GeneReference(hgnc_file)
    mane_map = mane_nm_map(mane_file)
    clinvar_nm = clinvar_nm_counts(clinvar_file)

    rows: list[dict] = []
    for gene in sorted({v["symbol"] for v in hgnc_ref.lookup.values() if not v["symbol"].startswith("MT-")}):
        print(f"Processing transcript map for {gene}...")
        gene_clinvar_counts = clinvar_nm.get(gene, {})
        mane_nm = mane_map.get(gene, "")
        clinvar_nm_choice = dominant_clinvar_nm(gene_clinvar_counts)
        preferred_nm = mane_nm or clinvar_nm_choice
        selection_rule = "mane_select" if mane_nm else "clinvar_dominant_nm" if clinvar_nm_choice else "no_preferred_nm"
        enst_id, enst_seq = enst_protein(gene)
        if not preferred_nm or not enst_seq:
            rows.append({
                "gene": gene,
                "preferred_nm": preferred_nm,
                "mane_nm": mane_nm,
                "clinvar_dominant_nm": clinvar_nm_choice,
                "toga_enst": enst_id,
                "map_status": "missing_nm_or_toga",
                "identity_fraction": None,
                "coverage_fraction": None,
                "mapped_positions": 0,
                "nm_protein_length": None,
                "enst_protein_length": len(enst_seq) if enst_seq else None,
                "map_json": "{}",
                "selection_rule": selection_rule,
                "clinvar_nm_counts_json": json.dumps(gene_clinvar_counts, sort_keys=True),
            })
            continue

        nm_seq = nm_protein(preferred_nm)
        time.sleep(0.35)
        if not nm_seq:
            rows.append({
                "gene": gene,
                "preferred_nm": preferred_nm,
                "mane_nm": mane_nm,
                "clinvar_dominant_nm": clinvar_nm_choice,
                "toga_enst": enst_id,
                "map_status": "nm_fetch_failed",
                "identity_fraction": None,
                "coverage_fraction": None,
                "mapped_positions": 0,
                "nm_protein_length": None,
                "enst_protein_length": len(enst_seq) if enst_seq else None,
                "map_json": "{}",
                "selection_rule": selection_rule,
                "clinvar_nm_counts_json": json.dumps(gene_clinvar_counts, sort_keys=True),
            })
            continue

        if nm_seq == enst_seq:
            pos_map: dict[int, int | None] = {i: i for i in range(1, len(nm_seq) + 1)}
            status = "identity"
            identity_fraction = 1.0
            coverage_fraction = 1.0
            mapped_positions = len(nm_seq)
        else:
            pos_map = build_pos_map(nm_seq, enst_seq)
            mapped = sum(
                1 for nm_pos, enst_pos in pos_map.items()
                if enst_pos is not None and nm_seq[nm_pos - 1] == enst_seq[enst_pos - 1]
            )
            identity_fraction = mapped / len(nm_seq) if nm_seq else 0.0
            mapped_positions = sum(1 for enst_pos in pos_map.values() if enst_pos is not None)
            coverage_fraction = mapped_positions / len(nm_seq) if nm_seq else 0.0
            status = "mapped"

        rows.append({
            "gene": gene,
            "preferred_nm": preferred_nm,
            "mane_nm": mane_nm,
            "clinvar_dominant_nm": clinvar_nm_choice,
            "toga_enst": enst_id,
            "map_status": status,
            "identity_fraction": identity_fraction,
            "coverage_fraction": coverage_fraction,
            "mapped_positions": mapped_positions,
            "nm_protein_length": len(nm_seq),
            "enst_protein_length": len(enst_seq),
            "map_json": json.dumps({str(k): v for k, v in pos_map.items()}, sort_keys=True),
            "selection_rule": selection_rule,
            "clinvar_nm_counts_json": json.dumps(gene_clinvar_counts, sort_keys=True),
        })

    pd.DataFrame(rows).to_parquet(TRANSCRIPT_MAPS_PARQUET, index=False)

    compat = {}
    for row in rows:
        compat[row["gene"]] = {
            "nm": row["preferred_nm"],
            "mane_nm": row["mane_nm"],
            "clinvar_dominant_nm": row["clinvar_dominant_nm"],
            "enst": row["toga_enst"],
            "type": row["map_status"],
            "identity_fraction": row["identity_fraction"],
            "coverage_fraction": row["coverage_fraction"],
            "mapped_positions": row["mapped_positions"],
            "nm_protein_length": row["nm_protein_length"],
            "enst_protein_length": row["enst_protein_length"],
            "selection_rule": row["selection_rule"],
            "clinvar_nm_counts": json.loads(row["clinvar_nm_counts_json"]),
            "map": json.loads(row["map_json"]),
        }
    with open(TRANSCRIPT_MAPS_JSON, "w", encoding="utf-8") as handle:
        json.dump(compat, handle, indent=2)
    with open(COMPAT_TRANSCRIPT_MAPS_JSON, "w", encoding="utf-8") as handle:
        json.dump(compat, handle, indent=2)
    print(f"Wrote transcript maps: {TRANSCRIPT_MAPS_PARQUET}")


if __name__ == "__main__":
    main()
