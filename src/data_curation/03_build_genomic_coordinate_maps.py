#!/usr/bin/env python3
"""
Build genomic coordinate rescue maps for target nuclear genes.

These maps act as the authoritative rescue layer for defined TOGA and MANE
transcript discordance cases.
"""

from pathlib import Path
import json
import sys
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

DATA_DIR = ROOT / "data"
DERIVED_CURATED_DIR = DATA_DIR / "derived" / "curated"
LEGACY_REFERENCE_DIR = DATA_DIR / "reference"
GENOMIC_MAPS_PARQUET = DERIVED_CURATED_DIR / "genomic_coordinate_maps.parquet"
GENOMIC_MAPS_JSON = DERIVED_CURATED_DIR / "genomic_coordinate_maps.json"
COMPAT_GENOMIC_MAPS_JSON = LEGACY_REFERENCE_DIR / "genomic_coordinate_maps.json"

ENSEMBL_REST = "https://rest.ensembl.org"
TOGA_ENST_GENES = {
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


def ensure_layout() -> None:
    for path in (DERIVED_CURATED_DIR, LEGACY_REFERENCE_DIR):
        path.mkdir(parents=True, exist_ok=True)


def fetch_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.loads(response.read())


def build_genomic_cds_map(enst_id: str) -> tuple[dict, dict] | tuple[None, None]:
    url = f"{ENSEMBL_REST}/lookup/id/{enst_id}?expand=1&content-type=application/json"
    try:
        data = fetch_json(url)
    except Exception:
        return None, None

    strand = data.get("strand")
    exons = data.get("Exon", [])
    translation = data.get("Translation")
    if strand is None or not exons or not translation:
        return None, None

    cds_low = translation["start"]
    cds_high = translation["end"]
    prot_len = translation["length"]
    exons_sorted = sorted(exons, key=lambda item: item["start"])

    positions: dict[int, dict] = {}
    cds_pos = 0
    if strand == 1:
        exon_iter = exons_sorted
        pos_fn = lambda exon: range(exon["start"], exon["end"] + 1)
    else:
        exon_iter = reversed(exons_sorted)
        pos_fn = lambda exon: range(exon["end"], exon["start"] - 1, -1)

    for exon in exon_iter:
        for genomic_pos in pos_fn(exon):
            if cds_low <= genomic_pos <= cds_high:
                cds_pos += 1
                positions[genomic_pos] = {"cds_pos": cds_pos, "aa_pos": (cds_pos - 1) // 3 + 1}

    meta = {
        "enst": enst_id,
        "strand": strand,
        "cds_genomic_start": cds_high if strand == -1 else cds_low,
        "cds_genomic_end": cds_low if strand == -1 else cds_high,
        "prot_len": prot_len,
        "n_cds_positions": len(positions),
        "map_status": "built",
    }
    return meta, positions


def main() -> None:
    import pandas as pd

    ensure_layout()
    rows: list[dict] = []
    compat: dict = {}
    for gene, enst_id in sorted(TOGA_ENST_GENES.items()):
        meta, positions = build_genomic_cds_map(enst_id)
        time.sleep(0.2)
        if meta is None:
            rows.append({"gene": gene, "enst": enst_id, "map_status": "fetch_failed", "map_json": "{}"})
            compat[gene] = {"enst": enst_id, "map_status": "fetch_failed", "map": {}}
            continue
        rows.append({**meta, "gene": gene, "map_json": json.dumps({str(k): v for k, v in positions.items()}, sort_keys=True)})
        compat[gene] = {**meta, "map": {str(k): v for k, v in positions.items()}}

    pd.DataFrame(rows).to_parquet(GENOMIC_MAPS_PARQUET, index=False)
    with open(GENOMIC_MAPS_JSON, "w", encoding="utf-8") as handle:
        json.dump(compat, handle, indent=2)
    with open(COMPAT_GENOMIC_MAPS_JSON, "w", encoding="utf-8") as handle:
        json.dump(compat, handle, indent=2)
    print(f"Wrote genomic maps: {GENOMIC_MAPS_PARQUET}")


if __name__ == "__main__":
    main()
