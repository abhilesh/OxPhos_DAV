"""
src/data_prep/00g_build_genomic_coordinate_maps.py

Builds GRCh38 genomic position → TOGA ENST CDS/AA position maps for the 13
nuclear OXPHOS genes where the TOGA ENST differs from the MANE Select transcript.

Why this is needed
──────────────────
TOGA (Zoonomia project) was built with GENCODE v34 (Ensembl 100, 2020). MANE Select
was standardised later. For 13 genes, TOGA chose a longer isoform that predates MANE.
ClinVar records annotated against *older* NM_ versions often used the same longer
isoform that TOGA has — not the current MANE-matched NM_.

The NM_→ENST global position maps built by 00f fail for these genes because they
compare the *current* NM_ (shorter, MANE-matched) with the *TOGA* ENST (longer),
producing a systematic offset that sends every variant to the wrong alignment column.

Strategy B: genomic coordinate anchor
──────────────────────────────────────
GRCh38 genomic positions are stable across NM_ versions. For each TOGA ENST:
  1. Fetch exon structure + CDS boundaries from Ensembl REST.
  2. Walk exons in transcript order, assigning CDS position to every exonic base
     that lies within the CDS (between Translation.start and Translation.end).
  3. Store {grch38_pos → {cds_pos, aa_pos}} for each coding base.

At classification time, look up var["genomic_pos"] directly in this map to obtain
the correct TOGA CDS and AA positions, bypassing the NM_ c. coordinate entirely.

Watchouts implemented
─────────────────────
• Off-by-one: Translation.start / end are 1-based inclusive genomic coords from the
  Ensembl REST API. Verified against ref allele in the classification script.
• Strand: minus-strand genes traverse each exon from end→start, and exons are
  visited in reverse genomic order (5'→3' of the transcript).
• Intronic / UTR variants: positions outside the CDS window or in introns are simply
  absent from the map → classify script logs GENOMIC_POS_NOT_IN_ENST.
• Ensembl version: current REST API (v112) used; ENST IDs are stable. The ref-allele
  check in the classify script catches any coordinate drift from assembly patches.

Run from project root inside the Docker container:
    python src/data_prep/00g_build_genomic_coordinate_maps.py
"""

import json
import time
import urllib.request
from pathlib import Path

ROOT     = Path(__file__).resolve().parents[2]
OUT_FILE = ROOT / "data" / "reference" / "genomic_coordinate_maps.json"

ENSEMBL_REST = "https://rest.ensembl.org"

# The 13 TOGA ENSTs that differ from MANE Select and cause position mismatches.
# Determined by comparing transcript_position_maps.json identity_fraction < 100%
# OR by 100% identity but with a systematic N-terminal offset (NDUFB1, ATP5MC2).
TOGA_ENST_GENES = {
    "ATP5MC2":  "ENST00000673498",
    "ATP5MF":   "ENST00000449683",
    "ATP5PF":   "ENST00000400099",
    "COX5A":    "ENST00000568783",
    "COXFA4L2": "ENST00000556732",
    "NDUFA10":  "ENST00000307300",
    "NDUFA11":  "ENST00000418389",
    "NDUFA13":  "ENST00000428459",
    "NDUFB1":   "ENST00000617122",
    "NDUFS6":   "ENST00000469176",
    "NDUFS7":   "ENST00000414651",
    "NDUFV2":   "ENST00000400033",
    "UQCRB":    "ENST00000523920",
}


def _fetch_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())


def build_genomic_cds_map(enst_id: str) -> tuple[dict, dict] | tuple[None, None]:
    """
    Fetches transcript data from Ensembl REST and returns:
        (meta, positions)

    meta      : {enst, strand, cds_genomic_start (ATG), cds_genomic_end (stop), prot_len}
    positions : {grch38_pos (int): {"cds_pos": int, "aa_pos": int}}

    Returns (None, None) on failure.
    """
    url = (f"{ENSEMBL_REST}/lookup/id/{enst_id}"
           "?expand=1&content-type=application/json")
    try:
        data = _fetch_json(url)
    except Exception as exc:
        print(f"FETCH ERROR: {exc}")
        return None, None

    strand = data.get("strand")
    exons  = data.get("Exon", [])
    transl = data.get("Translation")

    if strand is None or not exons or not transl:
        print(f"INCOMPLETE REST RESPONSE (strand={strand}, "
              f"exons={len(exons)}, transl={bool(transl)})")
        return None, None

    # Ensembl REST convention (both strands):
    #   Translation.start = lower  genomic coordinate of the CDS
    #   Translation.end   = higher genomic coordinate of the CDS
    #   → For + strand: start = ATG,  end = last base of stop codon
    #   → For − strand: end   = ATG,  start = last base of stop codon
    cds_low  = transl["start"]   # always the smaller coordinate
    cds_high = transl["end"]     # always the larger coordinate
    prot_len = transl["length"]  # amino acids (stop codon excluded)

    # Sort exons by ascending genomic coordinate
    exons_sorted = sorted(exons, key=lambda e: e["start"])

    positions: dict[int, dict] = {}
    cds_pos = 0

    if strand == 1:
        # Plus strand: walk exons low→high (5'→3')
        for exon in exons_sorted:
            for gpos in range(exon["start"], exon["end"] + 1):
                if cds_low <= gpos <= cds_high:
                    cds_pos += 1
                    aa_pos = (cds_pos - 1) // 3 + 1
                    positions[gpos] = {"cds_pos": cds_pos, "aa_pos": aa_pos}
    else:
        # Minus strand: walk exons high→low (5'→3')
        for exon in reversed(exons_sorted):
            for gpos in range(exon["end"], exon["start"] - 1, -1):
                if cds_low <= gpos <= cds_high:
                    cds_pos += 1
                    aa_pos = (cds_pos - 1) // 3 + 1
                    positions[gpos] = {"cds_pos": cds_pos, "aa_pos": aa_pos}

    # Sanity check: CDS length = prot_len * 3 codons + 1 stop codon * 3 bases
    expected_cds = prot_len * 3 + 3
    if abs(len(positions) - expected_cds) > 3:
        print(f"\n  WARNING: CDS length mismatch — "
              f"got {len(positions)} positions, expected {expected_cds} "
              f"({prot_len} aa). Check exon/translation boundaries.")

    meta = {
        "enst":              enst_id,
        "strand":            strand,
        "cds_genomic_start": cds_high if strand == -1 else cds_low,   # ATG
        "cds_genomic_end":   cds_low  if strand == -1 else cds_high,  # stop end
        "prot_len":          prot_len,
        "n_cds_positions":   len(positions),
    }
    return meta, positions


def main():
    print("Building GRCh38 genomic → TOGA ENST CDS/AA position maps\n")

    # Load existing entries to allow idempotent re-runs
    existing: dict = {}
    if OUT_FILE.exists():
        with open(OUT_FILE) as f:
            existing = json.load(f)
        print(f"Loaded {len(existing)} existing entries from {OUT_FILE}\n")

    results: dict = dict(existing)

    for gene, enst_id in sorted(TOGA_ENST_GENES.items()):
        if gene in existing and existing[gene].get("enst") == enst_id:
            print(f"  [SKIP]  {gene:<14} {enst_id} — already built")
            continue

        print(f"  {gene:<14} {enst_id} ... ", end="", flush=True)
        meta, positions = build_genomic_cds_map(enst_id)
        time.sleep(0.2)   # Ensembl rate limit: ≤15 req/s

        if meta is None:
            print("FAILED")
            continue

        results[gene] = {
            **meta,
            # JSON keys must be strings; values are {"cds_pos": int, "aa_pos": int}
            "map": {str(k): v for k, v in positions.items()},
        }
        print(f"built  strand={meta['strand']:+d}  "
              f"prot_len={meta['prot_len']}aa  "
              f"cds_positions={meta['n_cds_positions']}")

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_FILE, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved {len(results)} gene maps → {OUT_FILE}")


if __name__ == "__main__":
    main()
