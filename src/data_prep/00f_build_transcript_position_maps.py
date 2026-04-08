"""
src/data_prep/00f_build_transcript_position_maps.py

Builds NM_ → ENST amino-acid position maps for nuclear OXPHOS genes where the
ClinVar NM_ transcript and the TOGA ENST transcript encode different isoforms.

Why this is needed
──────────────────
TOGA publishes one alignment per gene; the ENST it chose may differ from the
NM_ used by ClinVar.  When the two isoforms have different CDS starts, internal
insertions, or entirely different N-termini, every c. coordinate from ClinVar
lands at the wrong position in the TOGA alignment.  The local ±10-residue anchor
strategy in AlignmentParser.find_sequence_anchor fails when the offset is large
(e.g. NDUFS1: 15 aa) or when the shared region is absent in that window.

Approach
────────
1. For each nuclear gene: fetch the NM_ protein sequence from NCBI Entrez and
   extract the ENST protein from the TOGA AA alignment (already on disk).
2. If the two sequences are identical → record as "identity" (no remapping needed).
3. If they differ → run a global pairwise alignment (BLOSUM62, affine gap penalties)
   and build a complete nm_aa_pos → enst_aa_pos map.  Positions that align to a
   gap in the ENST (inserted residues unique to NM_) map to None.
4. Save the result as data/reference/transcript_position_maps.json.

AlignmentParser reads this file at class initialisation and uses it in
check_compensation as a direct lookup, replacing find_sequence_anchor for genes
that have a map entry.

Run from project root inside the Docker container:
    python src/data_prep/00f_build_transcript_position_maps.py
"""

import json
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

from utils.parsers import GeneReference
from utils.utils import get_latest

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
TOGA_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_aa"
CURATED_DIR = DATA_DIR / "annotations" / "curated"
REF_DIR = DATA_DIR / "reference"
OUT_FILE = REF_DIR / "transcript_position_maps.json"

Entrez.email = "oxphos_dav@analysis.local"

# ── Aligner ────────────────────────────────────────────────────────────────────
# BLOSUM62 with affine gap penalties; designed for closely-related isoforms.
_aligner = PairwiseAligner()
_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
_aligner.mode = "global"
_aligner.open_gap_score = -10
_aligner.extend_gap_score = -0.5


# ── Sequence helpers ───────────────────────────────────────────────────────────

def _enst_protein(gene: str) -> tuple[str, str]:
    """
    Returns (enst_id, protein_seq) from the TOGA AA alignment for the gene.
    Gaps in the alignment are stripped.  Returns ("", "") if no file or no
    human reference sequence.
    """
    aa_fasta = TOGA_AA_DIR / f"{gene}_aa_alignment.fasta"
    if not aa_fasta.exists():
        return "", ""
    for rec in SeqIO.parse(aa_fasta, "fasta"):
        if rec.id.startswith("Homo_sapiens"):
            seq = str(rec.seq).replace("-", "").replace("*", "")
            m = re.search(r"(ENST\d+)", rec.id)
            enst_id = m.group(1) if m else "UNKNOWN"
            return enst_id, seq
    return "", ""


def _nm_protein(nm_id: str) -> str:
    """
    Fetches the CDS translation for an NM_ accession from NCBI.
    Returns the protein sequence string, or "" on failure.
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=nm_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        for feat in record.features:
            if feat.type == "CDS":
                aa = feat.qualifiers.get("translation", [""])[0]
                if aa:
                    return aa
    except Exception as e:
        print(f"    Entrez error for {nm_id}: {e}")
    return ""


# ── Position map builder ───────────────────────────────────────────────────────

def _build_pos_map(nm_seq: str, enst_seq: str) -> dict[int, int | None]:
    """
    Global pairwise alignment of nm_seq vs enst_seq.

    Returns {nm_aa_pos (1-indexed): enst_aa_pos (1-indexed) | None}.

    A position maps to None when it aligns to a gap in the ENST — meaning that
    residue exists only in the NM_ isoform and has no counterpart in TOGA.
    """
    alignments = _aligner.align(nm_seq, enst_seq)
    best = next(iter(alignments))   # highest-score alignment

    # Reconstruct aligned strings from the alignment object
    aligned_nm   = str(best[0])
    aligned_enst = str(best[1])

    pos_map: dict[int, int | None] = {}
    nm_pos = enst_pos = 0

    for nm_char, enst_char in zip(aligned_nm, aligned_enst):
        nm_gap   = (nm_char   == "-")
        enst_gap = (enst_char == "-")

        if not nm_gap:
            nm_pos += 1
        if not enst_gap:
            enst_pos += 1

        if not nm_gap:
            pos_map[nm_pos] = None if enst_gap else enst_pos

    return pos_map


def _identity_fraction(nm_seq: str, enst_seq: str, pos_map: dict) -> float:
    """Fraction of NM_ positions that map to the same residue in ENST."""
    matched = sum(
        1 for nm_pos, enst_pos in pos_map.items()
        if enst_pos is not None
        and nm_pos <= len(nm_seq)
        and enst_pos <= len(enst_seq)
        and nm_seq[nm_pos - 1] == enst_seq[enst_pos - 1]
    )
    return matched / len(nm_seq) if nm_seq else 0.0


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("Building NM_ → ENST amino-acid position maps for nuclear OXPHOS genes\n")

    # Load gene list
    hgnc_file = get_latest(DATA_DIR, "Canonical_OXPHOS_Subunits_HGNC*.csv")
    hgnc_ref  = GeneReference(hgnc_file)
    target_genes = sorted({
        data["symbol"]
        for data in hgnc_ref.lookup.values()
        if data.get("symbol") and not data["symbol"].startswith("MT-")
    })

    # Load ClinVar NM_ per gene from curated data
    nm_per_gene: dict[str, str] = {}
    nuc_curated = CURATED_DIR / "nucDNA_annotations_v2.json"
    if nuc_curated.exists():
        nm_counts: dict[str, Counter] = defaultdict(Counter)
        with open(nuc_curated) as f:
            for v in json.load(f):
                gene = v.get("locus", "")
                tx   = v.get("transcript_id", "").split(".")[0]
                if gene and tx.startswith("NM_"):
                    nm_counts[gene][tx] += 1
        nm_per_gene = {g: c.most_common(1)[0][0] for g, c in nm_counts.items()}
    print(f"NM_ transcripts found for {len(nm_per_gene)} genes in curated data.\n")

    # Load existing maps (to allow re-runs without re-fetching)
    existing: dict = {}
    if OUT_FILE.exists():
        with open(OUT_FILE) as f:
            existing = json.load(f)
        print(f"Loaded {len(existing)} existing map entries from {OUT_FILE}\n")

    results: dict = dict(existing)
    n_identity = n_mapped = n_skipped = n_failed = 0

    for gene in target_genes:
        nm_id = nm_per_gene.get(gene)
        if not nm_id:
            continue  # no ClinVar variants for this gene

        # Skip if already processed and NM_ hasn't changed
        if gene in existing and existing[gene].get("nm") == nm_id:
            n_skipped += 1
            continue

        enst_id, enst_seq = _enst_protein(gene)
        if not enst_seq:
            print(f"  [SKIP]   {gene:<15} no TOGA alignment")
            continue

        print(f"  {gene:<15} {nm_id:<16} fetching from NCBI...", end=" ", flush=True)
        nm_seq = _nm_protein(nm_id)
        time.sleep(0.35)  # NCBI rate limit: ~3 req/s

        if not nm_seq:
            print("FETCH FAILED")
            n_failed += 1
            continue

        if nm_seq == enst_seq:
            print(f"identical ({len(nm_seq)} aa)")
            results[gene] = {
                "nm": nm_id, "enst": enst_id,
                "type": "identity",
                "nm_len": len(nm_seq), "enst_len": len(enst_seq),
                "map": {},
            }
            n_identity += 1
            continue

        # Sequences differ — build global alignment map
        pos_map  = _build_pos_map(nm_seq, enst_seq)
        id_frac  = _identity_fraction(nm_seq, enst_seq, pos_map)
        n_mapped_pos = sum(1 for v in pos_map.values() if v is not None)
        n_gap_pos    = sum(1 for v in pos_map.values() if v is None)

        print(
            f"mapped  nm={len(nm_seq)}aa  enst={len(enst_seq)}aa  "
            f"id={id_frac:.1%}  "
            f"mapped={n_mapped_pos}  gaps={n_gap_pos}"
        )
        results[gene] = {
            "nm": nm_id, "enst": enst_id,
            "type": "mapped",
            "nm_len": len(nm_seq), "enst_len": len(enst_seq),
            "identity_fraction": round(id_frac, 4),
            "map": {str(k): v for k, v in pos_map.items()},  # JSON keys must be str
        }
        n_mapped += 1

    # Save
    REF_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_FILE, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n{'='*55}")
    print("POSITION MAP SUMMARY")
    print(f"{'='*55}")
    print(f"  Identity (no remapping needed) : {n_identity + n_skipped}")
    print(f"  Mapped (alignment built)       : {n_mapped}")
    print(f"  Fetch failures                 : {n_failed}")
    print(f"  Total genes in map file        : {len(results)}")
    print(f"\nSaved → {OUT_FILE}")

    # Print the mapped genes as a quick review
    mapped_genes = [(g, d) for g, d in results.items() if d["type"] == "mapped"]
    if mapped_genes:
        print(f"\nGenes with non-identity maps ({len(mapped_genes)}):")
        print(f"  {'Gene':<15} {'NM_':<16} {'ENST':<20} {'id%':<8} {'nm_len':<8} {'enst_len'}")
        for g, d in sorted(mapped_genes):
            print(
                f"  {g:<15} {d['nm']:<16} {d['enst']:<20} "
                f"{d.get('identity_fraction',0):.1%}    "
                f"{d['nm_len']:<8} {d['enst_len']}"
            )


if __name__ == "__main__":
    main()
