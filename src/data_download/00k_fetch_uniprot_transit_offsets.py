#!/usr/bin/env python3
"""
src/data_download/00k_fetch_uniprot_transit_offsets.py

For every nucDNA OXPHOS gene that has residue-anchoring failures in the
structure mapping stage (large_offset_candidate / mature_offset_candidate),
this script:

  1. Queries UniProt Swiss-Prot for the human reviewed entry and extracts the
     "Transit peptide" (or "Signal peptide") cleavage position.

  2. Validates the annotation by pairwise-aligning the TOGA human reference
     sequence (pre-protein) against the UniProt canonical sequence to confirm
     they represent the same protein variant.

  3. For validated genes (≥ 90% full-sequence identity, offset ≤ 100 aa),
     compares the observed offset distribution from the failure audit against
     the UniProt annotated offset to flag agreement or disagreement.

  4. Writes results to:
       data/reference/uniprot_transit_offsets.tsv   — annotated offset table
       data/reference/uniprot_transit_alignment_report.txt  — verbose details

  5. Appends validated entries to:
       data/reference/structural_anchor_exception_registry.tsv
     (existing entries are not overwritten; genes already present are skipped)

Run from project root:
  python src/data_download/00k_fetch_uniprot_transit_offsets.py

Requires: biopython, pandas (available in the oxphos_dav conda environment)
"""

import csv
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from utils.uniprot_transit_peptide import fetch_transit_annotation  # noqa: E402

TOGA_AA_DIR = ROOT / "data" / "alignments" / "toga_hg38_aa"
FAILURE_AUDIT = ROOT / "results" / "structural" / "structure_mapping_failure_audit.csv"
OUT_TSV = ROOT / "data" / "reference" / "uniprot_transit_offsets.tsv"
REPORT_TXT = ROOT / "data" / "reference" / "uniprot_transit_alignment_report.txt"
REGISTRY_TSV = ROOT / "data" / "reference" / "structural_anchor_exception_registry.tsv"

# Genes to query: all nucDNA genes that show up as mature_offset_candidate in
# the failure audit.  mtDNA genes (MT-*) are excluded — their "offsets" are
# unrelated to import transit peptides.
CANDIDATE_GENES = [
    # Complex I — matrix arm (iron-sulfur cluster assembly)
    "NDUFV1", "NDUFV2", "NDUFV3",
    "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS6", "NDUFS7", "NDUFS8",
    "NDUFAB1",
    # Complex I — peripheral arm
    "NDUFA2", "NDUFA6", "NDUFA9", "NDUFA10", "NDUFA13",
    # Complex I — membrane arm accessory
    "NDUFB2", "NDUFB3", "NDUFB5", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB11",
    "NDUFC1",
    # Complex II (SDH — already in registry, included for completeness)
    "SDHA", "SDHB", "SDHC", "SDHD",
    # Complex III
    "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRB", "UQCRH", "CYC1",
    # Complex IV
    "COX4I1", "COX4I2", "COX5A", "COX5B",
    "COX6A1", "COX6A2", "COX6B1", "COX6C",
    "COX7A1", "COX7A2", "COX7B", "COX7C", "COX8A",
    "NDUFA4",
    # Complex V
    "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
    "ATP5MC1", "ATP5MC2", "ATP5MC3",
    "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO",
]

# Identity threshold for accepting TOGA ↔ UniProt alignment as same isoform
SEQ_IDENTITY_THRESHOLD = 0.90

# Max reasonable transit peptide length (biologically motivated: longest known
# mammalian MTS is ~100 aa; flag anything beyond this for manual review)
MAX_TRANSIT_LEN = 100

# Buffer aa added to the UniProt offset when writing to the exception registry
# (provides slack for the sliding-window anchor search)
REGISTRY_BUFFER = 5


# ── Alignment helpers ──────────────────────────────────────────────────────────

def _make_aligner(mode: str = "global") -> PairwiseAligner:
    aln = PairwiseAligner()
    aln.mode = mode
    aln.match_score = 2
    aln.mismatch_score = -1
    aln.open_gap_score = -4
    aln.extend_gap_score = -0.5
    return aln


def global_identity(seq_a: str, seq_b: str) -> float:
    """Fraction of aligned columns that are identical (gaps excluded)."""
    if not seq_a or not seq_b:
        return 0.0
    aln = _make_aligner("global")
    try:
        best = next(iter(aln.align(seq_a, seq_b)))
    except StopIteration:
        return 0.0
    a_aln, b_aln = str(best[0]), str(best[1])
    matches = sum(a == b for a, b in zip(a_aln, b_aln) if a != "-" and b != "-")
    compared = sum(1 for a, b in zip(a_aln, b_aln) if a != "-" and b != "-")
    return matches / compared if compared else 0.0


def local_identity_at_start(query: str, target: str) -> float:
    """
    Local alignment of query against target.  Reports identity over the aligned
    region.  Used to verify that the mature sequence (query) aligns to the
    start of the PDB / TOGA sequence (target).
    """
    if not query or not target:
        return 0.0
    aln = _make_aligner("local")
    try:
        best = next(iter(aln.align(query, target)))
    except StopIteration:
        return 0.0
    a_aln, b_aln = str(best[0]), str(best[1])
    matches = sum(a == b for a, b in zip(a_aln, b_aln) if a != "-" and b != "-")
    compared = sum(1 for a, b in zip(a_aln, b_aln) if a != "-" and b != "-")
    return matches / compared if compared else 0.0


# ── TOGA sequence loading ──────────────────────────────────────────────────────

def load_toga_human_seqs(toga_dir: Path) -> dict[str, str]:
    """Return {gene: gap-free human reference sequence} from TOGA AA FASTAs."""
    seqs = {}
    if not toga_dir.exists():
        return seqs
    for fasta in sorted(toga_dir.glob("*_aa_alignment.fasta")):
        gene = fasta.name.replace("_aa_alignment.fasta", "")
        for rec in SeqIO.parse(fasta, "fasta"):
            if rec.id.startswith("Homo_sapiens|9606"):
                seqs[gene] = str(rec.seq).replace("-", "")
                break
    return seqs


# ── Failure audit summary ──────────────────────────────────────────────────────

def load_observed_offsets(audit_path: Path) -> dict[str, dict]:
    """
    Return {gene: {median_offset, min_offset, max_offset, n_rows}} from the
    structure_mapping_failure_audit.csv for mature_offset_candidate rows.
    """
    if not audit_path.exists():
        return {}
    df = pd.read_csv(audit_path)
    df = df[df["status_category"] == "mature_offset_candidate"].copy()
    df["offset_val"] = (
        df["status"]
        .str.extract(r"large_offset_candidate_([+-]\d+)")
        .squeeze()
        .astype(float)
    )
    df = df.dropna(subset=["offset_val"])
    result = {}
    for gene, grp in df.groupby("interpreted_gene"):
        vals = grp["offset_val"]
        result[gene] = {
            "n_rows": int(vals.count()),
            "median_offset": float(vals.median()),
            "min_offset": float(vals.min()),
            "max_offset": float(vals.max()),
        }
    return result


# ── Registry helpers ───────────────────────────────────────────────────────────

def load_registry(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    registry = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = (row.get("gene") or "").strip()
            if gene:
                registry[gene] = row
    return registry


def append_registry_entry(path: Path, gene: str, max_offset: int, rationale: str,
                           evidence: str, notes: str) -> None:
    """Append a new entry to the anchor exception registry TSV."""
    fieldnames = [
        "gene", "allow_extended_anchor", "max_offset",
        "rationale", "evidence_status", "notes",
    ]
    write_header = not path.exists() or path.stat().st_size == 0
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        if write_header:
            w.writeheader()
        w.writerow({
            "gene": gene,
            "allow_extended_anchor": "True",
            "max_offset": max_offset,
            "rationale": rationale,
            "evidence_status": evidence,
            "notes": notes,
        })


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("Loading TOGA human reference sequences...")
    toga_seqs = load_toga_human_seqs(TOGA_AA_DIR)
    print(f"  {len(toga_seqs)} genes found in TOGA alignment dir.")

    print("Loading observed failure offsets from structure mapping audit...")
    observed = load_observed_offsets(FAILURE_AUDIT)
    print(f"  {len(observed)} genes with mature_offset_candidate rows.")

    print("Loading existing anchor exception registry...")
    existing_registry = load_registry(REGISTRY_TSV)
    print(f"  {len(existing_registry)} existing entries: {list(existing_registry.keys())}")

    output_rows = []
    report_lines = []
    new_registry_genes = []

    total = len(CANDIDATE_GENES)
    for idx, gene in enumerate(CANDIDATE_GENES, 1):
        print(f"  [{idx:3d}/{total}] {gene:<16}", end=" ", flush=True)
        time.sleep(0.3)  # polite rate-limiting for UniProt API

        annotation = fetch_transit_annotation(gene)
        if annotation is None:
            print("no UniProt entry found")
            output_rows.append({
                "gene": gene, "uniprot_accession": "", "uniprot_gene_name": "",
                "feature_type": "", "transit_end": "", "mature_start": "",
                "uniprot_offset": "",
                "toga_seq_len": len(toga_seqs.get(gene, "")),
                "toga_uniprot_full_identity": "",
                "mature_start_local_identity": "",
                "validated": "no_entry",
                "observed_n_rows": observed.get(gene, {}).get("n_rows", 0),
                "observed_median_offset": observed.get(gene, {}).get("median_offset", ""),
                "observed_min_offset": observed.get(gene, {}).get("min_offset", ""),
                "observed_max_offset": observed.get(gene, {}).get("max_offset", ""),
                "offset_audit_agreement": "",
                "notes": "no reviewed UniProt entry",
            })
            continue

        accession = annotation["accession"]
        transit_end = annotation["transit_end"]
        mature_start = annotation["mature_start"]
        uniprot_seq = annotation["sequence"]
        evidence = annotation["evidence"] or ""
        feature_type = annotation["feature_type"] or ""
        toga_seq = toga_seqs.get(gene, "")

        # ── Validate: TOGA vs UniProt full sequence identity ──────────────────
        full_id = global_identity(toga_seq, uniprot_seq) if toga_seq else 0.0

        # ── Validate: mature sequence aligns to TOGA ──────────────────────────
        mature_id = 0.0
        if transit_end and toga_seq:
            mature_seq_uniprot = uniprot_seq[transit_end:]
            mature_id = local_identity_at_start(mature_seq_uniprot, toga_seq)

        # ── Agreement with observed offset distribution ────────────────────────
        obs = observed.get(gene, {})
        offset_agreement = ""
        if transit_end and obs:
            obs_med = obs.get("median_offset", 0)
            delta = abs(obs_med - transit_end)
            if delta <= 5:
                offset_agreement = f"agree_delta{delta:.0f}"
            elif delta <= 15:
                offset_agreement = f"close_delta{delta:.0f}"
            else:
                offset_agreement = f"mismatch_delta{delta:.0f}"
        elif transit_end and not obs:
            offset_agreement = "no_audit_data"
        elif not transit_end:
            offset_agreement = "no_transit_annotated"

        # ── Validation decision ───────────────────────────────────────────────
        if not transit_end:
            validated = "no_transit"
        elif transit_end > MAX_TRANSIT_LEN:
            validated = f"transit_too_long_{transit_end}"
        elif full_id < SEQ_IDENTITY_THRESHOLD:
            validated = f"seq_mismatch_{full_id:.2f}"
        elif "mismatch" in offset_agreement:
            validated = f"offset_mismatch_obs{obs.get('median_offset', '?'):.0f}_uni{transit_end}"
        elif mature_id < 0.80:
            validated = f"mature_alignment_weak_{mature_id:.2f}"
        else:
            validated = "validated"

        print(
            f"accession={accession}  transit_end={transit_end}  "
            f"full_id={full_id:.2f}  mature_id={mature_id:.2f}  "
            f"agreement={offset_agreement}  → {validated}"
        )

        output_rows.append({
            "gene": gene,
            "uniprot_accession": accession,
            "uniprot_gene_name": annotation["gene_name"],
            "feature_type": feature_type,
            "transit_end": transit_end if transit_end else "",
            "mature_start": mature_start if transit_end else "",
            "uniprot_offset": transit_end if transit_end else "",
            "toga_seq_len": len(toga_seq),
            "toga_uniprot_full_identity": f"{full_id:.4f}" if toga_seq else "",
            "mature_start_local_identity": f"{mature_id:.4f}" if transit_end and toga_seq else "",
            "validated": validated,
            "observed_n_rows": obs.get("n_rows", 0),
            "observed_median_offset": obs.get("median_offset", ""),
            "observed_min_offset": obs.get("min_offset", ""),
            "observed_max_offset": obs.get("max_offset", ""),
            "offset_audit_agreement": offset_agreement,
            "notes": f"evidence={evidence}" if evidence else "",
        })

        report_lines.append(
            f"\n{'='*70}\n"
            f"Gene: {gene}  ({accession})\n"
            f"Feature: {feature_type}\n"
            f"Transit end: {transit_end}  Mature start: {mature_start}\n"
            f"UniProt evidence: {evidence}\n"
            f"TOGA seq len: {len(toga_seq)}  UniProt seq len: {len(uniprot_seq)}\n"
            f"Full-sequence identity (TOGA vs UniProt): {full_id:.4f}\n"
            f"Mature-sequence local identity: {mature_id:.4f}\n"
            f"Observed failure offsets (n={obs.get('n_rows',0)}): "
            f"median={obs.get('median_offset','?')}  "
            f"range=[{obs.get('min_offset','?')}, {obs.get('max_offset','?')}]\n"
            f"Offset agreement: {offset_agreement}\n"
            f"Validated: {validated}\n"
        )

        # ── Update registry for validated genes ───────────────────────────────
        if validated == "validated" and gene not in existing_registry:
            # max_offset must cover the transit peptide AND any PDB N-terminal
            # resolution gap.  Use the larger of: UniProt transit_end+buffer or
            # the maximum observed positive offset from the failure audit+buffer.
            obs_max_pos = obs.get("max_offset", 0)
            obs_max_pos = int(obs_max_pos) if obs_max_pos and float(obs_max_pos) > 0 else 0
            max_offset_for_registry = max(transit_end + REGISTRY_BUFFER,
                                          obs_max_pos + REGISTRY_BUFFER)
            new_registry_genes.append(gene)
            append_registry_entry(
                REGISTRY_TSV,
                gene=gene,
                max_offset=max_offset_for_registry,
                rationale="validated_transit_peptide_offset",
                evidence=f"uniprot_{accession}_transit_end_{transit_end}",
                notes=(
                    f"UniProt {accession} annotates transit peptide end at "
                    f"position {transit_end}; validated by TOGA sequence alignment "
                    f"(full_id={full_id:.2f}, mature_id={mature_id:.2f}); "
                    f"offset_agreement={offset_agreement}"
                ),
            )
            print(f"    → Added to anchor registry (max_offset={max_offset_for_registry})")
        elif validated == "validated" and gene in existing_registry:
            print(f"    → Already in registry, skipped.")

    # ── Write output TSV ───────────────────────────────────────────────────────
    fieldnames = [
        "gene", "uniprot_accession", "uniprot_gene_name", "feature_type",
        "transit_end", "mature_start", "uniprot_offset",
        "toga_seq_len", "toga_uniprot_full_identity", "mature_start_local_identity",
        "validated",
        "observed_n_rows", "observed_median_offset", "observed_min_offset", "observed_max_offset",
        "offset_audit_agreement", "notes",
    ]
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(output_rows)
    print(f"\nOffset table written: {OUT_TSV}  ({len(output_rows)} rows)")

    with open(REPORT_TXT, "w") as fh:
        fh.write("UniProt Transit Peptide Annotation Report\n")
        fh.write("=" * 70 + "\n")
        fh.write(f"Genes queried: {len(CANDIDATE_GENES)}\n")
        validated_list = [r["gene"] for r in output_rows if r["validated"] == "validated"]
        fh.write(f"Validated: {len(validated_list)}\n")
        fh.write(f"New registry entries added: {len(new_registry_genes)}\n")
        fh.write("\nNew registry entries:\n")
        for g in new_registry_genes:
            fh.write(f"  {g}\n")
        fh.write("\n")
        fh.writelines(report_lines)
    print(f"Alignment report written: {REPORT_TXT}")

    print(f"\nNew anchor registry entries added: {len(new_registry_genes)}")
    if new_registry_genes:
        print("  " + ", ".join(new_registry_genes))

    print("\nDone.")


if __name__ == "__main__":
    main()
