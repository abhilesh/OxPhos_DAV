#!/usr/bin/env python3
"""
src/structural/01_find_compensating_partners.py

For every non-discarded c-DAR (Tier A, B, C), test whether the structurally
contacting residues show co-evolutionary enrichment in the species that carry
the disease amino acid.

Logic:
  - cDAR species   : those listed in var["lineages_with_disease_allele"] — they carry the
                     human-pathogenic AA as wild-type.
  - Background     : all other species in the alignment that have a readable AA
                     at the DAR position (not gap / X / ! / *).
  - For each contact residue: one-sided Fisher's exact test per alternative AA,
    asking "is this AA enriched in cDAR-harboring species?"
  - Benjamini-Hochberg FDR applied across all tests per DAR.

Structural incompatibility flag:
  A contact pair (dar_alt_aa, contact_human_aa) is flagged as
  "likely_incompatible" when:
    - contact_class is hbond or electrostatic  (charge/H-bond network sensitive)
    - sensitivity ≥ 0.5 AND specificity ≥ 0.5  (strong concordance)
  These are candidates where the disease AA physically clashes with the human
  contact AA, and compensating species resolve it by changing the contact residue.

Concordance metrics per candidate pair:
  sensitivity  = a / (a + c)   fraction of cDAR species carrying the alt contact AA
  specificity  = a / (a + b)   fraction of alt-contact-AA species that are cDAR
  odds_ratio                   from Fisher's exact

Output:
  results/structural/compensatory_partners.csv   — all significant pairs (FDR ≤ 0.10)
  results/structural/concordance_summary.csv     — per-cDAR enrichment counts by tier

Run from project root:
  python src/structural/01_find_compensating_partners.py
"""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from Bio import SeqIO
from scipy.stats import fisher_exact

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parents[2]
CURATED_DIR = ROOT / "data" / "annotations" / "curated"
MT_JSON     = CURATED_DIR / "cdav_classifications_mtDNA.json"
NUC_JSON    = CURATED_DIR / "cdav_classifications_nucDNA.json"
CONTACTS    = ROOT / "results" / "structural" / "dar_contacts_cbcb8A.csv"
TOGA_AA_DIR = ROOT / "data" / "alignments" / "toga_hg38_aa"
MT_AA_DIR   = ROOT / "data" / "alignments" / "mtdna_aa"
OUT_DIR     = ROOT / "results" / "structural"

_MASK            = {"-", "X", "!", "*"}
_INCOMP_CLASSES  = {"hbond", "electrostatic"}   # contact types sensitive to AA identity
_INCOMP_MIN_SENS = 0.5
_INCOMP_MIN_SPEC = 0.5

# TOGA filename → canonical gene name (mirrors 00_map_davs_to_structure.py)
_TOGA_TO_CANONICAL = {"COXFA4": "NDUFA4"}

# Tier display order for summary table
_TIER_ORDER = {"Tier A": 0, "Tier B": 1, "Tier C": 2}


# ── Alignment loading ─────────────────────────────────────────────────────────

def load_alignments(genes: set) -> dict[str, dict[str, str]]:
    """
    Load per-gene AA alignments for the requested gene set.
    Returns {gene: {species_shortname: gapped_sequence}}.
    Species shortname = header.split("|")[0].
    """
    alns: dict[str, dict[str, str]] = {}
    for aln_dir in (TOGA_AA_DIR, MT_AA_DIR):
        if not aln_dir.exists():
            continue
        for fasta in aln_dir.glob("*_aa_alignment.fasta"):
            gene = fasta.name.replace("_aa_alignment.fasta", "")
            gene = _TOGA_TO_CANONICAL.get(gene, gene)
            if gene not in genes:
                continue
            seqs = {}
            for rec in SeqIO.parse(fasta, "fasta"):
                sp = rec.id.split("|")[0]
                seqs[sp] = str(rec.seq).upper()
            alns[gene] = seqs
    return alns


def build_pos_to_col(ref_seq: str) -> dict[int, int]:
    """Maps 1-based ungapped biological position → 0-based alignment column."""
    pos_map: dict[int, int] = {}
    pos = 0
    for col, ch in enumerate(ref_seq):
        if ch not in _MASK:
            pos += 1
            pos_map[pos] = col
    return pos_map


def get_ref_seq(aln: dict[str, str]) -> tuple[str, str]:
    """Return (species_name, gapped_sequence) for Homo_sapiens."""
    for sp, seq in aln.items():
        if sp == "Homo_sapiens":
            return sp, seq
    raise ValueError("Human reference not found in alignment.")


# ── BH FDR ────────────────────────────────────────────────────────────────────

def bh_fdr(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values."""
    n = len(p_values)
    if n == 0:
        return []
    order     = np.argsort(p_values)
    adj       = np.array(p_values, dtype=float)
    for rank, idx in enumerate(order, 1):
        adj[idx] = min(1.0, p_values[idx] * n / rank)
    for i in range(n - 2, -1, -1):
        adj[order[i]] = min(adj[order[i]], adj[order[i + 1]])
    return adj.tolist()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # ── Load all non-discarded c-DARs ─────────────────────────────────────────
    cdars: dict[str, dict] = {}
    for json_path in (MT_JSON, NUC_JSON):
        if not json_path.exists():
            continue
        for var in json.load(open(json_path)):
            tier = var.get("tier", "")
            if "Discarded" in tier:
                continue
            if not var.get("is_cdav_amino_acid"):
                continue
            cdars[var["ann_id"]] = var

    by_tier = Counter(v["tier"] for v in cdars.values())
    print(f"Non-discarded c-DARs: {len(cdars)}")
    for t in sorted(by_tier, key=lambda x: _TIER_ORDER.get(x, 99)):
        print(f"  {t}: {by_tier[t]}")

    if not cdars:
        print("No c-DARs found.")
        return

    # ── Load contacts for relevant DARs ───────────────────────────────────────
    contacts: list[dict] = []
    with open(CONTACTS) as f:
        for row in csv.DictReader(f):
            if row["ann_id"] in cdars:
                contacts.append(row)
    print(f"Contacts for these DARs: {len(contacts)}")

    # ── Identify all genes needed and load alignments ─────────────────────────
    needed_genes = {
        g for row in contacts
        for g in (row["dar_locus"], row["contact_gene"])
        if not g.startswith("Unknown(")
    }
    alns = load_alignments(needed_genes)
    print(f"Alignments loaded for {len(alns)} genes.")

    # ── Per-DAR Fisher's exact tests ──────────────────────────────────────────
    by_dar: dict[str, list] = defaultdict(list)
    for row in contacts:
        by_dar[row["ann_id"]].append(row)

    all_candidates: list[dict] = []
    summary_rows:   list[dict] = []

    for ann_id, dar_contacts in by_dar.items():
        var          = cdars[ann_id]
        dar_gene     = var["locus"].split("/")[0]
        tier         = var["tier"]
        cdar_spp     = set(var.get("lineages_with_disease_allele", []))
        dar_aa_coord = int(dar_contacts[0]["dar_aa_coord"])

        if dar_gene not in alns:
            continue

        dar_aln = alns[dar_gene]
        try:
            _, dar_ref_gapped = get_ref_seq(dar_aln)
        except ValueError:
            continue
        dar_col_map = build_pos_to_col(dar_ref_gapped)
        if dar_aa_coord not in dar_col_map:
            continue
        dar_col = dar_col_map[dar_aa_coord]

        readable_spp = {
            sp for sp, seq in dar_aln.items()
            if sp != "Homo_sapiens"
            and len(seq) > dar_col
            and seq[dar_col] not in _MASK
        }
        cdar_in_aln = cdar_spp & readable_spp
        bg_spp      = readable_spp - cdar_spp

        if not cdar_in_aln:
            continue

        # De-duplicate by (contact_gene, contact_refseq_pos)
        seen_contacts: set = set()
        raw_tests: list[dict] = []

        for row in dar_contacts:
            contact_gene = row["contact_gene"]
            if contact_gene.startswith("Unknown(") or contact_gene not in alns:
                continue
            try:
                contact_pos = int(row["contact_refseq_pos"])
            except (ValueError, TypeError):
                continue

            key = (contact_gene, contact_pos)
            if key in seen_contacts:
                continue
            seen_contacts.add(key)

            c_aln = alns[contact_gene]
            try:
                _, c_ref_gapped = get_ref_seq(c_aln)
            except ValueError:
                continue
            c_col_map = build_pos_to_col(c_ref_gapped)
            if contact_pos not in c_col_map:
                continue
            c_col = c_col_map[contact_pos]

            human_contact_aa = row["contact_aa"]
            contact_class    = row["contact_class"]
            contact_type     = row["contact_type"]

            cdar_contact_aas = [
                c_aln[sp][c_col]
                for sp in cdar_in_aln
                if sp in c_aln
                and len(c_aln[sp]) > c_col
                and c_aln[sp][c_col] not in _MASK
            ]
            bg_contact_aas = [
                c_aln[sp][c_col]
                for sp in bg_spp
                if sp in c_aln
                and len(c_aln[sp]) > c_col
                and c_aln[sp][c_col] not in _MASK
            ]

            if not cdar_contact_aas:
                continue

            for alt_aa in {aa for aa in cdar_contact_aas if aa != human_contact_aa}:
                a = cdar_contact_aas.count(alt_aa)
                c = len(cdar_contact_aas) - a
                b = bg_contact_aas.count(alt_aa)
                d = len(bg_contact_aas) - b

                if a == 0:
                    continue

                _, p_val    = fisher_exact([[a, b], [c, d]], alternative="greater")
                sensitivity = a / (a + c) if (a + c) > 0 else 0.0
                specificity = a / (a + b) if (a + b) > 0 else 0.0

                incompatible = (
                    contact_class in _INCOMP_CLASSES
                    and sensitivity >= _INCOMP_MIN_SENS
                    and specificity >= _INCOMP_MIN_SPEC
                )

                raw_tests.append({
                    "ann_id":             ann_id,
                    "tier":               tier,
                    "dar_gene":           dar_gene,
                    "dar_genome":         var.get("genome", ""),
                    "dar_aa_coord":       dar_aa_coord,
                    "dar_ref_aa":         var.get("ref_aa", ""),
                    "dar_alt_aa":         var.get("alt_aa", ""),
                    "contact_gene":       contact_gene,
                    "contact_refseq_pos": contact_pos,
                    "contact_human_aa":   human_contact_aa,
                    "contact_alt_aa":     alt_aa,
                    "contact_class":      contact_class,
                    "contact_type":       contact_type,
                    "n_cdar_spp":         len(cdar_in_aln),
                    "n_bg_spp":           len(bg_spp),
                    "n_cdar_with_alt":    a,
                    "n_bg_with_alt":      b,
                    "sensitivity":        sensitivity,
                    "specificity":        specificity,
                    "p_value":            p_val,
                    "fdr":                None,
                    "likely_incompatible": incompatible,
                })

        if not raw_tests:
            continue

        fdrs = bh_fdr([t["p_value"] for t in raw_tests])
        for t, fdr in zip(raw_tests, fdrs):
            t["fdr"] = fdr

        all_candidates.extend(raw_tests)

        sig = [t for t in raw_tests if t["fdr"] <= 0.10]
        summary_rows.append({
            "ann_id":            ann_id,
            "tier":              tier,
            "dar_gene":          dar_gene,
            "dar_genome":        var.get("genome", ""),
            "dar_aa_coord":      dar_aa_coord,
            "dar_ref_aa":        var.get("ref_aa", ""),
            "dar_alt_aa":        var.get("alt_aa", ""),
            "n_cdar_spp":        len(cdar_in_aln),
            "n_bg_spp":          len(bg_spp),
            "n_contacts_tested": len(raw_tests),
            "n_sig_contacts":    len(sig),
            "n_incompatible":    sum(1 for t in sig if t["likely_incompatible"]),
            "intra_genomic_sig": sum(1 for t in sig if t["contact_type"] != "mt-nuc"),
            "inter_genomic_sig": sum(1 for t in sig if t["contact_type"] == "mt-nuc"),
            "sig_contact_genes": ",".join(sorted({t["contact_gene"] for t in sig})),
        })

    # ── Write outputs ──────────────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sig_candidates = sorted(
        [t for t in all_candidates if t.get("fdr") is not None and t["fdr"] <= 0.10],
        key=lambda x: (_TIER_ORDER.get(x["tier"], 99), x["ann_id"], x["fdr"]),
    )

    pair_fields = [
        "ann_id", "tier", "dar_gene", "dar_genome",
        "dar_aa_coord", "dar_ref_aa", "dar_alt_aa",
        "contact_gene", "contact_refseq_pos", "contact_human_aa", "contact_alt_aa",
        "contact_class", "contact_type", "likely_incompatible",
        "n_cdar_spp", "n_bg_spp", "n_cdar_with_alt", "n_bg_with_alt",
        "sensitivity", "specificity", "p_value", "fdr",
    ]
    with open(OUT_DIR / "compensatory_partners.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=pair_fields)
        w.writeheader()
        for row in sig_candidates:
            row["p_value"]     = f"{row['p_value']:.3e}"
            row["fdr"]         = f"{row['fdr']:.3e}"
            row["sensitivity"] = f"{row['sensitivity']:.3f}"
            row["specificity"] = f"{row['specificity']:.3f}"
            w.writerow(row)

    summary_fields = [
        "ann_id", "tier", "dar_gene", "dar_genome",
        "dar_aa_coord", "dar_ref_aa", "dar_alt_aa",
        "n_cdar_spp", "n_bg_spp",
        "n_contacts_tested", "n_sig_contacts",
        "n_incompatible", "intra_genomic_sig", "inter_genomic_sig",
        "sig_contact_genes",
    ]
    with open(OUT_DIR / "concordance_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary_fields)
        w.writeheader()
        w.writerows(
            sorted(summary_rows,
                   key=lambda r: (_TIER_ORDER.get(r["tier"], 99), -r["n_sig_contacts"]))
        )

    # ── Summary ───────────────���────────────────────────────────────────────────
    sig_by_tier    = Counter(t["tier"]         for t in sig_candidates)
    intra_by_tier  = Counter(t["tier"] for t in sig_candidates if t["contact_type"] != "mt-nuc")
    inter_by_tier  = Counter(t["tier"] for t in sig_candidates if t["contact_type"] == "mt-nuc")
    incomp_by_tier = Counter(t["tier"] for t in sig_candidates if t["likely_incompatible"] == True)

    print(f"\n{'='*65}")
    print("COMPENSATORY PARTNER ANALYSIS — ALL TIERS")
    print(f"{'='*65}")
    print(f"{'Tier':<10} {'cDARs':>7} {'sig pairs':>10} {'intra':>7} {'inter':>7} {'incompatible':>13}")
    print("-" * 60)
    for t in ("Tier A", "Tier B", "Tier C"):
        n_cdars = by_tier.get(t, 0)
        print(
            f"  {t:<8} {n_cdars:>7} {sig_by_tier[t]:>10} "
            f"{intra_by_tier[t]:>7} {inter_by_tier[t]:>7} {incomp_by_tier[t]:>13}"
        )
    print("-" * 60)
    print(
        f"  {'Total':<8} {len(cdars):>7} {len(sig_candidates):>10} "
        f"{sum(intra_by_tier.values()):>7} {sum(inter_by_tier.values()):>7} "
        f"{sum(incomp_by_tier.values()):>13}"
    )

    total_sig = len(sig_candidates)
    if total_sig:
        intra_frac = sum(intra_by_tier.values()) / total_sig * 100
        inter_frac = sum(inter_by_tier.values()) / total_sig * 100
        print(f"\nCompensation is {intra_frac:.0f}% intra-genomic, {inter_frac:.0f}% inter-genomic (mt↔nuc)")

    n_with_sig = sum(1 for r in summary_rows if r["n_sig_contacts"] > 0)
    print(f"\nDARs with ≥1 significant contact: {n_with_sig} / {len(summary_rows)}")
    print(f"\nSignificant pairs   : {OUT_DIR / 'compensatory_partners.csv'}")
    print(f"Concordance summary : {OUT_DIR / 'concordance_summary.csv'}")


if __name__ == "__main__":
    main()