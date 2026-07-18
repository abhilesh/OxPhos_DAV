#!/usr/bin/env python3
"""
src/mutagenesis/00_prioritize_pairs.py

Score all 694 significant compensatory pairs for mutagenesis prioritization.

Scoring rubric (max 17 points):
  Phylogenetic significance : min(pagel_fdr, branch_cooccur_fdr) ≤ 0.01 → 3; ≤ 0.05 → 2; ≤ 0.10 → 1
  Contact class             : electrostatic → 3; hbond → 2; hydrophobic → 1; vdw → 0
  Physicochemical type      : charge_reversal → 3; charge_rescue → 2; volume_swap → 1; else → 0
  Epistatic pre-adaptation  : contact_first → 2; co_occurring → 1; else → 0
  Likely incompatible       : True → 2
  Convergence               : n_dar_gain_branches ≥ 10 → 2; ≥ 2 → 1; 1 → 0
  Confidence tier           : high_confidence → 1
  Contact proximity         : Cβ-Cβ distance < 5 Å → 1

Outputs:
  results/mutagenesis/prioritized_pairs.csv  — all pairs scored and ranked
  results/mutagenesis/top_targets.csv        — top 50 by priority score
  results/mutagenesis/top_targets.md         — human-readable table

Usage:
  docker run --rm -v $(pwd):/app oxphos_dav_analysis conda run -n oxphos_dav \\
      python src/mutagenesis/00_prioritize_pairs.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

# ─── Paths ────────────────────────────────────────────────────────────────────
PARTNERS_CSV   = ROOT / "results" / "structural" / "compensatory_partners.csv"
CONTACTS_CSV   = ROOT / "results" / "structural" / "dar_contacts_cbcb8A.csv"
STRUCTURES_DIR = ROOT / "data" / "structures"
OUT_DIR        = ROOT / "results" / "mutagenesis"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ─── Scoring maps ─────────────────────────────────────────────────────────────
CONTACT_CLASS_SCORE = {
    "electrostatic": 3,
    "hbond":         2,
    "hydrophobic":   1,
    "vdw":           0,
}

PHYSICO_SCORE = {
    "charge_reversal": 3,
    "charge_rescue":   2,
    "volume_swap":     1,
    "polarity_swap":   1,
    "same_direction":  0,
    "unclassified":    0,
}

TIMING_SCORE = {
    "contact_first": 2,
    "co_occurring":  1,
    "contact_after": 0,
    "no_contact_change": 0,
}

# ─── Scoring function ──────────────────────────────────────────────────────────
def phylo_score(row: pd.Series) -> int:
    fdrs = [row.get("pagel_fdr"), row.get("branch_cooccur_fdr")]
    fdrs = [float(v) for v in fdrs if pd.notna(v)]
    if not fdrs:
        # Fall back to fisher_fdr if neither Pagel nor branch is available
        fisher = row.get("fisher_fdr")
        if pd.notna(fisher):
            fdrs = [float(fisher)]
    if not fdrs:
        return 0
    best = min(fdrs)
    if best <= 0.01:
        return 3
    if best <= 0.05:
        return 2
    if best <= 0.10:
        return 1
    return 0


def score_pair(row: pd.Series, cbcb_dist: float) -> int:
    s = phylo_score(row)

    cc = str(row.get("contact_class", "")).lower().strip()
    s += CONTACT_CLASS_SCORE.get(cc, 0)

    pt = str(row.get("physicochemical_type", "")).lower().strip()
    s += PHYSICO_SCORE.get(pt, 0)

    dt = str(row.get("dominant_timing", "")).lower().strip()
    s += TIMING_SCORE.get(dt, 0)

    incompatible = row.get("likely_incompatible")
    if str(incompatible).lower() in ("true", "1", "yes"):
        s += 2

    n_gains = row.get("n_dar_gain_branches")
    if pd.notna(n_gains):
        n = int(float(n_gains))
        if n >= 10:
            s += 2
        elif n >= 2:
            s += 1

    tier = str(row.get("confidence_tier", "")).lower().strip()
    if tier == "high_confidence":
        s += 1

    if pd.notna(cbcb_dist) and cbcb_dist < 5.0:
        s += 1

    return s


# ─── Cβ distance from CIF ─────────────────────────────────────────────────────
_structure_cache: dict = {}


def load_structure(pdb_id: str):
    if pdb_id in _structure_cache:
        return _structure_cache[pdb_id]
    try:
        from Bio.PDB import MMCIFParser
    except ImportError:
        print("WARNING: biopython not available — skipping distance calculation")
        _structure_cache[pdb_id] = None
        return None
    cif_path = STRUCTURES_DIR / f"{pdb_id}.cif"
    if not cif_path.exists():
        _structure_cache[pdb_id] = None
        return None
    parser = MMCIFParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure(pdb_id, str(cif_path))
    _structure_cache[pdb_id] = structure
    return structure


def get_cb_vector(structure, chain_id: str, resnum):
    """Return Cβ vector (Cα for Gly) for a given chain and residue number."""
    try:
        model = structure[0]
        if chain_id not in model:
            return None
        chain = model[chain_id]
        resnum_int = int(float(resnum))
        for residue in chain:
            het, seq, icode = residue.id
            if het == " " and seq == resnum_int:
                atom_name = "CA" if residue.get_resname().strip() == "GLY" else "CB"
                if atom_name in residue:
                    return residue[atom_name].get_vector()
    except Exception:
        pass
    return None


def cbcb_dist(pdb_id: str, dar_chain: str, dar_res, contact_chain: str, contact_res) -> float:
    structure = load_structure(pdb_id)
    if structure is None:
        return np.nan
    v1 = get_cb_vector(structure, str(dar_chain).strip(), dar_res)
    v2 = get_cb_vector(structure, str(contact_chain).strip(), contact_res)
    if v1 is None or v2 is None:
        return np.nan
    return float((v1 - v2).norm())


# ─── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    print("Loading compensatory partners...")
    partners = pd.read_csv(PARTNERS_CSV, dtype={"ann_id": str})
    print(f"  {len(partners)} pairs")

    print("Loading structural contacts for coordinate join...")
    contacts = pd.read_csv(CONTACTS_CSV, dtype={"ann_id": str}, low_memory=False)

    # Normalise column name: dar_locus → dar_gene
    if "dar_locus" in contacts.columns and "dar_gene" not in contacts.columns:
        contacts = contacts.rename(columns={"dar_locus": "dar_gene"})

    coord_cols = [
        "dar_gene", "dar_aa_coord", "dar_alt_aa", "contact_gene", "contact_refseq_pos",
        "pdb_id", "dar_chain", "dar_struct_res", "contact_chain", "contact_resnum",
    ]
    contacts_sub = contacts[[c for c in coord_cols if c in contacts.columns]].drop_duplicates()

    # Ensure join columns are consistently typed: float → int → str (handles "45.0" == "45")
    for df in (partners, contacts_sub):
        for col in ("dar_aa_coord", "contact_refseq_pos"):
            if col in df.columns:
                df[col] = df[col].apply(
                    lambda v: str(int(float(v))) if pd.notna(v) else str(v)
                )

    # Join on physical residue positions only — ann_id formats differ between files
    join_keys = ["dar_gene", "dar_aa_coord", "dar_alt_aa", "contact_gene", "contact_refseq_pos"]
    merged = partners.merge(contacts_sub, on=join_keys, how="left")
    print(f"  After join: {len(merged)} rows (may include multiple PDB entries per pair)")

    # Pre-load structures
    pdb_ids = [p for p in merged["pdb_id"].dropna().unique() if str(p) != "nan"]
    print(f"Loading {len(pdb_ids)} CIF structures: {sorted(pdb_ids)}")
    for pid in sorted(pdb_ids):
        print(f"  {pid}...", end=" ", flush=True)
        s = load_structure(str(pid))
        print("OK" if s is not None else "NOT FOUND")

    # Compute Cβ distances
    print("Computing Cβ-Cβ distances...")
    dists = []
    for _, row in merged.iterrows():
        pid = row.get("pdb_id")
        if pd.isna(pid):
            dists.append(np.nan)
            continue
        d = cbcb_dist(
            str(pid),
            row["dar_chain"],
            row["dar_struct_res"],
            row["contact_chain"],
            row["contact_resnum"],
        )
        dists.append(d)
    merged["cbcb_dist_A"] = dists
    n_valid = sum(pd.notna(d) for d in dists)
    print(f"  {n_valid}/{len(dists)} distances computed successfully")

    # Per pair, keep the PDB entry with the smallest Cβ distance
    pair_key = ["ann_id", "dar_gene", "dar_aa_coord", "dar_alt_aa",
                "contact_gene", "contact_refseq_pos", "contact_alt_aa"]
    pair_key = [c for c in pair_key if c in merged.columns]

    dedup = (
        merged
        .sort_values("cbcb_dist_A", na_position="last")
        .groupby(pair_key, sort=False)
        .first()
        .reset_index()
    )
    print(f"  After deduplication (best distance per pair): {len(dedup)} pairs")

    # Score
    print("Scoring all pairs...")
    scores = [score_pair(row, row.get("cbcb_dist_A", np.nan)) for _, row in dedup.iterrows()]
    dedup["priority_score"] = scores
    dedup = dedup.sort_values("priority_score", ascending=False).reset_index(drop=True)
    dedup["priority_rank"] = range(1, len(dedup) + 1)

    # Add intraprotein flag
    dedup["is_intraprotein"] = dedup["dar_gene"] == dedup["contact_gene"]

    # ─── Save outputs ─────────────────────────────────────────────────────────
    out_all = OUT_DIR / "prioritized_pairs.csv"
    dedup.to_csv(out_all, index=False)
    print(f"\nSaved all {len(dedup)} scored pairs → {out_all}")

    top50 = dedup.head(50).copy()
    out_top = OUT_DIR / "top_targets.csv"
    top50.to_csv(out_top, index=False)
    print(f"Saved top 50 targets → {out_top}")

    # Markdown table
    md_display_cols = [
        "priority_rank", "priority_score",
        "dar_gene", "dar_aa_coord", "dar_ref_aa", "dar_alt_aa",
        "contact_gene", "contact_refseq_pos", "contact_human_aa", "contact_alt_aa",
        "contact_class", "contact_type", "is_intraprotein",
        "physicochemical_type", "dominant_timing",
        "likely_incompatible", "n_dar_gain_branches",
        "confidence_tier", "cbcb_dist_A", "pdb_id",
        "pagel_fdr", "branch_cooccur_fdr",
    ]
    md_display_cols = [c for c in md_display_cols if c in top50.columns]

    def fmt_cell(c: str, v) -> str:
        if pd.isna(v):
            return ""
        if c == "cbcb_dist_A":
            return f"{float(v):.2f}"
        if c in ("pagel_fdr", "branch_cooccur_fdr", "fisher_fdr"):
            return f"{float(v):.2e}"
        return str(v)

    score_dist = dedup["priority_score"].describe()
    lines = [
        "# Top Mutagenesis Targets",
        "",
        f"Ranked {len(dedup)} significant compensatory pairs by priority score (max 17 points).",
        "",
        f"Score distribution — mean: {score_dist['mean']:.1f}, "
        f"median: {score_dist['50%']:.1f}, "
        f"max: {score_dist['max']:.0f}, "
        f"min: {score_dist['min']:.0f}",
        "",
        f"Intraprotein pairs: {dedup['is_intraprotein'].sum()} "
        f"({100*dedup['is_intraprotein'].mean():.1f}%)",
        f"Interprotein pairs: {(~dedup['is_intraprotein']).sum()}",
        "",
        "| " + " | ".join(md_display_cols) + " |",
        "| " + " | ".join(["---"] * len(md_display_cols)) + " |",
    ]
    for _, row in top50.iterrows():
        cells = [fmt_cell(c, row.get(c)) for c in md_display_cols]
        lines.append("| " + " | ".join(cells) + " |")

    out_md = OUT_DIR / "top_targets.md"
    out_md.write_text("\n".join(lines) + "\n")
    print(f"Saved markdown table → {out_md}")

    # ─── Summary stats ────────────────────────────────────────────────────────
    print("\n── Priority score distribution (all pairs) ─────────────────")
    print(dedup["priority_score"].value_counts().sort_index(ascending=False).to_string())

    print("\n── Top 10 pairs ─────────────────────────────────────────────")
    preview_cols = [
        "priority_rank", "priority_score", "dar_gene", "dar_aa_coord",
        "dar_alt_aa", "contact_gene", "contact_refseq_pos", "contact_alt_aa",
        "contact_class", "physicochemical_type", "dominant_timing",
        "n_dar_gain_branches", "cbcb_dist_A",
    ]
    preview_cols = [c for c in preview_cols if c in dedup.columns]
    print(dedup.head(10)[preview_cols].to_string(index=False))

    n_intra = dedup["is_intraprotein"].sum()
    print(f"\nIntraprotein: {n_intra}  Interprotein: {len(dedup)-n_intra}")
    n_dist = dedup["cbcb_dist_A"].notna().sum()
    print(f"Cβ distances computed for {n_dist}/{len(dedup)} pairs")
    n_close = (dedup["cbcb_dist_A"] < 5.0).sum()
    print(f"Pairs with Cβ-Cβ < 5 Å: {n_close}")


if __name__ == "__main__":
    main()
