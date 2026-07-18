#!/usr/bin/env python3
"""
src/mutagenesis/03_compile_targets.py

Merge prioritization scores, FoldX ΔΔG, and MI/APC scores into a final
composite ranking. Generates the mutagenesis_report.md for experimental use.

Composite score additions:
  +2 if ddg_rescue_stab > 0.5 AND ddg_stab_dar > 0.5  (FoldX confirms rescue)
  +1 if mi_percentile > 75                              (sequence co-variation)

Experiment recommendation logic:
  intraprotein + contact_class ∈ {electrostatic, hbond}:
    → "thermal shift assay + site-directed mutagenesis"
  interprotein (mt-nuc) + ddg_rescue_bind available:
    → "co-immunoprecipitation + mutagenesis + binding assay"
  interprotein (nuc-nuc or mt-mt):
    → "size-exclusion chromatography + mutagenesis"
  other:
    → "site-directed mutagenesis + functional assay"

Inputs (all optional; script works if only prioritized_pairs.csv is present):
  results/mutagenesis/prioritized_pairs.csv  — required
  results/mutagenesis/foldx_ddg.csv          — optional
  results/mutagenesis/mi_scores.csv          — optional

Outputs:
  results/mutagenesis/final_targets.csv
  results/mutagenesis/mutagenesis_report.md

Usage:
  docker run --rm -v $(pwd):/app oxphos_dav_analysis conda run -n oxphos_dav \\
      python src/mutagenesis/03_compile_targets.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

# ─── Paths ────────────────────────────────────────────────────────────────────
IN_PRIORITY  = ROOT / "results" / "mutagenesis" / "prioritized_pairs.csv"
IN_FOLDX     = ROOT / "results" / "mutagenesis" / "foldx_ddg.csv"
IN_MI        = ROOT / "results" / "mutagenesis" / "mi_scores.csv"
OUT_DIR      = ROOT / "results" / "mutagenesis"

# ─── Join key ─────────────────────────────────────────────────────────────────
JOIN_KEYS = ["ann_id", "dar_gene", "dar_aa_coord", "contact_gene", "contact_refseq_pos"]


def load_optional(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        print(f"  Optional file not found, skipping: {path.name}")
        return None
    df = pd.read_csv(path, dtype={"ann_id": str})
    print(f"  Loaded {path.name}: {len(df)} rows")
    return df


# ─── Experiment recommendation ────────────────────────────────────────────────
def recommend_experiment(row: pd.Series) -> str:
    ct = str(row.get("contact_type", "")).lower()
    cc = str(row.get("contact_class", "")).lower()
    is_intra = str(row.get("is_intraprotein", "True")).lower() == "true"
    ddg_rescue_bind = row.get("ddg_rescue_bind")
    has_bind = pd.notna(ddg_rescue_bind)

    if ct == "mt-nuc":
        if has_bind:
            return "Co-IP pulldown + SDM + binding assay (mCSM-PPI2 validation recommended)"
        return "Co-IP pulldown + SDM"
    if is_intra and cc in ("electrostatic", "hbond"):
        return "Thermal shift assay (DSF/nanoDSF) + SDM"
    if not is_intra:
        return "SEC-MALS + SDM + thermal shift"
    return "SDM + functional assay (CI/III/IV/V activity)"


# ─── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    if not IN_PRIORITY.exists():
        sys.exit(
            f"\nERROR: {IN_PRIORITY} not found.\n"
            "  Run 00_prioritize_pairs.py first."
        )

    print("Loading inputs...")
    priority = pd.read_csv(IN_PRIORITY, dtype={"ann_id": str})
    print(f"  prioritized_pairs.csv: {len(priority)} rows")

    foldx = load_optional(IN_FOLDX)
    mi    = load_optional(IN_MI)

    # Normalise join key types
    for df in [priority, foldx, mi]:
        if df is None:
            continue
        for col in ("dar_aa_coord", "contact_refseq_pos"):
            if col in df.columns:
                df[col] = df[col].astype(str)

    # Merge optional analyses
    final = priority.copy()
    if foldx is not None:
        foldx_cols = JOIN_KEYS + [
            c for c in foldx.columns
            if c not in JOIN_KEYS and c not in final.columns
        ]
        final = final.merge(foldx[foldx_cols], on=JOIN_KEYS, how="left")

    if mi is not None:
        mi_cols = JOIN_KEYS + [
            c for c in mi.columns
            if c not in JOIN_KEYS and c not in final.columns
        ]
        final = final.merge(mi[mi_cols], on=JOIN_KEYS, how="left")

    # ─── Composite score ──────────────────────────────────────────────────────
    composite = final["priority_score"].fillna(0).astype(float)

    # FoldX bonus: +2 if both rescue AND destabilisation exceed 1.0 kcal/mol
    # (1.0 kcal/mol is the defensible threshold given FoldX RMSE of ~1.25 kcal/mol;
    #  0.5 kcal/mol is reported for context but is within FoldX run-to-run noise)
    if "ddg_rescue_stab" in final.columns and "ddg_stab_dar" in final.columns:
        foldx_bonus = (
            (final["ddg_rescue_stab"] > 1.0) & (final["ddg_stab_dar"] > 1.0)
        ).astype(float) * 2
        composite += foldx_bonus

        # FoldX tier classification
        conditions = [
            (final["ddg_rescue_stab"] > 1.0) & (final["ddg_stab_dar"] > 1.0),
            (final["ddg_rescue_stab"] > 0.5) | (final["ddg_stab_dar"] > 0.5),
        ]
        choices = ["strong", "borderline"]
        final["foldx_tier"] = np.select(conditions, choices, default="not_supported")
        final.loc[final["ddg_rescue_stab"].isna(), "foldx_tier"] = "not_computed"
    else:
        final["ddg_rescue_stab"] = np.nan
        final["ddg_stab_dar"]    = np.nan
        final["foldx_tier"]      = "not_computed"

    # Epistasis significance gate: require |epistasis| > 2 × ddg_epistasis_sd
    if "ddg_epistasis" in final.columns and "ddg_epistasis_sd" in final.columns:
        final["ddg_epistasis_sig"] = (
            final["ddg_epistasis_sd"].notna() &
            (final["ddg_epistasis_sd"] > 0) &
            (final["ddg_epistasis"].abs() > 2 * final["ddg_epistasis_sd"])
        )
    else:
        final["ddg_epistasis_sig"] = False

    # MI/DCA bonus: prefer dca_di_percentile if available (phylogenetically corrected),
    # fall back to mi_percentile (confounded, retained for transparency only)
    if "dca_di_percentile" in final.columns:
        # Exclude pairs flagged as unreliable (insufficient M_eff)
        dca_note = final.get("dca_note", pd.Series("", index=final.index)).fillna("")
        mi_bonus = (
            (final["dca_di_percentile"] > 75) &
            (dca_note != "insufficient_meff")
        ).fillna(False).astype(float)
        composite += mi_bonus
    elif "mi_percentile" in final.columns:
        mi_bonus = (final["mi_percentile"] > 75).fillna(False).astype(float)
        composite += mi_bonus
    else:
        final["mi_percentile"] = np.nan

    final["composite_score"] = composite.astype(int)
    final = final.sort_values("composite_score", ascending=False).reset_index(drop=True)
    final["composite_rank"] = range(1, len(final) + 1)

    # ─── Deduplicate on physical positions ────────────────────────────────────
    # Multiple clinical variants (different ann_id) can share the same physical
    # residue pair (same AA change via different nucleotide mutations). For
    # mutagenesis prioritization, the physical pair is what matters — scores
    # are identical for duplicates. Retain count of clinical variants.
    phys_key = [
        "dar_gene", "dar_aa_coord", "dar_alt_aa",
        "contact_gene", "contact_refseq_pos", "contact_alt_aa",
    ]
    phys_key = [c for c in phys_key if c in final.columns]

    # Count clinical variants per physical pair
    n_var = final.groupby(phys_key, sort=False).size().rename("n_clinical_variants").reset_index()
    final = final.drop_duplicates(subset=phys_key, keep="first")
    final = final.merge(n_var, on=phys_key, how="left")

    # Re-rank after dedup
    final = final.sort_values("composite_score", ascending=False).reset_index(drop=True)
    final["composite_rank"] = range(1, len(final) + 1)

    print(f"  After physical deduplication: {len(final)} unique residue pairs")

    # Experiment recommendation
    final["recommended_experiment"] = final.apply(recommend_experiment, axis=1)

    # ─── Save final CSV ───────────────────────────────────────────────────────
    out_csv = OUT_DIR / "final_targets.csv"
    final.to_csv(out_csv, index=False)
    print(f"\nSaved final_targets.csv: {len(final)} unique pairs → {out_csv}")

    # ─── Report ───────────────────────────────────────────────────────────────
    top20 = final.head(20)

    report_cols = [
        "composite_rank", "composite_score", "priority_score",
        "n_clinical_variants",
        "dar_gene", "dar_aa_coord", "dar_ref_aa", "dar_alt_aa",
        "contact_gene", "contact_refseq_pos", "contact_human_aa", "contact_alt_aa",
        "contact_class", "contact_type", "is_intraprotein",
        "physicochemical_type", "dominant_timing",
        "n_dar_gain_branches", "confidence_tier",
        "cbcb_dist_A", "pdb_id",
        "pagel_fdr", "branch_cooccur_fdr",
        "ddg_stab_dar", "ddg_stab_dar_sd", "ddg_rescue_stab", "foldx_tier",
        "ddg_epistasis", "ddg_epistasis_sd", "ddg_epistasis_sig",
        "dca_di", "dca_di_percentile",      # plmDCA (phylogenetically corrected; primary)
        "mi_apc", "mi_percentile",           # APC-MI (phylogenetically confounded; context only)
        "recommended_experiment",
    ]
    report_cols = [c for c in report_cols if c in top20.columns]

    def fmt(col: str, v) -> str:
        if pd.isna(v):
            return "—"
        if col in ("cbcb_dist_A",):
            return f"{float(v):.2f} Å"
        if col in ("pagel_fdr", "branch_cooccur_fdr", "fisher_fdr"):
            return f"{float(v):.2e}"
        if col in ("ddg_stab_dar", "ddg_rescue_stab", "ddg_epistasis",
                   "ddg_bind_dar", "ddg_rescue_bind"):
            return f"{float(v):+.2f}"
        if col in ("ddg_stab_dar_sd", "ddg_epistasis_sd"):
            return f"±{float(v):.2f}"
        if col == "ddg_epistasis_sig":
            return "yes" if v else "no"
        if col == "foldx_tier":
            return str(v)
        if col in ("dca_di",):
            return f"{float(v):.4f}"
        if col in ("dca_di_percentile", "mi_percentile"):
            return f"{float(v):.1f}%"
        if col in ("mi_apc",):
            return f"{float(v):.4f}"
        return str(v)

    has_foldx = "ddg_rescue_stab" in top20.columns and top20["ddg_rescue_stab"].notna().any()
    has_mi    = "mi_apc" in top20.columns and top20["mi_apc"].notna().any()

    lines = [
        "# Mutagenesis Target Report",
        "",
        f"**Top 20 compensatory pairs** ranked by composite score "
        f"(max {17 + 2 + 1} points: 17 priority + 2 FoldX + 1 DCA/MI).",
        "",
        "## Data availability",
        f"- Priority scoring (geometric + phylogenetic): {len(final)} pairs scored",
        f"- FoldX ΔΔG: {'available' if has_foldx else 'not computed (run 01_foldx_ddg.py)'}",
        f"- MI/APC: {'available' if has_mi else 'not computed (run 02_mi_analysis.py)'}",
        "",
        "## Scoring rubric",
        "| Criterion | Points |",
        "|-----------|--------|",
        "| Phylogenetic FDR ≤ 0.01 | 3 |",
        "| Contact class: electrostatic | 3 |",
        "| Physicochemical: charge_reversal | 3 |",
        "| Epistatic pre-adaptation: contact_first | 2 |",
        "| Likely incompatible | 2 |",
        "| Convergence ≥ 10 origins | 2 |",
        "| FoldX ΔΔG_rescue > 1.0 kcal/mol AND ΔΔG_DAR > 1.0 kcal/mol | +2 |",
        "| plmDCA DI > 75th percentile (or MI_APC if DCA not run) | +1 |",
        "| Confidence tier: high | 1 |",
        "| Cβ-Cβ distance < 5 Å | 1 |",
        "",
        "## Top 20 targets",
        "",
        "| " + " | ".join(report_cols) + " |",
        "| " + " | ".join(["---"] * len(report_cols)) + " |",
    ]
    for _, row in top20.iterrows():
        cells = [fmt(c, row.get(c)) for c in report_cols]
        lines.append("| " + " | ".join(cells) + " |")

    lines += [
        "",
        "## Interpretive notes",
        "",
        "**Stability (intraprotein pairs)**:",
        "- foldx_tier='strong': ΔΔG_DAR > 1.0 AND ΔΔG_rescue > 1.0 kcal/mol (reported threshold)",
        "- foldx_tier='borderline': one or both values 0.5–1.0 kcal/mol (within FoldX noise floor; suggestive only)",
        "- ΔΔG_epistasis reported where ddg_epistasis_sig=yes: |epistasis| > 2 × propagated SD",
        "- FoldX RMSE vs experiment: ~1.25 kcal/mol; membrane/cryo-EM structures have higher noise",
        "",
        "**Binding (interprotein pairs)**:",
        "- ΔΔG_rescue_bind > 1.0 kcal/mol → compensatory contact restores subunit-subunit affinity",
        "- For top interprotein pairs, validate with mCSM-PPI2: https://biosig.unimelb.edu.au/mcsm_ppi2/",
        "",
        "**plmDCA / MI**:",
        "- dca_di_percentile > 75 → direct coupling score (plmDCA) in top quartile; phylogenetically corrected",
        "- mi_percentile > 75 → APC-MI in top quartile (phylogenetically confounded; shown for context only)",
        "- High DCA + phylogenetic significance + FoldX rescue (strong tier) = strongest mutagenesis candidate",
        "",
        "**Complex-specific notes**:",
        "- 9I4I (Complex I, 4 Å resolution): FoldX predictions less precise than other structures",
        "- 8GS8 (Complex II, 2.5 Å), 9HZL (CIII), 9I6F (CIV), 8H9S (CV): high-confidence structures",
    ]

    out_md = OUT_DIR / "mutagenesis_report.md"
    out_md.write_text("\n".join(lines) + "\n")
    print(f"Saved mutagenesis_report.md → {out_md}")

    # ─── Console summary ──────────────────────────────────────────────────────
    print("\n── Top 20 pairs (composite score) ───────────────────────────")
    preview = [
        "composite_rank", "composite_score",
        "dar_gene", "dar_aa_coord", "dar_alt_aa",
        "contact_gene", "contact_refseq_pos", "contact_alt_aa",
        "contact_class", "contact_type",
    ]
    preview = [c for c in preview if c in top20.columns]
    print(top20[preview].to_string(index=False))

    print(f"\nComposite score distribution (all {len(final)} pairs):")
    print(final["composite_score"].value_counts().sort_index(ascending=False).head(10).to_string())


if __name__ == "__main__":
    main()
