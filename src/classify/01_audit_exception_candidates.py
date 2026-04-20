#!/usr/bin/env python3
"""
Export focused audit tables for unresolved and warning-bearing classified rows.
"""

from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from utils.exception_registry import load_exception_registry, match_exception_entry


DATA_DIR = ROOT / "data"
DERIVED_CLASSIFIED_DIR = DATA_DIR / "derived" / "classified"
DERIVED_REFERENCE_DIR = DATA_DIR / "derived" / "reference"
REFERENCE_DIR = DATA_DIR / "reference"

CLASSIFIED_MASTER_PARQUET = DERIVED_CLASSIFIED_DIR / "variants_master_classified.parquet"
AUDIT_ROWS_TSV = DERIVED_CLASSIFIED_DIR / "exception_candidate_rows.tsv"
AUDIT_SUMMARY_TSV = DERIVED_CLASSIFIED_DIR / "exception_candidate_summary.tsv"
AUDIT_METADATA_JSON = DERIVED_CLASSIFIED_DIR / "exception_candidate_audit_metadata.json"

EXCEPTION_REGISTRY = DERIVED_REFERENCE_DIR / "variant_exception_registry.tsv"
if not EXCEPTION_REGISTRY.exists():
    EXCEPTION_REGISTRY = REFERENCE_DIR / "variant_exception_registry.tsv"


def ensure_layout() -> None:
    DERIVED_CLASSIFIED_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    import pandas as pd

    ensure_layout()
    df = pd.read_parquet(CLASSIFIED_MASTER_PARQUET)
    registry = load_exception_registry(EXCEPTION_REGISTRY)

    nuc = df[df["genome"] == "nucDNA"].copy()
    unresolved = nuc[nuc["classification_status"] == "unresolved"].copy()
    warning = nuc[
        (nuc["classification_status"] == "classified") &
        (nuc["mismatch_code"] == "REF_ALLELE_MISMATCH")
    ].copy()
    audit = pd.concat([unresolved, warning], ignore_index=True)

    exception_scope = []
    exception_class = []
    exception_code = []
    exception_decision = []
    manual_review_status = []
    for row in audit.to_dict(orient="records"):
        entry = match_exception_entry(
            registry,
            gene=row.get("classification_gene", ""),
            variant_id=row.get("variant_id", ""),
        ) or {}
        exception_scope.append(entry.get("scope", row.get("exception_scope", "")))
        exception_class.append(entry.get("exception_class", row.get("exception_class", "")))
        exception_code.append(entry.get("exception_code", row.get("exception_code", "")))
        exception_decision.append(entry.get("decision", row.get("exception_decision", "")))
        manual_review_status.append(entry.get("manual_review_status", row.get("manual_review_status", "")))

    audit["audit_exception_scope"] = exception_scope
    audit["audit_exception_class"] = exception_class
    audit["audit_exception_code"] = exception_code
    audit["audit_exception_decision"] = exception_decision
    audit["audit_manual_review_status"] = manual_review_status

    keep_cols = [
        "variant_id",
        "source_variant_group_id",
        "classification_gene",
        "aa_change",
        "hgvs_c",
        "hgvs_p",
        "genomic_pos",
        "classification_status",
        "classification_coordinate_method",
        "classification_coordinate_status",
        "mismatch_code",
        "mismatch_reason",
        "ref_allele_match",
        "classification_transcript_map_identity",
        "classification_transcript_map_coverage",
        "exception_scope",
        "exception_class",
        "exception_code",
        "exception_decision",
        "manual_review_status",
        "audit_exception_scope",
        "audit_exception_class",
        "audit_exception_code",
        "audit_exception_decision",
        "audit_manual_review_status",
    ]
    audit[keep_cols].to_csv(AUDIT_ROWS_TSV, sep="\t", index=False)

    summary = (
        audit.groupby(
            [
                "classification_gene",
                "classification_status",
                "mismatch_code",
                "audit_exception_class",
                "audit_exception_decision",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="n_rows")
        .sort_values(["n_rows", "classification_gene"], ascending=[False, True])
    )
    summary.to_csv(AUDIT_SUMMARY_TSV, sep="\t", index=False)

    metadata = {
        "classified_input": str(CLASSIFIED_MASTER_PARQUET.relative_to(ROOT)),
        "exception_registry": str(EXCEPTION_REGISTRY.relative_to(ROOT)) if EXCEPTION_REGISTRY.exists() else None,
        "unresolved_nuc_rows": int(len(unresolved)),
        "classified_warning_rows": int(len(warning)),
        "audit_rows": int(len(audit)),
        "outputs": [
            str(AUDIT_ROWS_TSV.relative_to(ROOT)),
            str(AUDIT_SUMMARY_TSV.relative_to(ROOT)),
        ],
    }
    with open(AUDIT_METADATA_JSON, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    print(f"Wrote exception audit rows: {AUDIT_ROWS_TSV}")
    print(f"Wrote exception audit summary: {AUDIT_SUMMARY_TSV}")
    print(f"Wrote exception audit metadata: {AUDIT_METADATA_JSON}")


if __name__ == "__main__":
    main()
