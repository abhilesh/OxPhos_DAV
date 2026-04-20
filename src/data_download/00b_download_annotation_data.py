#!/usr/bin/env python3
"""
Download and register disease-associated variant annotation resources.

Primary sources:
- MITOMAP: https://www.mitomap.org/MITOMAP/MutationsCoding
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/docs/data_file_download/
- MitImpact: https://mitimpact.css-mendel.it/
- PhyloTree: https://www.phylotree.org/
- MANE: https://www.ncbi.nlm.nih.gov/refseq/MANE/
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from utils.download_io import (
    fetch_stream_to_file,
    safe_validate,
    validate_gzip_tsv,
    validate_json,
    validate_tsv,
    validate_zip,
)
from utils.download_paths import (
    LEGACY_ANNOTATIONS_RAW_DIR,
    LEGACY_REFERENCE_DIR,
    ensure_layout,
    latest_existing,
    load_manifest,
    record_file,
    sync_compat_file,
    upsert_record,
    write_manifest,
)
from utils.download_resources import RESOURCE_CONFIG

VALIDATORS = {
    "gzip_tsv": validate_gzip_tsv,
    "json": validate_json,
    "tsv": validate_tsv,
    "zip": validate_zip,
}


def download_if_missing(resource_name: str) -> Path:
    cfg = RESOURCE_CONFIG[resource_name]
    patterns = [(cfg["target"].parent, cfg["pattern"])]
    patterns.extend((cfg["target"].parent, pattern) for pattern in cfg.get("fallback_patterns", []))
    existing = latest_existing(patterns)
    target = existing or cfg["target"]
    if target.exists():
        return target

    if resource_name == "myvariant_dbnsfp_gnomad":
        raise RuntimeError(
            "MyVariant acquisition is no longer a blind gene-batch fallback. "
            "Use the existing dated file for this phase or implement curated variant identity queries later."
        )

    print(f"Downloading {resource_name} to {target}...")
    fetch_stream_to_file(cfg["url"], target)
    return target


def register_download(resource_name: str, target: Path) -> None:
    cfg = RESOURCE_CONFIG[resource_name]
    validator = VALIDATORS[cfg["file_kind"]]
    validation = safe_validate(target, validator)
    records = upsert_record(
        load_manifest(),
        record_file(
            resource_name,
            cfg["url"],
            target,
            cfg["file_kind"],
            version_hint=cfg.get("version_hint", ""),
            active_for_pipeline=True,
            validation_status=validation,
        ),
    )
    write_manifest(records)

    if target.parent == LEGACY_ANNOTATIONS_RAW_DIR or target.parent == LEGACY_REFERENCE_DIR:
        return
    if "annotations" in str(target.parent):
        sync_compat_file(target, LEGACY_ANNOTATIONS_RAW_DIR / target.name)
    else:
        sync_compat_file(target, LEGACY_REFERENCE_DIR / target.name)


def main() -> None:
    ensure_layout()
    ordered = [
        "clinvar_variant_summary",
        "mitomap_coding_variants",
        "mitimpact_db",
        "phylotree_build17",
        "mane_grch38_summary",
        "toga_overview_hg38",
        "myvariant_dbnsfp_gnomad",
    ]
    for resource_name in ordered:
        cfg = RESOURCE_CONFIG[resource_name]
        patterns = [(cfg["target"].parent, cfg["pattern"])]
        patterns.extend((cfg["target"].parent, pattern) for pattern in cfg.get("fallback_patterns", []))
        existing = latest_existing(patterns)
        if existing is None and resource_name == "myvariant_dbnsfp_gnomad":
            print("Skipping fresh MyVariant fetch; an existing dated file is required until identity-based retrieval is implemented.")
            continue
        target = existing or download_if_missing(resource_name)
        register_download(resource_name, target)
        print(f"Registered {resource_name}: {target}")


if __name__ == "__main__":
    main()
