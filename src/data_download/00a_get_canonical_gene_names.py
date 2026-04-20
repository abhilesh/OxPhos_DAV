#!/usr/bin/env python3
"""
Download canonical gene names for OXPHOS genes from HGNC.
Source: https://www.genenames.org/data/genegroup/#!/group/639
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from utils.download_io import fetch_bytes, safe_validate, validate_tsv
from utils.download_paths import (
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


def main() -> None:
    ensure_layout()
    cfg = RESOURCE_CONFIG["hgnc_oxphos_gene_list"]
    existing = latest_existing([(cfg["target"].parent, cfg["pattern"])])
    target = existing or cfg["target"]

    if not target.exists():
        print(f"Downloading HGNC OXPHOS gene list to {target}...")
        target.write_bytes(fetch_bytes(cfg["url"]))

    validation = safe_validate(target, validate_tsv)
    manifest = upsert_record(
        load_manifest(),
        record_file(
            "hgnc_oxphos_gene_list",
            cfg["url"],
            target,
            cfg["file_kind"],
            active_for_pipeline=True,
            validation_status=validation,
        ),
    )
    write_manifest(manifest)
    sync_compat_file(target, LEGACY_REFERENCE_DIR / target.name)
    print(f"Using HGNC OXPHOS gene list: {target}")


if __name__ == "__main__":
    main()
