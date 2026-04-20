#!/usr/bin/env python3
"""
Validate local downloads and probe remote resource availability.

Validation is non-destructive: current dated inputs remain the active analysis
files unless a later step explicitly replaces them.
"""

from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from utils.download_io import remote_probe, safe_validate, validate_gzip_tsv, validate_json, validate_tsv, validate_zip
from utils.download_paths import VALIDATION_REPORT_JSON, ensure_layout, load_manifest, resolve_manifest_path, write_manifest

VALIDATORS = {
    "gzip_tsv": validate_gzip_tsv,
    "json": validate_json,
    "tsv": validate_tsv,
    "zip": validate_zip,
}


def remote_status_suffix(probe: dict) -> str:
    if probe.get("reachable"):
        return ""
    status = probe.get("status")
    if status:
        return f"remote_http_{status}"
    return "remote_unreachable"


def main() -> None:
    ensure_layout()
    records = load_manifest()
    updated = []
    report = []
    for record in records:
        path = resolve_manifest_path(record.local_path)
        validator = VALIDATORS.get(record.file_kind)
        local_status = safe_validate(path, validator) if validator else "unknown_validator"
        if record.resource_name == "myvariant_dbnsfp_gnomad" and local_status == "ok":
            local_status = "ok;completeness_legacy_unchecked"
        probe = remote_probe(record.source_url)
        remote_suffix = remote_status_suffix(probe)
        validation_status = local_status if not remote_suffix else f"{local_status};{remote_suffix}"
        record.validation_status = validation_status
        updated.append(record)
        report.append({
            "resource_name": record.resource_name,
            "local_path": record.local_path,
            "download_date": record.download_date,
            "local_validation": local_status,
            "remote_probe": probe,
            "active_for_pipeline": record.active_for_pipeline,
        })

    write_manifest(updated)
    with open(VALIDATION_REPORT_JSON, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    print(f"Wrote validation report: {VALIDATION_REPORT_JSON}")


if __name__ == "__main__":
    main()
