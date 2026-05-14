#!/usr/bin/env python3
"""
Verify frozen pipeline artifacts against data/reference/pipeline_audit_registry.tsv.

This checker is read-only: it hashes existing files and compares byte sizes and
row counts recorded at audit time. Parquet row-count validation is enabled when
pyarrow is available; hashes are always checked.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = ROOT / "data" / "reference" / "pipeline_audit_registry.tsv"


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def line_count_minus_header(path: Path) -> int:
    with open(path, "rb") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def jsonl_count(path: Path) -> int:
    with open(path, "rb") as handle:
        return sum(1 for _ in handle)


def parquet_count(path: Path) -> int | None:
    try:
        import pyarrow.parquet as pq
    except ImportError:
        return None
    return pq.ParquetFile(path).metadata.num_rows


def observed_row_count(path: Path) -> int | None:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return parquet_count(path)
    if suffix == ".jsonl":
        return jsonl_count(path)
    if suffix in {".csv", ".tsv"}:
        return line_count_minus_header(path)
    return None


def load_registry(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check_registry(registry_path: Path, stages: set[str] | None, strict_row_counts: bool) -> int:
    rows = load_registry(registry_path)
    failures: list[str] = []
    warnings: list[str] = []
    checked = 0

    for row in rows:
        if row.get("frozen", "").lower() != "true":
            continue
        if stages and row.get("stage_id") not in stages:
            continue

        checked += 1
        artifact = ROOT / row["artifact_path"]
        label = f"{row['stage_id']} {row['artifact_path']}"

        if not artifact.exists():
            failures.append(f"{label}: missing artifact")
            continue

        expected_size = int(row["byte_size"])
        actual_size = artifact.stat().st_size
        if actual_size != expected_size:
            failures.append(f"{label}: byte_size {actual_size} != {expected_size}")

        actual_sha = sha256_file(artifact)
        if actual_sha != row["sha256"]:
            failures.append(f"{label}: sha256 {actual_sha} != {row['sha256']}")

        expected_rows = row.get("row_count", "NA")
        if expected_rows not in {"", "NA"}:
            actual_rows = observed_row_count(artifact)
            if actual_rows is None:
                msg = f"{label}: row_count not checked; install pyarrow for parquet row counts"
                if strict_row_counts:
                    failures.append(msg)
                else:
                    warnings.append(msg)
            elif actual_rows != int(expected_rows):
                failures.append(f"{label}: row_count {actual_rows} != {expected_rows}")

    print(f"Checked frozen artifacts: {checked}")
    for warning in warnings:
        print(f"WARNING: {warning}")
    if failures:
        print("Freeze check failed:")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print("Freeze check passed.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--registry",
        type=Path,
        default=DEFAULT_REGISTRY,
        help="Path to pipeline audit registry TSV.",
    )
    parser.add_argument(
        "--stage",
        action="append",
        choices=["0_download", "1_curation", "2_classification"],
        help="Restrict verification to one stage. Can be passed multiple times.",
    )
    parser.add_argument(
        "--strict-row-counts",
        action="store_true",
        help="Fail if row counts cannot be checked because optional dependencies are missing.",
    )
    args = parser.parse_args()
    stages = set(args.stage) if args.stage else None
    return check_registry(args.registry, stages, args.strict_row_counts)


if __name__ == "__main__":
    sys.exit(main())
