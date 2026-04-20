from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Optional

TODAY = date.today().isoformat()
DATE_RE = re.compile(r"(\d{4}-\d{2}-\d{2})")

ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
RAW_ANNOTATIONS_DIR = RAW_DIR / "annotations"
RAW_REFERENCE_DIR = RAW_DIR / "reference"
DERIVED_DIR = DATA_DIR / "derived"
DERIVED_REFERENCE_DIR = DERIVED_DIR / "reference"

# Compatibility paths kept alive during transition.
LEGACY_ANNOTATIONS_RAW_DIR = DATA_DIR / "annotations" / "raw"
LEGACY_REFERENCE_DIR = DATA_DIR / "reference"

DOWNLOAD_MANIFEST_JSON = DERIVED_REFERENCE_DIR / "download_manifest.json"
DOWNLOAD_MANIFEST_CSV = DERIVED_REFERENCE_DIR / "download_manifest.csv"
DOWNLOAD_MANIFEST_PARQUET = DERIVED_REFERENCE_DIR / "download_manifest.parquet"
VALIDATION_REPORT_JSON = DERIVED_REFERENCE_DIR / "download_validation_report.json"


@dataclass
class DownloadRecord:
    resource_name: str
    source_url: str
    download_date: str
    local_path: str
    file_kind: str
    version_hint: str
    checksum: str
    active_for_pipeline: bool
    validation_status: str


def ensure_layout() -> None:
    for path in (
        RAW_ANNOTATIONS_DIR,
        RAW_REFERENCE_DIR,
        DERIVED_REFERENCE_DIR,
        LEGACY_ANNOTATIONS_RAW_DIR,
        LEGACY_REFERENCE_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def parse_date_from_name(path: Path) -> str:
    match = DATE_RE.search(path.name)
    return match.group(1) if match else ""


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sorted_existing(pattern: str, directory: Path) -> list[Path]:
    return sorted(directory.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)


def latest_existing(patterns: Iterable[tuple[Path, str]]) -> Optional[Path]:
    matches: list[Path] = []
    for directory, pattern in patterns:
        matches.extend(sorted_existing(pattern, directory))
    return matches[0] if matches else None


def write_manifest(records: list[DownloadRecord]) -> None:
    ensure_layout()
    by_resource: dict[str, list[DownloadRecord]] = {}
    for record in records:
        by_resource.setdefault(record.resource_name, []).append(record)

    normalized: list[DownloadRecord] = []
    for resource_records in by_resource.values():
        resource_records.sort(key=lambda rec: (bool(rec.download_date), rec.download_date, rec.local_path))
        active_idx = len(resource_records) - 1
        for idx, record in enumerate(resource_records):
            record.active_for_pipeline = idx == active_idx
            normalized.append(record)

    normalized.sort(key=lambda rec: (rec.resource_name, rec.download_date, rec.local_path))
    data = [asdict(r) for r in normalized]
    with open(DOWNLOAD_MANIFEST_JSON, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)

    fieldnames = list(DownloadRecord.__dataclass_fields__.keys())
    with open(DOWNLOAD_MANIFEST_CSV, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

    try:
        import pandas as pd

        pd.DataFrame(data).to_parquet(DOWNLOAD_MANIFEST_PARQUET, index=False)
    except Exception:
        # Parquet is produced when the runtime includes pyarrow/fastparquet.
        pass


def load_manifest() -> list[DownloadRecord]:
    if not DOWNLOAD_MANIFEST_JSON.exists():
        return []
    with open(DOWNLOAD_MANIFEST_JSON, encoding="utf-8") as handle:
        raw = json.load(handle)
    return [DownloadRecord(**row) for row in raw]


def upsert_record(records: list[DownloadRecord], new_record: DownloadRecord) -> list[DownloadRecord]:
    kept: list[DownloadRecord] = []
    for rec in records:
        if rec.resource_name == new_record.resource_name and rec.local_path == new_record.local_path:
            continue
        if new_record.active_for_pipeline and rec.resource_name == new_record.resource_name:
            rec.active_for_pipeline = False
        kept.append(rec)
    kept.append(new_record)
    kept.sort(key=lambda rec: (rec.resource_name, rec.download_date, rec.local_path))
    return kept


def record_file(
    resource_name: str,
    source_url: str,
    path: Path,
    file_kind: str,
    version_hint: str = "",
    active_for_pipeline: bool = True,
    validation_status: str = "not_validated",
) -> DownloadRecord:
    return DownloadRecord(
        resource_name=resource_name,
        source_url=source_url,
        download_date=parse_date_from_name(path),
        local_path=str(path.relative_to(ROOT)),
        file_kind=file_kind,
        version_hint=version_hint,
        checksum=sha256_file(path),
        active_for_pipeline=active_for_pipeline,
        validation_status=validation_status,
    )


def resolve_manifest_path(local_path: str) -> Path:
    return ROOT / local_path


def sync_compat_file(src: Path, legacy_dest: Path) -> None:
    legacy_dest.parent.mkdir(parents=True, exist_ok=True)
    if legacy_dest.resolve() == src.resolve():
        return
    shutil.copy2(src, legacy_dest)
