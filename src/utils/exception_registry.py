from pathlib import Path
import csv


def load_exception_registry(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def normalize_registry_entry(entry: dict) -> dict:
    normalized = {k: (v.strip() if isinstance(v, str) else v) for k, v in entry.items()}
    return normalized


def match_exception_entry(
    entries: list[dict],
    *,
    gene: str = "",
    variant_id: str = "",
) -> dict | None:
    gene = str(gene or "").strip()
    variant_id = str(variant_id or "").strip()

    normalized = [normalize_registry_entry(entry) for entry in entries]

    for entry in normalized:
        if entry.get("scope") == "variant" and variant_id and entry.get("variant_id") == variant_id:
            return entry

    for entry in normalized:
        if entry.get("scope") == "gene" and gene and entry.get("gene") == gene:
            return entry

    return None
