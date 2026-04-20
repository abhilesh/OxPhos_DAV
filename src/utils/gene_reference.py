from __future__ import annotations

import csv
from pathlib import Path


class GeneReference:
    def __init__(self, hgnc_csv: Path):
        self.lookup: dict[str, dict] = {}
        self.mt_genes: list[tuple[str, int, int]] = []
        self._load_hgnc(hgnc_csv)

    def _load_hgnc(self, hgnc_csv: Path) -> None:
        with open(hgnc_csv, "r", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                symbol = row["Approved symbol"].strip()
                entry = {
                    "name": row.get("Approved name", ""),
                    "group": row.get("Group name", ""),
                    "symbol": symbol,
                    "primary_symbol": symbol,
                }
                self.lookup[symbol] = entry
                for prev in str(row.get("Previous symbols", "")).split(","):
                    prev = prev.strip()
                    if prev:
                        self.lookup[prev] = entry

    def get_gene_data(self, symbol: str) -> dict:
        return self.lookup.get(symbol, {})

    def load_coordinates(self, coord_tsv: Path) -> None:
        if not coord_tsv.exists():
            return
        with open(coord_tsv, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                symbol = row["gene"]
                if symbol not in self.lookup:
                    continue
                raw_start = int(row["start"])
                raw_end = int(row["end"])
                coord_start = min(raw_start, raw_end)
                coord_end = max(raw_start, raw_end)
                self.lookup[symbol].update(
                    {
                        "chr": row.get("chr", ""),
                        "start": coord_start,
                        "end": coord_end,
                        "strand": row["strand"],
                    }
                )
                if row.get("chr") in ("MT", "chrM", "") or symbol.startswith("MT-"):
                    self.mt_genes.append((symbol, coord_start, coord_end))

    def get_mt_locus_by_position(self, pos: int) -> str:
        genes = [gene for gene, start, end in self.mt_genes if start <= pos <= end]
        return "/".join(genes) if genes else "Non-OXPHOS"
