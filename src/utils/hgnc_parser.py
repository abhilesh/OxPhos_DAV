import csv
from pathlib import Path


class GeneReference:
    def __init__(self, file_path: Path):
        self.lookup = self._load_hgnc_reference(file_path)

    def _load_hgnc_reference(self, file_path: Path) -> dict:
        """Loads HGNC data into a dictionary for O(1) retrieval."""
        gene_reference = {}

        with open(file_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")

            for row in reader:
                primary_symbol = row.get("Approved symbol", "").strip()
                if not primary_symbol:
                    continue

                gene_data = {
                    "hgnc_id": row.get("HGNC ID", ""),
                    "primary_symbol": primary_symbol,
                    "group": row.get("Group name", "Unknown OXPHOS Complex"),
                    "ensembl_id": row.get("Ensembl gene ID", ""),
                }

                gene_reference[primary_symbol] = gene_data

                aliases = row.get("Alias symbols", "").split(",")
                for alias in aliases:
                    alias = alias.strip()
                    if alias and alias not in gene_reference:
                        gene_reference[alias] = gene_data

                prev_symbols = row.get("Previous symbols", "").split(",")
                for prev in prev_symbols:
                    prev = prev.strip()
                    if prev and prev not in gene_reference:
                        gene_reference[prev] = gene_data

        return gene_reference

    def get_gene_data(self, symbol: str) -> dict:
        """Returns the gene data dictionary if found, otherwise None."""
        return self.lookup.get(str(symbol).strip())

    def is_target(self, symbol: str) -> bool:
        """Quick boolean check if a gene symbol exists in the reference."""
        return str(symbol).strip() in self.lookup
