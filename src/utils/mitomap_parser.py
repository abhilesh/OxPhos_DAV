import re
import pandas as pd
from pathlib import Path


class MitomapParser:
    def __init__(self, hgnc_reference=None):
        self.hgnc_reference = hgnc_reference
        self.mt_genes = [
            ("MT-ND1", 3307, 4262),
            ("MT-ND2", 4470, 5511),
            ("MT-CO1", 5904, 7445),
            ("MT-CO2", 7586, 8269),
            ("MT-ATP8", 8366, 8572),
            ("MT-ATP6", 8527, 9207),
            ("MT-CO3", 9207, 9990),
            ("MT-ND3", 10059, 10404),
            ("MT-ND4L", 10470, 10766),
            ("MT-ND4", 10760, 12137),
            ("MT-ND5", 12337, 14148),
            ("MT-ND6", 14149, 14673),
            ("MT-CYB", 14747, 15887),
        ]

    def _get_locus(self, pos: int) -> str:
        genes = [name for name, start, end in self.mt_genes if start <= pos <= end]
        return "/".join(genes) if genes else "Non-OXPHOS"

    def parse(self, file_path: Path) -> list:
        # Using windows-1252 to handle legacy bytes like 0xa0
        df = pd.read_csv(
            file_path, sep="\t", on_bad_lines="skip", encoding="windows-1252"
        )
        df.columns = [c.strip().lower() for c in df.columns]

        clean_variants = []

        for _, row in df.iterrows():
            try:
                pos = int(row.get("pos", 0))
            except ValueError:
                continue

            ref = str(row.get("ref", "")).strip().upper()
            alt = str(row.get("alt", "")).strip().upper()
            aachange = str(row.get("aachange", "")).strip()

            locus = self._get_locus(pos)

            # Exclude non-OXPHOS variants and those not in HGNC reference (if provided)
            if locus == "Non-OXPHOS":
                continue

            if self.hgnc_reference:
                locus_genes = locus.split("/")
                if not any(self.hgnc_reference.is_target(g) for g in locus_genes):
                    continue

            # Strict parsing for point mutations and valid amino acid changes
            if len(ref) > 1 or len(alt) > 1:
                continue
            if (
                not aachange
                or "noncoding" in aachange.lower()
                or "frameshift" in aachange.lower()
            ):
                continue
            if "*" in aachange or "Ter" in aachange or "Stop" in aachange:
                continue

            match = re.match(r"^([a-zA-Z]+)(\d+)([a-zA-Z]+)$", aachange)
            if not match:
                continue

            is_synonymous = match.group(1).upper() == match.group(3).upper()

            clean_variants.append(
                {
                    "genome": "mtDNA",
                    "locus": locus,
                    "nt_change": f"m.{pos}{ref}>{alt}",
                    "aa_change": aachange,
                    "is_synonymous": is_synonymous,
                    "disease": str(row.get("disease", "")).strip(),
                    "clinical_status": str(row.get("status", "")).strip(),
                    "rCRS_pos": pos,
                    "ref": ref,
                    "alt": alt,
                }
            )

        return clean_variants
