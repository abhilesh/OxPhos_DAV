import csv
import io
import json
import re
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import pandas as pd


# ==========================================
# 1. CORE DATA MODEL
# ==========================================


@dataclass
class VariantAnnotation:
    """Standardized representation for both mtDNA and nucDNA annotations."""

    # Core parsed attributes
    ann_id: str  # e.g., 'm.8573G>A' or ClinVar AlleleID
    locus: str  # e.g., 'MT-ND5' or 'SDHA'
    nc_change: str  # e.g., 'm.8573G>A' or 'c.327G>C'
    aa_change: str  # e.g., 'L109F'
    is_synonymous: bool  # True if synonymous, False if missense
    disease: str  # Phenotype or disease string
    genome: str  # 'mtDNA' or 'nucDNA'
    reference_assembly: str  # e.g., 'rCRS' or 'GRCh38'
    clinical_status: str  # Original status string from the database

    # Database-specific metadata
    ref_nt: str
    alt_nt: str
    genomic_pos: int
    clinvar_stars: Optional[int] = None
    clinvar_review_status: Optional[str] = None

    # Supplementary metadata
    is_haplogroup_marker: Optional[bool] = None
    mitoclass: Optional[str] = None
    population_frequency: Optional[float] = None

    # Downstream classification attributes
    tier: Optional[str] = None
    pathogenic_score: Optional[float] = None

    def to_dict(self):
        return asdict(self)


# ==========================================
# 2. BASE FILTERS & REFERENCES
# ==========================================


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
        return self.lookup.get(str(symbol).strip())

    def is_target(self, symbol: str) -> bool:
        return str(symbol).strip() in self.lookup


# ==========================================
# 3. MITOCHONDRIAL PARSERS
# ==========================================


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


class PhylotreeParser:
    def __init__(self, zip_path: Path):
        self.markers = self._load(zip_path)

    def _load(self, zip_path: Path) -> dict:
        markers = {}
        with zipfile.ZipFile(zip_path) as zf:
            htm_names = [
                n
                for n in zf.namelist()
                if n.lower().endswith((".htm", ".html")) and "__MACOSX" not in n
            ]
            if not htm_names:
                return markers

            with zf.open(htm_names[0]) as fh:
                content = fh.read().decode("windows-1252", errors="ignore")
                text = re.sub(r"<[^>]+>", " ", content)
                text = re.sub(r"[()]", " ", text)
                text = text.replace("\udca0", " ").replace("\xa0", " ")

                for tok in text.split():
                    tok = tok.strip()
                    if not tok:
                        continue

                    if re.fullmatch(r"\d+\.\d+[a-zA-Z]?", tok.lstrip("@")):
                        continue
                    if re.fullmatch(r"\d+d", tok.lstrip("@")):
                        continue

                    m = re.fullmatch(r"(@?)(\d+)([A-Za-z]?)!?", tok)
                    if m:
                        is_back = m.group(1) == "@"
                        pos = int(m.group(2))
                        explicit_alt = m.group(3).upper() if m.group(3) else None

                        if not is_back:
                            if pos not in markers:
                                markers[pos] = set()
                            markers[pos].add(explicit_alt)
        return markers

    def is_haplogroup(self, pos: int, ref: str, alt: str) -> bool:
        if pos not in self.markers:
            return False
        transition = {"A": "G", "G": "A", "C": "T", "T": "C"}.get(ref.upper())
        for explicit_alt in self.markers[pos]:
            if explicit_alt == alt or (explicit_alt is None and alt == transition):
                return True
        return False


class MitimpactParser:
    def __init__(self, zip_path: Path):
        self.lookup = self._load(zip_path)

    def _load(self, zip_path: Path) -> dict:
        lookup = {}
        with zipfile.ZipFile(zip_path) as zf:
            data_file = next(
                f for f in zf.namelist() if f.endswith((".txt", ".tsv", ".csv"))
            )
            with zf.open(data_file) as fh:
                reader = csv.DictReader(
                    io.TextIOWrapper(fh, encoding="utf-8", errors="ignore"),
                    delimiter="\t",
                )
                for row in reader:
                    try:
                        pos = int(row["Start"])
                        ref = str(row.get("Ref", "")).strip().upper()
                        alt = str(row.get("Alt", "")).strip().upper()
                        key = (pos, ref, alt)

                        score_str = str(row.get("APOGEE2_score", "")).strip()
                        score = float(score_str) if score_str else None

                        mitoclass = str(row.get("Mitoclass1", "")).strip().lower()

                        if key in lookup:
                            existing_score, existing_class = lookup[key]
                            combined_class = (
                                f"{existing_class} | {mitoclass}"
                                if mitoclass and mitoclass not in existing_class
                                else existing_class
                            )
                            lookup[key] = (score, combined_class)
                        else:
                            lookup[key] = (score, mitoclass)
                    except (ValueError, KeyError):
                        continue
        return lookup

    def get_metrics(self, pos: int, ref: str, alt: str) -> tuple:
        return self.lookup.get((pos, ref, alt), (None, ""))


# ==========================================
# 4. NUCLEAR PARSERS
# ==========================================


class ClinvarParser:
    three_to_one = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
    }

    def __init__(self, hgnc_reference=None):
        self.hgnc_reference = hgnc_reference

    def _parse_protein(self, name_field: str):
        m = re.search(
            r"p\.([A-Z][a-z]{0,2}|\*)(\d+)([A-Z][a-z]{0,2}|\*)", str(name_field)
        )
        if not m:
            return None, False

        wt_raw = m.group(1)
        mut_raw = m.group(3)

        stop_flags = ["Ter", "*", "X"]
        if (
            mut_raw in stop_flags
            or wt_raw in stop_flags
            or "Ter" in mut_raw
            or "Ter" in wt_raw
        ):
            return None, False

        wt = self.three_to_one.get(wt_raw, wt_raw)
        mut = self.three_to_one.get(mut_raw, mut_raw)

        return f"{wt}{m.group(2)}{mut}", wt == mut

    def parse(self, file_path: Path) -> list:
        usecols = [
            "Assembly",
            "Chromosome",
            "Start",
            "ReferenceAllele",
            "AlternateAllele",
            "GeneSymbol",
            "ClinicalSignificance",
            "Type",
            "Name",
            "PhenotypeList",
            "ReviewStatus",
            "VariationID",
        ]

        clean_variants = []
        chunk_iterator = pd.read_csv(
            file_path,
            sep="\t",
            compression="gzip",
            usecols=usecols,
            chunksize=100000,
            low_memory=False,
        )

        for chunk in chunk_iterator:
            mask = (
                (chunk["Assembly"] == "GRCh38")
                & (chunk["Type"] == "single nucleotide variant")
                & (chunk["Chromosome"] != "MT")
            )
            filtered_chunk = chunk[mask]

            for _, row in filtered_chunk.iterrows():
                gene = str(row["GeneSymbol"]).strip()

                if self.hgnc_reference and not self.hgnc_reference.is_target(gene):
                    continue

                name = str(row.get("Name", "")).strip()
                aa_change, is_syn = self._parse_protein(name)

                if not aa_change:
                    continue

                nc_match = re.search(r"(c\.[0-9]+[A-Z]>[A-Z])", name)
                nc_change = (
                    nc_match.group(1)
                    if nc_match
                    else f"g.{row['Start']}{row['ReferenceAllele']}>{row['AlternateAllele']}"
                )

                review = str(row.get("ReviewStatus", "")).lower()
                if "practice guideline" in review:
                    stars = 4
                elif "expert panel" in review:
                    stars = 3
                elif "multiple submitters" in review and "no conflicts" in review:
                    stars = 2
                elif "single submitter" in review or "conflicting" in review:
                    stars = 1
                else:
                    stars = 0

                clean_variants.append(
                    {
                        "genome": "nucDNA",
                        "chromosome": str(row["Chromosome"]),
                        "locus": gene,
                        "allele_id": str(row["VariationID"]),
                        "nt_change": nc_change,
                        "aa_change": aa_change,
                        "is_synonymous": is_syn,
                        "disease": str(row.get("PhenotypeList", "")).strip(),
                        "clinical_status": str(
                            row.get("ClinicalSignificance", "")
                        ).strip(),
                        "review_status": str(row.get("ReviewStatus", "")).strip(),
                        "stars": stars,
                        "grch38_pos": row["Start"],
                        "ref": row["ReferenceAllele"],
                        "alt": row["AlternateAllele"],
                    }
                )

        df_clean = pd.DataFrame(clean_variants).drop_duplicates(subset=["allele_id"])
        return df_clean.to_dict("records")


class MyVariantParser:
    def __init__(self, json_path: Path):
        self.lookup = self._load(json_path)

    def _load(self, json_path: Path) -> dict:
        lookup = {}
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)

        for item in data:
            if item.get("notfound"):
                continue

            vid = item.get("_id")
            if not vid:
                continue

            revel = 0.0
            dbnsfp = item.get("dbnsfp", {})
            if dbnsfp:
                revel_score = dbnsfp.get("revel", {}).get("score", 0.0)
                if isinstance(revel_score, list):
                    revel_score = revel_score[0]
                try:
                    revel = float(revel_score)
                except (ValueError, TypeError):
                    pass

            af = 0.0
            gnomad = item.get("gnomad_exome", item.get("gnomad_genome", {}))
            if gnomad:
                af_val = gnomad.get("af", 0.0)
                if isinstance(af_val, list):
                    af_val = af_val[0]
                try:
                    af = float(af_val)
                except (ValueError, TypeError):
                    pass

            lookup[vid] = {"revel": revel, "gnomad_af": af}
        return lookup

    def get_metrics(self, chrom: str, pos: int, ref: str, alt: str) -> tuple:
        vid = f"chr{chrom}:g.{pos}{ref}>{alt}"
        metrics = self.lookup.get(vid, {})
        return metrics.get("revel", 0.0), metrics.get("gnomad_af", 0.0)
