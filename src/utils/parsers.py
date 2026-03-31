import csv
import gzip
import io
import json
import re
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
import pandas as pd


@dataclass
class VariantAnnotation:
    ann_id: str
    locus: str
    nc_change: str
    aa_change: str
    is_synonymous: bool
    disease: str
    genome: str
    reference_assembly: str
    clinical_status: str
    ref_nt: str
    alt_nt: str
    ref_aa: str = ""
    alt_aa: str = ""
    genomic_pos: int = 0
    clinvar_stars: int = 0
    clinvar_review_status: str = ""
    pathogenic_score: float = 0.0
    population_frequency: float = 0.0
    mitoclass: str = ""
    is_haplogroup_marker: bool = False
    tier: str = "Unassigned"

    def to_dict(self):
        return asdict(self)


class GeneReference:
    def __init__(self, hgnc_csv: Path):
        self.lookup = {}
        self.mt_genes = []
        self._load_hgnc(hgnc_csv)

    def _load_hgnc(self, hgnc_csv: Path):
        with open(hgnc_csv, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sym = row["Approved symbol"].strip()
                entry = {
                    "name": row["Approved name"],
                    "group": row["Group name"],
                    "symbol": sym,
                }
                self.lookup[sym] = entry
                # Also index previous symbols so ClinVar variants annotated
                # with old gene names (e.g. NDUFA4 -> COXFA4) are accepted
                for prev in str(row.get("Previous symbols", "")).split(","):
                    prev = prev.strip()
                    if prev:
                        self.lookup[prev] = entry

    def load_coordinates(self, coord_tsv: Path):
        if not coord_tsv.exists():
            return
        with open(coord_tsv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                symbol = row["gene"]
                if symbol in self.lookup:
                    raw_start, raw_end = int(row["start"]), int(row["end"])
                    # Minus-strand genes are stored as start > end in the coordinate
                    # file; normalise to (min, max) for range checks
                    coord_start = min(raw_start, raw_end)
                    coord_end = max(raw_start, raw_end)
                    self.lookup[symbol].update(
                        {
                            "chr": row["chr"],
                            "start": coord_start,
                            "end": coord_end,
                            "strand": row["strand"],
                        }
                    )
                    if row["chr"] in ["MT", "chrM"] or symbol.startswith("MT-"):
                        self.mt_genes.append(
                            (symbol, coord_start, coord_end)
                        )

    def get_mt_locus_by_position(self, pos: int) -> str:
        genes = [gene for gene, start, end in self.mt_genes if start <= pos <= end]
        return "/".join(genes) if genes else "Non-OXPHOS"

    def get_gene_data(self, symbol: str) -> dict:
        """Retrieves coordinate and strand data for a specific gene locus."""
        return self.lookup.get(symbol, {})


class MitomapParser:
    def __init__(self, hgnc_reference: GeneReference):
        self.hgnc_reference = hgnc_reference

    def parse(self, file_path: Path) -> list:
        df = pd.read_csv(
            file_path, sep="\t", on_bad_lines="skip", encoding="windows-1252"
        )
        df.columns = [c.strip().lower() for c in df.columns]
        clean_variants = []
        for _, row in df.iterrows():
            try:
                pos = int(row.get("pos", 0))
            except:
                continue
            locus = self.hgnc_reference.get_mt_locus_by_position(pos)
            if locus == "Non-OXPHOS":
                continue

            ref, alt = str(row.get("ref", "")).upper(), str(row.get("alt", "")).upper()
            aachange = str(row.get("aachange", "")).strip()
            if not aachange or any(
                x in aachange.lower() for x in ["noncoding", "frameshift", "*", "ter"]
            ):
                continue

            match = re.match(r"^([a-zA-Z]+)(\d+)([a-zA-Z]+)$", aachange)
            if not match:
                continue

            # Strand correction for MT-ND6
            if 14149 <= pos <= 14673:
                comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
                c_ref, c_alt = comp.get(ref, ref), comp.get(alt, alt)
            else:
                c_ref, c_alt = ref, alt

            clean_variants.append(
                {
                    "genome": "mtDNA",
                    "locus": locus,
                    "nt_change": f"m.{pos}{ref}>{alt}",
                    "aa_change": aachange,
                    "ref_aa": match.group(1).upper(),
                    "alt_aa": match.group(3).upper(),
                    "is_synonymous": match.group(1).upper() == match.group(3).upper(),
                    "disease": str(row.get("disease", "")).strip(),
                    "clinical_status": str(row.get("status", "")).strip(),
                    "rCRS_pos": pos,
                    "ref": c_ref,
                    "alt": c_alt,
                    "genomic_ref": ref,
                    "genomic_alt": alt,
                }
            )
        return clean_variants


class ClinvarParser:
    def __init__(self, hgnc_reference=None):
        self.hgnc_reference = hgnc_reference
        self.three_to_one = {
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

    def _parse_protein(self, name):
        m = re.search(r"p\.([A-Z][a-z]{2}|\*)(\d+)([A-Z][a-z]{2}|\*)", str(name))
        if not m or "*" in m.groups():
            return None, False
        wt, mut = self.three_to_one.get(m.group(1), "X"), self.three_to_one.get(
            m.group(3), "X"
        )
        return f"{wt}{m.group(2)}{mut}", wt == mut

    def parse(self, file_path: Path) -> list:
        clean = []
        for chunk in pd.read_csv(
            file_path, sep="\t", compression="gzip", chunksize=50000, low_memory=False
        ):
            mask = (
                (chunk["Assembly"] == "GRCh38")
                & (chunk["Type"] == "single nucleotide variant")
                & (chunk["Chromosome"] != "MT")
            )
            for _, row in chunk[mask].iterrows():
                gene = str(row["GeneSymbol"]).strip()
                if self.hgnc_reference:
                    if gene not in self.hgnc_reference.lookup:
                        continue
                    # Resolve to canonical approved symbol (handles renamed genes)
                    gene = self.hgnc_reference.lookup[gene].get("symbol", gene)
                aa_change, is_syn = self._parse_protein(row["Name"])
                if not aa_change or "X" in aa_change:
                    continue
                _aa_m = re.match(r"^([A-Z])(\d+)([A-Z])$", aa_change)
                ref_aa = _aa_m.group(1) if _aa_m else ""
                alt_aa = _aa_m.group(3) if _aa_m else ""

                review = str(row.get("ReviewStatus", "")).strip()
                r = review.lower()
                if "guideline" in r:
                    stars = 4
                elif "expert" in r:
                    stars = 3
                elif "multiple" in r and "no conflicts" in r:
                    stars = 2
                elif "single submitter" in r or "conflicting" in r:
                    stars = 1
                else:
                    stars = 0

                clean.append(
                    {
                        "genome": "nucDNA",
                        "chromosome": str(row["Chromosome"]),
                        "locus": gene,
                        "allele_id": str(row["VariationID"]),
                        "nt_change": str(row["Name"]),
                        "aa_change": aa_change,
                        "ref_aa": ref_aa,
                        "alt_aa": alt_aa,
                        "is_synonymous": is_syn,
                        "disease": str(row.get("PhenotypeList", "")),
                        "clinical_status": str(row.get("ClinicalSignificance", "")),
                        "review_status": review,
                        "stars": stars,
                        "grch38_pos": int(row["Start"]),
                        "ref": str(row["ReferenceAlleleVCF"]),
                        "alt": str(row["AlternateAlleleVCF"]),
                    }
                )
        return clean


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


class MyVariantParser:
    def __init__(self, json_path: Path):
        with open(json_path, "r") as f:
            self.lookup = {i["_id"]: i for i in json.load(f) if not i.get("notfound")}

    def _extract_float(self, value):
        if isinstance(value, list):
            return self._extract_float(value[0]) if value else 0.0
        if isinstance(value, dict):
            for k in ["score", "af", "value"]:
                if k in value:
                    return self._extract_float(value[k])
            return 0.0
        try:
            return float(value)
        except:
            return 0.0

    def get_metrics(self, chrom, pos, ref, alt):
        item = self.lookup.get(f"chr{chrom}:g.{pos}{ref}>{alt}", {})
        revel = self._extract_float(
            item.get("dbnsfp", {}).get("revel", {}).get("score", 0.0)
        )
        af = self._extract_float(
            item.get("gnomad_exome", item.get("gnomad_genome", {})).get("af", 0.0)
        )
        return revel, af
