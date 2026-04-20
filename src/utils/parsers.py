"""
src/utils/parsers.py

All parsers used by the OXPHOS curation and classification pipeline.

  GeneReference     — HGNC gene metadata, coordinate lookup, MT-locus resolution
  MitomapParser     — MITOMAP CodingVariants TSV → rich evidence dicts
  ClinvarParser     — ClinVar variant_summary.txt.gz → rich evidence dicts
  MyVariantParser   — pre-downloaded MyVariant JSON → full metrics dict
  PhylotreeParser   — PhyloTree ZIP → haplogroup marker lookup
  MitimpactParser   — MitImpact ZIP → APOGEE2 / Mitoclass lookup

Each parser returns plain dicts; 01_curate_variants.py converts them to
VariantRecord objects.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
import zipfile
from pathlib import Path
from typing import Optional

import pandas as pd


# ── Helper: scalar extraction from nested MyVariant JSON ─────────────────────

def _first(val, default=None):
    """Return first element if list, val itself otherwise, or default if falsy."""
    if val is None:
        return default
    if isinstance(val, list):
        return val[0] if val else default
    return val


def _float(val, default=None) -> Optional[float]:
    try:
        return float(_first(val, default))
    except (TypeError, ValueError):
        return default


def _int(val, default=None) -> Optional[int]:
    try:
        v = _first(val, default)
        return int(v) if v is not None else default
    except (TypeError, ValueError):
        return default


# ── GeneReference ─────────────────────────────────────────────────────────────

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
                    "primary_symbol": sym,
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


# ── MitomapParser ─────────────────────────────────────────────────────────────

class MitomapParser:
    """
    Parse MITOMAP CodingVariants TSV.

    New fields vs. original:
      - source_db_version: extracted from file name (YYYY-MM-DD)
      - source_record_id: MITOMAP internal id
      - mitomap_plasmy: homo / hetero / homo/hetero
      - mitomap_pubmed_count: number of PubMed IDs
    """

    def __init__(self, hgnc_reference):
        self.hgnc_reference = hgnc_reference

    @staticmethod
    def _extract_version(file_path: Path) -> str:
        m = re.search(r"(\d{4}-\d{2}-\d{2})", file_path.name)
        return m.group(1) if m else ""

    @staticmethod
    def _plasmy(homo: str, hetero: str) -> str:
        h1 = str(homo).strip().lower() not in ("", "no", "nr", "nan")
        h2 = str(hetero).strip().lower() not in ("", "no", "nr", "nan")
        if h1 and h2:
            return "homo/hetero"
        if h1:
            return "homo"
        if h2:
            return "hetero"
        return ""

    @staticmethod
    def _pubmed_count(pubmed_field: str) -> int:
        raw = str(pubmed_field).strip()
        if not raw or raw.lower() in ("nan", ""):
            return 0
        return len([x for x in raw.split(",") if x.strip()])

    def parse(self, file_path: Path) -> list[dict]:
        version = self._extract_version(file_path)
        df = pd.read_csv(
            file_path, sep="\t", on_bad_lines="skip", encoding="windows-1252"
        )
        df.columns = [c.strip().lower() for c in df.columns]

        clean = []
        for _, row in df.iterrows():
            try:
                pos = int(row.get("pos", 0))
            except (ValueError, TypeError):
                continue

            locus = self.hgnc_reference.get_mt_locus_by_position(pos)
            if locus == "Non-OXPHOS":
                continue

            ref = str(row.get("ref", "")).upper().strip()
            alt = str(row.get("alt", "")).upper().strip()
            # Skip multi-nucleotide variants (e.g. m.8993TG>CA) — single-base
            # injection in check_compensation cannot handle simultaneous changes
            if len(ref) != 1 or len(alt) != 1:
                continue
            aachange = str(row.get("aachange", "")).strip()
            if not aachange or any(
                x in aachange.lower() for x in ["noncoding", "frameshift", "*", "ter"]
            ):
                continue

            match = re.match(r"^([a-zA-Z]+)(\d+)([a-zA-Z]+)$", aachange)
            if not match:
                continue

            # Strand correction for MT-ND6 (minus-strand, rCRS 14149-14673)
            comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
            if 14149 <= pos <= 14673:
                c_ref, c_alt = comp.get(ref, ref), comp.get(alt, alt)
            else:
                c_ref, c_alt = ref, alt

            disease = str(row.get("disease", "")).strip()
            status  = str(row.get("status", "")).strip()

            clean.append({
                "source_db":          "MITOMAP",
                "source_db_version":  version,
                "source_record_id":   str(row.get("id", "")),
                "genome":             "mtDNA",
                "locus":              locus,
                "nt_change":          f"m.{pos}{ref}>{alt}",
                "hgvs_c":             f"m.{pos}{ref}>{alt}",
                "hgvs_p":             "",
                "aa_change":          aachange,
                "ref_aa":             match.group(1).upper(),
                "alt_aa":             match.group(3).upper(),
                "is_synonymous":      match.group(1).upper() == match.group(3).upper(),
                "disease":            disease,
                "clinical_status":    status,
                "rCRS_pos":           pos,
                "ref":                c_ref,
                "alt":                c_alt,
                "genomic_ref":        ref,
                "genomic_alt":        alt,
                "mitomap_plasmy":     self._plasmy(
                    row.get("homoplasmy", ""),
                    row.get("heteroplasmy", ""),
                ),
                "mitomap_pubmed_count": self._pubmed_count(row.get("pubmed_ids", "")),
            })
        return clean


# ── ClinvarParser ─────────────────────────────────────────────────────────────

class ClinvarParser:
    """
    Parse ClinVar variant_summary.txt.gz.

    New fields vs. original:
      - source_db_version, source_record_id
      - hgvs_c, hgvs_p (extracted from Name field)
      - transcript_id (from NM_... in Name field)
      - clinvar_submitters_n, clinvar_last_evaluated, clinvar_conflicting
      - disease_terms: pipe-split list
      - genomic_ref / genomic_alt kept separate from CDS ref/alt
    """

    def __init__(self, hgnc_reference=None):
        self.hgnc_reference = hgnc_reference
        self._3to1 = {
            "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q",
            "Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K",
            "Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
            "Tyr":"Y","Val":"V",
        }

    @staticmethod
    def _extract_version(file_path: Path) -> str:
        m = re.search(r"(\d{4}-\d{2}-\d{2})", file_path.name)
        return m.group(1) if m else ""

    def _parse_protein(self, name: str):
        m = re.search(r"p\.([A-Z][a-z]{2}|\*)(\d+)([A-Z][a-z]{2}|\*)", str(name))
        if not m or "*" in m.groups():
            return None, False
        wt  = self._3to1.get(m.group(1), "X")
        mut = self._3to1.get(m.group(3), "X")
        return f"{wt}{m.group(2)}{mut}", wt == mut

    @staticmethod
    def _extract_hgvs_c(name: str) -> tuple[str, str, str, str]:
        """
        Extract (hgvs_c, transcript_id, cds_ref, cds_alt) from the Name field.
        Name format: NM_001369.4(NDUFS1):c.1057C>T (p.Arg353Cys)
        Returns ("NM_001369.4:c.1057C>T", "NM_001369.4", "C", "T") or ("","","","")
        """
        m_tx = re.search(r"(NM_[\d.]+|NR_[\d.]+|NP_[\d.]+)", str(name))
        m_c  = re.search(r"c\.\d+([ACGT])>([ACGT])", str(name))
        transcript_id = m_tx.group(1) if m_tx else ""
        if m_c:
            cds_ref = m_c.group(1)
            cds_alt = m_c.group(2)
            hgvs_c  = f"{transcript_id}:{m_c.group(0)}" if transcript_id else m_c.group(0)
        else:
            cds_ref = cds_alt = hgvs_c = ""
        return hgvs_c, transcript_id, cds_ref, cds_alt

    @staticmethod
    def _extract_hgvs_p(name: str) -> str:
        m = re.search(r"(p\.[A-Za-z0-9]+)", str(name))
        return m.group(1) if m else ""

    @staticmethod
    def _stars(review: str) -> int:
        r = review.lower()
        if "guideline" in r:        return 4
        if "expert" in r:           return 3
        if "multiple" in r and "no conflicts" in r: return 2
        if "single submitter" in r or "conflicting" in r: return 1
        return 0

    def parse(self, file_path: Path) -> list[dict]:
        version = self._extract_version(file_path)
        clean = []

        for chunk in pd.read_csv(
            file_path, sep="\t", compression="gzip",
            chunksize=50000, low_memory=False,
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
                    gene = self.hgnc_reference.lookup[gene].get("symbol", gene)

                aa_change, is_syn = self._parse_protein(row["Name"])
                if not aa_change or "X" in aa_change:
                    continue

                _aa_m  = re.match(r"^([A-Z])(\d+)([A-Z])$", aa_change)
                ref_aa = _aa_m.group(1) if _aa_m else ""
                alt_aa = _aa_m.group(3) if _aa_m else ""

                name        = str(row["Name"])
                hgvs_c, transcript_id, cds_ref, cds_alt = self._extract_hgvs_c(name)
                hgvs_p      = self._extract_hgvs_p(name)
                review      = str(row.get("ReviewStatus", "")).strip()
                status      = str(row.get("ClinicalSignificance", "")).strip()
                phenotypes  = str(row.get("PhenotypeList", "")).strip()
                disease_terms = [
                    t.strip() for t in phenotypes.split("|") if t.strip()
                ] if phenotypes else []

                clean.append({
                    "source_db":            "ClinVar",
                    "source_db_version":    version,
                    "source_record_id":     str(row["VariationID"]),
                    "genome":               "nucDNA",
                    "chromosome":           str(row["Chromosome"]),
                    "locus":                gene,
                    "allele_id":            str(row["VariationID"]),
                    "nt_change":            name,
                    "hgvs_c":               hgvs_c,
                    "hgvs_p":               hgvs_p,
                    "transcript_id":        transcript_id,
                    "aa_change":            aa_change,
                    "ref_aa":               ref_aa,
                    "alt_aa":               alt_aa,
                    "is_synonymous":        is_syn,
                    "disease":              phenotypes,
                    "disease_terms":        disease_terms,
                    "clinical_status":      status,
                    "review_status":        review,
                    "stars":                self._stars(review),
                    "clinvar_submitters_n": _int(row.get("NumberSubmitters"), 0),
                    "clinvar_last_evaluated": str(row.get("LastEvaluated", "")).strip(),
                    "clinvar_conflicting":  "conflicting" in status.lower(),
                    "grch38_pos":           int(row["Start"]),
                    "ref":                  cds_ref or str(row["ReferenceAlleleVCF"]),
                    "alt":                  cds_alt or str(row["AlternateAlleleVCF"]),
                    "genomic_ref":          str(row["ReferenceAlleleVCF"]),
                    "genomic_alt":          str(row["AlternateAlleleVCF"]),
                })
        return clean


# ── MyVariantParser ───────────────────────────────────────────────────────────

class MyVariantParser:
    """
    Expanded extraction from pre-downloaded MyVariant JSON.

    Returns a full metrics dict instead of just (revel, af).
    """

    _POP_ORDER = [
        "af_afr", "af_amr", "af_asj", "af_eas", "af_fin", "af_nfe", "af_oth",
    ]
    _POP_LABEL = {
        "af_afr": "afr", "af_amr": "amr", "af_asj": "asj",
        "af_eas": "eas", "af_fin": "fin", "af_nfe": "nfe", "af_oth": "oth",
    }

    def __init__(self, json_path: Path):
        with open(json_path) as f:
            self._lookup = {
                r["_id"]: r for r in json.load(f) if not r.get("notfound")
            }

    def _gnomad(self, rec: dict) -> dict:
        """Extract global + popmax AF, AC, AN, nhomalt from gnomad_genome or exome."""
        g = rec.get("gnomad_genome") or rec.get("gnomad_exome") or {}
        if not g:
            return {}

        af_dict  = g.get("af", {})
        ac_dict  = g.get("ac", {})
        an_dict  = g.get("an", {})
        hom_dict = g.get("hom", {})

        af_global = _float(af_dict.get("af") if isinstance(af_dict, dict) else af_dict)
        ac_global = _int(ac_dict.get("ac")   if isinstance(ac_dict, dict) else ac_dict)
        an_global = _int(an_dict.get("an")   if isinstance(an_dict, dict) else an_dict)
        nhomalt   = _int(hom_dict.get("hom") if isinstance(hom_dict, dict) else hom_dict)

        # Popmax: highest AF across continental populations
        popmax_af, popmax_pop = None, None
        if isinstance(af_dict, dict):
            for pop_key in self._POP_ORDER:
                paf = _float(af_dict.get(pop_key))
                if paf is not None and (popmax_af is None or paf > popmax_af):
                    popmax_af  = paf
                    popmax_pop = self._POP_LABEL.get(pop_key, pop_key)

        return {
            "gnomad_af_global":  af_global,
            "gnomad_af_popmax":  popmax_af,
            "gnomad_popmax_pop": popmax_pop,
            "gnomad_ac":         ac_global,
            "gnomad_an":         an_global,
            "gnomad_nhomalt":    nhomalt,
        }

    def get_all_metrics(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> dict:
        key = f"chr{chrom}:g.{pos}{ref}>{alt}"
        rec = self._lookup.get(key, {})
        d   = rec.get("dbnsfp", {})

        # AlphaMissense
        am = d.get("alphamissense", {})
        am_score = _float(_first(am.get("score") if isinstance(am, dict) else None))
        am_pred  = _first(am.get("pred") if isinstance(am, dict) else None)
        am_class = None
        if am_pred:
            lp = am_pred.lower()
            am_class = (
                "likely_pathogenic" if lp in ("lp", "p", "pathogenic", "likely pathogenic")
                else "likely_benign" if lp in ("lb", "b", "benign", "likely benign")
                else "ambiguous"
            )

        # ESM1b (lower = more deleterious)
        esm1b = d.get("esm1b", {})
        esm1b_score = _float(_first(
            esm1b.get("score") if isinstance(esm1b, dict) else esm1b
        ))

        result = {
            "revel_score":         _float(_first(d.get("revel", {}).get("score") if isinstance(d.get("revel"), dict) else d.get("revel"))),
            "phylop_100vert":      _float(d.get("phylop", {}).get("100way_vertebrate", {}).get("score") if isinstance(d.get("phylop"), dict) else None),
            "gerp_rs":             _float(d.get("gerp++", {}).get("rs") if isinstance(d.get("gerp++"), dict) else None),
            "alphamissense_score": am_score,
            "alphamissense_class": am_class,
            "esm1b_score":         esm1b_score,
            "mpc_score":           _float(d.get("mpc", {}).get("score") if isinstance(d.get("mpc"), dict) else None),
        }
        result.update(self._gnomad(rec))
        return result

    # Backward-compat shim for old call sites
    def get_metrics(self, chrom, pos, ref, alt):
        m = self.get_all_metrics(chrom, pos, ref, alt)
        return m.get("revel_score") or 0.0, m.get("gnomad_af_global") or 0.0


# ── PhylotreeParser ───────────────────────────────────────────────────────────

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


# ── MitimpactParser ───────────────────────────────────────────────────────────

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
