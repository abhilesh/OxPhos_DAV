#!/usr/bin/env python3
"""
Classify curated OXPHOS variants as cDAV or uDAV at the amino-acid and
nucleotide levels.

This stage is filter-late:
- every curated row is retained
- ineligible rows are marked as skipped by policy
- unresolved rows are retained with explicit mismatch metadata
- classified rows carry both AA-level and NT-level cDAV/uDAV calls
"""

from pathlib import Path
import csv
import hashlib
import json
import logging
import math
import re
import sys
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from Bio.Seq import Seq

from utils.alignment_parser import AlignmentParser


# ── Paths ────────────────────────────────────────────────────────────────────

DATA_DIR = ROOT / "data"
DERIVED_CURATED_DIR = DATA_DIR / "derived" / "curated"
DERIVED_CLASSIFIED_DIR = DATA_DIR / "derived" / "classified"
COMPAT_CURATED_DIR = DATA_DIR / "annotations" / "curated"
REFERENCE_DIR = DATA_DIR / "reference"
DERIVED_REFERENCE_DIR = DATA_DIR / "derived" / "reference"
LOG_DIR = DATA_DIR / "logs"

CURATED_MASTER_PARQUET = DERIVED_CURATED_DIR / "variants_master_curated.parquet"
CURATED_MASTER_JSONL = DERIVED_CURATED_DIR / "variants_master_curated.jsonl"
TRANSCRIPT_MAPS_JSON = DERIVED_CURATED_DIR / "transcript_position_maps.json"
GENOMIC_MAPS_JSON = DERIVED_CURATED_DIR / "genomic_coordinate_maps.json"

CLASSIFIED_MASTER_PARQUET = DERIVED_CLASSIFIED_DIR / "variants_master_classified.parquet"
CLASSIFIED_MASTER_JSONL = DERIVED_CLASSIFIED_DIR / "variants_master_classified.jsonl"
CLASSIFIED_ALL_PARQUET = DERIVED_CLASSIFIED_DIR / "classified_all.parquet"
CLASSIFIED_CLEAN_PARQUET = DERIVED_CLASSIFIED_DIR / "classified_clean.parquet"
CLASSIFIED_WARNING_PARQUET = DERIVED_CLASSIFIED_DIR / "classified_warning.parquet"
CLASSIFICATION_RUN_METADATA_JSON = DERIVED_CLASSIFIED_DIR / "classification_run_metadata.json"
CLASSIFICATION_MISMATCH_JSONL = DERIVED_CLASSIFIED_DIR / "classification_mismatch_log.jsonl"
CLASSIFICATION_QC_SUMMARY_JSON = DERIVED_CLASSIFIED_DIR / "classification_qc_summary.json"
CLASSIFICATION_QC_CHECKS_TSV = DERIVED_CLASSIFIED_DIR / "classification_qc_checks.tsv"

COMPAT_MT_JSON = COMPAT_CURATED_DIR / "cdav_classifications_mtDNA.json"
COMPAT_NUC_JSON = COMPAT_CURATED_DIR / "cdav_classifications_nucDNA.json"
OUT_MISMATCH_JSON = LOG_DIR / "cdav_mismatches.json"
OUT_LOG_FILE = LOG_DIR / "cdav_classification.log"
EXCEPTION_REGISTRY = DERIVED_REFERENCE_DIR / "variant_exception_registry.tsv"

MT_COORD_FILE = DERIVED_REFERENCE_DIR / "mtdna_gene_coordinates.tsv"
if not MT_COORD_FILE.exists():
    MT_COORD_FILE = REFERENCE_DIR / "mtdna_gene_coordinates.tsv"

TOGA_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_aa"
TOGA_NT_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"
MT_AA_DIR = DATA_DIR / "alignments" / "mtdna_aa"
MT_NT_DIR = DATA_DIR / "alignments" / "mtdna_codon"


# ── Constants ────────────────────────────────────────────────────────────────

REF_ALLELE_MISMATCH = "REF_ALLELE_MISMATCH"
TRANSCRIPT_MISMATCH = "TRANSCRIPT_MISMATCH"
ANCHOR_NOT_FOUND = "ANCHOR_NOT_FOUND"
POSITION_NOT_IN_ENST = "POSITION_NOT_IN_ENST"
GENOMIC_POS_NOT_IN_ENST = "GENOMIC_POS_NOT_IN_ENST"
CODON_EXTRACTION_FAILURE = "CODON_EXTRACTION_FAILURE"
NO_ALIGNMENT = "NO_ALIGNMENT"
COORD_PARSE_FAILURE = "COORD_PARSE_FAILURE"

_MT_ALIAS = {
    "MT-COX1": "MT-CO1",
    "MT-COX2": "MT-CO2",
    "MT-COX3": "MT-CO3",
    "MT-CYTB": "MT-CYB",
}

_MISMATCH_DESCRIPTIONS = {
    REF_ALLELE_MISMATCH:
        "reference base in the curated row does not match the human CDS base at the corrected alignment position",
    TRANSCRIPT_MISMATCH:
        "corrected mutant codon translates to an amino acid different from the curated amino-acid consequence",
    ANCHOR_NOT_FOUND:
        "wild-type amino acid was not found within the fallback anchor window in the alignment",
    POSITION_NOT_IN_ENST:
        "the NM residue maps to a gap in the TOGA ENST isoform",
    GENOMIC_POS_NOT_IN_ENST:
        "the genomic coordinate does not fall within the CDS of the TOGA ENST rescue map",
    CODON_EXTRACTION_FAILURE:
        "the corrected coding position could not be converted into a valid codon in the nucleotide alignment",
    NO_ALIGNMENT:
        "required AA or codon alignment file is absent for this interpreted gene",
    COORD_PARSE_FAILURE:
        "the classifier could not parse a valid amino-acid position and residue change from curated fields",
}


# ── Setup ────────────────────────────────────────────────────────────────────

def ensure_layout() -> None:
    for path in (DERIVED_CLASSIFIED_DIR, COMPAT_CURATED_DIR, LOG_DIR):
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> logging.Logger:
    ensure_layout()
    logger = logging.getLogger("cdav_classify")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(fmt)
    logger.addHandler(console)

    file_handler = logging.FileHandler(OUT_LOG_FILE, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)

    return logger


def repo_rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


# ── Loaders ──────────────────────────────────────────────────────────────────

def load_curated_rows():
    import pandas as pd

    if not CURATED_MASTER_PARQUET.exists():
        raise FileNotFoundError(f"Missing curated master table: {CURATED_MASTER_PARQUET}")
    df = pd.read_parquet(CURATED_MASTER_PARQUET)
    return [normalize_json_value(row) for row in df.to_dict(orient="records")]


def load_transcript_maps(logger: logging.Logger) -> dict[str, dict]:
    if not TRANSCRIPT_MAPS_JSON.exists():
        logger.warning("Transcript maps not found at %s", TRANSCRIPT_MAPS_JSON)
        return {}
    with open(TRANSCRIPT_MAPS_JSON, encoding="utf-8") as handle:
        raw = json.load(handle)
    tx_maps: dict[str, dict] = {}
    for gene, entry in raw.items():
        if entry.get("type") in {"mapped", "identity"}:
            tx_maps[gene] = {
                "map": {int(k): v for k, v in entry.get("map", {}).items()},
                "identity_fraction": entry.get("identity_fraction"),
                "coverage_fraction": entry.get("coverage_fraction"),
                "type": entry.get("type"),
                "nm": entry.get("nm"),
                "enst": entry.get("enst"),
            }
    logger.info("Loaded transcript maps for %d genes.", len(tx_maps))
    return tx_maps


def load_genomic_maps(logger: logging.Logger) -> dict[str, dict[int, dict]]:
    if not GENOMIC_MAPS_JSON.exists():
        logger.warning("Genomic rescue maps not found at %s", GENOMIC_MAPS_JSON)
        return {}
    with open(GENOMIC_MAPS_JSON, encoding="utf-8") as handle:
        raw = json.load(handle)
    genomic_maps: dict[str, dict[int, dict]] = {}
    for gene, entry in raw.items():
        genomic_maps[gene] = {int(k): v for k, v in entry.get("map", {}).items()}
    logger.info("Loaded genomic rescue maps for %d genes.", len(genomic_maps))
    return genomic_maps


def load_mt_coords() -> dict[str, tuple[int, int, str]]:
    coords: dict[str, tuple[int, int, str]] = {}
    if not MT_COORD_FILE.exists():
        return coords
    with open(MT_COORD_FILE, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            start = int(row["start"])
            end = int(row["end"])
            entry = (min(start, end), max(start, end), row["strand"])
            coords[row["gene"]] = entry
            if row["gene"] in _MT_ALIAS:
                coords[_MT_ALIAS[row["gene"]]] = entry
    return coords


MT_COORDS = load_mt_coords()


# ── Helpers ──────────────────────────────────────────────────────────────────

def mismatch_record(mtype: str, genome: str, locus: str, var: dict, detail: str) -> dict:
    return {
        "type": mtype,
        "genome": genome,
        "locus": locus,
        "variant_id": var.get("variant_id", var.get("ann_id", "")),
        "source_variant_group_id": var.get("source_variant_group_id", ""),
        "tier": var.get("tier", ""),
        "aa_change": var.get("aa_change", ""),
        "detail": detail,
    }


def normalize_json_value(value):
    if isinstance(value, dict):
        return {str(k): normalize_json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [normalize_json_value(v) for v in value]

    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        converted = tolist()
        if converted is not value:
            return normalize_json_value(converted)

    item = getattr(value, "item", None)
    if callable(item):
        try:
            converted = item()
        except Exception:
            converted = value
        if converted is not value:
            return normalize_json_value(converted)

    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def write_jsonl(rows: list[dict], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(normalize_json_value(row)) + "\n")


def genomic_to_cds(genomic_pos: int, locus: str) -> int | None:
    if locus not in MT_COORDS:
        return None
    start, end, strand = MT_COORDS[locus]
    return (genomic_pos - start + 1) if strand == "+" else (end - genomic_pos + 1)


def parse_aa_position_from_text(text: str) -> int | None:
    if not text:
        return None
    match = re.search(r"p\.[A-Za-z*]+(\d+)[A-Za-z*]+", text)
    if match:
        return int(match.group(1))
    match = re.search(r"[A-Za-z*]+(\d+)[A-Za-z*]+", text)
    if match:
        return int(match.group(1))
    return None


def parse_variant_coordinates(var: dict) -> tuple[int | None, str | None, str | None, int | None]:
    wt_aa = normalize_json_value(var.get("frame_specific_ref_aa"))
    if not wt_aa:
        wt_aa = normalize_json_value(var.get("ref_aa"))
    mut_aa = normalize_json_value(var.get("frame_specific_alt_aa"))
    if not mut_aa:
        mut_aa = normalize_json_value(var.get("alt_aa"))

    aa_pos = None
    if var.get("frame_specific_codon_index") not in (None, ""):
        try:
            aa_pos = int(var["frame_specific_codon_index"])
        except (TypeError, ValueError):
            aa_pos = None
    if aa_pos is None:
        for field in ("frame_specific_hgvs_p", "hgvs_p", "aa_change"):
            aa_pos = parse_aa_position_from_text(str(var.get(field, "") or ""))
            if aa_pos is not None:
                break

    nt_pos = None
    if var.get("genome") == "nucDNA":
        nc_text = str(var.get("hgvs_c") or var.get("nc_change") or "")
        match = re.search(r"c\.(\d+)", nc_text)
        if match:
            nt_pos = int(match.group(1))
    else:
        genomic_pos = var.get("genomic_pos")
        locus = var.get("interpreted_gene") or var.get("locus", "")
        if genomic_pos not in (None, "") and locus:
            nt_pos = genomic_to_cds(int(genomic_pos), locus)

    return aa_pos, wt_aa, mut_aa, nt_pos


def get_alignment_paths(locus: str, genome: str) -> tuple[Path, Path]:
    if genome == "nucDNA":
        return (
            TOGA_AA_DIR / f"{locus}_aa_alignment.fasta",
            TOGA_NT_DIR / f"{locus}_codon_alignment.fasta",
        )
    return (
        MT_AA_DIR / f"{locus}_aa_alignment.fasta",
        MT_NT_DIR / f"{locus}_codon_alignment.fasta",
    )


def classification_defaults(var: dict) -> dict:
    classification_gene = var.get("interpreted_gene") or var.get("locus", "")
    return {
        "classification_status": None,
        "classification_eligible": bool(var.get("eligible_core_comparative_pipeline")),
        "classification_exclusion_reason": None,
        "classification_gene": classification_gene,
        "classification_gene_source": "interpreted_gene" if var.get("interpreted_gene") else "locus",
        "classification_coordinate_method": None,
        "classification_coordinate_status": None,
        "classification_transcript_map_identity": None,
        "classification_transcript_map_coverage": None,
        "exception_applied": bool(var.get("exception_class")),
        "classification_exception_action": var.get("exception_decision"),
        "classification_used_overlap_frame": bool(var.get("is_overlap")),
        "classification_basis": None,
        "is_cdav_amino_acid": None,
        "is_cdav_nucleotide": None,
        "is_udav_amino_acid": None,
        "is_udav_nucleotide": None,
        "n_species_aligned": None,
        "n_species_with_disease_allele": None,
        "lineages_with_disease_allele": [],
        "n_species_with_disease_codon": None,
        "lineages_with_disease_codon": [],
        "cdav_support_level": None,
        "ref_allele_match": None,
        "mismatch_reason": None,
        "mismatch_code": None,
        "alignment_file_aa": None,
        "alignment_file_nt": None,
        "alignment_source": "TOGA" if var.get("genome") == "nucDNA" else "mtDNA_MACSE",
    }


def support_level(n_species: int | None) -> str | None:
    if n_species is None:
        return None
    if n_species == 0:
        return "none"
    if n_species == 1:
        return "single_species"
    return "multi_species"


def prefer_genomic_map(tx_entry: dict | None) -> bool:
    if not tx_entry:
        return True
    identity = tx_entry.get("identity_fraction")
    coverage = tx_entry.get("coverage_fraction")
    if identity is None or coverage is None:
        return True
    return identity < 0.90 or coverage < 0.95


def set_skipped_by_policy(var: dict) -> dict:
    var["classification_status"] = "skipped_by_policy"
    var["classification_eligible"] = False
    var["classification_exclusion_reason"] = var.get("core_pipeline_exclusion_reason") or "core_pipeline_ineligible"
    var["classification_coordinate_status"] = "not_attempted"
    var["classification_basis"] = "unresolved"
    return var


def set_unresolved(var: dict, reason: str, code: str, detail: str, mismatches: list[dict], logger: logging.Logger) -> dict:
    var["classification_status"] = "unresolved"
    var["classification_exclusion_reason"] = reason
    var["classification_basis"] = "unresolved"
    var["mismatch_reason"] = detail
    var["mismatch_code"] = code
    mismatches.append(mismatch_record(code, var["genome"], var["classification_gene"], var, detail))
    logger.warning("[%s] %s | %s | %s | %s", var["genome"], code, var["classification_gene"], var.get("variant_id"), detail)
    return var


def process_variants(
    variants: list[dict],
    loaded_alignments: dict,
    mismatches: list[dict],
    logger: logging.Logger,
    tx_maps: dict[str, dict],
    genomic_maps: dict[str, dict[int, dict]],
) -> tuple[list[dict], dict[str, dict[str, int]]]:
    enriched: list[dict] = []
    summary: dict[str, dict[str, int]] = {
        "mtDNA": defaultdict(int),
        "nucDNA": defaultdict(int),
    }
    missing_loci: dict[tuple[str, str], int] = {}

    total = len(variants)
    for idx, raw_var in enumerate(variants):
        if idx % 500 == 0:
            logger.info("Processed %d / %d variants...", idx, total)

        var = dict(raw_var)
        var = normalize_json_value(var)
        var.update(classification_defaults(var))
        genome = var["genome"]
        locus = var["classification_gene"]
        tx_entry = tx_maps.get(locus)
        if tx_entry:
            var["classification_transcript_map_identity"] = tx_entry.get("identity_fraction")
            var["classification_transcript_map_coverage"] = tx_entry.get("coverage_fraction")
        summary[genome]["rows_total"] += 1

        if not var.get("eligible_core_comparative_pipeline", False):
            enriched.append(set_skipped_by_policy(var))
            summary[genome]["skipped_by_policy"] += 1
            continue

        aa_pos, wt_aa, mut_aa, nt_pos = parse_variant_coordinates(var)
        if not aa_pos or not wt_aa or not mut_aa or not nt_pos:
            detail = (
                f"aa_pos={aa_pos} wt_aa={wt_aa} mut_aa={mut_aa} nt_pos={nt_pos} "
                f"variant_id={var.get('variant_id')}"
            )
            var["classification_coordinate_status"] = "parse_failed"
            enriched.append(set_unresolved(var, "coordinate_parse_failure", COORD_PARSE_FAILURE, detail, mismatches, logger))
            summary[genome]["unresolved"] += 1
            continue

        used_genomic_map = False
        genomic_map_preferred = genome == "nucDNA" and locus in genomic_maps and prefer_genomic_map(tx_entry)
        if genomic_map_preferred:
            gpos = var.get("genomic_pos")
            if gpos is not None:
                gmap_entry = genomic_maps[locus].get(int(gpos))
                if gmap_entry is None:
                    detail = f"locus={locus} genomic_pos={gpos} not in TOGA ENST CDS rescue map"
                    var["classification_coordinate_method"] = "genomic_map"
                    var["classification_coordinate_status"] = "genomic_pos_not_in_enst"
                    enriched.append(set_unresolved(var, "genomic_pos_not_in_enst", GENOMIC_POS_NOT_IN_ENST, detail, mismatches, logger))
                    summary[genome]["unresolved"] += 1
                    continue
                aa_pos = gmap_entry["aa_pos"]
                nt_pos = gmap_entry["cds_pos"]
                used_genomic_map = True

        alignment_key = (genome, locus)
        if alignment_key not in loaded_alignments:
            aa_path, nt_path = get_alignment_paths(locus, genome)
            if not aa_path.exists() or not nt_path.exists():
                loaded_alignments[alignment_key] = None
            else:
                tx_for_locus = None if used_genomic_map else tx_entry.get("map") if tx_entry else None
                loaded_alignments[alignment_key] = AlignmentParser(aa_path, nt_path, genome, tx_pos_map=tx_for_locus)

        parser = loaded_alignments[alignment_key]
        aa_path, nt_path = get_alignment_paths(locus, genome)
        var["alignment_file_aa"] = repo_rel(aa_path) if aa_path.exists() else None
        var["alignment_file_nt"] = repo_rel(nt_path) if nt_path.exists() else None

        if parser is None:
            missing_loci[(genome, locus)] = missing_loci.get((genome, locus), 0) + 1
            var["classification_coordinate_status"] = "missing_alignment"
            enriched.append(set_unresolved(var, "missing_alignment", NO_ALIGNMENT, f"locus={locus}", mismatches, logger))
            summary[genome]["unresolved"] += 1
            continue

        var["classification_coordinate_method"] = (
            "genomic_map" if used_genomic_map else
            "transcript_map" if parser.tx_pos_map is not None else
            "anchor_fallback"
        )
        var["classification_coordinate_status"] = "attempted"

        ref_nt = str(normalize_json_value(var.get("coding_ref_nt")) or normalize_json_value(var.get("ref_nt")) or "").upper()
        alt_nt = str(normalize_json_value(var.get("coding_alt_nt")) or normalize_json_value(var.get("alt_nt")) or "").upper()
        comp = parser.check_compensation(aa_pos, wt_aa, mut_aa, nt_pos, alt_nt)
        var["n_species_aligned"] = max(len(parser.aa_alignment) - 1, 0)

        if not comp["anchor_found"]:
            if comp.get("position_not_in_enst"):
                reason = "position_not_in_enst"
                code = POSITION_NOT_IN_ENST
                detail = f"locus={locus} aa_pos={aa_pos} wt_aa={wt_aa} aa_change={var.get('aa_change')}"
            else:
                reason = "anchor_not_found"
                code = ANCHOR_NOT_FOUND
                detail = f"locus={locus} aa_pos={aa_pos} wt_aa={wt_aa} aa_change={var.get('aa_change')}"
            var["classification_coordinate_status"] = reason
            enriched.append(set_unresolved(var, reason, code, detail, mismatches, logger))
            summary[genome]["unresolved"] += 1
            continue

        ref_base = comp["ref_base_found"]
        ref_ok = ref_base not in ("POS_NOT_IN_MAP",) and ref_base == ref_nt
        var["ref_allele_match"] = bool(ref_ok)
        ref_mismatch_detail = None
        if not ref_ok:
            ref_mismatch_detail = (
                f"locus={locus} raw_cds={nt_pos} corrected_cds={comp['corrected_nt_pos']} "
                f"expected={ref_nt} found={ref_base}"
            )
            mismatches.append(mismatch_record(REF_ALLELE_MISMATCH, genome, locus, var, ref_mismatch_detail))
            logger.warning("[%s] %s | %s | %s | %s", genome, REF_ALLELE_MISMATCH, locus, var.get("variant_id"), ref_mismatch_detail)

        mut_codon = comp["mut_codon"]
        if mut_codon is None:
            detail = f"locus={locus} corrected_cds={comp['corrected_nt_pos']} alt_nt={alt_nt}"
            var["classification_coordinate_status"] = "codon_extraction_failure"
            enriched.append(set_unresolved(var, "codon_extraction_failure", CODON_EXTRACTION_FAILURE, detail, mismatches, logger))
            summary[genome]["unresolved"] += 1
            continue

        if genome == "nucDNA":
            try:
                translated_aa = str(Seq(mut_codon).translate(table=1))
            except Exception:
                translated_aa = "ERR"
            expected_aa = mut_aa
            if translated_aa != expected_aa:
                detail = (
                    f"locus={locus} raw_cds={nt_pos} corrected_cds={comp['corrected_nt_pos']} "
                    f"codon={mut_codon} translated={translated_aa} expected={expected_aa} toga_tx={parser.transcript_id}"
                )
                var["classification_coordinate_status"] = "transcript_mismatch"
                enriched.append(set_unresolved(var, "transcript_mismatch", TRANSCRIPT_MISMATCH, detail, mismatches, logger))
                summary[genome]["unresolved"] += 1
                continue

        aa_cdav = bool(comp["aa_cdav"])
        nt_cdav = bool(comp["nt_cdav"])
        aa_species = comp["aa_species"]
        nt_species = comp["nt_species"]

        var.update({
            "classification_status": "classified",
            "classification_eligible": True,
            "classification_exclusion_reason": None,
            "classification_coordinate_status": "classified",
            "classification_basis": (
                "nt_and_aa" if nt_cdav and aa_cdav else
                "nt_only" if nt_cdav else
                "aa_only" if aa_cdav else
                "no_disease_allele_detected"
            ),
            "is_cdav_amino_acid": aa_cdav,
            "is_cdav_nucleotide": nt_cdav,
            "is_udav_amino_acid": not aa_cdav,
            "is_udav_nucleotide": not nt_cdav,
            "n_species_with_disease_allele": len(aa_species),
            "lineages_with_disease_allele": aa_species,
            "n_species_with_disease_codon": len(nt_species),
            "lineages_with_disease_codon": nt_species,
            "cdav_support_level": support_level(len(aa_species)),
            "mismatch_reason": ref_mismatch_detail,
            "mismatch_code": None if ref_ok else REF_ALLELE_MISMATCH,
        })
        enriched.append(var)

        summary[genome]["classified"] += 1
        if aa_cdav:
            summary[genome]["aa_cdav"] += 1
        if nt_cdav:
            summary[genome]["nt_cdav"] += 1
        if var["is_udav_amino_acid"]:
            summary[genome]["aa_udav"] += 1
        if var["is_udav_nucleotide"]:
            summary[genome]["nt_udav"] += 1

    for (genome, locus), count in sorted(missing_loci.items()):
        detail = f"locus={locus} affected_variants={count}"
        mismatches.append({
            "type": NO_ALIGNMENT,
            "genome": genome,
            "locus": locus,
            "variant_id": f"({count} variants)",
            "source_variant_group_id": "",
            "tier": "",
            "aa_change": "",
            "detail": detail,
        })

    return enriched, summary


def warning_reasons(row: dict) -> list[str]:
    reasons: list[str] = []
    if row.get("exception_applied"):
        reasons.append("exception_applied")
    if row.get("mismatch_code") == REF_ALLELE_MISMATCH:
        reasons.append(REF_ALLELE_MISMATCH)
    return reasons


def apply_classification_subsets(rows: list[dict]) -> tuple[list[dict], list[dict], list[dict], dict]:
    classified_all: list[dict] = []
    classified_clean: list[dict] = []
    classified_warning: list[dict] = []

    for row in rows:
        if row.get("classification_status") != "classified":
            row["classification_subset"] = None
            row["classification_warning_reason"] = None
            row["is_clean_classified"] = None
            continue

        reasons = warning_reasons(row)
        classified_all.append(row)
        if reasons:
            row["classification_subset"] = "classified_warning"
            row["classification_warning_reason"] = "|".join(reasons)
            row["is_clean_classified"] = False
            classified_warning.append(row)
        else:
            row["classification_subset"] = "classified_clean"
            row["classification_warning_reason"] = None
            row["is_clean_classified"] = True
            classified_clean.append(row)

    counts = {
        "classified_all": len(classified_all),
        "classified_clean": len(classified_clean),
        "classified_warning": len(classified_warning),
    }
    return classified_all, classified_clean, classified_warning, counts


def qc_result(check_name: str, failed_rows: list[dict], failure_rule: str) -> dict:
    example_variant_ids: list[str] = []
    for row in failed_rows[:10]:
        variant_id = row.get("variant_id") or row.get("ann_id") or ""
        if variant_id:
            example_variant_ids.append(str(variant_id))
    return {
        "check_name": check_name,
        "status": "pass" if not failed_rows else "fail",
        "n_failed": len(failed_rows),
        "failure_rule": failure_rule,
        "example_variant_ids": example_variant_ids,
    }


def build_qc_checks(curated_rows: list[dict], enriched_rows: list[dict]) -> list[dict]:
    valid_basis = {"nt_and_aa", "nt_only", "aa_only", "no_disease_allele_detected"}
    curated_overlap_mt = {
        row.get("variant_id")
        for row in curated_rows
        if row.get("genome") == "mtDNA" and row.get("derived_from_overlap_duplication")
    }
    classified_overlap_mt = {
        row.get("variant_id")
        for row in enriched_rows
        if row.get("genome") == "mtDNA" and row.get("derived_from_overlap_duplication")
    }

    checks = [
        qc_result(
            "no_clinvar_mitochondrial_row_eligible_for_nucdna_branch",
            [
                row for row in enriched_rows
                if row.get("source_db") == "ClinVar"
                and row.get("encoded_by") == "mitochondrial"
                and row.get("eligible_core_comparative_pipeline")
            ],
            "ClinVar mitochondrial-encoded rows must remain ineligible for the nuclear comparative branch.",
        ),
        qc_result(
            "no_overlap_derived_mtdna_row_lost_during_classification",
            [
                {"variant_id": variant_id}
                for variant_id in sorted(curated_overlap_mt - classified_overlap_mt)
            ],
            "Every curated overlap-derived mtDNA row must appear in the classified master table.",
        ),
        qc_result(
            "no_unresolved_row_counted_as_udav",
            [
                row for row in enriched_rows
                if row.get("classification_status") == "unresolved"
                and (row.get("is_udav_amino_acid") is True or row.get("is_udav_nucleotide") is True)
            ],
            "Rows with classification_status == unresolved must not be counted as uDAV.",
        ),
        qc_result(
            "every_classified_row_has_valid_classification_basis",
            [
                row for row in enriched_rows
                if row.get("classification_status") == "classified"
                and row.get("classification_basis") not in valid_basis
            ],
            "Classified rows must carry one of nt_and_aa, nt_only, aa_only, no_disease_allele_detected.",
        ),
        qc_result(
            "every_classified_row_has_resolved_alignment_source",
            [
                row for row in enriched_rows
                if row.get("classification_status") == "classified"
                and not row.get("alignment_source")
            ],
            "Classified rows must record a resolved alignment_source.",
        ),
        qc_result(
            "every_classified_row_has_resolved_coordinate_method",
            [
                row for row in enriched_rows
                if row.get("classification_status") == "classified"
                and not row.get("classification_coordinate_method")
            ],
            "Classified rows must record a resolved classification_coordinate_method.",
        ),
    ]
    return checks


# ── Reporting ────────────────────────────────────────────────────────────────

def print_mismatch_summary(mismatches: list[dict], logger: logging.Logger) -> None:
    by_type: dict[str, list[dict]] = defaultdict(list)
    for item in mismatches:
        by_type[item["type"]].append(item)

    logger.info("")
    logger.info("=" * 60)
    logger.info("MISMATCH / INCONSISTENCY SUMMARY")
    logger.info("=" * 60)
    logger.info("  Total mismatches logged : %d", len(mismatches))
    logger.info("")
    for mtype in [
        REF_ALLELE_MISMATCH,
        TRANSCRIPT_MISMATCH,
        ANCHOR_NOT_FOUND,
        POSITION_NOT_IN_ENST,
        GENOMIC_POS_NOT_IN_ENST,
        CODON_EXTRACTION_FAILURE,
        NO_ALIGNMENT,
        COORD_PARSE_FAILURE,
    ]:
        count = len(by_type.get(mtype, []))
        if count:
            logger.info("  %-28s : %d", mtype, count)
            logger.info("    %s", _MISMATCH_DESCRIPTIONS[mtype])
    logger.info("=" * 60)


def print_summary(summary: dict[str, dict[str, int]], logger: logging.Logger) -> None:
    logger.info("")
    logger.info("=" * 60)
    logger.info("cDAV CLASSIFICATION SUMMARY")
    logger.info("=" * 60)
    for genome in ("mtDNA", "nucDNA"):
        stats = summary[genome]
        logger.info("")
        logger.info("%s:", genome)
        logger.info("  Total rows              : %d", stats.get("rows_total", 0))
        logger.info("  Classified              : %d", stats.get("classified", 0))
        logger.info("  Unresolved              : %d", stats.get("unresolved", 0))
        logger.info("  Skipped by policy       : %d", stats.get("skipped_by_policy", 0))
        logger.info("  AA-level cDAV           : %d", stats.get("aa_cdav", 0))
        logger.info("  NT-level cDAV           : %d", stats.get("nt_cdav", 0))
        logger.info("  AA-level uDAV           : %d", stats.get("aa_udav", 0))
        logger.info("  NT-level uDAV           : %d", stats.get("nt_udav", 0))


# ── Entry point ──────────────────────────────────────────────────────────────

def main() -> None:
    import pandas as pd

    logger = setup_logging()
    logger.info("Initializing cDAV/uDAV classification from the curated master table")

    if not EXCEPTION_REGISTRY.exists():
        raise FileNotFoundError(
            f"Missing required exception registry stage contract: {EXCEPTION_REGISTRY}"
        )
    exception_registry_metadata = {
        "path": repo_rel(EXCEPTION_REGISTRY),
        "sha256": sha256_file(EXCEPTION_REGISTRY),
    }

    curated_rows = load_curated_rows()
    logger.info("Loaded %d curated rows from %s", len(curated_rows), CURATED_MASTER_PARQUET)

    tx_maps = load_transcript_maps(logger)
    genomic_maps = load_genomic_maps(logger)

    loaded_alignments: dict[str, AlignmentParser | None] = {}
    mismatches: list[dict] = []

    enriched_rows, summary = process_variants(
        curated_rows,
        loaded_alignments=loaded_alignments,
        mismatches=mismatches,
        logger=logger,
        tx_maps=tx_maps,
        genomic_maps=genomic_maps,
    )
    classified_all, classified_clean, classified_warning, subset_counts = apply_classification_subsets(enriched_rows)
    qc_checks = build_qc_checks(curated_rows, enriched_rows)
    qc_summary = {
        "exception_registry": exception_registry_metadata,
        "n_checks": len(qc_checks),
        "n_failed_checks": sum(1 for item in qc_checks if item["status"] != "pass"),
        "subset_counts": subset_counts,
        "checks": qc_checks,
    }

    df = pd.DataFrame(enriched_rows)
    df.to_parquet(CLASSIFIED_MASTER_PARQUET, index=False)
    pd.DataFrame(classified_all).to_parquet(CLASSIFIED_ALL_PARQUET, index=False)
    pd.DataFrame(classified_clean).to_parquet(CLASSIFIED_CLEAN_PARQUET, index=False)
    pd.DataFrame(classified_warning).to_parquet(CLASSIFIED_WARNING_PARQUET, index=False)
    write_jsonl(enriched_rows, CLASSIFIED_MASTER_JSONL)
    write_jsonl(mismatches, CLASSIFICATION_MISMATCH_JSONL)
    pd.DataFrame(qc_checks).to_csv(CLASSIFICATION_QC_CHECKS_TSV, sep="\t", index=False)
    with open(CLASSIFICATION_QC_SUMMARY_JSON, "w", encoding="utf-8") as handle:
        json.dump(normalize_json_value(qc_summary), handle, indent=2)

    mt_rows = [row for row in enriched_rows if row.get("genome") == "mtDNA"]
    nuc_rows = [row for row in enriched_rows if row.get("genome") == "nucDNA"]
    with open(COMPAT_MT_JSON, "w", encoding="utf-8") as handle:
        json.dump(normalize_json_value(mt_rows), handle, indent=2)
    with open(COMPAT_NUC_JSON, "w", encoding="utf-8") as handle:
        json.dump(normalize_json_value(nuc_rows), handle, indent=2)
    with open(OUT_MISMATCH_JSON, "w", encoding="utf-8") as handle:
        json.dump(normalize_json_value(mismatches), handle, indent=2)

    metadata = {
        "input_rows": len(curated_rows),
        "output_rows": len(enriched_rows),
        "canonical_input": repo_rel(CURATED_MASTER_PARQUET),
        "canonical_output": repo_rel(CLASSIFIED_MASTER_PARQUET),
        "subset_outputs": [
            repo_rel(CLASSIFIED_ALL_PARQUET),
            repo_rel(CLASSIFIED_CLEAN_PARQUET),
            repo_rel(CLASSIFIED_WARNING_PARQUET),
        ],
        "compatibility_outputs": [repo_rel(COMPAT_MT_JSON), repo_rel(COMPAT_NUC_JSON)],
        "mismatch_log": repo_rel(CLASSIFICATION_MISMATCH_JSONL),
        "qc_outputs": {
            "summary_json": repo_rel(CLASSIFICATION_QC_SUMMARY_JSON),
            "checks_tsv": repo_rel(CLASSIFICATION_QC_CHECKS_TSV),
        },
        "exception_registry": exception_registry_metadata,
        "summary": {genome: dict(stats) for genome, stats in summary.items()},
        "subset_counts": subset_counts,
    }
    with open(CLASSIFICATION_RUN_METADATA_JSON, "w", encoding="utf-8") as handle:
        json.dump(normalize_json_value(metadata), handle, indent=2)

    logger.info("Wrote classified master table: %s", CLASSIFIED_MASTER_PARQUET)
    logger.info("Wrote classified subset tables: %s, %s, %s", CLASSIFIED_ALL_PARQUET, CLASSIFIED_CLEAN_PARQUET, CLASSIFIED_WARNING_PARQUET)
    logger.info("Wrote classified JSONL view: %s", CLASSIFIED_MASTER_JSONL)
    logger.info("Wrote classification metadata: %s", CLASSIFICATION_RUN_METADATA_JSON)
    logger.info("Wrote classification QC summary: %s", CLASSIFICATION_QC_SUMMARY_JSON)
    print_summary(summary, logger)
    print_mismatch_summary(mismatches, logger)


if __name__ == "__main__":
    main()
