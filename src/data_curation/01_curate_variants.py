#!/usr/bin/env python3
"""
Curate disease-associated OXPHOS variants into the canonical master table.

This is the filter-late curation layer. Downloaded source rows are retained as
metadata-bearing records whenever they can be represented safely.
"""

from pathlib import Path
import csv
import gzip
import json
import re
import sys
from dataclasses import asdict
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from utils.gene_reference import GeneReference
from utils.exception_registry import load_exception_registry, match_exception_entry
from utils.mt_overlap import load_human_mt_sequences, load_mtdna_gene_coords, parse_simple_aa_change
from utils.variant_record import VariantRecord

DATA_DIR = ROOT / "data"
RAW_ANNOTATIONS_DIR = DATA_DIR / "raw" / "annotations"
RAW_REFERENCE_DIR = DATA_DIR / "raw" / "reference"
COMPAT_CURATED_DIR = DATA_DIR / "annotations" / "curated"
DERIVED_CURATED_DIR = DATA_DIR / "derived" / "curated"
DERIVED_REFERENCE_DIR = DATA_DIR / "derived" / "reference"
DERIVED_RESULTS_DIR = DATA_DIR / "derived" / "results"
MT_CODON_DIR = DATA_DIR / "alignments" / "mtdna_codon"
EXCEPTION_REGISTRY = DERIVED_REFERENCE_DIR / "variant_exception_registry.tsv"
if not EXCEPTION_REGISTRY.exists():
    EXCEPTION_REGISTRY = DATA_DIR / "reference" / "variant_exception_registry.tsv"

MASTER_CURATED_PARQUET = DERIVED_CURATED_DIR / "variants_master_curated.parquet"
MASTER_CURATED_JSONL = DERIVED_CURATED_DIR / "variants_master_curated.jsonl"
EXCLUDED_INVENTORY_PARQUET = DERIVED_CURATED_DIR / "variants_excluded_inventory.parquet"
CURATION_RUN_METADATA_JSON = DERIVED_CURATED_DIR / "curation_run_metadata.json"
CURATION_CROSS_SOURCE_GROUPS_TSV = DERIVED_CURATED_DIR / "cross_source_duplicate_groups.tsv"
COMPAT_MT_JSON = COMPAT_CURATED_DIR / "mtDNA_annotations.json"
COMPAT_NUC_JSON = COMPAT_CURATED_DIR / "nucDNA_annotations.json"
CLINVAR_MITOMAP_OVERLAP_SUMMARY_JSON = DERIVED_RESULTS_DIR / "clinvar_mitomap_cross_source_overlap_summary.json"
CLINVAR_MT_PRESENT_TSV = DERIVED_RESULTS_DIR / "clinvar_mt_variants_present_in_mitomap.tsv"
CLINVAR_MT_ABSENT_TSV = DERIVED_RESULTS_DIR / "clinvar_mt_variants_absent_from_mitomap.tsv"


def ensure_layout() -> None:
    for path in (DERIVED_CURATED_DIR, COMPAT_CURATED_DIR, DERIVED_RESULTS_DIR):
        path.mkdir(parents=True, exist_ok=True)


def repo_rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def latest(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)
    if not matches:
        raise FileNotFoundError(f"Missing required file matching {pattern} in {directory}")
    return matches[0]


def extract_date(path: Path) -> str:
    match = re.search(r"(\d{4}-\d{2}-\d{2})", path.name)
    return match.group(1) if match else ""


def dedupe_terms(raw: str) -> list[str]:
    seen: list[str] = []
    for token in re.split(r"[|;]", raw):
        token = token.strip()
        if token and token not in seen:
            seen.append(token)
    return seen


def variant_flags(ref: str, alt: str, aa_change: str) -> dict:
    ref = str(ref).strip().upper()
    alt = str(alt).strip().upper()
    aa_lower = str(aa_change).strip().lower()
    is_snv = len(ref) == 1 and len(alt) == 1
    is_mnv = len(ref) == len(alt) and len(ref) > 1
    is_indel = len(ref) != len(alt)
    is_frameshift = "frameshift" in aa_lower
    is_nonsense = "*" in aa_change or "ter" in aa_lower
    is_coding = "noncoding" not in aa_lower and bool(aa_change)
    parsed = parse_simple_aa_change(aa_change) if aa_change else None
    is_synonymous = parsed is not None and parsed[0] == parsed[2]
    is_missense = parsed is not None and parsed[0] != parsed[2] and parsed[2] != "*"
    variant_class = "snv" if is_snv else "mnv" if is_mnv else "indel" if is_indel else "other"
    return {
        "variant_class": variant_class,
        "is_snv": is_snv,
        "is_mnv": is_mnv,
        "is_indel": is_indel,
        "is_coding": is_coding,
        "is_missense": is_missense,
        "is_synonymous": is_synonymous,
        "is_nonsense": is_nonsense,
        "is_frameshift": is_frameshift,
        "is_splice_related": False,
        "parsed_aa_change": parsed,
    }


def core_eligibility(flags: dict, parse_failure_reason: str = "") -> tuple[bool, str]:
    if not flags["is_coding"]:
        return False, "noncoding"
    if not flags["is_snv"]:
        return False, "non_snv"
    if flags["is_frameshift"]:
        return False, "frameshift"
    if flags["is_nonsense"]:
        return False, "nonsense"
    if parse_failure_reason:
        return False, parse_failure_reason
    if not flags["is_missense"]:
        return False, "not_missense"
    return True, ""


def plasmy(homo: str, hetero: str) -> str:
    has_homo = str(homo).strip().lower() not in ("", "no", "nr", "nan")
    has_hetero = str(hetero).strip().lower() not in ("", "no", "nr", "nan")
    if has_homo and has_hetero:
        return "homo/hetero"
    if has_homo:
        return "homo"
    if has_hetero:
        return "hetero"
    return ""


def pubmed_count(pubmed_field: str) -> int:
    raw = str(pubmed_field).strip()
    if not raw or raw.lower() == "nan":
        return 0
    return len([x for x in raw.split(",") if x.strip()])


def iter_mitomap_rows(file_path: Path) -> list[dict]:
    import pandas as pd

    df = pd.read_csv(file_path, sep="\t", on_bad_lines="skip", encoding="windows-1252")
    df.columns = [c.strip().lower() for c in df.columns]
    return df.to_dict(orient="records")


THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}


def parse_clinvar_protein(name: str) -> tuple[str, str, str] | None:
    match = re.search(r"p\.([A-Z][a-z]{2}|\*)(\d+)([A-Z][a-z]{2}|\*)", str(name))
    if not match:
        return None
    wt = THREE_TO_ONE.get(match.group(1), "X")
    alt = THREE_TO_ONE.get(match.group(3), "X")
    return wt, match.group(2), alt


def extract_hgvs_c(name: str) -> tuple[str, str, str, str]:
    transcript_match = re.search(r"(NM_[\d.]+|NR_[\d.]+|NP_[\d.]+)", str(name))
    cds_match = re.search(r"c\.\d+([ACGT])>([ACGT])", str(name))
    transcript_id = transcript_match.group(1) if transcript_match else ""
    if cds_match:
        cds_ref = cds_match.group(1)
        cds_alt = cds_match.group(2)
        hgvs_c = f"{transcript_id}:{cds_match.group(0)}" if transcript_id else cds_match.group(0)
    else:
        cds_ref = cds_alt = hgvs_c = ""
    return hgvs_c, transcript_id, cds_ref, cds_alt


def extract_hgvs_p(name: str) -> str:
    match = re.search(r"(p\.[A-Za-z0-9*]+)", str(name))
    return match.group(1) if match else ""


def stars(review: str) -> int:
    review = review.lower()
    if "guideline" in review:
        return 4
    if "expert" in review:
        return 3
    if "multiple" in review and "no conflicts" in review:
        return 2
    if "single submitter" in review or "conflicting" in review:
        return 1
    return 0


def iter_clinvar_rows(file_path: Path):
    import pandas as pd

    for chunk in pd.read_csv(file_path, sep="\t", compression="gzip", chunksize=50000, low_memory=False):
        for row in chunk.to_dict(orient="records"):
            yield row


def base_mt_record(row: dict, version: str, hgnc_ref: GeneReference, mt_sequences: dict) -> tuple[list[VariantRecord], dict | None]:
    try:
        pos = int(row.get("pos", 0))
    except (TypeError, ValueError):
        return [], {
            "source_db": "MITOMAP",
            "source_record_id": str(row.get("id", "")),
            "raw_source_file": "MITOMAP",
            "exclusion_reason": "invalid_position",
        }

    loci_raw = hgnc_ref.get_mt_locus_by_position(pos)
    if loci_raw == "Non-OXPHOS":
        return [], {
            "source_db": "MITOMAP",
            "source_record_id": str(row.get("id", "")),
            "raw_source_file": "MITOMAP",
            "exclusion_reason": "non_oxphos",
            "rCRS_pos": pos,
        }

    genomic_ref = str(row.get("ref", "")).strip().upper()
    genomic_alt = str(row.get("alt", "")).strip().upper()
    aa_change_raw = str(row.get("aachange", "")).strip()
    flags = variant_flags(genomic_ref, genomic_alt, aa_change_raw)

    parse_status = "parsed"
    parse_failure_reason = ""
    if not flags["is_coding"]:
        parse_status = "parsed_noncore"
        parse_failure_reason = "noncoding"
    elif flags["is_frameshift"]:
        parse_status = "parsed_noncore"
        parse_failure_reason = "frameshift"
    elif flags["is_nonsense"]:
        parse_status = "parsed_noncore"
        parse_failure_reason = "nonsense"
    elif flags["parsed_aa_change"] is None and aa_change_raw:
        parse_status = "parsed_noncore"
        parse_failure_reason = "unparsable_aa"

    overlap_genes = loci_raw.split("/")
    is_overlap = len(overlap_genes) > 1
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    coding_ref = complement.get(genomic_ref, genomic_ref) if 14149 <= pos <= 14673 else genomic_ref
    coding_alt = complement.get(genomic_alt, genomic_alt) if 14149 <= pos <= 14673 else genomic_alt
    shared_group_id = f"MITOMAP:{row.get('id', '') or f'm.{pos}{genomic_ref}>{genomic_alt}'}"

    outputs: list[VariantRecord] = []
    for gene in overlap_genes if is_overlap else [overlap_genes[0]]:
        frame_data = {
            "coordinate_status": "source_derived",
            "codon_index": None,
            "ref_aa": flags["parsed_aa_change"][0] if flags["parsed_aa_change"] else "",
            "alt_aa": flags["parsed_aa_change"][2] if flags["parsed_aa_change"] else "",
            "hgvs_p": f"p.{aa_change_raw}" if aa_change_raw else "",
            "aa_change": aa_change_raw,
            "is_synonymous": flags["is_synonymous"],
            "is_missense": flags["is_missense"],
        }
        if is_overlap and gene in mt_sequences and coding_alt:
            frame_data = mt_sequences[gene].consequence_for_variant(pos, coding_alt)

        eligible, exclusion_reason = core_eligibility(
            {
                **flags,
                "is_missense": bool(frame_data.get("is_missense")),
                "is_synonymous": bool(frame_data.get("is_synonymous")),
            },
            "" if frame_data["coordinate_status"] == "resolved" or not is_overlap else "unresolved_overlap_translation",
        )

        record = VariantRecord(
            ann_id=f"m.{pos}{genomic_ref}>{genomic_alt}",
            variant_id=f"{shared_group_id}:{gene}",
            source_variant_group_id=shared_group_id,
            raw_record_id=str(row.get("id", "")),
            source_db="MITOMAP",
            source_db_version=version,
            source_record_id=str(row.get("id", "")),
            genome="mtDNA",
            reference_assembly="rCRS",
            locus=gene,
            loci_raw=loci_raw,
            interpreted_gene=gene,
            nc_change=f"m.{pos}{genomic_ref}>{genomic_alt}",
            aa_change=frame_data["aa_change"],
            ref_nt=coding_ref,
            alt_nt=coding_alt,
            ref_aa=frame_data["ref_aa"],
            alt_aa=frame_data["alt_aa"],
            genomic_pos=pos,
            is_synonymous=bool(frame_data["is_synonymous"]),
            is_overlap=is_overlap,
            overlap_group=loci_raw if is_overlap else None,
            overlap_genes=overlap_genes if is_overlap else [],
            overlap_role="primary_frame_specific_row" if is_overlap else None,
            derived_from_overlap_duplication=is_overlap,
            shared_source_variant_group_id=shared_group_id if is_overlap else None,
            frame_specific_hgvs_p=frame_data["hgvs_p"],
            frame_specific_ref_aa=frame_data["ref_aa"],
            frame_specific_alt_aa=frame_data["alt_aa"],
            frame_specific_is_synonymous=frame_data["is_synonymous"],
            frame_specific_is_missense=frame_data["is_missense"],
            frame_specific_codon_index=frame_data["codon_index"],
            frame_specific_coordinate_status=frame_data["coordinate_status"],
            coding_ref_nt=coding_ref,
            coding_alt_nt=coding_alt,
            genomic_ref_nt=genomic_ref,
            genomic_alt_nt=genomic_alt,
            hgvs_c=f"m.{pos}{genomic_ref}>{genomic_alt}",
            hgvs_p=frame_data["hgvs_p"],
            coordinate_resolution_method="mt_overlap_frame_translation" if is_overlap else "source",
            coordinate_resolution_status=frame_data["coordinate_status"],
            reference_allele_representation="coding_strand_for_nd6_else_genomic",
            variant_class=flags["variant_class"],
            parse_status=parse_status,
            parse_failure_reason=parse_failure_reason,
            curation_status="curated",
            curation_notes="overlap_duplicated" if is_overlap else "",
            eligible_core_comparative_pipeline=eligible,
            core_pipeline_exclusion_reason=exclusion_reason,
            disease=str(row.get("disease", "")).strip(),
            disease_terms=dedupe_terms(str(row.get("disease", "")).strip()),
            clinical_status=str(row.get("status", "")).strip(),
            mitomap_plasmy=plasmy(row.get("homoplasmy", ""), row.get("heteroplasmy", "")),
            mitomap_pubmed_count=pubmed_count(row.get("pubmed_ids", "")),
        )
        record.populate_gene_context()
        record.populate_substitution_properties()
        outputs.append(record)
    return outputs, None


def curate_clinvar_row(row: dict, version: str, hgnc_ref: GeneReference) -> tuple[VariantRecord | None, dict | None]:
    gene_raw = str(row.get("GeneSymbol", "")).strip()
    assembly = str(row.get("Assembly", "")).strip()
    var_type = str(row.get("Type", "")).strip().lower()
    chromosome = str(row.get("Chromosome", "")).strip()
    name = str(row.get("Name", "")).strip()
    hgvs_c, transcript_id, cds_ref, cds_alt = extract_hgvs_c(name)
    hgvs_p = extract_hgvs_p(name)
    protein = parse_clinvar_protein(name)

    if gene_raw not in hgnc_ref.lookup:
        return None, {
            "source_db": "ClinVar",
            "source_record_id": str(row.get("VariationID", "")),
            "raw_source_file": "ClinVar",
            "exclusion_reason": "non_oxphos",
            "gene": gene_raw,
        }

    gene = hgnc_ref.lookup[gene_raw]["symbol"]
    genomic_ref = str(row.get("ReferenceAlleleVCF", "")).strip().upper()
    genomic_alt = str(row.get("AlternateAlleleVCF", "")).strip().upper()
    flags = variant_flags(cds_ref or genomic_ref, cds_alt or genomic_alt, "")
    flags["is_snv"] = var_type == "single nucleotide variant"
    flags["is_mnv"] = "multi" in var_type
    flags["is_indel"] = "deletion" in var_type or "insertion" in var_type
    flags["variant_class"] = "snv" if flags["is_snv"] else "non_snv"
    flags["is_splice_related"] = "splice" in name.lower()

    parse_status = "parsed"
    parse_failure_reason = ""
    aa_change = ""
    ref_aa = ""
    alt_aa = ""
    is_syn = False
    is_missense = False
    if protein is None:
        parse_status = "parsed_noncore"
        parse_failure_reason = "unparsable_aa"
    else:
        ref_aa, aa_pos, alt_aa = protein
        aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
        is_syn = ref_aa == alt_aa
        is_missense = ref_aa != alt_aa and alt_aa != "*" and ref_aa != "X" and alt_aa != "X"

    eligible, exclusion_reason = core_eligibility(
        {
            **flags,
            "is_coding": protein is not None,
            "is_missense": is_missense,
            "is_synonymous": is_syn,
            "is_nonsense": alt_aa == "*" if alt_aa else False,
            "is_frameshift": False,
        },
        parse_failure_reason if assembly == "GRCh38" and chromosome != "MT" else "",
    )
    if assembly != "GRCh38":
        eligible = False
        exclusion_reason = "non_grch38"
    elif chromosome == "MT":
        eligible = False
        exclusion_reason = "mt_in_clinvar_nuclear_branch"
    elif not flags["is_snv"]:
        eligible = False
        exclusion_reason = "non_snv"

    record = VariantRecord(
        ann_id=str(row.get("VariationID", "")),
        variant_id=f"ClinVar:{row.get('VariationID', '')}:{gene}",
        source_variant_group_id=f"ClinVar:{row.get('VariationID', '')}",
        raw_record_id=str(row.get("VariationID", "")),
        source_db="ClinVar",
        source_db_version=version,
        source_record_id=str(row.get("VariationID", "")),
        genome="nucDNA",
        reference_assembly=assembly or "GRCh38",
        locus=gene,
        loci_raw=gene_raw,
        interpreted_gene=gene,
        nc_change=name,
        aa_change=aa_change,
        ref_nt=cds_ref or genomic_ref,
        alt_nt=cds_alt or genomic_alt,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        genomic_pos=int(row.get("Start")) if row.get("Start") not in (None, "") else 0,
        is_synonymous=is_syn,
        coding_ref_nt=cds_ref or genomic_ref,
        coding_alt_nt=cds_alt or genomic_alt,
        genomic_ref_nt=genomic_ref,
        genomic_alt_nt=genomic_alt,
        hgvs_c=hgvs_c,
        hgvs_p=hgvs_p,
        transcript_id=transcript_id,
        coordinate_resolution_method="source",
        coordinate_resolution_status="source",
        reference_allele_representation="cds_if_present_else_genomic",
        variant_class=flags["variant_class"],
        parse_status=parse_status,
        parse_failure_reason=parse_failure_reason,
        curation_status="curated",
        curation_notes="",
        eligible_core_comparative_pipeline=eligible,
        core_pipeline_exclusion_reason=exclusion_reason,
        disease=str(row.get("PhenotypeList", "")).strip(),
        disease_terms=dedupe_terms(str(row.get("PhenotypeList", "")).strip()),
        clinical_status=str(row.get("ClinicalSignificance", "")).strip(),
        clinvar_stars=stars(str(row.get("ReviewStatus", "")).strip()),
        clinvar_review_status=str(row.get("ReviewStatus", "")).strip(),
        clinvar_submitters_n=int(row.get("NumberSubmitters") or 0),
        clinvar_last_evaluated=str(row.get("LastEvaluated", "")).strip(),
        clinvar_conflicting="conflicting" in str(row.get("ClinicalSignificance", "")).lower(),
    )
    record.populate_gene_context()
    record.populate_substitution_properties()
    return record, None


def write_jsonl(records: list[dict], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for row in records:
            handle.write(json.dumps(row) + "\n")


def write_compat_json(records: list[dict]) -> None:
    mt = [row for row in records if row["genome"] == "mtDNA"]
    nuc = [row for row in records if row["genome"] == "nucDNA"]
    with open(COMPAT_MT_JSON, "w", encoding="utf-8") as handle:
        json.dump(mt, handle, indent=2)
    with open(COMPAT_NUC_JSON, "w", encoding="utf-8") as handle:
        json.dump(nuc, handle, indent=2)


def apply_exception_metadata(record: VariantRecord, registry: list[dict]) -> None:
    entry = match_exception_entry(
        registry,
        gene=record.interpreted_gene or record.locus,
        variant_id=record.variant_id or "",
    )
    if not entry:
        return

    record.exception_scope = entry.get("scope") or None
    record.exception_class = entry.get("exception_class") or None
    record.exception_code = entry.get("exception_code") or None
    record.exception_decision = entry.get("decision") or None
    record.manual_review_status = entry.get("manual_review_status") or None
    record.replacement_nm = entry.get("replacement_nm") or None
    record.replacement_enst = entry.get("replacement_enst") or None
    record.rescue_method = entry.get("rescue_method") or None
    record.exception_rationale = entry.get("rationale") or None
    record.exception_evidence_source = entry.get("evidence_source") or None
    record.exception_reviewed_by = entry.get("reviewed_by") or None
    record.exception_review_date = entry.get("review_date") or None
    record.exception_notes = entry.get("notes") or None


def allele_consistency(record: dict) -> tuple[bool | None, bool | None]:
    genomic_ref = str(record.get("genomic_ref_nt") or "").strip().upper()
    genomic_alt = str(record.get("genomic_alt_nt") or "").strip().upper()
    coding_ref = str(record.get("ref_nt") or "").strip().upper()
    coding_alt = str(record.get("alt_nt") or "").strip().upper()

    ref_equal = None if not genomic_ref or not coding_ref else genomic_ref == coding_ref
    alt_equal = None if not genomic_alt or not coding_alt else genomic_alt == coding_alt
    return ref_equal, alt_equal


def cross_source_match_key(record: dict) -> str | None:
    if str(record.get("encoded_by") or "") != "mitochondrial":
        return None

    gene = str(record.get("interpreted_gene") or record.get("locus") or "").strip()
    pos = record.get("genomic_pos")
    genomic_ref = str(record.get("genomic_ref_nt") or "").strip().upper()
    genomic_alt = str(record.get("genomic_alt_nt") or "").strip().upper()
    if not gene or pos in (None, "") or not genomic_ref or not genomic_alt:
        return None
    return f"{gene}|{int(pos)}|{genomic_ref}|{genomic_alt}"


def attach_cross_source_duplicate_linkage(
    curated_records: list[dict],
) -> tuple[list[dict], list[dict], list[dict], dict]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    present_rows: list[dict] = []
    absent_rows: list[dict] = []

    for record in curated_records:
        ref_equal, alt_equal = allele_consistency(record)
        record["genomic_equals_ref_nt"] = ref_equal
        record["genomic_equals_alt_nt"] = alt_equal
        record["cross_source_match_key_nt"] = cross_source_match_key(record)
        if record["cross_source_match_key_nt"]:
            grouped[record["cross_source_match_key_nt"]].append(record)

    cross_source_groups: list[dict] = []
    shared_keys: set[str] = set()
    for key, members in grouped.items():
        sources = sorted({str(item.get("source_db") or "") for item in members if item.get("source_db")})
        genes = sorted({str(item.get("interpreted_gene") or item.get("locus") or "") for item in members})
        variant_ids = sorted({str(item.get("variant_id") or "") for item in members if item.get("variant_id")})
        partner_sources = ",".join(sources)
        partner_variant_ids = sorted(variant_ids)

        gene, pos, ref_nt, alt_nt = key.split("|")
        group_id = f"XSDUP:MT:{gene}:{pos}:{ref_nt}:{alt_nt}"
        is_shared = "ClinVar" in sources and "MITOMAP" in sources
        if is_shared:
            shared_keys.add(key)

        cross_source_groups.append({
            "cross_source_duplicate_group_id": group_id,
            "cross_source_match_key_nt": key,
            "interpreted_gene": gene,
            "genomic_pos": int(pos),
            "genomic_ref_nt": ref_nt,
            "genomic_alt_nt": alt_nt,
            "sources_present": sources,
            "n_rows": len(members),
            "n_variant_ids": len(partner_variant_ids),
            "n_source_dbs": len(sources),
            "is_shared_mt_clinvar_mitomap": is_shared,
            "partner_variant_ids": partner_variant_ids,
        })

        for record in members:
            source_db = str(record.get("source_db") or "")
            if is_shared:
                status = "shared_mt_clinvar_mitomap"
            elif source_db == "ClinVar":
                status = "clinvar_only_mt"
            elif source_db == "MITOMAP":
                status = "mitomap_only_mt"
            else:
                status = "not_applicable"

            record["cross_source_duplicate_group_id"] = group_id
            record["cross_source_duplicate_status"] = status
            record["matched_in_clinvar"] = "ClinVar" in sources
            record["matched_in_mitomap"] = "MITOMAP" in sources
            record["cross_source_partner_count"] = len(partner_variant_ids)
            record["cross_source_partner_variant_ids"] = partner_variant_ids
            record["cross_source_partner_sources"] = sources

            if source_db == "ClinVar":
                out_row = {
                    "variant_id": record.get("variant_id"),
                    "raw_record_id": record.get("raw_record_id"),
                    "source_record_id": record.get("source_record_id"),
                    "interpreted_gene": record.get("interpreted_gene"),
                    "genomic_pos": record.get("genomic_pos"),
                    "genomic_ref_nt": record.get("genomic_ref_nt"),
                    "genomic_alt_nt": record.get("genomic_alt_nt"),
                    "hgvs_c": record.get("hgvs_c"),
                    "hgvs_p": record.get("hgvs_p"),
                    "clinical_status": record.get("clinical_status"),
                    "core_pipeline_exclusion_reason": record.get("core_pipeline_exclusion_reason"),
                    "cross_source_match_key_nt": key,
                    "cross_source_duplicate_group_id": group_id,
                }
                if is_shared:
                    present_rows.append(out_row)
                else:
                    absent_rows.append(out_row)

    for record in curated_records:
        if record.get("cross_source_match_key_nt"):
            continue
        record["cross_source_duplicate_status"] = "not_applicable"
        record["matched_in_clinvar"] = bool(record.get("source_db") == "ClinVar")
        record["matched_in_mitomap"] = bool(record.get("source_db") == "MITOMAP")
        record["cross_source_partner_count"] = 0
        record["cross_source_partner_variant_ids"] = []
        record["cross_source_partner_sources"] = []

    summary = {
        "clinvar_mt_rows_total": int(sum(
            1 for row in curated_records
            if row.get("source_db") == "ClinVar" and row.get("encoded_by") == "mitochondrial"
        )),
        "clinvar_mt_unique_variant_ids_total": int(len({
            row.get("variant_id") for row in curated_records
            if row.get("source_db") == "ClinVar" and row.get("encoded_by") == "mitochondrial"
        })),
        "clinvar_mt_unique_raw_record_ids_total": int(len({
            row.get("raw_record_id") for row in curated_records
            if row.get("source_db") == "ClinVar" and row.get("encoded_by") == "mitochondrial"
        })),
        "mitomap_rows_total": int(sum(1 for row in curated_records if row.get("source_db") == "MITOMAP")),
        "mitomap_unique_variant_ids_total": int(len({
            row.get("variant_id") for row in curated_records if row.get("source_db") == "MITOMAP"
        })),
        "shared_unique_nt_events": int(len(shared_keys)),
        "clinvar_mt_rows_present_in_mitomap": int(len(present_rows)),
        "clinvar_mt_unique_variant_ids_present_in_mitomap": int(len({row["variant_id"] for row in present_rows})),
        "clinvar_mt_unique_raw_record_ids_present_in_mitomap": int(len({row["raw_record_id"] for row in present_rows})),
        "clinvar_mt_rows_absent_from_mitomap": int(len(absent_rows)),
        "clinvar_mt_unique_variant_ids_absent_from_mitomap": int(len({row["variant_id"] for row in absent_rows})),
        "clinvar_mt_unique_raw_record_ids_absent_from_mitomap": int(len({row["raw_record_id"] for row in absent_rows})),
    }
    return cross_source_groups, present_rows, absent_rows, summary


def main() -> None:
    import pandas as pd

    ensure_layout()
    hgnc_file = latest(RAW_REFERENCE_DIR, "Canonical_OXPHOS_Subunits_HGNC_*.csv")
    mitomap_file = latest(RAW_ANNOTATIONS_DIR, "MITOMAP_CodingVariants_*.tsv")
    clinvar_file = latest(RAW_ANNOTATIONS_DIR, "ClinVar_VariantSummary_*.txt.gz")
    mt_coords_file = DERIVED_REFERENCE_DIR / "mtdna_gene_coordinates.tsv"

    hgnc_ref = GeneReference(hgnc_file)
    hgnc_ref.load_coordinates(mt_coords_file)
    mt_coords = load_mtdna_gene_coords(mt_coords_file)
    mt_sequences = load_human_mt_sequences(MT_CODON_DIR, mt_coords)
    exception_registry = load_exception_registry(EXCEPTION_REGISTRY)

    curated_records: list[dict] = []
    excluded_rows: list[dict] = []

    for row in iter_mitomap_rows(mitomap_file):
        records, excluded = base_mt_record(row, extract_date(mitomap_file), hgnc_ref, mt_sequences)
        for record in records:
            apply_exception_metadata(record, exception_registry)
            curated_records.append(asdict(record))
        if excluded:
            excluded_rows.append(excluded)

    for row in iter_clinvar_rows(clinvar_file):
        record, excluded = curate_clinvar_row(row, extract_date(clinvar_file), hgnc_ref)
        if record is not None:
            apply_exception_metadata(record, exception_registry)
            curated_records.append(asdict(record))
        if excluded:
            excluded_rows.append(excluded)

    cross_source_groups, present_rows, absent_rows, overlap_summary = attach_cross_source_duplicate_linkage(curated_records)

    curated_df = pd.DataFrame(curated_records)
    excluded_df = pd.DataFrame(excluded_rows)
    cross_source_df = pd.DataFrame(cross_source_groups)
    curated_df.to_parquet(MASTER_CURATED_PARQUET, index=False)
    if not excluded_df.empty:
        excluded_df.to_parquet(EXCLUDED_INVENTORY_PARQUET, index=False)
    if not cross_source_df.empty:
        cross_source_df.to_csv(CURATION_CROSS_SOURCE_GROUPS_TSV, sep="\t", index=False)
    write_jsonl(curated_records, MASTER_CURATED_JSONL)
    write_compat_json(curated_records)
    present_df = pd.DataFrame(present_rows)
    absent_df = pd.DataFrame(absent_rows)
    if not present_df.empty:
        present_df = present_df.sort_values(["interpreted_gene", "genomic_pos", "variant_id"])
    if not absent_df.empty:
        absent_df = absent_df.sort_values(["interpreted_gene", "genomic_pos", "variant_id"])
    present_df.to_csv(CLINVAR_MT_PRESENT_TSV, sep="\t", index=False)
    absent_df.to_csv(CLINVAR_MT_ABSENT_TSV, sep="\t", index=False)
    with open(CLINVAR_MITOMAP_OVERLAP_SUMMARY_JSON, "w", encoding="utf-8") as handle:
        json.dump(overlap_summary, handle, indent=2)

    metadata = {
        "master_curated_rows": int(len(curated_df)),
        "excluded_inventory_rows": int(len(excluded_df)),
        "mt_overlap_rows": int(sum(1 for row in curated_records if row.get("is_overlap"))),
        "cross_source_duplicate_groups": int(len(cross_source_groups)),
        "cross_source_overlap_summary": overlap_summary,
        "canonical_output": repo_rel(MASTER_CURATED_PARQUET),
        "compatibility_outputs": [repo_rel(COMPAT_MT_JSON), repo_rel(COMPAT_NUC_JSON)],
        "cross_source_outputs": [
            repo_rel(CURATION_CROSS_SOURCE_GROUPS_TSV),
            repo_rel(CLINVAR_MITOMAP_OVERLAP_SUMMARY_JSON),
            repo_rel(CLINVAR_MT_PRESENT_TSV),
            repo_rel(CLINVAR_MT_ABSENT_TSV),
        ],
        "source_files": {
            "hgnc": repo_rel(hgnc_file),
            "mitomap": repo_rel(mitomap_file),
            "clinvar": repo_rel(clinvar_file),
        },
        "exception_registry": repo_rel(EXCEPTION_REGISTRY) if EXCEPTION_REGISTRY.exists() else None,
    }
    with open(CURATION_RUN_METADATA_JSON, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    print(f"Wrote curated master table: {MASTER_CURATED_PARQUET}")
    print(f"Wrote curated JSONL view: {MASTER_CURATED_JSONL}")
    print(f"Wrote curation metadata: {CURATION_RUN_METADATA_JSON}")


if __name__ == "__main__":
    main()
