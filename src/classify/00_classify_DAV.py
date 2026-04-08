import csv
import json
import logging
import re
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
from utils.alignment_parser import AlignmentParser

# ==== Configuration ====
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
CURATED_DIR = DATA_DIR / "annotations" / "curated"
LOG_DIR = DATA_DIR / "logs"

# Input: v2 curated records from 02_curate_variants.py
MT_CURATED = CURATED_DIR / "mtDNA_annotations_v2.json"
NUC_CURATED = CURATED_DIR / "nucDNA_annotations_v2.json"

# Alignment Directories
TOGA_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_aa"
TOGA_NT_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"
MT_AA_DIR = DATA_DIR / "alignments" / "mtdna_aa"
MT_NT_DIR = DATA_DIR / "alignments" / "mtdna_codon"

# Output
OUT_JSON_MT = CURATED_DIR / "cdav_classifications_mtDNA.json"
OUT_JSON_NUC = CURATED_DIR / "cdav_classifications_nucDNA.json"
OUT_MISMATCH_JSON = LOG_DIR / "cdav_mismatches.json"
OUT_LOG_FILE = LOG_DIR / "cdav_classification.log"

MT_COORD_FILE    = DATA_DIR / "reference" / "mtdna_gene_coordinates.tsv"
TX_MAP_FILE      = DATA_DIR / "reference" / "transcript_position_maps.json"
GENOMIC_MAP_FILE = DATA_DIR / "reference" / "genomic_coordinate_maps.json"

_MT_ALIAS = {
    "MT-COX1": "MT-CO1",
    "MT-COX2": "MT-CO2",
    "MT-COX3": "MT-CO3",
    "MT-CYTB": "MT-CYB",
}

# ── Mismatch type constants ───────────────────────────────────────────────────

# Used as the "type" field in each mismatch record and in the per-type summary.
REF_ALLELE_MISMATCH      = "REF_ALLELE_MISMATCH"
TRANSCRIPT_MISMATCH      = "TRANSCRIPT_MISMATCH"
ANCHOR_NOT_FOUND         = "ANCHOR_NOT_FOUND"
POSITION_NOT_IN_ENST     = "POSITION_NOT_IN_ENST"
GENOMIC_POS_NOT_IN_ENST  = "GENOMIC_POS_NOT_IN_ENST"
CODON_EXTRACTION_FAILURE = "CODON_EXTRACTION_FAILURE"
NO_ALIGNMENT             = "NO_ALIGNMENT"
COORD_PARSE_FAILURE      = "COORD_PARSE_FAILURE"

_MISMATCH_DESCRIPTIONS = {
    REF_ALLELE_MISMATCH:
        "ref_nt in curated variant does not match the human CDS base at that "
        "position in the alignment — likely isoform sequence divergence at this site",
    TRANSCRIPT_MISMATCH:
        "mutant codon translates to a different AA than annotated; the TOGA "
        "ENST and the ClinVar NM_ have divergent sequence at this position "
        "(variant is unclassifiable via TOGA)",
    ANCHOR_NOT_FOUND:
        "wildtype AA not found within ±10 residues of the annotated position in "
        "the alignment — possible isoform offset beyond the search window",
    POSITION_NOT_IN_ENST:
        "NM_ residue at this position maps to a gap in the TOGA ENST isoform; "
        "the residue exists only in the NM_ transcript and has no ENST counterpart",
    GENOMIC_POS_NOT_IN_ENST:
        "GRCh38 genomic position is not within any CDS exon of the TOGA ENST; "
        "the variant falls in an intron or UTR of the TOGA transcript",
    CODON_EXTRACTION_FAILURE:
        "CDS position absent from the alignment nucleotide map; codon could not "
        "be constructed (alignment may be truncated)",
    NO_ALIGNMENT:
        "AA and/or NT alignment FASTA files not found for this locus",
    COORD_PARSE_FAILURE:
        "could not parse aa_pos / wt_aa / mut_aa from the aa_change string",
}


# ── Logging setup ─────────────────────────────────────────────────────────────

def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("cdav_classify")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    # Console: INFO and above (progress + warnings)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    # File: DEBUG and above (every mismatch detail)
    fh = logging.FileHandler(OUT_LOG_FILE, mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# ── Mismatch record builder ───────────────────────────────────────────────────

def _mismatch(
    mtype: str, genome: str, locus: str, var: dict, detail: str
) -> dict:
    """Returns a structured mismatch record suitable for JSON output."""
    return {
        "type":       mtype,
        "genome":     genome,
        "locus":      locus,
        "variant_id": var.get("ann_id", ""),
        "tier":       var.get("tier", ""),
        "aa_change":  var.get("aa_change", ""),
        "detail":     detail,
    }


# ── Gene-coordinate helpers ───────────────────────────────────────────────────

def _load_mt_coords() -> dict:
    coords = {}
    if not MT_COORD_FILE.exists():
        return coords
    with open(MT_COORD_FILE) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            start, end = int(row["start"]), int(row["end"])
            entry = (min(start, end), max(start, end), row["strand"])
            coords[row["gene"]] = entry
            if row["gene"] in _MT_ALIAS:
                coords[_MT_ALIAS[row["gene"]]] = entry
    return coords


MT_COORDS = _load_mt_coords()


def _genomic_to_cds(genomic_pos: int, locus: str) -> int | None:
    if locus not in MT_COORDS:
        return None
    start, end, strand = MT_COORDS[locus]
    return (genomic_pos - start + 1) if strand == "+" else (end - genomic_pos + 1)


# ── Variant coordinate parsing ────────────────────────────────────────────────

def parse_variant_coordinates(variant: dict) -> tuple:
    """Returns (aa_pos, wt_aa, mut_aa, nt_pos) or (None, None, None, None)."""
    aa_str   = variant.get("aa_change", "")
    aa_match = re.search(r"([a-zA-Z]+)(\d+)([a-zA-Z]+)", aa_str)
    if not aa_match:
        return None, None, None, None

    wt_aa  = aa_match.group(1)
    aa_pos = int(aa_match.group(2))
    mut_aa = aa_match.group(3)

    if variant["genome"] == "nucDNA":
        nc_str   = variant.get("hgvs_c") or variant.get("nc_change", "")
        nc_match = re.search(r"c\.(\d+)", nc_str)
        if nc_match:
            return aa_pos, wt_aa, mut_aa, int(nc_match.group(1))
    else:
        genomic_pos = variant.get("genomic_pos")
        locus = variant.get("locus", "").split("/")[0]
        if genomic_pos and locus:
            cds_pos = _genomic_to_cds(genomic_pos, locus)
            if cds_pos:
                return aa_pos, wt_aa, mut_aa, cds_pos

    return aa_pos, wt_aa, mut_aa, (aa_pos * 3) - 2


# ── Alignment helpers ─────────────────────────────────────────────────────────

def get_alignment_paths(locus: str, genome: str) -> tuple:
    if genome == "nucDNA":
        return (
            TOGA_AA_DIR / f"{locus}_aa_alignment.fasta",
            TOGA_NT_DIR / f"{locus}_codon_alignment.fasta",
        )
    return (
        MT_AA_DIR / f"{locus}_aa_alignment.fasta",
        MT_NT_DIR / f"{locus}_codon_alignment.fasta",
    )


def _check_ref_allele_raw(
    parser: AlignmentParser, nt_pos: int, expected_ref: str
) -> tuple[bool, str]:
    """Checks ref allele at a raw (pre-anchor) NT position.

    Only used for logging pre-correction mismatches.  The authoritative check
    (using anchor-corrected position) is done inside check_compensation.
    Returns (matches: bool, found_base: str).
    """
    if nt_pos not in parser.nt_map:
        return False, "POS_NOT_IN_MAP"
    col_idx = parser.nt_map[nt_pos]
    found = parser.nt_alignment[parser.ref_header][col_idx]
    return found == expected_ref.upper(), found


# ── Core processing ─────────────────���─────────────────────────────────────────

def process_variants(
    variants: list,
    genome: str,
    loaded_alignments: dict,
    mismatches: list,
    logger: logging.Logger,
    tx_maps: dict | None = None,
    genomic_maps: dict | None = None,
) -> tuple:
    """
    Classifies non-discarded, non-synonymous variants as cDAV or uDAV.

    cDAV: mutant AA found in ≥1 non-human species (AA-level);
          strictest subset uses the exact mutant codon (NT-level).
    uDAV: mutant AA absent from all aligned species.
    Undetermined: alignment absent or coordinates unresolvable.

    All mismatches / inconsistencies are appended to `mismatches` and logged.
    Returns (enriched_list, global_stats, tier_stats, ref_mismatch_count).
    """
    enriched = []
    ref_mismatch = 0
    missing_loci: dict[str, int] = {}   # locus → count of affected variants
    tier_stats   = defaultdict(lambda: {"Total": 0, "aa_cDAV": 0, "nt_cDAV": 0})
    global_stats = {"Total": 0, "aa_cDAV": 0, "nt_cDAV": 0}

    total = len(variants)
    for idx, var in enumerate(variants):
        if idx % 500 == 0:
            logger.info("[%s] Processed %d / %d...", genome, idx, total)

        tier  = var.get("tier", "Discarded")
        locus = var.get("locus", "").split("/")[0]

        if "Discarded" in tier or var["is_synonymous"]:
            continue

        # ── Coordinate parsing ────────────────────────────────────────────────
        aa_pos, wt_aa, mut_aa, nt_pos = parse_variant_coordinates(var)
        if not aa_pos:
            detail = f"aa_change='{var.get('aa_change', '')}'"
            logger.warning("[%s] %s | %s | %s | %s",
                           genome, COORD_PARSE_FAILURE, locus,
                           var.get("ann_id", ""), detail)
            mismatches.append(_mismatch(COORD_PARSE_FAILURE, genome, locus, var, detail))
            continue

        # ── Strategy B: genomic coordinate override ───────────────────────────
        # For the 13 genes where TOGA used a non-MANE isoform, look up the
        # GRCh38 genomic position directly in the TOGA exon structure to get
        # the correct TOGA CDS and AA positions.  This bypasses the NM_→ENST
        # position map (which compared the *current* NM_ against TOGA and was
        # wrong because ClinVar used an older NM_ version matching TOGA).
        using_genomic_map = False
        if genome == "nucDNA" and genomic_maps and locus in genomic_maps:
            gpos = var.get("genomic_pos")
            if gpos is not None:
                gmap_entry = genomic_maps[locus].get(int(gpos))
                if gmap_entry is not None:
                    aa_pos = gmap_entry["aa_pos"]
                    nt_pos = gmap_entry["cds_pos"]
                    using_genomic_map = True
                else:
                    # Variant genomic position not in any CDS exon of TOGA ENST
                    detail = (
                        f"locus={locus}  genomic_pos={gpos}  "
                        f"not in TOGA ENST CDS exons (intronic or UTR)"
                    )
                    logger.warning("[%s] %s | %s | %s | %s",
                                   genome, GENOMIC_POS_NOT_IN_ENST, locus,
                                   var.get("ann_id", ""), detail)
                    mismatches.append(
                        _mismatch(GENOMIC_POS_NOT_IN_ENST, genome, locus, var, detail)
                    )
                    var.update({
                        "is_cdav_amino_acid":            None,
                        "is_cdav_nucleotide":            None,
                        "n_species_aligned":             None,
                        "n_species_with_disease_allele": None,
                        "lineages_with_disease_allele":  [],
                        "ref_allele_match":              None,
                        "mismatch_reason":               "Genomic position not in TOGA ENST",
                    })
                    enriched.append(var)
                    global_stats["Total"] += 1
                    tier_stats[tier]["Total"] += 1
                    continue

        # ── Alignment loading ─────────────────────────────────────────────────
        if locus not in loaded_alignments:
            aa_path, nt_path = get_alignment_paths(locus, genome)
            if not aa_path.exists() or not nt_path.exists():
                loaded_alignments[locus] = None
            else:
                # Genomic-map genes: pass tx_pos_map=None — genomic override
                # gives TOGA-space coordinates directly; no NM_→ENST remap needed.
                # Other genes: use the global NM_→ENST position map if available.
                tx_for_locus = (
                    None if (genomic_maps and locus in genomic_maps)
                    else (tx_maps or {}).get(locus)
                )
                loaded_alignments[locus] = AlignmentParser(
                    aa_path, nt_path, genome, tx_pos_map=tx_for_locus,
                )

        parser = loaded_alignments[locus]
        if parser is None:
            missing_loci[locus] = missing_loci.get(locus, 0) + 1
            var.update({
                "is_cdav_amino_acid":        None,
                "is_cdav_nucleotide":        None,
                "n_species_aligned":         None,
                "n_species_with_disease_allele": None,
                "lineages_with_disease_allele":  [],
                "ref_allele_match":          None,
                "mismatch_reason":           "No alignment file",
            })
            enriched.append(var)
            global_stats["Total"] += 1
            tier_stats[tier]["Total"] += 1
            continue

        n_species_aligned = len(parser.aa_alignment) - 1
        ref_nt = var["ref_nt"]
        alt_nt = var["alt_nt"]

        # ── Single call: anchor correction + codon construction + species scan ─
        # Codon is built at the anchor-corrected NT position, not at the raw
        # ClinVar c. coordinate — this eliminates isoform-offset errors where
        # the TOGA ENST and the ClinVar NM_ have different CDS starts.
        comp = parser.check_compensation(aa_pos, wt_aa, mut_aa, nt_pos, alt_nt)
        mut_codon = comp["mut_codon"]

        # ── Anchor-not-found ──────────────────────────────────────────────────
        if not comp["anchor_found"]:
            if comp.get("position_not_in_enst"):
                mtype  = POSITION_NOT_IN_ENST
                detail = (
                    f"locus={locus}  aa_pos={aa_pos}  wt_aa={wt_aa}  "
                    f"aa_change={var.get('aa_change')}  "
                    f"(NM_ residue absent from ENST isoform)"
                )
                mismatch_reason = "Position not in ENST isoform"
            else:
                mtype  = ANCHOR_NOT_FOUND
                detail = (
                    f"locus={locus}  aa_pos={aa_pos}  wt_aa={wt_aa}  "
                    f"aa_change={var.get('aa_change')}"
                )
                mismatch_reason = "Anchor not found"
            logger.warning("[%s] %s | %s | %s | %s",
                           genome, mtype, locus,
                           var.get("ann_id", ""), detail)
            mismatches.append(_mismatch(mtype, genome, locus, var, detail))
            var.update({
                "is_cdav_amino_acid":            None,
                "is_cdav_nucleotide":            None,
                "n_species_aligned":             n_species_aligned,
                "n_species_with_disease_allele": None,
                "lineages_with_disease_allele":  [],
                "ref_allele_match":              None,
                "mismatch_reason":               mismatch_reason,
            })
            enriched.append(var)
            global_stats["Total"] += 1
            tier_stats[tier]["Total"] += 1
            continue

        # ── Ref-allele check (at anchor-corrected position) ───────────────────
        ref_base = comp["ref_base_found"]
        ref_ok = ref_base not in ("POS_NOT_IN_MAP",) and ref_base == ref_nt.upper()
        if not ref_ok:
            ref_mismatch += 1
            var["ref_allele_match"] = False
            var["mismatch_reason"] = (
                f"Expected {ref_nt} at corrected CDS pos "
                f"{comp['corrected_nt_pos']} (raw: {nt_pos}), found {ref_base}"
            )
            detail = (
                f"locus={locus}  aa_change={var.get('aa_change')}  "
                f"raw_cds={nt_pos}  corrected_cds={comp['corrected_nt_pos']}  "
                f"expected={ref_nt}  found={ref_base}"
            )
            logger.warning("[%s] %s | %s | %s | %s",
                           genome, REF_ALLELE_MISMATCH, locus,
                           var.get("ann_id", ""), detail)
            mismatches.append(_mismatch(REF_ALLELE_MISMATCH, genome, locus, var, detail))
        else:
            var["ref_allele_match"] = True
            var["mismatch_reason"]  = None

        # ── Codon extraction failure ──────────────────────────────────────────
        if mut_codon is None:
            detail = (
                f"locus={locus}  corrected_cds={comp['corrected_nt_pos']}  "
                f"alt_nt={alt_nt}  aa_change={var.get('aa_change')}"
            )
            logger.warning("[%s] %s | %s | %s | %s",
                           genome, CODON_EXTRACTION_FAILURE, locus,
                           var.get("ann_id", ""), detail)
            mismatches.append(_mismatch(CODON_EXTRACTION_FAILURE, genome, locus, var, detail))
            var.update({
                "is_cdav_amino_acid":            None,
                "is_cdav_nucleotide":            None,
                "n_species_aligned":             n_species_aligned,
                "n_species_with_disease_allele": None,
                "lineages_with_disease_allele":  [],
            })
            enriched.append(var)
            global_stats["Total"] += 1
            tier_stats[tier]["Total"] += 1
            continue

        # ── Transcript validation (nucDNA only) ───────────────────────────────
        # Uses anchor-corrected codon; genuine splice-isoform mismatches still fire.
        if genome == "nucDNA":
            try:
                translated_aa = str(Seq(mut_codon).translate(table=1))
            except Exception:
                translated_aa = "ERR"

            expected_aa = var.get("alt_aa") or mut_aa
            if translated_aa != expected_aa:
                ref_mismatch += 1
                var["ref_allele_match"] = False
                var["mismatch_reason"] = (
                    f"Codon {mut_codon} → {translated_aa}, "
                    f"expected {expected_aa} (TOGA: {parser.transcript_id})"
                )
                detail = (
                    f"locus={locus}  aa_change={var.get('aa_change')}  "
                    f"raw_cds={nt_pos}  corrected_cds={comp['corrected_nt_pos']}  "
                    f"codon={mut_codon}  translated={translated_aa}  "
                    f"expected={expected_aa}  toga_tx={parser.transcript_id}"
                )
                logger.warning("[%s] %s | %s | %s | %s",
                               genome, TRANSCRIPT_MISMATCH, locus,
                               var.get("ann_id", ""), detail)
                mismatches.append(_mismatch(TRANSCRIPT_MISMATCH, genome, locus, var, detail))
                var.update({
                    "is_cdav_amino_acid":            None,
                    "is_cdav_nucleotide":            None,
                    "n_species_aligned":             n_species_aligned,
                    "n_species_with_disease_allele": None,
                    "lineages_with_disease_allele":  [],
                })
                enriched.append(var)
                global_stats["Total"] += 1
                tier_stats[tier]["Total"] += 1
                continue

        # ── cDAV classification ───────────────────────────────────────────────
        var.update({
            "is_cdav_amino_acid":        comp["aa_cdar"],
            "is_cdav_nucleotide":        comp["nt_cdar"],
            "n_species_aligned":         n_species_aligned,
            "n_species_with_disease_allele": len(comp["aa_species"]),
            "lineages_with_disease_allele":  comp["aa_species"],
        })
        enriched.append(var)

        global_stats["Total"] += 1
        tier_stats[tier]["Total"] += 1
        if comp["aa_cdar"]:
            global_stats["aa_cDAV"] += 1
            tier_stats[tier]["aa_cDAV"] += 1
        if comp["nt_cdar"]:
            global_stats["nt_cDAV"] += 1
            tier_stats[tier]["nt_cDAV"] += 1

    # Log missing loci once per locus (not per variant) to reduce noise
    for missing_locus, count in sorted(missing_loci.items()):
        detail = f"locus={missing_locus}  affected_variants={count}"
        logger.warning("[%s] %s | %s | — | %s", genome, NO_ALIGNMENT, missing_locus, detail)
        # Add one representative record to the mismatch file
        mismatches.append({
            "type":       NO_ALIGNMENT,
            "genome":     genome,
            "locus":      missing_locus,
            "variant_id": f"({count} variants)",
            "tier":       "",
            "aa_change":  "",
            "detail":     detail,
        })

    return enriched, global_stats, dict(tier_stats), ref_mismatch


# ── Output helpers ────────────────────────────────────────────────────────────

def print_mismatch_summary(mismatches: list, logger: logging.Logger) -> None:
    """Logs a per-type breakdown of all mismatches to console and file."""
    by_type: dict[str, list] = defaultdict(list)
    for m in mismatches:
        by_type[m["type"]].append(m)

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
        n = len(by_type.get(mtype, []))
        if n:
            logger.info("  %-28s : %d", mtype, n)
            logger.info("    %s", _MISMATCH_DESCRIPTIONS[mtype])

    # Per-genome breakdown
    for genome in ("mtDNA", "nucDNA"):
        g_items = [m for m in mismatches if m["genome"] == genome]
        if not g_items:
            continue
        logger.info("")
        logger.info("  %s (%d total):", genome, len(g_items))
        g_by_type: dict[str, list] = defaultdict(list)
        for m in g_items:
            g_by_type[m["type"]].append(m)
        for mtype, items in sorted(g_by_type.items(), key=lambda x: -len(x[1])):
            # Show up to 5 example variant IDs inline
            examples = ", ".join(m["variant_id"] for m in items[:5])
            tail = f" … (+{len(items)-5} more)" if len(items) > 5 else ""
            logger.info("    %-28s %3d  e.g. %s%s", mtype, len(items), examples, tail)

    logger.info("=" * 60)
    logger.info("Full details → %s", OUT_MISMATCH_JSON)
    logger.info("Full log     → %s", OUT_LOG_FILE)


def print_genome_summary(
    genome: str,
    enriched: list,
    global_stats: dict,
    tier_stats: dict,
    ref_mismatch: int,
    logger: logging.Logger,
) -> None:
    n_undetermined = sum(1 for v in enriched if v.get("is_cdav_amino_acid") is None)
    classified     = global_stats["Total"] - n_undetermined

    logger.info("")
    logger.info("%s:", genome)
    logger.info("  Analyzed variants        : %d", global_stats["Total"])
    logger.info("  Ref-allele / TX mismatches: %d", ref_mismatch)
    logger.info("  Undetermined (no alignment): %d", n_undetermined)
    if classified > 0:
        logger.info("  AA-level cDAVs           : %d (%.1f%%)",
                    global_stats["aa_cDAV"], global_stats["aa_cDAV"] / classified * 100)
        logger.info("  NT-level cDAVs           : %d (%.1f%%)",
                    global_stats["nt_cDAV"], global_stats["nt_cDAV"] / classified * 100)

    logger.info("")
    logger.info("  %-12s %7s %9s %9s", "Tier", "Total", "aa_cDAV", "nt_cDAV")
    logger.info("  %s", "-" * 40)
    for tier in sorted(tier_stats.keys()):
        s = tier_stats[tier]
        aa_pct = f"({s['aa_cDAV']/s['Total']*100:.0f}%)" if s["Total"] else ""
        nt_pct = f"({s['nt_cDAV']/s['Total']*100:.0f}%)" if s["Total"] else ""
        logger.info("  %-12s %7d %5d %-5s %5d %s",
                    tier, s["Total"], s["aa_cDAV"], aa_pct, s["nt_cDAV"], nt_pct)
    logger.info("-" * 50)


# ── Entry point ───────────────────────────────────────────────────────────────

def _load_genomic_maps(logger: logging.Logger) -> dict:
    """
    Loads genomic_coordinate_maps.json.

    Returns {gene: {grch38_pos (int): {"cds_pos": int, "aa_pos": int}}} for
    the 13 TOGA genes whose ENST differs from the MANE Select transcript.

    These maps let the classify script look up the correct TOGA CDS/AA position
    directly from the ClinVar GRCh38 genomic coordinate, bypassing the NM_→ENST
    position remapping entirely (Strategy B).
    """
    if not GENOMIC_MAP_FILE.exists():
        logger.warning(
            "genomic_coordinate_maps.json not found at %s — "
            "Strategy B unavailable; isoform-mismatch variants will remain "
            "undetermined. Run src/data_prep/00g_build_genomic_coordinate_maps.py.",
            GENOMIC_MAP_FILE,
        )
        return {}

    with open(GENOMIC_MAP_FILE) as f:
        raw = json.load(f)

    genomic_maps: dict[str, dict[int, dict]] = {}
    for gene, entry in raw.items():
        genomic_maps[gene] = {int(k): v for k, v in entry["map"].items()}

    logger.info(
        "Loaded genomic coordinate maps for %d genes (Strategy B).", len(genomic_maps)
    )
    return genomic_maps


def _load_tx_maps(logger: logging.Logger) -> dict:
    """
    Loads transcript_position_maps.json.

    Returns {gene: {nm_aa_pos (int): enst_aa_pos (int|None)}} for genes whose
    type is "mapped".  Identity genes and missing genes return no entry (the
    AlignmentParser fallback anchor search handles them).
    """
    if not TX_MAP_FILE.exists():
        logger.warning(
            "transcript_position_maps.json not found at %s — "
            "falling back to sequence-anchor strategy for all genes. "
            "Run src/data_prep/00f_build_transcript_position_maps.py to generate it.",
            TX_MAP_FILE,
        )
        return {}

    with open(TX_MAP_FILE) as f:
        raw = json.load(f)

    tx_maps: dict[str, dict[int, int | None]] = {}
    n_identity = n_mapped = 0
    for gene, entry in raw.items():
        if entry.get("type") == "mapped":
            tx_maps[gene] = {int(k): v for k, v in entry["map"].items()}
            n_mapped += 1
        else:
            n_identity += 1

    logger.info(
        "Loaded transcript position maps: %d genes remapped, %d identity (no remap needed).",
        n_mapped, n_identity,
    )
    return tx_maps


def main():
    logger = setup_logging()
    logger.info("Initializing cDAV Classification Engine")

    with open(MT_CURATED) as f:
        mt_variants = json.load(f)
    with open(NUC_CURATED) as f:
        nuc_variants = json.load(f)

    logger.info("Loaded %d mtDNA variants, %d nucDNA variants.",
                len(mt_variants), len(nuc_variants))

    tx_maps      = _load_tx_maps(logger)
    genomic_maps = _load_genomic_maps(logger)

    mismatches: list[dict] = []
    loaded_alignments: dict = {}

    logger.info("")
    logger.info("Processing mtDNA variants...")
    mt_enriched, mt_global, mt_tier, mt_mismatch = process_variants(
        mt_variants, "mtDNA", loaded_alignments, mismatches, logger,
        tx_maps=tx_maps, genomic_maps=genomic_maps,
    )

    logger.info("")
    logger.info("Processing nucDNA variants...")
    nuc_enriched, nuc_global, nuc_tier, nuc_mismatch = process_variants(
        nuc_variants, "nucDNA", loaded_alignments, mismatches, logger,
        tx_maps=tx_maps, genomic_maps=genomic_maps,
    )

    # Save classification outputs
    CURATED_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON_MT, "w", encoding="utf-8") as f:
        json.dump(mt_enriched, f, indent=2)
    with open(OUT_JSON_NUC, "w", encoding="utf-8") as f:
        json.dump(nuc_enriched, f, indent=2)

    logger.info("")
    logger.info("Saved %d mtDNA records  → %s", len(mt_enriched), OUT_JSON_MT)
    logger.info("Saved %d nucDNA records → %s", len(nuc_enriched), OUT_JSON_NUC)

    # Save structured mismatch report
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_MISMATCH_JSON, "w", encoding="utf-8") as f:
        json.dump(mismatches, f, indent=2)

    # Print summaries
    logger.info("")
    logger.info("=" * 50)
    logger.info("cDAV CLASSIFICATION SUMMARY")
    logger.info("=" * 50)
    print_genome_summary("mtDNA",  mt_enriched,  mt_global,  mt_tier,  mt_mismatch,  logger)
    print_genome_summary("nucDNA", nuc_enriched, nuc_global, nuc_tier, nuc_mismatch, logger)

    print_mismatch_summary(mismatches, logger)


if __name__ == "__main__":
    main()
