import re
import csv
import json
import zipfile
import io
import pandas as pd
from pathlib import Path
from utils.VariantAnnotation import VariantAnnotation


class MitoVariantAnnotator:
    def __init__(
        self,
        mitimpact_zip: Path,
        phylotree_zip: Path,
        hgnc_reference=None,
        pathogenic_threshold=0.716,
    ):
        self.pathogenic_threshold = pathogenic_threshold
        self.mitimpact_lookup = self._load_mitimpact(mitimpact_zip)
        self.phylotree_markers = self._load_phylotree(phylotree_zip)
        self.hgnc_reference = hgnc_reference  # Inject the GeneReference object

        # Standard rCRS coordinates for MT-OXPHOS genes
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

    def _get_locus(self, pos):
        """Map standard rCRS coordinate to mitochondrial gene, natively allowing overlaps."""
        genes = [name for name, start, end in self.mt_genes if start <= pos <= end]
        return "/".join(genes) if genes else "Non-OXPHOS"

    def _transition(self, nt):
        """Return the natural transition partner of a nucleotide."""
        return {"A": "G", "G": "A", "C": "T", "T": "C"}.get(nt.upper())

    def _load_mitimpact(self, zip_path):
        lookup = {}
        row_count = 0
        duplicate_count = 0
        skipped_count = 0

        with zipfile.ZipFile(zip_path) as zf:
            data_file = next(
                f
                for f in zf.namelist()
                if f.endswith(".txt") or f.endswith(".tsv") or f.endswith(".csv")
            )
            with zf.open(data_file) as fh:
                reader = csv.DictReader(
                    io.TextIOWrapper(fh, encoding="utf-8", errors="ignore"),
                    delimiter="\t",
                )
                for row in reader:
                    row_count += 1
                    try:
                        key = (
                            int(row["Start"]),
                            str(row["Ref"]).strip().upper(),
                            str(row["Alt"]).strip().upper(),
                        )
                        score_str = str(row.get("APOGEE2_score", "")).strip()
                        score = float(score_str) if score_str else 0.0
                        mitoclass = str(row.get("Mitoclass1", "")).strip()
                        gene = str(row.get("Gene_symbol", "")).strip()

                        if key in lookup:
                            duplicate_count += 1
                            existing_score, existing_class, existing_genes = lookup[key]

                            # Merge gene names
                            combined_genes = (
                                f"{existing_genes}/{gene}"
                                if gene not in existing_genes
                                else existing_genes
                            )
                            # Merge Mitoclass if different reading frames yield different predictions
                            combined_class = (
                                f"{existing_class} | {mitoclass}"
                                if mitoclass and mitoclass not in existing_class
                                else existing_class
                            )

                            lookup[key] = (score, combined_class, combined_genes)
                        else:
                            lookup[key] = (score, mitoclass, gene)

                    except (ValueError, KeyError, TypeError):
                        skipped_count += 1

        print(f"MitImpact File Stats: {row_count} total data rows.")
        print(
            f"  -> Deduplicated {duplicate_count} overlapping gene entries by merging them."
        )
        if skipped_count > 0:
            print(f"  -> Skipped {skipped_count} malformed rows.")
        print(f"Loaded {len(lookup)} unique MitImpact loci...")

        return lookup

    def _load_phylotree(self, zip_path):
        """Parse PhyloTree exactly using the HTML table logic."""
        markers = {}
        with zipfile.ZipFile(zip_path) as zf:
            for filename in zf.namelist():
                if "__MACOSX" in filename:
                    continue
                if filename.endswith((".htm", ".html")):
                    with zf.open(filename) as fh:
                        # Use windows-1252 as defined in original PhyloTree export
                        content = fh.read().decode("windows-1252", errors="ignore")

                        # 1. Strip all HTML tags
                        text = re.sub(r"<[^>]+>", " ", content)
                        # 2. Strip uncertainty brackets
                        text = re.sub(r"[()]", " ", text)

                        # 3. Parse tokens
                        for tok in text.split():
                            tok = tok.strip()
                            if not tok or "." in tok or "d" in tok.lower():
                                continue  # Skip indels

                            # Parse optional back-mut (@), pos, optional alt, optional (!)
                            m = re.fullmatch(r"(@?)(\d+)([A-Za-z]?)!?", tok)
                            if m:
                                is_back = m.group(1) == "@"
                                pos = int(m.group(2))
                                explicit_alt = (
                                    m.group(3).upper() if m.group(3) else None
                                )

                                # We don't discard based on back-mutations
                                if not is_back:
                                    if pos not in markers:
                                        markers[pos] = set()
                                    markers[pos].add(explicit_alt)

        print(
            f"Loaded haplogroup definitions for {len(markers)} mitochondrial positions..."
        )
        return markers

    def curate(self, raw_df: pd.DataFrame) -> list[VariantAnnotation]:
        curated_objects = []
        for _, row in raw_df.iterrows():
            try:
                pos = int(row.get("pos", 0))
            except ValueError:
                continue

            ref = str(row.get("ref", "")).strip().upper()
            alt = str(row.get("alt", "")).strip().upper()
            aachange = str(row.get("aachange", "")).strip()

            locus = self._get_locus(pos)

            # 1. Structural Filters
            if locus == "Non-OXPHOS":
                continue

            # 2. HGNC Validation
            # Handle overlapping genes (e.g., MT-ATP8/MT-ATP6) by checking if at least one is a target
            if self.hgnc_reference:
                locus_genes = locus.split("/")
                if not any(self.hgnc_reference.is_target(g) for g in locus_genes):
                    continue

            # 3. Exclude noncoding, indels, and stop codons
            if (
                not aachange
                or "noncoding" in aachange.lower()
                or "syn" in aachange.lower()
                and len(aachange) < 5
            ):
                # Note: MITOMAP sometimes puts literal "syn" in the column instead of the amino acid.
                # We drop those if we cannot parse the actual codon change.
                continue
            if len(ref) > 1 or len(alt) > 1 or "frameshift" in aachange.lower():
                continue
            if "*" in aachange or "Ter" in aachange or "Stop" in aachange:
                continue

            # 4. Strict Parsing for Missense and Synonymous
            match = re.match(r"^([a-zA-Z]+)(\d+)([a-zA-Z]+)$", aachange)
            if not match:
                # Explicitly drop variants that fail the regex parser (garbage strings)
                continue

            # We no longer discard the variant if group(1) == group(3).
            # Synonymous variants are now retained and processed.

            apogee_score, mitoclass, _ = self.mitimpact_lookup.get(
                (pos, ref, alt), (0.0, "", "")
            )
            status = str(row.get("status", "")).strip()

            # --- DYNAMIC HAPLOGROUP EVALUATION ---
            is_haplogroup = False
            if pos in self.phylotree_markers:
                for explicit_alt in self.phylotree_markers[pos]:
                    if explicit_alt == alt:
                        is_haplogroup = True
                        break
                    elif explicit_alt is None:  # Implicit transition
                        if alt == self._transition(ref):
                            is_haplogroup = True
                            break

            # --- UNIFIED 3-TIER SCALE FOR MTDNA ---
            if is_haplogroup or "benign" in status.lower() or "[lb]" in status.lower():
                tier = "Discarded"
            elif "Cfrm [P]" in status or "Cfrm [LP]" in status:
                tier = "Tier 1"
            elif (
                apogee_score >= self.pathogenic_threshold
                and "pathogenic" in str(mitoclass).lower()
            ):
                tier = "Tier 2"
            else:
                tier = "Tier 3"

            curated_objects.append(
                VariantAnnotation(
                    ann_id=f"m.{pos}{ref}>{alt}",
                    locus=locus,
                    nc_change=f"{pos}{ref}>{alt}",
                    aa_change=aachange,
                    disease=str(row.get("disease", "Unknown")).strip(),
                    genome="mtDNA",
                    tier=tier,
                    pathogenic_score=apogee_score,
                )
            )

        return curated_objects


class NucVariantAnnotator:
    # Standard ClinVar rating guidelines (evaluates evidence quality)
    STAR_MAP = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, single submitter": 1,
        "criteria provided, conflicting classifications": 1,
        "no assertion criteria provided": 0,
        "no assertion provided": 0,
        "no interpretation for the single variant": 0,
    }

    # Standard amino acid mapping
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

    def __init__(self, myvariant_json: Path, min_clinvar_stars=0):
        self.min_stars = min_clinvar_stars
        self.annotation_lookup = self._load_myvariant(myvariant_json)

    def _load_myvariant(self, json_path):
        """Extracts REVEL (dbNSFP) and gnomAD frequencies safely from MyVariant data."""
        lookup = {}
        try:
            with open(json_path, "r") as f:
                data = json.load(f)
                hits = data.get("hits", data) if isinstance(data, dict) else data

                for item in hits:
                    revel = 0.0
                    af = 0.0

                    # 1. Pathogenicity (REVEL from dbNSFP)
                    dbnsfp = item.get("dbnsfp", {})
                    if dbnsfp:
                        r_score = dbnsfp.get("revel", {}).get("score", 0.0)
                        if isinstance(r_score, list):
                            r_score = r_score[0]
                        try:
                            revel = float(r_score)
                        except Exception:
                            pass

                    # 2. Population Frequency (gnomAD)
                    gnomad = item.get("gnomad_exome", item.get("gnomad_genome", {}))
                    if gnomad:
                        af_val = gnomad.get("af", 0.0)
                        if isinstance(af_val, list):
                            af_val = af_val[0]
                        try:
                            af = float(af_val)
                        except Exception:
                            pass

                    # 3. Establish links
                    clinvar = item.get("clinvar", {})
                    allele_id = str(clinvar.get("allele_id", ""))
                    query_id = str(item.get("query", ""))

                    metrics = {"revel": revel, "af": af}
                    if allele_id:
                        lookup[allele_id] = metrics
                    if query_id:
                        lookup[query_id] = metrics
        except Exception:
            pass
        return lookup

    def _parse_protein(self, name_field):
        m = re.search(
            r"p\.([A-Z][a-z]{0,2}|\*)(\d+)([A-Z][a-z]{0,2}|\*)", str(name_field)
        )
        if not m:
            return ""

        wt_raw = m.group(1)
        mut_raw = m.group(3)

        # Strict Missense Filter: Reject Stop-gains / Nonsense AND Stop-Losses
        stop_flags = ["Ter", "*", "X"]
        if (
            mut_raw in stop_flags
            or wt_raw in stop_flags
            or "Ter" in mut_raw
            or "Ter" in wt_raw
        ):
            return ""

        wt = self.three_to_one.get(wt_raw, wt_raw)
        mut = self.three_to_one.get(mut_raw, mut_raw)

        # Reject Synonymous variants
        if wt == mut:
            return ""

        return f"{wt}{m.group(2)}{mut}"

    def curate(self, raw_df: pd.DataFrame) -> list[VariantAnnotation]:
        curated_objects = []
        for _, row in raw_df.iterrows():
            rev_stat = str(row.get("ReviewStatus", "")).lower()
            stars = self.STAR_MAP.get(rev_stat, 0)

            # Evidence Threshold
            if stars < self.min_stars:
                continue

            name = str(row.get("Name", "")).strip()
            aa_change = self._parse_protein(name)
            if not aa_change:
                continue

            nc_match = re.search(r"(c\.[0-9]+[A-Z]>[A-Z])", name)
            nc_change = nc_match.group(1) if nc_match else name
            gene = str(row.get("GeneSymbol", "")).strip()
            allele_id_str = str(row.get("AlleleID", name))

            # Fetch metrics from JSON
            metrics = self.annotation_lookup.get(
                allele_id_str, {"revel": 0.0, "af": 0.0}
            )
            gnomad_af = metrics["af"]
            revel_score = metrics["revel"]

            # --- UNIFIED 3-TIER SCALE FOR NUCDNA WITH POPULATION SIEVE ---
            if gnomad_af > 0.01:
                tier = "Discarded"  # Fails the gnomAD Mendelian population screen
            elif stars >= 2:
                tier = "Tier 1"
            elif stars == 1:
                tier = "Tier 2"
            else:
                tier = "Tier 3"

            curated_objects.append(
                VariantAnnotation(
                    ann_id=allele_id_str,
                    locus=gene,
                    nc_change=nc_change,
                    aa_change=aa_change,
                    disease=str(row.get("PhenotypeList", "Unknown")).strip(),
                    genome="nucDNA",
                    tier=tier,
                    pathogenic_score=revel_score,
                )
            )
        return curated_objects
