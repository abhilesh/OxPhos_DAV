# Resolved Nuances and Fixes in the Current Pipeline

This note captures issues that appear to have been addressed in the current codebase, based on the implementations in `src/data_prep`, `src/align`, `src/classify`, `src/structural`, and `src/utils`.

It is intentionally separate from unresolved risks and review findings.

## 1. Historical gene-name aliases are now handled explicitly

Problem:
Older source records and external resources may use outdated symbols such as `COXFA4`, while the rest of the pipeline uses the current HGNC symbol `NDUFA4`.

How this was fixed:
- `GeneReference` now indexes both approved symbols and previous symbols from the HGNC table.
- TOGA filenames are also remapped to canonical HGNC symbols in downstream scripts.

Where:
- `src/utils/parsers.py`
- `src/structural/00_map_davs_to_structure.py`
- `src/structural/01_find_compensating_partners.py`

Why this matters:
- Prevents silent record loss during parsing.
- Prevents alignments and structural mappings from being skipped just because the data source used an older symbol.

## 2. MT-ND6 minus-strand handling has been corrected at parse time

Problem:
MT-ND6 is encoded on the mitochondrial minus strand, so raw genomic alleles from MITOMAP cannot be used directly as coding-strand alleles.

How this was fixed:
- The MITOMAP parser complements `ref` and `alt` for positions in the MT-ND6 interval while preserving the original genomic alleles separately.

Where:
- `src/utils/parsers.py`

Why this matters:
- Prevents codon injection and amino-acid interpretation from being wrong for all ND6 variants.
- Preserves both coding-strand and genomic-strand representations for later use.

## 3. CDS-relative and genomic alleles are now kept separate

Problem:
For both mtDNA and nucDNA, downstream steps need to distinguish between coding-sequence alleles and raw genomic alleles.

How this was fixed:
- Parsers now retain `ref` / `alt` alongside `genomic_ref` / `genomic_alt`.
- Variant records preserve both layers of representation.

Where:
- `src/utils/parsers.py`
- `src/data_prep/01_curate_variants.py`

Why this matters:
- Avoids mixing strand-corrected CDS alleles with raw genomic alleles.
- Makes downstream checks against gnomAD/MyVariant and mtDNA resources much safer.

## 4. Transcript mismatches between ClinVar NM_ and TOGA ENST are no longer handled only by local guesswork

Problem:
ClinVar protein and coding coordinates are often attached to an NM_ transcript that does not match the ENST used in the TOGA alignment.

How this was fixed:
- A direct NM_ to ENST amino-acid position map can be built in `00f_build_transcript_position_maps.py`.
- The classifier uses this map when available instead of relying only on local anchor search.

Where:
- `src/data_prep/00f_build_transcript_position_maps.py`
- `src/utils/alignment_parser.py`
- `src/classify/00_classify_DAV.py`

Why this matters:
- Reduces widespread off-by-offset errors in nucDNA cDAV calls.
- Prevents systematic misclassification when TOGA and ClinVar disagree on isoform choice.

## 5. A second rescue strategy exists for the hardest TOGA transcript mismatches

Problem:
For some genes, the current MANE-matched NM_ transcript is not actually the isoform that matches the TOGA ENST, so a pure NM_ to ENST map can still be wrong.

How this was fixed:
- `00g_build_genomic_coordinate_maps.py` builds a genomic-position to TOGA-CDS/AA map directly from ENST exon structures.
- The classifier uses genomic position as an override for those genes before attempting ordinary transcript remapping.

Where:
- `src/data_prep/00g_build_genomic_coordinate_maps.py`
- `src/classify/00_classify_DAV.py`

Why this matters:
- This is the strongest fix for transcript-version drift because genomic coordinates are more stable than transcript-relative coordinates.

## 6. cDAV codon construction is now done at the corrected alignment position

Problem:
If the annotated ClinVar coordinate is shifted relative to the TOGA alignment, building the mutant codon at the raw coordinate produces the wrong mutant amino acid.

How this was fixed:
- `AlignmentParser.check_compensation()` first resolves the corrected amino-acid position, then shifts the nucleotide position accordingly, then constructs the mutant codon.

Where:
- `src/utils/alignment_parser.py`
- `src/classify/00_classify_DAV.py`

Why this matters:
- Prevents false transcript mismatches caused only by isoform offsets.
- Makes AA-level and NT-level cDAV calls consistent with the corrected alignment coordinate.

## 7. The classifier now records undetermined cases instead of forcing a call

Problem:
Variants with missing alignments, transcript mismatch, anchor failure, or absent ENST coverage should not be treated as uDAVs.

How this was fixed:
- The classifier assigns `None`-like unresolved states for cDAV fields and stores a structured mismatch reason.
- Mismatch categories include anchor failure, transcript mismatch, missing alignment, codon extraction failure, and genomic position not in ENST.

Where:
- `src/classify/00_classify_DAV.py`

Why this matters:
- Prevents unresolved variants from contaminating the uncompensated set.
- Preserves a clear audit trail for later debugging.

## 8. Non-biological sequence states are filtered out of species-level cDAV calls

Problem:
Masked or invalid states such as gaps, frameshift markers, stops, and unknown residues should not count as support for or against compensation.

How this was fixed:
- `AlignmentParser.check_compensation()` excludes species carrying `-`, `!`, `*`, or `X` at the relevant amino-acid or codon positions.

Where:
- `src/utils/alignment_parser.py`

Why this matters:
- Prevents damaged alignments or codon-disrupted species from inflating or suppressing cDAV counts.

## 9. Cross-genome species matching was redesigned around TaxID rather than name heuristics

Problem:
Species-name matching between TOGA and mtDNA sources is brittle because names may differ even when the taxonomy is the same.

How this was fixed:
- `00d_taxid_species_map.py` now extracts TaxIDs directly from sanitized FASTA headers and builds the overlap from TaxID.

Where:
- `src/data_prep/00d_taxid_species_map.py`
- `src/data_prep/02_sanitize_all_alignments.py`

Why this matters:
- Greatly reduces false mismatches between mtDNA and nucDNA species sets.
- Makes cross-genome overlap reproducible and taxonomically stable.

## 10. Structural mapping no longer depends only on chain labels or exact sequence identity

Problem:
PDB chain naming and sequence content are often inconsistent with gene naming and reference isoforms.

How this was fixed:
- Chain-to-gene assignment first uses the RCSB API, then falls back to sequence-based assignment.
- Residue mapping uses global alignment after chain assignment.
- Final residue identity is checked against the expected reference amino acid.

Where:
- `src/structural/00_map_davs_to_structure.py`

Why this matters:
- Prevents many silent chain-mapping failures.
- Adds a final amino-acid sanity check before a DAR is declared structurally mapped.

## 11. Isoform proxy handling exists for a few structurally absent subunits

Problem:
Some preferred structures contain one human isoform but not another closely related isoform used in the variant dataset.

How this was fixed:
- A small explicit proxy table maps certain tissue-specific isoforms to their structural equivalents.
- The mapping status is tagged so proxy use is visible downstream.

Where:
- `src/structural/00_map_davs_to_structure.py`

Why this matters:
- Keeps otherwise mappable variants in the analysis while preserving traceability.

## 12. Structural contact extraction is now explicit and reproducible

Problem:
Earlier structural-contact logic can easily become ad hoc if the distance rule and contact classes are not clearly encoded.

How this was fixed:
- The structural script explicitly defines:
  - the residue used for contact geometry (`Cβ`, or `Cα` for glycine),
  - the 8 Å neighborhood rule,
  - the contact class priority order (`hbond`, `electrostatic`, `hydrophobic`, `vdw`).

Where:
- `src/structural/00_map_davs_to_structure.py`

Why this matters:
- Makes downstream structural contact calls consistent and auditable.

## 13. Fisher significance is no longer presented as the only structural partner signal

Problem:
Species-level Fisher tests are vulnerable to phylogenetic clustering and should not be the only basis for “compensatory partner” calls.

How this was fixed:
- The structural partner script now treats Fisher as a retained comparison column and adds phylogenetically aware columns for Pagel and branch co-occurrence analyses.
- It also writes all tested pairs, rather than only the initially filtered subset.

Where:
- `src/structural/01_find_compensating_partners.py`

Why this matters:
- Prevents circular filtering and preserves the full tested landscape for re-analysis.

## 14. The pipeline now has explicit invariant-based checks around cDAV output

Problem:
Classification outputs are easy to corrupt silently if AA-level and NT-level support or mismatch logic drift over time.

How this was fixed:
- The test suite includes invariant checks for AA/NT consistency, discarded/synonymous exclusion, and direct alignment-backed verification of cDAV claims.

Where:
- `tests/test_cdav_invariants.py`
- `tests/test_DAV_classification.py`
- `tests/test_parsers.py`

Why this matters:
- Gives you regression coverage for the most failure-prone classification logic.

