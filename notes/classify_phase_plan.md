# Phase 3 Plan: `classify` Refactor Against the Canonical Curated Master Table

This note defines the next refactor target for `src/classify/00_classify_DAV.py` so the classification stage is consistent with:

- the current `data_download` and `data_curation` refactor
- the filter-late rule in `notes/pipeline_consistency_and_manuscript_prep.md`
- the resolved and unresolved issues recorded in `notes/resolved_nuances_and_fixes.md` and `notes/unresolved_issues_ranked.md`

It is written against the current curated schema in:

- `data/derived/curated/variants_master_curated.parquet`
- `data/derived/curated/variants_master_curated.jsonl`

and against the current implementation state in:

- `src/data_curation/01_curate_variants.py`
- `src/data_curation/02_build_transcript_position_maps.py`
- `src/data_curation/03_build_genomic_coordinate_maps.py`
- `src/classify/00_classify_DAV.py`

---

## 1. Validation Status Before Classification Refactor

The `data_curation` layer is in a usable state for planning the classify refactor.

What is validated:

- the canonical curated outputs exist under `data/derived/curated/`
- the flattened `src/data_curation/*.py` scripts compile cleanly
- the current curated JSONL confirms the expected filter-late fields are present, including:
  - `variant_id`
  - `source_variant_group_id`
  - `loci_raw`
  - `interpreted_gene`
  - `is_overlap`
  - `eligible_core_comparative_pipeline`
  - `core_pipeline_exclusion_reason`
  - `parse_status`
  - `frame_specific_*`
  - `coordinate_resolution_*`
  - `tier`
- overlap-aware mtDNA rows are present in the curated outputs

Important caveat:

- this was validated from existing curated outputs plus script compile checks
- I did not rerun the full `data_curation` stage in the current host environment because the host Python is missing some analysis dependencies used by parts of the pipeline

That is sufficient for planning and code refactor design, but the classify refactor should still include a clean end-to-end rerun check in the container once implemented.

---

## 2. Current Classify-Stage Mismatches With the Refactored Pipeline

The current `src/classify/00_classify_DAV.py` is still aligned to the old pipeline contract, not the current curation contract.

### 2.1 It still reads compatibility JSONs instead of the canonical curated master table

Current behavior:

- reads:
  - `data/annotations/curated/mtDNA_annotations.json`
  - `data/annotations/curated/nucDNA_annotations.json`

Why this is now wrong:

- those files are compatibility views
- the canonical source of truth is now `data/derived/curated/variants_master_curated.parquet`

Required fix:

- classification must read the canonical curated master table
- JSON compatibility views may be kept only as optional export targets, not as the primary input contract

### 2.2 It still drops records by policy instead of retaining them with metadata

Current behavior:

- skips records when:
  - `"Discarded" in tier`
  - `var["is_synonymous"]`

Why this is now wrong:

- the refactored pipeline is explicitly filter-late
- exclusion by absence is no longer acceptable for core stage outputs

Required fix:

- keep every curated row in the classification output
- annotate ineligibility with metadata fields such as:
  - `classification_eligible`
  - `classification_exclusion_reason`
  - `classification_status`

### 2.3 It still collapses overlap loci using `locus.split("/")[0]`

Current behavior:

- uses `var.get("locus", "").split("/")[0]`

Why this is now wrong:

- the curation layer already emits frame-specific overlap rows
- the authoritative gene for classification is now `interpreted_gene`, not `locus.split("/")[0]`

Required fix:

- always use:
  - `interpreted_gene` for gene identity
  - `loci_raw` only for provenance
  - `is_overlap` and `overlap_group` for overlap-aware annotation

### 2.4 It still assumes the old coordinate parsing model

Current behavior:

- infers amino-acid and coding coordinates from:
  - `aa_change`
  - `hgvs_c`
  - `locus`

Why this is now incomplete:

- the curated schema now provides richer and more explicit fields:
  - `frame_specific_hgvs_p`
  - `frame_specific_ref_aa`
  - `frame_specific_alt_aa`
  - `frame_specific_is_synonymous`
  - `frame_specific_is_missense`
  - `coordinate_resolution_status`
  - `coding_ref_nt`
  - `coding_alt_nt`
  - `genomic_ref_nt`
  - `genomic_alt_nt`

Required fix:

- the classifier must prefer explicit curated fields over reparsing free-text fields
- reparsing should only be a fallback with explicit failure logging

### 2.5 It still writes old-style classified JSON outputs only

Current behavior:

- writes:
  - `cdav_classifications_mtDNA.json`
  - `cdav_classifications_nucDNA.json`

Why this is now incomplete:

- the current pipeline should produce a canonical classified master table
- genome-specific JSON files can remain as compatibility views only

Required fix:

- write canonical classified outputs under `data/derived/classified/`

---

## 3. Resolved Upstream Issues the Classifier Must Respect

These issues are already addressed upstream and should not be re-broken during classify refactor.

### 3.1 Overlap-aware mtDNA curation already exists

Implication for classify:

- do not collapse `MT-ATP6/MT-ATP8`
- classify each curated row independently
- preserve the shared source event relationship via:
  - `source_variant_group_id`
  - `shared_source_variant_group_id`

### 3.2 ND6 strand correction already exists in curation

Implication for classify:

- use `coding_ref_nt` and `coding_alt_nt` where classification is in CDS space
- do not try to re-derive ND6 strand logic from genomic alleles

### 3.3 Transcript and genomic rescue maps are already first-class curated support products

Implication for classify:

- continue using transcript and genomic map rescue logic
- but load them from the current derived curated/reference outputs, not the old `data/reference/*.json` assumptions as the primary contract

### 3.4 Unresolved mapping categories already have a strong conceptual model

Implication for classify:

- continue preserving unresolved states
- never coerce unresolved rows into the uDAV set

---

## 4. Unresolved Issues the Classifier Must Explicitly Solve

These are the main unresolved items from `notes/unresolved_issues_ranked.md` that are classify-relevant.

### 4.1 Overlap-locus collapse

Status:

- unresolved in classify

Fix requirement:

- replace all `split("/")[0]` logic with `interpreted_gene`
- use overlap-aware fields already emitted by curation

### 4.2 Filter-late violation by skipping rows

Status:

- unresolved in classify

Fix requirement:

- keep all rows
- annotate policy ineligibility explicitly

### 4.3 Fallback amino-acid anchoring remains ambiguous

Status:

- still unresolved conceptually

Fix requirement:

- preserve the current anchor method for compatibility
- add metadata that distinguishes:
  - direct transcript-map resolution
  - genomic-map resolution
  - anchor-based fallback
  - ambiguous anchor
  - anchor failure

Strong recommendation:

- do not silently accept a non-unique fallback anchor
- ambiguous anchors should be unresolved, not forced

### 4.4 Residual transcript/alignment mismatch remains a legitimate unresolved category

Status:

- acceptable if explicitly recorded

Fix requirement:

- maintain and expand explicit unresolved categories
- ensure downstream analyses can exclude unresolved rows without conflating them with uDAVs

---

## 5. Canonical Input Contract for the Refactored Classifier

The refactored classifier should use the following as canonical inputs.

### 5.1 Required inputs

- `data/derived/curated/variants_master_curated.parquet`
- `data/derived/curated/transcript_position_maps.parquet` or normalized JSON export
- `data/derived/curated/genomic_coordinate_maps.parquet` or normalized JSON export
- `data/derived/curated/alignment_sanitation_manifest.parquet`
- alignment FASTAs under:
  - `data/alignments/toga_hg38_aa/`
  - `data/alignments/toga_hg38_codon/`
  - `data/alignments/mtdna_aa/`
  - `data/alignments/mtdna_codon/`
- `data/reference/mtdna_gene_coordinates.tsv` only if needed as a compatibility source

### 5.2 Input row policy

Every row from the curated master table should enter classification.

Classification should decide, per row:

- eligible for comparative classification
- unresolved after attempting classification
- classified as cDAV / uDAV
- skipped by policy but retained

---

## 6. Canonical Output Contract for the Refactored Classifier

The classifier should write canonical outputs under:

- `data/derived/classified/`

Required outputs:

- `variants_master_classified.parquet`
- `variants_master_classified.jsonl`
- `classification_run_metadata.json`
- `classification_mismatch_log.jsonl`
- optional compatibility exports:
  - `data/annotations/curated/cdav_classifications_mtDNA.json`
  - `data/annotations/curated/cdav_classifications_nucDNA.json`

The canonical analytical product should be:

- one row in, one row out relative to the curated master table

No row should disappear between curation and classification.

---

## 7. Required Classification Fields

These fields should be added without overwriting curated provenance.

### 7.1 Eligibility and status

- `classification_status`
  - `classified`
  - `unresolved`
  - `skipped_by_policy`
- `classification_eligible`
- `classification_exclusion_reason`
  - `core_pipeline_ineligible`
  - `synonymous`
  - `nonsense`
  - `frameshift`
  - `non_snv`
  - `noncoding`
  - `unparsable_aa`
  - `missing_alignment`
  - `position_not_in_enst`
  - `anchor_ambiguous`
  - `anchor_not_found`
  - `codon_extraction_failure`
  - `transcript_mismatch`
  - `genomic_pos_not_in_enst`
  - `ref_allele_mismatch`

### 7.2 Coordinate-resolution provenance

- `classification_gene`
  - should equal `interpreted_gene`
- `classification_gene_source`
  - `interpreted_gene`
- `classification_coordinate_method`
  - `frame_specific_curated`
  - `genomic_map`
  - `transcript_map`
  - `anchor_fallback`
- `classification_coordinate_status`
- `classification_used_overlap_frame`

### 7.3 cDAV/uDAV calls

- `is_udav_amino_acid`
- `is_udav_nucleotide`
- `is_cdav_amino_acid`
- `is_cdav_nucleotide`
- `classification_basis`
  - `aa_only`
  - `nt_and_aa`
  - `unresolved`

### 7.4 Species support and support strength

- `n_species_aligned`
- `n_species_with_disease_allele`
- `lineages_with_disease_allele`
- `n_species_with_disease_codon`
- `lineages_with_disease_codon`
- `cdav_support_level`
  - `none`
  - `single_species`
  - `multi_species`
  - `multi_origin` if later derivable phylogenetically

### 7.5 Logging and audit fields

- `ref_allele_match`
- `mismatch_reason`
- `mismatch_code`
- `alignment_file_aa`
- `alignment_file_nt`
- `alignment_source`
  - `TOGA`
  - `mtDNA_MACSE`

---

## 8. Exact Logic Changes Required in `00_classify_DAV.py`

### 8.1 Input loader

Replace:

- `MT_CURATED`
- `NUC_CURATED`

With:

- load canonical classified input from `variants_master_curated.parquet`
- optionally export by genome later

### 8.2 Record identity

Replace:

- `ann_id` as the main classifier identifier

With:

- `variant_id` as the canonical row identifier
- preserve `source_variant_group_id` for grouping overlap-derived rows

### 8.3 Gene identity

Replace:

- `locus.split("/")[0]`

With:

- `interpreted_gene`

### 8.4 Eligibility gate

Replace:

- `if "Discarded" in tier or var["is_synonymous"]: continue`

With:

- classify every row into one of:
  - `classified`
  - `unresolved`
  - `skipped_by_policy`

Base rule:

- if `eligible_core_comparative_pipeline` is `False`, keep the row and mark:
  - `classification_status = skipped_by_policy`
  - `classification_eligible = False`
  - `classification_exclusion_reason = core_pipeline_exclusion_reason`

### 8.5 Coordinate parsing

Prefer, in order:

1. `frame_specific_ref_aa`
2. `frame_specific_alt_aa`
3. `frame_specific_hgvs_p`
4. `hgvs_p`
5. `aa_change`

Prefer nucleotide fields:

1. `coding_ref_nt`
2. `coding_alt_nt`
3. `ref_nt`
4. `alt_nt`

### 8.6 mtDNA coordinate conversion

Use:

- `interpreted_gene`
- `genomic_pos`
- current mtDNA coordinate table

Do not use:

- `locus.split("/")[0]`

### 8.7 Output layer

Append classification fields onto the curated row rather than creating an unrelated reduced object.

---

## 9. Recommended Classify Refactor Order

### Step 1. Convert input/output contract

- read canonical curated master table
- write canonical classified master table
- keep compatibility JSON exports only as derived outputs

### Step 2. Replace gene-identity and overlap handling

- remove `split("/")[0]`
- use `interpreted_gene`
- add overlap-aware metadata to classification outputs

### Step 3. Replace skip-by-absence logic

- retain ineligible and unresolved rows
- introduce `classification_status`, `classification_eligible`, and `classification_exclusion_reason`

### Step 4. Upgrade coordinate parsing

- use explicit curated frame-specific fields first
- keep free-text reparsing only as fallback

### Step 5. Improve fallback anchoring transparency

- distinguish transcript-map, genomic-map, and anchor-fallback methods
- record ambiguous-anchor states explicitly

### Step 6. Update downstream compatibility expectations

- structural scripts should consume `variants_master_classified.parquet`
- compatibility JSONs should remain temporary only if needed

---

## 10. Tests Required Before Accepting the Classify Refactor

### 10.1 Master-table preservation

- number of rows in `variants_master_classified.parquet` equals number of rows in `variants_master_curated.parquet`

### 10.2 Overlap handling

- overlap-derived mtDNA rows remain duplicated after classification
- `MT-ATP6` and `MT-ATP8` interpretations are classified independently
- no classify code path uses `split("/")[0]`

### 10.3 Policy retention

- rows with `eligible_core_comparative_pipeline = False` are retained
- such rows receive explicit `classification_status = skipped_by_policy`

### 10.4 Resolved vs unresolved separation

- unresolved rows do not become uDAVs
- cDAV and uDAV are never assigned when coordinate resolution fails

### 10.5 Coordinate-method provenance

- each classified row records whether the call came from:
  - genomic-map rescue
  - transcript-map rescue
  - anchor fallback

### 10.6 Schema consistency

- `classification_gene` always equals `interpreted_gene`
- `variant_id` survives unchanged through classification
- canonical outputs can regenerate legacy JSON exports without semantic loss

---

## 11. Practical Decision

The classify refactor should be treated as:

- a schema migration
- a filter-late correction
- an overlap-handling correction

not just as a local patch to `00_classify_DAV.py`.

The key invariant for this phase is:

> Every curated row survives classification as a row, with either a valid call or an explicit status explaining why it could not be classified.

