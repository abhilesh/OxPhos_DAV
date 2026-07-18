## Scope

This note plans three closely related implementation steps before downstream structural and phylogenetic work:

1. add stage-level QC summaries
2. formalize `classified_all` / `classified_clean` / `classified_warning`
3. lock the exception registry as a stage contract
4. add cross-source duplicate linkage between mitochondrial `ClinVar` and `MITOMAP`

The goal is to keep the pipeline metadata-complete and filter-late while making downstream analytical subsets explicit and reproducible.

## Implementation Status

Implemented on `2026-04-20` in the current codebase.

Current outputs now written by the pipeline:

- `data/derived/curated/cross_source_duplicate_groups.tsv`
- `data/derived/results/clinvar_mitomap_cross_source_overlap_summary.json`
- `data/derived/results/clinvar_mt_variants_present_in_mitomap.tsv`
- `data/derived/results/clinvar_mt_variants_absent_from_mitomap.tsv`
- `data/derived/classified/classified_all.parquet`
- `data/derived/classified/classified_clean.parquet`
- `data/derived/classified/classified_warning.parquet`
- `data/derived/classified/classification_qc_summary.json`
- `data/derived/classified/classification_qc_checks.tsv`

Current stage-level QC status:

- all `6` required QC checks pass
- exception registry provenance is recorded by path and `sha256`
- classified subset counts:
  - `classified_all = 6809`
  - `classified_clean = 6340`
  - `classified_warning = 469`

Schema nuance locked in by this implementation:

- `classification_basis` must distinguish:
  - `nt_and_aa`
  - `nt_only`
  - `aa_only`
  - `no_disease_allele_detected`
- even in a missense-focused pipeline, `nt_only` remains a valid theoretical and
  schema-level outcome because cross-species nucleotide recurrence does not
  guarantee amino-acid recurrence in the non-human codon context
- in the current audited dataset there are no `nt_only` classified rows
- mitochondrial-encoded `ClinVar` rows remain retained for traceability but are
  not eligible for the active nucDNA comparative branch

## Current Audited State

Cross-source mitochondrial overlap was audited from:

- `data/derived/curated/variants_master_curated.parquet`

Current overlap definition:

- nucleotide-level match on:
  - `interpreted_gene`
  - `genomic_pos`
  - `genomic_ref_nt` if present, otherwise `ref_nt`
  - `genomic_alt_nt` if present, otherwise `alt_nt`

Current counts:

- mitochondrial `ClinVar` rows retained for traceability: `4582`
- unique mitochondrial `ClinVar` variant IDs: `2291`
- `MITOMAP` rows: `388`
- shared unique mitochondrial nucleotide events between `ClinVar` and `MITOMAP`: `234`
- mitochondrial `ClinVar` rows present in `MITOMAP`: `468`
- unique mitochondrial `ClinVar` variant IDs present in `MITOMAP`: `234`
- mitochondrial `ClinVar` rows absent from `MITOMAP`: `4114`
- unique mitochondrial `ClinVar` variant IDs absent from `MITOMAP`: `2057`

Important interpretation:

- the biologically meaningful overlap is `234` unique shared nucleotide events
- the `468` overlapping `ClinVar` rows reflect duplicated curation interpretations per source record, not `468` distinct biological events
- this is expected under Scope-B traceability and should be represented explicitly as cross-source linkage, not destructive deduplication

Explicit variant inventories written from the audit:

- [clinvar_mitomap_cross_source_overlap_summary.json](/Users/ad2347/Documents/OxPhos_DAV/data/derived/results/clinvar_mitomap_cross_source_overlap_summary.json)
- [clinvar_mt_variants_present_in_mitomap.tsv](/Users/ad2347/Documents/OxPhos_DAV/data/derived/results/clinvar_mt_variants_present_in_mitomap.tsv)
- [clinvar_mt_variants_absent_from_mitomap.tsv](/Users/ad2347/Documents/OxPhos_DAV/data/derived/results/clinvar_mt_variants_absent_from_mitomap.tsv)

These two TSVs are the explicit record-level answers to:

- which mitochondrial `ClinVar` variants are also represented in `MITOMAP`
- which mitochondrial `ClinVar` variants are not represented in `MITOMAP`

## Implementation Plan

### 1. Add stage-level QC summaries

Add a deterministic QC summary writer at the end of classification.

Canonical outputs:

- `data/derived/classified/classification_qc_summary.json`
- `data/derived/classified/classification_qc_checks.tsv`
- optional compatibility section appended to `notes/running_results.md`

Required checks:

- no mitochondrial-encoded `ClinVar` row is eligible for the nuclear comparative branch
- no overlap-derived mtDNA row is lost during classification
- no `unresolved` row is counted as `uDAV`
- every classified row has a non-null, valid `classification_basis`
- every classified row has a resolved `alignment_source`
- every classified row has a resolved coordinate method or a documented mtDNA/source exception

Recommended implementation details:

- define each QC check as a named function returning:
  - `check_name`
  - `status`
  - `n_failed`
  - `failure_rule`
  - `example_variant_ids`
- write failures even when `n_failed = 0` so the report is schema-stable
- treat QC as stage metadata, not just console logging

Expected code location:

- primary logic in `src/classify/00_classify_DAV.py`
- low-level summary helper in `src/utils/`

### 2. Formalize classified subsets

The classification stage should emit explicit downstream subsets rather than relying on ad hoc filtering later.

Canonical subset products:

- `classified_all.parquet`
- `classified_clean.parquet`
- `classified_warning.parquet`

Definitions:

- `classified_all`
  - all rows with `classification_status == "classified"`
- `classified_warning`
  - classified rows with any of:
    - `exception_applied == True`
    - `mismatch_code == "REF_ALLELE_MISMATCH"`
    - future explicit warning flags
- `classified_clean`
  - classified rows excluding `classified_warning`

Rules:

- no row should disappear from the canonical classified master table
- these are derived views, not replacements for the master table
- counts for all three subsets should be reported in QC metadata

Recommended extension fields:

- `classification_subset`
- `classification_warning_reason`
- `is_clean_classified`

### 3. Preserve the exception registry as a locked stage contract

The exception registry is now part of classification correctness and should be treated as an input artifact, not informal patch logic.

Required contract rules:

- classification behavior must read exceptions only from:
  - `data/derived/reference/variant_exception_registry.tsv`
- code should fail loudly if exception fields are referenced but the registry is missing
- no inline hard-coded rescue exceptions should remain in classification logic
- every exception-bearing classified row must retain:
  - `exception_scope`
  - `exception_class`
  - `exception_code`
  - `exception_decision`
  - `classification_exception_action`

Recommended additional metadata:

- `exception_registry_version`
- `exception_registry_sha256`
- `exception_registry_loaded_at_runtime`

Recommended enforcement:

- write a small registry validation step before classification begins
- write registry provenance into `classification_qc_summary.json`

### 4. Add cross-source duplicate linkage

This should be implemented in `data_curation`, because it is source provenance over curated biological events.

The objective is not to collapse source records. The objective is to link records that represent the same biological mitochondrial nucleotide event across `ClinVar` and `MITOMAP`.

Canonical linkage fields to add to curated rows:

- `cross_source_match_key_nt`
- `cross_source_duplicate_group_id`
- `cross_source_duplicate_status`
- `matched_in_clinvar`
- `matched_in_mitomap`
- `cross_source_partner_count`
- `cross_source_partner_variant_ids`
- `cross_source_partner_sources`

Recommended status values:

- `shared_mt_clinvar_mitomap`
- `clinvar_only_mt`
- `mitomap_only_mt`
- `not_applicable`

Recommended duplicate-group rule:

- only assign duplicate groups for mitochondrial rows
- use the nucleotide key:
  - `interpreted_gene`
  - `genomic_pos`
  - normalized reference allele
  - normalized alternate allele
- assign a stable group ID such as:
  - `XSDUP:MT:<gene>:<pos>:<ref>:<alt>`

Important restrictions:

- do not merge or drop rows
- do not infer amino-acid duplication from nucleotide duplication
- for `MT-ATP6/MT-ATP8` overlap rows, linkage remains frame-specific through `interpreted_gene`
- if later desired, an additional event-level overlap abstraction can be added above the frame-specific layer

Recommended outputs:

- curated master table updated with cross-source linkage fields
- `data/derived/curated/cross_source_duplicate_groups.tsv`
- `data/derived/results/clinvar_mitomap_cross_source_overlap_summary.json`

### 5. Explicitly preserve the “ClinVar not in MITOMAP” inventory

This is not just a planning note. It should remain a stable derived result because it answers a manuscript-relevant provenance question.

Keep as generated artifacts:

- `data/derived/results/clinvar_mt_variants_present_in_mitomap.tsv`
- `data/derived/results/clinvar_mt_variants_absent_from_mitomap.tsv`

The absent file is the explicit inventory of mitochondrial `ClinVar` variants not represented in `MITOMAP` under the current nucleotide-level matching rule.

## Implementation Order

1. add cross-source linkage fields during curation
2. rerun curation
3. update classification to emit QC summaries and subset views
4. rerun classification
5. append stable summary metrics to `notes/running_results.md`
6. document the stage contract and exception policy in `notes/curation_and_classification_issues_and_fixes.md`

This order is preferable because classification QC should consume the final curated schema rather than a temporary one.

## Test Plan

### Cross-source linkage

- every mitochondrial `MITOMAP` row gets one of:
  - `shared_mt_clinvar_mitomap`
  - `mitomap_only_mt`
- every mitochondrial `ClinVar` row gets one of:
  - `shared_mt_clinvar_mitomap`
  - `clinvar_only_mt`
- the number of unique shared duplicate groups equals `234` for the current audited dataset
- overlap linkage never changes row counts in the master curated table

### Stage-level QC

- QC report writes even when all checks pass
- QC failures include example `variant_id` values
- `classification_basis` is valid for every classified row
- no `is_udav_*` flag is set on unresolved rows

### Classified subset views

- `classified_all`, `classified_clean`, and `classified_warning` can be regenerated deterministically from the master classified table
- `classified_clean` is a proper subset of `classified_all`
- `classified_warning` captures all exception-bearing classified rows

### Exception contract

- if the registry changes, provenance fields in the QC summary change accordingly
- no exception-bearing row is classified without the registry metadata being propagated

## Manual Review Notes

The following items remain manual scientific review tasks and should not be conflated with this implementation pass:

- auditing `Unassigned` complex genes
- revisiting whether some warning rows should remain in `classified_clean`
- deciding whether downstream manuscript tables should report:
  - source-level counts
  - biological-event counts
  - or both

The correct answer is probably both, but they must stay clearly labeled.
