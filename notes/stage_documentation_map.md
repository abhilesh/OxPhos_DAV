# Stage Documentation Map

This note defines the current documentation map for the OXPHOS DAV project.

Its purpose is to make clear:

- which notes are the current authoritative records
- which notes are supplementary deep dives
- which notes are planning/history rather than final stage records
- which note is the running results log

## Status labels

Use these labels when interpreting the notes:

- `Authoritative`: current best source of truth for a stage or topic
- `Supplementary`: important supporting detail, but not the main stage record
- `Historical / Planning`: useful background, but not the final record of what
  was implemented
- `Results Log`: cumulative output summary, not the main methods/decision record

## 1. Project-level notes

### [pipeline_consistency_and_manuscript_prep.md](/Users/ad2347/Documents/OxPhos_DAV/notes/pipeline_consistency_and_manuscript_prep.md)

Status:

- `Authoritative`

Role:

- overall project organization plan
- Scope-B traceability rule
- repo/manuscript preparation guidance
- stage organization logic

Use this for:

- top-level project contract and structure

## 2. Upstream pipeline: `data_download` / reference layer

### [resolved_nuances_and_fixes.md](/Users/ad2347/Documents/OxPhos_DAV/notes/resolved_nuances_and_fixes.md)

Status:

- `Authoritative`

Role:

- resolved upstream pipeline issues
- reference/download-layer fixes
- general resolved implementation nuances across the early pipeline

Use this for:

- current resolved-state record for `data_download` and shared upstream
  infrastructure

### [unresolved_issues_ranked.md](/Users/ad2347/Documents/OxPhos_DAV/notes/unresolved_issues_ranked.md)

Status:

- `Authoritative`

Role:

- ranked unresolved issues across the pipeline
- prioritization of what still needs work

Use this for:

- current unresolved-state record
- action prioritization

## 3. `data_curation` and `classify`

### [curation_and_classification_issues_and_fixes.md](/Users/ad2347/Documents/OxPhos_DAV/notes/curation_and_classification_issues_and_fixes.md)

Status:

- `Authoritative`

Role:

- main final record for the `data_curation` stage
- main final record for the `classify` stage
- documents problems, fixes, transcript-map rebuild, unresolved subsets,
  warning-bearing rows, and exception logic

Use this for:

- the implemented curation/classification contract
- stage-specific issue history and fixes

### [classify_phase_plan.md](/Users/ad2347/Documents/OxPhos_DAV/notes/classify_phase_plan.md)

Status:

- `Historical / Planning`

Role:

- design/refactor plan for classification

Use this for:

- understanding the intended classify refactor logic

Do not use this as the only record of current classify behavior.

### [stage_qc_and_cross_source_duplicate_linkage_plan.md](/Users/ad2347/Documents/OxPhos_DAV/notes/stage_qc_and_cross_source_duplicate_linkage_plan.md)

Status:

- `Supplementary`

Role:

- focused note for stage-level QC
- cross-source duplicate linkage logic

Use this for:

- QC and cross-database linkage details

## 4. Structural stage

### [structural_phase_plan.md](/Users/ad2347/Documents/OxPhos_DAV/notes/structural_phase_plan.md)

Status:

- `Authoritative`

Role:

- main structural-stage design and current-state note
- structural source validation
- structure panel logic
- current structural bottlenecks and policy decisions

Use this for:

- the current structural-stage contract and strategy

### [structural_mapping_major_problems_report.md](/Users/ad2347/Documents/OxPhos_DAV/notes/structural_mapping_major_problems_report.md)

Status:

- `Authoritative`

Role:

- concise problem report for the structural mapping stage
- summarizes major failure classes and affected genes

Use this for:

- quick structural-stage status review

### [structural_anchoring_refactor_notes.md](/Users/ad2347/Documents/OxPhos_DAV/notes/structural_anchoring_refactor_notes.md)

Status:

- `Supplementary`

Role:

- details of the anchoring refactor
- explains direct, diagnostic, and registry-gated offset logic

Use this for:

- anchoring-method provenance

### [structural_high_burden_gene_audit_plan.md](/Users/ad2347/Documents/OxPhos_DAV/notes/structural_high_burden_gene_audit_plan.md)

Status:

- `Supplementary`

Role:

- formal high-burden gene audit framework

Use this for:

- targeted structural gene-audit methodology

### [non_sdh_large_offset_audit.md](/Users/ad2347/Documents/OxPhos_DAV/notes/non_sdh_large_offset_audit.md)

Status:

- `Supplementary`

Role:

- focused audit of non-`SDH` large-offset candidates

Use this for:

- rationale for keeping the anchor registry conservative

### [gene_specific_pipeline_issues_by_stage.md](/Users/ad2347/Documents/OxPhos_DAV/notes/gene_specific_pipeline_issues_by_stage.md)

Status:

- `Supplementary`

Role:

- cross-stage gene-by-gene issue map

Use this for:

- following difficult genes across download/curation/classify/structural stages

## 5. Results

### [running_results.md](/Users/ad2347/Documents/OxPhos_DAV/notes/running_results.md)

Status:

- `Results Log`

Role:

- cumulative running results summary
- source inventory snapshot
- cDAV metrics
- structural metrics

Use this for:

- current output summaries and reported counts

Do not use this as the sole authoritative methods note for any stage.

## 6. Practical reading order

If you want the shortest path to the current pipeline state, read in this order:

1. [pipeline_consistency_and_manuscript_prep.md](/Users/ad2347/Documents/OxPhos_DAV/notes/pipeline_consistency_and_manuscript_prep.md)
2. [resolved_nuances_and_fixes.md](/Users/ad2347/Documents/OxPhos_DAV/notes/resolved_nuances_and_fixes.md)
3. [unresolved_issues_ranked.md](/Users/ad2347/Documents/OxPhos_DAV/notes/unresolved_issues_ranked.md)
4. [curation_and_classification_issues_and_fixes.md](/Users/ad2347/Documents/OxPhos_DAV/notes/curation_and_classification_issues_and_fixes.md)
5. [structural_phase_plan.md](/Users/ad2347/Documents/OxPhos_DAV/notes/structural_phase_plan.md)
6. [structural_mapping_major_problems_report.md](/Users/ad2347/Documents/OxPhos_DAV/notes/structural_mapping_major_problems_report.md)
7. [running_results.md](/Users/ad2347/Documents/OxPhos_DAV/notes/running_results.md)

## 7. Current recommendation

For manuscript-grade documentation, the current best stage mapping is:

- `data_download` / upstream reference layer:
  - [resolved_nuances_and_fixes.md](/Users/ad2347/Documents/OxPhos_DAV/notes/resolved_nuances_and_fixes.md)
  - [unresolved_issues_ranked.md](/Users/ad2347/Documents/OxPhos_DAV/notes/unresolved_issues_ranked.md)
- `data_curation`:
  - [curation_and_classification_issues_and_fixes.md](/Users/ad2347/Documents/OxPhos_DAV/notes/curation_and_classification_issues_and_fixes.md)
- `classify`:
  - [curation_and_classification_issues_and_fixes.md](/Users/ad2347/Documents/OxPhos_DAV/notes/curation_and_classification_issues_and_fixes.md)
- `structural`:
  - [structural_phase_plan.md](/Users/ad2347/Documents/OxPhos_DAV/notes/structural_phase_plan.md)
  - [structural_mapping_major_problems_report.md](/Users/ad2347/Documents/OxPhos_DAV/notes/structural_mapping_major_problems_report.md)
- results:
  - [running_results.md](/Users/ad2347/Documents/OxPhos_DAV/notes/running_results.md)
