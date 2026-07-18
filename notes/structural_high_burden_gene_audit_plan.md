# Structural High-Burden Gene Audit Plan

This document formalizes the next audit step for the structural mapping stage.

## Why this audit is needed

For nucDNA rows, there are three different coordinate systems in play:

1. the transcript/protein consequence carried by the ClinVar-derived curated row
2. the human TOGA ENST sequence used to define the comparative alignment space
3. the protein chain and residue numbering present in the selected structural
   model

These are not guaranteed to be the same transcript product or residue
numbering system.

The current pipeline already partially accounts for this by:

- reading the curated transcript-position maps built during `data_curation`
- remapping nucDNA amino-acid coordinates from the ClinVar/curation space into
  the TOGA ENST alignment space before structural anchoring
- aligning the human TOGA reference sequence to each structural chain

What it has not yet made explicit enough is:

- which ClinVar transcript ID was present on each row
- which preferred NM transcript and TOGA ENST were used for reconciliation
- whether the structural attempt used an identity map, a remapped coordinate, or
  had no valid transcript projection
- how large the coordinate shift was between ClinVar space and TOGA space

That makes it harder to tell whether a structural failure is due to:

- transcript mismatch
- mature-protein offset
- unresolved structure segment
- absent chain/model
- isoform proxy limitations

## Audit goals

The high-burden audit should answer:

1. Which failures are likely transcript-reconciliation problems?
2. Which failures are likely mature-protein/processed-protein offset problems?
3. Which failures are mainly due to unresolved density or incomplete model
   coverage?
4. Which failures are due to chain absence rather than coordinate mismatch?
5. Which genes are poor structural candidates in the current human panel and
   should be treated as sensitivity classes rather than rescued automatically?

## High-priority gene set

The first pass should focus on genes with the largest residual failure burden.

### Residue-anchoring priority genes

- `NDUFV3`
- `SDHA`
- `SDHB`
- `SDHC`
- `SDHD`
- `NDUFS3`
- `NDUFS4`
- `NDUFS6`
- `NDUFS7`
- `NDUFA9`
- `NDUFA10`
- `NDUFAB1`
- `COX6A1`
- `ATP5F1A`
- `ATP5F1B`
- `ATP5F1D`
- `ATP5MC2`
- `ATP5MC3`

### Chain/model-gap priority genes

- `CYC1`
- `UQCRC1`
- `UQCRC2`
- `COX4I1`
- `COX5A`
- `COX6A1`
- `COX6B1`
- `UQCRFS1`
- `UQCRB`
- `UQCRQ`
- `UQCRH`

### Isoform-proxy priority genes

- `COX4I2`
- `COX6A2`
- `COX7A1`

## Required audit outputs

The structural stage should write explicit transcript-reconciliation metadata so
the high-burden gene review can be done from output tables rather than code
inspection alone.

Required per-row fields:

- `clinvar_transcript_id`
- `preferred_nm`
- `toga_enst`
- `transcript_map_type`
- `transcript_map_identity`
- `transcript_map_coverage`
- `transcript_reconciliation_status`
- `transcript_coord_delta`

Recommended statuses:

- `mt_native`
- `no_transcript_map`
- `transcript_identity`
- `transcript_remapped`
- `position_not_in_enst`

Required derived audit products:

- `results/structural/structure_transcript_reconciliation_audit.csv`
- `results/structural/structure_mapping_failure_audit.csv`

## Expected failure classes

The audit should separate at least these classes:

### 1. Transcript / sequence-model mismatch

Typical signature:

- nonzero transcript coordinate shift or `position_not_in_enst`
- persistent structural anchoring failure after transcript remapping

### 2. Mature-protein offset candidate

Typical signature:

- transcript reconciliation succeeds
- chain exists
- residue anchors only if a larger-than-default shift is allowed
- failures cluster near termini

### 3. Unresolved structure segment

Typical signature:

- transcript reconciliation succeeds
- chain exists
- local mapping region is plausible
- residue itself is not built in the structure

### 4. Chain/model absence

Typical signature:

- no chain for the gene in the selected model
- often concentrated in specific accessory subunits and specific PDB entries

### 5. Isoform proxy limitation

Typical signature:

- direct chain absent
- proxy gene available or partially available
- mapping confidence depends on isoform substitution rather than direct
  representation

## Implementation plan

### Step 1. Make transcript reconciliation explicit

Update the structural mapper so each nucDNA row records:

- original ClinVar transcript ID from the classified row
- preferred NM and TOGA ENST from the transcript map
- transcript map quality values
- whether the coordinate was unchanged, remapped, or not projectable
- the raw-to-remapped amino-acid coordinate difference

### Step 2. Emit transcript audit outputs

Write an explicit transcript reconciliation audit table that summarizes failures
and mappings by:

- gene
- reconciliation status
- transcript map type
- map quality
- mapping outcome

### Step 3. Use the new audit to guide rescue rules

Only after transcript reconciliation is explicit should we add further rescue
logic such as:

- wider bounded anchor windows for mature-protein offsets
- local sequence-window rescue
- model-specific preferred anchoring rules

## Interpretation rule

For nucDNA, the structural mapper is not using the ClinVar transcript directly
as the structure coordinate system.

It is using:

- ClinVar consequence position as the starting coordinate
- curation-derived NM→ENST mapping to project into TOGA human alignment space
- human TOGA reference sequence to structural chain alignment to project into
  structure space

This is scientifically defensible, but only if the reconciliation metadata is
recorded explicitly and audited gene by gene.
