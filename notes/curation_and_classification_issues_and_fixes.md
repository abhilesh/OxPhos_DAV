# Curation And Classification Issues, Problems, And Fixes

This note records the main issues encountered during the `data_curation` and
`classify` refactor, how they affected the pipeline, and how they were solved.

It is intended to document the transition from the older JSON- and
drop-filter-based workflow to the current metadata-first, canonical Parquet
workflow.

---

## 1. Curation-stage issues and fixes

### 1.1 Early-dropping records instead of retaining them with metadata

**Problem**

The older curation logic dropped many records during parsing or normalization,
including rows that were:

- outside the core missense-comparative branch
- noncoding or otherwise not eligible for downstream codon-aware analysis
- difficult to interpret from source fields alone

This violated the intended Scope-B, filter-late contract.

**Why this mattered**

- The absence of a row became overloaded to mean "excluded", "unparsable", or
  "never seen".
- Downstream denominators could become unstable.
- Traceability from raw source records to curated interpretations was weakened.

**Fix**

The curation layer was refactored so that:

- the canonical output is `data/derived/curated/variants_master_curated.parquet`
- rows are retained with metadata rather than dropped
- eligibility is expressed explicitly using fields such as:
  - `eligible_core_comparative_pipeline`
  - `core_pipeline_exclusion_reason`
  - `parse_status`
  - `curation_status`

**Result**

- Curation is now metadata-first and filter-late.
- Exclusion is represented explicitly rather than by row absence.

---

### 1.2 `MT-ATP6/MT-ATP8` overlap loci were biologically collapsed

**Problem**

Composite mitochondrial loci such as `MT-ATP6/MT-ATP8` were at risk of being
collapsed to one gene interpretation, typically by taking the first token.

**Why this mattered**

For overlap nucleotides, a single mtDNA variant can have two distinct
frame-specific consequences:

- missense in one gene
- synonymous or unresolved in the other

Collapsing the locus discards valid biology and breaks compensation analyses.

**Fix**

Overlap-aware curation was implemented in
[01_curate_variants.py](/Users/ad2347/Documents/OxPhos_DAV/src/data_curation/01_curate_variants.py):

- overlap loci are detected from the raw locus field
- overlap rows are duplicated into frame-specific interpretations
- each interpretation gets its own curated row
- overlap metadata is recorded explicitly:
  - `interpreted_gene`
  - `is_overlap`
  - `overlap_group`
  - `overlap_genes`
  - `overlap_role`
  - `derived_from_overlap_duplication`
  - `shared_source_variant_group_id`
- frame-specific amino-acid fields are emitted separately

**Result**

- `MT-ATP6` and `MT-ATP8` are now treated as separate biological
  interpretations of the same nucleotide event.
- Downstream classification can operate on the correct interpreted frame.

---

### 1.3 Gene identity was too loosely defined

**Problem**

Older code paths sometimes treated `locus` as the authoritative gene identity.
For overlap rows and some downstream logic, this could still imply ambiguous
gene assignment.

**Why this mattered**

- Complex assignment
- encoded-genome assignment
- overlap handling
- downstream classification logic

all need a single authoritative interpreted gene.

**Fix**

The curated schema was updated so that:

- `interpreted_gene` is the authoritative biological gene assignment
- `loci_raw` remains source provenance
- gene-context population in
  [variant_record.py](/Users/ad2347/Documents/OxPhos_DAV/src/utils/variant_record.py)
  now uses `interpreted_gene` and no longer falls back to collapsing composite
  loci via `split("/")[0]`

**Result**

- Curation and classification now agree on the same gene identity contract.

---

### 1.4 Transcript-position maps were conceptually correct but technically broken

**Problem**

The transcript-map builder in
[02_build_transcript_position_maps.py](/Users/ad2347/Documents/OxPhos_DAV/src/data_curation/02_build_transcript_position_maps.py)
was producing implausible NM-to-ENST mappings for many nuclear genes.

Observed symptoms included:

- extremely low `identity_fraction` values for genes known to be near-identical
  between the chosen NM transcript and the TOGA Homo sapiens ENST sequence
- large stretches of positions mapped to `None`
- many downstream `POSITION_NOT_IN_ENST` classification failures

**Root cause**

The TOGA amino-acid FASTA parser assumed one-line sequences by reading a header
line and then only the next sequence line. For multi-line FASTA entries this
truncated the human TOGA sequence badly.

This made otherwise valid NM-to-ENST comparisons appear mismatched.

**Fix**

The transcript-map builder was changed to:

- parse TOGA FASTA files correctly using `Bio.SeqIO.parse`
- read the full human sequence from the alignment file
- keep the TOGA ENST identifier from the human header
- compute map-quality metadata:
  - `identity_fraction`
  - `coverage_fraction`
  - `mapped_positions`
  - `nm_protein_length`
  - `enst_protein_length`
- record both:
  - `mane_nm`
  - `clinvar_dominant_nm`

**Result**

After rebuilding the transcript maps in Docker:

- all 87 nuclear genes produced mapped transcript entries
- the previously problematic genes such as `SDHA`, `SDHC`, `SDHD`, `NDUFS1`,
  `NDUFV1`, `COX4I2`, and `NDUFA6` now have biologically coherent maps
- most of the spurious nuclear unresolved burden disappeared downstream

---

### 1.5 Transcript selection needed explicit provenance

**Problem**

Using only one implicit transcript preference rule made it difficult to audit
why a given nuclear gene used a particular NM transcript.

**Fix**

The transcript-map output now records:

- `nm`
- `mane_nm`
- `clinvar_dominant_nm`
- `selection_rule`

**Result**

- Transcript choice is now auditable.
- Cases where MANE and dominant ClinVar usage differ can be identified and
  reviewed explicitly.

---

### 1.6 Canonical curation outputs needed to replace older compatibility-only views

**Problem**

The older workflow relied heavily on JSON compatibility outputs rather than a
single canonical curated analytical product.

**Fix**

The curation stage now writes canonical outputs under `data/derived/curated/`,
including:

- `variants_master_curated.parquet`
- `variants_master_curated.jsonl`
- `transcript_position_maps.parquet`
- `genomic_coordinate_maps.parquet`
- `alignment_sanitation_manifest.parquet`
- `curation_run_metadata.json`

Compatibility JSONs remain as derived exports only.

**Result**

- The curation stage now has a stable canonical contract for downstream stages.

---

## 2. Classification-stage issues and fixes

### 2.1 Classification still used older inputs instead of the curated master table

**Problem**

The older classifier was aligned to previous JSON-based outputs and not to the
new curated master-table contract.

**Fix**

The classifier was rewritten in
[00_classify_DAV.py](/Users/ad2347/Documents/OxPhos_DAV/src/classify/00_classify_DAV.py)
to:

- read from `data/derived/curated/variants_master_curated.parquet`
- write canonical outputs under `data/derived/classified/`
- retain compatibility JSONs only as derived views

**Result**

- Classification is now driven by the canonical curated schema.

---

### 2.2 Classification still violated the filter-late rule

**Problem**

The older classification logic skipped some records by policy rather than
retaining them as explicit states.

**Fix**

The new classifier retains all curated rows and marks them using:

- `classification_status`
- `classification_eligible`
- `classification_exclusion_reason`

Rows can now be:

- `classified`
- `unresolved`
- `skipped_by_policy`

**Result**

- Downstream analyses can keep full denominators and exclude rows explicitly.

---

### 2.3 Overlap loci were still at risk of being collapsed during classification

**Problem**

The older classify logic still relied on patterns such as
`locus.split("/")[0]`, which would collapse overlap-aware curation rows back to
one gene.

**Fix**

Classification was updated to:

- use `interpreted_gene` as the authoritative classification gene
- preserve overlap metadata
- record whether an overlap-derived row used a frame-specific interpretation via
  `classification_used_overlap_frame`

**Result**

- Overlap-derived `MT-ATP6` and `MT-ATP8` rows are classified independently.

---

### 2.4 Parquet-loaded values caused downstream JSON export failures

**Problem**

The first full classify run completed biologically but failed when exporting
JSON/JSONL because Parquet-loaded rows contained:

- `numpy.ndarray` values in list-like columns
- `NaN` values in optional frame-specific fields

This caused errors such as:

- `Object of type ndarray is not JSON serializable`

**Fix**

The classifier now normalizes loaded values before inference and export:

- NumPy arrays are converted to Python lists
- NumPy scalars are converted to native Python scalars
- `NaN` is converted to `None`
- optional amino-acid fields are normalized before use

**Result**

- Canonical Parquet, JSONL, and compatibility JSON outputs all write cleanly.

---

### 2.5 Nuclear mismatch categories were initially dominated by false failures

**Problem**

Before fixing transcript-map generation, nuclear classification produced a very
large unresolved burden:

- `POSITION_NOT_IN_ENST`: `3468`
- `GENOMIC_POS_NOT_IN_ENST`: `69`
- `ANCHOR_NOT_FOUND`: `13`
- `TRANSCRIPT_MISMATCH`: `9`
- `REF_ALLELE_MISMATCH`: `8`

The dominant genes were:

- `SDHA`
- `SDHD`
- `SDHC`
- `NDUFS1`
- `NDUFV1`

**Interpretation**

This looked at first like a widespread biological transcript discordance
problem, but in practice much of it was a technical artifact caused by broken
transcript maps.

**Fix**

After correcting transcript-map generation and rebuilding the maps, classification
was rerun in Docker.

**Result**

Nuclear classification improved substantially:

- classified rows increased from `2995` to `6486`
- unresolved rows dropped from `3559` to `68`
- AA-level nuclear cDAVs increased by `899`
- NT-level nuclear cDAVs increased by `760`

This was the single largest correctness gain in the classify refactor.

---

### 2.6 Genomic rescue was being applied too aggressively

**Problem**

In an intermediate classify version, any gene present in the genomic-rescue map
set was routed through genomic rescue preferentially, even when the transcript
map was already high quality.

That could manufacture `GENOMIC_POS_NOT_IN_ENST` failures unnecessarily.

**Fix**

The classifier now records transcript-map quality and only prefers genomic
rescue when transcript-map quality is poor:

- low identity
- low coverage
- or missing transcript map

Otherwise, it uses the transcript map directly.

**Result**

- Classification method selection is now data-driven rather than hard-coded by
  gene membership in the rescue set.
- Remaining unresolved categories are more likely to represent real residual
  discordance.

---

### 2.7 Mismatch categories were refined into meaningful unresolved states

**Problem**

Without explicit mismatch categories, unresolved rows are hard to interpret and
easy to confuse with uDAVs.

**Fix**

The classifier now logs structured mismatch categories, including:

- `REF_ALLELE_MISMATCH`
- `TRANSCRIPT_MISMATCH`
- `POSITION_NOT_IN_ENST`
- `GENOMIC_POS_NOT_IN_ENST`
- `ANCHOR_NOT_FOUND`
- `CODON_EXTRACTION_FAILURE`
- `NO_ALIGNMENT`
- `COORD_PARSE_FAILURE`

These are written to:

- `data/derived/classified/classification_mismatch_log.jsonl`

**Result**

- Unresolved rows are now auditable and can be excluded without conflating them
  with uDAVs.

---

## 3. Current classification state after the fixes

After the transcript-map rebuild and classify refactor, the current outputs are
written successfully in Docker and reflect the corrected curation/classify
contract.

### 3.1 Current classified counts

**mtDNA**

- total rows: `388`
- classified: `323`
- unresolved: `0`
- skipped by policy: `65`

**nucDNA**

- total rows: `33391`
- classified: `6486`
- unresolved: `68`
- skipped by policy: `26837`

### 3.2 Current cDAV counts

**mtDNA**

- AA-level cDAVs: `167`
- NT-level cDAVs: `149`

**nucDNA**

- AA-level cDAVs: `2046`
- NT-level cDAVs: `1745`

### 3.3 Current remaining nuclear mismatch categories

At the current state, remaining nuclear mismatches are:

- `GENOMIC_POS_NOT_IN_ENST`: `45`
- `POSITION_NOT_IN_ENST`: `2`
- `REF_ALLELE_MISMATCH`: `19`
- `TRANSCRIPT_MISMATCH`: `21`

Main residual genes:

- `NDUFS6`
- `NDUFA13`
- `NDUFA11`
- `UQCRB`
- `NDUFS7`
- `COXFA4L2`
- `ATP5PF`
- `NDUFV2`
- `NDUFB1`
- `NDUFA10`
- `ATP5MC2`

These now look like genuine residual transcript/genomic discordance or
reference/consequence conflicts rather than a pipeline-wide map-construction
failure.

---

## 4. Remaining residual issues

The major solved problems were technical and structural. The remaining issues
are now narrower and more biologically specific.

### 4.1 Low-identity transcript-map genes still need targeted review

Some genes still have weaker NM-to-ENST concordance, for example:

- `NDUFS6`
- `NDUFA13`
- `NDUFA11`
- `UQCRB`
- `NDUFS7`
- `COXFA4L2`
- `ATP5PF`
- `COX5A`

These are the best candidates for manual review or more tailored rescue logic.

### 4.2 Some rows still show true reference/consequence disagreement

The remaining `REF_ALLELE_MISMATCH` and `TRANSCRIPT_MISMATCH` rows, especially
in:

- `NDUFV2`
- `NDUFB1`
- `NDUFA10`
- `ATP5MC2`

should be reviewed row-by-row because they may reflect:

- transcript-version differences
- ClinVar consequence inconsistencies
- rescue-map coordinate drift
- genuine incompatibility between selected NM and TOGA ENST models

---

## 5. Summary

The most important solved issues were:

- moving curation to a canonical filter-late master table
- implementing overlap-aware mtDNA curation
- replacing composite-locus collapse with `interpreted_gene`
- fixing transcript-map generation for nuclear genes
- aligning classification with the canonical curated schema
- retaining all rows with explicit classification status
- normalizing Parquet-loaded values for stable export
- reducing false unresolved nuclear classifications dramatically

The pipeline is now much closer to a scientifically defensible curation and
classification framework for both mtDNA and nucDNA cDAV/uDAV analysis.

---

## 6. Exception-handling implementation

An explicit exception-handling framework was added so that residual problematic
rows are handled deterministically rather than by hidden one-off code paths.

### 6.1 Exception registry

A registry file is now maintained at:

- `data/derived/reference/variant_exception_registry.tsv`
- compatibility copy: `data/reference/variant_exception_registry.tsv`

The registry records:

- scope (`gene` or `variant`)
- target gene or variant
- exception class
- exception code
- decision
- manual review status
- optional replacement transcript/model fields
- rationale and review notes

### 6.2 Current implemented exception classes

- `transcript_model_incompatible`
- `gene_specific_translation_conflict`
- `classified_with_warning`

### 6.3 Current implemented decisions

- `keep_unresolved`
- `classify_with_warning`

### 6.4 Pipeline integration

The exception registry is now loaded during curation and its metadata is
attached to curated rows using fields such as:

- `exception_scope`
- `exception_class`
- `exception_code`
- `exception_decision`
- `manual_review_status`
- `replacement_nm`
- `replacement_enst`
- `rescue_method`
- `exception_rationale`
- `exception_notes`

Classification now carries these fields through and also annotates:

- `exception_applied`
- `classification_exception_action`

### 6.5 Exception audit export

A dedicated audit entrypoint now exists:

- [01_audit_exception_candidates.py](/Users/ad2347/Documents/OxPhos_DAV/src/classify/01_audit_exception_candidates.py)

It exports:

- `data/derived/classified/exception_candidate_rows.tsv`
- `data/derived/classified/exception_candidate_summary.tsv`
- `data/derived/classified/exception_candidate_audit_metadata.json`

These tables combine:

- unresolved nuclear rows
- classified rows carrying `REF_ALLELE_MISMATCH`

and overlay the current exception-registry assignments for manual review.

---

## 7. Focused audit of the remaining nuclear unresolved and mismatch-heavy sets

This section summarizes the focused audit of the remaining unresolved nuclear
rows plus the residual `REF_ALLELE_MISMATCH` and `TRANSCRIPT_MISMATCH` sets
after the transcript-map rebuild and classify refactor.

### 6.1 Residual unresolved set

Current unresolved nuclear rows: `68`

Breakdown:

- `GENOMIC_POS_NOT_IN_ENST`: `45`
- `TRANSCRIPT_MISMATCH`: `21`
- `POSITION_NOT_IN_ENST`: `2`

Gene-level breakdown:

- `NDUFS6`: `12` `GENOMIC_POS_NOT_IN_ENST`
- `NDUFA13`: `10` `GENOMIC_POS_NOT_IN_ENST`
- `NDUFV2`: `10` `TRANSCRIPT_MISMATCH`
- `NDUFA11`: `7` `GENOMIC_POS_NOT_IN_ENST`
- `UQCRB`: `5` `GENOMIC_POS_NOT_IN_ENST`
- `NDUFS7`: `5` `GENOMIC_POS_NOT_IN_ENST`
- `NDUFB1`: `5` `TRANSCRIPT_MISMATCH`
- `NDUFA10`: `3` `TRANSCRIPT_MISMATCH`
- `ATP5MC2`: `3` `TRANSCRIPT_MISMATCH`
- `COXFA4L2`: `3` `GENOMIC_POS_NOT_IN_ENST`
- `NDUFA10`: `2` `POSITION_NOT_IN_ENST`
- `ATP5PF`: `2` `GENOMIC_POS_NOT_IN_ENST`
- `COX5A`: `1` `GENOMIC_POS_NOT_IN_ENST`

### 6.2 Audit interpretation by category

#### A. Low-identity genomic-rescue genes

These genes still route through genomic rescue because their transcript-map
quality is lower than the threshold used for transcript-first classification:

- `NDUFS6`
- `NDUFA13`
- `NDUFA11`
- `UQCRB`
- `NDUFS7`
- `COXFA4L2`
- `ATP5PF`
- `COX5A`

Representative transcript-map quality values:

- `NDUFS6`: identity `0.871`, coverage `1.0`
- `NDUFA13`: identity `0.792`, coverage `0.854`
- `NDUFA11`: identity `0.816`, coverage `1.0`
- `UQCRB`: identity `0.838`, coverage `1.0`
- `NDUFS7`: identity `0.887`, coverage `0.953`
- `COXFA4L2`: identity `0.828`, coverage `0.966`
- `ATP5PF`: identity `0.898`, coverage `1.0`
- `COX5A`: identity `0.787`, coverage `0.960`

**Pattern**

For these genes, the unresolved rows are almost entirely
`GENOMIC_POS_NOT_IN_ENST`, not generalized parsing or anchor failures.

This means:

- the pipeline is finding the right source row
- the genomic coordinate is being interpreted consistently
- but the ClinVar genomic position still does not land in the CDS model used by
  the selected TOGA ENST rescue map

**Interpretation**

These residual failures are most consistent with genuine transcript-model
discordance rather than a generic pipeline bug.

This is especially likely for:

- `NDUFA13`
- `NDUFA11`
- `UQCRB`
- `COXFA4L2`
- `COX5A`

where the identity or coverage values indicate only moderate compatibility
between the selected NM protein and the TOGA ENST protein.

**Recommended next step**

- review the exact TOGA human ENST and ClinVar/MANE transcript model for these
  genes
- confirm whether an alternate human ENST in TOGA is biologically closer to the
  MANE/ClinVar protein
- if no better ENST exists, keep these rows unresolved and document them as
  transcript-model incompatibilities

#### B. `NDUFV2` transcript mismatches

`NDUFV2` contributes `10` unresolved rows, all `TRANSCRIPT_MISMATCH`.

Transcript-map quality is relatively high:

- identity `0.936`
- coverage `1.0`

Representative failures include:

- `A6G`: corrected codon translates to `R`, expected `G`
- `R10L`: corrected codon translates to `I`, expected `L`
- `H17Q`: corrected codon translates to `*`, expected `Q`

Several also carry ref-allele disagreement in the warning log.

**Pattern**

- the position mapping is stable
- the mismatch happens at the codon/consequence level
- several events produce stop or unrelated amino-acid outcomes

**Interpretation**

This does not look like broad positional drift. It looks more like one of:

- ClinVar consequence annotation attached to a transcript that is not
  equivalent to the TOGA human coding model in the local N-terminal region
- incorrect reference nucleotide representation for these rows in the curated
  fields
- a transcript-version or exon-choice discrepancy not resolved by the current
  NM-to-ENST protein map alone

**Audit conclusion**

`NDUFV2` should be treated as a manual review target, not automatically
converted into classified cDAV/uDAV calls.

#### C. `NDUFB1` transcript mismatches

`NDUFB1` contributes `5` unresolved rows, all `TRANSCRIPT_MISMATCH`.

Transcript-map quality is nominally perfect:

- identity `1.0`
- coverage `1.0`

However, the mapped coordinates show a large N-terminal offset:

- raw CDS positions such as `10`, `16`, `19`, `29`, `40`
- corrected CDS positions around `151` to `181`

Representative failures include:

- `H6Y`: corrected codon translates to `*`, expected `Y`
- `W4R`: corrected codon translates to `I`, expected `R`
- `P10L`: corrected codon translates to `V`, expected `L`

**Interpretation**

This pattern is too coherent to be random noise. The most likely explanation is
that the selected TOGA ENST carries a long N-terminal extension or alternate
start context relative to the NM transcript used by ClinVar, and the current
protein-based position map is not enough to reconcile the coding frame
interpretation for these N-terminal substitutions.

**Audit conclusion**

`NDUFB1` should be treated as a transcript-model/translation-context exception.
It likely needs a gene-specific rule or manual exclusion from the automatic
comparative branch.

#### D. `NDUFA10` edge cases

`NDUFA10` contributes:

- `2` `POSITION_NOT_IN_ENST`
- `3` `TRANSCRIPT_MISMATCH`
- `1` classified row with `REF_ALLELE_MISMATCH`

Transcript-map quality is fairly strong:

- identity `0.927`
- coverage `0.972`

The two `POSITION_NOT_IN_ENST` rows are both at amino-acid position `217`:

- `R217Q`
- `R217W`

The `TRANSCRIPT_MISMATCH` rows are much farther downstream:

- `R337C`
- `K350M`
- `K350Q`

**Interpretation**

This suggests `NDUFA10` contains at least two distinct residual problems:

- a local region near residue `217` that is absent or gap-mapped in the TOGA
  ENST model
- a more distal region where transcript-projected codons no longer match the
  curated consequence

**Audit conclusion**

`NDUFA10` is not a global map failure, but a region-specific discordance case.
It merits a manual transcript/model comparison around residues `217` and
`337-350`.

#### E. `ATP5MC2` and `ATP5MF` mismatches despite high map quality

`ATP5MC2` contributes `3` `TRANSCRIPT_MISMATCH` rows with:

- identity `1.0`
- coverage `1.0`

Representative events:

- `R34G`
- `S29Y`
- `C35F`

All translate to incorrect codons after correction, despite a high-quality map.

`ATP5MF` contributes `1` classified `REF_ALLELE_MISMATCH` row:

- `C7F`

with:

- identity `0.936`
- coverage `1.0`

**Interpretation**

For `ATP5MC2`, this looks less like transcript-map failure and more like an
issue in the source allele/consequence representation, local CDS framing, or
reference-base normalization for these specific rows.

For `ATP5MF`, the row is still classifiable but carries a reference-base
warning, so it should be treated as lower-confidence rather than unresolved.

### 6.3 Residual `REF_ALLELE_MISMATCH` set

After the final classify rerun, the classified nuclear rows carrying
`REF_ALLELE_MISMATCH` are only:

- `ClinVar:897590:NDUFA10` (`R337H`)
- `ClinVar:3972444:ATP5MF` (`C7F`)

These rows were still classifiable, but the human aligned reference base did
not match the curated expected base at the corrected CDS position.

**Interpretation**

These should remain flagged in downstream analyses as classified-but-warning
rows, not treated identically to cleanly matched classified rows.

### 6.4 Audit summary

The remaining unresolved and mismatch-heavy nuclear cases fall into three broad
classes:

1. **Low-identity or incomplete transcript-model concordance**
   affecting genomic-rescue genes such as `NDUFA13`, `NDUFA11`, `UQCRB`,
   `NDUFS7`, `COXFA4L2`, `ATP5PF`, and `COX5A`.

2. **Gene-specific transcript/consequence discordance despite decent mapping**
   especially `NDUFV2`, `NDUFB1`, and `NDUFA10`.

3. **A small number of classified but warning-bearing rows**
   such as `NDUFA10 R337H` and `ATP5MF C7F`.

These no longer look like broad pipeline failures. They now look like the kind
of narrow residual cases that should be handled by targeted exception review.

---

## 8. Mitigations and safeguards used in the pipeline

This section records the main mitigations used to reduce false classification
while preserving traceability.

### 8.1 Filter-late retention

Mitigation:

- rows are retained with explicit status rather than dropped

Purpose:

- preserves denominators
- separates ineligibility from unresolved biology

### 8.2 Overlap-aware mtDNA duplication

Mitigation:

- overlap loci such as `MT-ATP6/MT-ATP8` are duplicated into frame-specific
  rows

Purpose:

- avoids collapsing distinct biological interpretations

### 8.3 Transcript-map-first coordinate rescue

Mitigation:

- classification first uses deterministic NM-to-ENST transcript maps when they
  are available and of acceptable quality

Purpose:

- avoids local heuristic coordinate guessing when explicit mapping exists

### 8.4 Genomic rescue for poor transcript-concordance genes

Mitigation:

- genomic rescue is used preferentially only when transcript-map quality is
  poor or insufficient

Purpose:

- provides a deterministic rescue path for transcript-discordant genes
- avoids forcing all genes through genomic rescue unnecessarily

### 8.5 Anchor fallback

Mitigation:

- when no transcript map is available, the classifier can fall back to a local
  amino-acid anchor search in the TOGA or mtDNA alignment
- this is implemented in
  [alignment_parser.py](/Users/ad2347/Documents/OxPhos_DAV/src/utils/alignment_parser.py)
  using a strict coordinate check followed by a local `+/-10` residue search
  for the expected wild-type amino acid

Purpose:

- rescues coordinate offsets when a direct map does not exist
- avoids assuming the source amino-acid position is directly aligned in the
  comparative reference

Guardrails:

- if the expected wild-type residue is not found, the row remains unresolved
- if corrected coding coordinates cannot produce a codon, the row remains
  unresolved
- if the corrected codon does not reproduce the curated consequence, the row is
  labeled `TRANSCRIPT_MISMATCH`

### 8.6 Explicit mismatch categories

Mitigation:

- unresolved states are emitted as explicit mismatch codes rather than folded
  into uDAV calls

Purpose:

- prevents unresolved rows from being misinterpreted as uncompensated
- supports targeted follow-up review

### 8.7 Classified-with-warning retention

Mitigation:

- rows with `REF_ALLELE_MISMATCH` can remain classified if the broader
  classification succeeds

Purpose:

- preserves potentially useful calls while making their warning state explicit

---

## 9. Manual review steps required

The remaining exception classes cannot be resolved safely by automation alone.
The following manual review steps should be followed and recorded in the
exception registry.

### 9.1 Required workflow for each reviewed gene or variant

1. Open the exception candidate row in
   `data/derived/classified/exception_candidate_rows.tsv`.
2. Confirm the current registry assignment in
   `data/derived/reference/variant_exception_registry.tsv`.
3. Compare the ClinVar `HGVS c.` and `HGVS p.` annotation against:
   - the selected NM transcript
   - the MANE transcript if different
   - the TOGA human ENST used by the alignment
4. Check whether the TOGA Homo sapiens amino-acid sequence is biologically
   equivalent to the NM/MANE protein in the local region of interest.
5. If genomic rescue is involved, confirm whether the ClinVar genomic position
   should lie inside the CDS of the selected TOGA human ENST.
6. Decide one of:
   - keep unresolved
   - classify with warning
   - adopt alternate transcript/model
7. Record the decision in the exception registry with:
   - `manual_review_status`
   - `replacement_nm` or `replacement_enst` if used
   - `rationale`
   - `reviewed_by`
   - `review_date`
8. Rerun curation and classification after registry updates.

### 9.2 Priority manual review list

Highest priority:

- `NDUFV2`
- `NDUFB1`
- `NDUFA10`

Second priority:

- `NDUFA13`
- `NDUFA11`
- `UQCRB`
- `NDUFS6`
- `NDUFS7`

Third priority:

- `COXFA4L2`
- `ATP5PF`
- `COX5A`
- `ATP5MC2`
- `ATP5MF`

### 9.3 Manual decisions that should not be made without evidence

Do not:

- widen anchor fallback heuristics without recording the change
- convert unresolved rows to uDAV by default
- overwrite curated amino-acid consequences with alignment-derived consequences
- add silent gene-specific overrides in code without a registry entry
- discard rows instead of recording an explicit decision

---

## 10. Manuscript-oriented summary

The curation and classification framework was redesigned to enforce a
metadata-first, filter-late analysis strategy in which all downloaded
disease-variant source records are retained in inventory and only later marked
for eligibility within the comparative cDAV/uDAV branch. A major conceptual
correction was the adoption of frame-specific interpretation for overlapping
mtDNA loci, especially `MT-ATP6/MT-ATP8`, so that a single nucleotide event can
be represented as two distinct biological interpretations when appropriate.

The largest technical correction concerned the nuclear transcript-reconciliation
layer. During audit, a substantial fraction of initially unresolved nucDNA
variants was traced not to genuine biological incompatibility, but to an error
in the transcript-position map builder that truncated multi-line TOGA human
protein sequences. After rebuilding the NM-to-ENST maps correctly and aligning
classification to the canonical curated master table, nuclear unresolved rows
fell from `3559` to `68`, while classified nuclear rows increased from `2995`
to `6486`. This substantially increased the recoverable nucDNA cDAV set and
showed that most early unresolved calls reflected transcript-mapping failure
rather than true absence of interpretable comparative evidence.

After these corrections, the remaining unresolved nuclear variants are
concentrated in a small set of genes with either weak NM-to-TOGA concordance or
gene-specific transcript/consequence discordance, notably `NDUFS6`, `NDUFA13`,
`NDUFA11`, `UQCRB`, `NDUFS7`, `NDUFV2`, `NDUFB1`, and `NDUFA10`. These cases
now represent explicit, reviewable incompatibilities between clinical variant
annotations and the comparative reference framework, rather than hidden
pipeline-level inconsistency. The current framework therefore supports
comparative classification of mtDNA and nucDNA cDAVs with substantially improved
traceability, overlap awareness, and transcript-model rigor.

---

## 11. Schema nuances added during stage-level QC and cross-source linkage

### 11.1 `classification_basis` and missense-only eligibility

The classification schema now distinguishes:

- `nt_and_aa`
- `nt_only`
- `aa_only`
- `no_disease_allele_detected`

An implementation bug had previously allowed any nucleotide-positive row to be
labeled `nt_and_aa`, even if amino-acid support was absent. This was corrected
when stage-level QC was added.

Important scientific interpretation:

- the core comparative branch still remains missense-focused
- rows are only eligible for the core branch if the curated consequence is
  missense
- this does **not** imply that amino-acid and nucleotide cDAV support must
  always coincide in non-human species

Why `nt_only` is still needed in the schema:

- a human disease variant can be missense in the human reference background
- the same disease nucleotide can appear in another species at the homologous
  nucleotide position
- but due to codon-context differences in that species, the resulting amino-acid
  state need not match the human disease amino-acid state

So `nt_only` is a valid comparative outcome even under a missense-focused
pipeline. It is not a synonym/silent-mutation category. It is a discordance
between nucleotide-level and amino-acid-level recurrence across species.

Current audited state:

- there are currently no `nt_only` classified rows in the dataset
- the schema still needs to represent that category correctly so future data or
  edge cases are not mislabeled

### 11.2 ClinVar mitochondrial rows are retained but not used for core cDAV/uDAV classification

The curated inventory retains mitochondrial-encoded `ClinVar` rows for
traceability under the Scope-B policy. However:

- they are **not** used as eligible rows in the nucDNA comparative branch
- they are marked ineligible with
  `core_pipeline_exclusion_reason = mt_in_clinvar_nuclear_branch`
- stage-level QC now explicitly verifies that no mitochondrial-encoded
  `ClinVar` row is eligible for the nuclear comparative branch

This means:

- mitochondrial `ClinVar` rows are present in the curated and classified master
  tables as inventory rows
- but they are not part of the active cDAV/uDAV comparative set used for the
  mtDNA or nucDNA branch outputs
- they remain useful for provenance, cross-database overlap auditing, and
  manuscript reporting

### 11.3 Database-specific nuance: mtDNA variants occur in ClinVar as well as MITOMAP

The disease-variant sources are not cleanly partitioned by genome:

- `MITOMAP` is the primary mtDNA disease-variant source for the active mtDNA
  comparative branch
- `ClinVar` also contains mitochondrial-encoded rows

These mitochondrial `ClinVar` rows are therefore handled as:

- retained source inventory rows
- excluded from the active nuclear comparative branch
- linked against `MITOMAP` using genomic mitochondrial nucleotide identity

Current audited cross-source counts:

- mitochondrial `ClinVar` rows retained: `4582`
- unique mitochondrial `ClinVar` variant IDs: `2291`
- unique mitochondrial `ClinVar` variant IDs also present in `MITOMAP`: `234`
- unique mitochondrial `ClinVar` variant IDs absent from `MITOMAP`: `2057`

These counts are now materialized in:

- `data/derived/results/clinvar_mitomap_cross_source_overlap_summary.json`
- `data/derived/results/clinvar_mt_variants_present_in_mitomap.tsv`
- `data/derived/results/clinvar_mt_variants_absent_from_mitomap.tsv`
