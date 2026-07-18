# Structural Anchoring Refactor Notes

This note records the anchoring refactor applied to
`src/structural/00_map_davs_to_structure.py`.

## Why the refactor was needed

Before this change, most structure-space failures that occurred after chain
assignment and transcript reconciliation collapsed into a single status:

- `position_not_anchored`

That status mixed several biologically and technically different situations:

- a small local mismatch where the residue simply was not found in the default
  search window
- a larger offset that could reflect mature-protein processing or numbering
  differences
- a gap or unresolved segment in the built structure
- a local sequence-model conflict
- a reference position that was not represented in the chain-position map

Because all of these ended up in the same bucket, the audit could show which
genes were problematic, but not why.

## What was changed

The mapper now uses a staged anchoring approach.

### Stage 1. Default anchor search

The script first tries the existing conservative anchor search:

- direct coordinate
- then `±10` residues around the target position

If successful, the row is labeled as:

- `ok`
- or `isoform_offset_*` when the best anchor is shifted

### Stage 2. Extended anchor search

If the default search fails, the script performs a wider bounded search:

- `±80` residues

This search is not used silently. If it finds a candidate, the failure is now
classified as:

- `large_offset_candidate_<delta>`

rather than being hidden inside generic anchoring failure.

This is intended to identify genes that may need a mature-protein offset rule
or a model-specific rescue, not to auto-accept weak mappings.

## Conservative policy update

An initial permissive run showed that allowing the extended window globally
greatly increased the mapped set, but also introduced a large number of
large-offset rescues that were not yet biologically justified.

That behavior was not retained as the final policy.

Instead, extended-window rescue is now gated by an explicit registry:

- `data/reference/structural_anchor_exception_registry.tsv`

Current enabled genes:

- `SDHA`
- `SDHB`
- `SDHC`
- `SDHD`

These were chosen as the clearest first-pass mature-offset candidates from the
audit.

All other genes now behave conservatively:

- the script still records large-offset candidates diagnostically
- but those candidates are not converted into successful mappings unless the
  gene is explicitly enabled in the registry

This keeps the current structural output split into:

- high-confidence direct/default-window mappings
- explicit registry-approved extended-offset rescues
- unresolved large-offset candidates

### Stage 3. Structured failure labeling

If no anchor is found, the script now distinguishes several failure types:

- `unresolved_structure_segment_candidate`
- `possible_large_offset_or_gap`
- `unmapped_reference_position`
- `reference_sequence_conflict`
- `anchor_window_exhausted`

These all remain failure states, but they are now separated explicitly in the
audit outputs.

## New metadata emitted by the structural stage

The structural outputs now carry:

- `anchor_method`
- `anchor_failure_detail`

and the structural transcript audit now also records:

- `clinvar_transcript_id`
- `preferred_nm`
- `toga_enst`
- `transcript_map_type`
- `transcript_map_identity`
- `transcript_map_coverage`
- `transcript_reconciliation_status`
- `transcript_coord_delta`

## New interpretation

For nucDNA rows, the structural stage now exposes three separate layers of
reconciliation:

1. ClinVar / curated transcript consequence space
2. TOGA ENST human alignment space
3. structure-chain residue space

The anchoring refactor makes it easier to tell whether a failure is primarily:

- transcript / ENST projection failure
- structure-chain absence
- unresolved structure segment
- large-offset candidate
- local sequence conflict

## Intended next use

This refactor is not the end of structural rescue. It is a diagnostic step that
prepares the next pass.

The immediate next analyses enabled by this change are:

- quantify how many failures are true large-offset candidates
- identify genes with likely mature-protein processing offsets
- separate unresolved-structure failures from model-absence failures
- decide where a cautious anchor rescue rule is scientifically defensible
