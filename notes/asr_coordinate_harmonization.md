# ASR Coordinate Harmonization

**Date:** 2026-05-18
**Scope:** Stage 4/phylogenetic ASR consumers that use `data/phylo/ancestral_state_maps.json`.

## Problem

IQTree `.state` files report ancestral states using **1-based alignment-site IDs** from the gapped IQTree FASTA. Downstream biological outputs use **human ungapped protein positions**:

- `dar_aa_coord` in structural and partner tables
- `contact_refseq_pos` in structural and partner tables
- ASR lookups in branch co-occurrence, phylogenetic timing, Pagel prep, and PyEvolve condition sets

The old ASR parser mixed coordinate systems:

- internal/root states and branch changes were keyed by IQTree alignment columns
- leaf states were derived by counting ungapped residues per species

This was unsafe for gapped alignments because species-specific ungapped position 8 is not necessarily the same alignment column as human protein position 8. For example, `SDHC` human protein position 8 maps to IQTree alignment column 19, not alignment column 8.

## Decision

Use **human ungapped protein position** as the active ASR coordinate system. This matches the coordinate system used by `dar_aa_coord` and `contact_refseq_pos`.

Alignment columns where Homo sapiens has a gap, `X`, or stop are excluded from the active ASR map because they do not correspond to a human protein residue and cannot represent a DAV/contact position.

## Implementation

Updated `src/phylo/01_parse_ancestral_states.py` to harmonize IQTree alignment sites through the Homo sapiens gapped alignment from:

```text
data/phylo/iqtree_jobs/{gene}/{gene}.fasta
```

For each gene, the parser now builds:

- `protein_pos_to_alignment_site`
- `alignment_site_to_protein_pos`
- `n_alignment_sites`
- `n_human_protein_positions`
- `n_human_gap_alignment_sites_dropped`
- `coordinate_system = human_protein_position`

The parser then:

1. Reads IQTree `.state` files using IQTree alignment-site IDs.
2. Converts internal-node and root states from alignment-site IDs to human protein positions.
3. Drops IQTree alignment sites where Homo sapiens has no residue.
4. Reads leaf states from the same human-mapped alignment columns, not by species-specific ungapped residue counts.
5. Builds `branches` and `branches_lc` after coordinate conversion, so branch changes are keyed by human protein position.
6. Preserves existing active field names for downstream compatibility:
   - `root_states`
   - `leaf_states`
   - `branches`
   - `branches_lc`

Added compact audit output:

```text
data/phylo/asr_coordinate_harmonization_audit.tsv
```

## Regression Checks

Real-gene spot checks after regeneration:

| Gene | Human protein position | IQTree alignment site | Human leaf AA |
|---|---:|---:|---|
| SDHC | 8 | 19 | H |
| SDHC | 56 | 74 | I |
| SDHC | 149 | 168 | Y |
| MT-ATP6 | 192 | 195 | I |
| MT-ND5 | 398 | 422 | T |

Synthetic unit tests were added in:

```text
tests/test_asr_coordinate_harmonization.py
```

They verify that:

- human protein position 8 can map to alignment column 19
- internal/root states at alignment column 19 appear as `root_states["8"]`
- leaf states are read from the human alignment column, not species ungapped residue 8
- branch changes are keyed by human protein position after conversion

Focused tests passed:

```text
tests/test_asr_coordinate_harmonization.py
tests/test_phylo_timing.py
```

## Outputs Regenerated Locally

The following outputs were regenerated after harmonization:

```text
data/phylo/ancestral_state_maps.json
data/phylo/asr_coordinate_harmonization_audit.tsv
results/structural/all_tested_pairs.csv
results/structural/concordance_summary.csv
data/phylo/pagel_jobs/
scripts/slurm_pagel_array.sh
results/phylo/timing_annotations.csv
results/phylo/contact_first_revised_test.csv
results/phylo/multi_origin_binomial.csv
```

The local structural partner rerun skipped local Pagel because the local Pagel runtime/tree integration was not available in that mode. Therefore:

- `all_tested_pairs.csv` has refreshed Fisher and branch-cooccurrence fields
- `compensatory_partners.csv` is not final until Pagel is rerun on HPC and merged

## Downstream Reruns Required

Treat previous Pagel, PyEvolve, timing interpretations, and combined phylogenetic summaries as stale until these are rerun from the harmonized ASR map.

Required next reruns:

1. Transfer regenerated `data/phylo/pagel_jobs/`, `src/`, `setup.py`, and `scripts/slurm_pagel_array.sh` to `lambda`.
2. Rerun Pagel SLURM jobs.
3. Copy Pagel results back to `data/phylo/pagel_results/`.
4. Run `src/phylo/merge_pagel_results.py`.
5. Rerun full PyEvolve conditional permissiveness jobs on HPC.
6. Merge PyEvolve chunks.
7. Rerun compartment, mito-nuclear, and combined summary/report scripts that depend on Pagel or PyEvolve significance.

DCA/APC-MI does not need recomputation solely because of ASR harmonization. It uses its own alignment mapping. However, any summary combining DCA/APC-MI with revised phylogenetic significance must be regenerated.

## Scientific Interpretation

After this fix, ASR positions used in branch co-occurrence, timing, Pagel prep, and PyEvolve condition sets refer to the same biological coordinate as the structural and variant tables: **human ungapped protein position**.

This makes the following comparisons valid:

- `dar_aa_coord` ↔ `root_states` / `leaf_states` / `branches`
- `contact_refseq_pos` ↔ `root_states` / `leaf_states` / `branches`
- branch-level DAV gain/loss events ↔ contact-site state changes
- timing categories such as `contact_first`, `co_occurring`, and `contact_after`

Previous interpretations that assumed ASR site numbers were already RefSeq/human protein positions should be considered superseded.
