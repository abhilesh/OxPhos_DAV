# Project Contract

## Goal
Study compensated and uncompensated disease-associated variation in OXPHOS genes, using human disease-causing mutations to investigate functional and evolutionary mechanisms across mitochondrial and nuclear genomes.

## Scientific Scope
- Core analysis scope: missense, disease-causing mutations in OXPHOS genes
- Traceability scope: retain all downloaded source records in inventory, even when later excluded from the core comparative branch

## Current Phase Boundary
The current work covers `data_download` and `data_curation`:
- `data_download`: downloading, inventorying, validating, and normalizing raw external inputs
- `data_curation`: generating canonical filter-late curated products, transcript/genomic support maps, and alignment sanitation metadata
- preserving current downloaded files and existing expensive alignments as the active analysis inputs

## Non-Negotiable Rules
- Do not recompute MACSE or other expensive alignment steps unless explicitly approved
- Keep current dated download files as the source of truth during this phase
- Do not modify or regenerate expensive existing alignments; only annotate or sanitize them non-destructively
- Record exclusions, uncertainty, and validation results as metadata rather than silently dropping context
- Treat overlap loci, especially `MT-ATP6/MT-ATP8`, as frame-specific duplicated interpretations rather than collapsing to one gene

## Data Layout
- `data/raw/`: dated downloaded third-party assets
- `data/derived/reference/`: generated canonical reference products
- `data/derived/curated/`: canonical Parquet curation products, support-map products, and sanitation manifests
- legacy compatibility paths may remain live during transition

## Traceability Rule
All downloaded source records should be traceable through manifest metadata, validation metadata, and stable resource naming, even if they are not part of the core missense comparative branch later.

## Curation Rule
The canonical curated product is a Parquet master table with filter-late metadata. Exclusions, unresolved states, transcript/genomic rescue status, and overlap-specific frame interpretations must be recorded explicitly rather than encoded by record omission.
