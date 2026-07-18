# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Summary

This project identifies **compensated disease-associated residues (c-DARs)** in mammalian OXPHOS genes by finding non-human species that naturally harbor human-pathogenic amino acids. It combines clinical variant databases (MITOMAP for mtDNA, ClinVar for nucDNA), cross-species sequence alignments (TOGA), structural data from cryo-EM, and phylogenetic analysis.

- **AA-level c-DAR**: Any non-human species harboring the human-pathogenic amino acid.
- **NT-level c-DAR**: Stricter — the species uses the exact same codon as the human mutation.

## Environment

All analysis runs inside the Docker container `oxphos_dav_analysis`. Build and enter it:

```bash
docker build -t oxphos_dav_analysis .
docker run -it --rm -v $(pwd):/app oxphos_dav_analysis
```

The devcontainer (`.devcontainer/devcontainer.json`) automates this for VS Code. Inside the container, the `oxphos_dav` conda environment is active and the package is installed as editable (`pip install -e .`).

## Running Scripts

All pipeline scripts are run from the **project root** with `python src/<stage>/<script>.py`. `pytest.ini` sets `pythonpath = src` so tests resolve imports without the container.

## Tests

```bash
# Run all tests
pytest

# Run a single test file
pytest tests/test_cdav_invariants.py -v

# Run a single test
pytest tests/test_cdav_invariants.py::test_name -v
```

Tests in `tests/` validate pipeline invariants against the live Parquet outputs in `data/derived/`. Tests skip (not fail) when required Parquet files are missing. Legacy tests are in `archive/legacy_tests/` and are not run by default.

## Pipeline Stages (in order)

### 1. `src/data_download/` — Download raw data
Scripts `00a`–`00j` download annotations (MITOMAP, ClinVar, MitImpact, dbSNFP, gnomAD), reference data (TOGA, MANE, gene coords, aaindex, PhyloTree), and PDB structures. `00i_validate_downloads.py` checks integrity.

### 2. `src/data_curation/` — Curate and map variants
- `01_curate_variants.py` — produces `data/derived/curated/variants_master_curated.parquet` (the canonical variant table). Sources: MITOMAP + ClinVar. Cross-source duplicate groups are linked.
- `02_build_transcript_position_maps.py` — NM_ → ENST position maps (JSON).
- `03_build_genomic_coordinate_maps.py` — genomic coordinate lookup tables.
- `04_sanitize_all_alignments.py` — validates TOGA alignment FASTA files.

### 3. `src/align/` — Translate sequences
- `00_align_translate_mtDNA.py` — codon alignments for mtDNA genes.
- `00_translate_nucDNA.py` — protein alignments for nuclear genes.

### 4. `src/classify/` — Classify DAVs as cDAV or uDAV
- `00_classify_DAV.py` — reads curated Parquet, applies `AlignmentParser` per gene, writes `data/derived/classified/variants_master_classified.parquet` plus QC splits (`classified_clean`, `classified_warning`, `classified_all`). Filter-late: all rows are retained; ineligible rows are marked as skipped.
- `01_audit_exception_candidates.py` — identifies rows for the exception registry.

### 5. `src/structural/` — Map variants to cryo-EM structures
- `00_map_davs_to_structure.py` — maps classified variants onto PDB/CIF structures. Uses local alignment for chain→gene assignment and global alignment for position mapping. Contacts are defined as Cβ–Cβ ≤ 8 Å (Cα for Gly) and classified as hbond / electrostatic / hydrophobic / vdw. Outputs to `results/structural/`.
- `01_find_compensating_partners.py` — tests structural contacts for co-evolutionary enrichment using Fisher's exact, Pagel's discrete (R phytools), and branch co-occurrence (IQTree ancestral states). Outputs to `results/structural/compensatory_partners.csv`.

### 6. `src/phylo/` — Phylogenetic analysis
HPC-bound scripts (`slurm_*.sh` in `scripts/`). Uses IQTree for ancestral state reconstruction and VertLife mammal tree for Pagel's discrete test. Pagel jobs run via R (`pagel_discrete.R`). Results land in `results/phylo/`.

### 7. `src/mutagenesis/` — Mutagenesis target prioritization
Ranks 693 c-DAR/compensatory partner pairs by FoldX ΔΔG, mutual information (pydca/evcouplings), and physicochemical scores. Outputs to `results/mutagenesis/`.

## Key Data Flow

```
data/raw/annotations/        ← downloaded source files (MITOMAP, ClinVar, etc.)
data/raw/reference/          ← downloaded reference data (TOGA, MANE, etc.)
data/alignments/             ← TOGA FASTA alignments per gene (AA + NT)
data/derived/curated/        ← variants_master_curated.parquet (post-curation)
data/derived/classified/     ← variants_master_classified.parquet (post-classify)
results/structural/          ← structure mapping + compensatory partner outputs
results/phylo/               ← phylogenetic test outputs
results/mutagenesis/         ← mutagenesis prioritization outputs
```

## Key Utility Modules (`src/utils/`)

- `variant_record.py` — `VariantRecord` dataclass; central evidence schema with physicochemical properties (BLOSUM62, Miyata distance, KD hydrophobicity, volume). Fields for structural/cross-species data are `None` at curation time and populated downstream.
- `alignment_parser.py` — `AlignmentParser` wraps TOGA AA+NT FASTA per gene; exposes `check_compensation()` which returns AA-level and NT-level cDAV calls. Uses transcript position maps (`tx_pos_map`) when available; falls back to sequence anchoring.
- `exception_registry.py` — loads `data/reference/variant_exception_registry.tsv` and `data/reference/structural_anchor_exception_registry.tsv`. Exceptions can be scoped to `variant` or `gene`; variant-scope takes priority.
- `gene_reference.py` — canonical OXPHOS gene list and genome assignment (mtDNA vs nucDNA).
- `mt_overlap.py` — handles mtDNA gene overlaps (e.g., MT-ATP8/MT-ATP6).

## Reference Files (`data/reference/`)

| File | Purpose |
|------|---------|
| `Canonical_OXPHOS_Subunits_HGNC_*.csv` | Canonical gene list |
| `structure_model_manifest.tsv` | PDB ID → gene → chain mapping for structural stage |
| `structural_anchor_exception_registry.tsv` | Manual offsets for isoform mismatches in structure mapping |
| `variant_exception_registry.tsv` | Per-variant or per-gene classification overrides |
| `transcript_position_maps.json` / `genomic_coordinate_maps.json` | Built by curation stage; used by classify |
| `TOGA_overview_table_hg38_*.tsv` | TOGA ortholog quality scores |
| `aaindex_properties_*.json` | Physicochemical AA property matrices (loaded by `variant_record.py`) |

## Notes Directory

`notes/` contains the authoritative pipeline documentation. Read in this order for context:
1. `pipeline_consistency_and_manuscript_prep.md` — overall project contract
2. `resolved_nuances_and_fixes.md` — resolved upstream issues
3. `unresolved_issues_ranked.md` — current open issues (prioritized)
4. `curation_and_classification_issues_and_fixes.md` — curation/classify contract
5. `structural_phase_plan.md` + `structural_mapping_major_problems_report.md` — structural stage
6. `running_results.md` — current output counts and metrics
7. `final_analysis_filtering_decisions.md` — filtering flags that must be applied at the final analysis layer (not intermediate scripts)

## Obsidian Wiki Project Tracker

When running `/wiki-update` or syncing project knowledge, follow the owner-specific policy in:

`/Users/ad2347/Documents/Obsidian_Wiki/AGENTS.md`

In particular, update the project-management layer:
- `projects/<ProjectName>/status.md`
- `projects/<ProjectName>/tasks.md`
- `projects/<ProjectName>/tasks-done.md`

Do not create per-thread task pages unless explicitly requested. Merge durable knowledge into existing project pages, move completed/superseded tasks to `tasks-done.md`, refresh QMD after writing, and leave unrelated dirty files alone.
