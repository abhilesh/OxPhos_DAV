# Pipeline Consistency, Filter-Late Audit, and Manuscript-Ready Repo Plan

This document has three goals:

1. Audit whether the current pipeline excludes variants before analysis
2. Define a metadata-first, filter-late implementation target
3. Recommend a repository structure suitable for manuscript submission

The guiding principle is:

> Keep all in-scope variants in the working dataset as long as possible.  
> Record exclusions, mismatches, uncertainty, and filter criteria as metadata fields.  
> Apply biological or quality filters only in final analysis views and manuscript tables.

---

## 1. Executive Summary

The current codebase is **not yet fully filter-late**.

There are three layers of early exclusion:

- `Source parsing exclusions`: variants are dropped before they ever become `VariantRecord`s
- `Classification exclusions`: some curated variants are skipped entirely instead of being retained with unresolved status
- `Downstream analysis exclusions`: some analyses only load non-discarded cDAVs rather than all variants with metadata

At the moment:

- `01_curate_variants.py` is **not metadata-complete for all raw records**, because its parsers already drop several record classes
- `00_classify_DAV.py` is **not fully filter-late**, because synonymous and `Discarded` variants are skipped rather than retained with classification metadata
- `00_map_davs_to_structure.py` is **closest to filter-late**, because it loads all records from the classification outputs, but it still silently drops records with unparsable amino-acid coordinates
- `01_find_compensating_partners.py` is intentionally **cDAV-only**, which is fine as a derived analysis step, but it should not be the first point where records disappear

If your goal is to run the analysis on all variants first and defer filtering to the end, the main redesign should be:

- curate everything in scope
- classify everything in scope
- map everything in scope
- annotate every failure as metadata
- derive filtered analysis tables only at the final summarization stage

---

## 2. Filter Audit by Pipeline Stage

## 2.1 `src/data_prep/01_curate_variants.py`

### Current behavior

`01_curate_variants.py` itself mostly assembles metadata and does not add much filtering logic.  
The bigger issue is that the parser classes it depends on already exclude records before curation.

### Upstream exclusions happening before curation

#### MITOMAP parser

In [src/utils/parsers.py](/Users/ad2347/Documents/OxPhos_DAV/src/utils/parsers.py:158), the parser drops:

- all mtDNA variants outside OXPHOS loci
- multi-nucleotide variants
- noncoding variants
- frameshifts
- stop-gain / `Ter` variants
- variants with unparseable amino-acid strings

This means the curated mtDNA dataset is already restricted to **single-nucleotide coding substitutions with parseable AA changes**.

#### ClinVar parser

In [src/utils/parsers.py](/Users/ad2347/Documents/OxPhos_DAV/src/utils/parsers.py:300), the parser drops:

- all non-GRCh38 entries
- all non-SNV entries
- all MT chromosome entries
- all genes not in the HGNC lookup
- all entries without parseable protein changes
- all entries whose protein change parses to `X`

So the curated nucDNA dataset is already restricted to **GRCh38 nuclear SNVs with parseable amino-acid substitutions in target genes**.

### Interpretation

This is acceptable only if your declared scope is:

- OXPHOS-associated coding SNVs
- with interpretable amino-acid consequences

It is **not** acceptable if “all variants” means:

- all ClinVar OXPHOS variants
- all MITOMAP coding variants
- indels
- nonsense
- frameshift
- splice
- noncoding

### Recommendation

Decide explicitly between these two scopes:

### Scope A: amino-acid comparative pipeline

If this project is fundamentally about amino-acid and codon-aware compensation, define the core scope as:

- coding OXPHOS variants
- single-nucleotide substitutions
- amino-acid interpretable variants

Then keep excluded classes in a separate raw inventory table with metadata fields:

- `excluded_from_core_pipeline = True`
- `exclusion_stage = parser`
- `exclusion_reason = non_snv | frameshift | nonsense | noncoding | unparsable_aa | non_oxphos`

### Scope B: all clinical variants inventory

If you want true whole-dataset bookkeeping, then parsers should not drop these records entirely.  
They should emit them into a broader master table and mark whether they are eligible for the codon-aware comparative branch.

---

## 2.2 `src/classify/00_classify_DAV.py`

### Current behavior

The classifier excludes some records before classification:

In [src/classify/00_classify_DAV.py](/Users/ad2347/Documents/OxPhos_DAV/src/classify/00_classify_DAV.py:244), it does:

- `if "Discarded" in tier or var["is_synonymous"]: continue`

So these records never enter the classification output JSONs at all.

### Why this conflicts with your goal

If you want all filtering deferred until the last steps, then:

- synonymous variants should still receive classification-stage metadata
- discarded variants should still survive as rows with explicit reasons

Right now they disappear from the main classified dataset.

### Good part of the current implementation

For many other failures, the script behaves the way you want:

- missing alignment -> record retained with `None` classification fields
- anchor failure -> retained
- transcript mismatch -> retained
- codon extraction failure -> retained
- genomic position not in ENST -> retained

That is already metadata-first behavior for unresolved-but-in-scope cases.

### Recommendation

Refactor classification so that it never drops records after curation.

Instead of skipping synonymous and discarded variants, keep them and annotate:

- `classification_eligible = False`
- `classification_exclusion_reason = synonymous | discarded_tier`
- `is_cdav_amino_acid = None`
- `is_cdav_nucleotide = None`

Also add:

- `classification_status = classified | unresolved | skipped_by_policy`

This would make the classification outputs your true master analysis-ready tables.

---

## 2.3 `src/structural/00_map_davs_to_structure.py`

### Current behavior

This script is relatively close to your desired model.

It loads all records from the classification JSONs, not just cDAVs:
- [src/structural/00_map_davs_to_structure.py](/Users/ad2347/Documents/OxPhos_DAV/src/structural/00_map_davs_to_structure.py:487)

That means both cDAVs and uDAVs can be structurally mapped.

### Remaining early exclusion

It still silently skips records when:

- `aa_change` cannot be parsed into an amino-acid coordinate
- no complex ID is found
- no structure file exists

The first case is the one that most clearly violates metadata-first behavior:

- [src/structural/00_map_davs_to_structure.py](/Users/ad2347/Documents/OxPhos_DAV/src/structural/00_map_davs_to_structure.py:495)
- [src/structural/00_map_davs_to_structure.py](/Users/ad2347/Documents/OxPhos_DAV/src/structural/00_map_davs_to_structure.py:497)

### Recommendation

Even when `aa_coord` cannot be parsed, emit a mapping-status row:

- `mapping_status = aa_coord_unparseable`
- `structure_mapping_eligible = False`

The structural mapping table should contain one row per variant, even for failed mappings.

---

## 2.4 `src/structural/01_find_compensating_partners.py`

### Current behavior

This script explicitly loads only:

- non-discarded variants
- AA-level cDAVs

See [src/structural/01_find_compensating_partners.py](/Users/ad2347/Documents/OxPhos_DAV/src/structural/01_find_compensating_partners.py:356)

### Interpretation

This is acceptable only if this script is understood as a **derived cDAV-only analysis layer**.

It should not be treated as the place where the master dataset is defined.

### Recommendation

Keep this as a derived analysis, but:

- rename its inputs conceptually as “eligible cDAV subset”
- ensure the master dataset already exists upstream with all variants and metadata
- consider later adding a matched uDAV structural-neighbor table for direct cDAV vs uDAV comparisons

---

## 2.5 `src/data_prep/00f_build_transcript_position_maps.py`

### Current behavior

This script builds transcript maps only for genes that already have ClinVar-supported NM_ usage in curated data:

- [src/data_prep/00f_build_transcript_position_maps.py](/Users/ad2347/Documents/OxPhos_DAV/src/data_prep/00f_build_transcript_position_maps.py:191)

If no ClinVar NM_ is seen for a gene, the gene is skipped.

### Why this matters

This is not variant filtering directly, but it means coordinate support infrastructure is built only for the currently observed variant set, not for all target genes.

### Recommendation

Build transcript-position maps for all targeted nuclear OXPHOS genes, not only those already represented in current ClinVar records.

That makes the pipeline more stable for future updates and avoids data-dependent infrastructure.

---

## 2.6 `src/data_prep/00c_download_zoonomia_TOGA_alignments.py`

### Current behavior

Transcript preference is partially data-dependent:

- it uses the curated nucDNA annotations to identify the most common ClinVar NM_ per gene
- [src/data_prep/00c_download_zoonomia_TOGA_alignments.py](/Users/ad2347/Documents/OxPhos_DAV/src/data_prep/00c_download_zoonomia_TOGA_alignments.py:139)

### Why this matters

This is not immediate variant filtering, but it does make alignment selection depend on the current disease dataset.

### Recommendation

For manuscript-grade reproducibility, define transcript selection deterministically:

- primary rule: MANE Select when available
- explicit exception table for known TOGA/MANE mismatches
- current ClinVar NM_ usage retained as metadata, not as the main selection driver

---

## 3. What Needs to Change for a True Filter-Late Pipeline

## 3.1 Recommended data model

You should have one master row per variant that persists through the pipeline.

Suggested metadata groups:

### Identity and provenance
- `variant_id`
- `source_db`
- `source_record_id`
- `genome`
- `reference_assembly`
- `raw_source_file`
- `source_db_version`

### Core variant description
- `locus`
- `loci_raw`
- `overlap_locus_flag`
- `variant_class`
- `is_snv`
- `is_mnv`
- `is_indel`
- `is_coding`
- `is_missense`
- `is_synonymous`
- `is_nonsense`
- `is_frameshift`
- `is_splice_related`

### Scope flags
- `eligible_core_comparative_pipeline`
- `core_pipeline_exclusion_reason`
- `classification_eligible`
- `structure_mapping_eligible`
- `partner_analysis_eligible`

### Coordinate resolution
- `hgvs_c`
- `hgvs_p`
- `transcript_id`
- `genomic_pos`
- `coordinate_resolution_method = raw | tx_map | genomic_map | anchor | unresolved`
- `coordinate_resolution_status`

### Classification outputs
- `classification_status = classified | unresolved | skipped_by_policy`
- `is_cdav_amino_acid`
- `is_cdav_nucleotide`
- `n_species_aligned`
- `n_species_with_disease_allele`
- `lineages_with_disease_allele`
- `ref_allele_match`
- `mismatch_reason`

### Structural outputs
- `structure_mapping_status`
- `pdb_id`
- `pdb_chain`
- `pdb_resnum`
- `proxy_mapping_used`
- `n_structural_contacts`

### Filtering metadata for late-stage analysis
- `tier`
- `is_sdh`
- `is_haplogroup_marker`
- `clinvar_stars`
- `mitomap_pubmed_count`
- `gnomad_af_global`
- `alphamissense_score`
- `apogee2_score`
- any final analysis inclusion flags should be derived, not destructive

---

## 3.2 Recommended implementation changes

### Curation stage

Do not silently drop in-scope but non-core records.

Create:
- `variants_master_raw.jsonl` or `.parquet`
- `variants_master_curated.parquet`

Every raw record should either:
- become a curated record, or
- become a curated record with `eligible_core_comparative_pipeline = False`

### Classification stage

Never skip rows after loading curated variants.

Instead:
- annotate each row with status fields
- preserve unresolved and ineligible rows
- write one full table, not only a biologically filtered subset

### Structural stage

Emit a mapping-status row for every variant loaded from the classification table.

### Derived analysis stage

Only here should you create filtered views such as:
- `cdav_only`
- `non_discarded`
- `tier_1_2_only`
- `sdh_excluded`

These should be views or derived files, not replacements for the master dataset.

---

## 4. File-by-File Implementation Checklist

## Must refactor now

### `src/utils/parsers.py`
- Add exclusion metadata instead of dropping non-core records entirely, if you want full inventory coverage
- At minimum, produce a sidecar exclusions table during parsing

### `src/data_prep/01_curate_variants.py`
- Write a master curated table with both eligible and ineligible records
- Prefer columnar format such as Parquet in addition to JSON

### `src/classify/00_classify_DAV.py`
- Remove hard skip for `Discarded` and `is_synonymous`
- Convert policy skips into metadata
- Preserve one row per input variant

### `src/structural/00_map_davs_to_structure.py`
- Emit failure rows for unparsable AA coordinates
- Preserve one row per classified variant

## Strongly recommended next

### `src/data_prep/00f_build_transcript_position_maps.py`
- Build maps for all target genes, not only currently represented genes

### `src/data_prep/00c_download_zoonomia_TOGA_alignments.py`
- Make transcript selection rules deterministic and explicit

### `src/structural/01_find_compensating_partners.py`
- Keep current cDAV-only logic, but label it clearly as a derived view
- Consider future uDAV-matched comparison tables

---

## 5. Recommended Repository Organization for Manuscript Submission

Your current repo is workable for development, but for manuscript submission I would move toward a clearer separation of:

- immutable inputs
- reproducible intermediate products
- final analysis tables
- manuscript artifacts

## 5.1 Suggested top-level layout

```text
OxPhos_DAV/
  README.md
  environment/
    Dockerfile
    environment.yml
    requirements.txt
  config/
    pipeline.yaml
    paths.yaml
    thresholds.yaml
    manuscript_cohorts.yaml
  data/
    raw/
      annotations/
      alignments/
      structures/
    external/
      reference/
    derived/
      curated/
      classified/
      structural/
      phylo/
      mutagenesis/
  src/
    oxphos_dav/
      data_prep/
      align/
      classify/
      structural/
      phylo/
      mutagenesis/
      utils/
      schemas/
  workflows/
    Snakefile or justfile
    rules/
  results/
    figures/
    tables/
    supplementary/
  manuscript/
    main/
    supplementary/
    cover_letter/
  docs/
    methods/
    decisions/
    audits/
  tests/
```

## 5.2 Data product philosophy

Use three explicit classes of outputs:

### `data/raw`
- downloaded, immutable inputs
- never edited manually

### `data/derived`
- reproducible intermediate tables
- each stage writes versioned or well-named outputs
- these should include master tables, not only filtered subsets

### `results`
- final tables and figures used in the paper
- always derived from `data/derived`

This separation helps a lot for review, lab handoff, and future revisions.

## 5.3 Prefer tables over nested JSON for master products

For manuscript preparation, I strongly recommend:

- Parquet for master datasets
- CSV only for export tables
- JSON only for special nested structures such as ancestral-state maps

Suggested master tables:

- `data/derived/curated/variants_master.parquet`
- `data/derived/classified/variants_classified.parquet`
- `data/derived/structural/variants_structural_map.parquet`
- `data/derived/structural/variant_contact_pairs.parquet`

Then create analysis-specific views:

- `cdav_view.parquet`
- `udav_view.parquet`
- `tiered_variant_view.parquet`
- `mt_vs_nuc_matched_view.parquet`

## 5.4 Add explicit schema definitions

Create schema documents or dataclasses for:

- curated variant master table
- classified variant table
- structural mapping table
- partner table

This will help reviewers and future you understand what each file means.

## 5.5 Add a true workflow entry point

Right now the repo is script-driven. For manuscript submission, it will be much stronger if you have:

- a `Makefile`, `justfile`, or `Snakemake` workflow
- named stages
- one command to regenerate core outputs

For example:

```text
make curate
make classify
make structural
make phylo
make manuscript-tables
```

## 5.6 Separate exploratory notes from stable method docs

Right now `notes/` mixes:

- methodological caveats
- scientific interpretation
- development reasoning

For manuscript readiness, split this into:

- `docs/methods/` for stable methods docs
- `docs/audits/` for filter audits and validation
- `notes/` for active exploratory thinking only

## 5.7 Freeze manuscript analysis cohorts in config files

Instead of hard-coding tiers or exclusions in scripts, define them in config:

- `config/manuscript_cohorts.yaml`

This could include:

- full dataset
- non-SDH dataset
- tier 1-2 dataset
- matched-stringency mt vs nuc dataset

That makes the manuscript analysis both transparent and rerunnable.

---

## 6. Recommended Near-Term Work Plan

### Phase 1: make the pipeline metadata-first
- stop dropping synonymous and discarded variants at classification
- preserve one row per curated variant throughout classification and structural mapping
- create explicit exclusion/status fields

### Phase 2: stabilize infrastructure
- fix overlap-locus handling
- fix path/file consistency
- fix mtDNA downloader/schema mismatches
- make transcript-selection rules deterministic

### Phase 3: manuscript-ready outputs
- create master Parquet tables
- create derived cohort views
- add workflow entry points
- move stable docs into manuscript-oriented structure

---

## 7. Bottom Line

If your requirement is:

> “Run the analysis on all variants and defer filtering till the last steps”

then the current pipeline is not there yet.

The most important changes are:

- stop destructive skipping in classification
- preserve excluded/ineligible records as rows with metadata
- distinguish master datasets from filtered analysis views
- organize the repo around reproducible data products rather than only stage-specific scripts

With those changes, the pipeline will be much more scientifically defensible and much easier to present in a manuscript.

