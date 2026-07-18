# Running Results

This document records the important analysis results accumulated so far for the OXPHOS DAV study. It will be updated as additional stages of the pipeline are completed.

## 1. Current scope

- Current implemented stages: `data_download`, `data_curation`, and `classify`.
- Current targeted OXPHOS gene set from HGNC: `100` genes total (`13` mitochondrial, `87` nuclear).
- Current curated DAV inventory: `33779` curated interpretation rows.
- Current disease-variant sources entering curation: `MITOMAP` for mtDNA and `ClinVar` for nucDNA.

## 2. Data sources used so far

The current download and reference layer includes the following active source inputs:

| Resource | Download date | Local file | Source URL | Validation |
| --- | --- | --- | --- | --- |
| clinvar_variant_summary | 2026-03-26 | data/raw/annotations/ClinVar_VariantSummary_2026-03-26.txt.gz | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz | ok |
| hgnc_oxphos_gene_list | 2026-03-25 | data/raw/reference/Canonical_OXPHOS_Subunits_HGNC_2026-03-25.csv | https://www.genenames.org/cgi-bin/genegroup/download?id=639&type=branch | ok |
| mane_grch38_summary | 2026-04-18 | data/raw/reference/MANE_GRCh38_v1.5_2026-04-18.txt.gz | https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz | ok |
| mitimpact_db | 2026-03-26 | data/raw/annotations/MitImpact_db_2026-03-26.txt.zip | https://mitimpact.css-mendel.it/cdn/MitImpact_db_3.1.3.txt.zip | ok |
| mitomap_coding_variants | 2026-03-26 | data/raw/annotations/MITOMAP_CodingVariants_2026-03-26.tsv | https://www.mitomap.org/cgi-bin/disease.cgi | ok;remote_http_403 |
| myvariant_dbnsfp_gnomad | 2026-03-26 | data/raw/annotations/MyVariant_dbNSFP_gnomAD_2026-03-26.json | https://myvariant.info/v1/query | ok;completeness_legacy_unchecked |
| phylotree_build17 | 2026-03-26 | data/raw/reference/PhyloTree_build_17_2026-03-26.zip | https://www.phylotree.org/builds/mtDNA_tree_Build_17%20-%20rCRS-oriented%20version.zip | ok |
| toga_overview_hg38 | 2026-04-07 | data/raw/reference/TOGA_overview_table_hg38_2026-04-07.tsv | https://genome.senckenberg.de/download/TOGA/human_hg38_reference/overview.table.tsv | ok |

Additional comparative and reference assets currently present in the repo and used by the pipeline include:

- TOGA codon-aware and amino-acid alignments for nuclear genes under `data/alignments/toga_hg38_aa/` and `data/alignments/toga_hg38_codon/`.
- Existing mtDNA codon-aware and amino-acid alignments under `data/alignments/mtdna_codon/` and `data/alignments/mtdna_aa/`.
- Canonical transcript and genomic rescue maps under `data/derived/curated/`.
- Exception registry and focused audit products under `data/derived/reference/` and `data/derived/classified/`.
- Human OXPHOS structure models downloaded from the RCSB PDB under `data/structures/`, using a manifest-backed structure panel that includes primary and validation models.

Current structure panel sources:

| Complex | Role | PDB | URL |
| --- | --- | --- | --- |
| CI | primary | 9I4I | https://www.rcsb.org/structure/9I4I |
| CI | validation | 9TI4 | https://www.rcsb.org/structure/9TI4 |
| CI | reference | 5XTH | https://www.rcsb.org/structure/5XTH |
| CII | primary | 8GS8 | https://www.rcsb.org/structure/8GS8 |
| CIII | primary | 9HZL | https://www.rcsb.org/structure/9HZL |
| CIII | reference | 5XTE | https://www.rcsb.org/structure/5XTE |
| CIII | validation | 5XTH | https://www.rcsb.org/structure/5XTH |
| CIV | primary | 9I7U | https://www.rcsb.org/structure/9I7U |
| CIV | validation | 9I6F | https://www.rcsb.org/structure/9I6F |
| CIV | reference | 5Z62 | https://www.rcsb.org/structure/5Z62 |
| CIV | validation | 5XTH | https://www.rcsb.org/structure/5XTH |
| CV | primary | 8H9S | https://www.rcsb.org/structure/8H9S |
| CV | validation | 8H9T | https://www.rcsb.org/structure/8H9T |
| CV | validation | 8H9U | https://www.rcsb.org/structure/8H9U |

## 3. Downloaded genes and variant inventory

- HGNC targeted gene list: `100` OXPHOS genes total (`13` mtDNA-encoded, `87` nucDNA-encoded).
- Curated mtDNA interpretation rows from MITOMAP: `388` rows spanning `9` genes.
- Curated nucDNA interpretation rows from ClinVar: `33391` rows spanning `100` genes.
- Unique source-variant groups: `378` for mtDNA and `16779` for nucDNA.
- mtDNA curated rows are interpretation-level rows and include overlap-derived duplications where biologically required, especially for `MT-ATP6/MT-ATP8`.

## 4. Current cDAV classification results

The current classifier identifies compensated and uncompensated DAVs at both the amino-acid and nucleotide levels for mtDNA and nucDNA.

| Genome | Total rows | Classified | Eligible | Unresolved | Skipped by policy | AA-level cDAVs | NT-level cDAVs | AA cDAV % of classified | NT cDAV % of classified |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mtDNA | 388 | 323 | 323 | 0 | 65 | 167 | 149 | 51.7 | 46.13 |
| nucDNA | 33391 | 6486 | 6554 | 68 | 26837 | 2046 | 1745 | 31.54 | 26.9 |

Main current cDAV proportions using classified rows as the denominator:

- `mtDNA`: AA-level cDAVs `167/323` = `51.7%`; NT-level cDAVs `149/323` = `46.13%`.
- `nucDNA`: AA-level cDAVs `2046/6486` = `31.54%`; NT-level cDAVs `1745/6486` = `26.9%`.

At the current state of the pipeline, this corresponds to:

- `mtDNA`: `167` AA-level cDAVs and `149` NT-level cDAVs among `323` classified rows.
- `nucDNA`: `2046` AA-level cDAVs and `1745` NT-level cDAVs among `6486` classified rows.

### 4.1 Classified subsets and QC status

The classification stage now writes explicit derived subsets and stage-level QC outputs under `data/derived/classified/`.

Current subset counts:

- `classified_all`: `6809`
- `classified_clean`: `6340`
- `classified_warning`: `469`

Current warning composition:

- `467` rows are warning-classified because `exception_applied = True`
- `2` rows are warning-classified because both `exception_applied = True` and `mismatch_code = REF_ALLELE_MISMATCH`

Current stage-level QC status:

- all `6` stage-level QC checks pass
- no mitochondrial-encoded `ClinVar` row is eligible for the nuclear comparative branch
- no overlap-derived mtDNA row is lost during classification
- no unresolved row is counted as `uDAV`
- every classified row has a valid `classification_basis`
- every classified row has a resolved `alignment_source`
- every classified row has a resolved `classification_coordinate_method`

Relevant outputs:

- `data/derived/classified/classified_all.parquet`
- `data/derived/classified/classified_clean.parquet`
- `data/derived/classified/classified_warning.parquet`
- `data/derived/classified/classification_qc_summary.json`
- `data/derived/classified/classification_qc_checks.tsv`

## 5. cDAV metrics by OXPHOS complex

| genome | complex_id | classified_rows | aa_cdav | nt_cdav | aa_cdav_pct | nt_cdav_pct |
| --- | --- | --- | --- | --- | --- | --- |
| mtDNA | CI | 224 | 119 | 107 | 53.12 | 47.77 |
| mtDNA | CV | 99 | 48 | 42 | 48.48 | 42.42 |
| nucDNA | CI | 2018 | 732 | 642 | 36.27 | 31.81 |
| nucDNA | CII | 3101 | 812 | 672 | 26.19 | 21.67 |
| nucDNA | CIII | 402 | 135 | 116 | 33.58 | 28.86 |
| nucDNA | CIV | 319 | 160 | 140 | 50.16 | 43.89 |
| nucDNA | CV | 501 | 142 | 117 | 28.34 | 23.35 |
| nucDNA | Unassigned | 145 | 65 | 58 | 44.83 | 40.0 |

A tab-delimited copy of the complex-level breakdown is available at [dav_metrics_by_complex.tsv](/Users/ad2347/Documents/OxPhos_DAV/data/derived/results/dav_metrics_by_complex.tsv).

## 5A. Current structural mapping state

The structural mapping stage now consumes the classified Parquet master table and a manifest-backed panel of primary and validation human OXPHOS structures.

Current structural mapping outputs:

- `results/structural/dar_structure_map.csv`
- `results/structural/dar_contacts_cbcb8A.csv`
- `results/structural/dar_mito_nuc_contacts.csv`
- `results/structural/structure_model_summary.csv`

Primary versus validation usage:

- primary structures provide the default structural mapping and default contact set
- validation structures are used to test whether a mapping or contact is preserved across an alternative model or state
- structural support is therefore interpretable as:
  - primary-only support
  - validation-only support
  - multi-model support

Current verified mapping summary from the latest local structural outputs
(`results/structural/`, refreshed on 2026-05-14):

- classified rows loaded: `33779`
- unique source-variant summaries in `structure_model_summary.csv`: `17167`
- mapping rows written: `37549`
- mapping-eligible variants: `6673`
- unique variants with at least one successful structural map: `6326`
- true mapping rate: `94.8%`
- mapped model rows by status category: `9771`
  - `mapped_direct`: `8321`
  - `mapped_with_extended_offset_rescue`: `1130`
  - `mapped_with_isoform_offset`: `320`
- unique variants with at least one successful structural map:
  - `6326` by `structure_model_summary.csv` support class
  - `6326` by mapped-row status categories in `dar_structure_map.csv`
- canonical contact rows written: `64968`
- validation contact rows written: `29231`
- union contact rows written: `94199`
- canonical mito-nuclear contact rows written: `1208`

Gene-symbol boundary:

- `interpreted_gene` remains the HGNC/classification analytical key
- `structure_gene`, `dar_structure_gene`, and `contact_structure_gene` record
  PDB/UniProt resource aliases
- COXFA4 rows are therefore retained as `COXFA4` analytically and as `NDUFA4`
  structurally

Current per-variant structural support distribution:

- `10841` variant summaries have `0` mapped models
- `3454` variant summaries have primary-model-only support
- `2872` variant summaries have both primary and validation support
- `0` variant summaries have validation-only support

### 5A.1 Structural coverage caveat

The structural stage is now explicitly multi-model, but structural support is
still uneven across the classified variant set.

Why:

- a variant only counts as mapped in a given model if the script can identify
  the chain, align the human reference sequence, anchor the residue position,
  and verify the expected human reference amino acid in that structure
- this favors the subset of complexes, subunits, and residue positions that are
  well resolved and sequence-tractable in the available human models

Interpretation:

- successful mappings remain concentrated in the models that currently support
  robust residue anchoring
- many nucDNA rows remain structurally unmapped because current structures do
  not yet support confident residue-level mapping for those subunits or
  positions
- absence of structural support should therefore not be interpreted as absence
  of a structural mechanism

Practical consequence:

- structural enrichment results should always be interpreted relative to the
  structurally mappable subset
- cDAV versus uDAV structural comparisons require explicit structure-coverage
  denominators by complex and gene

### 5A.2 Structural coverage denominators

Structure-coverage summaries are now exported to:

- `data/derived/results/structural_mapping_coverage_by_complex.tsv`
- `data/derived/results/structural_mapping_coverage_by_gene.tsv`

Current per-variant support-class summary:

- `primary_and_validation`: `2872`
- `primary_only`: `3454`
- `validation_only`: `0`
- `unmapped`: `10841`

Current structural audit outputs:

- `results/structural/structure_mapping_failure_audit.csv`

Current structural row-status distribution from
`results/structural/dar_structure_map.csv`:

- `structural_ineligible`: `27106`
- `mapped_direct`: `8321`
- `mapped_with_extended_offset_rescue`: `1130`
- `mapped_with_isoform_offset`: `320`
- `secondary_not_attempted`: `306`
- `mature_offset_candidate`: `198`
- `residue_anchoring_failure`: `151`
- `chain_assignment_or_model_gap`: `17`

Interpretation of these categories:

- `structural_ineligible` rows are not structural mapping failures; they are
  rows outside the active structural branch, mostly skipped/non-classified
  source interpretations retained under the filter-late contract
- `secondary_not_attempted` rows are deliberate non-attempts after a variant
  has already been handled through the primary structural path
- `mapped_direct`, `mapped_with_extended_offset_rescue`, and
  `mapped_with_isoform_offset` are successful mappings
- `mature_offset_candidate` rows are plausible large-offset mappings that are
  deliberately kept diagnostic unless covered by the conservative anchor
  exception registry
- `residue_anchoring_failure` rows are current hard mapping failures where the
  expected residue cannot be confidently anchored in the selected model
- `chain_assignment_or_model_gap` rows are retained as explicit structural
  coverage gaps, not silently dropped

Current residual failure cases:

- `mature_offset_candidate`: `198` model rows
- `residue_anchoring_failure`: `151` model rows
- `chain_assignment_or_model_gap`: `17` model rows

The compact failure audit table currently contains `179` grouped rows:

- `mature_offset_candidate`: `150` audit groups
- `residue_anchoring_failure`: `28` audit groups
- `chain_assignment_or_model_gap`: `1` audit group

Current anchor-policy state after the conservative registry-gated rerun:

- extended-offset rescue is enabled only for `SDHA`, `SDHB`, `SDHC`, and `SDHD`
- non-`SDH` large-offset candidates are still recorded, but remain unresolved by
  policy
- this keeps the structural mapped set conservative while preserving the audit
  evidence for future manual review
- current registry-approved rescued-offset support contributes `1130` model rows and `751`
  variant summaries with at least one extended-offset structural mapping

Current complex-level structural coverage from the refreshed denominator export:

| complex_id | genome | variants_total | variants_mapped | mapping_pct |
| --- | --- | --- | --- | --- |
| CI | mtDNA | 224 | 224 | 100.0 |
| CV | mtDNA | 99 | 93 | 93.94 |
| CI | nucDNA | 2018 | 1670 | 82.76 |
| CII | nucDNA | 3101 | 2627 | 84.71 |
| CIII | nucDNA | 402 | 356 | 88.56 |
| CIV | nucDNA | 319 | 210 | 65.83 |
| CV | nucDNA | 501 | 389 | 77.64 |
| Unassigned | nucDNA | 145 | 0 | 0.0 |

Important interpretation:

- the earlier mtDNA-dominant mapping state was a methodological artifact, not a
  stable biological result
- applying transcript-position remapping and fixing partial chain-assignment
  fallback recovered substantial nucDNA structural coverage
- expanding the corrected structure panel with `5Z62` and `8H9U` increased
  mapped model rows and contact coverage further, especially for `CIV` and `CV`
- the remaining structural gaps are now concentrated in specific failure
  classes rather than hidden inside a generic unmapped majority
- downstream structural enrichment analyses should use the current
  `structure_model_summary.csv` and `structure_mapping_failure_audit.csv`
  outputs rather than the older denominator snapshot alone
- the current structural mapping is complete enough for denominator-aware
  structural analysis, but not for treating every classified variant as equally
  structurally testable

## 5B. Current compensating-partner state

The compensating-partner stage was rerun after the COXFA4/NDUFA4 structural
alias-boundary fix.

Current outputs:

- `results/structural/all_tested_pairs.csv`: `9068` rows
- `results/structural/compensatory_partners.csv`: `199` rows
- `results/structural/concordance_summary.csv`: `1447` rows

COXFA4/NDUFA4 handling:

- `dar_gene` and `contact_gene` use the HGNC/classification key `COXFA4`
- `dar_structure_gene` and `contact_structure_gene` preserve the PDB/UniProt
  resource alias `NDUFA4`
- COXFA4 cDAVs now contribute `19` tested DAR-contact rows
- no duplicate primary `(variant_id, dar_gene, dar_aa_coord, contact_gene,
  contact_refseq_pos, contact_alt_aa)` tuples were detected

## 6. Supplementary information

### 6.1 Important fixes and nuanced issues resolved so far

- The pipeline was converted to a metadata-first, filter-late framework: records are retained with explicit eligibility and exclusion metadata rather than being dropped during parsing or classification.
- Overlapping mtDNA loci, especially MT-ATP6/MT-ATP8, are duplicated into frame-specific interpretations, so one nucleotide event can yield two independent curated and classified rows.
- Gene identity was standardized around interpreted_gene, replacing older logic that collapsed composite loci such as locus.split("/")[0].
- The nuclear transcript-position map builder was corrected to parse full multi-line TOGA human protein FASTA sequences. This removed a major source of false POSITION_NOT_IN_ENST calls.
- Classification now reads the canonical curated Parquet master table, retains all rows, and writes canonical classified outputs under data/derived/classified/.
- Parquet-to-JSON export was stabilized by normalizing NumPy arrays and NaN values before writing compatibility JSON and JSONL outputs.
- Coordinate rescue is now ordered more defensibly: deterministic transcript maps are preferred when map quality is sufficient, genomic rescue is used for lower-concordance genes, and anchor fallback remains available only as a guarded last-resort heuristic.
- Residual problematic genes and rows are now handled via an explicit exception registry plus exported audit tables instead of hidden one-off code paths.

### 6.2 Important mitigations currently in use

- Anchor fallback: if no transcript map exists, the classifier can search within +/-10 amino-acid positions for the expected wild-type residue in the aligned human sequence before projecting codon coordinates.
- Genomic rescue: for lower-concordance genes, genomic coordinates are projected directly into the CDS space of the selected TOGA ENST model.
- Mismatch categories are explicit and not conflated with uDAV calls: REF_ALLELE_MISMATCH, TRANSCRIPT_MISMATCH, POSITION_NOT_IN_ENST, GENOMIC_POS_NOT_IN_ENST, ANCHOR_NOT_FOUND, CODON_EXTRACTION_FAILURE, NO_ALIGNMENT, and COORD_PARSE_FAILURE.
- Classified-with-warning rows are preserved but flagged when a cDAV/uDAV classification succeeds despite a reference-base disagreement.

### 6.3 DAV metrics by gene

The full gene-level breakdown is also exported as [dav_metrics_by_gene.tsv](/Users/ad2347/Documents/OxPhos_DAV/data/derived/results/dav_metrics_by_gene.tsv). The current gene-level supplementary tables are reproduced below.

#### mtDNA genes

| gene | complex_id | classified_rows | aa_cdav | nt_cdav | aa_cdav_pct | nt_cdav_pct |
| --- | --- | --- | --- | --- | --- | --- |
| MT-ND1 | CI | 74 | 36 | 31 | 48.65 | 41.89 |
| MT-ND2 | CI | 21 | 16 | 16 | 76.19 | 76.19 |
| MT-ND3 | CI | 12 | 6 | 5 | 50.0 | 41.67 |
| MT-ND4 | CI | 23 | 15 | 12 | 65.22 | 52.17 |
| MT-ND4L | CI | 6 | 2 | 2 | 33.33 | 33.33 |
| MT-ND5 | CI | 58 | 28 | 27 | 48.28 | 46.55 |
| MT-ND6 | CI | 30 | 16 | 14 | 53.33 | 46.67 |
| MT-ATP6 | CV | 83 | 36 | 30 | 43.37 | 36.14 |
| MT-ATP8 | CV | 16 | 12 | 12 | 75.0 | 75.0 |

#### nucDNA genes

| gene | complex_id | classified_rows | aa_cdav | nt_cdav | aa_cdav_pct | nt_cdav_pct |
| --- | --- | --- | --- | --- | --- | --- |
| NDUFA1 | CI | 27 | 16 | 12 | 59.26 | 44.44 |
| NDUFA10 | CI | 103 | 48 | 40 | 46.6 | 38.83 |
| NDUFA11 | CI | 43 | 21 | 17 | 48.84 | 39.53 |
| NDUFA12 | CI | 42 | 20 | 16 | 47.62 | 38.1 |
| NDUFA13 | CI | 44 | 5 | 3 | 11.36 | 6.82 |
| NDUFA2 | CI | 34 | 16 | 13 | 47.06 | 38.24 |
| NDUFA3 | CI | 14 | 11 | 9 | 78.57 | 64.29 |
| NDUFA5 | CI | 12 | 7 | 7 | 58.33 | 58.33 |
| NDUFA6 | CI | 40 | 20 | 16 | 50.0 | 40.0 |
| NDUFA7 | CI | 22 | 5 | 5 | 22.73 | 22.73 |
| NDUFA8 | CI | 33 | 6 | 6 | 18.18 | 18.18 |
| NDUFA9 | CI | 146 | 72 | 70 | 49.32 | 47.95 |
| NDUFAB1 | CI | 33 | 21 | 17 | 63.64 | 51.52 |
| NDUFB1 | CI | 7 | 4 | 4 | 57.14 | 57.14 |
| NDUFB10 | CI | 75 | 36 | 32 | 48.0 | 42.67 |
| NDUFB11 | CI | 46 | 20 | 18 | 43.48 | 39.13 |
| NDUFB2 | CI | 19 | 12 | 11 | 63.16 | 57.89 |
| NDUFB3 | CI | 26 | 14 | 12 | 53.85 | 46.15 |
| NDUFB4 | CI | 12 | 8 | 8 | 66.67 | 66.67 |
| NDUFB5 | CI | 38 | 21 | 17 | 55.26 | 44.74 |
| NDUFB6 | CI | 16 | 7 | 7 | 43.75 | 43.75 |
| NDUFB7 | CI | 26 | 10 | 9 | 38.46 | 34.62 |
| NDUFB8 | CI | 56 | 21 | 21 | 37.5 | 37.5 |
| NDUFB9 | CI | 59 | 16 | 15 | 27.12 | 25.42 |
| NDUFC1 | CI | 15 | 8 | 7 | 53.33 | 46.67 |
| NDUFC2 | CI | 3 | 1 | 1 | 33.33 | 33.33 |
| NDUFS1 | CI | 207 | 45 | 41 | 21.74 | 19.81 |
| NDUFS2 | CI | 114 | 23 | 17 | 20.18 | 14.91 |
| NDUFS3 | CI | 87 | 32 | 27 | 36.78 | 31.03 |
| NDUFS4 | CI | 57 | 16 | 14 | 28.07 | 24.56 |
| NDUFS5 | CI | 15 | 9 | 9 | 60.0 | 60.0 |
| NDUFS6 | CI | 46 | 18 | 16 | 39.13 | 34.78 |
| NDUFS7 | CI | 67 | 24 | 21 | 35.82 | 31.34 |
| NDUFS8 | CI | 84 | 16 | 13 | 19.05 | 15.48 |
| NDUFV1 | CI | 201 | 22 | 19 | 10.95 | 9.45 |
| NDUFV2 | CI | 60 | 14 | 13 | 23.33 | 21.67 |
| NDUFV3 | CI | 89 | 67 | 59 | 75.28 | 66.29 |
| SDHA | CII | 1569 | 307 | 251 | 19.57 | 16.0 |
| SDHB | CII | 743 | 150 | 126 | 20.19 | 16.96 |
| SDHC | CII | 383 | 166 | 140 | 43.34 | 36.55 |
| SDHD | CII | 406 | 189 | 155 | 46.55 | 38.18 |
| CYC1 | CIII | 99 | 25 | 23 | 25.25 | 23.23 |
| UQCR10 | CIII | 7 | 1 | 1 | 14.29 | 14.29 |
| UQCR11 | CIII | 11 | 7 | 7 | 63.64 | 63.64 |
| UQCRB | CIII | 27 | 11 | 7 | 40.74 | 25.93 |
| UQCRC1 | CIII | 94 | 29 | 28 | 30.85 | 29.79 |
| UQCRC2 | CIII | 91 | 34 | 28 | 37.36 | 30.77 |
| UQCRFS1 | CIII | 33 | 16 | 13 | 48.48 | 39.39 |
| UQCRH | CIII | 15 | 6 | 5 | 40.0 | 33.33 |
| UQCRQ | CIII | 25 | 6 | 4 | 24.0 | 16.0 |
| COX4I1 | CIV | 35 | 16 | 15 | 45.71 | 42.86 |
| COX4I2 | CIV | 49 | 31 | 25 | 63.27 | 51.02 |
| COX5A | CIV | 26 | 16 | 15 | 61.54 | 57.69 |
| COX5B | CIV | 20 | 11 | 10 | 55.0 | 50.0 |
| COX6A1 | CIV | 53 | 23 | 21 | 43.4 | 39.62 |
| COX6A2 | CIV | 27 | 4 | 3 | 14.81 | 11.11 |
| COX6B1 | CIV | 25 | 15 | 13 | 60.0 | 52.0 |
| COX6C | CIV | 15 | 7 | 5 | 46.67 | 33.33 |
| COX7A1 | CIV | 16 | 4 | 3 | 25.0 | 18.75 |
| COX7A2 | CIV | 12 | 7 | 6 | 58.33 | 50.0 |
| COX7B | CIV | 20 | 13 | 12 | 65.0 | 60.0 |
| COX7C | CIV | 6 | 4 | 4 | 66.67 | 66.67 |
| COX8A | CIV | 15 | 9 | 8 | 60.0 | 53.33 |
| ATP5F1A | CV | 140 | 18 | 16 | 12.86 | 11.43 |
| ATP5F1B | CV | 50 | 13 | 9 | 26.0 | 18.0 |
| ATP5F1C | CV | 39 | 14 | 10 | 35.9 | 25.64 |
| ATP5F1D | CV | 85 | 23 | 18 | 27.06 | 21.18 |
| ATP5F1E | CV | 12 | 1 | 1 | 8.33 | 8.33 |
| ATP5MC1 | CV | 19 | 15 | 15 | 78.95 | 78.95 |
| ATP5MC2 | CV | 18 | 11 | 9 | 61.11 | 50.0 |
| ATP5MC3 | CV | 23 | 8 | 8 | 34.78 | 34.78 |
| ATP5PB | CV | 46 | 10 | 7 | 21.74 | 15.22 |
| ATP5PD | CV | 23 | 9 | 8 | 39.13 | 34.78 |
| ATP5PF | CV | 13 | 4 | 3 | 30.77 | 23.08 |
| ATP5PO | CV | 33 | 16 | 13 | 48.48 | 39.39 |
| ATP5IF1 | Unassigned | 14 | 4 | 3 | 28.57 | 21.43 |
| ATP5ME | Unassigned | 18 | 3 | 3 | 16.67 | 16.67 |
| ATP5MF | Unassigned | 11 | 5 | 5 | 45.45 | 45.45 |
| ATP5MG | Unassigned | 18 | 7 | 5 | 38.89 | 27.78 |
| ATP5MJ | Unassigned | 3 | 2 | 1 | 66.67 | 33.33 |
| ATP5MK | Unassigned | 9 | 9 | 7 | 100.0 | 77.78 |
| COX6B2 | Unassigned | 12 | 6 | 6 | 50.0 | 50.0 |
| COX7B2 | Unassigned | 12 | 11 | 10 | 91.67 | 83.33 |
| COX8C | Unassigned | 15 | 11 | 11 | 73.33 | 73.33 |
| COXFA4 | Unassigned | 17 | 3 | 3 | 17.65 | 17.65 |
| COXFA4L2 | Unassigned | 14 | 2 | 2 | 14.29 | 14.29 |
| COXFA4L3 | Unassigned | 2 | 2 | 2 | 100.0 | 100.0 |

### 6.4 Current residual classification caveats

- Residual unresolved nuclear rows remain concentrated in a small set of genes with transcript-model incompatibility or gene-specific transcript/consequence discordance, especially `NDUFS6`, `NDUFA13`, `NDUFA11`, `UQCRB`, `NDUFS7`, `NDUFV2`, `NDUFB1`, and `NDUFA10`.
- These are now tracked explicitly via the exception registry rather than being hidden in implicit script logic.
- Focused review tables are available at `data/derived/classified/exception_candidate_rows.tsv` and `data/derived/classified/exception_candidate_summary.tsv`.

### 6.5 Current schema and database-handling nuances

- `classification_basis` now distinguishes `nt_and_aa`, `nt_only`, `aa_only`, and `no_disease_allele_detected`.
- This correction was necessary because nucleotide-positive rows were previously at risk of being labeled `nt_and_aa` even when amino-acid support was absent.
- The pipeline remains missense-focused in eligibility: only missense curated rows enter the core comparative branch.
- However, `nt_only` remains a valid schema state in principle because a human missense disease nucleotide can recur at the homologous nucleotide site in another species without reproducing the same amino-acid state in that species.
- In the current audited dataset there are no `nt_only` classified rows, but the schema now records the distinction correctly.

- Mitochondrial `ClinVar` rows are retained for traceability but are not used as eligible rows in the nuclear comparative branch.
- These rows are explicitly marked ineligible with `core_pipeline_exclusion_reason = mt_in_clinvar_nuclear_branch`.
- Stage-level QC now verifies that no mitochondrial-encoded `ClinVar` row becomes eligible for the nuclear branch.

- Cross-source mitochondrial overlap between `ClinVar` and `MITOMAP` is now recorded explicitly.
- Current audited overlap counts:
  - mitochondrial `ClinVar` rows retained: `4582`
  - unique mitochondrial `ClinVar` variant IDs: `2291`
  - unique mitochondrial `ClinVar` variant IDs also represented in `MITOMAP`: `234`
  - unique mitochondrial `ClinVar` variant IDs not represented in `MITOMAP`: `2057`
- Record-level inventories are written to:
  - `data/derived/results/clinvar_mt_variants_present_in_mitomap.tsv`
  - `data/derived/results/clinvar_mt_variants_absent_from_mitomap.tsv`
