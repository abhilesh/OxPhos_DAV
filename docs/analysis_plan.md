# OxPhos DAV — Analysis Plan

**Project:** Compensated disease-associated variants (c-DARs) in OXPHOS  
**Last updated:** 2026-05-13

---

## 1. Scientific question

Identify positions in OXPHOS subunits where a human-pathogenic amino acid substitution is naturally tolerated in non-human mammals. The hypothesis is that such tolerance is mechanistically enabled by compensating changes at structurally neighbouring residues that co-evolved with the variant site.

Two levels of evidence for compensation:

- **AA-level cDAV**: ≥1 non-human species carries the human disease amino acid at the orthologous position. Discovery-level; one-species and multi-species cases should be stratified.
- **NT-level cDAV**: ≥1 non-human species uses the exact same codon as the human pathogenic allele. More stringent; requires the same nucleotide change.

The core question: are cDAV structural contact neighbours enriched for co-evolutionary signal (Pagel's discrete test, branch co-occurrence, mutual information) compared to uDAV neighbours?

---

## 2. Data sources

| Resource | Use | Local path |
|---|---|---|
| ClinVar variant_summary | nucDNA pathogenic/likely-pathogenic variants in OXPHOS genes | `data/raw/annotations/ClinVar_VariantSummary_2026-03-26.txt.gz` |
| MITOMAP coding variants | mtDNA disease-associated variants | `data/raw/annotations/MITOMAP_CodingVariants_2026-03-26.tsv` |
| TOGA hg38 alignments | Codon-aware and AA-level cross-species alignments for 241 mammals | `data/alignments/toga_hg38_aa/`, `data/alignments/toga_hg38_codon/` |
| mtDNA alignments | mtDNA codon-aware and AA-level species alignments | `data/alignments/mtdna_aa/`, `data/alignments/mtdna_codon/` |
| HGNC canonical OXPHOS gene list | Gene universe (100 genes: 13 mtDNA, 87 nucDNA) | `data/reference/Canonical_OXPHOS_Subunits_HGNC_2026-03-25.csv` |
| MANE GRCh38 v1.5 | Canonical transcript selection for nucDNA genes | `data/reference/MANE_GRCh38_v1.5_2026-04-18.txt.gz` |
| TOGA overview table | Ortholog quality scores per gene/species | `data/reference/TOGA_overview_table_hg38_2026-04-07.tsv` |
| PDB/CIF structures | Cryo-EM structures for OXPHOS CI–CV | `data/structures/` (manifest-backed panel) |
| PhyloTree build 17 | mtDNA haplogroup tree | `data/reference/PhyloTree_build_17_2026-03-26.zip` |
| VertLife mammal tree | Timetree for Pagel's discrete test | Used by `src/phylo/` |
| MitImpact / dbNSFP / gnomAD | Variant severity and population frequency annotations | `data/raw/annotations/` (gnomAD via MyVariant; completeness legacy-unchecked) |
| aaindex | Physicochemical amino acid property matrices | `data/reference/aaindex_properties_2026-04-07.json` |

---

## 3. Main pipeline

### Stage 0 — Data download
**Scripts:** `src/data_download/00a_download_mitomap.py` through `00j_download_structures.py`  
**Validation:** `00i_validate_downloads.py`

Downloads and validates all raw inputs. Outputs are immutable once downloaded. Key outputs:
- Raw annotation files in `data/raw/annotations/`
- TOGA FASTA alignments in `data/alignments/`
- PDB/CIF structure files in `data/structures/`

---

### Stage 1 — Data curation
**Scripts:** `src/data_curation/01_curate_variants.py`, `02_build_transcript_position_maps.py`, `03_build_genomic_coordinate_maps.py`, `04_sanitize_all_alignments.py`

Standardises variant records from MITOMAP and ClinVar into a canonical schema. Key design: **filter-late** — all in-scope records are retained with eligibility metadata rather than dropped during parsing.

Key outputs:
- `data/derived/curated/variants_master_curated.parquet` — canonical variant master table (33,779 rows: 388 mtDNA + 33,391 nucDNA)
- `data/reference/transcript_position_maps.json` — NM_→ENST AA position maps
- `data/reference/genomic_coordinate_maps.json` — genomic-position→CDS/AA maps

---

### Stage 2 — cDAV classification
**Scripts:** `src/classify/00_classify_DAV.py`, `src/align/00_align_translate_mtDNA.py`, `src/align/00_translate_nucDNA.py`

Reads curated parquet, runs `AlignmentParser.check_compensation()` per variant, and classifies as cDAV / uDAV / unresolved. Coordinate rescue cascade: transcript-position map → genomic-coordinate map → anchor fallback (±10 aa).

Key outputs:
- `data/derived/classified/variants_master_classified.parquet` — all rows with classification metadata
- `data/derived/classified/classified_clean.parquet` — no exceptions or ref-allele mismatches (6,340 rows)
- `data/derived/classified/classified_warning.parquet` — exception-applied or REF_ALLELE_MISMATCH rows (469 rows)
- `data/derived/classified/classification_qc_summary.json` + `classification_qc_checks.tsv`

Current counts (2026-05-13):

| Genome | Classified | AA-cDAV | NT-cDAV | AA-cDAV% |
|---|---|---|---|---|
| mtDNA | 323 | 167 | 149 | 51.7% |
| nucDNA | 6,486 | 2,046 | 1,745 | 31.5% |

---

### Stage 3 — Structural mapping
**Script:** `src/structural/00_map_davs_to_structure.py`

Maps classified variants onto primary and validation cryo-EM structures. Chain-to-gene assignment uses RCSB API then sequence fallback. Position mapping uses global alignment; final AA identity check required.

Structure panel (`data/reference/structure_model_manifest.tsv`):

| Complex | Primary | Validation | Notes |
|---|---|---|---|
| CI | 9I4I (2.63 Å) | 9TI4 (2.66 Å) | |
| CII | 8GS8 (2.86 Å) | — | |
| CIII | 9HZL (2.52 Å) | — | |
| CIV | 9I7U (3.15 Å, NDUFA4-bound) | 9I6F (2.95 Å, assembly intermediate) | Mature state kept as primary despite lower resolution |
| CV | 8H9S (2.53 Å, state 1) | 8H9U (2.61 Å, state 3a, p2); 8H9T (2.77 Å, state 2, p3) | Priority corrected 2026-05-13 |

Contact definition: Cβ–Cβ ≤ 8 Å (Cα for Gly); classified as hbond / electrostatic / hydrophobic / vdw.

Key outputs (`results/structural/`):
- `dar_structure_map.csv` — one row per variant × model (37,549 rows)
- `dar_contacts_canonical.csv` — contacts from primary models (64,968 rows)
- `dar_contacts_validation.csv` — contacts from validation models (29,231 rows)
- `dar_contacts_cbcb8A.csv` — union with `contact_source` column (94,199 rows)
- `dar_mito_nuc_contacts.csv` — mito-nuclear interface contacts (1,208 rows)
- `structure_model_summary.csv`, `structure_mapping_failure_audit.csv`

Current mapping summary (2026-05-13):

| Category | Count |
|---|---|
| Classified rows loaded | 33,779 |
| skipped_by_policy (benign/VUS, not attempted) | 26,902 |
| Mapping-eligible | 6,673 |
| Mapped (unique variants) | 6,326 |
| Eligible but failed | 347 |
| **True mapping rate** | **94.8%** |

Failure breakdown: 198 `mature_offset_candidate` model rows, 151 `residue_anchoring_failure` model rows, and 17 `chain_assignment_or_model_gap` model rows. `structure_mapping_failure_audit.csv` groups these into 179 audit rows.

Gene-symbol boundary: Stage 3 preserves `interpreted_gene` as the HGNC/classification key and uses `structure_gene` for PDB/UniProt aliases. For the renamed CIV subunit, outputs retain `COXFA4` analytically and record `NDUFA4` as structural provenance.

---

### Stage 4 — Compensating partners
**Script:** `src/structural/01_find_compensating_partners.py`

Tests whether structural contact neighbours of cDAV positions are enriched for co-evolutionary signatures. Input: `dar_contacts_cbcb8A.csv`. Uses three complementary tests:
- Fisher's exact test (species co-occurrence; retained for comparison, not primary evidence)
- Pagel's discrete test via R `phytools` (phylogenetically corrected co-evolution)
- Branch co-occurrence using IQTree ancestral states

Outputs:
- `all_tested_pairs.csv` — all tested DAR-contact pairs (9,068 rows)
- `compensatory_partners.csv` — significant derived view (199 rows)
- `concordance_summary.csv` — per-cDAV summary (1,447 rows)

**Status:** Rerun after the COXFA4/NDUFA4 boundary fix. COXFA4 cDAVs now enter the partner analysis without duplicate COXFA4/NDUFA4 rows.

---

### Stage 5 — Phylogenetic timing
**Scripts:** `src/phylo/` + HPC Slurm scripts in `scripts/`

Uses IQTree ancestral state reconstruction on the VertLife mammal tree to test whether compensating substitutions at partner residues preceded or followed the appearance of the disease allele at the variant site.

**Status:** HPC-bound, not yet run.

---

## 4. Supplementary analyses

### 4a — FoldX ΔΔG (structural stability)
**Scripts:** `src/mutagenesis/`  
**Purpose:** Estimate how much the pathogenic variant destabilises the subunit structure (ΔΔG_variant) and whether the compensating partner substitution partially rescues stability (ΔΔG_double_mutant vs ΔΔG_variant).  
**Method:** FoldX `RepairPDB` on the primary structure, then `BuildModel` for single and double mutants. BLOSUM62, Miyata distance, KD hydrophobicity, and volume delta from `src/utils/variant_record.py` complement the structural stability score.

### 4b — Physicochemical exchangeability
**Scripts:** `src/utils/variant_record.py`, `src/mutagenesis/`  
**Purpose:** Score the severity of each pathogenic substitution and the compensating substitution on independent physicochemical scales (BLOSUM62 log-odds, Miyata distance, KD hydrophobicity delta, volume delta).  
**Use:** Rank cDAV/partner pairs and test whether compensating partners tend to match or buffer the physicochemical impact of the disease allele.

### 4c — Shannon entropy and direct coupling analysis (DCA)
**Scripts:** `src/mutagenesis/`  
**Purpose:** Detect co-evolutionary signal between the variant site and its structural contact neighbours in the TOGA mammalian alignment.  
**Method:** Shannon entropy per alignment column (site conservation); pairwise mutual information; direct coupling analysis via pydca/evcouplings for residue-pair co-evolution scores.  
**Use:** Independent, sequence-only evidence for compensation that does not depend on structural contact geometry.

---

## 5. Key design principles

| Principle | Implementation |
|---|---|
| Filter-late | Variants retained through curation and classification with eligibility fields; filtering only at derived analysis views |
| Metadata-first | All failures, exclusions, and uncertainty recorded as column values; no silent row-drops after curation |
| Two-level evidence | AA-level (discovery) and NT-level (stringent) cDAV calls kept separate; NT-only is a valid schema state |
| Multi-model structural panel | Primary models for headline results; validation models for reproducibility support classes |
| Exception registry | Transit-peptide offset rescue gated on UniProt evidence; per-variant overrides explicit and auditable |
| Coordinate rescue cascade | Deterministic tx-map → genomic-map preferred; anchor fallback only as a last resort |
| Cross-genome via TaxID | mtDNA/nucDNA species overlap built from TaxID, not name strings |

---

## 6. Repo layout (key paths)

```
data/
  raw/annotations/          ← downloaded source variant files
  raw/reference/            ← TOGA, MANE, HGNC, PhyloTree, etc.
  alignments/               ← TOGA + mtDNA FASTA alignments
  structures/               ← PDB/CIF files
  derived/curated/          ← variants_master_curated.parquet
  derived/classified/       ← variants_master_classified.parquet + QC
  reference/                ← exception registries, manifest, position maps

src/
  data_download/            ← 00a–00j download scripts
  data_curation/            ← 01–04 curation scripts
  align/                    ← translation scripts
  classify/                 ← 00_classify_DAV.py, 01_audit_exception_candidates.py
  structural/               ← 00_map_davs_to_structure.py, 01_find_compensating_partners.py
  phylo/                    ← phylogenetic analysis scripts
  mutagenesis/              ← FoldX, DCA, exchangeability scripts
  utils/                    ← variant_record.py, alignment_parser.py, gene_reference.py, etc.

results/
  structural/               ← structure map, contacts, partner CSVs
  phylo/                    ← ancestral state outputs
  mutagenesis/              ← FoldX and DCA outputs

notes/                      ← authoritative pipeline documentation
docs/                       ← stable reference documents (this file)
tests/                      ← invariant tests against live parquet outputs
```
