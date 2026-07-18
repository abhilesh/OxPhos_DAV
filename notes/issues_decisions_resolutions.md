# Issues, Decisions, and Resolutions

This document is a consolidated ledger of problems discovered, decisions made, and resolutions applied across all pipeline stages completed to date. It is intended as a durable memory for downstream work — particularly for stages not yet run (compensating partners, phylogenetics, mutagenesis) where these decisions propagate.

Each entry follows the format: **Issue → Decision → Resolution → Status**.

**Coverage:** Stages 0–3 (data download, curation, classification, structural mapping) through 2026-05-13.

---

## Stage 0 — Data download

### D0-1: MITOMAP HTTP 403 on programmatic access
**Issue:** MITOMAP's disease variant endpoint returns HTTP 403 when accessed programmatically, blocking automated download.  
**Decision:** Download manually from https://www.mitomap.org/cgi-bin/disease.cgi and commit the TSV to `data/raw/annotations/`.  
**Resolution:** File at `data/raw/annotations/MITOMAP_CodingVariants_2026-03-26.tsv`. Download validator marks status as `ok;remote_http_403`.  
**Status:** Accepted limitation. Re-download on data refresh requires manual step.

### D0-2: MyVariant / dbNSFP / gnomAD completeness not guaranteed
**Issue:** MyVariant fetch used gene-level batch queries with `size=1000`; no check that all variants for a gene were returned (gene could have >1000 ClinVar entries).  
**Decision:** Do not rely on MyVariant-derived gnomAD/dbNSFP fields for primary manuscript claims until identity-based re-acquisition is implemented.  
**Resolution:** Current file marked `ok;completeness_legacy_unchecked`. Programmatic re-acquisition deferred.  
**Status:** Open. Downstream stages should not use gnomAD/dbNSFP annotation as a primary filter without refreshing this dataset via per-variant ID queries.

### D0-3: Cross-genome species matching redesigned from name heuristics to TaxID
**Issue:** Matching species between TOGA (nucDNA) and mtDNA alignment sources via name strings was brittle — the same species could appear under different name formats.  
**Decision:** Build species overlap from NCBI TaxID extracted from FASTA headers.  
**Resolution:** `src/data_download/00e_build_taxid_species_map.py` extracts TaxIDs; `src/data_curation/04_sanitize_all_alignments.py` normalises FASTA headers. Cross-genome overlap is TaxID-keyed.  
**Status:** Resolved. TaxID is the authoritative cross-genome species identifier throughout the pipeline.

---

## Stage 1 — Data curation

### D1-1: Early-drop records violated filter-late contract
**Issue:** Older curation parsers dropped non-core records (noncoding, frameshift, nonsense, non-SNV, unparseable AA) before they ever became rows in the master table. Absence of a row was overloaded to mean "excluded", "unparsable", or "never seen".  
**Decision:** Adopt filter-late/Scope-A contract: the pipeline targets coding SNVs with interpretable AA changes; records outside that scope are retained with `eligible_core_comparative_pipeline = False` and an explicit exclusion reason rather than being silently dropped.  
**Resolution:** Curation layer refactored. Master table is `data/derived/curated/variants_master_curated.parquet`. Eligibility expressed via `eligible_core_comparative_pipeline`, `core_pipeline_exclusion_reason`, `parse_status`, and `curation_status` fields.  
**Status:** Resolved.

### D1-2: MT-ATP6/MT-ATP8 overlap loci collapsed to first gene
**Issue:** Composite mitochondrial loci (`MT-ATP6/MT-ATP8`) were reduced to the first token by `.split("/")[0]`, discarding the second gene's frame-specific interpretation. One nucleotide event can have two distinct coding consequences.  
**Decision:** Implement overlap-aware duplication: each overlap nucleotide event produces two independent curated rows, one per interpreted gene.  
**Resolution:** `src/data_curation/01_curate_variants.py` detects overlap loci and emits duplicate rows with `interpreted_gene`, `is_overlap`, `overlap_group`, `overlap_genes`, `overlap_role`, `derived_from_overlap_duplication`, `shared_source_variant_group_id` fields. Classification QC verifies no overlap row is lost.  
**Status:** Resolved for curation, classification, and structural mapping. Caution: legacy phylogenetic scripts (`src/phylo/`) may still contain `.split("/")[0]` patterns — audit before running overlap-aware phylogenetic analyses.

### D1-3: Gene identity too loosely defined (`locus` vs `interpreted_gene`)
**Issue:** Older code paths used `locus` as the authoritative gene name. For overlap loci, `locus` is a composite string (`MT-ATP6/MT-ATP8`) rather than the gene whose frame the interpretation belongs to.  
**Decision:** `interpreted_gene` is the authoritative single-gene identity field throughout the pipeline.  
**Resolution:** Classification, structural mapping, and compensating-partners scripts all read `interpreted_gene` (falling back to `classification_gene`, then `locus`) for gene lookups.  
**Status:** Resolved.

### D1-4: MT-ND6 minus-strand handling
**Issue:** MT-ND6 is encoded on the mitochondrial minus strand; raw MITOMAP alleles are genomic-strand alleles and cannot be used directly as coding-strand alleles for codon construction.  
**Decision:** Complement `ref` and `alt` for positions in the MT-ND6 interval at parse time.  
**Resolution:** MITOMAP parser applies reverse-complement for MT-ND6 coding-strand alleles; original genomic alleles stored separately as `genomic_ref`/`genomic_alt`.  
**Status:** Resolved.

### D1-5: CDS and genomic alleles kept separate
**Issue:** Downstream steps (gnomAD lookup, mtDNA resource cross-referencing) need to distinguish strand-corrected CDS alleles from raw genomic alleles.  
**Decision:** Store both representations: `ref`/`alt` (coding strand) and `genomic_ref`/`genomic_alt` (genomic strand).  
**Resolution:** Implemented in parsers and propagated through `VariantRecord`. Both fields are present in the curated master table.  
**Status:** Resolved.

### D1-6: Cross-source mtDNA overlap between MITOMAP and ClinVar
**Issue:** The same mtDNA variant can appear in both MITOMAP and ClinVar; without explicit linkage, they appear as independent records downstream.  
**Decision:** Retain all records but link cross-source duplicates; exclude mt-ClinVar rows from the nuclear comparative branch (they are not eligible for TOGA-based classification).  
**Resolution:** Cross-source linkage written during curation. mt-ClinVar rows retained with `core_pipeline_exclusion_reason = mt_in_clinvar_nuclear_branch`. Current overlap counts documented in `notes/running_results.md` §6.5.  
**Status:** Resolved.

---

## Stage 2 — Classification

### D2-1: Coordinate rescue cascade ordering
**Issue:** Without a deterministic priority for coordinate sources, classification could silently use the weakest source (anchor fallback) for genes where a better source exists.  
**Decision:** Prefer transcript-position map (NM_→ENST) → genomic-coordinate map → anchor fallback (last resort only).  
**Resolution:** `src/classify/00_classify_DAV.py` applies this cascade; `coordinate_resolution_method` field records which source was used (`tx_map | genomic_map | anchor`). Anchor fallback applies only when no map coverage exists.  
**Status:** Resolved.

### D2-2: Transcript-position map parser truncated human protein sequences
**Issue:** The transcript-position map builder was reading only the first line of multi-line FASTA sequences, producing truncated human reference proteins and causing false `POSITION_NOT_IN_ENST` failures for positions in the second or later FASTA lines.  
**Decision:** Fix the FASTA parser to accumulate all continuation lines.  
**Resolution:** `src/data_curation/02_build_transcript_position_maps.py` corrected to parse full multi-line FASTA. A large class of false POSITION_NOT_IN_ENST calls was eliminated.  
**Status:** Resolved.

### D2-3: Anchor fallback ambiguity
**Issue:** The ±10 aa anchor search can find multiple occurrences of the same WT amino acid within the window, causing the anchor to resolve to the wrong position.  
**Decision:** Record ambiguous anchor cases as unresolved rather than forcing a mapping; flag for future short-peptide-context anchor improvement.  
**Resolution:** Ambiguous anchors produce `mismatch_code = ANCHOR_NOT_FOUND` or `classification_status = unresolved` rather than a forced cDAV/uDAV call. Current unresolved nucDNA rows concentrated in NDUFS6, NDUFA13, NDUFA11, UQCRB, NDUFS7, NDUFV2, NDUFB1, NDUFA10.  
**Status:** Partially resolved. The safe-fallback behaviour is implemented. The underlying ambiguity in some genes remains; tracked via exception registry candidate tables.

### D2-4: Exception registry for hard cases
**Issue:** A small number of genes and variants require manual offset or classification overrides that cannot be resolved algorithmically (transcript-version drift, unusual isoform choice, validated transit-peptide offsets).  
**Decision:** Maintain an explicit exception registry rather than encoding one-off logic in scripts.  
**Resolution:** `data/reference/variant_exception_registry.tsv` for per-variant/per-gene classification overrides. `data/reference/structural_anchor_exception_registry.tsv` for transit-peptide offset rescue, gated on UniProt evidence and TOGA alignment validation. Exception-applied rows tagged `classification_warning = True`.  
**Status:** Resolved. Registry is the authoritative place for all manual overrides.

### D2-5: Skipped-by-policy rows must be retained, not dropped
**Issue:** Benign and VUS ClinVar variants (and synonymous variants) were previously skipped entirely by the classifier, making them invisible in all downstream outputs. This caused the structural mapping "ineligible" count to be misinterpreted as a structural failure.  
**Decision:** Retain all rows from the curated master table; benign/VUS/synonymous variants receive `classification_status = "skipped_by_policy"` and `structure_mapping_eligible = False`. They are never attempted for structural mapping.  
**Resolution:** Implemented in `src/classify/00_classify_DAV.py`. The 26,902 skipped_by_policy rows are visible in the master classified parquet and accounted for explicitly in the structural mapping summary.  
**Status:** Resolved. Critical for correct structural mapping-rate reporting (see D3-3).

### D2-6: `classification_basis` distinction
**Issue:** `nt_and_aa` was being assigned even when AA-level support was absent, overcounting AA-level evidence.  
**Decision:** Use four explicit values: `nt_and_aa`, `nt_only`, `aa_only`, `no_disease_allele_detected`.  
**Resolution:** Corrected in `src/classify/00_classify_DAV.py`. Current audited dataset has no `nt_only` rows (correct; would require same codon but different AA, which is rare), but the schema correctly supports it.  
**Status:** Resolved.

### D2-7: Non-biological sequence states excluded from species counts
**Issue:** Gaps (`-`), frameshift markers (`!`), stops (`*`), and unknown residues (`X`) in TOGA alignments should not count as species support for or against compensation.  
**Decision:** Exclude these states from the species-level cDAV count in `AlignmentParser.check_compensation()`.  
**Resolution:** Implemented in `src/utils/alignment_parser.py`. All cDAV and uDAV species counts are based on biologically valid residue states only.  
**Status:** Resolved.

---

## Stage 3 — Structural mapping

### D3-1: Structure manifest audit — CV priority corrected (2026-05-13)
**Issue:** ATP synthase (CV) validation structures were listed with an incorrect priority order: 8H9T (2.77 Å, state 2) was at priority 2 and 8H9U (2.61 Å, state 3a) was at priority 3. Higher resolution should win in contact deduplication, so 8H9U should be priority 2.  
**Decision:** Swap priorities: 8H9U → priority 2 (higher resolution), 8H9T → priority 3.  
**Resolution:** `data/reference/structure_model_manifest.tsv` updated. Notes columns updated to document the rationale.  
**Status:** Resolved.

### D3-2: CIV primary structure — mature state over higher resolution (2026-05-13)
**Issue:** 9I6F (HIGD2A-bound assembly intermediate, 2.95 Å) has higher resolution than 9I7U (NDUFA4-bound mature state, 3.15 Å). The question was whether resolution or physiological state should determine the primary.  
**Decision:** Retain 9I7U as primary. The mature NDUFA4-bound state is physiologically relevant for interpreting pathogenic variants in functioning Complex IV; the assembly intermediate represents a transient state not representative of the complex in which pathogenic variants exert their effect.  
**Resolution:** `data/reference/structure_model_manifest.tsv` unchanged for CIV. Notes column documents the rationale.  
**Status:** Resolved and documented.

### D3-3: COXFA4/NDUFA4 resource-alias boundary — 17 variants silently dropped (2026-05-14)
**Issue:** The frozen HGNC table uses approved symbol `COXFA4` and previous symbol `NDUFA4`, while PDB/UniProt structural resources still commonly annotate the same CIV subunit as `NDUFA4`. The structural mapper originally needed the structural alias to find CIV chains, but overwriting `interpreted_gene` with `NDUFA4` made Stage 3 outputs inconsistent with the frozen Stage 2 analytical key.  
**Decision:** Keep `interpreted_gene = COXFA4` as the curated/classified analytical key and add explicit structural provenance fields for resource aliases.  
**Resolution:** `src/structural/00_map_davs_to_structure.py` now preserves `interpreted_gene = COXFA4` and writes `structure_gene = NDUFA4`; contact outputs write `dar_locus = COXFA4`, `dar_structure_gene = NDUFA4`, `contact_gene = COXFA4`, and `contact_structure_gene = NDUFA4` where applicable. After rerun: 17 COXFA4 classified variants map structurally; `dar_contacts_cbcb8A.csv` contains 140 COXFA4 DAR contact rows with NDUFA4 structural provenance, and `dar_mito_nuc_contacts.csv` contains 19 COXFA4 DAR mito-nuclear contact rows.  
**Status:** Resolved.

### D3-4: Mapping-rate misreporting — skipped_by_policy conflated with failures (2026-05-13)
**Issue:** The structural mapping summary printed "Classified rows loaded: N" vs "Mapped model rows: M", implying a mapping rate of M/N. With 26,902 skipped_by_policy rows, the naive ratio produced a misleadingly low apparent mapping rate.  
**Decision:** Report separate counts for skipped_by_policy (never attempted), mapping-eligible, mapped unique variants, eligible-but-failed, and the true mapping rate (mapped/eligible).  
**Resolution:** Summary block in `src/structural/00_map_davs_to_structure.py` replaced. True mapping rate: 6326/6673 = **94.8%**.  
**Status:** Resolved.

### D3-5: Residue anchoring design
**Issue:** PDB chain naming and sequence content are often inconsistent with gene naming and reference isoforms. Simple chain-label or exact-sequence lookups fail silently.  
**Decision:** Use a two-step approach: (1) RCSB API for chain-to-gene assignment, with local sequence-similarity fallback; (2) global alignment for position mapping; (3) final AA identity check against the expected reference amino acid.  
**Resolution:** Implemented in `src/structural/00_map_davs_to_structure.py`. All three steps are required before a variant is declared `mapped_direct`.  
**Status:** Resolved.

### D3-6: Multi-model contact deduplication and contact_source tracking
**Issue:** With multiple structures per complex (primary + validation), contacts from different models needed to be reconciled without double-counting.  
**Decision:** Primary model contacts → `dar_contacts_canonical.csv`; validation model contacts → `dar_contacts_validation.csv`; union with `contact_source` field → `dar_contacts_cbcb8A.csv`. Downstream analyses should use `contact_source` to distinguish primary-only, validation-only, and multi-model-supported contacts.  
**Resolution:** Implemented. Three separate CSVs written. `contact_source` column in union file is the handle for support-class stratification.  
**Status:** Resolved.

### D3-7: Isoform proxy mappings
**Issue:** Some genes in the variant set (COX4I2, COX6A2, COX7A1) are absent from the selected structures. They can be mapped through a closely related isoform (COX4I1, COX6A1, COX7A2 respectively), but this is a structural assumption.  
**Decision:** Allow proxy mappings but tag them explicitly; treat as a sensitivity class rather than primary evidence.  
**Resolution:** Proxy mappings tagged `mapped_with_isoform_offset` and `proxy_mapping_used = True` in `dar_structure_map.csv`. Downstream analyses should stratify direct vs proxy-based mappings.  
**Status:** Resolved at the labelling level. Current `dar_structure_map.csv` has zero `other` status-category rows.

### D3-8: No-panel-coverage genes (accepted limitation)
**Issue:** 12 genes in the canonical OXPHOS variant set have classified variants but are absent from all active cryo-EM chains (peripheral or tissue-specific subunits not resolved in current structures).  
**Decision:** Accept as a structural panel limitation. Do not add these genes to GENE_COMPLEX. Disclose in manuscript.  
**Genes:** ATP5ME (18 variants), ATP5MG (18), COXFA4L2 (14), COX8C (15), ATP5IF1 (14), COX6B2 (12), COX7B2 (12), ATP5MF (11), ATP5MK (9), ATP5MJ (3), COXFA4L3 (2). Note: COXFA4 (17 variants) is structurally represented through the NDUFA4 resource alias and is no longer a panel gap (resolved by D3-3).  
**Status:** Accepted limitation, documented in `notes/structural_mapping_eligibility_gaps.md`.

### D3-9: Duplicate-protomer handling (open)
**Issue:** For symmetric or repeated assemblies, the current mapper stores one chain per gene (the last chain assigned wins). Contact neighbourhoods may differ slightly by local resolution or assembly context across protomers.  
**Decision:** Defer full protomer-aware mapping; note as a caveat for protomer-specific contact claims.  
**Status:** Open. Affects manuscript interpretation if making claims about specific protomer contacts. Gene-level structural mapping is not affected.

### D3-10: Residual mapping failures
**Issue:** 347 mapping-eligible variants remain unmapped after the registry-gated pipeline.  
**Classes:**
- `mature_offset_candidate` (198 model rows; 150 grouped failure-audit rows): Large positive AA offsets consistent with unresolved transit-peptide cleavage. Conservative policy: only SDHA/SDHB/SDHC/SDHD and validated-by-UniProt genes are in the extended-offset rescue registry. Remaining candidates are preserved as diagnostic rows.
- `residue_anchoring_failure` (151 model rows; 28 grouped failure-audit rows): Position present in the alignment but absent from the PDB chain (unresolved loop or terminal region). Per-variant, not a panel gap.
**Decision:** Maintain conservative registry-gated policy. Do not rescue large offsets without UniProt transit-peptide evidence.  
**Status:** Open. Future work: AlphaFold supplement for residue-anchoring failures; registry expansion for additional validated transit-peptide genes.

---

## Open issues propagating to downstream stages

| ID | Issue | Affects | Priority |
|---|---|---|---|
| D0-2 | MyVariant/gnomAD completeness legacy-unchecked | Any analysis using population frequency filters | Fix before manuscript |
| D1-2 (residual) | Legacy `locus.split("/")[0]` patterns in phylo scripts | Phylogenetic timing of overlap-locus variants | Fix before running phylo |
| D2-3 (residual) | Anchor ambiguity in ~8 genes | ~68 unresolved nucDNA rows | Document as limitation |
| D3-9 | Duplicate-protomer handling | Protomer-specific contact claims | Document as limitation |
| D3-10 | 198 mature_offset_candidate + 151 anchoring-failure model rows | Structural coverage for transit-peptide genes | AlphaFold supplement or registry expansion |
