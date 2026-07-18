# Ranked Unresolved Issues for the Current Pipeline

This document turns the currently unresolved issues into a practical action list.

## Current status update, 2026-05-12

This note is preserved as a historical issue ledger. Older problem statements
are not deleted; instead, each item is marked with its current status.

Current interpretation:

- the original curation/classification rerun blockers have mostly been resolved
  in the current `src/data_download`, `src/data_curation`, `src/classify`, and
  `src/structural` layout
- structural mapping now runs and produces substantial coverage, but it is not a
  fully neutral denominator across genes or complexes
- the main remaining manuscript-grade issues are structural eligibility,
  residue anchoring edge cases, proxy/isoform sensitivity, and downstream
  phylogenetic/partner-analysis assumptions

Current structural mapping status from the latest local outputs:

- `results/structural/dar_structure_map.csv`: `37532` mapping rows
- mapped model rows by status category:
  - `mapped_direct`: `8298`
  - `mapped_with_extended_offset_rescue`: `443`
  - `mapped_with_isoform_offset`: `311`
- non-mapped structural diagnostic rows:
  - `mature_offset_candidate`: `505`
  - `residue_anchoring_failure`: `154`
  - `other`: `6`
- `results/structural/dar_contacts_cbcb8A.csv`: `92894` contact rows
- `results/structural/dar_mito_nuc_contacts.csv`: `2101` mito-nuclear contact
  rows

Conclusion:

- structural mapping is complete enough for denominator-aware structural
  analysis
- it is not complete enough to treat all classified variants as equally
  structurally testable
- no glaring stage-breaking exception is visible in the latest outputs, but the
  structural caveats below remain active for manuscript interpretation

The categories are:
- `Must Fix Before Rerun`: likely to break the pipeline or systematically misclassify variants if you regenerate results
- `Should Fix Before Manuscript`: unlikely to block execution, but important for validity, interpretation, or reviewer scrutiny
- `Acceptable as Documented Limitations`: real caveats that can be carried if clearly disclosed and sensitivity-checked
- `Resolved / Historical`: older issue retained for provenance; no longer an
  active blocker in the current implementation

## Must Fix Before Rerun

No current curation/classification rerun blocker is known from this list. The
items that were originally in this section are retained below and marked with
their current status.

### 1. Overlapping mtDNA loci are still collapsed to the first gene

Current status:

- `Resolved / Historical` for curation, classification, and structural mapping.
- Current curation duplicates overlap loci into frame-specific interpretations
  with `interpreted_gene`, overlap metadata, and shared source-event fields.
- Current classification QC verifies that no overlap-derived mtDNA row is lost.
- Current structural mapping reads `interpreted_gene` rather than treating
  composite source locus as the authoritative gene.
- Residual caution: legacy/phylo scripts still contain `locus.split("/")[0]`
  patterns and should be reviewed before reusing those paths for overlap-aware
  analyses.

Problem:
Several scripts still reduce composite loci such as `MT-ATP6/MT-ATP8` to the first token using `.split("/")[0]`.

Where:
- `src/classify/00_classify_DAV.py`
- `src/utils/variant_record.py`
- `src/structural/00_map_davs_to_structure.py`
- likely anywhere else using `locus.split("/")[0]` as authoritative gene identity

Why this is critical:
- The wrong gene can be used for CDS mapping, complex assignment, gene context, and structural mapping.
- Overlap loci are biologically special cases, so forcing them into a single-gene interpretation can create systematic errors rather than random noise.

What to do:
- Introduce explicit handling for overlapping mtDNA loci.
- Store a primary locus and a secondary/overlap locus separately, or a structured `loci` field.
- Require downstream steps to resolve which reading frame/gene the amino-acid change belongs to instead of defaulting to the first token.

### 2. `00c_download_mtDNA_NCBI_seqs.py` expects `primary_symbol`, but `GeneReference` does not provide it

Current status:

- `Resolved / Historical`.
- Current `src/utils/gene_reference.py` provides `primary_symbol`.
- Current `src/data_download/00c_download_mtdna_ncbi_seqs.py` uses
  `gene_data.get("symbol") or gene_data.get("primary_symbol")`.
- The old `src/data_prep` path named below has been superseded by
  `src/data_download`.

Problem:
The mtDNA raw sequence download script uses `data["primary_symbol"]`, but `GeneReference` currently stores `symbol`.

Where:
- `src/data_prep/00c_download_mtDNA_NCBI_seqs.py`
- `src/utils/parsers.py`

Why this is critical:
- The script is not reliable for rebuilding the mtDNA raw CDS set from scratch.
- That breaks end-to-end reproducibility of the mtDNA branch of the pipeline.

What to do:
- Standardize on `symbol` everywhere, or add `primary_symbol` consistently in `GeneReference`.
- Re-test the mtDNA download step from a clean state.

### 3. The new TaxID mapping output is incompatible with `02_sanitize_all_alignments.py`

Current status:

- `Resolved / Historical` for the current stage layout.
- TaxID species overlap is now built by
  `src/data_download/00e_build_taxid_species_map.py`.
- Alignment sanitation is now handled by
  `src/data_curation/04_sanitize_all_alignments.py`, which scans current FASTA
  headers directly and writes `alignment_sanitation_manifest.parquet` plus
  `mt_accession_header_normalization.parquet`.
- The old `src/data_prep/00d_taxid_species_map.py` and
  `src/data_prep/02_sanitize_all_alignments.py` references are legacy names.

Problem:
The taxid overlap builder now writes one schema, but the sanitizer still expects the older accession-oriented schema.

Where:
- `src/data_prep/00d_taxid_species_map.py`
- `src/data_prep/02_sanitize_all_alignments.py`

Why this is critical:
- mtDNA header sanitation cannot be reliably regenerated from the current mapping file.
- This creates hidden dependence on previously generated files and undermines rerun reproducibility.

What to do:
- Decide which file is the source of truth:
  - a species overlap table for analysis, or
  - an accession-to-species remapping table for header sanitation.
- Either restore the older accession-based mapping file for sanitation, or update the sanitizer to use the current schema and another source for accession lookup.

### 4. Path conventions are inconsistent across rerunnable scripts

Current status:

- `Partly Resolved / Partly Historical`.
- The structural path issue named below is resolved:
  `src/structural/00_map_davs_to_structure.py` now uses
  `Path(__file__).resolve().parents[2]`.
- The old `src/data_prep` file references are legacy names in the current tree.
- Remaining caution: before a clean full rerun, re-check every active
  `src/data_download`, `src/data_curation`, `src/align`, `src/classify`,
  `src/structural`, and `src/phylo` script against the current emitted file
  names.

Problem:
Some scripts still assume obsolete filenames or locations, and one structural script assumes Docker-only root paths.

Where:
- `src/structural/00_map_davs_to_structure.py` uses `Path("/app")`
- `src/data_prep/00_fetch_reference_gene_coords.py` writes `reference_gene_coordinates{date}.tsv`
- `src/data_prep/01_curate_variants.py` expects `nucdna_gene_coordinates.tsv`
- `src/data_prep/00b_download_annotation_data.py` looks for the HGNC file under `data/`, while the HGNC download script writes to `data/reference/`

Why this is critical:
- These are classic rerun failures: the pipeline may work only because files already exist in the current workspace.

What to do:
- Standardize all paths and filenames.
- Replace hard-coded `/app` with `Path(__file__).resolve().parents[2]`.
- Make every downstream script consume the exact file emitted by the corresponding upstream step.

### 5. MyVariant download completeness is not guaranteed

Current status:

- `Resolved as an execution-safety issue; still a data-source limitation`.
- Current `src/data_download/00b_download_annotation_data.py` no longer does a
  blind broad MyVariant gene-batch fetch.
- Fresh MyVariant acquisition is explicitly skipped unless/until
  identity-based retrieval is implemented.
- The existing MyVariant-derived resource remains marked with
  `ok;completeness_legacy_unchecked`, so manuscript claims should not depend on
  complete MyVariant/dbNSFP/gnomAD coverage without a renewed acquisition pass.

Problem:
The MyVariant fetch batches genes and requests up to `size=1000`, but there is no check that all relevant variants were returned.

Where:
- `src/data_prep/00b_download_annotation_data.py`

Why this is critical:
- Missing MyVariant records will silently reduce coverage of gnomAD/dbNSFP-based annotations.
- This can distort severity filters and downstream stratification.

What to do:
- Add pagination or query by variant IDs rather than broad gene queries.
- At minimum, log total hits per query and fail when returned records equal the size cap.

## Should Fix Before Manuscript

### 6. Fallback amino-acid anchoring is still ambiguous in some cases

Current status:

- `Still open / narrowed`.
- Classification-stage QC currently passes, but the same general anchoring
  problem remains relevant to structural residue anchoring.
- Latest structural outputs contain `505` `mature_offset_candidate` rows and
  `154` `residue_anchoring_failure` rows in `dar_structure_map.csv`; these are
  not stage-breaking, but they must be preserved as structural coverage limits.

Problem:
When transcript/genomic maps are unavailable, the classifier falls back to finding the expected WT amino acid within a `±10 aa` window.

Where:
- `src/utils/alignment_parser.py`

Why it matters:
- If the same residue occurs more than once in the window, the anchor can resolve to the wrong position.
- This affects mtDNA and any nucDNA cases not covered by transcript remapping.

What to do:
- Use short peptide-context anchoring rather than a single residue.
- Record whether an anchor was unique or ambiguous.
- Consider classifying ambiguous anchors as unresolved instead of forcing a mapping.

### 7. Structural mapping collapses multiple copies of the same subunit to one chain

Current status:

- `Still open / manuscript robustness`.
- The current structural mapper produces substantial contact output, but
  duplicate protomer handling should still be audited before relying on
  protomer-specific contact claims.
- This is not a glaring exception for gene-level structural mapping, but it is
  a real caveat for contact-neighborhood interpretation.

Problem:
The structure cache stores one `gene_to_chain` value per gene, so duplicate protomers overwrite each other.

Where:
- `src/structural/00_map_davs_to_structure.py`

Why it matters:
- In symmetric or repeated assemblies, you are effectively choosing one arbitrary copy.
- Contact neighborhoods may differ slightly by local resolution, missing residues, or assembly context.

What to do:
- Store `gene -> [chain_ids]` rather than a single chain.
- Either map all copies and collapse later, or explicitly select the best-resolved copy by rule.

### 8. Isoform proxy mapping is useful but still a strong assumption

Current status:

- `Still open / sensitivity class`.
- The current mapper tags proxy-derived mappings with `proxy_mapping_used` and
  status categories such as `mapped_with_isoform_offset`.
- This is acceptable for exploratory structural mapping, but proxy mappings
  should remain separate from direct mappings in headline analyses.

Problem:
Some gene products are structurally mapped through a proxy isoform.

Where:
- `src/structural/00_map_davs_to_structure.py`

Why it matters:
- The assumption that numbering and local structure are interchangeable is not always safe, especially near termini or isoform-specific segments.

What to do:
- Keep the proxy system, but annotate proxy-derived mappings separately in all downstream summaries.
- Consider excluding proxy-based mappings from primary headline counts and treating them as sensitivity analyses.

### 9. One preferred PDB per complex is a simplification that can change contact calls

Current status:

- `Partly Resolved / still open for robustness`.
- The structural stage now uses a manifest-backed multi-model panel with
  primary, validation, and reference roles.
- The remaining issue is not a single-PDB-only design; it is that some models,
  especially respirasome/reference contexts, are not chain-complete substitutes
  for standalone complex structures.
- Multi-model support classes should be used in downstream structural claims.

Problem:
The structural pipeline chooses a single preferred structure for each complex.

Where:
- `src/structural/00_map_davs_to_structure.py`

Why it matters:
- Contact networks differ by conformational state, assembly state, resolution, bound cofactors, and missing loops.
- Some residues unresolved in one structure may be resolved in another.

What to do:
- For high-priority cases, validate contacts against at least one alternative structure or state.
- Add structure-quality metadata to mapped DARs.

### 10. Cross-dataset harmonization is still asymmetric between mtDNA and nucDNA

Current status:

- `Still open / interpretation`.
- Cross-source mtDNA ClinVar/MITOMAP overlap is now explicitly linked, and
  mitochondrial ClinVar rows are retained but excluded from the active nuclear
  comparative branch.
- The broader mtDNA-vs-nucDNA evidence asymmetry remains and should be handled
  by stratified reporting.

Problem:
The two branches of the pipeline rely on different source databases, transcript conventions, and population resources.

Where:
- throughout `src/data_prep`, `src/classify`, and `src/structural`

Why it matters:
- Even if execution is correct, biological comparisons between mtDNA and nucDNA can still reflect curation asymmetry rather than true biology.

What to do:
- Keep mtDNA vs nucDNA comparisons stratified by evidence quality.
- Explicitly define which comparisons are discovery-level and which are stringency-matched.

### 11. Translation masking is conservative but may reduce usable signal

Current status:

- `Still open / acceptable if quantified`.
- No evidence from the current notes indicates this has been fully quantified.
- Keep as a pre-manuscript sensitivity/statistics item.

Problem:
Any codon containing `N`, `X`, or `-` is translated to `X`.

Where:
- `src/align/00_translate_nucDNA.py`

Why it matters:
- This is safer than overcalling residues, but it can shrink usable species support near local alignment damage.

What to do:
- Keep this behavior, but quantify the fraction of alignment positions/species lost to masking.
- Report whether cDAV/uDAV calls are enriched among poorly resolved codon columns.

## Acceptable as Documented Limitations

### 12. The current cDAV definition is still discovery-oriented rather than stringent

Problem:
AA-level cDAV currently means the human disease amino acid is seen in at least one non-human species.

Why this can remain a limitation:
- It is a defensible discovery definition if you later stratify by support.
- It becomes problematic only when treated as equally strong evidence across one-species and many-species cases.

How to document it:
- Distinguish one-species, multi-species, and multi-origin cDAVs in all main analyses.

### 13. Structural contacts do not capture all compensation mechanisms

Problem:
The pipeline focuses on direct contact neighbors defined from static structures.

Why this can remain a limitation:
- Compensation may occur through second-shell residues, assembly effects, allostery, lipids, cofactors, or expression differences.

How to document it:
- Treat the structural partner analysis as one mechanistic layer, not an exhaustive definition of compensation.

### 14. Some variants will remain unresolved because reference, transcript, and alignment systems are not perfectly commensurate

Problem:
Even with transcript maps and genomic maps, some variants will not map cleanly into the chosen alignment framework.

Why this can remain a limitation:
- This is normal in multi-resource comparative genomics.

How to document it:
- Preserve unresolved categories explicitly.
- Never force unresolved cases into the uDAV set.

## Suggested Priority Order

### Immediate rerun blockers
- No active blocker from the original list is currently known for curation,
  classification, or structural mapping.
- Before a full clean rerun, verify active script path contracts end to end.
- Before reusing phylogenetic preparation scripts, replace or audit remaining
  legacy `locus.split("/")[0]` overlap handling.
- Treat MyVariant/dbNSFP/gnomAD as a legacy-completeness-limited annotation
  layer unless identity-based reacquisition is implemented.

### Next for analytical robustness
- Build a structural eligibility registry by gene/model/support class.
- Stratify structural analyses by `primary_only`, `primary_and_validation`,
  direct mapping, extended-offset rescue, and isoform-proxy support.
- Improve or explicitly freeze the mature-offset and residue-anchoring policy.
- Audit duplicate-protomer handling before making protomer-specific contact
  claims.
- Validate high-priority contacts against alternative structures or model roles.

### Can remain as stated limitations
- Discovery-level cDAV definition
- Static direct-contact interpretation
- Residual unresolved transcript/alignment cases
- Incomplete structural testability for genes absent from the selected human
  structure panel

## Practical Rule of Thumb

Before regenerating core outputs:
- all `data_prep` stages should run cleanly from scratch
- classification should not silently coerce overlap loci or ambiguous anchors
- structural mapping should run with repo-relative paths and deterministic chain-selection logic

Before writing the manuscript:
- unresolved but non-fatal issues should be either sensitivity-tested or explicitly disclosed
