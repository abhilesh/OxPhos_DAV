# OXPHOS DAV — Obsidian-Wiki Notes Missing from Logseq DB

- generated_at:: 2026-09-01T00:00:00+01:00
- project:: OXPHOS_DAV
- suggested_import_target:: [[Compensatory Variant Analysis]]
- suggested_workstream_targets:: [[Data curation]], [[Structural mapping pipeline]], [[Compensatory evolution pipeline]]
- source_graph:: Abhilesh_Knowledge_Graph
- source_wiki:: /Users/ad2347/Documents/Obsidian_Wiki/projects/OXPHOS_DAV/
- import_mode:: manual only
- database_policy:: The Logseq DB graph was inspected with read-only CLI commands and was not modified.
- timestamp_policy:: Dates below are source dates. Exact times are included only where the source page recorded them; no times have been invented.
- currentness_policy:: Later verification-ledger entries supersede older prose and earlier Logseq result snapshots.

## Audit conclusion

- timestamp:: 2026-09-01
  - finding:: The Logseq page [[Compensatory Variant Analysis]] is an early project outline, created 2026-05-11T23:32:01.580Z and updated 2026-05-21T09:10:12.796Z. It contains data-source notes, six structure models, three checkpoints, and links to three workstreams, but no current project status or verification ledger.
  - finding:: The Logseq workstream pages are almost empty: [[Data curation]] was created 2026-05-21T09:08:17.108Z, [[Structural mapping pipeline]] 2026-05-21T09:10:12.224Z, and [[Compensatory evolution pipeline]] 2026-05-21T09:05:14.876Z. Their useful content lives in linked reference blocks rather than in the workstream pages themselves.
  - finding:: The latest substantive Logseq reference found was dated 2026-07-23 and described a Klink et al. 2021 paper. The Obsidian wiki contains project updates through 2026-08-21, so the Logseq DB is missing the corrected late-July/August state.
  - finding:: The largest missing layer is not another method description; it is the current distinction between verified outputs, inherited claims, known-stale values, corrected candidate lists, and open work.

## Chronological missing-note stream

- timestamp:: 2026-07-08
  - note:: The durable project overview replaced the old manual-note reading order. The active project is /Users/ad2347/Documents/OxPhos_DAV; the pipeline runs from source acquisition and curation through alignment, classification, structural mapping, phylogenetics, and mutagenesis prioritization.
  - note:: The project asks whether human-pathogenic OXPHOS residues are naturally tolerated elsewhere in mammals and whether nearby structural residues provide co-evolutionary context. AA-level and exact-codon NT-level cDAV calls must remain separate.
  - note:: The live repository and current compiled wiki pages govern current state; historical manual-note folders are reference material only.
  - source:: projects/OXPHOS_DAV/OXPHOS_DAV.md; projects/OXPHOS_DAV/references/manual-notes-archive-plan.md

- timestamp:: 2026-07-09
  - note:: Legacy OXPHOS_DAR directories were archived to /Users/ad2347/Documents/Archives/oxphos_dar_legacy_2026-07-09.tar.zst with a checksum. The active OXPHOS_DAV project directory was left untouched.
  - note:: This separates the older OXPHOS_DAR clinical-tier pipeline from the active OXPHOS_DAV AA/NT cDAV pipeline; their numbers must not be merged.
  - source:: projects/OXPHOS_DAV/references/legacy-data-directory-archive-2026-07-09.md; projects/OXPHOS_DAV/references/cdar-reanalysis-2025-stratified-counts.md

- timestamp:: 2026-07-16
  - note:: Sixty-two historical OXPHOS_DAV notes were audited and their durable technical content was compiled into the wiki. Four adjacent research-program notes were excluded from the archive clearance.
  - note:: The audit identified a stale-number firewall: old notes contain incompatible pair totals, partner totals, Pagel totals, Pyvolve totals, contact totals, and timing distributions. Historical clinical Tier A/B/C/D labels also collide with later compensation-support tier names.
  - note:: The current structural contract, curation/classification contract, phylogenetic interpretation, filtering flags, gene caution list, and deferred-analysis backlog were created or consolidated to preserve the durable content.
  - source:: projects/OXPHOS_DAV/references/legacy-manual-notes-audit-2026-07-16.md; projects/OXPHOS_DAV/concepts/curation-and-classification-contract.md; projects/OXPHOS_DAV/concepts/structural-mapping-contract.md

- timestamp:: 2026-07-17
  - note:: Docker image 48a9069ef362 rebuilt successfully; core Python/R/IQ-TREE/MACSE smoke checks passed and the active tests passed 129/129 with pytest tests -q.
  - note:: Bare pytest recursively discovered archive/legacy_tests and failed on a removed utils.variant_annotators import. This was a test-discovery/configuration defect, not an active-suite failure; pytest.ini was subsequently changed but requires a fresh Docker verification.
  - note:: The cleared historical manual notes were confirmed archived. This was an environment and active-test validation, not a complete end-to-end freeze rerun.
  - source:: projects/OXPHOS_DAV/concepts/current-status-and-validation.md; projects/OXPHOS_DAV/references/manual-notes-archive-plan.md

- timestamp:: 2026-07-18
  - note:: Previously untracked repository content was committed across 6b3eb0a..beedd76; the P3 fix used _species_set() for list/string species-set representations. The malformed .gitignore rule .DS_Store.worktrees/ was corrected and local FoldX tooling was kept external/gitignored.
  - note:: The user explicitly chose to leave the abandoned feature/phylo-v2-overhaul branch/worktree untouched for this cleanup pass. That is a deferred merge/retirement decision, not adoption of the alternative.
  - source:: projects/OXPHOS_DAV/references/claude-cleanup-2026-07-18.md; projects/OXPHOS_DAV/concepts/decisions-and-caveats.md

- timestamp:: 2026-07-20
  - note:: Structural audit bugs P1–P9 were checked directly. P1/P2/P3/P4/P6/P7/P8/P9 were fixed or addressed; P5 was substantially addressed via the UniProt transit-peptide rescue. The final P7–P9 commit was f3b5384 and all 7 commits were pushed to origin/main.
  - note:: The structural-anchor registry expanded from SDHA/B/C/D to 31 UniProt-validated genes after 60 genes were queried. The 100-gene structure-panel registry was added, and the chain-assignment threshold became a named/logged 0.30 constant, although it remains empirically uncalibrated.
  - note:: Current model roles were clarified, especially CIV 9I7U as the primary mature NDUFA4-bound model and 9I6F as validation. Old tables that reverse those roles are superseded.
  - source:: projects/OXPHOS_DAV/references/claude-bug-verification-2026-07-20.md; projects/OXPHOS_DAV/concepts/structural-mapping-contract.md

- timestamp:: 2026-07-21
  - note:: Direct verification of curation, alignment, and classification found the current counts and contracts correct. The older REF_ALLELE_MISMATCH=19 claim is stale: current data has 2, and both are classified warnings rather than unresolved rows. Alignment masking is 0.50% overall and cDAV positions are modestly enriched on masked columns.
  - note:: Structural outputs were regenerated: 126 chain_assignment_or_model_gap rows split into 82 known gene-unrepresented panel gaps and 44 remaining failures; other mapping/contact totals were numerically unchanged.
  - note:: The corrected Pagel regeneration attempted all 4,227 pairs, obtained 896 valid p-values, and retained 47 significant pairs after global FDR correction. All four old named FoldX-era examples were absent from the corrected set; branch co-occurrence contributed zero pairs.
  - note:: A full phylo code sweep found a pre-existing 43-chunk Pagel pathway and a duplicate newer chunking implementation. The older merge path does not recompute branch_cooccur_fdr. A stale hardcoded 0.632 within-genome baseline was also identified for later correction.
  - note:: Re-running 00_prioritize_pairs.py exposed a silent tier1_pairs merge gap: up to 8/17 rubric points were zeroed. This explains why the old ranking and FoldX candidate list could not be trusted.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/concepts/decisions-and-caveats.md

- timestamp:: 2026-07-29
  - note:: All ten dependent phylogenetic outputs were re-run in dependency order against corrected structural/partner results and directly verified. The corrected outputs include 314/912 Pyvolve significance, 98.3%/94.3% mito-nuclear permissive timing, the rate-controlled nuclear-versus-mtDNA result, compartment analysis, revised contact-first, multi-origin analysis, and tier extraction.
  - note:: Tier extraction produced 14 Tier-1A, 33 Tier-1B, and 39 Tier-1C rows. The former raw Pagel top candidates are Tier-1B rather than Tier-1A; the full 14-pair Tier-1A list is recorded below.
  - note:: The corrected within-genome baseline is 0.551181. Historical 0.632 and 44-Pagel-pair statements must remain labelled as superseded.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/concepts/unified-compensation-scoring-alternatives.md

- timestamp:: 2026-08-20
  - note:: A number-free project-goals-and-analysis-arms page was added to clarify the project map. The canonical-results page was corrected to distinguish the 14-pair Tier-1A priority set from the wider 33-pair Tier-1B Pagel-significant set.
  - note:: Some other wiki prose still contains older mutagenesis or result wording. The later ledger entries and live result files take precedence over those stale passages.
  - source:: projects/OXPHOS_DAV/concepts/project-goals-and-analysis-arms.md; projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md

- timestamp:: 2026-08-21
  - note:: Origin coverage was rechecked on the current 47-pair set: 28 pairs have n_cdav_spp≥10 and all 28 resolve via timing_annotations_v2.csv. The older v1 timing file disagrees for 8/28 and must not supply current origin counts.
  - note:: The mutagenesis ranking code was fixed to merge the current v2 timing fields and verified by a live rerun. All 47 pairs scored; the maximum is 11/17, led by SDHD:132L↔104K. FoldX remains unrun for the current set and lambda is not staged.
  - note:: The current task queue and verification ledger were updated with these results, but the old 63-of-80 origin statement remains in a stale closing paragraph and must not be re-imported as current.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/concepts/phylogenetic-evidence-and-interpretation.md; projects/OXPHOS_DAV/tasks.md

## What is missing by Logseq area

- timestamp:: 2026-09-01
  - area:: Project orientation and governance
  - missing:: A current restatement of the two claims: c-DAR existence and the narrower, site-specific mechanism claim.
  - missing:: The AA-level versus exact-codon NT-level definitions, and the warning that NT support must be a subset of AA support.
  - missing:: The pipeline-verification-ledger rule: check current code and data before trusting a note or result number.
  - missing:: The current interpretation that the dominant signal is constitutive structural permissiveness under background purifying selection, not a universal reactive-rescue process.
  - wiki_sources:: projects/OXPHOS_DAV/concepts/project-goals-and-analysis-arms.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md

- timestamp:: 2026-09-01
  - area:: Data curation and classification
  - overlap:: Logseq already contains the historical coordinate-resolution logic: NM_ to ENST mapping, genomic-map fallback, allele orientation for MT-ND6, codon reconstruction, alignment masks, and explicit mismatch states.
  - missing:: The current filter-late state model: every curated row remains visible as classified, unresolved, or skipped_by_policy; unresolved is missing or discordant evidence, not a negative biological result.
  - missing:: The analytical-key rule that interpreted_gene must remain separate from structure/resource aliases such as NDUFA4/COXFA4.
  - missing:: Current curation and classification counts, clean/warning split, the corrected REF_ALLELE_MISMATCH count, and the measured alignment-masking audit.
  - missing:: The warning that clinical evidence tiers and compensation-support tiers are separate axes and must not share an unqualified Tier A/B/C/D label.
  - wiki_sources:: projects/OXPHOS_DAV/concepts/curation-and-classification-contract.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-09-01
  - area:: Alignment validation
  - overlap:: Logseq records codon-aware MACSE translation for mtDNA and TOGA-based protein translation for nuclear genes.
  - missing:: Current code/data verification: nuclear translation uses genetic-code table 1; mtDNA uses MACSE vertebrate-mitochondrial code table 2 via -gc_def 2.
  - missing:: The measured masking result: 0.50% of nuclear codon positions masked overall, 0.13%–1.68% per gene, versus 25.06% structural alignment gaps; cDAV positions are modestly enriched on masked columns (0.107% versus 0.054%, Mann–Whitney p=2.4e-12), but only 9.6% of classified rows touch any masked column.
  - wiki_sources:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-09-01
  - area:: Structural mapping
  - overlap:: Logseq records the three coordinate spaces, transit-peptide/mature-chain offsets, a structure manifest, eligibility ideas, and an older 94.8% mapping-rate snapshot.
  - missing:: The current active model policy: CI 9I4I; CII 8GS8; CIII 9HZL; CIV 9I7U primary with 9I6F validation; CV 8H9S primary with 8H9U and 8H9T validation priorities. Legacy 5XTH, 5XTE, and 5Z62 are reference/provenance models, not active denominator models.
  - missing:: The current UniProt transit-peptide rescue state: 60 genes queried, 31 validated, and 31 entries in the anchor-rescue registry; older SDH-only coverage is superseded.
  - missing:: The 100-gene structure-panel eligibility registry: 82 structurally_testable_direct, 5 structurally_testable_proxy, 9 structurally_unrepresented, and 4 no_eligible_variants (MT-CO1, MT-CO2, MT-CO3, MT-CYB).
  - missing:: The denominator rules: distinguish deliberate ineligibility from eligible-but-failed mapping; distinguish direct, registry-rescued, and proxy mappings; report primary/validation model and single/multi-model support.
  - missing:: The P1–P9 structural bug audit, including global FDR correction, cross-model contact deduplication, the named 0.30 chain-assignment threshold, and the known-gap versus unexpected-failure split.
  - missing:: Current structural output totals: 37,677 model-level mapping rows, 93,693 Cβ/Cβ ≤8 Å contact rows, 1,208 mito-nuclear contacts, and 9,375 tested DAR-contact pairs.
  - wiki_sources:: projects/OXPHOS_DAV/concepts/structural-mapping-contract.md; projects/OXPHOS_DAV/references/gene-specific-caution-list.md; projects/OXPHOS_DAV/references/claude-bug-verification-2026-07-20.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-09-01
  - area:: Compensatory-evolution and phylogenetic results
  - overlap:: Logseq has early method notes for Fisher screening, Pagel, Pyvolve, APC-MI, DCA, branch timing, and an abandoned SIMMAP/corHMM alternative. It also has a May result summary.
  - missing:: The corrected Pagel result: 4,227 pairs attempted, 896 valid p-values, and 47 significant pairs at global Pagel FDR ≤0.10. The earlier 462-pair result was inflated by the pre-fix per-DAR FDR scope.
  - missing:: Branch co-occurrence contributes 0 significant pairs globally in the corrected run; the original branch counts were retracted as a coordinate artifact.
  - missing:: The corrected Pyvolve result: 314/912 pairs significant at p<0.05, with 267/912 at p<0.01; the historical 310/912 is an older stochastic snapshot.
  - missing:: Corrected timing and directionality: within-genome baseline 0.551181; mito-nuclear timing is 98.3% permissive in the mt-DAR direction and 94.3% in the nuclear-DAR direction, with 918/937 combined events permissive and no supported directional asymmetry after Pyvolve calibration.
  - missing:: The 2026-07-29 result inventory for all ten dependent phylo outputs, including compartment analysis, rate-controlled Pyvolve, revised contact-first, multi-origin, and tier extraction.
  - missing:: The v1/v2 origin-count warning: among the current 47 significant pairs, 28 have n_cdav_spp≥10 and all 28 resolve using timing_annotations_v2.csv; the older timing_annotations.csv disagrees for 8/28, sometimes by up to sixfold, and must not supply current origin counts.
  - wiki_sources:: projects/OXPHOS_DAV/concepts/phylogenetic-evidence-and-interpretation.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md

- timestamp:: 2026-09-01
  - area:: Tiering and mutagenesis prioritization
  - missing:: The current tier extraction: 14 Tier-1A pairs, 33 Tier-1B pairs, and 39 Tier-1C pairs in an 86-row tier file.
  - missing:: Tier-1A is specifically Pagel-significant plus contact_first or permissive_background timing. It excludes co_adaptation, constitutively_permissive, and contact_after by the current implementation, even though some excluded categories may be scientifically interesting.
  - missing:: The corrected prioritization rerun: all 47 current pairs scored; maximum score 11/17; 5 pairs are under 5 Å; 42 are intraprotein and 5 interprotein; 43/47 have timing/physicochemical matches.
  - missing:: The current top scored pair is SDHD:132L↔104K, score 11/17, hbond, 3.83 Å, charge_rescue, co_occurring timing. Other score-10 examples are SDHD:132L↔130T, SDHB:228A↔227R, ATP5MF:78M↔82, and COX5A:55C↔59.
  - missing:: FoldX has not been run on the current 47-pair set. The historical 50/693 FoldX results cover the superseded ranking and are not current evidence.
  - wiki_sources:: projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/tasks.md

- timestamp:: 2026-09-01
  - area:: Project maintenance and open work
  - missing:: The current task queue: figure-generation scripts; reconciliation or retirement of duplicate Pagel chunking paths; FoldX for the current 47 pairs; extension and execution of the freeze registry; Docker pytest rerun; calibration of the 0.30 chain threshold; decision on the abandoned phylo-v2 worktree; and review of deferred extensions.
  - missing:: The current infrastructure boundary: FoldX is an external, gitignored dependency; lambda lacks the FoldX binary, current prioritized_pairs.csv, and data/structures/*.cif, and unattended SSH key access is blocked.
  - missing:: The legacy archive boundary: old OXPHOS_DAR directories were archived separately on 2026-07-09; the active OXPHOS_DAV project was not modified.
  - missing:: The project-wide manuscript rule: every quoted number needs a value, denominator, filter, source file, generating script, and analysis date.
  - wiki_sources:: projects/OXPHOS_DAV/tasks.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/references/legacy-data-directory-archive-2026-07-09.md; projects/OXPHOS_DAV/concepts/decisions-and-caveats.md

## Current project definition and analysis map

- timestamp:: 2026-08-20
  - definition:: A c-DAR/cDAV is a human-pathogenic OXPHOS amino-acid substitution that is naturally tolerated in at least one non-human mammal at the homologous position.
  - definition:: AA-level c-DAR means any non-human species carries the human disease amino acid.
  - definition:: NT-level c-DAR means a non-human species carries the exact same codon as the human mutation; this is stricter than AA-level support.
  - goal:: Establish whether c-DARs exist at meaningful scale, then test whether structural contact neighbours show co-evolutionary support for tolerance.
  - interpretation:: Existence is well supported. Mechanistic support is real but site-specific, strongest in mtDNA and in defined nuclear subsets; it is not a universal property of every cDAV-contact pair.
  - interpretation:: The dominant timing result supports constitutive permissiveness under background purifying selection, not a generic claim that a compensatory change arose reactively after a disease mutation.
  - source:: projects/OXPHOS_DAV/concepts/project-goals-and-analysis-arms.md; projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md

- timestamp:: 2026-08-20
  - analysis_arm:: Curation merges MITOMAP and ClinVar and maps transcript/genomic coordinates. Question: is there a sufficiently inventoried variant set?
  - analysis_arm:: Classification applies TOGA cross-species alignments and calls AA-cDAV and NT-cDAV. Question: which variants qualify?
  - analysis_arm:: Structural mapping places classified variants on cryo-EM models and identifies Cβ-Cβ ≤8 Å neighbours. Question: which c-DARs have a testable physical context?
  - analysis_arm:: Compensatory-partner testing applies Fisher screening, Pagel's discrete test, and branch co-occurrence. Question: is neighbour co-evolution enriched beyond chance?
  - analysis_arm:: Timing uses IQ-TREE ancestral states to classify contact_first, permissive_background, co_adaptation, constitutively_permissive, and contact_after events. Question: was tolerance pre-existing or reactive?
  - analysis_arm:: Pyvolve supplies a gene-specific neutral permutation null. Question: does permissiveness exceed tree-aware neutral expectation?
  - analysis_arm:: MI/APC-MI/plmDCA provide sequence-level coupling evidence. Question: does an independent sequence signal agree?
  - analysis_arm:: Tiering combines Pagel and timing evidence. Question: which pairs are the most defensible next candidates?
  - analysis_arm:: Mutagenesis prioritization combines geometry, physicochemical rescue, phylogenetic timing, and DCA/MI; FoldX is a downstream stability estimate. Question: which pairs should be tested experimentally?
  - source:: projects/OXPHOS_DAV/concepts/project-goals-and-analysis-arms.md

## Current verified result inventory

- timestamp:: 2026-08-21
  - curation:: 33,779 retained interpretation rows across 100 canonical genes: 13 mtDNA and 87 nuclear genes. An excluded-record inventory contains 8,875,040 rows.
  - classification_mtDNA:: 388 total; 323 classified; 0 unresolved; 65 skipped_by_policy; 167 AA-cDAV; 149 NT-cDAV.
  - classification_nucDNA:: 33,391 total; 6,486 classified; 68 unresolved; 26,837 skipped_by_policy; 2,046 AA-cDAV; 1,745 NT-cDAV.
  - classification_total:: 6,809 classified rows: 6,340 clean and 469 warning. All curated rows remain represented through explicit state fields.
  - classification_correction:: REF_ALLELE_MISMATCH is 2 in current data, not the stale 19 in the older issue note; both are classified with classified_warning rather than unresolved.
  - alignment:: Nuclear X-masking is 0.50% of codon positions overall, with a 0.13%–1.68% gene range; structural gaps are 25.06% and are a distinct phenomenon.
  - structure:: 37,677 model-level mapping rows; 93,693 Cβ/Cβ-8 Å contact rows; 1,208 mito-nuclear contact rows.
  - structure_status:: Regenerated structural mapping split 126 chain_assignment_or_model_gap rows into 82 gene_unrepresented_in_panel rows and 44 remaining unexpected/other failures; other mapped categories were numerically unchanged.
  - partners:: 9,375 DAR-contact pairs tested; 47 retained significant pairs after corrected global Pagel FDR; all 47 are Pagel-driven and branch_cooccur_fdr contributes 0.
  - pyvolve:: 912 conditional-permissiveness rows; 314/912=34.4% significant at perm_p<0.05; 267/912=29.3% at perm_p<0.01; median OR 1.73.
  - timing_mito_nuclear:: mt-DAR direction 26 pairs and 867 branch events, 852/867=98.3% permissive, p=2.2e-194 versus corrected baseline 0.551181; nuclear-DAR direction 20 pairs and 70 events, 66/70=94.3%, p=5.4e-13; combined 918/937=98.0%, p=3.1e-205.
  - compartments:: 41,373-pair universe: 809 with Pyvolve, 40,527 with APC-MI, and 47 with Pagel. The mito-nuclear APC-MI enrichment remains significant at p=2.83e-10; 166 mito-nuclear pairs are significant in at least one test.
  - rate_controlled:: Nuclear pairs retain higher permissive signal after rate control: Mann–Whitney U=69043, p=1.6e-05, median OR 1.736 nuclear versus 1.056 mtDNA; CMH χ²=7.45, p=6.3e-03.
  - tiers:: results/phylo/tier1_pairs.csv has 86 rows: 14 Tier-1A, 33 Tier-1B, and 39 Tier-1C.
  - origins:: Of the current 47 significant pairs, 28 have n_cdav_spp≥10 and all 28 have current origin counts in timing_annotations_v2.csv. The old v1 file disagrees for 8/28 and must not be used for current origin counts.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md

## Tier-1A current priority set

- timestamp:: 2026-07-29
  - definition:: Tier-1A is the current Pagel-significant plus permissive/informative timing set used for mutagenesis priority. It is not a claim that every pair has many independent origins or high timing confidence.
  - pair:: ATP5MF:78→M ↔ ATP5MF:80; pagel_fdr=9.3e-6; contact=nuc-nuc; timing=permissive_background; timing_confidence=low; independent_origins=5.
  - pair:: ATP5MF:93→Q ↔ ATP5MF:96; pagel_fdr=6.6e-3; contact=nuc-nuc; timing=permissive_background; timing_confidence=low; independent_origins=5.
  - pair:: SDHC:149→H ↔ SDHC:145; pagel_fdr=4.4e-2; contact=nuc-nuc; timing=permissive_background; timing_confidence=low; independent_origins=2.
  - pair:: SDHD:132→L ↔ SDHD:135; pagel_fdr=4.4e-2; contact=nuc-nuc; timing=permissive_background; timing_confidence=high; independent_origins=1.
  - pair:: SDHD:108→E ↔ SDHD:111; pagel_fdr=4.4e-2; contact=nuc-nuc; timing=permissive_background; timing_confidence=high; independent_origins=1.
  - pair:: SDHA:565→S ↔ SDHA:502; pagel_fdr=4.4e-2; contact=nuc-nuc; timing=contact_first; timing_confidence=high; independent_origins=1.
  - pair:: NDUFA9:23→T ↔ NDUFA9:253; pagel_fdr=5.1e-2; contact=nuc-nuc; timing=contact_first; timing_confidence=high; independent_origins=3.
  - pair:: SDHD:45→R ↔ SDHD:80; pagel_fdr=7.6e-2; contact=nuc-nuc; timing=permissive_background; timing_confidence=high; independent_origins=1.
  - pair:: SDHD:74→I ↔ SDHD:76; pagel_fdr=7.6e-2; contact=nuc-nuc; timing=contact_first; timing_confidence=low; independent_origins=4.
  - pair:: MT-ND2:122→A ↔ MT-ND2:123; pagel_fdr=8.2e-2; contact=mt-mt; timing=contact_first; timing_confidence=high; independent_origins=2.
  - pair:: UQCRB:13→R ↔ UQCRB:10; pagel_fdr=9.5e-2; contact=nuc-nuc; timing=permissive_background; timing_confidence=high; independent_origins=1.
  - pair:: SDHC:150→K ↔ SDHC:151; pagel_fdr=9.9e-2; contact=nuc-nuc; timing=permissive_background; timing_confidence=high; independent_origins=1.
  - pair:: SDHC:151→T ↔ SDHC:150; pagel_fdr=9.9e-2; contact=nuc-nuc; timing=permissive_background; timing_confidence=high; independent_origins=1.
  - pair:: SDHD:86→T ↔ SDHD:91; pagel_fdr=9.9e-2; contact=nuc-nuc; timing=permissive_background; timing_confidence=high; independent_origins=1.
  - caveat:: The tier definition deliberately excludes co_adaptation, constitutively_permissive, and contact_after. Those categories can still be scientifically important and should not be described as biologically weaker merely because they are Tier-1B.
  - caveat:: Several Tier-1A pairs have one independent origin or low timing confidence; single-origin pseudoreplication remains a caveat.
  - caveat:: FoldX has not been run on these 14 pairs.
  - source:: projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md; results/phylo/tier1_pairs.csv

## Tier-1B examples and why the old Logseq shortlist is unsafe

- timestamp:: 2026-07-29
  - definition:: Tier-1B is Pagel-significant but timing is absent or in a category excluded from the current Tier-1A definition. The raw Pagel top list is not automatically the mutagenesis-priority list.
  - pair:: ATP5MF:78V→M ↔ ATP5MF:82; pagel_fdr=8.7e-7; contact=hbond; contact_class=nuc-nuc.
  - pair:: ATP5MF:78V→M ↔ ATP5MF:80; pagel_fdr=9.3e-6; contact=vdw; contact_class=nuc-nuc.
  - pair:: ATP5MF:93R→Q ↔ ATP5MF:96; pagel_fdr=6.6e-3; contact=vdw; contact_class=nuc-nuc.
  - pair:: COX5A:55R→C ↔ COX5A:59; pagel_fdr=6.6e-3; contact=hbond; contact_class=nuc-nuc.
  - pair:: SDHB:228E→A ↔ SDHB:227; pagel_fdr=1.2e-2; contact=hbond; contact_class=nuc-nuc.
  - pair:: SDHD:108G→E ↔ SDHD:104; pagel_fdr=1.2e-2; contact=hbond; contact_class=nuc-nuc.
  - pair:: SDHA:6G→A ↔ SDHA:254; pagel_fdr=1.2e-2; contact=vdw; contact_class=nuc-nuc.
  - pair:: MT-ATP6:177A→T ↔ MT-ATP6:176; pagel_fdr=4.3e-2; contact=hbond; contact_class=mt-mt.
  - warning:: The former named examples MT-ATP6:192↔189, SDHC:146↔144, SDHC:46↔47, and COX6C:24↔25 are absent from the corrected 47-pair set. They are historical examples only.
  - warning:: The earlier Logseq result note reported 44 Pagel pairs, 310/912 Pyvolve pairs, 15 Tier-1A and 38 Tier-1C pairs, and a 50/693 FoldX pilot. Those are superseded snapshots, not current values.
  - source:: projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; Logseq block/page 34179; Logseq block/page 34171

## Curation and classification contract

- timestamp:: 2026-07-21
  - rule:: Keep all curated rows through classification and expose explicit state fields: classified, unresolved, or skipped_by_policy.
  - rule:: Treat unresolved as missing or discordant evidence, never as uDAV.
  - rule:: Use interpreted_gene as the analytical gene identity; preserve aliases and structure/resource names in separate provenance fields.
  - rule:: AA-level cDAV requires a readable non-human sequence carrying the disease amino acid. NT-level cDAV requires the exact disease codon.
  - rule:: classification_basis distinguishes nt_and_aa, nt_only, aa_only, and no_disease_allele_detected. The current classified snapshot has no nt_only rows, but the schema must not assume that state is impossible.
  - rule:: Select human TaxID 9606 explicitly; count biological positions on the ungapped human reference; reject partial-gap codons; mask -, !, *, and X; handle IUPAC ambiguity explicitly.
  - rule:: Use vertebrate mitochondrial code table 2 for mtDNA and standard table 1 for nuclear genes.
  - rule:: Complement minus-strand alleles, especially MT-ND6, before codon reconstruction while retaining CDS-relative and genomic alleles separately.
  - rule:: Preserve MITOMAP/ClinVar cross-source linkage rather than destructively deduplicating mitochondrial records.
  - historical_debugging:: The mismatch cascade was 862→274→104→69; a multiline FASTA parser fix raised later nuclear classified rows from 2,995 to 6,486 and reduced unresolved rows from 3,559 to 68. These are provenance facts, not current counts to substitute for the current snapshot.
  - source:: projects/OXPHOS_DAV/concepts/curation-and-classification-contract.md

## Structural mapping contract and gene cautions

- timestamp:: 2026-07-20
  - eligibility:: A row is structurally ineligible if it is skipped_by_policy, lacks complex/subunit assignment, lacks a usable TOGA/reference sequence, lacks an active model, or lacks a usable amino-acid coordinate.
  - reporting:: Any failure after eligibility must be reported as eligible-but-failed with a specific reason, not merged with deliberate ineligibility.
  - coordinate_policy:: Chain-to-gene assignment and position mapping must be sequence-backed; direct/±10 anchoring is default; wider ±80 searches are diagnostic unless a registry-approved offset exists.
  - contact_policy:: Default contact shell is Cβ–Cβ≤8 Å, using Cα for glycine. Contact labels are heuristic and must not be treated as biochemical proof.
  - sensitivity:: Report direct versus registry-rescued versus proxy mapping; primary versus validation model; single versus multi-model support; intra-chain versus inter-chain contacts; and 6 Å sensitivity where an 8 Å conclusion depends on the shell.
  - registry:: structure_panel_eligibility.tsv covers all 100 canonical genes: 82 direct, 5 proxy, 9 unrepresented, 4 no eligible variants.
  - no_eligible_variants:: MT-CO1, MT-CO2, MT-CO3, and MT-CYB are data-availability cases, not evidence of absent structural representation.
  - caution_genes:: High-value manual review list: NDUFV3, SDHB, NDUFS3, NDUFS4, CYC1, UQCRC1, UQCRC2, COX6A1, COX4I1, COX6B1, UQCRFS1, UQCRB, UQCRQ, UQCRH, COX5A, NDUFAB1.
  - source:: projects/OXPHOS_DAV/concepts/structural-mapping-contract.md; projects/OXPHOS_DAV/references/gene-specific-caution-list.md

## Phylogenetic method safeguards

- timestamp:: 2026-08-21
  - Fisher:: Species-level Fisher association is descriptive/screening evidence only because species are phylogenetically pseudoreplicated.
  - Pagel:: Pagel is the main phylogenetically rigorous correlated-binary-trait axis, but low or degenerate variation can make it unidentifiable or cause long runtimes; use guards, timeouts, and checkpointing.
  - branch_cooccurrence:: Matched-branch co-change is invalid if internal-node IDs from independent IQ-TREE runs are compared directly; the corrected current global result is 0 significant pairs.
  - Pyvolve:: Conditional permutation uses gene-specific neutral sequence evolution and is the appropriate tree-aware null for permissiveness; Q-matrix fallback and trait-prevalence dependence remain caveats.
  - MI:: Raw MI is ancestry- and indirect-correlation-confounded. APC-MI corrects positional averages but does not remove all phylogenetic or indirect effects. plmDCA needs large Meff relative to sequence length.
  - nuclear_DCA:: A global nuclear DCA null is primarily an Meff/power limitation, not evidence against nuclear co-evolution; signal appears in Meff≥10 strata, with 13/73 genes showing signal in that stratum.
  - directionality:: Do not describe mito-nuclear timing as one genome anticipating the other. The defensible interpretation is strong constitutive compatibility at conserved interfaces.
  - origins:: Species count is not independent-origin count. Use timing_annotations_v2.csv; for the current 47-pair population all 28 pairs with n_cdav_spp≥10 are resolved, but single-origin cases remain a pseudoreplication risk.
  - source:: projects/OXPHOS_DAV/concepts/phylogenetic-evidence-and-interpretation.md

## Filtering, manuscript, and interpretation flags

- timestamp:: 2026-07-20
  - F1:: MyVariant batch size is 1,000 and truncation verification remains a source-acquisition caveat; do not make primary gnomAD/dbNSFP claims without reacquisition.
  - F2:: ClinVar contains duplicated mtDNA records across GRCh37/38; use MITOMAP for mtDNA counts and ClinVar for nucDNA counts.
  - F3:: Structural panel gaps are now partly superseded by the 2026-07-20 eligibility registry; report known gaps separately from failures.
  - F4:: Eligible mapping failures must be reported separately from intentional ineligibility.
  - F5:: Fisher results are pseudoreplicated by species inheritance.
  - F6:: Pagel has many NA cases for MT-MT or degenerate traits; use branch/Pyvolve where appropriate and do not compare raw Pagel values across compartments without care.
  - F7:: A very large contact_first odds ratio may be driven by one origin; validate with independent-origin and Pyvolve analyses.
  - F8:: Tier labels are for final analysis prioritization, not an early filter on the pipeline.
  - F9:: COXFA4/NDUFA4 is an alias-display issue; retain COXFA4 as the analytical key.
  - F10:: likely_incompatible is a mechanistic flag, not a significance result.
  - source:: projects/OXPHOS_DAV/references/final-analysis-filtering-flags.md; projects/OXPHOS_DAV/concepts/decisions-and-caveats.md

## Open work and deferred ideas

- timestamp:: 2026-08-21
  - open:: Add or locate reproducible figure-generation scripts and map result tables to manuscript figures.
  - open:: Reconcile or retire the duplicate Pagel chunking infrastructure. The pre-existing prepare_pagel_hpc.py pathway creates 43 chunks at chunk size 100, but its merge step does not recompute branch_cooccur_fdr; the newer self-contained path and the older path should not be silently conflated.
  - open:: Run FoldX ΔΔG for the current 47-pair set, beginning with the corrected prioritized_pairs.csv ranking. Lambda still needs the FoldX binary, current prioritized_pairs.csv, and structure CIFs staged; SSH key access is also an infrastructure issue.
  - open:: Extend pipeline_audit_registry.tsv to structural-partner, phylo, and mutagenesis outputs and rerun tools/check_pipeline_freeze.py.
  - open:: Run pytest tests -q and bare pytest inside the current Docker image to verify the pytest.ini/P3/P7–P9 fixes together.
  - open:: Empirically calibrate the named 0.30 chain-assignment identity threshold with a labeled correct/incorrect dataset or manual review of divergent short subunits.
  - open:: Decide whether to retire or merge the abandoned feature/phylo-v2-overhaul branch/worktree; it contains a SIMMAP+corHMM alternative with only small smoke tests and was never adopted by main.
  - open:: Complete or explicitly scope out the remaining ~643 of the old 693-pair FoldX plan; this is separate from rerunning FoldX on the current 47 pairs.
  - deferred:: Six-Å structural sensitivity; solvent accessibility, functional annotation, entropy, and alignment-confidence covariates; matched nulls; mutation-rate controls; LOEUF/pLI; carefully designed dN/dS/BUSTED/RELAX analyses; 3–5-position compensatory networks; and disease-only/contact-only/double-mutant experimental designs.
  - promotion_rule:: Do not promote a deferred analysis without a current script or assay plan, frozen input cohort, explicit expected output, and a claim-to-output ledger entry.
  - source:: projects/OXPHOS_DAV/tasks.md; projects/OXPHOS_DAV/concepts/deferred-analysis-and-validation-backlog.md

## Historical notes that must remain labelled as historical

- timestamp:: 2026-07-16
  - legacy_counts:: Old manual notes contain incompatible snapshots: 8,087 versus 9,375 tested pairs; 693/819 versus 462 versus 47 retained partners; 37 versus 44 versus 47 Pagel pairs; 3,673/232 versus 912/310 or 314 Pyvolve results; and several incompatible timing distributions.
  - legacy_semantics:: Older NT-cDAR usage sometimes meant an amino-acid match rather than exact codon identity. The current NT definition is exact codon identity.
  - legacy_candidates:: April provisional examples were MT-ATP6 I192T/T189A; MT-ATP6 I192T/S188T; NDUFA9 M220I/S216Y; MT-ATP6 T28S/MT-ATP8 F24; and SDHC I46T/G47S. They are not current candidates.
  - archive:: Sixty-two historical OXPHOS_DAV notes were audited and durable content compiled into the wiki; the cleared old manual notes were archived on 2026-07-17. Four adjacent research-program notes were excluded from that clearance.
  - archive_boundary:: The legacy OXPHOS_DAR data directories were archived on 2026-07-09 to /Users/ad2347/Documents/Archives/oxphos_dar_legacy_2026-07-09.tar.zst with checksum; the active OXPHOS_DAV directory remained untouched.
  - source:: projects/OXPHOS_DAV/references/legacy-manual-notes-audit-2026-07-16.md; projects/OXPHOS_DAV/references/provisional-april-2026-results-and-targets.md; projects/OXPHOS_DAV/references/legacy-data-directory-archive-2026-07-09.md

- timestamp:: 2026-07-18
  - repository_state:: Previously untracked project content was committed across six commits 6b3eb0a..beedd76; the final P7–P9 fix commit f3b5384 was pushed on 2026-07-20.
  - cleanup:: The malformed .gitignore line .DS_Store.worktrees/ was split into working ignore rules, with local settings, .superpowers, and tools/foldx also ignored.
  - branch_decision:: The user chose to leave the phylo-v2 branch/worktree untouched for the cleanup pass; this is a deferred decision, not evidence that the alternative is live.
  - testing:: Docker image 48a9069ef362 passed the active suite at 129 tests on 2026-07-17, but the bare pytest discovery issue was only later addressed in pytest.ini and still requires a fresh Docker verification.
  - source:: projects/OXPHOS_DAV/references/claude-cleanup-2026-07-18.md; projects/OXPHOS_DAV/concepts/current-status-and-validation.md

## Source-page timestamp index

- timestamp:: 2026-08-20T17:37:00+01:00
  - source_page:: projects/OXPHOS_DAV/OXPHOS_DAV.md
  - page_created:: 2026-07-08T18:35:00+01:00
  - page_updated:: 2026-08-20T17:37:00+01:00

- timestamp:: 2026-08-21T13:00:00+01:00
  - source_page:: projects/OXPHOS_DAV/tasks.md
  - page_created:: 2026-07-15
  - page_updated:: 2026-08-21T13:00:00+01:00

- timestamp:: 2026-08-21T13:00:00+01:00
  - source_page:: projects/OXPHOS_DAV/tasks-done.md
  - page_created:: 2026-07-15
  - page_updated:: 2026-08-21T13:00:00+01:00

- timestamp:: 2026-08-21T13:30:00+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md
  - page_created:: 2026-07-21T00:00:00+01:00
  - page_updated:: 2026-08-21T13:30:00+01:00

- timestamp:: 2026-08-21T13:00:00+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/phylogenetic-evidence-and-interpretation.md
  - page_created:: 2026-07-16T12:25:47+01:00
  - page_updated:: 2026-08-21T13:00:00+01:00

- timestamp:: 2026-08-20T17:30:00+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md
  - page_created:: 2026-07-20T22:00:00+01:00
  - page_updated:: 2026-08-20T17:30:00+01:00

- timestamp:: 2026-08-20T17:00:00+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/current-status-and-validation.md
  - page_created:: 2026-07-08T18:35:00+01:00
  - page_updated:: 2026-08-20T17:00:00+01:00

- timestamp:: 2026-08-20T17:37:00+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/project-goals-and-analysis-arms.md
  - page_created:: 2026-08-20T17:37:00+01:00
  - page_updated:: 2026-08-20T17:37:00+01:00

- timestamp:: 2026-07-21T19:30:00+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/decisions-and-caveats.md
  - page_created:: 2026-07-08T18:35:00+01:00
  - page_updated:: 2026-07-21T19:30:00+01:00

- timestamp:: 2026-07-20T21:00:00+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/curation-and-classification-contract.md
  - page_created:: 2026-07-16T12:25:47+01:00
  - page_updated:: 2026-07-20T21:00:00+01:00

- timestamp:: 2026-07-20T21:00:00+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/structural-mapping-contract.md
  - page_created:: 2026-07-16T12:25:47+01:00
  - page_updated:: 2026-07-20T21:00:00+01:00

- timestamp:: 2026-07-29T09:45:14+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/unified-compensation-scoring-alternatives.md
  - page_created:: 2026-07-29T09:45:14+01:00
  - page_updated:: 2026-07-29T09:45:14+01:00

- timestamp:: 2026-07-29T00:00:00+01:00
  - source_page:: projects/OXPHOS_DAV/concepts/uniprot-offset-validation-methodology.md
  - page_created:: 2026-07-29T00:00:00+01:00
  - page_updated:: 2026-07-29T00:00:00+01:00

- timestamp:: 2026-07-08T19:20:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/claude-audit-2026-07-08.md
  - page_created:: 2026-07-08T19:20:00+01:00
  - page_updated:: 2026-07-08T19:20:00+01:00

- timestamp:: 2026-07-20T20:25:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/claude-bug-verification-2026-07-20.md
  - page_created:: 2026-07-20T20:00:00+01:00
  - page_updated:: 2026-07-20T20:25:00+01:00

- timestamp:: 2026-07-18T14:18:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/claude-cleanup-2026-07-18.md
  - page_created:: 2026-07-18T14:18:00+01:00
  - page_updated:: 2026-07-18T14:18:00+01:00

- timestamp:: 2026-07-20T21:15:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/final-analysis-filtering-flags.md
  - page_created:: 2026-07-20T21:15:00+01:00
  - page_updated:: 2026-07-20T21:15:00+01:00

- timestamp:: 2026-07-20T21:15:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/gene-specific-caution-list.md
  - page_created:: 2026-07-20T21:15:00+01:00
  - page_updated:: 2026-07-20T21:15:00+01:00

- timestamp:: 2026-07-16T12:00:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/legacy-data-directory-archive-2026-07-09.md
  - page_created:: 2026-07-15T12:00:00+01:00
  - page_updated:: 2026-07-15T12:00:00+01:00

- timestamp:: 2026-07-17T11:58:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/legacy-manual-notes-audit-2026-07-16.md
  - page_created:: 2026-07-16T12:25:47+01:00
  - page_updated:: 2026-07-17T11:58:00+01:00

- timestamp:: 2026-07-17T11:58:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/manual-notes-archive-plan.md
  - page_created:: 2026-07-08T18:35:00+01:00
  - page_updated:: 2026-07-17T11:58:00+01:00

- timestamp:: 2026-07-20T21:15:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/phylo-v2-overhaul-design.md
  - page_created:: 2026-07-20T21:15:00+01:00
  - page_updated:: 2026-07-20T21:15:00+01:00

- timestamp:: 2026-07-20T21:15:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/project-contract-and-scope-history.md
  - page_created:: 2026-07-20T21:15:00+01:00
  - page_updated:: 2026-07-20T21:15:00+01:00

- timestamp:: 2026-07-16T12:25:47+01:00
  - source_page:: projects/OXPHOS_DAV/references/provisional-april-2026-results-and-targets.md
  - page_created:: 2026-07-16T12:25:47+01:00
  - page_updated:: 2026-07-16T12:25:47+01:00

- timestamp:: 2026-07-29T00:00:00+01:00
  - source_page:: projects/OXPHOS_DAV/references/cdar-reanalysis-2025-stratified-counts.md
  - page_created:: 2026-07-29T00:00:00+01:00
  - page_updated:: 2026-07-29T00:00:00+01:00

- timestamp:: 2026-07-20T20:25:00+01:00
  - source_page:: projects/OXPHOS_DAV/skills/canonical-pipeline-runbook.md
  - page_created:: 2026-07-08T18:35:00+01:00
  - page_updated:: 2026-07-20T20:25:00+01:00

## Manual import note

- timestamp:: 2026-09-01
  - instruction:: Import this file manually into the Logseq DB graph, or split its dated top-level sections into the three existing workstream pages and the main [[Compensatory Variant Analysis]] page.
  - instruction:: Keep the historical section labelled historical and keep the current result inventory near the top of the main page.
  - instruction:: Do not overwrite the existing Logseq DB page or its backing data automatically.
  - verification:: This Markdown file is a derived import aid; it is not a replacement for the Obsidian wiki or the live repository.
