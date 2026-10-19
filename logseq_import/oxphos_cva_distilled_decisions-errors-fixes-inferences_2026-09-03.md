# OXPHOS DAV — Distilled Decisions, Errors, Fixes, and Inferences

- generated_at:: 2026-09-03T00:00:00+01:00
- evidence_cutoff:: 2026-08-21
- project:: OXPHOS_DAV
- source_graph:: Abhilesh_Knowledge_Graph
- suggested_import_target:: [[Compensatory Variant Analysis]]
- suggested_workstream_targets:: [[Data curation]], [[Structural mapping pipeline]], [[Compensatory evolution pipeline]]
- source_wiki:: /Users/ad2347/Documents/Obsidian_Wiki/projects/OXPHOS_DAV/
- companion_to:: /Users/ad2347/Documents/OxPhos_DAV/logseq_import/oxphos_cva_missing_wiki_notes_2026-09-01.md
- import_mode:: manual only
- database_policy:: The Logseq DB and its backing data were not modified.
- currentness_policy:: These notes summarize the latest verified wiki state; a future rerun may change result values.

## Important decisions

- timestamp:: 2026-07-08
  - decision:: Treat the live OXPHOS_DAV repository plus the verification ledger as the operational source of truth. Historical manual notes are context, not current instructions.
  - decision:: Record new project findings in the Obsidian wiki and verification ledger, not in repo notes/docs that can silently drift from the wiki.
  - rationale:: The project had already accumulated prose whose numbers outlived the code and output state.
  - source:: projects/OXPHOS_DAV/OXPHOS_DAV.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-09
  - decision:: Keep the active OXPHOS_DAV project separate from the legacy OXPHOS_DAR pipeline and archive the old directories without modifying active data.
  - decision:: Do not merge OXPHOS_DAR clinical-tier results with OXPHOS_DAV AA/NT cDAV results.
  - rationale:: The projects use different filters, definitions, and denominators.
  - source:: projects/OXPHOS_DAV/references/legacy-data-directory-archive-2026-07-09.md; projects/OXPHOS_DAV/references/cdar-reanalysis-2025-stratified-counts.md

- timestamp:: 2026-07-16
  - decision:: Use a filter-late pipeline: retain curated records and label them classified, unresolved, or skipped_by_policy rather than silently dropping them.
  - decision:: Treat unresolved as missing or discordant evidence, never as uDAV.
  - decision:: Keep interpreted_gene separate from structure/resource aliases such as COXFA4/NDUFA4.
  - decision:: Keep clinical evidence, compensation support, and structural evidence as separate axes; never use an unqualified Tier A/B/C/D label.
  - source:: projects/OXPHOS_DAV/concepts/curation-and-classification-contract.md; projects/OXPHOS_DAV/concepts/decisions-and-caveats.md

- timestamp:: 2026-07-16
  - decision:: Define AA-level cDAV as a non-human species carrying the disease amino acid, and NT-level cDAV as a non-human species carrying the exact disease codon.
  - decision:: Preserve AA and NT results separately; NT support should be a subset of AA support.
  - rationale:: An amino-acid match can arise through a different codon, so NT-level evidence is stricter.
  - source:: projects/OXPHOS_DAV/concepts/curation-and-classification-contract.md

- timestamp:: 2026-07-20
  - decision:: Use sequence-backed coordinate resolution across clinical transcript, TOGA alignment, and structure-chain spaces.
  - decision:: Use direct/±10 anchoring by default; accept wider offsets only through an explicit registry.
  - decision:: Treat direct mapping, registry-rescued mapping, proxy mapping, model role, and mapping failure as separate evidence states.
  - source:: projects/OXPHOS_DAV/concepts/structural-mapping-contract.md

- timestamp:: 2026-07-20
  - decision:: Use the active structure manifest rather than historical structure tables: CI 9I4I, CII 8GS8, CIII 9HZL, CIV 9I7U primary/9I6F validation, and CV 8H9S primary with 8H9U/8H9T validation priorities.
  - decision:: Use Cβ–Cβ≤8 Å contacts, Cα for glycine, and report a 6 Å sensitivity analysis when conclusions depend on the shell.
  - source:: projects/OXPHOS_DAV/concepts/structural-mapping-contract.md

- timestamp:: 2026-07-21
  - decision:: Use globally corrected Pagel FDR as the main phylogenetically rigorous partner-test axis.
  - decision:: Use species-level Fisher only as descriptive screening because species are phylogenetically pseudoreplicated.
  - decision:: Treat branch co-occurrence cautiously because internal node IDs from independent IQ-TREE runs are not directly comparable.
  - decision:: Use Pyvolve as a tree-aware neutral permutation control and interpret DCA/MI according to Meff and ancestry limitations.
  - source:: projects/OXPHOS_DAV/concepts/phylogenetic-evidence-and-interpretation.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-18
  - decision:: Leave the abandoned feature/phylo-v2-overhaul branch/worktree untouched for now; do not treat the SIMMAP+corHMM alternative as the adopted pipeline.
  - rationale:: Main continues to use Pagel, branch timing, and IQ-TREE ancestral-state analysis; the alternative had only small smoke tests.
  - source:: projects/OXPHOS_DAV/references/phylo-v2-overhaul-design.md; projects/OXPHOS_DAV/concepts/decisions-and-caveats.md

- timestamp:: 2026-07-29
  - decision:: Keep the newer self-contained Pagel chunking and the older prepare_pagel_hpc.py pathway under deliberate review rather than silently deleting or conflating either.
  - rationale:: The older path is real and produces 43 chunks, but its merge step does not recompute branch_cooccur_fdr.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/tasks.md

## Errors encountered and fixes implemented

- timestamp:: 2026-07-17
  - error:: Bare pytest recursively collected archive/legacy_tests and failed on the removed utils.variant_annotators import.
  - fix:: Restrict pytest discovery to tests/ and exclude archive/.worktrees/.devcontainer in pytest.ini.
  - status:: The active suite previously passed 129/129; the combined bare-pytest and active-suite behavior still requires a fresh Docker verification.
  - source:: projects/OXPHOS_DAV/references/claude-cleanup-2026-07-18.md; projects/OXPHOS_DAV/concepts/current-status-and-validation.md

- timestamp:: 2026-07-18
  - error:: P3 used lineages_with_disease_allele without guarding list/string representations.
  - fix:: Add _species_set() coercion at both structural-partner call sites.
  - follow_up:: A numpy-array representation later exposed a third case; that was fixed during the Pagel regeneration.
  - source:: projects/OXPHOS_DAV/references/claude-cleanup-2026-07-18.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-20
  - error:: Structural audit P1 mutated status with proxy text; P2 ignored the registry max offset; P4 allowed duplicate contacts across models; P5/P6 mishandled transcript-to-chain offsets and restricted rescue to SDH; P7 applied FDR per DAR while testing a global threshold; P8 hid an uncalibrated 0.30 chain identity cutoff; P9 had no empirical panel-eligibility registry.
  - fix:: Remove status mutation, pass registry offsets, deduplicate physical contacts, add UniProt transit-peptide rescue, globally correct FDR, name/log the identity threshold, and add the 100-gene eligibility registry.
  - status:: All P1–P9 were fixed or substantially addressed; P8's threshold value remains uncalibrated.
  - source:: projects/OXPHOS_DAV/references/claude-bug-verification-2026-07-20.md; projects/OXPHOS_DAV/concepts/structural-mapping-contract.md

- timestamp:: 2026-07-20
  - error:: Initial UniProt offset work was blocked by network restrictions and showed that several genes have multiple conflicting offset candidates, making naive gene-wide rescue unsafe.
  - fix:: Use the later completed validation pass: 60 genes queried, 31 validated by UniProt transit-peptide and TOGA cross-checks, with the 31 validated entries added to the rescue registry.
  - source:: projects/OXPHOS_DAV/concepts/uniprot-offset-validation-methodology.md; projects/OXPHOS_DAV/concepts/structural-mapping-contract.md

- timestamp:: 2026-07-21
  - error:: A stale classification note reported REF_ALLELE_MISMATCH=19 and treated those cases as unresolved.
  - fix:: Recheck current Parquet output and correct the claim to 2; both current rows are classified warnings, not unresolved rows.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/concepts/curation-and-classification-contract.md

- timestamp:: 2026-07-21
  - error:: Pagel/partner outputs were stale and the earlier per-DAR FDR scope inflated the retained set; the first rerun attempt also hit a numpy-array species-set bug and a wall-time limit.
  - fix:: Add real Pagel chunking, fix the array representation, rerun all 4,227 pairs, and apply global FDR.
  - result:: 896 valid p-values and 47 significant pairs, replacing the stale 462-pair set; all four old named examples disappeared, and branch co-occurrence contributed zero pairs.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-21
  - error:: phylo sibling scripts still contained unguarded species-set access and locus.split("/")[0] overlap-gene fallbacks.
  - fix:: Patch both latent code paths; the fallback had not fired for the 2,213 classified cDAVs, so this was a risk closure rather than a known output correction.
  - source:: projects/OXPHOS_DAV/concepts/decisions-and-caveats.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-21
  - error:: 06b_temporal_ordering_mitonuc.py hardcoded a 0.632 within-genome baseline based on the old 44-pair result.
  - fix:: Recompute the baseline after the corrected phylo rerun; current value is 0.551181.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-21
  - error:: 00_prioritize_pairs.py never merged tier1_pairs.csv, silently assigning zero to up to 8 of 17 rubric points and producing an untrustworthy ranking.
  - fix:: Merge timing_annotations_v2.csv for current timing/origin/confidence fields and use the older timing file only for deterministic physicochemical_type where no v2 equivalent exists; switch to timing_confidence="high".
  - result:: The live rerun scored all 47 current pairs and raised the maximum score from 7 to 11.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-29
  - error:: Eight dependent phylo files and two always-run conditional-permissiveness outputs were stale after the structural correction; the Pyvolve array did not generate the latter two files.
  - fix:: Rerun all ten outputs in dependency order, including revised contact-first and multi-origin analyses.
  - result:: Current outputs are directly verified rather than merely traced as rerun.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-18
  - error:: The repository had previously untracked content and a malformed .gitignore rule joining .DS_Store and .worktrees.
  - fix:: Commit the project content in six commits, repair the ignore rules, and push the final P7–P9 fix commit f3b5384 to origin/main on 2026-07-20.
  - source:: projects/OXPHOS_DAV/references/claude-cleanup-2026-07-18.md

- timestamp:: 2026-08-21
  - error:: Lambda was not ready for FoldX: no FoldX binary, current prioritized_pairs.csv, or structure CIFs were staged; SSH key authentication was also rejected outside an existing ControlMaster session.
  - fix:: None yet; stage the three required inputs and resolve unattended SSH access before sbatch.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/tasks.md

## Important inferences

- timestamp:: 2026-08-20
  - inference:: The existence claim is strong: cDAVs occur at substantial rates among classified pathogenic residues, higher in mtDNA than nuclear DNA and higher at AA than exact-codon NT level.
  - evidence:: mtDNA 167/323 AA-cDAV and 149/323 NT-cDAV; nucDNA 2,046/6,486 AA-cDAV and 1,745/6,486 NT-cDAV.
  - source:: projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md

- timestamp:: 2026-08-20
  - inference:: The mechanism claim is real but site-specific, not universal. The corrected Pagel set is 47 pairs, not hundreds of broadly supported contacts.
  - evidence:: Pagel found 47 significant pairs from 896 valid p-values after 4,227 attempts; SDHD/SDHC contribute about 40%, but non-SDH pairs remain across CI, CIII, CIV, and CV.
  - source:: projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-29
  - inference:: Timing primarily supports constitutive structural permissiveness under background purifying selection, not reactive compensation after the disease allele.
  - evidence:: Within-genome timing is based on corrected baseline 0.551181; mito-nuclear events are 98.3% permissive in the mt-DAR direction and 94.3% in the nuclear-DAR direction, with no supported directional asymmetry after Pyvolve calibration.
  - source:: projects/OXPHOS_DAV/concepts/phylogenetic-evidence-and-interpretation.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-29
  - inference:: The Pyvolve enrichment is independently supportive of a non-random permissiveness signal, but it does not prove organismal fitness compensation.
  - evidence:: 314/912 pairs are significant at perm_p<0.05 versus an approximately 5% null expectation; odds ratios still depend on trait prevalence.
  - source:: projects/OXPHOS_DAV/concepts/phylogenetic-evidence-and-interpretation.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-21
  - inference:: The nuclear DCA null is mainly a power/Meff problem, not evidence that nuclear co-evolution is absent.
  - evidence:: Signal appears in Meff≥10 strata; mtDNA DCA-DI enrichment replicates across 8/9 genes, while most nuclear alignments are underpowered.
  - source:: projects/OXPHOS_DAV/concepts/phylogenetic-evidence-and-interpretation.md

- timestamp:: 2026-08-21
  - inference:: Species count is not independent-origin count and is unsafe as a proxy for phylogenetic replication.
  - evidence:: All 28 current significant pairs with n_cdav_spp≥10 resolve in timing_annotations_v2.csv, but the old v1 file disagrees for 8/28, up to sixfold; single-origin cases remain common.
  - source:: projects/OXPHOS_DAV/concepts/phylogenetic-evidence-and-interpretation.md; projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

- timestamp:: 2026-07-20
  - inference:: A structural non-mapping result is not biological evidence that compensation is absent.
  - evidence:: Mapping failures mix deliberate policy exclusions, absent model representation, isoform incompatibility, unresolved residues, and genuine eligible failures; the new registry separates these categories.
  - source:: projects/OXPHOS_DAV/concepts/structural-mapping-contract.md

- timestamp:: 2026-08-20
  - inference:: Tier-1A is a framework-specific priority class, not a universal ranking of biological interest.
  - evidence:: It includes Pagel-significant pairs with contact_first or permissive_background timing but excludes co_adaptation, constitutively_permissive, and contact_after; some Tier-1B pairs may be more direct compensation hypotheses.
  - source:: projects/OXPHOS_DAV/concepts/canonical-results-and-goal-interpretation.md

- timestamp:: 2026-08-21
  - inference:: The current mutagenesis ranking is computational prioritization, not experimental validation.
  - evidence:: All 47 current pairs are scored, maximum 11/17, but FoldX has not been run on the current set and no wet-lab rescue has been demonstrated.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/concepts/decisions-and-caveats.md

- timestamp:: 2026-08-21
  - inference:: The most important reproducibility risk is temporal drift between code, outputs, and prose—not necessarily a flawed analysis method.
  - evidence:: Direct checks uncovered stale rankings, stale candidate examples, a silent merge gap, duplicated chunking paths, and stale summaries even after several items had been described as fixed.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md

## Current compact status

- timestamp:: 2026-08-21
  - status:: 33,779 curated rows across 100 canonical genes; 6,809 classified rows; 6,340 clean and 469 warning.
  - status:: Structural outputs: 37,677 model mappings, 93,693 contacts, 1,208 mito-nuclear contacts, and 9,375 tested pairs.
  - status:: Phylo outputs: 47 corrected Pagel-significant pairs; 314/912 Pyvolve-significant pairs; 14 Tier-1A, 33 Tier-1B, and 39 Tier-1C rows.
  - status:: Mutagenesis: 47/47 pairs scored; top SDHD:132L↔104K at 11/17; FoldX current-set run still open.
  - status:: Remaining major work: FoldX staging/run, freeze-registry coverage, figure scripts, duplicate Pagel-path decision, Docker test rerun, threshold calibration, and phylo-v2 branch decision.
  - source:: projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md; projects/OXPHOS_DAV/tasks.md

## Import instruction

- timestamp:: 2026-09-03
  - instruction:: This is the distilled companion to the comprehensive 2026-09-01 import bundle. Import it as a compact decision-and-results layer, then use the longer bundle for detailed candidate tables and page timestamps.
  - instruction:: Keep each timestamp:: block as a separate Logseq block if you want to query decisions, errors, fixes, and inferences independently.
  - instruction:: Do not automatically overwrite the existing Logseq DB page or workstream pages.
