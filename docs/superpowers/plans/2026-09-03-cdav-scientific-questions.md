# cDAV Scientific Questions Research Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the current cDAV framework into a reproducible set of bias-controlled analyses that distinguish evolutionary permissiveness, structural mechanism, recurrence, higher-order interactions, and translational usefulness.

**Architecture:** Build one analysis contract at the DAR-contact-pair level, then run pre-specified analyses in three layers: (1) prevalence and null calibration, (2) structural mechanism and evolutionary timing, and (3) network/predictive extensions. Generate manuscript figures only from versioned analysis tables, with every headline number traceable to a source file and filter definition.

**Tech Stack:** Existing Docker image `oxphos_dav_analysis`; Python with pandas, scipy, statsmodels, scikit-learn, matplotlib/seaborn; existing R/phytools Pagel workflow; pytest; optional FreeSASA/DSSP and FoldX only after their inputs are staged and recorded.

**Spec:** The additional scientific questions identified in the current project discussion, interpreted against the verified pipeline state in `/Users/ad2347/Documents/Obsidian_Wiki/projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md`.

## Global Constraints

- Read the pipeline verification ledger before trusting any result; use later verified ledger rows over stale prose in older status and figure-plan documents.
- Use `timing_annotations_v2.csv` for all timing and origin claims; the older timing file is not a source of truth for those claims.
- Treat raw and acquired data as immutable; write only derived tables to `results/analysis/` and figures to a separate generated-output directory.
- Keep Logseq DB strictly read-only.
- Record new findings, decisions, and audit results in the Obsidian project wiki, especially the verification ledger and project status/tasks pages.
- Pre-specify primary contrasts, denominators, matching variables, and multiple-testing correction before inspecting final p-values.
- Do not interpret a missing test as a biological null when alignment depth, structural coverage, or variant ascertainment is inadequate.

## Current Evidence Baseline

The implementation should begin from the current verified state rather than historical values:

- 100 canonical OXPHOS genes: 13 mtDNA and 87 nuclear-encoded.
- 6,809 classified variants: 6,340 clean and 469 warnings.
- 2,213 classified AA-level cDAVs and 1,894 classified NT-level cDAVs.
- 41,373 DAR-contact pairs in the compartment universe; 9,375 tested DAR-contact pairs in the structural/phylogenetic analyses.
- Corrected Pagel regeneration: 47 significant pairs at global FDR ≤ 0.10.
- Current Pyvolve output: 314/912 pairs significant at permutation p < 0.05.
- Current extended timing output: 98.3% permissive for mtDNA-DAR directions and 94.3% for nuclear-DAR directions.
- Current Tier-1 extraction: 86 rows, comprising 14 Tier-1A, 33 Tier-1B, and 39 Tier-1C pairs.
- FoldX has not yet been run on the corrected 47-pair set.

## Planned File Structure

The following files are the intended implementation surface; they are not created by this planning task.

- Create `src/analysis/00_build_analysis_contract.py`: assemble the pair-level analysis table, provenance manifest, denominators, and filter audit.
- Create `src/analysis/01_null_and_bias_controls.py`: matched nulls, label permutations, sensitivity analyses, and corrected effect-size tables.
- Create `src/analysis/02_codon_and_structural_mechanisms.py`: AA-versus-NT cDAV analysis and structural-location/physicochemical mechanism tests.
- Create `src/analysis/03_origins_and_networks.py`: recurrent origins, timing summaries, mito-nuclear directionality, and higher-order contact networks.
- Create `src/analysis/04_predictive_validation.py`: leakage-controlled candidate scoring and held-out validation when an independent outcome is available.
- Create `src/analysis/05_make_manuscript_figures.py`: deterministic generation of the main and supplementary figures from derived tables only.
- Create `tests/test_analysis_contract.py` and `tests/test_analysis_invariants.py`: schema, denominator, no-leakage, and statistical-output invariants.
- Create derived outputs under `results/analysis/` and rendered figures under `reports/figures/`; never write derived results into raw-data directories.

## Research Sequence

### Task 1: Freeze the analysis contract

**Files:**
- Create: `src/analysis/00_build_analysis_contract.py`
- Create: `tests/test_analysis_contract.py`
- Create: `results/analysis/analysis_contract.json`
- Create: `results/analysis/dar_contact_analysis_table.parquet`

**Purpose:** Establish one auditable row definition before adding new tests.

- [ ] Define the primary unit as one DAR position paired with one structural contact position, with stable keys containing gene, residue, partner gene, partner residue, structure/model, genome compartment, and complex.
- [ ] Join only current inputs: classified variant parquet, canonical structural contacts, `tier1_pairs.csv`, `conditional_permissiveness.csv`, `timing_annotations_v2.csv`, `timing_mitonuclear_extended.csv`, `pagel_by_complex.csv`, `pyvolve_rate_controlled.csv`, and current DCA/MI outputs.
- [ ] Record cDAV definition separately for AA-level and NT-level status; never collapse them into one label.
- [ ] Add alignment depth, species count, structure coverage, contact distance, contact class, intraprotein/interprotein status, mito-nuclear status, SDH status, and all analysis-specific missingness indicators.
- [ ] Store source file paths, source modification times, row counts, and filter expressions in `analysis_contract.json`.
- [ ] Test unique pair keys, expected current denominators, v2 timing provenance, absence of raw-data writes, and one-to-many joins that could duplicate pairs.

**Gate:** No downstream analysis proceeds until the contract reproduces the current 912 Pyvolve rows, 47 current Pagel-significant pairs, and the current Tier-1 counts without manual editing.

### Task 2: Quantify and control ascertainment and evolutionary-rate bias

**Files:**
- Create: `src/analysis/01_null_and_bias_controls.py`
- Create: `results/analysis/null_model_effects.csv`
- Create: `results/analysis/sensitivity_analysis.csv`
- Modify: `tests/test_analysis_invariants.py`

**Purpose:** Answer whether the apparent cDAV/contact signal exceeds what can be explained by conservation, mutation opportunity, alignment quality, cDAV prevalence, or disease ascertainment.

- [ ] Define primary outcomes before testing: Pyvolve significance, log odds ratio, Pagel FDR, and timing class.
- [ ] Construct matched controls within gene, genome compartment, complex, alignment-depth bin, `n_cdav_spp` bin, and residue conservation bin wherever the data support matching.
- [ ] Run cDAV/uDAV label permutations within gene and analysis stratum; preserve the observed number of cDAV positions and contact topology.
- [ ] Run contact versus non-contact comparisons using the same residue and gene eligibility rules.
- [ ] Repeat the primary analyses after SDH exclusion, leave-one-complex-out analysis, and exclusion of warning-classified variants.
- [ ] Report effect sizes and confidence intervals, not only p-values; apply Benjamini–Hochberg correction separately for each pre-specified test family.
- [ ] Include a negative-control analysis using pairs that cannot be structurally contacted under the active model policy.

**Gate:** A result is called robust only if its direction and practical effect remain stable across the matched null, SDH exclusion, warning exclusion, and leave-one-complex-out analyses. Otherwise classify it as exploratory.

### Task 3: Test codon accessibility versus amino-acid permissiveness

**Files:**
- Create: `src/analysis/02_codon_and_structural_mechanisms.py`
- Create: `results/analysis/aa_vs_nt_cdar.csv`
- Create: `results/analysis/codon_accessibility_model.csv`
- Create: `tests/test_analysis_invariants.py`

**Purpose:** Determine whether cDAV observations primarily reflect convergence on a tolerated amino acid or reuse of the exact nucleotide codon.

- [ ] Compare AA-level and NT-level cDAV rates by genome, complex, gene, variant class, codon degeneracy, and nucleotide distance from the human allele.
- [ ] Model NT-level cDAV conditional on AA-level cDAV, codon degeneracy, transition/transversion opportunity, GC content, alignment depth, and species sampling.
- [ ] Use gene-stratified or mixed-effects models so highly represented genes do not dominate the comparison.
- [ ] Report the number of variants for which NT-level classification is unavailable separately from genuine NT-level non-cDAV outcomes.
- [ ] Interpret exact-codon recurrence as evidence about mutational accessibility and historical reuse, not automatically as stronger functional compensation.

**Gate:** The manuscript must state whether NT-level cDAV is a mechanistic subset of AA-level cDAV or a distinct signal after codon opportunity is controlled.

### Task 4: Localize the structural mechanism of compensation

**Files:**
- Modify: `src/structural/00_map_davs_to_structure.py` only if a missing metric is needed and after a separate audit
- Create: `results/analysis/structural_mechanism_effects.csv`
- Create: `results/analysis/structure_coverage_by_gene.csv`
- Extend: `src/analysis/02_codon_and_structural_mechanisms.py`

**Purpose:** Test whether cDAV-associated contacts occupy specific structural environments and whether rescue is chemically interpretable.

- [ ] Compare cDAV-associated and matched uDAV-associated residues by Cβ distance, contact class, intraprotein/interprotein status, interface status, and complex.
- [ ] Quantify whether significant Pagel/Pyvolve pairs are enriched at interfaces, buried positions, or high-contact-degree residues.
- [ ] Use the current hbond, electrostatic, hydrophobic, and van der Waals labels to test physicochemical rescue categories.
- [ ] Add solvent accessibility only as a derived metric from the active structure models, with model ID and calculation method recorded per residue; keep it supplementary if coverage is incomplete.
- [ ] Separate direct residue contacts from proxy mappings and report structure-panel coverage explicitly.
- [ ] Test whether distance/contact-class effects survive matching on conservation and alignment quality.

**Gate:** Do not claim a preferred chemical mechanism unless the effect survives contact-distance, gene, and conservation matching and is reproducible across at least two complexes or clearly labelled as a case study.

### Task 5: Resolve evolutionary timing, recurrence, and directionality

**Files:**
- Create: `src/analysis/03_origins_and_networks.py`
- Create: `results/analysis/origin_recurrence.csv`
- Create: `results/analysis/timing_directionality.csv`
- Create: `tests/test_analysis_invariants.py`

**Purpose:** Extend the existing pre-adaptation result into a formal comparison of constitutive permissiveness, reactive rescue, recurrent origins, and mito-nuclear directionality.

- [ ] Use `timing_annotations_v2.csv` and `timing_mitonuclear_extended.csv` exclusively for timing claims.
- [ ] Summarize timing by genome direction, complex, SDH/non-SDH status, and intraprotein/interprotein class.
- [ ] Compare contact-first, permissive-background, constitutively-permissive, co-adaptive, and contact-after events using branch-event denominators rather than pair counts alone.
- [ ] Quantify independent origins per pair and test whether recurrent origins exceed a gene-specific null using `multi_origin_binomial.csv` or its current reproducible equivalent.
- [ ] Test mito-nuclear symmetry with matched Pyvolve outcomes and confidence intervals; do not use low-powered Fisher asymmetry as the primary directionality test.
- [ ] Perform a sensitivity analysis excluding residues with near-invariant contact states, because extreme contact-first rates can reflect structural conservation rather than active tracking.

**Gate:** The primary conclusion must distinguish “pre-existing permissiveness” from “active compensatory rescue”; these are separate biological mechanisms.

### Task 6: Build and test higher-order compensatory networks

**Files:**
- Extend: `src/analysis/03_origins_and_networks.py`
- Create: `results/analysis/compensatory_network_nodes.csv`
- Create: `results/analysis/compensatory_network_edges.csv`
- Create: `results/analysis/network_permutation_tests.csv`

**Purpose:** Determine whether pairwise signals form residue-level modules or are isolated events.

- [ ] Build graphs from active structural contacts, retaining edge type, distance, genome compartment, and evidence source.
- [ ] Overlay current Pagel-significant, Pyvolve-significant, Tier-1, and timing-supported edges without treating different evidence types as interchangeable.
- [ ] Identify connected components, residue hubs, cross-subunit bridges, and three-to-five-residue motifs.
- [ ] Compare observed module size and centrality against gene- and degree-preserving permutations.
- [ ] Require replication across independent evidence layers before calling a module biologically supported.
- [ ] Use network results for hypothesis generation and supplementary figures unless independent validation shows predictive value.

**Gate:** Do not infer epistasis solely from graph connectivity; network claims require a null model and, for strong claims, experimental or independent evolutionary support.

### Task 7: Evaluate predictive and translational utility

**Files:**
- Create: `src/analysis/04_predictive_validation.py`
- Create: `results/analysis/candidate_scores_current.csv`
- Create: `results/analysis/held_out_validation.csv`
- Extend: `src/mutagenesis/00_prioritize_pairs.py` only after confirming the current 47-pair merge inputs

**Purpose:** Test whether the framework can prioritize useful experimental targets or improve interpretation of untested variants.

- [ ] Define the prediction target before model fitting: current high-confidence pair status, independent functional rescue, or an external clinical outcome. Do not use a post hoc composite label as its own validation target.
- [ ] Use gene- or complex-held-out splits to prevent closely related pairs from appearing in both training and test sets.
- [ ] Compare simple baselines against models using cDAV status, timing, Pagel/Pyvolve evidence, contact geometry, physicochemical class, and alignment depth.
- [ ] Report calibration, precision among top-ranked candidates, and confidence intervals; do not report only ranking scores.
- [ ] Stage FoldX for the corrected current 47-pair set, record the exact structure/model and FoldX version, and keep the historical ~643/693 plan separate.
- [ ] Treat ESM2 and clinical phenotype integration as extensions requiring explicit data and validation, not as current conclusions.

**Gate:** A translational claim requires held-out or independent validation. Without it, report a ranked experimental hypothesis set only.

### Task 8: Generate the manuscript evidence package

**Files:**
- Create: `src/analysis/05_make_manuscript_figures.py`
- Create: `reports/figures/figure_1_cdar_landscape.pdf`
- Create: `reports/figures/figure_2_alignment_and_dca.pdf`
- Create: `reports/figures/figure_3_pyvolve_and_nulls.pdf`
- Create: `reports/figures/figure_4_pagel_and_complexes.pdf`
- Create: `reports/figures/figure_5_timing_and_mito_nuclear.pdf`
- Create: `reports/figures/figure_6_targets_and_mechanism.pdf`
- Create: `results/analysis/figure_manifest.tsv`

**Main-figure narrative:**

- [ ] **Figure 1:** establish scale and cDAV prevalence: AA/NT cDAV rates by genome, complex-level distribution, and gene-level coverage.
- [ ] **Figure 2:** show the measurement boundary: mtDNA DCA enrichment alongside nuclear effective-sequence-depth stratification. The nuclear global DCA result must be labelled as power-limited, not biological absence.
- [ ] **Figure 3:** show Pyvolve enrichment, effect-size distributions, rate-controlled mtDNA versus nuclear comparison, and matched-null sensitivity.
- [ ] **Figure 4:** show corrected Pagel results across complexes, SDH versus non-SDH stratification, and the distribution of evidence types. Move a full contact-class heatmap to supplement unless it has a significant, interpretable result.
- [ ] **Figure 5:** show timing classes, branch-event denominators, and bidirectional mito-nuclear permissiveness. Make the constitutive-permissiveness versus reactive-rescue distinction visually explicit.
- [ ] **Figure 6:** show the current Tier-1 evidence matrix, a ranked candidate plot, and one or two structural snapshots. Do not show FoldX values until the corrected 47-pair analysis has been run.
- [ ] Put full DCA distributions, Meff distributions, all Pagel attempts including NA results, complete Tier-1B/C tables, network permutations, and sensitivity analyses in supplementary figures/tables.
- [ ] Generate every panel from named columns in derived tables; prohibit hand-entered counts in plotting code.
- [ ] Write `figure_manifest.tsv` with panel ID, input table, row filter, statistical test, multiple-testing family, and headline caption claim.

**Gate:** A figure is manuscript-ready only when its numbers can be regenerated in one Docker command and its caption can identify the exact denominator and statistical test.

### Task 9: Validate, document, and promote conclusions

**Files:**
- Modify: `data/reference/pipeline_audit_registry.tsv`
- Modify: `/Users/ad2347/Documents/Obsidian_Wiki/projects/OXPHOS_DAV/concepts/pipeline-verification-ledger.md`
- Modify: `/Users/ad2347/Documents/Obsidian_Wiki/projects/OXPHOS_DAV/status.md`
- Modify: `/Users/ad2347/Documents/Obsidian_Wiki/projects/OXPHOS_DAV/tasks.md`
- Modify: `/Users/ad2347/Documents/Obsidian_Wiki/projects/OXPHOS_DAV/tasks-done.md`

- [ ] Run active Docker tests and bare pytest after the current pytest-discovery configuration is verified.
- [ ] Extend the freeze registry to structural, phylogenetic, and mutagenesis outputs, then rerun the freeze checker.
- [ ] Record each new result, audit, or changed number in the wiki ledger before updating manuscript prose.
- [ ] Mark stale historical values in the project status and figure-plan references so old values such as 44 Pagel pairs or 310 Pyvolve pairs cannot be mistaken for current results.
- [ ] Record unresolved limitations: uncalibrated 0.30 chain-assignment threshold, incomplete raw-download re-verification, duplicate Pagel chunk pathways, unrun FoldX, deferred ESM2, and the status of the phylo-v2 worktree.

## Recommended Execution Order

1. Tasks 1–2: freeze the pair-level contract and establish robust nulls.
2. Tasks 3–5: answer codon, structural, complex, and timing questions.
3. Task 8: generate a provisional manuscript figure set from the robust results.
4. Task 6: add network analyses only if the pairwise evidence supports coherent modules.
5. Task 7: run FoldX and predictive validation for the final experimental target set.
6. Task 9: perform final reproducibility and wiki-ledger promotion.

## Expected Scientific Decision Points

- If cDAV/contact effects survive matched nulls, the framework supports genuine evolutionary permissiveness rather than simple conservation bias.
- If NT-level enrichment disappears after codon-opportunity control, the dominant signal is amino-acid permissiveness rather than exact-codon reuse.
- If structural effects remain after conservation matching, specific physical environments are implicated; otherwise the mechanism should be described as broad structural constraint.
- If contact-first events dominate but contact states are nearly invariant, describe the result as constitutive compatibility, not active anticipation or rescue.
- If network modules exceed degree-preserving nulls and replicate across evidence layers, elevate them to mechanistic hypotheses; otherwise retain pairwise claims.
- If held-out prediction is not better than simple baselines, use the framework for candidate prioritization rather than clinical prediction.

