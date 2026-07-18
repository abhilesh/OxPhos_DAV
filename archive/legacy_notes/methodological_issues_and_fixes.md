# Methodological Issues and Planned Fixes

Three major methodological issues have been identified in the OxPhos cDAV analysis, plus five additional blocking issues required for publication readiness. This document describes each issue precisely, explains why it matters, details the correct fix, and is explicit about which fixes fully resolve the problem versus which remain sensitivity analyses.

---

## Issue 1: Circular Logic in the OR = 24.88 for Derived cDAV Contact-First Enrichment

### What the analysis currently claims

The `contact_first` enrichment test reports: "Derived cDAVs: OR = 24.88 (95% CI 11.35–54.51), p = 8.8 × 10⁻⁴¹. 96.7% of significant pairs had contact_first timing." This is presented as evidence for epistatic pre-adaptation — that compensating contacts structurally pre-adapt the environment before the pathogenic residue appears.

### The circular logic

For a **derived** cDAV on branch X (appeared on a specific recent branch):

1. `contact_first = True` means the contact alt-AA was at the **parent node** of branch X before the cDAV appeared.
2. All descendants of branch X inherit the contact alt-AA → `n_cdav_with_alt ≈ n_cdav_spp`.
3. This inheritance drives **both** significance tests simultaneously:
   - `branch_cooccur_fdr`: the cDAV clade and contact clade are the same clade by ancestry → trivially significant.
   - `fisher_fdr`: n_cdav_with_alt ≈ n_cdav_spp forces the Fisher OR high regardless of functional coupling.
4. The existing test ("is contact_first enriched in significant pairs?") is therefore testing a near-tautology for derived cDAVs.

**Switching from `branch_cooccur_fdr` to `fisher_fdr` as the significance criterion does NOT fix this.** Both tests are mechanically driven by the same n_cdav_with_alt ≈ n_cdav_spp inheritance pattern, just via different routes.

### Why Bernoulli sampling also fails

An initial fix proposal — drawing n_cdav_spp Bernoulli(p_contact) samples to generate the null — is also wrong. The cDAV species are a clade: they are phylogenetically correlated, so their contact states are correlated through shared ancestry. Bernoulli draws assume each cDAV species is an independent draw from the background frequency, which they are not. This makes the null OR distribution too narrow (underestimates clade-level variance), giving anti-conservative p-values in dense clades.

### The correct fix: tree-aware simulation via pyvolve

The null must preserve phylogenetic non-independence. For each derived cDAV on branch X with contact_first = True:

1. Load the per-gene Newick tree (`data/phylo/iqtree_jobs/{gene}/{gene}_tree.nwk`).
2. Fit an empirical substitution model for the contact column: WAG rate matrix with empirical amino acid frequencies as stationary distribution.
3. Simulate null contact evolution on the full tree for 10,000 replicates using `pyvolve.Evolver` — no coupling to the cDAV character.
4. For each simulation, count how many cDAV species (descendants of branch X) carry the contact alt-AA. Compute the simulated Fisher OR.
5. Empirical p-value: fraction of simulations where simulated OR ≥ observed OR.

**p_contact must use non-cDAV species only**: the background frequency should be estimated from `n_bg_with_alt / n_bg_spp` (non-cDAV species), not from the pooled estimate. Pooling includes the cDAV species being tested and biases p_contact upward when the cDAV clade is already contact-enriched, deflating the null toward "no enrichment."

This null correctly reflects the expected contact state enrichment in a cDAV clade due to shared ancestry alone. If the real OR falls within the null distribution, the signal is attributable to inheritance. If it consistently exceeds the null, there is genuine enrichment beyond what clade structure predicts.

### Multi-origin test — properly defined baseline

For DARs with `n_dar_gain_branches ≥ 2`: compare the fraction of cDAV gain events on contact+ branches to the genome-wide fraction of branches with the contact alt-AA at the parent node, **weighted by branch length**. Use `ancestral_state_maps.json` for the per-node contact state and per-gene tree branch lengths. A binomial test per DAR gives a principled p-value. This subset excludes single-origin DARs (which drive the OR = 24.88), so the result may differ substantially and should be reported separately.

### Sensitivity analysis (secondary, not primary)

Re-run the contact_first enrichment test using only `fisher_fdr` as the significance criterion (excluding `branch_cooccur_fdr`). This reduces but does not eliminate the circular logic — the Fisher OR is still mechanically inflated by inheritance. Report as a sensitivity check labelled explicitly as partially confounded.

### Pre-registered decision thresholds

Before running the tree-aware test, record these interpretation criteria:
- OR_conditional > 5 with majority perm_p < 0.05 → pre-adaptation signal is genuine; keep headline claim
- 2 < OR_conditional ≤ 5 → real but weaker; qualify abstract
- OR_conditional ≤ 2 or majority perm_p > 0.05 → retract pre-adaptation claim from abstract

### Reporting expectation

The OR will fall substantially — likely from ~25 to 2–5. If it holds above 5 under tree-aware simulation, the result is genuinely strong. If it falls below 2, the current abstract framing of epistatic pre-adaptation is unsupported.

### Implementation

New script: `src/phylo/04_conditional_permissiveness.py`

New outputs:
- `results/phylo/conditional_permissiveness.csv` — per-pair: `perm_p`, `null_median_or`, `observed_or`, `contact_first`, `n_dar_gain_branches`, ASR confidence stratum
- `results/phylo/multi_origin_binomial.csv` — per-DAR binomial test with branch-length-weighted expected rate
- `results/phylo/contact_first_revised_test.csv` — fisher-only sensitivity analysis (labelled as such)

---

## Issue 2: Phylogenetically-Confounded MI — Why the Standard Permutation Fails and What to Do Instead

### What the analysis currently claims

MI_APC scores with `mi_percentile > 75` are used as supporting evidence for functional coupling in the composite mutagenesis score.

### The problem

APC corrects for entropic background MI but not for phylogenetic non-independence. In a 289-species mammalian alignment, two positions both conserved within a single clade (e.g., a primate-specific substitution in protein A and an unrelated primate-specific substitution in protein B) will have elevated MI_APC. This is shared ancestry, not functional coupling.

### Why random column permutation is wrong

Shuffling column j randomly across species destroys its phylogenetic structure, making the null MI distribution unrealistically low. The observed MI (which carries phylogenetic inflation from both columns) almost always beats this artificially low null — generating p-values that appear significant but mean nothing about functional coupling.

Concrete example: two primate-specific substitutions in unrelated OxPhos proteins will have high MI_APC (both fixed in primates, absent elsewhere). After shuffling column j, primate-specificity is destroyed → null MI ≈ 0 → observed MI appears hugely significant. The resulting p-value is meaningless.

This is worse than no correction because it produces numbers that look more rigorous without addressing the underlying problem.

### The correct fixes

**Primary: plmDCA via pydca (mandatory — commit to plmDCA, not MF-DCA)**

Use `pydca.plmdca.PlmDCA`, not `pydca.meanfield_dca.MeanFieldDCA`. Mean-field DCA is the older approximation; plmDCA has been the field standard for ~10 years and gives substantially better discrimination between direct and indirect couplings. Two key mechanisms handle phylogenetic bias:

1. **Sequence reweighting** (seqid=0.8): sequences with >80% pairwise identity are downweighted, reducing over-representation of dense clades.
2. **Regularized graphical model**: finds direct couplings while marginalizing over indirect couplings through third positions, which includes many phylogenetically-driven transitive correlations.

**Interprotein DCA requires careful verification** (non-trivial):
- Species row pairing must be verified: row N in dar_gene MSA must correspond to row N in contact_gene MSA from the same species.
- M_eff check required: plmDCA requires M_eff > L (effective sequences > concatenated alignment length) for reliable results. With ~289 species and seqid=0.8 reweighting, M_eff is typically 50–150 — potentially below the concatenated length for longer protein pairs. Flag pairs where M_eff < L as unreliable.
- mt-nuc rate asymmetry: mtDNA evolves ~10× faster. This asymmetry can bias DCA for the 38 mt-nuc inter-protein pairs. Flag these for additional scrutiny.

**Cross-check: EVcouplings for top 50 pairs**

For the top 50 pairs, compute DCA scores independently using EVcouplings. Pairs scoring highly in both implementations are the strongest candidates. Pairs where the two implementations disagree are flagged as uncertain. Report the agreement statistics.

**Secondary: tree-aware null via pyvolve (top 50 pairs)**

For the top 50 pairs by plmDCA score, validate with a proper tree-aware null. Simulate null evolution of column j on the per-gene tree using a fitted WAG + empirical frequency model — NOT raw frequency tables alone. A proper substitution model (rate matrix + rate heterogeneity) is required to produce a biologically realistic null; using only empirical frequencies produces incorrect substitution dynamics.

```python
import pyvolve
tree = pyvolve.read_tree(file=str(tree_path))
# Empirical AA frequencies in pyvolve ordering
freqs = compute_empirical_aa_freqs(col_j_data)  # 20-element vector
model = pyvolve.Model("WAG", {"state_freqs": freqs})  # WAG matrix + custom stationary
partition = pyvolve.Partition(models=model, size=1)
for _ in range(1000):
    evolver = pyvolve.Evolver(tree=tree, partitions=partition)
    evolver()
    sim_col_j = extract_column_from_simulation(evolver)
    null_mi_vals.append(compute_mi_apc(col_i, sim_col_j, bg_stats))
tree_pvalue = np.mean(np.array(null_mi_vals) >= observed_mi_apc)
```

Empirical p-value: fraction of simulations with MI_APC ≥ observed. This null preserves clade structure, so inflation due to shared ancestry appears in both null and observed distributions.

### Reporting expectation

After plmDCA, expect the number of pairs with significant coevolution evidence to fall substantially from the 176/693 currently meeting `mi_percentile > 75`. The DCA-significant pairs with M_eff adequacy represent the genuinely robust signals. MI_APC scores should be retained in output for transparency but not used as a primary criterion.

### Implementation

Updated script: `src/mutagenesis/02_mi_analysis.py`

New output columns in `results/mutagenesis/mi_scores.csv`:
- `dca_di` — plmDCA direct information score
- `dca_di_percentile` — rank among all column pairs in that alignment
- `dca_meff` — effective sequence count after reweighting
- `dca_note` — flags: `insufficient_meff`, `mt_nuc_rate_asymmetry_warning`, `evcouplings_disagrees`
- `dca_tree_pvalue` — pyvolve null p-value (top 50 only; NaN otherwise)
- `evcouplings_di` — EVcouplings DI score (top 50 only; NaN otherwise)

Dockerfile additions: `pip install pydca pyvolve evcouplings`

---

## Issue 3: FoldX Threshold at Noise Floor + Epistasis SD

### What the analysis currently claims

"9 of 50 pairs satisfy both criteria for genuine compensatory rescue (ΔΔG_DAR > 0.5 AND ΔΔG_rescue > 0.5 kcal/mol). The NDUFA9 M220I/S216Y pair exhibits the strongest epistasis signal (ΔΔG_epistasis = −1.22 kcal/mol)."

### The problems

**1. Threshold at noise floor**: FoldX RMSE vs experiment is ~1.25 kcal/mol. Run-to-run SD for standard soluble proteins is 0.2–0.6 kcal/mol; higher for membrane proteins and cryo-EM structures (the majority of our pairs). The 0.5 kcal/mol threshold cannot reliably distinguish real compensation from FoldX stochastic variation.

**2. SD not being parsed**: Column 1 of `Average_*_Repair.fxout` contains the SD across 3 replicates. The current parser reads only column 2 (mean ΔΔG) and discards the SD — we have the data but are not using it.

**3. Compounded epistasis variance**:
```
ΔΔG_epistasis = ΔΔG_double − (ΔΔG_DAR + ΔΔG_contact)
SD_epistasis  = sqrt(SD_DAR² + SD_contact² + SD_double²)
```
With typical per-condition SDs of 0.2–0.5 kcal/mol, SD_epistasis reaches 0.4–0.9 kcal/mol. The NDUFA9 M220I/S216Y epistasis of −1.22 kcal/mol may be only 1.5 SDs from zero — not the 2+ SDs its magnitude implies. The "synergistic stabilization" claim requires explicit verification.

### The fix

**1. Parse SD from existing fxout files** (no FoldX re-run needed): Update `parse_buildmodel_ddg()` to return `(ddg_mean, ddg_sd)` from columns 2 and 1 respectively. Add `--reparse` mode to walk `results/mutagenesis/pdb_converted/` and re-extract without new FoldX runs.

**2. Compute propagated epistasis SD**: `ddg_epistasis_sd = sqrt(SD_DAR² + SD_contact² + SD_double²)`

**3. Raise threshold to 1.0 kcal/mol**: Re-classify pairs:
- `foldx_tier = "strong"` — rescue > 1.0 AND DAR destabilization > 1.0 kcal/mol
- `foldx_tier = "borderline"` — rescue 0.5–1.0 (within FoldX error; suggestive, not confirmed)
- `foldx_tier = "not_supported"` — rescue < 0.5

**4. Epistasis gate**: require `|ΔΔG_epistasis| > 2 × ddg_epistasis_sd` for non-additive epistasis claims. If NDUFA9 fails this gate, remove it from the abstract. If it passes, it is the lead case study.

**5. Bootstrap CI (preferred over 2-SD threshold)**: With only 3 replicates, normality is questionable. If individual replicate ΔΔG files (`Dif_*_Repaired_1.fxout` etc.) are present in `pdb_converted/`, bootstrap the epistasis estimate by resampling replicates. More reliable than assuming Gaussian error with 2 degrees of freedom.

**6. More replicates for top case studies**: For the 5–10 pairs featured most prominently, re-run FoldX with `--numberOfRuns=10`. This tightens CIs substantially and resolves whether borderline pairs (rescue 0.5–1.0) are genuine or noise. Cost: ~1 hour total compute for 10 pairs.

### This fix is complete

Unlike Issues 1 and 2, Issue 3 requires no new methodology — just parsing existing output, propagating errors, and raising the threshold. All the data needed already exist in the fxout files.

### Reporting

Report both thresholds: "N pairs at 1.0 kcal/mol threshold; 9 pairs at 0.5 kcal/mol threshold (within FoldX error, shown for completeness)." For borderline pairs, use language: "borderline suggestive, within FoldX run-to-run error." Do not describe them as "rescue-confirmed."

---

## Five Additional Blocking Issues

### A. Tier Stratification (blocking for publication)

All 2,282 cDAVs are currently `Tier = Unassigned`. The dataset mixes confirmed pathogenic variants with VUSes of unknown significance. Tier stratification using MITOMAP status, ClinVar star rating, gnomAD AF, AlphaMissense, and dbNSFP is required before any clinical or mechanistic claim can be made. All key results must be reported stratified by tier.

Suggested tiers:
- **Tier 1**: MITOMAP Confirmed OR ClinVar ★★★★ pathogenic
- **Tier 2**: MITOMAP Reported OR ClinVar ★★★ pathogenic/likely pathogenic
- **Tier 3**: ClinVar ★★ with functional evidence OR gnomAD AF < 0.0001 AND AlphaMissense ≥ 0.8
- **Tier 4**: VUS/unclassified

### B. SDH Stratification (blocking for cross-genome comparison)

Complex II (SDH, ~39.7% of nucDNA cDAVs) causes disease through succinate accumulation, not OXPHOS dysfunction. The "nuclear > mitochondrial compensation" finding (10.8% vs 6.1%) is likely heavily driven by SDH. Re-run `03_comparative_analysis.py` with SDH excluded. Report with and without SDH side-by-side. If the difference disappears after exclusion, revise or retract the cross-genome claim.

### C. Matched Null for Compensation Rate Baseline

The 8.6% significance rate (693/8,087) has no comparison baseline. Generate N=5 matched control pairs per real pair from the same PDB structure, matched on: same complex, same contact class, same Cβ-Cβ distance bin, same secondary structure type. Run the same statistical tests on controls. Report the matched-null significance rate alongside 8.6% — the difference is the actual enrichment due to disease-variant co-evolution, not structural proximity.

### D. ASR Confidence Stratification

The `contact_first` timing depends on parent-node ASR posteriors in `ancestral_state_maps.json`. Apply contact_first analyses separately for:
- High-confidence ASR: posterior probability ≥ 0.80 for the contact alt-AA at the parent node
- Lower-confidence ASR: posterior probability 0.5–0.80

If the OR is higher in the high-confidence stratum and lower in the uncertain stratum, that pattern validates the contact_first signal. If they are similar, the signal may be driven by ASR uncertainty.

### E. Dataset Asymmetry Test for mt vs nuc Comparison

mtDNA (443 variants, heavily curated) vs nucDNA (6,554 variants, mixed quality) cannot be directly compared for compensation rates. Extract a matched-stringency nucDNA subset (ClinVar ★★★+ pathogenic, gnomAD AF < 0.001). Compare compensation rates between this matched subset and mtDNA. Report both the full comparison and the matched-stringency comparison. If the difference disappears in the matched subset, the apparent genomic difference reflects curation bias.

---

## Disclosure Standards

For any manuscript, the following must be stated explicitly:

**Issue 1**: Report the original OR = 24.88 with a note that it reflects clade-structure inheritance. Report the tree-aware conditional permutation result as the primary result. State the pre-registered interpretation thresholds. The fisher-only sensitivity analysis should appear as supplementary material.

**Issue 2**: State that APC-MI was used for initial pair identification, but plmDCA is the primary reported co-evolution metric. Note that interprotein pairs with M_eff < L are flagged as having unreliable DCA estimates. Retain MI_APC in supplementary data.

**Issue 3**: Report both thresholds (0.5 and 1.0 kcal/mol) with pair counts. Report epistasis SD explicitly. State whether bootstrap CI or propagated SD was used for epistasis confidence. Report re-run results for top 10 case studies at N=10 replicates.

**Partial fixes disclosed as partial are more credible than partial fixes presented as complete.** The conditional permutation and plmDCA approaches are genuinely more rigorous; they should be the primary analyses, with original results as historical context.

---

## Run Order

1. **Issue 3 first** (FoldX reparse): minutes; tells you immediately which pairs survive 1.0 kcal/mol and whether NDUFA9 survives 2-SD before committing to expensive validation.
2. **Issue 2 in parallel** (plmDCA overnight) + **Issues A, B, D, E** (SDH filter, tier join, ASR stratification, dataset asymmetry): run in parallel with Issue 3.
3. **Issue 1** (tree-aware permutation): most computationally intensive; run after Issues 2 and 3 clarify which pairs to focus on.
4. **Step 4**: re-compile targets with revised criteria after all fixes complete.
5. **Step 5**: one revision pass of all notes/reports.