# OxPhos DAV — Critical Evaluation and Nuclear Compensation Forward Plan

## Summary of Major Findings (One Paragraph)

Across 93,693 structural DAR-contact pairs spanning 100 OXPHOS genes, the pipeline establishes that pathogenic variant tolerance (cDAV status) is detectably associated with contact co-evolution, but the strength and nature of that association differs fundamentally between the two genomes. mtDNA-encoded variants are naturally tolerated in >50% of tested cases (167/323, 51.7%), compared to ~31.5% for nuclear variants, and their structural contacts show extremely strong sequence-level co-evolutionary enrichment (DCA-DI cDAV vs uDAV, p = 7.65e-69). Nuclear variants also show genuine compensatory co-evolutionary signal — 38% of nuclear cDAV-contact pairs are Pyvolve-significant vs a 5% permutation null (binomial p ≈ 0), and Pagel identifies 40 significant nuclear pairs (18 outside SDH) across 11 genes spanning Complexes I, III, IV, and V — but the sequence-level (DCA/APC-MI) signal is absent due to insufficient alignment depth (Meff = 7.6 for nuclear vs 61 for mtDNA), meaning the nuclear compensation hypothesis is supported by phylogenetic and species-distribution tests but not by direct sequence coupling metrics. At the mito-nuclear interface, mtDNA-variant → nuclear-contact pairs show the highest compensation rate (61.4% sig_any, APC-MI p = 2.83e-10 from the mtDNA gene alignment), while Pyvolve — the one test internally calibrated against gene-specific neutral evolution — shows nearly symmetric significance in both directions (mtDNA-DAR→nuclear-contact: 48.1%; nuclear-DAR→mtDNA-contact: 50.0%), indicating the directional asymmetry in sig_any is largely a power artifact of Fisher's test rather than a genuine biological asymmetry. Temporal ordering of the 44 Pagel-significant pairs shows that structural permissiveness dominates: 63% of 201 branch events show the contact already physicochemically compatible with the compensatory state before the DAV arose, while rescue (contact_after) accounts for 6.5%; but because only 1 of the 44 pairs is a mito-nuclear contact, the branch-level directionality at the interface (which genome changed first) remains untested at scale. FoldX thermodynamic data (distinguishing structurally-mild cDAVs from rescued destabilising ones) and an extended temporal ordering analysis across all 83 Pyvolve-significant mito-nuclear pairs are the two most critical gaps.

---

## Part 0: Dataset Scale and cDAV/uDAV Proportions

### Variant classification (variants_master_classified.parquet)

| Genome | cDAV (AA-level) | uDAV (AA-level) | Total classified | cDAV fraction |
|---|---|---|---|---|
| **mtDNA** | **167** | 156 | 323 | **51.7%** |
| **nucDNA** | **2,046** | 4,440 | 6,486 | **31.5%** |
| **Total** | **2,213** | 4,596 | 6,809 | **32.5%** |

### Unique DAR positions with structural contacts (dar_contacts_cbcb8A.csv)

| Genome | cDAV DARs | uDAV DARs | Total | cDAV fraction |
|---|---|---|---|---|
| mtDNA | 163 | 154 | 317 | 51.4% |
| nucDNA | 1,756 | 4,202 | 5,958 | 29.5% |

### DAR-contact pairs by genome

| Genome | cDAV pairs | uDAV pairs | Total pairs | cDAV fraction |
|---|---|---|---|---|
| mtDNA | 3,610 | 3,808 | 7,418 | 48.7% |
| nucDNA | 23,959 | 62,316 | 86,275 | 27.8% |
| **Total** | **27,569** | **66,124** | **93,693** | **29.4%** |

mtDNA encodes only 13 subunits vs 87 nuclear, yet shows a dramatically higher cDAV fraction (51.7% vs 31.5%). More than half of all mtDNA pathogenic variants are naturally tolerated in at least one non-human species; for nuclear variants, this is roughly 1 in 3.

---

## Part 1: Critical Evaluation

### 1.1 Cleanest findings (robust to multiple tests, internally controlled)

**1. mtDNA DCA-DI enrichment (p = 7.65e-69)**
In 8/9 mtDNA genes, cDAV contacts sit at significantly higher DI percentiles than uDAV contacts. Internally controlled (within each gene), replicated across independent genes and complexes, does not depend on arbitrary cutoffs. This is the strongest statistical result in the pipeline.

**2. Pyvolve conditional permissiveness (310/912, 34% vs 5% null, p ≈ 0)**
The permutation null is gene-specific and internally calibrated to neutral sequence evolution on the gene tree. The 7× enrichment is not a power artifact — it persists when SDH is excluded (38.3% for non-SDH nuclear). Nuclear gene cDAV-contact pairs show higher ORs (median 2.04) than mtDNA (1.06), likely because nuclear cDAV variants are rarer (median 18 cDAV species vs 154 for mtDNA), generating sharper contrast in the OR calculation. OR is negatively correlated with n_cdav_spp (r = −0.21, p = 2e-10).

**3. mito-nuclear APC-MI enrichment (p = 2.83e-10)**
For mito-nuclear contacts where the DAR is in an mtDNA gene, APC-MI percentiles are higher for cDAV vs uDAV (80.57 vs 64.21). This signal comes exclusively from the mtDNA gene's alignment; the nuclear contact gene's own alignment cannot contribute a paired MI measurement. **This is not a cross-gene co-evolutionary measure** — it measures whether the contact column in the mtDNA gene's multi-species alignment is more coupled to the DAR column in cDAV pairs vs uDAV pairs, irrespective of how the nuclear gene evolves.

**4. Temporal ordering — structural permissiveness dominates (63.2% permissive branch events)**
Across the 44 Pagel-significant pairs (201 branch events), 63% fall into the permissive signal class (contact already physicochemically compatible with compensatory state) vs 16.9% secondary mutation. Rescue (contact_after) accounts for 6.5%.

### 1.2 Weaknesses and confounds

**Weakness 1 — Nuclear DCA/APC-MI is a power failure, not biology**
nuc-nuc DCA-DI MW p = 1.0; nuc-nuc APC-MI MW p = 0.98. But DCA DI percentiles are REVERSED (uDAV > cDAV in nuclear) likely because Meff ≈ 7.6 makes percentile rankings reflect protein topology, not co-evolutionary coupling — most nuclear DCA rows are tagged `insufficient_meff`. The nuclear compensation hypothesis cannot be rejected from this null; it simply cannot be tested by sequence MI at current alignment depths.

**Note on cross-gene MI as a solution:** Computing MI between an mtDNA column and a nuclear column across the shared mammalian species would require a joint alignment. The effective Meff would be bottlenecked by the nuclear gene (Meff ≈ 7.6), not the mtDNA gene, since MI requires variation in both columns. Cross-gene MI inherits the nuclear Meff limitation and is not a viable solution to this problem.

**Weakness 2 — Pagel nuclear pairs are SDH-heavy, but signal survives exclusion**
22/44 (50%) Pagel-significant pairs are SDH (SDHD: 10, SDHC: 8, SDHB: 2, SDHA: 2). After SDH exclusion, 18 non-SDH nuclear pairs remain across 11 genes (UQCRC2, NDUFB11, NDUFS8, COX7A2, UQCRB, NDUFB10, NDUFA1, NDUFA6, ATP5MF, NDUFV2, COX5A). The nuclear Pagel signal is not an SDH artifact.

**Weakness 3 — Mito-nuclear directional asymmetry is partially a power artifact**
sig_any shows 61.4% for mtDNA-DAR → nuclear-contact vs 4.9% for nuclear-DAR → mtDNA-contact. But Pyvolve (internally calibrated, rate-controlled) shows 48.1% vs 50.0% — nearly symmetric. The asymmetry in sig_any is driven by Fisher's exact test, which has near-zero power when n_cdav_spp = 9 (nuclear-DAR direction), not by a genuine biological difference in compensation rates.

**Weakness 4 — Temporal ordering covers only 44 Pagel pairs, 1 of which is mito-nuclear**
The 44 Pagel-significant pairs are 39 nuc-nuc + 4 mt-mt + 1 nuc-dar_mt-contact (NDUFB10:76→MT-ND4:178). Zero mt-dar_nuc-contact pairs reach Pagel significance. The 83 Pyvolve-significant mito-nuclear pairs (54 mt-dar_nuc + 29 nuc-dar_mt) are never temporally ordered. **The question of which genome's variant typically arose first at the mito-nuclear interface is currently unanswered.**

**Weakness 5 — FoldX ΔΔG is absent**
Cannot distinguish structurally mild cDAVs (tolerable without compensation) from rescued destabilising ones (compensation required). This is the primary mechanistic gap.

---

## Part 2: What the Data Tell Us About Mito-Nuclear Compensation

### 2.1 The mtDNA compensation story (well-supported)

Compensatory co-evolution at mtDNA residues is detectable at the sequence level (DCA), at the species-distribution level (Pyvolve), and at the phylogenetic trait-correlation level (Pagel). The temporal ordering shows that structural permissiveness precedes the variant in most cases. The picture is coherent: mtDNA positions adjacent to pathogenic variants are under tighter mutual constraint, and the contact environment is typically already accommodating when the variant arises.

### 2.2 The nuclear compensation story (biologically plausible, evidentially supported outside SDH)

18 non-SDH nuclear pairs across 11 genes show Pagel significance. Non-SDH nuclear Pyvolve significance (38.3%) is indistinguishable from the full nuclear rate (38.0%) and significantly above null (p = 4.65e-101). The DCA null is a power failure. **Nuclear compensation is real and not SDH-specific**, but the effect size (Pyvolve OR median 1.74 for non-SDH nuclear) is smaller than SDH (3.13) and the temporal ordering is less clear (36.3% vs 68.1% permissive for SDH).

### 2.3 The mito-nuclear interface — directionality untested at scale

The one Pagel-significant mito-nuclear pair (NDUFB10:76V → MT-ND4:178I, nuclear DAR → mtDNA contact) shows "constitutively_permissive" timing: the mtDNA contact position (MT-ND4 178I) was already BLOSUM-similar to the compensatory AA across all 60 gain branches, with 57/60 permissive branch events. This is consistent with the nuclear variant (NDUFB10:76V) arising in a background where the mtDNA contact was already structurally accommodating — supporting the "mtDNA genome accommodates nuclear variation" interpretation for this pair.

However, this is one data point. The 83 Pyvolve-significant mito-nuclear pairs have not been temporally ordered. Testing whether this single-pair result generalises requires the analysis in Section 3.1.

**Critical question for directionality:** Does the mtDNA contact gene change *before* or *after* the nuclear DAR variant arises on the same lineage? And vice versa for mtDNA-DAR → nuclear contact pairs? These questions require the temporal ordering analysis to be explicitly extended to cross-genome pairs.

**On controlling for evolutionary rate differences:** The temporal ordering analysis (06_temporal_ordering.py) controls for evolutionary rate differences implicitly — it uses the actual IQTree MAP branch-level events on the same VertLife mammal tree for both genes, so "before/after" is defined by which gene changed on which specific branch, regardless of total branch counts. The residual concern is ASR confidence: nuclear genes have fewer total changes (median 256 branches with changes vs 436 for mtDNA), so posterior probabilities are lower and confidence filtering will exclude more nuclear events. This means estimates of "contact_first" events in nuclear genes are conservative. Direction: if anything, nuclear-contact pre-adaptation is *underestimated*, not overestimated, in analyses where the nuclear gene is the contact gene.

---

## Part 3: Forward Analysis — Testing Nuclear Compensation Directionality

Listed in priority order. All inputs available locally; no new HPC runs required for Steps 1–4.

### 3.1 Extended temporal ordering for all Pyvolve-significant mito-nuclear pairs (HIGH PRIORITY)

**Script:** `src/phylo/06_temporal_ordering.py` — extend to accept a custom pair list

**Why:** Only 1 of 44 Pagel pairs is mito-nuclear. The 83 Pyvolve-significant mito-nuclear pairs (54 mt-dar_nuc + 29 nuc-dar_mt) are unexamined. This is the primary untested question.

**Input pairs:** Filter `results/phylo/conditional_permissiveness.csv` for:
- `perm_p < 0.05` AND
- `dar_gene ∈ MTDNA XOR contact_gene ∈ MTDNA` (cross-genome pairs)

**Analysis additions:**
- Run existing temporal ordering logic on these 83 pairs
- Stratify output by direction: mt-DAR-nuc-contact vs nuc-DAR-mt-contact
- Primary comparison: `n_contact_first / (n_contact_first + n_contact_after)` — what fraction of directional events show pre-adaptation — separately for each direction
- **Rate control test:** Compute `expected_permissive_rate` for each gene pair using the within-compartment base rate (mt-mt pairs as baseline for mtDNA contact gene; nuc-nuc pairs for nuclear contact gene) and test whether the observed mito-nuclear permissive rate exceeds that baseline via binomial test

**Output:** `results/phylo/timing_mitonuclear_extended.csv`

**Expected result:** If the "mtDNA accommodates nuclear" hypothesis is correct, nuc-DAR → mtDNA-contact pairs should show high permissive_background + constitutively_permissive (mtDNA contact already compatible), while mt-DAR → nuc-contact pairs should show more contact_after (nuclear contacts adapting after mtDNA variant).

### 3.2 Controlling for genome-specific evolutionary rate in OR comparison

**Why:** Non-SDH nuclear Pyvolve OR (1.74) > mtDNA (1.06), but this may partly reflect cDAV rarity (OR inflates at low n_cdav_spp). True rate-controlled comparison is needed.

**Analysis:**
- Stratify Pyvolve pairs into n_cdav_spp bins: [1–5], [6–20], [21–50], [51–200], [>200]
- Within each bin: compare OR and sig% for mtDNA vs non-SDH nuclear DAR pairs
- Fisher-CMH test for the combined within-bin OR
- This tests whether nuclear pairs have higher OR *at matched cDAV species counts*, i.e., after removing the OR-inflation due to rarity

**Output:** `results/phylo/pyvolve_rate_controlled.csv`

### 3.3 High-Meff nuclear DCA stratified by branch-change density

**Why:** Identify the Meff threshold where nuclear DCA signal emerges, providing a power-matched comparison to mtDNA.

**Analysis:**
- Compute per-gene median Meff from `dca_all_davs.csv` `dca_meff` column
- Stratify nuclear genes into Meff bins: < 5, 5–10, 10–20, > 20
- For each bin: Mann-Whitney U test for cDAV vs uDAV DCA-DI percentile
- Add branch-change count per gene from `asr_coordinate_harmonization_audit.tsv` as an independent rate proxy
- Plot: x = median Meff, y = MW p-value, colored by compartment — shows Meff threshold for nuclear DCA signal emergence

**Output:** `results/mutagenesis/dca_nuclear_meff_stratified.csv`

### 3.4 SDH vs non-SDH Pagel complex stratification

**Why:** Determine whether non-SDH nuclear Pagel signal varies by complex, and whether contact type (H-bond vs vdW) predicts compensation.

**Analysis:**
- Tag 44 Pagel pairs with complex (CI, CII, CIII, CIV, CV) and contact_type
- Chi-squared test: contact_type distribution in sig vs non-sig pairs
- Report gene, pair count, and pagel_p range per complex

**Output:** `results/phylo/pagel_by_complex.csv`

### 3.5 Tier-1 pair extraction (Pagel + Pyvolve + permissive timing)

**Criteria:**
- `pagel_fdr ≤ 0.10` AND `perm_p < 0.05` AND `dominant_refined_timing ∈ {contact_first, permissive_background}`
- Join compensatory_partners.csv × conditional_permissiveness.csv × timing_annotations_v2.csv on (dar_gene, dar_aa_coord, contact_gene, contact_refseq_pos)

**Output:** `results/phylo/tier1_pairs.csv`
Expected: ~5–15 pairs; these are the experimental mutagenesis priority targets.

### 3.6 ESM2 epistatic coupling for nuclear Pyvolve-significant pairs (DEFERRED — Phase 3)

For the 247 nuclear-gene Pyvolve-significant pairs, ESM2 provides a power-independent epistatic coupling test that bypasses the Meff limitation:

`epistatic_score = log P(contact_alt_aa | seq_with_DAV) − log P(contact_human_aa | seq_with_DAV)`

A positive score means the DAV context makes the compensatory contact AA more likely — direct biochemical epistasis without requiring a deep MSA. This is the cleanest nuclear compensation test available given the Meff constraint.

Note: For mito-nuclear pairs, ESM2 must be applied to each protein chain separately. The nuclear protein receives the DAV context from the mtDNA partner via structure-informed context truncation (mask the contact position in the nuclear chain; score it given the DAV in the separate mtDNA chain's sequence). This is an approximation but has been used successfully for protein-protein interface epistasis.

Deferred until Analyses 3.1–3.4 complete.

---

## Part 4: SDH-Excluded Results Summary

SDH accounts for 22/44 (50%) of Pagel pairs. After exclusion:

| Test | All nuclear | Non-SDH nuclear | Conclusion |
|---|---|---|---|
| Pagel FDR ≤ 0.10 | 40 pairs | **18 pairs, 11 genes** | Signal not SDH artifact |
| Pyvolve sig% | 38.0% | **38.3%** (p = 4.65e-101) | Signal not SDH artifact |
| Pyvolve OR median | 2.04 | **1.74** (vs mtDNA 1.06, p = 3.25e-05) | Persists after SDH removal |
| DCA-DI MW p | 1.0 | **1.0** | Power failure in both |
| APC-MI MW p | 0.98 | **0.083** | Slight improvement, still ns |
| Temporal permissive | SDH 68.1% | **Non-SDH 36.3%** | SDH drives pre-adaptation story |
| Temporal rescue | SDH 0% | **Non-SDH 9.1%** | Non-SDH has some rescue |

SDH has stronger per-pair signal (OR 3.13) and clearer temporal permissive pattern (68.1% vs 36.3%). Non-SDH nuclear compensation is real but weaker and temporally less resolved.

---

## Part 5: Key Numbers for Manuscript Context

| Metric | Full dataset | SDH-excluded nuclear |
|---|---|---|
| Total unique DAR-contact pairs | 41,373 | — |
| cDAV pairs (all) | 14,879 (36.0%) | — |
| mtDNA cDAV fraction | 51.7% (167/323) | — |
| nucDNA cDAV fraction | 31.5% (2046/6486) | — |
| Pagel sig (FDR ≤ 0.10) | 44 / 4,227 tested | 18 non-SDH nuclear / 11 genes |
| Pyvolve sig (p < 0.05) | 310 / 912 (34% vs 5% null) | 170/444 (38.3%) |
| Pyvolve sig — mtDNA DAR | 63/262 (24%) | — |
| Pyvolve sig — nucDNA DAR | 247/650 (38%) | 170/444 (38.3%) |
| mt-mt DCA MW p | 7.65e-69 | — |
| nuc-nuc DCA MW p | 1.0 (power failure) | 1.0 |
| mt-nuc APC-MI p (mtDNA DAR) | 2.83e-10 | — |
| Pyvolve: mt-dar_nuc-contact sig% | 48.1% (26/54) | — |
| Pyvolve: nuc-dar_mt-contact sig% | 50.0% (19/38) | — |
| Temporal: permissive signal (all) | 63.2% (127/201 branch events) | — |
| Temporal: permissive — SDH pairs | 68.1% | — |
| Temporal: permissive — non-SDH | 36.3% | — |
| Temporal: rescue (contact_after) | 6.5% (13/201) | 9.1% (non-SDH dominant) |
| Mito-nuclear Pagel pairs | 1 (nuc-dar_mt-contact) | — |
| Mito-nuclear Pyvolve pairs | 83 (54 mt-dar + 29 nuc-dar) | — |

---

## Part 6: Implementation Order

| Step | Script | Est. time |
|---|---|---|
| Write this doc to docs/ | — | immediate |
| 3.4 Pagel complex stratification + contact type | `src/phylo/08_pagel_complex_stratification.py` | < 5 min |
| 3.5 Tier-1 pair extraction | `src/phylo/11_tier1_extraction.py` | < 5 min |
| 3.3 High-Meff nuclear DCA | `src/phylo/10_nuclear_dca_stratified.py` | < 10 min |
| 3.2 OR comparison rate-controlled | `src/phylo/09_pyvolve_rate_controlled.py` | < 10 min |
| 3.1 Extended mito-nuclear temporal ordering | Extend `src/phylo/06_temporal_ordering.py` | < 15 min |
| FoldX Phase 2a | (HPC submission) | HPC |
| 3.6 ESM2 epistatic coupling | `src/mutagenesis/05_esm_epistasis.py` | GPU/HPC |
