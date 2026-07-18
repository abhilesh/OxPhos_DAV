# OxPhos DAV — Results Interpretation and Figure Plan

**Date:** 2026-06-04  
**Status:** Analyses 3.1–3.5 complete; FoldX HPC and ESM2 coupling deferred

---

## Overview

Across 93,693 structural DAR-contact pairs spanning 100 OXPHOS genes, the pipeline establishes that pathogenic variant tolerance (cDAV status) is detectably associated with contact co-evolution, but the strength and nature of that association differs fundamentally between the two genomes. mtDNA-encoded variants are naturally tolerated in >50% of tested cases (167/323, 51.7%), compared to ~31.5% for nuclear variants, and their structural contacts show extremely strong sequence-level co-evolutionary enrichment (DCA-DI cDAV vs uDAV, p = 7.65e-69). Nuclear variants also show genuine compensatory co-evolutionary signal — 38% of nuclear cDAV-contact pairs are Pyvolve-significant vs a 5% permutation null (binomial p ≈ 0), and Pagel identifies 44 significant pairs including 18 outside SDH across 11 genes spanning Complexes I, III, IV, and V — but the sequence-level (DCA/APC-MI) signal is absent due to insufficient alignment depth (median nuclear Meff = 7.6 vs 61 for mtDNA), meaning the nuclear compensation hypothesis is supported by phylogenetic and species-distribution tests but not by direct sequence coupling. The rate-controlled OR comparison confirms this is not an artifact of cDAV rarity: nuclear non-SDH compensatory OR remains higher than mtDNA after n_cdav_spp bin stratification (CMH χ² = 6.6, p = 0.010). At the mito-nuclear interface, Pyvolve rates are nearly symmetric in both directions (48.1% mt-dar vs 50.0% nuc-dar) — rejecting the previously apparent directional asymmetry as a Fisher's test power artifact. Critically, temporal ordering of the 45 Pyvolve-significant mito-nuclear pairs reveals overwhelming pre-adaptation in both directions: 98.3% of mt-dar branch events show the nuclear contact already in a compatible state before the mtDNA variant arose (vs 63.2% within-genome baseline, p = 5.5e-145), and 94.1% of nuc-dar branch events are similarly permissive (p = 3.7e-9). This constitutive structural permissiveness at the mito-nuclear interface, combined with the DCA Meff threshold analysis (signal emerges at Meff ≥ 10 for nuclear genes), provides a mechanistically coherent picture: the interface is under exceptionally strong purifying selection to maintain structural compatibility, and the compensatory co-evolutionary signal reflects this conservation rather than dynamic sequence-level tracking.

---

## Results Interpretations

### 2.1 Dataset Scale and cDAV Proportions

mtDNA encodes 13 of 100 OXPHOS subunits (13%) but accounts for 51.7% of pathogenic variants that are naturally tolerated in non-human species (167/323 cDAVs vs 2,046/6,486 for nuclear, 31.5%). This 1.64-fold difference in cDAV fraction likely reflects clonal mtDNA inheritance (no recombination, maternal transmission), which allows mildly deleterious mitochondrial variants to fix in some lineages through genetic drift, whereas nuclear variants face bi-parental purifying selection. At the DAR-contact pair level, 93,693 pairs are available (7,418 mtDNA, 86,275 nuclear), with the mito-nuclear interface represented by 12 pairs in the Pagel/Fisher dataset and 45 Pyvolve-significant pairs in the conditional permissiveness dataset.

**What this tells us:** The elevated mtDNA cDAV fraction is the starting observation that motivates the entire analysis. It is not simply a reflection of mtDNA evolving faster — fast evolution would increase *all* variant types equally. Instead, a specific enrichment for pathogenic amino acids being tolerated in non-human species points to systematic compensatory co-evolution.

### 2.2 mtDNA Sequence-Level Co-evolution (DCA-DI)

DCA-DI percentile is significantly higher for cDAV contacts than uDAV contacts in 8/9 mtDNA genes (global MW p = 7.65e-69). This is the strongest, most internally controlled result in the pipeline: it is replicated independently across genes, within the same statistical framework, and does not depend on arbitrary cutoffs. The one non-significant gene (MT-ND4L) has the lowest Meff among all mtDNA genes and the fewest cDAV positions, consistent with a power limitation rather than a biological exception.

**Mechanistic interpretation:** cDAV positions sit at DCA-DI hot spots — positions with high direct evolutionary coupling to their structural contacts. The contact environment of a cDAV residue carries more information about the identity of that residue than the environment of a uDAV residue. This is precisely what co-evolution would predict: when a pathogenic AA is tolerated, it is because the surrounding residues have co-evolved to accommodate it.

### 2.3 Nuclear DCA: Power Failure, Not Biological Null

Nuclear DCA-DI shows no enrichment overall (MW p = 1.0), but the Meff stratification analysis identifies this as a power failure, not a biological null. Key evidence:

- Signal is absent at Meff ≤ 2 (25 genes, p = 1.0) and Meff 2–10 (35 genes, p = 1.0)
- Signal emerges at Meff 10–20 (8 genes, MW p = 0.0088) and persists at Meff > 20 (5 genes, p = 0.0065)
- Spearman r(Meff, −log10 p) = 0.37, p = 0.0014 across 73 nuclear genes
- Two genes show individually significant signals: NDUFS6 (Meff = 26.8, p = 9.8×10⁻⁴) and COXFA4 (Meff = 11.9, p = 0.026)
- 66% of nuclear genes (48/73) have Meff ≤ 5, below the detection threshold

The absence of a global nuclear DCA signal cannot be interpreted as absence of co-evolution. Nuclear OXPHOS genes in TOGA alignments average Meff_80 ≈ 7.6 (substitution-rate-adjusted effective sequences), compared to ≈61 for mtDNA. Cross-gene MI (between mtDNA and nuclear columns in a joint alignment) would not solve this: the joint Meff is bottlenecked by the nuclear gene, inheriting the same power failure.

**What this tells us:** The nuclear DCA null is a methodological limitation, not a biological finding. The correct statement is "DCA cannot be applied to nuclear OXPHOS genes at current alignment depth." The threshold of Meff ≥ 10 provides a concrete target: OXPHOS nuclear genes with Meff above this threshold (NDUFV3, ATP5MJ, NDUFA11, NDUFS6, COX6C, ATP5IF1, UQCRB, NDUFA1, NDUFB1, COXFA4, COX6B1, NDUFS5, NDUFA3) are the priority candidates for DCA analysis with improved alignments.

### 2.4 Pyvolve Conditional Permissiveness

310/912 tested DAR-contact pairs show Pyvolve significance (p < 0.05) compared to a 5% permutation null — a 7× enrichment. The permutation null is gene-specific and calibrated to neutral sequence evolution on the actual gene tree, making this result robust to gene-specific substitution rate differences. Nuclear sig% (38%) exceeds mtDNA (24%), but this reflects lower n_cdav_spp for nuclear variants (median 18 vs 154 for mtDNA), which inflates OR when cDAV species are rare. The Pyvolve test for excess significance remains p ≈ 0 after excluding SDH (non-SDH nuclear: 38.3%, p = 4.65×10⁻¹⁰¹).

Spearman r(n_cdav_spp, log_OR) = −0.21, p = 2×10⁻¹⁰ confirms the OR-inflation mechanism. Pairs with n_cdav_spp = 1–5 (nuclear only) show 34.2% significance at median OR = 1.97; pairs with n_cdav_spp > 200 (mtDNA) show 25% significance at median OR = 1.10.

**What this tells us:** The Pyvolve test detects whether the cDAV background makes the contact more likely to carry an alternative amino acid than the neutral expectation under the gene's own substitution model. This signal is present across all contact types (nuc-nuc, mt-mt, mt-nuc) and survives SDH exclusion, establishing that compensatory co-evolution is detectable by species-distribution criteria regardless of DCA power.

### 2.5 Rate-Controlled Pyvolve OR Comparison (New)

After stratifying by n_cdav_spp bins, non-SDH nuclear pairs show higher OR than mtDNA pairs in 4/5 bins. The Cochran-Mantel-Haenszel test, which controls for n_cdav_spp stratification, confirms the nuclear advantage: χ² = 6.6, p = 0.010. Mann-Whitney test (nuclear > mtDNA log-OR) p = 1.6×10⁻⁵.

| n_cdav_spp bin | mtDNA sig% (OR) | Non-SDH nuclear sig% (OR) |
|---|---|---|
| 1–5 | no pairs | 34.2% (1.97) |
| 6–20 | 26.7% (1.06) | 41.9% (1.96) |
| 21–50 | 57.9% (1.73) | 34.7% (1.49) |
| 51–200 | 17.5% (0.94) | 31.3% (1.44) |
| >200 | 25.0% (1.10) | 52.9% (1.71) |

The 21–50 bin is the only reversal (57.9% mtDNA vs 34.7% nuclear) but contains only 19 mtDNA pairs, and the CMH test corrects for this.

**What this tells us:** The higher nuclear compensatory OR is not purely a statistical artifact of cDAV rarity — it survives matched-count comparison. This suggests a real biological difference: nuclear cDAV contacts are more conditionally associated with alternative amino acids in cDAV species than in background species, and this effect is stronger than for mtDNA contacts. One interpretation is that nuclear cDAV positions are at tighter functional bottlenecks, such that compensatory changes are more specifically required (higher OR) even though cDAV species are fewer (lower power).

### 2.6 Pagel Phylogenetic Co-evolution

Pagel's discrete test identifies 44/4,227 tested pairs (1.0%) at FDR ≤ 0.10. The distribution by complex:

| Complex | mtDNA pairs | Nuclear pairs | SDH? | Genes |
|---|---|---|---|---|
| CI | 2 | 9 | No | MT-ND2, MT-ND6, NDUFA1, NDUFA6, NDUFB10, NDUFB11, NDUFS8, NDUFV2 |
| CII | 0 | 22 | Yes | SDHA, SDHB, SDHC, SDHD |
| CIII | 0 | 3 | No | UQCRB, UQCRC2 |
| CIV | 0 | 3 | No | COX5A, COX7A2 |
| CV | 2 | 3 | No | MT-ATP6, MT-ATP6, ATP5MF |

No contact-class enrichment is detected (χ² = 2.55, df = 3, p = 0.47). H-bond contacts are over-represented among significant pairs (70.5% vs 60.0% in non-significant pairs) but not significantly (Fisher p = 0.117). The compensation mechanism is not restricted to a specific bond type.

**What this tells us:** Pagel significance reflects phylogenetically correlated binary trait evolution (the DAR and the contact co-change along lineages). Coverage across all 4 non-SDH complexes (CI, CIII, CIV, CV) establishes that the signal is a general OXPHOS phenomenon, not an artifact of any single complex's architecture. The absence of contact-class specificity is biologically reasonable: structural compensation can occur through H-bonds (e.g., restoration of a lost H-bond), hydrophobic packing, and van der Waals interactions with similar co-evolutionary signatures.

### 2.7 Complex Stratification and Non-SDH Signal (New)

The 18 non-SDH nuclear Pagel-significant pairs span CI (9 pairs), CIII (3), CIV (3), CV (3), with Pagel p ranging from 9.7×10⁻¹⁰ (ATP5MF:78/93 within CV, the most significant non-SDH result) to 4.1×10⁻³. The CV pairs (ATP5MF) are the most significant non-SDH findings and likely reflect the tight structural coupling between ATPase subunits at the rotor-stator interface.

Comparing SDH to non-SDH nuclear:
- SDH: OR median 3.13, temporal permissive 68.1% — strong structural pre-adaptation
- Non-SDH: OR median 1.74, temporal permissive 36.3%, rescue 9.1% — weaker, more dynamic compensation

**What this tells us:** SDH compensation is qualitatively different from non-SDH nuclear compensation. SDH pairs are structurally pre-wired (contact already compatible, no change needed), while non-SDH nuclear pairs show genuine contact_after events (rescue), indicating secondary site mutations that arise after the DAV. This suggests two modes of compensation: constitutive permissiveness (structural conservation maintains compatibility) and adaptive compensation (secondary mutation rescues the DAV effect).

### 2.8 Temporal Ordering — Within-Genome Pagel Pairs

Across 40 testable Pagel pairs (201 branch events), the refined timing distribution:

| Category | Branch events | % |
|---|---|---|
| contact_first | 17 | 8.5% |
| permissive_background | 55 | 27.4% |
| constitutively_permissive | 55 | 27.4% |
| co_adaptation | 34 | 16.9% |
| contact_after | 13 | 6.5% |
| no_contact_change (genuine) | 27 | 13.4% |

Permissive signal (contact_first + permissive_background + constitutively_permissive) = 63.2% of all branch events. Rescue (contact_after) = 6.5%. SDH drives the permissive story (68.1%) while non-SDH nuclear shows more genuine rescue (9.1% contact_after vs 0% for SDH). This suggests SDH compensation is largely pre-wired while non-SDH nuclear compensation involves occasional secondary site mutation.

**What this tells us:** The dominant mechanism is structural pre-adaptation — pathogenic variants arise in environments where the surrounding contacts are already physically accommodating. This is consistent with the "permissive structural background" model and argues against the view that compensatory mutations must arise after the DAV to rescue it. Only 6.5% of events are clear rescues; the majority of "compensated" states were already present before the DAV arose.

### 2.9 Temporal Ordering — Mito-Nuclear Interface (New)

**This is the most striking new result.** The 45 Pyvolve-significant mito-nuclear pairs, stratified by direction:

**mt-dar direction (mtDNA DAR → nuclear contact gene), 26 pairs, 867 branch events:**

| Category | Events | % |
|---|---|---|
| contact_first | 843 | 97.2% |
| permissive_background | 3 | 0.3% |
| constitutively_permissive | 6 | 0.7% |
| co_adaptation | 2 | 0.2% |
| contact_after | 0 | 0.0% |
| no_contact_change | 5 | 0.6% |

Permissive signal: 98.3%. Rescue: 0%. Far above the within-genome baseline of 63.2% (binomial p = 5.5×10⁻¹⁴⁵). Every single pair has contact_first as its dominant timing class (26/26).

**nuc-dar direction (nuclear DAR → mtDNA contact gene), 19 pairs, 68 branch events:**

| Category | Events | % |
|---|---|---|
| contact_first | 36 | 52.9% |
| permissive_background | 8 | 11.8% |
| constitutively_permissive | 20 | 29.4% |
| co_adaptation | 1 | 1.5% |
| contact_after | 1 | 1.5% |
| no_contact_change | 2 | 2.9% |

Permissive signal: 94.1%. Rescue: 1.5%. Also significantly above the baseline (p = 3.7×10⁻⁹).

**Interpretation:** Both directions show overwhelming structural pre-adaptation. When an mtDNA DAR variant arises on a phylogenetic branch, the nuclear partner residue was already in the compensatory amino acid state on 97.2% of those branches. This is not a temporal ordering result in the traditional sense — it does not mean the nuclear gene "anticipated" the mtDNA change. Rather, it means nuclear residues at the mito-nuclear interface are under exceptionally strong purifying selection to maintain structural compatibility, so they are almost never in an incompatible state when an mtDNA change occurs.

**Critical caveat:** The 97.2% contact_first rate partly reflects structural conservation at the nuclear interface (these residues rarely change at all on any branch), not just co-evolutionary tracking. The Pyvolve test confirms that cDAV mtDNA positions are enriched for co-occurring alternative AAs in the nuclear gene — this conditions on cDAV vs uDAV status and shows genuine co-variation. The temporal ordering tells us the mechanism is constitutive permissiveness, not active rescue.

### 2.10 Mito-Nuclear Directional Symmetry

The full evidence picture for directionality:

| Test | mt-dar (mtDNA DAR → nuc contact) | nuc-dar (nuc DAR → mt contact) |
|---|---|---|
| sig_any (Fisher p < 0.05 OR branch p < 0.05) | 61.4% | 4.9% |
| Pyvolve sig% (perm_p < 0.05) | 48.1% | 50.0% |
| Temporal permissive rate | 98.3% | 94.1% |

The sig_any asymmetry (61.4% vs 4.9%) is a power artifact: nuclear-DAR pairs have median n_cdav_spp = 9, giving Fisher's exact test near-zero power. Pyvolve, calibrated per-gene, shows symmetric rates (48.1% vs 50.0%). Temporal ordering shows similarly high permissive rates in both directions. **There is no statistical support for directional asymmetry in the compensatory co-evolution signal at the mito-nuclear interface.**

**What this tells us:** The "mtDNA tolerates nuclear variation" and "nuclear accommodates mtDNA variation" hypotheses are not distinguishable with current data. Both directions show the same pattern: structural pre-adaptation dominates. The mito-nuclear interface as a whole is under strong bidirectional purifying selection to maintain compatibility.

### 2.11 Tier-1 Pair Identification

**Tier-1A (15 pairs, 5 non-SDH):** Pagel-significant + permissive timing — primary single-subunit mutagenesis targets.

Non-SDH priority targets:
| Pair | Pagel FDR | Timing | BLOSUM | Notes |
|---|---|---|---|---|
| ATP5MF:78M→80 | 9×10⁻⁶ | permissive_background | 2.0 | CV; strongest non-SDH Pagel signal |
| ATP5MF:93Q→96 | 6.6×10⁻³ | permissive_background | 1.0 | CV; second CV pair |
| NDUFV2:155L→157 | 4.1×10⁻² | permissive_background | 3.0 | CI |
| MT-ND2:122A→123 | 8.1×10⁻² | contact_first | 8.0 | CI mt-mt; strongest BLOSUM |
| UQCRB:13R→10 | 9.5×10⁻² | permissive_background | 2.0 | CIII |

**Tier-1C (38 mito-nuclear pairs):** Pyvolve-significant + contact_first timing — cross-interface targets.

High-confidence targets:
- NDUFS2:6S/20L ↔ MT-ND2:305 (OR = 1.3–2.6, high confidence) — CI interface
- COX6B1:67I ↔ MT-CO2:117 (OR = 9.6, high confidence) — CIV interface
- COX4I2:147T ↔ MT-CO2:191 (OR = 4.3, high confidence) — CIV interface
- COX5B:10A/24V ↔ MT-CO3:154 (OR = 1.7–2.9, high confidence) — CIV interface
- MT-ATP8:39L ↔ ATP5IF1:44/47/50 (OR = 1.6–1.7, high confidence) — CV interface
- MT-ATP8:39L ↔ ATP5PD:85/90 (OR = 1.7–1.8, high confidence) — CV interface

---

## Figure Plan

### Overall Narrative Structure for 6 Figures

```
Fig 1: Scale → What fraction of pathogenic variants are compensated, and where?
Fig 2: DCA → mtDNA has sequence co-evolution signal; nuclear can't be tested (Meff)
Fig 3: Pyvolve → Both genomes show species-distribution co-evolution; nuclear OR genuine
Fig 4: Pagel → Phylogenetic signal spans all complexes; not contact-class specific
Fig 5: Timing → Pre-adaptation dominates within-genome; mito-nuclear interface even stronger
Fig 6: Targets → Tier-1 pairs for experimental validation
```

---

### Figure 1 — Dataset Overview and cDAV Landscape

**Purpose:** Establish scale; the central observation is the mtDNA/nuclear cDAV asymmetry.

**1A — Horizontal stacked bar (2 bars, mtDNA and nucDNA)**
- x-axis: fraction 0–1; bars split cDAV (dark) vs uDAV (light)
- Label counts inside bars: "167 cDAV / 156 uDAV" and "2,046 cDAV / 4,440 uDAV"
- *Message:* mtDNA 51.7% vs nuclear 31.5%
- Data: `data/derived/classified/variants_master_classified.parquet`

**1B — Horizontal grouped bar by complex (CI–CV + mito-nuclear)**
- y-axis: complex; bars for n_cDAV_genes and n_cDAV_DARs per complex, colored by complex
- Separate bars or facets for mtDNA-encoded and nuclear-encoded subunits within each complex
- *Message:* cDAV positions span all complexes; mtDNA subunits of CI and CV are overrepresented
- Data: `data/derived/classified/variants_master_classified.parquet`

**1C — Bubble chart: 100 genes as dots**
- x = n_cDAV positions per gene, y = n_structural_contacts per gene, size = n_testable_Pagel_pairs, color = complex
- Label outliers: SDHC/SDHD (many Pagel pairs), MT-ATP8 (many mito-nuclear pairs), MT-ND5
- *Message:* No single gene dominates; the signal is distributed
- Data: `results/structural/compensatory_partners.csv`

**1D — Stacked bar: contact pair type × cDAV status**
- 3 bars (mt-mt, nuc-nuc, mt-nuc); each bar split cDAV/uDAV
- *Message:* mt-nuc pairs are a small fraction but the question of whether cDAV enrichment is symmetric across directions is addressed in later figures
- Data: `results/phylo/conditional_permissiveness.csv`

---

### Figure 2 — Sequence-Level Co-evolution (DCA)

**Purpose:** Establish mtDNA DCA signal; show Meff as the limiting factor for nuclear.

**2A — mtDNA DCA-DI box/violin plot (9 genes)**
- x: genes ordered by MW p-value (most significant first)
- y: DCA-DI percentile
- Two violins per gene (cDAV dark blue, uDAV light gray), with swarm overplot
- Significance markers (*** p<0.001, ** p<0.01, * p<0.05, ns) above each pair
- *Message:* 8/9 genes show significant enrichment; the signal is consistent across genes
- Data: `results/mutagenesis/dca_all_davs.csv` filtered `dar_genome == "mtDNA"`

**2B — Nuclear Meff vs DCA signal scatter**
- x: per-gene median Meff (log scale, x-axis), y: −log10(MW p-value) for cDAV > uDAV
- One dot per nuclear gene (n = 73); color by complex; size proportional to n_cDAV_pairs
- Horizontal dashed line at p = 0.05; vertical dashed line at Meff = 10
- Label genes with Meff > 10 (NDUFS6, COXFA4, etc.)
- Annotate: Spearman r = 0.37, p = 0.0014
- *Message:* Signal emerges at Meff ≥ 10; most nuclear genes fall below threshold
- Data: `results/mutagenesis/dca_nuclear_meff_stratified.csv` + per-gene MW computation

**2C — mito-nuclear APC-MI boxplot**
- Two grouped boxes: mt-mt pairs (left) and mt-dar mito-nuclear pairs (right)
- Each group: cDAV vs uDAV APC-MI percentile
- *Message:* The mtDNA gene alignment captures mito-nuclear co-variation (p = 2.83e-10 for mt-dar direction); this is the source of the mito-nuclear MI signal
- Data: `results/mutagenesis/dca_all_davs.csv` filtered by pair type and `inter_protein`

---

### Figure 3 — Species-Distribution Co-evolution (Pyvolve)

**Purpose:** Establish Pyvolve enrichment; show rate-controlled OR comparison; confirm non-SDH signal.

**3A — Pyvolve sig% by contact class (grouped bar)**
- x: contact class (mt-mt, nuc-nuc all, nuc-nuc nonSDH, mt-dar, nuc-dar)
- y: % pairs with perm_p < 0.05
- Error bars: 95% binomial CI (Wilson)
- Horizontal dashed line at 5% null
- *Message:* All classes significantly exceed 5% null; nuclear rate matches or exceeds mtDNA after SDH exclusion
- Data: `results/phylo/conditional_permissiveness.csv`

**3B — OR distribution by contact class (violin + swarm)**
- Same 5 classes on x-axis; y: observed OR (log scale)
- Median lines, swarm overplot for small-n classes (mt-nuc)
- Annotate median OR per class
- *Message:* Nuclear ORs are higher (median 1.74–2.04 vs mtDNA 1.06); some outliers at mito-nuclear interface (OR > 9)
- Data: `results/phylo/conditional_permissiveness.csv`

**3C — Rate-controlled OR comparison (line plot)**
- x: n_cdav_spp bins (1–5, 6–20, 21–50, 51–200, >200)
- y: sig% (left) and median OR (right, separate axis or subplot)
- Two lines: mtDNA (dashed, orange) and non-SDH nuclear (solid, blue)
- Annotate: CMH χ² = 6.6, p = 0.010
- *Message:* Nuclear OR remains higher than mtDNA at matched cDAV species count
- Data: `results/phylo/pyvolve_rate_controlled.csv`

**3D — SDH exclusion robustness (paired bars)**
- Two paired bars: nuclear all vs non-SDH nuclear
- Metrics shown: Pyvolve sig%, Pyvolve OR median, Pagel sig count
- *Message:* Signal persists after removing SDH; SDH is a strong contributor but not the driver
- Data: computed from `results/phylo/conditional_permissiveness.csv` and `results/structural/compensatory_partners.csv`

---

### Figure 4 — Phylogenetic Co-evolution (Pagel)

**Purpose:** Show Pagel results; complex stratification; non-SDH nuclear signal.

**4A — Horizontal lollipop: 44 significant pairs**
- y: one row per pair, sorted by complex then pagel_fdr (most significant at top)
- x: −log10(pagel_fdr)
- Color: complex (CI=steelblue, CII/SDH=crimson, CIII=seagreen, CIV=darkorange, CV=mediumpurple)
- Shape: direction (mt-mt=circle, nuc-nuc=square, mt-nuc=diamond)
- *Message:* SDH (CII, crimson) has many significant pairs; non-SDH spans all other complexes
- Data: `results/structural/compensatory_partners.csv`

**4B — Stacked bar by complex**
- x: complex category (CI mtDNA, CI nucDNA, CII/SDH, CIII, CIV, CV mtDNA, CV nucDNA)
- y: count of significant pairs
- Stack by direction/contact type
- *Message:* Non-SDH nuclear signal is distributed across complexes, not concentrated
- Data: `results/phylo/pagel_by_complex.csv`

**4C — Contact-class heatmap: sig vs non-sig**
- 4 rows (hbond, vdW, hydrophobic, electrostatic), 2 columns (sig, non-sig)
- Fill: % of total in that column; annotate counts
- Annotation: χ² = 2.55, p = 0.47
- *Message:* No contact-class enrichment; compensation mechanism is not bond-type specific
- Data: `results/phylo/pagel_by_complex.csv`

---

### Figure 5 — Temporal Ordering

**Purpose:** Show within-genome timing (pre-adaptation dominates); reveal extreme pre-adaptation at mito-nuclear interface.

**5A — Stacked bar: Pagel pairs (3 columns)**
- Columns: All Pagel (n=40), SDH pairs (n=22), Non-SDH (n=18)
- Stacks: 6 refined timing classes (contact_first=deep teal, permissive_background=teal, constitutively_permissive=light teal, co_adaptation=coral, contact_after=red, no_contact_change=gray)
- y-axis: % of branch events
- Annotate permissive % above each bar
- Data: `results/phylo/timing_annotations_v2.csv`

**5B — Stacked bar: mito-nuclear pairs (2 columns)**
- Same color scheme as 5A
- Columns: mt-dar (n=26, 867 events) and nuc-dar (n=19, 68 events)
- Place this panel next to 5A for visual contrast
- The near-100% contact_first bar for mt-dar is the key visual
- Data: `results/phylo/timing_mitonuclear_extended.csv`

**5C — Forest plot: permissive rate comparison**
- y-axis: pair class (All Pagel, SDH, Non-SDH, mt-dar mito-nuc, nuc-dar mito-nuc)
- x-axis: permissive signal % (0–100%)
- Points with 95% binomial CI bars
- Vertical dashed reference line at 63.2% (Pagel-all baseline)
- *Message:* Mito-nuclear pairs (98.3% and 94.1%) far exceed the within-genome baseline (63.2%)
- Data: computed from timing CSVs

**5D — Schematic cartoon (Inkscape)**
- Simplified phylogeny branch showing contact_first scenario
- Parent node: DAR = human ref AA, contact = compensatory AA (already compatible)
- Child node: DAR = pathogenic AA (cDAV state), contact = compensatory AA (unchanged)
- Label timing classes as arrows on the branch
- *Message:* Explain the temporal ordering framework visually for readers unfamiliar with ancestral state reconstruction

---

### Figure 6 — Integrated Tier-1 Pairs and Structural Context

**Purpose:** Translate analysis results into actionable experimental targets.

**6A — Pagel pairs scatter: phylogenetic signal × temporal support**
- x: −log10(pagel_fdr), y: Fisher OR (from compensatory_partners.csv)
- Color: dominant_refined_timing (contact_first/permissive_background = teal; co_adaptation = coral; others = gray)
- Shape: complex
- Size: n_dar_gain_branches
- Tier-1A pairs labeled with gene names
- *Message:* Only pairs in upper-right + teal are the highest-confidence targets
- Data: `results/phylo/tier1_pairs.csv`

**6B — Tier-1A heatmap (15 pairs)**
- Rows: pairs labeled "DAR_gene:pos → contact_gene:pos"
- Columns: pagel_fdr (color scale), BLOSUM contact (color scale), permissive events (color scale), timing_confidence (categorical)
- Side labels: SDH vs non-SDH
- *Message:* 5 non-SDH pairs with high BLOSUM + high permissive counts are the mutagenesis priority
- Data: `results/phylo/tier1_pairs.csv`

**6C — Structural snapshot (PyMOL/ChimeraX)**
- Recommended pair: NDUFS2:6S/20L ↔ MT-ND2:305 in Complex I
- Structure: 9I4I or 9TI4
- Show: NDUFS2 chain with positions 6 and 20 highlighted (spheres, blue); MT-ND2 position 305 highlighted (sphere, orange)
- Draw distance line between Cβ atoms to show contact
- Label with gene names and residue numbers
- *This panel to be generated manually in PyMOL or ChimeraX*

**6D — Evidence matrix (structured table)**
- 3 rows (mt-mt, nuc-nuc nonSDH, mt-nuc)
- 5 columns (DCA/Meff, Pyvolve sig, Rate-controlled OR, Pagel FDR, Temporal permissive)
- Fill: ● strong, ◐ partial, ○ absent (or color: green/yellow/red)
- Annotate key statistics in each cell
- *Message:* Summary of the entire evidence structure for the three genomic contexts

---

## Supplementary Figures

| Figure | Content | Data |
|---|---|---|
| S1 | Full DCA-DI violin for all 100 genes (Meff-labeled) | `dca_all_davs.csv` |
| S2 | Meff distribution for all 100 genes sorted (nuclear vs mtDNA) | `dca_all_davs.csv`, `asr_coordinate_harmonization_audit.tsv` |
| S3 | Pyvolve volcano: OR vs −log10(perm_p) for all 912 pairs, color by contact_type | `conditional_permissiveness.csv` |
| S4 | FoldX ΔΔG: cDAV vs uDAV DARs (exists at `foldx_ddg.csv`) | `results/mutagenesis/foldx_ddg.csv` |
| S5 | Manhattan-style Pagel p-values for all 4,227 tested pairs, colored by complex | `compensatory_partners.csv` |
| S6 | Full tier-1B/C table | `tier1_pairs.csv` |

---

## Key Numbers for Manuscript Text

| Metric | Value |
|---|---|
| Total classified variants | 6,809 (2,213 cDAV, 4,596 uDAV) |
| mtDNA cDAV fraction | 51.7% (167/323) |
| Nuclear cDAV fraction | 31.5% (2,046/6,486) |
| Total DAR-contact pairs | 93,693 |
| DCA enrichment (mtDNA) | MW p = 7.65×10⁻⁶⁹ |
| Nuclear DCA Meff threshold | Meff ≥ 10 → p < 0.01 |
| Genes with Meff ≥ 10 | 13/73 nuclear (18%) |
| Best nuclear DCA gene | NDUFS6 (Meff = 26.8, p = 9.8×10⁻⁴) |
| Pyvolve sig — all | 310/912 (34.0% vs 5% null) |
| Pyvolve sig — mtDNA DAR | 63/262 (24.0%) |
| Pyvolve sig — nuclear nonSDH | 170/444 (38.3%) |
| Pyvolve sig — mt-dar mito-nuc | 26/54 (48.1%) |
| Pyvolve sig — nuc-dar mito-nuc | 19/38 (50.0%) |
| Rate-controlled OR CMH | χ² = 6.6, p = 0.010 |
| Pagel sig (FDR ≤ 0.10) | 44/4,227 |
| Pagel non-SDH nuclear | 18 pairs, 11 genes, all complexes |
| Contact-class chi-sq | χ² = 2.55, p = 0.47 (non-significant) |
| Temporal permissive — all Pagel | 63.2% (127/201 events) |
| Temporal permissive — SDH | 68.1% |
| Temporal permissive — non-SDH | 36.3% |
| Temporal rescue (contact_after) | 6.5% (all); 9.1% (non-SDH) |
| Mito-nuclear pairs temporally ordered | 45 (26 mt-dar, 19 nuc-dar) |
| mt-dar permissive rate | 98.3% (vs baseline 63.2%, p = 5.5×10⁻¹⁴⁵) |
| nuc-dar permissive rate | 94.1% (vs baseline 63.2%, p = 3.7×10⁻⁹) |
| mt-dar rescue rate | 0.0% (0/867 events) |
| Tier-1A (Pagel + permissive timing) | 15 pairs (5 non-SDH) |
| Tier-1C (mito-nuc Pyvolve + contact_first) | 38 pairs |

---

## Data Sources per Figure

| Figure | Primary files |
|---|---|
| Fig 1 | `data/derived/classified/variants_master_classified.parquet`, `results/structural/dar_contacts_cbcb8A.csv`, `results/phylo/conditional_permissiveness.csv` |
| Fig 2 | `results/mutagenesis/dca_all_davs.csv`, `results/mutagenesis/dca_nuclear_meff_stratified.csv`, `data/phylo/asr_coordinate_harmonization_audit.tsv` |
| Fig 3 | `results/phylo/conditional_permissiveness.csv`, `results/phylo/pyvolve_rate_controlled.csv`, `results/structural/compensatory_partners.csv` |
| Fig 4 | `results/structural/compensatory_partners.csv`, `results/phylo/pagel_by_complex.csv` |
| Fig 5 | `results/phylo/timing_annotations_v2.csv`, `results/phylo/timing_mitonuclear_extended.csv` |
| Fig 6 | `results/phylo/tier1_pairs.csv`, `results/mutagenesis/foldx_ddg.csv` (structure: 9I4I or 9TI4) |

---

## Aesthetic Guide

**Color palette (consistent across all figures):**

```
Genomes:
  mtDNA:                  #E07B39  (warm orange)
  nucDNA:                 #4472C4  (steel blue)
  mito-nuclear interface: #D4A017  (gold/amber)

Complexes:
  CI:   #4472C4  (steel blue)
  CII:  #C0392B  (crimson — SDH)
  CIII: #27AE60  (sea green)
  CIV:  #E67E22  (dark orange)
  CV:   #8E44AD  (medium purple)

cDAV vs uDAV within a plot:
  cDAV: filled, darker shade of class color
  uDAV: open/hatched, 50% alpha

Refined timing classes:
  contact_first:             #1B7A6E  (deep teal)
  permissive_background:     #2EA89A  (teal)
  constitutively_permissive: #85C8C2  (light teal)
  co_adaptation:             #E87C4D  (coral/orange)
  contact_after:             #C0392B  (red)
  no_contact_change:         #BFBFBF  (gray)

Statistical annotation: *** p<0.001, ** p<0.01, * p<0.05, ns
```

**Typography:**
- All figures: 300 dpi, Helvetica or Arial
- Single-column width: 7 inches; double-column: 14 inches
- Axis labels: 8 pt; tick labels: 7 pt; figure panel letter: 10 pt bold
- Report exact p-values in figure legends, not just significance stars

**Panel letter convention:** A–D within each figure, upper-left bold (e.g., **A**, **B**)
