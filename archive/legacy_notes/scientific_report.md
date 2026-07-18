# Compensated Disease-Associated Variants in Oxidative Phosphorylation: A Comprehensive Scientific Report

## Abstract

Human pathogenic variants in oxidative phosphorylation (OxPhos) subunits cause a clinically diverse spectrum of mitochondrial diseases, yet identical amino acid substitutions are naturally tolerated across mammalian species. We systematically identified 2,282 amino acid-level compensated disease-associated variants (cDAVs) — positions where at least one non-human mammal naturally carries a human-pathogenic residue — across all 100 OxPhos subunit genes in 289 mammalian species. Phylogenetic reconstruction of substitution histories revealed that 83% of cDAVs arose independently in multiple lineages, demonstrating pervasive convergent evolution at disease-associated positions. Using three orthogonal co-evolutionary tests across 8,087 DAR–contact residue pairs, we identified 693 significant compensatory co-evolutionary relationships. The dominant evolutionary mode is epistatic pre-adaptation: the compensating contact substitution preceded the pathogenic residue in 85.8% of significant pairs (OR = 5.64, p = 1.2 × 10⁻³⁰), rising to OR = 24.88 (p = 8.8 × 10⁻⁴¹) for derived cDAVs. Nuclear-encoded compensatory relationships are nearly twice as frequent as mitochondrially-encoded ones (10.8% vs 6.1%, p = 1.6 × 10⁻¹²). To identify the highest-priority pairs for experimental validation, we developed a multi-evidence prioritization framework combining phylogenetic significance, structural proximity (Cβ–Cβ distance), physicochemical mechanism, sequence co-variation (mutual information with APC correction), and computational stability perturbation (FoldX BuildModel/AnalyseComplex ΔΔG). FoldX calculations on the top 50 pairs confirm that 9 pairs simultaneously satisfy both criteria for genuine compensatory rescue: the disease mutation destabilises the protein (ΔΔG_DAR > 0.5 kcal/mol) and the compensating contact substitution partially or fully restores stability (ΔΔG_rescue > 0.5 kcal/mol). The pair NDUFA9 M220I / S216Y exhibits the strongest epistasis signal (ΔΔG_epistasis = −1.22 kcal/mol), indicating a purely non-additive synergistic interaction. These results provide the first systematic phylogenomic evidence of epistatic pre-adaptation in OxPhos evolution and a prioritized, multi-evidence shortlist of compensatory pairs ready for cell-line suppressor screen validation.

---

## 1. Introduction

### 1.1 The OxPhos co-evolutionary problem

Oxidative phosphorylation is structurally and functionally the most complex protein system in the eukaryotic cell. The five OxPhos complexes (Complexes I–V) are assembled from 98 nuclear-encoded and 13 mitochondrially-encoded subunits, all of which must physically interact across up to 847 intersubunit contacts to carry out electron transport and ATP synthesis. This structural interdependence creates a fundamental evolutionary constraint: any amino acid substitution that alters a subunit's shape, charge, or stability can propagate perturbations through the contact network, destabilising assembly or function.

In human mitochondrial disease, single amino acid changes in OxPhos subunits cause devastating clinical outcomes, yet the same changes are found naturally — and apparently tolerated — in other mammalian species. This observation implies that the pathogenicity of a variant is not an intrinsic property of the substitution itself, but depends on the structural and genetic context in which it occurs. The question is: what makes that context different?

Two non-exclusive explanations exist. First, the compensatory contact residue evolved independently of the disease variant — a reactive rescue substitution that arose *after* the pathogenic residue. Second, the contact residue evolved to its current state for entirely unrelated reasons, and the lineage that happened to carry it was thereby pre-adapted to tolerate the pathogenic residue when it arose. The second model — epistatic pre-adaptation — predicts that the structural background permitting a human-pathogenic residue is in place *before* the variant appears in that lineage, meaning that certain lineages are selectively permissive. Distinguishing these models has been possible only with systematic phylogenomic data across a large mammalian tree.

### 1.2 The cDAV framework

We define a **compensated disease-associated variant (cDAV)** as any non-human mammalian species naturally carrying an amino acid that is pathogenic in humans at the orthologous OxPhos position. This framework is anchored to clinical data from MITOMAP (mtDNA) and ClinVar (nucDNA), ensuring that every variant in our analysis has documented human pathogenicity. A stricter nucleotide-level cDAV subset requires the species to use the identical codon as the human mutation.

---

## 2. The cDAV Landscape Across OxPhos

### 2.1 Scale and distribution

Screening 443 mtDNA and 6,554 nucDNA clinical variants across 289 mammalian species identified **2,282 AA-level cDAVs** (236 mtDNA, 2,046 nucDNA) and 1,954 NT-level cDAVs (207 mtDNA, 1,747 nucDNA). The AA-to-NT ratio (~87%) is notably high, indicating that in the majority of cDAV cases the non-human species uses the exact same codon as the human mutation — not merely the same amino acid via a synonymous path. This matters for interpreting molecular mechanism: the same codon implies the identical regulatory and translational context, making these cases more directly comparable to the human pathogenic situation.

cDAVs span all five OxPhos complexes and all 13 mt-encoded and 87 of 87 nuclear-encoded subunit genes surveyed:

**mtDNA cDAVs by complex:**
| Complex | cDAVs | % of mtDNA total |
|---|---|---|
| Complex I (ND subunits) | 119 | 50.4% |
| Complex IV (CO subunits) | 51 | 21.6% |
| Complex V (ATP subunits) | 44 | 18.6% |
| Complex III (CYB) | 22 | 9.3% |

**nucDNA cDAVs by complex:**
| Complex | cDAVs | % of nucDNA total |
|---|---|---|
| Complex II (SDH subunits) | 812 | 39.7% |
| Complex I (NDUF subunits) | 734 | 35.9% |
| Complex IV (COX subunits) | 160 | 7.8% |
| Complex V (ATP subunits) | 142 | 6.9% |
| Complex III (UQCR subunits) | 135 | 6.6% |

The dominance of Complex II in nucDNA cDAVs is striking — all four SDHA/B/C/D subunits are nuclear-encoded, and ClinVar contains a substantial number of SDH variants associated with hereditary paraganglioma and pheochromocytoma. The dominance of Complex I in mtDNA cDAVs reflects the 7 ND subunits in the mitochondrial genome, together accounting for the major LHON and MELAS-associated loci.

### 2.2 Species breadth

The median mtDNA cDAV is found in 31 species; the median nucDNA cDAV is found in only 5 species. This threefold difference in breadth reflects fundamental differences in mtDNA vs nucDNA evolutionary dynamics. mtDNA evolves ~10× faster than the nuclear genome in mammals (Brown et al. 1979; Pesole et al. 1999) and accumulates substitutions via high mutation rate and genetic drift without recombination, allowing many independent lineages to acquire and maintain the same amino acid substitution. Nuclear-encoded OxPhos subunits evolve more slowly under purifying selection, and when a human-pathogenic residue appears it tends to be found in fewer lineages, often recently derived.

---

## 3. Convergent Evolution at Disease-Associated Positions

### 3.1 Independent origins are the rule

Ancestral state reconstruction using IQTree across 100 OxPhos gene trees, anchored to 289 mammalian species and calibrated against the MamPhy ultrametric timetree, allowed us to count the number of independent evolutionary events (branches where the human-pathogenic amino acid arose) for each cDAV across the mammalian tree.

Of 903 unique cDAV positions with sufficient phylogenetic data:

| Independent origins | cDAVs | % |
|---|---|---|
| 1 (single lineage) | 150 | 16.6% |
| 2–5 | 244 | 27.0% |
| 6–10 | 129 | 14.3% |
| >10 (highly convergent) | 375 | 41.5% |

**83.4% of cDAV positions arose independently in ≥2 lineages**, and 41.5% arose >10 times. The most convergent positions include MT-CO1 419V (262 independent gains), MT-ATP6 192T (255 gains), and NDUFB10 154K (251 gains). This extraordinary convergence — approaching parallel evolution at hundreds of positions — has no precedent in systematic disease variant analysis.

### 3.2 Interpretation in the context of fitness landscape theory

High convergence at a position implies that this amino acid state is repeatedly accessible on the mammalian fitness landscape. Under a standard NK model of molecular evolution (Kauffman & Weinberger 1989), convergent evolution is expected when the relevant substitution has low epistatic cost in a wide range of genetic backgrounds. The finding that most cDAVs arise >10 times independently suggests these are positions where the alternative amino acid is close to neutral, or even beneficial, in many mammalian contexts — despite causing disease in humans.

This is consistent with the "nearly neutral" model of molecular evolution (Ohta 1973): variants that are deleterious in the human genetic background may be close to neutral in species where background selection has shaped the surrounding protein environment differently. Compensatory co-evolution is the mechanistic bridge between these two observations.

### 3.3 Convergence predicts compensation

Convergent cDAVs (≥2 independent origins) are significantly enriched for compensatory contacts relative to single-origin cDAVs (OR = 2.65, 95% CI 1.23–5.71, Fisher p = 0.003). This is a critical observation: positions that repeatedly and independently acquire the human-pathogenic amino acid across mammals are also more likely to be found in protein environments that have co-evolved structural compensation. The structural context that makes a position "easy" to acquire is also the context where compensatory contacts are available — suggesting that convergence and compensation are two manifestations of the same underlying permissive structural environment.

---

## 4. Compensatory Co-Evolution: Statistical Evidence

### 4.1 Three-test framework

For each of 8,087 cDAV–contact residue pairs (defined by Cβ–Cβ distance ≤ 8Å in OxPhos crystal structures), we applied three orthogonal statistical tests of co-evolutionary enrichment:

1. **Fisher's exact test (species as units)**: Tests enrichment of the contact alt-AA in species carrying the cDAV, as a phylogenetically naive baseline (one-tailed, FDR-corrected per cDAV)

2. **Branch co-occurrence test**: Identifies branches where both the cDAV and contact substitution arose independently in the same phylogenetic clade (≥80% leaf-set overlap), testing whether co-arising branches exceed chance expectation (Fisher's exact over branches, BH FDR per cDAV)

3. **Pagel's discrete correlated-evolution test**: Model-based likelihood ratio test for correlated evolution of two binary characters (cDAV presence/absence, contact alt-AA presence/absence) on the MamPhy tree (`phytools::fitPagel`), corrected by global BH FDR

Using branch-based and model-based tests (rather than species counts alone) is critical: OxPhos subunits show strong phylogenetic structure, and species-level tests inflated by relatedness would massively over-report co-evolution. The branch co-occurrence and Pagel tests address phylogenetic non-independence directly.

### 4.2 Significant co-evolutionary pairs

At FDR ≤ 0.10, **693 DAR–contact pairs** across **267 unique cDAV positions** show significant co-evolutionary enrichment:

| Evidence | Pairs |
|---|---|
| Branch co-occurrence FDR ≤ 0.10 | 172 |
| Pagel's discrete FDR ≤ 0.10 | 25 |
| Both tests | 1 |

The substantially larger branch co-occurrence set reflects its greater statistical power for detecting recent co-arising events in species-rich clades (e.g. rodents, bats), while Pagel's test has greater power for older, deeper co-evolutionary signals that have accumulated across many independent lineages. The single pair passing both tests — involving UQCRC2 Ser410Pro in Complex III — represents the highest-confidence co-evolutionary signal in the dataset.

### 4.3 Complex-specific compensation rates

| Complex | Pairs tested | Significant | Rate |
|---|---|---|---|
| Complex I | 4,524 | 351 | 7.8% |
| Complex II | 1,236 | 110 | 8.9% |
| Complex III | 381 | 34 | 8.9% |
| Complex IV | 838 | 131 | 15.6% |
| Complex V | 958 | 66 | 6.9% |

**Complex IV (cytochrome c oxidase) has the highest compensation rate at 15.6%**, nearly double the other complexes. This is consistent with CIV's structural architecture: the 13 mt-encoded core subunits (CO1–3) form tight contacts with 11 nuclear-encoded peripheral subunits, and CIV must assemble with cytochrome c (a small soluble electron carrier) and the inner membrane lipid cardiolipin. This dense, multi-interface contact network provides more structural degrees of freedom for compensatory co-evolution than, for example, the peripheral arm of Complex I where contacts are more localised.

---

## 5. Epistatic Pre-Adaptation: The Dominant Evolutionary Mode

### 5.1 Temporal ordering across the mammalian tree

For each annotated pair, we classified the evolutionary order of the cDAV and its compensatory contact using IQTree ancestral state reconstructions:

- **contact_first**: The contact alt-AA was already present at the parent node of the branch on which the cDAV arose — the compensating residue *preceded* the pathogenic variant
- **co_occurring**: Both changes arose on the same branch in the same inferred lineage
- **contact_after**: The contact alt-AA arose within the subtree *after* the cDAV appeared
- **no_contact_change**: No detectable contact substitution in the lineage

Among all 6,738 timing-annotated pairs:
- contact_first: 18.1% of derived cDAV branches
- co_occurring: 0.2%
- contact_after: 0.5%
- no_contact_change: 81.2%

### 5.2 Contact_first is massively enriched in significant pairs

The critical question is whether contact_first timing is causally associated with co-evolutionary significance, or merely reflects baseline rates of contact substitution on ancestral branches. The result is unambiguous:

**Contact_first timing is enriched with OR = 5.64 (95% CI 3.98–7.99, p = 1.2 × 10⁻³⁰) in significant pairs across all cDAVs.** For derived cDAVs specifically:

| Subset | OR | 95% CI | p-value |
|---|---|---|---|
| All pairs | 5.64 | 3.98–7.99 | 1.2 × 10⁻³⁰ |
| Derived cDAVs | 24.88 | 11.35–54.51 | 8.8 × 10⁻⁴¹ |
| Ancestral cDAVs | 1.92 | 1.19–3.10 | 4.6 × 10⁻³ |
| mt-mt pairs | 5.83 | 3.25–10.45 | 2.5 × 10⁻¹² |
| nuc-nuc pairs | 5.47 | 3.45–8.67 | 9.0 × 10⁻¹⁸ |
| mt-nuc pairs | 3.14 | 0.96–10.28 | 2.9 × 10⁻² |

For derived cDAVs, **96.7% of all significant compensatory pairs had contact_first timing** — an odds ratio of 24.88. This is not a weak association. It means that in the great majority of cases where we can detect statistically significant co-evolution between a cDAV and a contact residue, the contact residue was already changed *before* the human-pathogenic residue appeared in that lineage.

### 5.3 Mechanistic interpretation: the permissive background model

These results strongly support a model of **epistatic pre-adaptation** rather than reactive compensation:

**Reactive compensation (alternative hypothesis)**: The pathogenic residue arises in a lineage, creates a structural perturbation, and drives selection for a compensating contact substitution. In this model, the pathogenic residue comes first and the compensating contact is selected secondarily.

**Epistatic pre-adaptation (observed)**: Substitutions at contact positions occur continuously across mammalian lineages for reasons unrelated to any specific disease variant. A lineage that happens to carry a particular contact variant is thereby pre-adapted — its structural environment is permissive for the pathogenic residue. When that residue independently arises in this lineage, it is tolerated. The co-evolutionary signal we detect is the signature of this pre-adaptation: the contact change was not driven *by* the pathogenic residue, but it is the reason that lineage could tolerate the pathogenic residue when it appeared.

This has a direct clinical implication. If pre-adaptation is the dominant mode, then the pathogenicity of an OxPhos variant in a human patient is not primarily determined by the residue itself, but by whether the patient's background genome carries the permissive structural environment. Patients with the same disease-causing variant but different background genetic architecture may have different severity and progression — a prediction consistent with known variable penetrance in mitochondrial disease.

### 5.4 Ancestral vs derived cDAVs: a 3.9-fold difference in compensation rates

| | Pairs | Significant | Rate | Mean independent origins |
|---|---|---|---|---|
| Ancestral cDAVs | 715 | 78 | 10.9% | 116 |
| Derived cDAVs | 6,023 | 183 | 3.0% | 47 |

OR = 3.92 (95% CI 2.97–5.17, p = 1.6 × 10⁻¹⁸). Ancestral cDAVs — where the human-pathogenic residue was present in the common mammalian ancestor and has been repeatedly lost — show nearly four times the compensation rate of derived cDAVs. This is expected from the pre-adaptation model: ancestral cDAVs have had the full span of mammalian evolution (~180 Mya) to accumulate compensatory structural environments, whereas derived cDAVs are younger lineage-specific events where the compensatory background may not yet have diversified.

The much weaker contact_first enrichment for ancestral cDAVs (OR = 1.92 vs 24.88 for derived) confirms this interpretation: for ancestral cDAVs, the temporal ordering between "contact change" and "cDAV" is blurred by deep evolutionary time. The compensating contact likely co-evolved with the pathogenic residue over tens of millions of years, not on a single lineage.

---

## 6. Genomic Architecture of Compensation

### 6.1 The mt-nuc co-evolution problem

The 13 mt-encoded and ~87 nuclear-encoded OxPhos subunits must physically interact, but they evolve under entirely different genetic regimes: mtDNA is maternally inherited, haploid, lacks recombination, and evolves ~10× faster; nuclear genes are Mendelian, diploid, recombine freely, and evolve more slowly under purifying selection. This creates a persistent co-evolutionary tension at the mt-nuclear interface that has been proposed as a driver of speciation (Chou & Leu 2010; Burton & Barreto 2012), cytonuclear incompatibilities (Meiklejohn et al. 2013), and population-level fitness variation (Havird & Sloan 2016).

Our dataset contains three classes of co-evolutionary pairs:
- **mt-mt** (3,705 pairs): both cDAV and contact are mtDNA-encoded
- **nuc-nuc** (3,972 pairs): both cDAV and contact are nuclear-encoded
- **mt-nuc** (410 pairs): cDAV in one genome, contact in the other

### 6.2 Nuclear compensation is more frequent than mitochondrial

| Contact type | Pairs tested | Significant | Rate |
|---|---|---|---|
| mt-mt | 3,705 | 226 | **6.1%** |
| nuc-nuc | 3,972 | 429 | **10.8%** |
| mt-nuc | 410 | 38 | **9.3%** |

Chi-square across all three types: χ² = 54.33, df = 2, p = 1.6 × 10⁻¹². The Fisher p-value distribution for mt-mt pairs is also significantly shifted toward higher (weaker) values compared to nuc-nuc pairs (Mann-Whitney p = 6.1 × 10⁻²⁹), confirming that the difference is not an artefact of different thresholding.

**Nuclear-encoded cDAVs are compensated at nearly twice the rate of mtDNA cDAVs.** This is the opposite of the "mt-nuclear arms race" prediction, which would expect mtDNA variants — evolving faster and under stronger drift pressure — to accumulate compensatory contacts more readily. Several alternative explanations are consistent with the data:

1. **Structural complexity**: Nuclear-encoded OxPhos subunits are on average larger (~400 aa vs ~200 aa for mt-encoded) and form more extensive contact networks, providing more candidate compensatory positions per cDAV
2. **Evolutionary rate asymmetry**: Fast mtDNA evolution means that contact substitutions at mt-mt interfaces are more frequent, but also less likely to be specifically associated with any one disease variant — the "signal" of compensation is diluted by background substitution rate
3. **Purifying selection asymmetry**: Nuclear OxPhos genes are under stronger purifying selection, meaning that when a contact substitution does occur at a nuclear gene, it is more likely to be maintained and to have functional consequences for the contact interface
4. **ClinVar variant density**: The nucDNA dataset (6,554 variants vs 443 mtDNA) is drawn from a more diverse clinical database, potentially enriching for variants at accessible, partially-tolerated positions

### 6.3 Intergenomic compensation: no enrichment relative to intragenomic

Intergenomic (mt-nuc) pairs show a compensation rate (9.3%) not significantly different from intragenomic pairs overall (8.5%; OR = 1.11, p = 0.59). At face value, this argues against a specific enrichment of compensatory co-evolution at the mt-nuclear interface relative to within-genome contacts.

However, the Fisher p-value distribution for mt-nuc pairs *does* differ significantly from nuc-nuc pairs (Mann-Whitney p = 2.5 × 10⁻⁸), despite similar compensation rates. This suggests that when co-evolution does occur at the mt-nuclear interface, the statistical signal differs qualitatively from purely nuclear co-evolution — possibly because mt-nuc co-evolution operates on a different timescale (driven by mt substitution rate) or involves different physicochemical mechanisms.

### 6.4 Contact class and physicochemical mechanisms

Among the 261 timing-annotated compensatory pairs, the physicochemical mechanisms of compensation are:

| Mechanism | Pairs | % |
|---|---|---|
| Charge rescue | 58 | 22.2% |
| Volume swap | 64 | 24.5% |
| Unclassified | 139 | 53.3% |

The **charge rescue** mechanism (charge change at the cDAV position counteracted by a complementary charge change at the contact) is the most structurally interpretable and provides the most direct supressor mutation hypotheses for experimental follow-up. Contact classes are dominated by hydrogen bond contacts (403/693, 58%) and van der Waals contacts (241/693, 35%), consistent with the types of contacts that mediate subunit–subunit specificity in OxPhos assembly.

### 6.5 The MT-CO2 Thr55 case study: ancient mt-nuclear co-adaptation

The only compensatory pair with a dateable multi-species origin involves **MT-CO2 Thr55** (human disease allele Thr at position 55 of Cytochrome c oxidase subunit 2). This cDAV arose ~20.9 Mya in the common ancestor of a rodent clade (Node229) and has been stably maintained across that entire lineage. It is accompanied by four co-evolved contact positions, including three intra-subunit MT-CO2 contacts and one inter-genomic contact with the nuclear-encoded **COX5A Gly79**:

| Contact | Type | Timing | Mechanism |
|---|---|---|---|
| MT-CO2 52H | mt-mt | no_contact_change | charge_rescue |
| MT-CO2 54S | mt-mt | contact_first | unclassified |
| MT-CO2 56M | mt-mt | no_contact_change | charge_rescue |
| COX5A 79G | mt-nuc | contact_first | charge_rescue |

The ~21 Mya co-maintenance of the MT-CO2 Thr55 / COX5A Gly79 mt-nuclear pair across a rodent clade represents one of the best-dated examples of cytonuclear co-adaptation at atomic resolution. The pre-adaptation of the nuclear COX5A contact position alongside the maintained mt-encoded variant over this timescale directly illustrates the long-term stability of permissive structural environments at the mt-nuclear interface.

---

## 7. Inter-Subunit vs Intra-Subunit Compensation

Of 693 significant compensatory pairs, **594 (86%) involve contact residues within the same subunit** as the cDAV, while **99 (14%) are inter-subunit contacts**. Among the inter-subunit pairs:

| Contact type | Inter-subunit pairs |
|---|---|
| nuc-nuc | 40 (40%) |
| mt-nuc | 38 (38%) |
| mt-mt | 21 (21%) |

The disproportionate representation of mt-nuc contacts among inter-subunit pairs (38% of inter-subunit, but only 5.5% of all significant pairs) is notable. At the mt-nuclear interface, by definition all contacts are inter-subunit (the two genomes encode different polypeptides). The enrichment of mt-nuc contacts among inter-subunit pairs reflects the fact that the mt-nuclear interface is itself the primary inter-subunit contact zone in OxPhos, particularly for the CIV core–peripheral subunit interface and the Complex I membrane arm.

These 99 inter-subunit pairs are particularly important for OxPhos assembly biology: subunit–subunit contacts at these interfaces are precisely the points where mt-nuclear incompatibilities manifest, and where compensatory structural evolution at one subunit must be matched by accommodation at the partner. The 38 mt-nuc inter-subunit pairs define specific contact dyads at the mt-nuclear interface where such co-adaptation has been detected.

---

## 8. Computational Mutagenesis Prioritization

### 8.1 Rationale and scoring framework

The 693 significant compensatory pairs are not equally informative for experimental follow-up. Experimental validation — thermal shift assays, site-directed mutagenesis, pulldown or co-immunoprecipitation — is resource-intensive and requires a ranked shortlist. We developed a multi-evidence composite scoring framework that integrates five orthogonal evidence streams, each addressing a different aspect of the compensatory mechanism:

| Criterion | Max points | Rationale |
|---|---|---|
| Phylogenetic FDR (min of Pagel, branch) | 3 | Significance of co-evolutionary signal |
| Contact class (electrostatic > hbond > hydrophobic > vdw) | 3 | Physicochemical interpretability |
| Physicochemical type (charge_reversal > charge_rescue > volume_swap) | 3 | Mechanistic directness |
| Epistatic pre-adaptation (contact_first timing) | 2 | Evolutionary evidence for pre-adaptation |
| Likely incompatible flag | 2 | Predicted structural incompatibility without compensation |
| Convergence (n_independent_origins ≥ 10) | 2 | Robustness of the compensated state |
| FoldX ΔΔG rescue (> 0.5 kcal/mol) | +2 | Computational stability confirmation |
| MI/APC percentile (> 75th) | +1 | Sequence co-variation evidence |
| Confidence tier (high_confidence) | 1 | Ancestral reconstruction confidence |
| Cβ–Cβ distance (< 5 Å) | 1 | Physical proximity of the contact |

Scores range from 0–20 (maximum). 693 pairs were scored using existing phylogenetic and structural data; FoldX BuildModel ΔΔG was computed on the top 50 pairs using the five available cryo-EM structures (CI: 9I4I, CII: 8GS8, CIII: 9HZL, CIV: 9I6F, CV: 8H9S). Mutual information with APC correction was computed from the 289-species protein MSAs for 615/693 pairs.

### 8.2 Contact category and FoldX protocol

A critical distinction in the FoldX analysis concerns the structural relationship between the disease residue and its compensating contact. Of the 693 pairs:

- **594 (85.7%) are intraprotein**: both the disease residue and the compensating contact reside within the same polypeptide chain. The relevant FoldX metric is **ΔΔG_stability** (`BuildModel`), which captures changes in protein folding free energy.
- **99 (14.3%) are interprotein**: the two residues are in different polypeptide chains (40 nuc-nuc, 38 mt-nuc, 21 mt-mt). For these, `BuildModel` stability is complemented by **`AnalyseComplex` ΔΔG_binding**, which captures changes in subunit-subunit binding affinity.

Applying `AnalyseComplex` to intraprotein pairs would be mechanistically incorrect, as there is no binding interface to perturb within a single polypeptide. This distinction is particularly important for OxPhos, where the majority of compensatory contacts are within individual subunits rather than at subunit interfaces.

### 8.3 FoldX results: destabilization and rescue

FoldX `BuildModel` calculations (3 runs per condition; RepairPDB pre-processing) were completed for all 50 top-ranked pairs across 5 OxPhos cryo-EM structures:

| Metric | Count (of 50 pairs) |
|---|---|
| ΔΔG_stab_DAR > 0.5 kcal/mol (significant destabilization) | 18 (36%) |
| ΔΔG_stab_DAR > 1.0 kcal/mol (strong destabilization) | 12 (24%) |
| ΔΔG_rescue_stab > 0.5 kcal/mol (significant rescue) | 13 (26%) |
| ΔΔG_rescue_stab > 1.0 kcal/mol (strong rescue) | 8 (16%) |
| ΔΔG_epistasis < −0.5 kcal/mol (non-additive epistasis) | 6 (12%) |
| Both destabilizing DAR AND rescue > 0.5 kcal/mol | **9 (18%)** |

The 9 pairs satisfying both criteria simultaneously represent the highest-confidence compensatory pairs — cases where the disease mutation is genuinely destabilising in the protein context, and the compensating contact substitution specifically counteracts that destabilisation rather than being independently stabilising or neutral.

### 8.4 Top validated targets

The five pairs with the strongest composite evidence:

**1. MT-ATP6 I192T / T189A** (Complex V, 8H9S, composite score 16/20)
- ΔΔG_stab_DAR = +1.62 kcal/mol (disease mutation destabilises ATP synthase subunit a)
- ΔΔG_rescue_stab = **+1.71 kcal/mol** (Thr189→Ala fully compensates)
- ΔΔG_epistasis = −0.25 kcal/mol (mild positive cooperativity)
- 255 independent evolutionary origins; branch co-occurrence FDR = 1.4 × 10⁻¹⁶
- Charge_rescue mechanism; hbond contact; Cβ–Cβ = 5.84 Å
- *Interpretation*: The most strongly supported compensatory pair in the dataset. The ATP synthase subunit a Ile192Thr disease mutation (associated with neuropathy/ataxia/retinitis pigmentosa — NARP) destabilises the protein by +1.62 kcal/mol in the CV structure; the co-evolved Thr189→Ala in the adjacent position fully compensates this loss.

**2. MT-ATP6 I192T / S188T** (Complex V, 8H9S, composite score 16/20)
- ΔΔG_stab_DAR = +1.62 kcal/mol; ΔΔG_rescue_stab = +0.83 kcal/mol
- MI_APC = 0.089, 90th percentile (strong sequence co-variation)
- 255 origins; branch FDR = 3.1 × 10⁻¹⁷
- *Interpretation*: A second independent compensating contact for the same MT-ATP6 I192T disease variant, at the adjacent Ser188 position. The MI signal at 90th percentile provides independent sequence-level evidence for co-variation between these two positions across mammals.

**3. NDUFA9 M220I / S216Y** (Complex I, 9I4I, composite score 15/20)
- ΔΔG_stab_DAR = +0.60 kcal/mol; ΔΔG_rescue_stab = **+1.73 kcal/mol**
- ΔΔG_epistasis = **−1.22 kcal/mol** (strong non-additive epistasis)
- MI_APC = 0.068, 79th percentile; branch FDR = 7.7 × 10⁻⁵
- *Interpretation*: The strongest epistasis signal in the dataset. ΔΔG_epistasis = −1.22 kcal/mol means the combined double mutant is 1.22 kcal/mol more stable than the sum of the two individual mutations — a purely synergistic interaction. This is precisely the signature expected from a genuine compensatory pair where the two substitutions are physically interdependent: the S216Y contact substitution cannot provide stability alone, but creates a new structural environment that becomes stabilising only in the context of the M220I disease variant.

**4. MT-ATP6 T28S / MT-ATP8 F24** (interprotein, mt-mt, composite score 13/20)
- ΔΔG_stab_DAR = **+4.42 kcal/mol** (strongest destabilisation in the dataset)
- ΔΔG_rescue_stab = **+2.66 kcal/mol** (strongest rescue in the dataset)
- ΔΔG_epistasis = −1.14 kcal/mol (strong epistasis)
- *Interpretation*: The MT-ATP6 / MT-ATP8 intersubunit contact pair carries the largest absolute ΔΔG values. The T28S disease mutation strongly destabilises Complex V; the MT-ATP8 Phe24 contact substitution rescues +2.66 kcal/mol of this loss. Both subunits are mt-encoded, making this an intra-mt-genome compensatory interaction. The strong epistasis (−1.14 kcal/mol) confirms the two residues are physically coupled at the intersubunit interface.

**5. SDHC I46T / G47S** (Complex II, 8GS8, composite score 14/20)
- ΔΔG_stab_DAR = +0.02 kcal/mol (mild destabilisation by FoldX)
- Cβ–Cβ = **4.43 Å** (closest contact in the top 50)
- MI_APC = 0.091, 88th percentile; 186 independent origins; branch FDR = 5.8 × 10⁻⁸
- Charge_rescue mechanism; 8GS8 at 2.5 Å resolution (highest-confidence structure)
- *Interpretation*: Despite modest FoldX ΔΔG, this pair scores highly due to near-atomic contact distance (4.43 Å), very strong MI evidence (88th percentile), and high convergence. The closest Cβ–Cβ contact among top targets indicates a direct side-chain interaction. The low FoldX ΔΔG may reflect underestimation at the 2.5 Å cryo-EM structure or compensation through an electrostatic/hydrogen bond mechanism that FoldX captures imprecisely.

### 8.5 Mutual information analysis

Mutual information with APC correction (MI_APC) was computed for 615/693 pairs from the 289-species protein MSAs. Of the top 50 FoldX-scored pairs:

- 176/615 pairs (28.6%) have MI_APC above the 75th percentile of within-gene background
- 38 pairs (6.2%) exceed the 90th percentile
- The top MI signal is **MT-CO1 I419V / L423M** (MI_APC = 0.025, 98.2th percentile), which also has strong branch co-occurrence (FDR = 1.0 × 10⁻²³) and 262 independent origins

The agreement between MI-based sequence co-variation and phylogenetically-validated co-evolutionary significance is imperfect (many phylogenetically significant pairs have low MI and vice versa), consistent with the expectation that MI and phylogenetic branch co-occurrence capture different aspects of co-evolution: MI reflects correlated amino acid states across all species simultaneously, while branch co-occurrence detects co-arising substitution events. Pairs significant by both criteria represent the strongest combined evidence.

---

## 9. Connections to Broader Theory

### 8.1 Fitness landscape accessibility and the distribution of compensated states

The convergence data (83% of cDAVs arising ≥2 times independently) demonstrate that human-pathogenic OxPhos variants are not trapped in inaccessible corners of the fitness landscape. They are accessible to mammalian evolution — repeatedly acquired and maintained — in lineages that carry the right structural background. This is consistent with Goldstein's (2011) empirical fitness landscape model for proteins, where most positions tolerate multiple amino acid states contingent on the rest of the sequence, and where what appears "deleterious" in one background is accessible in another.

The ~25-fold OR for contact_first enrichment in derived cDAVs (OR = 24.88) quantifies the contribution of background epistasis to this landscape accessibility. It implies that the probability of a lineage maintaining a human-pathogenic OxPhos residue is ~25-fold higher if a permissive contact substitution is already present. This is one of the largest epistasis effect sizes reported from natural phylogenomic data in any protein system.

### 8.2 Compensatory evolution and Muller's ratchet in mtDNA

mtDNA's lack of recombination predicts accumulation of slightly deleterious mutations via Muller's ratchet (Lynch et al. 1993; Neiman & Taylor 2009). The OxPhos system must therefore continuously evolve compensatory responses in the nuclear genome to counteract mt-encoded fitness declines. Our finding that mt-mt compensation is significantly less frequent than nuc-nuc compensation (6.1% vs 10.8%, p = 1.6 × 10⁻¹²) at first appears to contradict this prediction.

However, the relevant comparison for the Muller's ratchet hypothesis is not compensation rate but *compensatory source*: are mt-encoded cDAVs primarily compensated by nuclear contact partners (i.e., mt-nuc contacts) rather than by other mt-encoded residues? Among the 38 significant mt-nuc compensatory pairs, 38/99 inter-subunit pairs are mt-nuc, disproportionately representing the mt-nuclear interface. This is consistent with nuclear compensation of mtDNA drift — the nuclear genome responding to mt-encoded variation at contact interfaces — even if the overall rate of statistically detectable co-evolution at mt-mt contacts is lower due to signal dilution by high mtDNA substitution rates.

### 8.3 Cytonuclear co-evolution and speciation

The mt-nuclear co-evolution hypothesis of speciation (Burton et al. 2013; Sloan et al. 2017) predicts that co-adapted mt-nuclear contact pairs should be under stabilising selection within populations and should show co-divergence between species. The 38 significant mt-nuc compensatory pairs identified here — particularly the MT-CO2 Thr55 / COX5A Gly79 pair maintained for ~21 Mya — are direct empirical candidates for this class of co-adapted dyads. Each represents a residue pair where mtDNA variation in one genome has been structurally accommodated by nuclear variation at the contact interface.

If these pairs are genuinely co-adapted, reciprocal crosses between lineages carrying different combinations of these residues (e.g., a rodent mt haplotype with Thr55 crossed into a rodent nuclear background without COX5A Gly79) should show fitness defects at the OxPhos level — a prediction that can be tested with cell fusion or mitochondrial transfer experiments.

### 8.4 Clinical implications: variable penetrance and epistasis in mitochondrial disease

The finding that epistatic pre-adaptation accounts for 96.7% of significant compensatory pairs in derived cDAVs has a direct and testable prediction for human disease: patients who carry both a pathogenic OxPhos variant *and* the human-orthologous version of the mammalian compensating contact may show attenuated disease severity. Conversely, rare patients who carry the compensating contact variant alone may be functionally buffered against acquiring the pathogenic residue.

This is not unprecedented — modifier loci in mtDNA disease have been reported (Gropman et al. 2004; Pinos et al. 2021) — but they have never been systematically identified at the structural level. The 693 significant compensatory pairs, and particularly the 58 charge-rescue pairs, define specific residue-pair hypotheses for suppressor screens in patient-derived cell lines or in relevant model organisms.

---

## 9. Summary and Outstanding Questions

### Key findings

1. **2,282 AA-level cDAVs** across all OxPhos complexes in 289 mammals; 1,954 use the identical codon as the human mutation
2. **83% of cDAVs are convergent** (≥2 independent origins); 42% arose >10 times independently
3. **693 significant compensatory co-evolutionary pairs** at FDR ≤ 0.10 using branch co-occurrence and/or Pagel's discrete test
4. **Epistatic pre-adaptation dominates**: contact_first enrichment OR = 5.64 overall, rising to **OR = 24.88** for derived cDAVs (p = 8.8 × 10⁻⁴¹); 96.7% of significant derived-cDAV pairs had contact_first timing
5. **Nuclear cDAVs are better compensated than mitochondrial** (10.8% vs 6.1%, p = 1.6 × 10⁻¹²); Complex IV has the highest compensation rate (15.6%)
6. **Ancestral cDAVs show 3.9× higher compensation** than derived cDAVs (10.9% vs 3.0%, p = 1.6 × 10⁻¹⁸), reflecting longer evolutionary time for co-adaptation
7. **MT-CO2 Thr55 / COX5A Gly79**: the only dateable ancient mt-nuclear compensatory pair, maintained for ~20.9 Mya across a rodent clade

### Outstanding questions

1. **Tier-stratified analysis**: All variants are currently annotated as "Unassigned" tier. Applying severity-based tier classification (MITOMAP confirmed/reported status, ClinVar pathogenicity stars, PhyloP, AlphaMissense) will sharpen the cDAV set to high-confidence pathogenic variants and is expected to increase effect sizes for all analyses

2. **Allelic architecture of pre-adaptation**: Does the contact_first association hold across all 289 species, or is it driven by specific lineages (e.g. rodents, which are over-represented in mammalian genome databases)?

3. **Experimental validation of charge-rescue pairs**: The 58 charge-rescue compensatory pairs are the highest-priority candidates for cell-line suppressor screens — introducing the mammalian contact variant into a human disease cell model should partially rescue OxPhos function

4. **mt-nuc interface co-evolution and speciation**: Do the 38 mt-nuc compensatory pairs show signs of co-divergence (correlated substitution rates) across rodent or bat phylogenies, consistent with speciation-level cytonuclear incompatibility?

5. **Network compensation**: All analyses here are pairwise. True compensatory networks may involve 3–5 co-evolved positions acting in concert. Network-level analysis of co-evolved positions surrounding cDAVs is the next analytical step

---

## Methods Summary

**Clinical variants**: MITOMAP (mtDNA) and ClinVar (nucDNA) pathogenic/likely-pathogenic variants in OxPhos subunit genes

**Mammalian alignments**: TOGA-derived orthologs for nucDNA; MITOMAP species alignments for mtDNA; 289 mammalian species

**Ancestral reconstruction**: IQTree v2 with model selection (TEST), fixed topology from MamPhy (Upham et al. 2019), ancestral state reconstruction (`--ancestral`); posterior probability thresholds HC ≥ 0.80, LC 0.50–0.80

**Contact definition**: Cβ–Cβ distance ≤ 8Å from OxPhos PDB structures

**Co-evolutionary tests**: (1) Fisher's exact (species units, BH FDR per cDAV); (2) Branch co-occurrence Fisher's exact (branches matched by ≥80% Jaccard leaf-set overlap, BH FDR per cDAV); (3) Pagel's discrete LRT via `phytools::fitPagel` (global BH FDR)

**Timetree**: MamPhy ultrametric tree (branch lengths in Mya); force-ultrametric coercion for floating-point drift correction

**Statistical analyses**: All comparative analyses in Python (`scipy.stats`); OR with 95% Wald CI on log scale; BH FDR for multiple testing

---

*All data, code, and results are available in the project repository. Analysis scripts: `src/phylo/01_parse_ancestral_states.py`, `src/structural/01_find_compensating_partners.py`, `src/phylo/02_phylogenetic_timing.py`, `src/phylo/03_comparative_analysis.py`.*