# Compensated Disease-Associated Variants in Oxidative Phosphorylation: Implications

## Background

Human pathogenic variants in oxidative phosphorylation (OxPhos) subunits cause a spectrum of mitochondrial diseases, yet the same amino acid substitutions are naturally tolerated in other mammalian species. These **compensated disease-associated variants (cDAVs)** — defined here at the amino acid level as any non-human species naturally harbouring the human-pathogenic residue — represent evolved solutions to protein dysfunction that human patients cannot access. A stricter nucleotide-level subset (cDAVs using the identical human-pathogenic codon) accounts for 1,954 of these events.

This analysis screened 6,997 clinical variants across 13 mitochondrially-encoded and 87 nuclear-encoded OxPhos genes in 289 mammalian species, identifying **2,282 AA-level cDAVs** (236 mtDNA, 2,046 nucDNA) from 443 mtDNA and 6,554 nucDNA variants. The central question is: what structural and evolutionary features allow these species to tolerate mutations that are pathogenic in humans?

---

## Key Findings

### 1. Compensatory co-evolution is pervasive across the OxPhos system

Of the 8,087 DAR–contact residue pairs tested for co-evolution, **693 pairs across 267 unique cDAV positions** show statistically significant co-occurrence of the cDAV and an altered contact residue (Fisher + branch co-occurrence FDR ≤ 0.1 or Pagel's discrete FDR ≤ 0.1). This represents a substantial fraction of the tested landscape and demonstrates that compensatory co-evolution is not an isolated phenomenon but a recurring feature of OxPhos evolution across mammals.

The breadth of this signal spans all five OxPhos complexes:

| Complex | Compensatory pairs |
|---|---|
| Complex I (NADH dehydrogenase) | 351 |
| Complex IV (Cytochrome c oxidase) | 131 |
| Complex II (Succinate dehydrogenase) | 110 |
| Complex III (bc1 complex) | 34 |
| Complex V (ATP synthase) | 66 |

Complex I dominates, consistent with it being the largest OxPhos complex (45 subunits), the most frequent source of mitochondrial disease, and the site with the highest density of subunit–subunit and subunit–cofactor contacts.

### 2. Two complementary lines of evidence converge on a small high-confidence set

- **172 pairs** pass branch co-occurrence FDR ≤ 0.1: the cDAV and contact substitution arose on the same phylogenetic branch more often than expected by chance, controlling for phylogenetic non-independence
- **25 pairs** pass Pagel's discrete correlated-evolution test (FDR ≤ 0.1): the two characters evolved in a correlated fashion across the tree under a model-based likelihood framework
- **1 pair** is supported by both tests independently

The difference in counts reflects the distinct statistical power profiles: branch co-occurrence is sensitive to recent co-arising events in species-dense clades, while Pagel's test has greater power when co-evolution has accumulated across many independent lineages. The Pagel-only set of 24 pairs therefore likely represents older, more deeply embedded compensatory relationships, while the branch-only set captures more recent or lineage-specific compensations.

### 3. Intra-subunit compensations predominate, but inter-subunit compensations reveal functional interfaces

**594/693 pairs (86%)** involve a contact residue within the same subunit as the cDAV. This is expected — packing contacts within a single polypeptide are the most proximal structural environment for a pathogenic substitution, and local rewiring of side-chain contacts is the most accessible evolutionary path.

However, **99 inter-subunit pairs (14%)** are particularly significant: these identify specific subunit–subunit interfaces at which charge and steric environments have co-evolved to accommodate the pathogenic residue. Given that OxPhos assembly defects are a leading cause of complex instability in mitochondrial disease, these interfaces are candidate loci for studying how disease-associated destabilisation can be buffered by co-evolved partner subunits.

Contact types among significant pairs reflect the structural diversity of compensations:
- Hydrogen bonds (403, 58%) — the most directional, most physically specific contact type
- Van der Waals contacts (241, 35%) — steric packing adjustments
- Hydrophobic contacts (43, 6%) — core packing changes
- Electrostatic contacts (6, <1%) — charge–charge interactions (likely underrepresented given the structural distance cutoff)

### 4. Two physicochemical mechanisms dominate the identifiable compensations

Among the 261 cDAV–contact pairs with full physicochemical annotation:

- **Charge rescue (58 pairs, 22%)**: The contact residue undergoes a charge change that counteracts the charge introduced or removed by the pathogenic DAR substitution. These represent the most mechanistically interpretable compensations — structural stability maintained through electrostatic neutralisation. Notably, 18/58 (31%) of charge rescue pairs involve ancestral cDAVs (the human-pathogenic residue present in the common ancestor), suggesting that charge-rescue co-evolution can be ancient enough to have become fixed across many mammalian lineages.

- **Volume swap (64 pairs, 25%)**: The contact residue undergoes a volume change that mirrors the steric change introduced by the DAR substitution — larger pathogenic residue paired with a smaller contact residue, or vice versa. These represent packing compensation: the cavity or clash introduced by the pathogenic substitution is geometrically offset by the contact partner.

- **Unclassified (139 pairs, 53%)**: Pairs where the physicochemical changes in DAR and contact residue do not fall into a simple biophysical category. These may involve indirect structural transmission (second-shell effects), co-evolved backbone geometry, or context-dependent interactions not captured by pairwise residue properties alone.

### 5. Evolutionary timing reveals pre-adaptation as a mechanism

The dominant timing classification among the 261 annotated pairs distinguishes two biologically distinct scenarios:

**No contact change (190/261, 73%)**: The cDAV arose in species that already carried the human reference residue at the contact position. In these cases, the contact residue was not the source of compensation — the cDAV is tolerated through other structural or functional mechanisms not captured in this pairwise contact analysis (e.g., global protein stability differences, altered protein expression levels, or compensations at non-contacting positions).

**Contact first (71/261, 27%)**: The compensating contact substitution is inferred to have arisen *before* the cDAV on the phylogenetic tree. This is the most biologically significant scenario: the contact variant pre-adapted the structural environment to be permissive of the human-pathogenic residue *before* that residue appeared. This temporal ordering implies genuine epistatic pre-adaptation — an evolutionary path where one substitution opens the fitness landscape for a second substitution that would otherwise be deleterious.

This 27% contact-first rate is likely a lower bound, as the analysis is limited to the single most parsimonious inferred gain branch per cDAV; parallel or convergent origins may inflate the apparent "no contact change" category.

### 6. The ancestral state of cDAVs reflects both ancient and lineage-specific events

Of the 261 annotated pairs:
- **78 (30%)** involve ancestral cDAVs — the human-pathogenic residue was present in the inferred common ancestor of the cDAV-bearing clade (reconstructed by IQTree at the root node of the cDAV subtree)
- **183 (70%)** involve derived cDAVs — the residue arose within a specific mammalian lineage after divergence

Ancestral cDAVs are particularly important: they imply the pathogenic substitution has been stably maintained across potentially tens of millions of years in multiple lineages, and that the compensatory structural environment is itself deeply conserved. These cases provide the strongest evolutionary evidence that the human-pathogenic residue is truly compatible with OxPhos function when the right structural context is present.

Derived cDAVs, by contrast, represent more recent acquisitions that may be lineage-specific adaptations (e.g., metabolic specialisations in hibernating mammals, diving mammals, or high-altitude species) rather than neutral variation.

---

## Comparative Phylogenetic Analyses

### 1. Independent origins: cDAV convergence is the rule, not the exception

Of 903 unique cDAV positions with timing annotations:

| Independent origins | cDAVs | % |
|---|---|---|
| 1 (single lineage) | 150 | 16.6% |
| 2–5 | 244 | 27.0% |
| 6–10 | 129 | 14.3% |
| >10 | 375 | 41.5% |

**83% of cDAVs arose independently in ≥2 lineages.** This extraordinary convergence demonstrates that these human-pathogenic residues are accessible on the mammalian fitness landscape and are repeatedly acquired across independent lineages. The median derived cDAV has 15 independent origins; the median ancestral cDAV has 115 independent loss events.

Convergent cDAVs (≥2 origins) are significantly enriched for compensatory partners relative to single-origin cDAVs (OR = 2.65, 95% CI 1.23–5.71, Fisher p = 0.003). This means that variants repeatedly acquired across mammals are more likely to be found in structural contexts that co-evolve compensating contacts — consistent with the interpretation that convergence selects for robust compensatory backgrounds.

### 2. Temporal ordering: epistatic pre-adaptation is the dominant mode

Contact_first timing (compensating contact arose before the cDAV on the same lineage) is massively enriched in significant co-evolutionary pairs:

- **All pairs**: OR = 5.64 (95% CI 3.98–7.99), p = 1.2 × 10⁻³⁰
- **Derived cDAVs**: OR = 24.88 (95% CI 11.35–54.51), p = 8.8 × 10⁻⁴¹
- **Ancestral cDAVs**: OR = 1.92 (95% CI 1.19–3.10), p = 0.005

**85.8% of all significant compensatory pairs had contact_first timing.** For derived cDAVs specifically, 96.7% of significant pairs had contact_first timing — an odds ratio of nearly 25. This is not statistical noise.

This finding reframes the mechanism of compensation: rather than the cDAV arising and then driving selection for a rescue substitution, the dominant mode appears to be **structural pre-adaptation** — the compensating contact substitution was already present in the lineage before the pathogenic residue appeared. The lineage that already carries a permissive structural background is more likely to tolerate, and therefore to maintain, the human-pathogenic residue.

The much weaker enrichment in ancestral cDAVs (OR = 1.92 vs 24.88 for derived) is expected: ancestral cDAVs were already present at the mammalian root, so the "contact_first" framing is less meaningful — the compensating contact co-evolved with a variant that predates the mammalian radiation.

The contact_first enrichment is consistent across mt-mt (OR = 5.83), nuc-nuc (OR = 5.47), and mt-nuc pairs (OR = 3.14), indicating that pre-adaptation is a universal mode across all genomic contexts.

### 3. Genomic architecture: nucDNA cDAVs are better compensated than mtDNA

| Contact type | Pairs | Significant | Rate |
|---|---|---|---|
| mt-mt | 3,705 | 226 | 6.1% |
| nuc-nuc | 3,972 | 429 | 10.8% |
| mt-nuc | 410 | 38 | 9.3% |

Chi-square across contact types: χ² = 54.33, df = 2, p = 1.6 × 10⁻¹². **Nuclear-encoded cDAVs (nuc-nuc) show nearly twice the compensation rate of mtDNA cDAVs (mt-mt):** 10.8% vs 6.1%. The Fisher p-value distributions also differ significantly (Mann-Whitney U, p = 6.1 × 10⁻²⁹), confirming a real difference in the strength of co-evolutionary signal.

This is consistent with the higher structural complexity and more extensive subunit–subunit contacts in nuclear-encoded OxPhos subunits, which provide more opportunities for compensatory contact evolution. mtDNA-encoded subunits are smaller (averaging ~200 aa vs ~400 aa), evolve faster (reducing signal-to-noise in co-evolution tests), and have fewer contact partners.

Intergenomic (mt-nuc) pairs show a compensation rate (9.3%) intermediate between mt-mt and nuc-nuc, and are not significantly different from intragenomic pairs overall (OR = 1.11, p = 0.59). However, their Fisher p distribution differs from nuc-nuc (p = 2.5 × 10⁻⁸), suggesting qualitatively different co-evolutionary dynamics at the mt-nuclear interface compared to purely nuclear contacts.

### 4. Ancestral cDAVs are compensated at 3.6× the rate of derived cDAVs

| | Pairs | Significant | Rate | Mean origins |
|---|---|---|---|---|
| Ancestral cDAVs | 715 | 78 | 10.9% | 116 |
| Derived cDAVs | 6,023 | 183 | 3.0% | 47 |

OR = 3.92 (95% CI 2.97–5.17), p = 1.6 × 10⁻¹⁸. **Ancestral cDAVs — where the human-pathogenic residue was present in the common mammalian ancestor — are nearly four times more likely to show significant compensatory co-evolution than derived cDAVs.** This is expected: ancestral cDAVs have had the entire mammalian radiation (~180 Mya) to accumulate compensatory contacts, while derived cDAVs are younger lineage-specific events.

The dramatic difference in mean independent origins (116 vs 47) reflects the fundamental difference in these categories: ancestral cDAVs are counted by loss events (species that reverted to the human reference), while derived cDAVs are counted by gain events. The high loss counts for ancestral cDAVs reflect widespread reversion to the "human" state across the mammalian tree.

---

## Broader Implications

### For understanding mitochondrial disease pathogenesis

The existence of 2,282 cDAVs across 289 mammalian species demonstrates that pathogenicity of OxPhos variants is not an intrinsic property of the substitution itself, but depends on the structural and genetic context in which it occurs. A variant that disrupts human Complex I assembly may be compatible with function in a rodent where adjacent contact residues have co-evolved to accommodate it.

This has direct implications for variant interpretation: pathogenicity prediction tools (REVEL, AlphaMissense, ESM-1b) trained on conservation metrics will systematically overpredict pathogenicity for positions where non-human mammals naturally carry the human disease allele. The 2,282 cDAV positions identified here should be treated with caution when using such tools.

### For compensatory mutation therapy

The 693 significant cDAV–contact pairs define a specific map of residue pairs where co-evolution has demonstrably compensated for the pathogenic substitution. The 58 charge-rescue pairs in particular provide testable hypotheses for intragenic suppressor mutations: introducing the mammalian contact variant into a human cell model of the corresponding mitochondrial disease could partially or fully restore protein function. This class of suppressor screen is currently limited by the lack of systematic identification of candidate compensatory positions — the present analysis provides exactly this catalogue for OxPhos.

### For understanding epistasis in protein evolution

The 71 contact-first pairs represent direct evidence of historical epistasis: one substitution rendering a second substitution tolerable. The physicochemical breakdown of these pairs — predominantly unclassified (56%) with significant volume swap (20%) and charge rescue (24%) contributions — suggests that pre-adaptation operates through multiple structural mechanisms. The preponderance of derived cDAVs in contact-first pairs (64/71, 90%) is consistent with recent co-evolution: the compensating contact arose in a specific lineage, and that same lineage subsequently acquired the human-pathogenic residue.

### Evolutionary age of compensatory events

Phylogenetic timing reveals that compensatory co-evolution is predominantly a recent, lineage-specific phenomenon. Of the 246 annotated partner pairs with age estimates:

- **242/246 (98%) are leaf-origin events** — the cDAV arose on a terminal branch in a single extant species, with no detectable shared ancestral origin across multiple lineages. These compensations represent private, species-specific adaptations that have not spread to related taxa.
- **4/246 (2%) have internal-node origins**, all clustering at **~20.9 Mya** at a single node (Node229) in the rodent phylogeny. These 4 pairs all involve the same position — **MT-CO2 Thr55** (human disease allele: Thr→X) — and include three intra-subunit contacts and one mt–nuclear contact with COX5A. The shared origin implies this cDAV arose once in the common ancestor of a rodent clade ~20.9 Mya and has been stably maintained alongside its compensating contact partners across that entire lineage. This is a rare example of ancient, fixed compensatory co-evolution within OxPhos.

The predominance of single-species (leaf-origin) events has an important implication: most mammalian cDAVs have not been "tested" by evolution across multiple independent genetic backgrounds. The stability of the compensated state is therefore known in only one genomic context, limiting the generalisability of derived suppressor hypotheses. The MT-CO2 55T case — maintained for ~21 million years across a rodent clade — provides substantially stronger evidence of robust compensation.

### For mt–nuclear co-evolution

**38 of the 693 compensatory pairs are mt–nuclear contact pairs** — cases where a mitochondrially-encoded cDAV is compensated by a nuclear-encoded contact partner, or vice versa. These are particularly significant because mt and nuclear genomes are inherited through different mechanisms (maternal vs. Mendelian), and incompatibilities between them underlie cytonuclear co-evolution, speciation barriers, and population-level fitness variation. Each of these 38 pairs represents a co-evolved mt–nuclear residue contact where the two genomes have jointly adapted to a human-pathogenic residue — a direct window into the molecular basis of cytonuclear epistasis.

The MT-CO2 Thr55 case is notable in this context: one of its four compensatory contacts (COX5A Gly79) is nuclear-encoded. The 20.9 Mya co-maintenance of this mt–nuclear pair implies that the two genomes have been co-evolving at this interface across the entire rodent clade, representing one of the best-dated examples of mt–nuclear co-adaptation at atomic resolution.

---

## Limitations and Caveats

1. **Tier assignment is currently Unassigned for all variants**: The severity-based tier classification (using MITOMAP, ClinVar, gnomAD, MitImpact, dbNSFP scores) has not yet been applied, meaning the analyses do not currently distinguish high-confidence pathogenic variants from those with uncertain clinical significance. Filtering to Tier 1–2 variants will substantially sharpen the cDAV set.

2. **Pagel's test coverage**: Only 769 of 3,736 Fisher-filtered pairs produced finite Pagel p-values. The majority of NA values reflect near-monomorphic trait distributions (cDAV present in very few species), which genuinely lack the statistical power for model-based correlated evolution tests. This is a biological constraint, not a methodological failure.

3. **Contact definition**: Compensatory contacts are defined by Cβ–Cβ distance ≤ 8Å in available PDB structures. This misses contacts transmitted through backbone geometry, water-mediated interactions, and allosteric propagation. The 53% unclassified physicochemical category partly reflects this limitation.

4. **Single-pair analysis**: The analysis tests each DAR–contact pair independently. True compensatory networks may involve multiple co-evolved positions acting in concert; pairwise tests will underestimate the full compensation.

5. **Inferred ancestral states**: IQTree posterior probability reconstructions at internal nodes carry uncertainty, particularly for deep nodes and positions under strong selection. High-confidence threshold (PP ≥ 0.80) is applied throughout, but misclassifications will occur.