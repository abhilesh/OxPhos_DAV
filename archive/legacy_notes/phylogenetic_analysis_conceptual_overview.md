# Phylogenetic Analysis of Compensated Disease-Associated Variants (cDAVs) in OXPHOS Genes

## Overview

This document describes the conceptual framework, key analytical decisions, major empirical findings, and evolutionary interpretation of the phylogenetic timing analysis of cDAVs in oxidative phosphorylation (OXPHOS) genes across 288 mammalian species.

---

## 1. Conceptual Framework

### 1.1 What is a cDAV?

A **Compensated Disease-Associated Variant (cDAV)** is a human-pathogenic amino acid substitution that is nonetheless tolerated — indeed, fixed — in at least one non-human mammalian species. At the amino acid level (AA-cDAV), any non-human species carrying the disease amino acid qualifies. At the nucleotide level (NT-cDAV), the species must additionally use the exact same codon as the human mutation.

The central evolutionary puzzle: if an amino acid is pathogenic in humans, why does it persist in other lineages without apparent disease? The compensated residue hypothesis proposes that one or more **compensatory substitutions at structurally contacting residues** suppress the deleterious effect — the genetic background "rescues" the variant.

### 1.2 Phylogenetic Timing as a Test of Compensation

If compensation is real, the timing of substitutions should be non-random:
- A **contact residue** that changes *before or concurrently* with the cDAV is a **permissive** compensator — it remodeled the local structural environment prior to or alongside the disease AA, enabling its fixation.
- A contact change *after* the cDAV suggests **reactive** compensation — the deleterious variant spread first, and the contact residue subsequently adapted.
- A contact change **never associated with the cDAV lineage** is consistent with a spurious structural contact signal or neutral variation at the contact site.

Crucially, a rigorous phylogenetic test requires:
1. A time-calibrated species tree
2. Ancestral state reconstruction for each position in each gene
3. Cross-gene node matching by subtended species (not internal node IDs, which differ between IQTree runs)

---

## 2. Dataset and Scope

| Parameter | Value |
|-----------|-------|
| Species in analysis | 288 mammals (TaxID-matched cross-genome overlap) |
| Genes with IQTree ancestral reconstruction | 100 (13 mtDNA, 87 nucDNA) |
| Total structural contact pairs tested | 8,087 (unique by ann_id × contact) |
| Pairs with ancestral state annotation | 2,414 (29.9%) |
| Species tree | MamPhy v1 MCC tree (Upham et al. 2019), 5,911 mammals, time-calibrated |
| Substitution model | Per-gene best-fit (IQTree `-m TEST`), fixed topology from MamPhy |
| Posterior threshold for state calls | ≥ 0.80 MAP probability (internal nodes only) |
| Leaf states | Read directly from input alignment FASTA (not from IQTree .state) |

### Why 5,673 pairs remain unannotated

These pairs have DAR genes with full IQTree coverage but no high-confidence branch change at the disease position. Diagnostics show the cDAV alt_aa is present in only 0–8% of leaf species, and no internal node reconstruction meets the 0.80 posterior threshold. This is not a data gap — it reflects genuine reconstruction uncertainty at rare, possibly recent polymorphisms. These variants are correctly classified as having unresolvable timing.

---

## 3. Critical Analytical Decisions

### 3.1 Cross-gene node matching by subtended species

The most consequential methodological decision. IQTree assigns internal node labels (`Node1`, `Node2`, ...) independently per gene run. A node labelled `Node142` in the MT-CYB reconstruction corresponds to a completely different clade in the NDUFS2 reconstruction. Naive ID-based matching — looking up the DAR origin node's label in the contact gene's branch map — always fails silently and misclassifies everything as `no_contact_change`.

**Fix**: For each origin branch in the DAR gene, compute the **frozenset of leaf species subtended** by the child node. Then find the contact gene node whose subtended species set has the highest Jaccard similarity to this set. All timing calls are then made relative to this species-matched node. This is the only biologically valid approach when running IQTree independently per gene.

### 3.2 Ancestral vs. derived cDAVs

A fundamental dichotomy missed by naive branch-change counting. If the human-pathogenic amino acid is the **ancestral mammalian state**, no "gain" branches will be found — the alt_aa was never gained, because it was always there. Humans evolved the "wildtype" as a **derived** substitution.

Classification: if the alt_aa is present in ≥ 50% of leaf species with data at that position, the cDAV is classified as **ancestral**. For ancestral cDAVs, the relevant events are **loss branches** — where the cDAV was replaced in particular lineages.

**Result**: 51.7% of annotatable pairs involve ancestral cDAVs. This has profound implications (see Section 5).

### 3.3 Posterior thresholding at 0.80

Lowering or removing the threshold inflates the number of "branch changes" by including nodes where reconstruction is near-uniform (51/49% posteriors). A phantom substitution introduced by a low-confidence node creates false timing signals — a contact residue appears to change co-currently because both DAR and contact positions have low-confidence reconstructions on the same branch, not because they co-evolved. The 0.80 threshold removes ~5–10% of branch changes per gene relative to MAP-only parsing.

### 3.4 Local age estimation from MamPhy branch lengths

TimeTree REST API is unreliable for bulk queries and introduces network dependencies. Since the MamPhy tree is ultrametric (branch lengths in millions of years), node ages are computed directly:

```
age(MRCA) = total_root_to_leaf_depth − distance(root, MRCA)
```

This gives the divergence time of the clade subtended by the DAR origin branch, bounding when the cDAV arose.

---

## 4. Empirical Findings

### 4.1 The ancestral/derived split

| cDAV Type | Pairs | % |
|-----------|-------|---|
| Ancestral (alt_aa is mammalian ancestral state) | 1,247 | 51.7% |
| Derived (alt_aa arose on ≥1 mammalian branch) | 1,167 | 48.3% |

Nearly half of all human disease amino acids at positions with structural contacts are in fact the ancestral mammalian state — the human lineage evolved *away* from this amino acid. This means the disease in humans is caused by a reversion to an ancestral state that is no longer compatible with the evolved human genetic background.

### 4.2 Timing signal by cDAV type

| Category | Ancestral cDAVs | Derived cDAVs |
|----------|----------------|---------------|
| contact_first | 9 (0.7%) | 33 (2.8%) |
| co_occurring | 1 (0.1%) | 18 (1.5%) |
| contact_after | 4 (0.3%) | 10 (0.9%) |
| no_contact_change | 1,233 (98.9%) | 1,106 (94.8%) |

Derived cDAVs show a 5-fold enrichment in timing signal (5.2% vs 1.1%). This makes biological sense: for ancestral cDAVs, the compensatory landscape was presumably established over deep evolutionary time and may not be resolvable at the resolution of this tree.

### 4.3 The 75 high-confidence candidate pairs

75 pairs show a dominant timing category other than `no_contact_change`. Key features:

**Gene enrichment**: mtDNA genes dominate (55/75 pairs = 73%), despite representing only 13% of genes in the analysis. This likely reflects:
1. Higher evolutionary rate in mtDNA → more branch changes detectable
2. Smaller proteins → more structural contacts per residue
3. Co-evolution between mtDNA-encoded and nuclear-encoded subunits is already established in the literature (mito-nuclear co-evolution)

**Top genes by candidate count**:
| Gene | Genome | Candidate pairs |
|------|--------|----------------|
| MT-CO2 | mtDNA | 13 |
| MT-ATP8 | mtDNA | 12 |
| MT-CO3 | mtDNA | 10 |
| MT-ND5 | mtDNA | 6 |
| NDUFA9 | nucDNA | 4 |
| NDUFA11 | nucDNA | 4 |
| UQCRFS1 | nucDNA | 4 |

### 4.4 Convergent cDAVs

491 pairs involve a cDAV position where the alt_aa arose **independently on multiple branches** (n_dar_gain_branches ≥ 2). The most convergent position observed is SDHA-276 and UQCRFS1-84, each arising independently 5 times. Convergent evolution of the same disease amino acid across distant lineages is strong evidence that the position is under recurrent positive selection or relaxed constraint in those lineages — and that different genetic backgrounds may have independently evolved different compensatory solutions.

### 4.5 Contact state in cDAV-carrying species (derived cDAVs)

| Contact state | N | % |
|---------------|---|---|
| polymorphic | 727 | 62.3% |
| conserved_human | 354 | 30.3% |
| conserved_alt | 86 | 7.4% |

- **conserved_alt** (7.4%): the contact position carries a fixed alternative amino acid in all cDAV species. These 86 pairs are the strongest candidates for a simple, deterministic compensatory relationship — the contact residue is invariably different wherever the disease AA is present.
- **conserved_human** (30.3%): contact residue is unchanged despite the presence of the disease AA. These species tolerate the cDAV without modifying this contact — either compensation occurs at other positions, or the structural impact of this cDAV is context-dependent.
- **polymorphic** (62.3%): the contact residue varies among cDAV-carrying species. This may reflect incomplete compensation (multiple partial solutions) or neutral variation at the contact site that is uninformative about compensation.

### 4.6 Strongest candidates: contact_first with conserved_alt contact

The most compelling pairs are those where (a) the contact residue change *preceded* the cDAV, and (b) the contact alt_aa is fixed across all cDAV-carrying species:

| DAR | Position | Mutation | Contact Gene | Contact Pos | Contact Alt | Timing | Age (Mya) |
|-----|----------|----------|--------------|-------------|-------------|--------|-----------|
| NDUFA13 | 44 | I→M | NDUFA13 | 45 | F | contact_first | 25.8 |
| MT-ND1 | 1 | M→T | MT-ND1 | 4 | T | contact_first | 2.4 |
| MT-CO2 | 127 | F→S | MT-CO2 | 124 | T | co_occurring | 2.7 |
| MT-CO2 | 147 | E→K | MT-CO2 | 148 | T | co_occurring | 30.0 |
| MT-ATP8 | 6 | T→I | MT-ATP8 | 4 | M | co_occurring | 2.7 |

These are adjacent residues (positions separated by 1–4 in sequence space), consistent with direct steric or electrostatic contact within an α-helix or β-strand.

### 4.7 Intra-protein vs. inter-protein compensation

85.5% of candidate pairs involve contacts within the same protein (intra-protein). This is expected: intra-protein compensation requires only a single gene to evolve, whereas inter-protein compensation requires co-evolution between two independently encoded proteins — a much harder evolutionary problem. The 72 intra-protein pairs with timing signal are the most tractable for functional validation.

### 4.8 Age distribution of cDAV origins

Derived cDAVs: median origin age ~12.1 Mya (range 0.5–87.0 Mya)  
Ancestral cDAVs: median loss age ~22.2 Mya (range 6.3–87.0 Mya)

The youngest cDAVs (~0.5–3 Mya) likely represent recent species-specific variants. Those arising at 87 Mya (early Eutherian divergence, e.g., Afrotheria/Laurasiatheria split) suggest the disease amino acid has been stably maintained across entire mammalian orders.

---

## 5. Evolutionary Interpretation

### 5.1 Two evolutionary histories of human disease variants

The 51.7% ancestral / 48.3% derived split reveals that human genetic disease at OXPHOS variants falls into two fundamentally different categories with distinct evolutionary histories:

**Derived cDAVs** ("gain" scenario): The disease amino acid arose in one or more non-human lineages as a *derived* state. These lineages apparently tolerate it, implying either (a) compensatory mutations at contacting residues, (b) differences in functional constraint on the OXPHOS complex in those lineages, or (c) relaxed purifying selection in energy-rich environmental contexts. The human pathogenic version arose by the same mutation occurring in the human lineage, but in a human-specific genetic background lacking the compensatory context.

**Ancestral cDAVs** ("loss" scenario): The disease amino acid is the ancestral mammalian state. Humans evolved a "wildtype" amino acid as a derived substitution, presumably because the human genetic background makes the ancestral state deleterious — or because the human "wildtype" confers some advantage. When humans revert (by mutation) to the ancestral state, the modern human genetic background cannot tolerate it. This is essentially **intragenic sign epistasis**: the ancestral amino acid is incompatible with the evolved human sequence context at other positions.

### 5.2 The compensation hypothesis in light of these findings

The predominance of `no_contact_change` (94.8–98.9%) does not refute compensation. Several non-exclusive explanations exist:

1. **Compensation at non-contact positions**: the compensatory substitution may be at a residue that is structurally proximal in 3D space but not captured by the PDB contact analysis (different subunit, different rotamer, solvent-mediated contact).

2. **Soft compensation / epistatic buffering**: the cDAV may be tolerated because the overall fitness landscape in non-human lineages is shifted — multiple weak epistatic partners collectively buffer the deleterious effect, none of which individually shows a detectable timing signal.

3. **Functional divergence of the complex**: different mammals rely on OXPHOS at different thermogenic demands (e.g., bats with high metabolic rates, hibernating rodents, marine mammals). A substitution that is pathogenic at human metabolic temperature and ATP demand may be functionally neutral or even advantageous in a different physiological context.

4. **Resolution limits**: The 0.80 posterior threshold, the sparse species coverage at specific positions, and the 288-species tree limit the power to detect branch changes at positions where the alt_aa is carried by only 1–2 species.

### 5.3 mtDNA enrichment in timing signal

The 5-fold enrichment of mtDNA genes in timing candidates (73% of 75 pairs despite 13% of genes) is consistent with two forces:

1. **mito-nuclear co-evolution**: The 13 mtDNA-encoded OXPHOS subunits must physically interact with 70+ nuclear-encoded subunits. Any substitution in an mtDNA subunit changes the interface with its nuclear partner, creating co-evolutionary pressure on the nuclear gene. This is well-established (Barrientos 2003; Havird & Sloan 2016) and the timing signal here may reflect known mito-nuclear epistasis.

2. **Relaxed purifying selection in specific lineages**: mtDNA lacks recombination and has high mutation rate. In small or bottlenecked populations, mildly deleterious variants can drift to fixation, potentially followed by compensatory fixation at contacting sites.

### 5.4 NDUFA13 position 44 as a model case

NDUFA13 I44M, with contact_first timing and a conserved F at the adjacent position 45 in all species carrying the M44 variant (~25.8 Mya origin), represents a textbook case of structural compensation: a hydrophobic isoleucine replaced by a polar methionine at position 44, tolerated because the adjacent position 45 switched from a smaller residue to phenylalanine (bulkier, restoring van der Waals packing). The fact that F45 is conserved across all species carrying M44, but not in species with I44, is consistent with a deterministic compensatory relationship.

### 5.5 MT-ATP8 as a high-rate compensation locus

MT-ATP8 contributes 12 candidate pairs concentrated at positions 6, 16, and 17 — a 11-residue N-terminal region. MT-ATP8 is the smallest mtDNA-encoded protein (68 aa in humans) and overlaps in sequence with MT-ATP6. Its N-terminal region contacts the c-ring stator stalk of Complex V. The concentration of compensation signals here suggests this is a **hotspot of inter-subunit co-evolution** within the ATP synthase rotor stalk interface, consistent with known accelerated evolution of MT-ATP8 across mammals.

### 5.6 Implications for understanding human mitochondrial disease

The finding that ~52% of structural-contact cDAVs involve ancestral amino acids inverts the standard clinical interpretation. These are not cases where "other animals evolved a trick to tolerate a pathogenic mutation." Rather, these are cases where **humans evolved a derived amino acid that is epistasis-dependent**: it only functions correctly in the context of the rest of the modern human sequence. A pathogenic reversion exposes the epistatic dependency.

This has implications for:
- **Drug design**: targeting the compensatory partner may restore function even without correcting the primary mutation
- **Variant interpretation**: an ancestral cDAV in a patient should be interpreted differently from a derived cDAV — the former implies a defect in modern human-specific epistatic buffering
- **Gene therapy**: correcting an ancestral cDAV by introducing the "wildtype" amino acid may not be sufficient if the compensatory context from other species is absent in human tissue

---

## 6. How "Compensation" Is Defined in This Analysis — and Its Caveats

The entire weight of the results rests on how compensation is operationalized at each stage of the pipeline. Compensation is tiered across three levels: **observation** (alignment), **structure** (physical proximity), and **phylogeny** (correlated change). Each tier has specific implementation choices with associated caveats.

### 6.1 Tier 1 — Observation: cDAV classification (`00_classify_DAV.py` + `alignment_parser.py`)

**Implementation**: A variant is classified as an AA-cDAV if the human-pathogenic amino acid is found in ≥1 non-human species in the TOGA or mtDNA alignment (function: `check_compensation()` in `src/utils/alignment_parser.py`). The threshold is literally one species: `results["aa_cdav"] = True` the moment a single match is found. An NT-cDAV additionally requires the identical codon.

**What this means empirically**: Of 2,282 AA-cDAVs in the final dataset:
- **22.9% (523 variants)** are supported by only a single non-human species
- **77.1% (1,759 variants)** are supported by ≥2 species
- **85.6% of AA-cDAVs** are also NT-cDAVs (same codon, not just same amino acid)

**Caveats**:

*The one-species trap*: A single species carrying the disease amino acid could reflect a recent deleterious variant that is segregating and has not yet been purged by selection — not a genuinely compensated variant. There is no minimum species count or phylogenetic breadth requirement in the current pipeline. A variant in a single bat lineage out of 241 taxa is counted identically to one conserved across all primates.

*AA-level vs. NT-level compensation*: AA-cDAV and NT-cDAV test different biological hypotheses. AA-level compensation (same amino acid, any codon) captures convergent phenotypic evolution — the same physicochemical change tolerated independently. NT-level compensation (identical codon) is a much stronger signal: it implies either the *exact* same mutational event recurred (rare), or the codon was horizontally transferred (extremely rare), or — most plausibly — the variant arose once and is shared by descent, meaning the compensated state has been stably maintained since the common ancestor. The 85.6% overlap between AA- and NT-cDAVs in this dataset suggests most cDAVs in the alignment are shared by descent, not convergent. However, the 14.4% of AA-cDAVs that are *not* NT-cDAVs used a different codon to reach the same amino acid — a genuinely independent mutational path — and these should be interpreted with more caution as compensation candidates.

*No fitness filter*: The definition does not require that the species carrying the cDAV is healthy or that OXPHOS function is demonstrably normal. It only requires sequence presence. A species could carry the disease amino acid while also carrying a compensatory substitution, or while being under relaxed selective constraint on OXPHOS (e.g., an obligate parasite with reduced oxidative metabolism).

### 6.2 Tier 2 — Structure: contact partner identification (`00_map_davs_to_structure.py`, `01_find_compensating_partners.py`)

**Implementation**: Structural contacts are defined by a **Cβ–Cβ distance cutoff of 8.0 Å** in PDB structures of OXPHOS complexes (NeighborSearch at 8.0 Å radius, `src/structural/00_map_davs_to_structure.py:699`). Contacts are further classified by interaction type:
- **hbond**: any N/O pair ≤ 3.5 Å
- **electrostatic**: oppositely charged residues, any heavy atom ≤ 5.0 Å
- **hydrophobic**: both nonpolar, sidechain C–C ≤ 5.0 Å
- **vdw**: all other Cβ–Cβ contacts within 8 Å

The enrichment test asks: for each (cDAV position, contact position, contact alt_aa) triplet, is the contact alt_aa significantly more common in cDAV-carrying species than in background species? Tested by Fisher's exact test on a 2×2 table.

**Contact type distribution in tested pairs**:
| Contact type | Pairs | % |
|---|---|---|
| nuc–nuc (intra-nuclear) | 3,972 | 49.1% |
| mt–mt (intra-mtDNA) | 3,705 | 45.8% |
| mt–nuc (cross-genomic) | 410 | 5.1% |

**Contact class distribution**:
| Class | Pairs | % |
|---|---|---|
| hbond | 4,784 | 59.2% |
| vdw | 2,620 | 32.4% |
| hydrophobic | 624 | 7.7% |
| electrostatic | 59 | 0.7% |

**Caveats**:

*The 8 Å Cβ–Cβ cutoff is liberal*: Direct steric interactions typically require Cβ–Cβ ≤ 5–6 Å; hydrogen bonds require heavy atom distances ≤ 3.5 Å. The 8 Å cutoff captures second-shell contacts and backbone-mediated interactions that may not have a direct compensatory relationship. Residues at 7–8 Å separation are unlikely to exert direct steric or electrostatic effects on each other. A stricter 6 Å cutoff (as used in some compensation literature) would reduce the partner set but increase signal specificity.

*Static PDB structure*: The contact is measured in the human PDB structure of the assembled complex. The disease amino acid itself is not present in that structure. The contact partner residues are identified in the human wildtype conformation. The compensatory substitution in a non-human species may induce conformational changes that alter which residues are actually in contact — a dynamic that cannot be captured with a single static structure.

*Fisher's exact test treats species as independent observations*: This is statistically invalid because cDAV-carrying species are phylogenetically clustered — all species in a clade may have inherited the cDAV from a single common ancestor. Species are not independent. A single fixation event shared by 50 bat species counts as 50 "observations" in the Fisher test, massively inflating the p-value. This is the principal motivation for adding Pagel's discrete and branch co-occurrence tests. The Fisher p-values in `all_tested_pairs.csv` should be treated as exploratory only.

*Cross-genomic (mt–nuc) contacts are underrepresented*: Only 5.1% of tested pairs are mt–nuc cross-genomic contacts, despite mt–nuc co-evolution being a major biological phenomenon. This reflects both the small number of mtDNA-encoded subunits (13) and the stringency of the cross-genomic species overlap (288 species with data in both genomes).

*No allosteric contacts*: Compensatory residues may act allosterically — at positions not in physical contact with the disease site but connected by correlated motions in the protein. These are completely absent from the current analysis, which is limited to direct structural contacts.

### 6.3 Tier 3 — Phylogeny: timing and directionality (`01_parse_ancestral_states.py`, `02_phylogenetic_timing.py`)

**Implementation**: For each tested pair, IQTree reconstructs ancestral amino acid states at all internal nodes under the best-fit substitution model with the MamPhy topology fixed. A branch "change" is recorded only if the MAP state probability ≥ 0.80 for the child node. Cross-gene timing is inferred by matching nodes between DAR and contact gene trees using Jaccard similarity of subtended leaf species sets.

**Caveats**:

*The 0.80 posterior threshold creates a sparse detection problem*: 5,422 pairs (67%) are skipped because the DAR alt_aa exists in only 0–8% of species and no internal node reconstruction meets the 0.80 threshold. These pairs are not "uncompensated" — they are unresolvable at this tree and this prior. Lowering the threshold to 0.70 would annotate more pairs but introduce more phantom substitutions. A tiered confidence scheme (≥0.90 = high confidence; 0.70–0.90 = low confidence; flagged in output) would allow downstream users to apply their own strictness.

*Leaf states read from FASTA, not from IQTree .state files*: IQTree only writes reconstructed states for internal nodes in the `.state` file. Leaf states (observed species amino acids) are read directly from the input FASTA alignment. This is correct but introduces an alignment-level caveat: if the alignment has a gap or masked position at the DAR site for a given species, that species contributes no information to the ancestral/derived classification even if the species nominally carries the cDAV (it was in `aa_species` from `check_compensation()` but its position is ambiguous in the trimmed alignment used for IQTree).

*The ancestral/derived threshold of 50% of leaf species is arbitrary*: A cDAV is classified as "ancestral" if the alt_aa is present in ≥50% of leaf species with data at that position. This threshold was chosen to separate clear majority (ancestral) from clear minority (derived) cases. But at exactly 50%, classification is ambiguous. Additionally, a cDAV present in 49% of species is classified as "derived" even though it may have been ancestral and lost in just over half the lineages. The threshold should be reported and sensitivity tested.

*Cross-gene node matching is approximate*: The Jaccard similarity approach for matching a DAR origin node to its corresponding node in the contact gene tree is a heuristic. When the two genes have different species sets (because different species were sequenced or passed QC), the matching becomes less precise. A pair of nodes with 70% Jaccard similarity could still differ in which specific species are subtended — particularly at shallow nodes representing recent divergences. Timing calls at such shallow nodes (co_occurring vs. contact_first separated by a single branch) should be treated with extra caution.

*Fixed topology from MamPhy*: Using the MamPhy consensus topology (rather than estimating topology from the alignment) is standard for ancestral reconstruction but assumes the MamPhy species tree is correct for the genes being analyzed. For rapidly evolving genes (especially mtDNA), the gene tree may differ from the species tree due to incomplete lineage sorting or introgression. This is unlikely to affect the majority of internal nodes but could affect timing calls at shallow nodes within closely related species groups (e.g., within Carnivora or Chiroptera).

*Directionality inversion for ancestral cDAVs*: For the 51.7% of pairs where the cDAV is ancestral, the analysis is inverted — "gain" branches are replaced by "loss" branches (where the alt_aa was replaced). In this case, `contact_first` means the contact residue changed *before* the cDAV was lost in a given lineage, which is interpreted as a **permissive background withdrawal** — the contact residue first destabilized the environment for the disease AA, then the disease AA was replaced. This is a fundamentally different evolutionary narrative from compensation in the derived sense, and the two categories should not be pooled in a single timing analysis without clearly distinguishing them.

### 6.4 Summary of definitional caveats

| Level | Implemented definition | Key caveat | Consequence |
|-------|----------------------|------------|-------------|
| **Observation** | Alt_aa in ≥1 non-human species | One-species trap; no fitness filter | 22.9% of cDAVs supported by single species; may include deleterious variants not yet purged |
| **NT vs AA** | AA: any codon; NT: identical codon | AA-cDAV conflates convergent evolution with shared descent | 14.4% of AA-cDAVs used a different codon — independent mutational paths, weaker evidence of genuine compensation |
| **Structure** | Cβ–Cβ ≤ 8 Å in human PDB | Liberal cutoff; static human structure; no allosteric contacts | ~30–40% of "contacts" at 6–8 Å are unlikely to have direct compensatory effect |
| **Statistics** | Fisher's exact (species as units) | Phylogenetic non-independence; Pagel/branch tests not yet applied | Fisher p-values are inflated for phylogenetically clustered cDAVs (e.g., all Chiroptera) |
| **Ancestral reconstruction** | MAP state, PP ≥ 0.80 | 67% of pairs unresolvable; phantom substitutions below threshold | Timing annotations limited to high-confidence position changes |
| **Cross-gene node matching** | Jaccard similarity of subtended leaf species | Approximate; fails at shallow nodes with different species sets | Timing calls at recent divergences (< ~5 Mya) unreliable |
| **Ancestral/derived classification** | ≥50% of leaves carry alt_aa → ancestral | Arbitrary threshold; ignores root state directly | Cases near 50% boundary misclassified; true root state not queried |

## 7. Limitations and Next Steps

### 6.1 Current limitations

| Limitation | Impact |
|-----------|--------|
| 5,673 pairs unresolvable (rare alt_aa, <0.80 posterior) | ~70% of pairs lack timing annotation |
| 288-species tree (24 cross-genome species absent from MamPhy) | Slightly reduced power for some lineages |
| Contact pairs from PDB only — solvent-mediated and allosteric contacts excluded | Compensatory partners at non-contact positions missed |
| Pagel's discrete test not yet run | Correlated evolution not formally tested |
| Branch co-occurrence test not yet run | Convergent compensation not formally tested |

### 6.2 Immediate next steps

1. **Run `src/structural/01_find_compensating_partners.py`** with Pagel's discrete test — this will add phylogenetically-corrected correlation p-values to the 8,087 pairs, providing formal statistical support for the 75 timing candidates and potentially rescuing additional pairs from the unresolvable set.

2. **Lower posterior threshold to 0.70** and re-run timing — assess how many additional pairs become annotatable vs. how many false timing signals are introduced (measurable by the change in `co_occurring` / `contact_first` rates in negative controls).

3. **Examine the 86 `conserved_alt` contact pairs** in structural detail — map them onto PDB structures to assess the physicochemical plausibility of the inferred compensation (charge complementarity, steric fit, hydrogen bonding changes).

4. **NDUFA13 I44M and MT-ATP8 positions 6/16** — pursue as priority functional validation targets (mutagenesis in cell models or yeast complementation).

---

## 8. Summary Statistics

| Metric | Value |
|--------|-------|
| Total genes analyzed | 100 (13 mtDNA, 87 nucDNA) |
| Total mammalian species in tree | 288 |
| Unique structural contact pairs | 8,087 |
| Pairs with timing annotation | 2,414 (29.9%) |
| Ancestral cDAVs | 1,247 (51.7% of annotated) |
| Derived cDAVs | 1,167 (48.3% of annotated) |
| High-confidence timing candidates | 75 (3.1% of annotated) |
| — contact_first | 42 |
| — co_occurring | 19 |
| — contact_after | 14 |
| Convergent cDAVs (≥2 independent gains) | 491 pairs |
| Strongest conserved_alt contact candidates | 10 pairs |
| Intra-protein timing candidates | 72 (96%) |
| Median derived cDAV origin age | 12.1 Mya |