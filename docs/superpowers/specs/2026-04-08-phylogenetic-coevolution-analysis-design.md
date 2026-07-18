# Design: Phylogenetic Co-evolution Analysis of cDAVs and Compensatory Partners

**Date:** 2026-04-08  
**Status:** Approved  
**Context:** OxPhos_DAV project — replaces Fisher's exact test in `01_find_compensating_partners.py` with phylogenetically-aware methods and adds branch-level timing annotation.

---

## Problem

The current compensatory partner analysis (`src/structural/01_find_compensating_partners.py`) uses Fisher's exact test to identify contact residues enriched in cDAV-carrying species. This is statistically invalid: cDAV-carrying species are phylogenetically clustered (e.g., all Chiroptera inherited the same disease AA once from a single ancestor), so they are **not independent observations**. The test inflates signal for traits that are clade-specific.

Additionally, the output `compensatory_partners.csv` is filtered to only Fisher-significant pairs (FDR ≤ 0.10) before writing, creating a **circular dependency**: the timing analysis would annotate pairs selected by the invalid test.

---

## Solution Overview

1. **Pagel's discrete model** (R `ape`/`phytools`): formal likelihood ratio test for correlated binary character evolution. Fits independent vs. dependent Markov transition models for two binary traits (cDAV presence, contact alt AA presence) on the phylogenetic tree.

2. **Branch co-occurrence enrichment** (IQTree `--ancestral`): infer ancestral states at all internal nodes, then count *branches* (not species) where the DAR AA and contact AA co-changed. Fisher's exact over branches is valid because branches are approximately independent evolutionary events.

3. **Phylogenetic timing annotation**: for each pair, determine on which mammalian branch(es) the cDAV arose, and whether the compensatory contact substitution occurred before, concurrently, or after in each lineage that carries it. Cross-validate with TimeTree divergence times.

4. **Decouple output from significance**: the partners script now writes all tested pairs to `all_tested_pairs.csv` with all three test columns. `compensatory_partners.csv` becomes a derived view filtered by the phylogenetically valid tests.

---

## Species Tree

**Primary**: VertLife (Upham et al. 2019) — 5,911 mammal species, fully ultrametric with calibrated divergence times. Download `MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.nwk` from vertlife.org/phylosubsets/ and place at `data/phylo/species_tree/mammaltree.nwk`.

**Cross-validation**: TimeTree REST API for divergence time estimates at key nodes identified during timing analysis.

**Per-test pruning**: the VertLife tree is pruned per test to the species intersection of the DAR gene alignment and contact gene alignment. Different gene pairs → different N. Tests with N < 20 are flagged `low_power=True`.

### Species scope by pair type

| Pair type | Species set |
|-----------|-------------|
| nuc–nuc intra-genomic | TOGA spp in DAR gene alignment ∩ contact gene alignment |
| mt–mt intra-genomic | mtDNA spp in both alignments |
| mt–nuc cross-genomic | Both alignments ∩ 405-species cross-genome overlap |

The 405-species cross-genome overlap comes from `data/reference/taxid_species_mapping.csv` (`Exact_TaxID_Match` rows).

---

## IQTree Step (HPC)

IQTree ancestral reconstruction is computationally expensive and runs on HPC (SLURM + conda IQTree). Scope: genes appearing in `all_tested_pairs.csv` only (DAR gene + contact gene per tested pair — approximately 20–40 genes).

**Command per gene:**
```bash
iqtree2 -s {gene}.fasta -te {gene}_tree.nwk -m TEST --ancestral \
        --prefix {gene} -nt 4 --quiet
```

---

## New Files

```
src/phylo/
  00_prep_iqtree_jobs.py       # Prune VertLife tree per gene, write HPC inputs
  01_parse_ancestral_states.py # Parse IQTree .state files → ancestral_state_maps.json
  02_phylogenetic_timing.py    # Annotate all pairs with timing + TimeTree ages
  pagel_discrete.R             # R helper: Pagel's discrete LRT

scripts/
  slurm_iqtree_array.sh        # SLURM job array template

data/phylo/
  species_tree/                # VertLife Newick
  iqtree_jobs/{gene}/          # Per-gene FASTA + pruned tree (transfer to HPC)
  ancestral_states/{gene}/     # IQTree .state + .treefile (copy back from HPC)
  cross_genome_species.txt     # 405-species overlap list
  ancestral_state_maps.json    # Parsed ancestral states (branch-level changes)

results/phylo/
  timing_annotations.csv       # Per-pair timing + TimeTree age bounds
```

---

## Modified Files

- `src/structural/01_find_compensating_partners.py` — output all tested pairs; add `pagel_p`/`pagel_fdr`, `branch_cooccur_p`/`branch_cooccur_fdr`, `n_species_in_test`, `low_power` columns; derive `compensatory_partners.csv` as filtered view
- `Dockerfile` — add `r-base`, `r-cran-ape`, `r-cran-phytools`

---

## Output Schema

### `results/structural/all_tested_pairs.csv` (new primary output)

All columns from the existing `compensatory_partners.csv`, plus:

| New field | Description |
|-----------|-------------|
| `fisher_p` | Renamed from `p_value` |
| `fisher_fdr` | Renamed from `fdr` |
| `n_species_in_test` | Species in the intersection for this pair |
| `low_power` | True if N < 20 |
| `pagel_p` | Pagel's discrete LRT p-value |
| `pagel_fdr` | BH-corrected across all pairs |
| `branch_cooccur_p` | Fisher p-value over branches |
| `branch_cooccur_fdr` | BH-corrected |
| `n_cooccur_branches` | Branches where both DAR and contact changed |
| `n_dar_only_branches` | Branches where only DAR changed |
| `n_contact_only_branches` | Branches where only contact changed |

### `results/structural/compensatory_partners.csv` (derived view)

Filtered to: `(pagel_fdr ≤ 0.10 OR branch_cooccur_fdr ≤ 0.10) AND low_power == False`.  
Falls back to Fisher-filtered if phylogenetic tests not yet run (graceful degradation).

### `results/phylo/timing_annotations.csv`

| Field | Description |
|-------|-------------|
| `ann_id`, `dar_gene`, `contact_gene`, `contact_alt_aa` | Pair identity |
| `n_dar_origin_branches` | Independent origins of the cDAV on the tree |
| `n_contact_first` | Lineages where contact changed before cDAV |
| `n_co_occurring` | Lineages where both changed on same branch |
| `n_contact_after` | Lineages where contact changed after cDAV |
| `n_no_contact_change` | cDAV lineages with no contact change |
| `dominant_timing` | Most common timing category |
| `dar_origin_node` | Tree node label of most common cDAV origin |
| `timetree_age_mya_min` | TimeTree lower bound divergence time (Mya) |
| `timetree_age_mya_max` | TimeTree upper bound divergence time (Mya) |

---

## Pipeline Order

1. (Docker) Run `src/structural/01_find_compensating_partners.py` → writes `all_tested_pairs.csv` with Pagel + branch columns (None until IQTree done), `compensatory_partners.csv` (Fisher fallback)
2. (Docker) Run `src/phylo/00_prep_iqtree_jobs.py` → `data/phylo/iqtree_jobs/` + `slurm_iqtree_array.sh`
3. (HPC) `rsync` jobs → HPC; `sbatch slurm_iqtree_array.sh`; `rsync` results back
4. (Docker) Run `src/phylo/01_parse_ancestral_states.py` → `ancestral_state_maps.json`
5. (Docker) Re-run `src/structural/01_find_compensating_partners.py` → fills branch co-occurrence columns; re-derives `compensatory_partners.csv`
6. (Docker) Run `src/phylo/02_phylogenetic_timing.py` → `timing_annotations.csv`

---

## Interpretation Notes

- Each gene pair has a different N (species in test intersection). Report `n_species_in_test` in all outputs to allow power-stratified interpretation.
- For mt–nuc cross-genomic pairs, N is bounded by the 405-species overlap. These tests will have lower power than intra-genomic pairs.
- Pagel's discrete assumes a continuous-time Markov model for each binary trait. It is sensitive to very rare traits (< 3 state changes); `low_power=True` (N < 20) partially guards against this but does not catch all edge cases.
- Multiple independent origins (high `n_dar_origin_branches`) strengthen the case for genuine parallel evolution and increase power for the branch co-occurrence test.
