# Final Analysis Filtering Decisions

This document tracks issues discovered during pipeline development that must be
accounted for when constructing the final filtered analysis views, manuscript
tables, and statistical tests. It is distinct from `issues_decisions_resolutions.md`
(which records resolved implementation decisions) — entries here are prospective
flags that must survive into the analysis layer.

Each entry states the issue, the correct filter or caveat to apply, and which
output files or analyses it affects.

---

## F1 — gnomAD / dbNSFP completeness not guaranteed

**Source:** IDR D0-2  
**Issue:** MyVariant batch queries used `size=1000` with no check for truncation. Genes
with >1000 ClinVar entries may have incomplete gnomAD/dbNSFP annotation rows.  
**Filter:** Do not use gnomAD allele frequency or dbNSFP scores as primary filters or
claims without first re-acquiring annotations per-variant via variant-ID-based queries.
Population frequency comparisons are currently unsupported as primary evidence.  
**Affects:** Any analysis stratifying variants by gnomAD AF or using dbNSFP functional
scores as primary evidence.

---

## F2 — ClinVar mtDNA variants are duplicate / ineligible rows

**Source:** Discovered during IQTree scope audit, 2026-05-14  
**Issue:** mtDNA variants in ClinVar appear twice in the download (once per genome
assembly: GRCh37 and GRCh38) because the mitochondrial reference (rCRS, NC_012920.1)
is unchanged between assemblies. All such rows are `skipped_by_policy`:
- GRCh37 copy → `non_grch38`
- GRCh38 copy → `mt_in_clinvar_nuclear_branch`

This is correct behaviour. MITOMAP is the authoritative mtDNA source. Affected genes
confirmed as 100% skipped_by_policy: MT-CO1 (476), MT-CO2 (260), MT-CO3 (340),
MT-CYB (666). These are not missing data — they are correctly excluded ClinVar
duplicates of MITOMAP-sourced variants.  
**Filter:** When reporting variant counts per gene or complex, use
`source == "MITOMAP"` for mtDNA genes and `source == "ClinVar"` for nucDNA genes.
Do not count ClinVar mtDNA rows in any denominators.  
**Affects:** Per-gene and per-complex variant count tables; supplementary tables
reporting source breakdown.

---

## F3 — 12 genes absent from the structural panel (accepted coverage gap)

**Source:** `notes/structural_mapping_eligibility_gaps.md` Problem 3  
**Issue:** The following genes have classified variants but no cryo-EM chain coverage
in any active structure. They are absent from all contact analyses.

| Gene | Complex | Classified variants | Reason |
|---|---|---|---|
| ATP5ME | CV | 18 | Peripheral CV subunit, not resolved |
| ATP5MG | CV | 18 | Peripheral CV subunit |
| COXFA4L2 | CIV | 17 | Tissue-specific CIV regulatory subunit |
| COX8C | CIV | 15 | Tissue-specific CIV isoform |
| ATP5IF1 | CV | 14 | CV inhibitory factor |
| COX6B2 | CIV | 12 | Testis-specific CIV isoform |
| COX7B2 | CIV | 12 | Testis-specific CIV isoform |
| ATP5MF | CV | 11 | Peripheral CV subunit |
| ATP5MK | CV | 9 | Peripheral CV subunit |
| ATP5MJ | CV | 3 | Peripheral CV subunit |
| COXFA4L3 | CIV | 2 | Tissue-specific CIV regulatory subunit |

**Filter:** Exclude these genes from all structural and compensatory partner analyses.
Disclose as "outside the structural panel" in the manuscript methods. Do not include
in mapping rate denominators.  
**Affects:** Structural mapping summaries; compensating partners; mutagenesis
prioritization; complex-level coverage tables.

---

## F4 — 347 eligible-but-failed structural mapping variants

**Source:** `docs/analysis_plan.md` Stage 3; `notes/structural_mapping_eligibility_gaps.md`  
**Issue:** 347 variants were mapping-eligible but failed:
- 505 `mature_offset_candidate` — transit-peptide offset not resolved in the exception
  registry (likely the true count is lower after deduplication; exact number to be
  confirmed on next mapping run)
- 154 `residue_anchoring_failure` — position absent from PDB chain

These are not bugs — they are genuine coordinate resolution failures. They are absent
from `dar_contacts_cbcb8A.csv` and therefore from all compensating partner analyses.  
**Filter:** Report separately as "mapping-eligible but structurally unresolved" in
the methods. Do not include in the mapped-variant denominator. If transit-peptide
offset rescue is extended, these can be partially recovered in a future run.  
**Affects:** Mapping rate reporting; per-complex coverage tables.

---

## F5 — Single-origin cDAVs inflate Fisher and branch co-occurrence

**Source:** Compensating partners design; `docs/analysis_plan.md` Stage 4  
**Issue:** A cDAV that arose once in a mammalian ancestor will be shared across an
entire clade. Fisher's exact (species as units) treats all clade members as independent
observations, inflating the test statistic. Branch co-occurrence is also inflated for
single-origin cDAVs because there is only one DAR branch and any contact change in
that clade appears co-occurring by definition.  
**Filter:**  
- Fisher p-values are **not valid** as primary evidence. Retained for comparison only.  
- For branch co-occurrence: stratify results by number of independent DAR origins
  (`n_dar_branches`). Single-origin pairs (n=1) require pyvolve conditional permutation
  as the primary arbiter before claiming co-evolution.  
- Pagel's discrete test is the preferred primary test for single-origin cDAVs (it
  models correlated evolution on the tree and is less sensitive to clade inflation).  
**Affects:** `compensatory_partners.csv` interpretation; all co-evolutionary enrichment
claims; manuscript figures showing OR or enrichment.

---

## F6 — Pagel MT-MT coverage bias (~50–70% NA rate)

**Source:** Compensating partners design; observed in current run  
**Issue:** mtDNA-encoded subunits are highly conserved across mammals. After tree
pruning to species where both the DAR and contact positions are readable, the contact
trait is often invariant (all 0 or all 1), causing `fitPagel` to return NA. MT-MT
pairs have a ~50–70% NA rate for Pagel p-values.  
**Filter:** For MT-MT pairs, do not use Pagel as the primary test. Use multi-origin
branch co-occurrence and pyvolve permutation instead. Cross-genome Pagel comparisons
(nucDNA vs mtDNA) are biased and should not be made without accounting for this NA
differential.  
**Affects:** Complex-level Pagel enrichment comparisons; any claim that MT-MT pairs
are less enriched than NUC-NUC pairs based on Pagel alone.

---

## F7 — contact_first odds ratio may be inflated (pyvolve required)

**Source:** Compensating partners output; plan file  
**Issue:** The observed contact_first OR of ~24.88 in the current run is high enough
to suspect inflation from single-origin cDAVs (F5 above). A contact alt AA that arose
before a single-origin DAR will always appear as contact_first regardless of whether
the two are causally linked.  
**Filter:** The contact_first OR must be validated against the pyvolve conditional
permutation null (simulate contact-column evolution on the gene tree without coupling
to the DAR; empirical p-value for whether observed OR exceeds the null). Do not report
contact_first as evidence of pre-adaptation until pyvolve permutation is complete
(`src/phylo/04_conditional_permissiveness.py`).  
**Affects:** Any claim of pre-adaptive compensation; Figure panels showing
contact_first enrichment.

---

## F8 — Tier assignment is deferred (all variants currently "Unassigned")

**Source:** IDR; pipeline design  
**Issue:** All variants in the current parquet have `tier = "Unassigned"`. Tier
assignment (Tier A: NT-level cDAV, multi-species; Tier B: AA-level, multi-species;
Tier C: AA-level, single-species) is deferred to the final analysis layer per the
filter-late contract.  
**Filter:** Apply tier assignment as a derived view at the final analysis step.
Do not use `tier` as a filter criterion in any intermediate script. When reporting
results by tier, derive tiers from `is_cdav_amino_acid`, `is_cdav_nucleotide`, and
`n_lineages_with_disease_allele` (or equivalent species count column) at analysis
time.  
**Affects:** All stratified results tables; Tier A/B/C breakdowns in manuscript.

---

## F9 — COXFA4 / NDUFA4 naming across pipeline layers

**Source:** `notes/structural_mapping_eligibility_gaps.md`; confirmed 2026-05-14  
**Issue:** This gene has three different names in use across pipeline layers:
- **Stage 0–2 (parquet):** `interpreted_gene = "COXFA4"` — correct per the frozen HGNC
  source (`Canonical_OXPHOS_Subunits_HGNC_2026-03-25.csv`, line 32), where COXFA4 is
  the approved symbol and NDUFA4 is a previous symbol. The parquet must not be patched.
- **Stage 3 structural outputs:** `dar_locus = "NDUFA4"` — the structural mapper
  remaps COXFA4→NDUFA4 internally for PDB/UniProt/chain lookups and writes the
  remapped name to the contacts CSV.
- **TOGA alignment filename:** `COXFA4_aa_alignment.fasta` — matches parquet.

**Resolution (Stage 4 onwards):** `01_find_compensating_partners.py` uses
`_ALIAS_EQUIV = {"COXFA4": "NDUFA4", "NDUFA4": "COXFA4"}`. `load_alignments` stores
the FASTA under both keys so structural-layer "NDUFA4" lookups and parquet-layer
"COXFA4" lookups both resolve. `ancestral_state_maps.json` is keyed as "COXFA4" and
`dar_gene` from the parquet is "COXFA4" — these are consistent. This is the correct
long-term pattern: each stage uses its own authoritative name, and explicit alias
bridges handle the transition.

**Filter:** In manuscript tables display as NDUFA4 (the widely-used name in the
literature). Query the parquet as `interpreted_gene == "COXFA4"`. Do not rename in
parquet or structural outputs — that would break the audit contract for those stages.  
**Affects:** Per-gene tables; CIV complex-level counts; figure labels.

---

## F10 — Incompatible contact pairs flag is a hypothesis, not a filter

**Source:** Compensating partners script (`likely_incompatible` column)  
**Issue:** Pairs flagged `likely_incompatible = True` (hbond or electrostatic contact
with sensitivity ≥ 0.50 and specificity ≥ 0.50) are candidates where the contact
change directly explains structural compensation — but this is a mechanistic hypothesis,
not a statistical result.  
**Filter:** Use `likely_incompatible` as a prioritisation flag for mutagenesis target
selection, not as an independent line of statistical evidence. Do not present it as a
significance criterion.  
**Affects:** Mutagenesis prioritisation rankings; supplementary tables of candidate
pairs.
