# Structural Mapping Eligibility Gaps

Identified during structural mapping audit, 2026-05-13.

## Problem 1: COXFA4/NDUFA4 gene symbol mismatch — 17 variants silently dropped

**Where:** `src/structural/00_map_davs_to_structure.py`, lines 745–746

**Why:** `_TOGA_TO_CANONICAL` remaps COXFA4→NDUFA4 for TOGA FASTA filenames (line 725)
but is never applied to `interpreted_gene` from the classified parquet (line 745).
ClinVar stores some variants under the old symbol COXFA4. When the script looks up
`GENE_COMPLEX.get("COXFA4")`, it returns None (GENE_COMPLEX has "NDUFA4"). The variant
also fails `gene in ref_seqs` because ref_seqs is keyed by "NDUFA4". Both conditions
make `structure_mapping_eligible = False`, so the variant is silently labelled
`structural_ineligible`.

**Impact:** 17 classified NDUFA4 variants are not mapped to CIV structures (9I7U / 9I6F)
and do not appear in `dar_contacts_cbcb8A.csv`, so they are also missing from the
compensating partners analysis in `01_find_compensating_partners.py`.

**Fix (one line):** After line 745 in `00_map_davs_to_structure.py`, add:
```python
gene = _TOGA_TO_CANONICAL.get(gene, gene)
```
This normalises `interpreted_gene` to the canonical HGNC name before all downstream
lookups (GENE_COMPLEX, ref_seqs, structure_mapping_eligible). The normalised name also
propagates into the output `interpreted_gene` column, making CSVs and the parquet
consistent with the canonical symbol used everywhere else.

**Downstream propagation:** After the fix, rerun `00_map_davs_to_structure.py` then
`01_find_compensating_partners.py`. No other scripts need changing — the alignment
loader in `01_find_compensating_partners.py` already has its own `_TOGA_TO_CANONICAL`
remapping (line 93) that converts `COXFA4_aa_alignment.fasta → "NDUFA4"`.

---

## Problem 2: Structural mapping rate is misreported

**Where:** `src/structural/00_map_davs_to_structure.py`, summary section, lines 1501–1512

**Why:** The summary prints "Classified rows loaded: N" vs "Mapped model rows: M", implying
an overall mapping rate of M/N. But ~16,844 of the ~17,006 rows have
`classification_status = "skipped_by_policy"` (benign/VUS variants not eligible for
structural analysis). These are always `structural_ineligible` and were never attempted.
Reporting M/N as the mapping rate conflates deliberate policy exclusions with genuine
mapping failures.

**Accurate denominator:** Only variants with `structure_mapping_eligible = True` were
attempted. Pre-fix counts (post CV-priority swap, pre COXFA4 fix):

| Category | Count |
|---|---|
| Total classified rows loaded | 17,006 |
| skipped_by_policy (benign/VUS, not attempted) | 16,844 |
| Mapping-eligible | 162 |
| Mapped (unique variants) | ~130 |
| Eligible but failed | ~32 |
| True mapping rate | ~80% |

**Fix:** Replace lines 1501–1512 in `00_map_davs_to_structure.py` to compute and print
separate counts for: skipped_by_policy, mapping-eligible, mapped unique variants,
eligible-but-failed, and true mapping rate (%).

---

## Problem 3: Genes with no structure panel coverage (documented gap)

The following 12 genes are in the canonical OXPHOS variant set with classified
variants, but are absent from GENE_COMPLEX and from all active structure chains:

| Gene | Complex | Classified variants | Reason absent |
|---|---|---|---|
| ATP5ME | CV | 18 | Peripheral CV subunit, not resolved in 8H9S/8H9T/8H9U |
| ATP5MG | CV | 18 | Peripheral CV subunit |
| COXFA4L2 | CIV | 17 | Tissue-specific CIV regulatory subunit |
| COX8C | CIV | 15 | Tissue-specific CIV isoform |
| ATP5IF1 | CV | 14 | CV inhibitory factor, not a structural subunit |
| COX6B2 | CIV | 12 | Testis-specific CIV isoform |
| COX7B2 | CIV | 12 | Testis-specific CIV isoform |
| ATP5MF | CV | 11 | Peripheral CV subunit |
| ATP5MK | CV | 9 | Peripheral CV subunit |
| ATP5MJ | CV | 3 | Peripheral CV subunit |
| COXFA4L3 | CIV | 2 | Tissue-specific CIV regulatory subunit |

**Status:** Accepted limitation. These are not added to GENE_COMPLEX because their
chains are absent from all active cryo-EM structures. Disclose in the manuscript as
outside the structural panel.

**Additional 63 variants** from genes that ARE in the panel fail due to coordinate
resolution failures (`unresolved` classification status or `position_not_in_enst`).
These are per-variant issues, not panel gaps.

---

## Problem 4: Downstream propagation to compensating partners

`01_find_compensating_partners.py` reads `dar_contacts_cbcb8A.csv` and uses
`dar_locus` and `contact_gene` for cross-species alignment lookups. Because NDUFA4
contacts are absent from that CSV (Problem 1), the compensating partner analysis
currently has no NDUFA4 entries. After fixing Problem 1 and rerunning both scripts,
NDUFA4 contacts will be present and will flow correctly.

The 12 no-panel-coverage genes (Problem 3) will remain absent from the contacts CSV
after the fix — this is correct behaviour and they should not appear in compensating
partner results.
