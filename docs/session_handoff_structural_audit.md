# Session Handoff: Structural Mapping Stage Audit

**Date:** 2026-05-12  
**Session:** Structural stage audit + CLAUDE.md improvement  
**Continuation from:** `/Users/ad2347/.claude/projects/-Users-ad2347-Documents-OxPhos-DAV/961ca287-93ab-4e8f-b735-bf2d35a1bb9e.jsonl`

---

## What Was Done

### 1. CLAUDE.md Rewrite

The existing `CLAUDE.md` was 5 lines (a project description paragraph + a Docker note). It was completely rewritten to include:

- Project summary (c-DAR / cDAV definitions, AA-level vs NT-level)
- Docker environment setup and devcontainer usage
- How to run pipeline scripts
- Test commands (all tests, single file, single test)
- All 7 pipeline stages with inputs, outputs, and key script names
- Key data flow diagram (`data/raw/` → `data/derived/` → `results/`)
- Key utility modules (`variant_record.py`, `alignment_parser.py`, `exception_registry.py`, `gene_reference.py`, `mt_overlap.py`)
- Reference files table
- Notes reading order (7 files, in order)

### 2. Structural Stage Audit

Comprehensive audit of `src/structural/00_map_davs_to_structure.py` (1509 lines) and `src/structural/01_find_compensating_partners.py` (862 lines), cross-referenced against all associated notes.

---

## Structural Stage: High-Level Flowcharts

### `00_map_davs_to_structure.py`

```
INPUT
  variants_master_classified.parquet
  structure_model_manifest.tsv
  structural_anchor_exception_registry.tsv
  PDB/CIF files (local cache)
  RCSB REST API
        │
        ▼
LOAD ACTIVE MODELS (active=True, sorted by complex + priority)
  Current active panel:
    CI:   9I4I (primary p1), 9TI4 (validation p2)
    CII:  8GS8 (primary p1)
    CIII: 9HZL (primary p1)
    CIV:  9I7U (primary p1), 9I6F (validation p2)
    CV:   8H9S (primary p1), 8H9T (p2), 8H9U (p3)
  Inactive: 5XTH, 5XTE, 5Z62
        │
        ▼
FOR EACH VARIANT (grouped by gene)
  ┌─ ELIGIBILITY GATE
  │   • Skip if gene not in structure panel
  │   • Apply isoform proxy map:
  │     COX4I2→COX4I1, COX6A2→COX6A1,
  │     COX7A1→COX7A2, ATP5MC2/3→ATP5MC1
  │   • Apply TOGA→canonical: COXFA4→NDUFA4
  │
  ├─ nucDNA COORD REMAP (remap_aa_coord_to_structure_space)
  │   NM_ position (ClinVar) → ENST position (TOGA)
  │   via transcript_position_maps.json
  │   ⚠ THREE-SPACE COORD GAP — main failure source for nucDNA
  │
  └─ FOR EACH MODEL (primary-first policy)
      ┌─ PRIMARY-FIRST GATE (post-May 8 2026)
      │   Validation models only attempted if ≥1 primary mapped
      │
      ├─ CHAIN ASSIGNMENT (fetch_chain_gene_map / assign_chain_to_gene)
      │   1. RCSB REST API → gene name per chain
      │   2. Fallback: PairwiseAligner local alignment, 30% threshold
      │   ⚠ 30% threshold uncalibrated for short homologous subunits
      │
      ├─ POSITION MAP (map_refseq_to_chain)
      │   Global alignment: canonical protein seq → PDB chain seq
      │   Builds {ref_pos → chain_pos} dict
      │
      ├─ RESIDUE ANCHORING (find_anchor) — 3 stages
      │   Stage 1: Direct match (window=0)
      │            → status: mapped_direct
      │   Stage 2: Sliding window ±10
      │            → status: mapped_with_isoform_offset (if offset)
      │   Stage 3: Diagnostic window ±80
      │            → status: mature_offset_candidate (label only, no rescue)
      │   Stage 4: Registry-gated extended rescue
      │            Only if allow_extended_anchor=True in exception registry
      │            Currently only SDHA/B/C/D (max_offset=80)
      │            → status: extended_offset_rescue
      │
      │   FAILURE LABELS (classify_anchor_failure):
      │     missing_from_pos_map / not_present_in_structure /
      │     outside_chain_range / mature_offset_candidate /
      │     anchoring_exhausted
      │   ⚠ BUG L983-984: always passes extended_window=80
      │     regardless of registry max_offset
      │
      ├─ CONTACT EXTRACTION
      │   Cβ–Cβ ≤ 8 Å (Cα for Gly)
      │   classify_contact: hbond > electrostatic > hydrophobic > vdw
      │
      └─ STATUS TAG MUTATION ⚠ BUG L1049
          status += f"(proxy={proxy_gene})" after classification
          Breaks status_category() startswith matching
          Proxy-mapped rows land in "other" instead of correct category

OUTPUT
  results/structural/structure_mapping_results.parquet
  results/structural/structure_contacts.parquet
  results/structural/mapping_summary.tsv
```

---

### `01_find_compensating_partners.py`

```
INPUT
  variants_master_classified.parquet  (AA-level cDAVs only)
  structure_contacts.parquet
  TOGA AA alignments per gene
  IQTree ancestral state maps (results/phylo/)
  VertLife mammal tree
        │
        ▼
LOAD cDAVs
  Filter: is_cdav_aa == True
  Build: {variant_id → set(species_with_disease_allele)}
  ⚠ L427: cdav_spp = set(var.get("lineages_with_disease_allele", []))
    Silent empty set if field is JSON string in parquet (type safety)
        │
        ▼
LOAD CONTACTS
  Join contacts to cDAVs on variant_id
  Group by (DAR × contact_residue × alt_AA)
        │
        ▼
PASS 1: PER (DAR × CONTACT × ALT_AA)
  ┌─ FISHER'S EXACT TEST
  │   2×2: species with/without disease allele × contact AA
  │   ⚠ INVALID as primary test (no phylogenetic correction)
  │   Kept for comparison; BH FDR computed per-DAR
  │   ⚠ FDR SCOPE ISSUE: per-DAR BH but globally applied threshold
  │
  └─ BRANCH CO-OCCURRENCE
      IQTree ancestral state maps → binary trait per branch
      Fisher on branch counts (with phylogenetic correction)
      BH FDR per-DAR
        │
        ▼
PASS 2: BATCH PAGEL'S DISCRETE (R subprocess)
  For each (DAR × contact) pair:
    Rscript pagel_discrete.R → p-value for correlated evolution
  ⚠ L543: rec_idx = len(all_records) + len(raw_tests)
    Fragile index but correct given append-only all_records
  BH FDR applied globally across all Pagel records
        │
        ▼
SIGNIFICANCE GATE (_get_sig)
  PRIMARY: (pagel_fdr ≤ 0.10 OR branch_cooccur_fdr ≤ 0.10)
           AND low_power == False
  FALLBACK: fisher_fdr ≤ 0.10
  (Fallback activates when phylo tests unavailable)
        │
        ▼
OUTPUT
  results/structural/compensatory_partners.csv
  results/structural/compensatory_partners_all.csv  (unfiltered)
```

---

## Identified Problems (Priority-Ranked)

### P1 — Proxy status tag mutation (correctness bug)

**File:** [src/structural/00_map_davs_to_structure.py](src/structural/00_map_davs_to_structure.py#L1049)  
**Line:** 1049

```python
# BUG: status is mutated after classify_anchor_failure() / status_category() calls
status += f"(proxy={proxy_gene})"
```

**Impact:** Proxy-mapped rows (COX4I2, COX6A2, COX7A1, ATP5MC2/3) land in the `"other"` bucket in all `status_category()` filtering instead of `"mapped_with_isoform_offset"` or `"mapped_direct"`. These rows are excluded from downstream counts and from contact extraction if the filter is applied before the tag mutation. Affects 412 rows in the current run (isoform_offset count).

**Fix:** Append proxy tag to a separate `status_detail` field, not to the primary `status` field.

---

### P2 — `classify_anchor_failure` always passes `extended_window=80`

**File:** [src/structural/00_map_davs_to_structure.py](src/structural/00_map_davs_to_structure.py#L983)  
**Lines:** 983-984

```python
# BUG: always passes 80 regardless of registry max_offset
classify_anchor_failure(..., diagnostic_anchor, 10, 80)
```

The function signature is `classify_anchor_failure(gene, pos, registry_entry, diagnostic_anchor, default_window, extended_window)`. The call hardcodes `extended_window=80` instead of reading `registry_entry.max_offset`. For non-SDH genes with no registry entry, this produces incorrect offset range labels in the failure diagnostics.

**Impact:** `mature_offset_candidate` labels may overstate the range of plausible mature-protein offsets for non-SDH genes. Minor diagnostic accuracy issue; does not affect which rows are rescued.

**Fix:**
```python
max_offset = registry_entry.max_offset if registry_entry else 10
classify_anchor_failure(..., diagnostic_anchor, 10, max_offset)
```

---

### P3 — `lineages_with_disease_allele` type safety

**File:** [src/structural/01_find_compensating_partners.py](src/structural/01_find_compensating_partners.py#L427)  
**Line:** 427

```python
cdav_spp = set(var.get("lineages_with_disease_allele", []))
```

If the parquet field is stored as a JSON string (`'["Homo_sapiens", ...]'`) rather than a list, `set(json_string)` produces a set of individual characters — silently wrong. The field type depends on how `variants_master_classified.parquet` was written.

**Fix:**
```python
raw = var.get("lineages_with_disease_allele", [])
cdav_spp = set(json.loads(raw) if isinstance(raw, str) else raw)
```

---

### P4 — Contact deduplication missing across models

**File:** [src/structural/01_find_compensating_partners.py](src/structural/01_find_compensating_partners.py)

The same (DAR × contact_residue) pair can appear from multiple models (e.g., 9I7U and 9I6F both contributing CIV contacts). `n_contacts_tested` in the partner analysis inflates when models overlap in coverage. No deduplication step exists before the Fisher/Pagel test loop.

**Impact:** Inflated contact counts; duplicate tests for the same biological contact. Minor effect on FDR given the per-DAR BH structure, but inflates the apparent contact burden.

**Fix:** Deduplicate on `(variant_id, contact_chain_resnum, contact_chain_id)` before the test loop, keeping the highest-confidence model's contact record.

---

### P5 — Three-space coordinate gap (structural root cause)

**File:** [src/structural/00_map_davs_to_structure.py](src/structural/00_map_davs_to_structure.py) — `remap_aa_coord_to_structure_space`

nucDNA variants traverse: **ClinVar NM_ position → TOGA ENST position → PDB chain residue position**. Each translation can fail independently:
- NM_→ENST: requires `transcript_position_maps.json` to cover the specific NM_ transcript and position
- ENST→PDB: requires the canonical protein sequence to align well to the PDB chain (global alignment)

**Current impact:** 5688/16779 nucDNA variants mapped (33.9%). The remaining 66.1% fail predominantly at the NM_→ENST step (missing transcript coverage) or the ENST→PDB step (isoform/transit peptide mismatch).

**Partial solution (not yet implemented):** UniProt transit peptide annotations (Solution B from `structural_mapping_diagnosis_and_solutions.md`) would improve the ENST→PDB alignment step for nuclear-encoded mitochondrial subunits.

---

### P6 — Conservative ±10 anchor window: non-SDH large-offset genes

**File:** [src/structural/00_map_davs_to_structure.py](src/structural/00_map_davs_to_structure.py) — `find_anchor`  
**Registry:** [data/reference/structural_anchor_exception_registry.tsv](data/reference/structural_anchor_exception_registry.tsv)

Only SDHA/B/C/D are granted extended rescue (max_offset=80). Other nuclear-encoded subunits with transit peptides (NDUFV1, NDUFS2, UQCRC1/2, COX5A, ATP5A1, etc.) have systematic mature-protein offsets of 20–80 aa. These are labelled `mature_offset_candidate` but not rescued.

Current count of `mature_offset_candidate` rows: ~526 (post-May 8 2026 rerun).

**Path forward:** Implement Solution B (UniProt mature annotations) before expanding the registry, to avoid rescuing wrong offsets. See `structural_mapping_diagnosis_and_solutions.md` for the full plan.

---

### P7 — Fisher FDR scope mismatch

**File:** [src/structural/01_find_compensating_partners.py](src/structural/01_find_compensating_partners.py#L591)  
**Lines:** 591-593

Fisher BH FDR is computed per-DAR (each DAR's test p-values corrected independently), but `_get_sig` applies a fixed 0.10 threshold. This is weaker control than a global BH pass over all Fisher tests simultaneously.

**Impact:** Minor inflation of Fisher false positives in the fallback path. The primary significance path (Pagel/branch) uses a global BH pass, so this only affects records where phylo tests are unavailable.

---

### P8 — Chain assignment 30% threshold uncalibrated

**File:** [src/structural/00_map_davs_to_structure.py](src/structural/00_map_davs_to_structure.py) — `assign_chain_to_gene`

The fallback local alignment threshold is 30% identity. OXPHOS has several short subunits (e.g., ND3, ND4L, ATP8, CYC1) where the best-match chain may be homologous at only 25–35% identity to a wrong gene. No calibration data exists to validate this threshold.

**Impact:** Potential chain misassignment for short or highly divergent subunits in the fallback path (RCSB API unavailable or returning no results).

---

### P9 — No structure-panel eligibility registry

There is no explicit list of which genes are expected to be absent from the panel (e.g., accessory subunits not present in any active PDB model). Genes absent from all models are indistinguishable from genes that failed chain assignment.

**Impact:** Audit statistics mix "expected absent" and "unexpected failure" categories. The `chain_assignment_or_model_gap` failure class (629 rows) is a mix of both.

---

## Current Mapping Summary (Post-May 8/9 2026 Rerun)

| Status | Count |
|---|---|
| mapped_direct | 10621 |
| extended_offset_rescue | 443 |
| isoform_offset | 412 |
| secondary_not_attempted | 1354 |
| mature_offset_candidate | 526 |
| chain_assignment_or_model_gap | 509 |
| residue_anchoring_failure | 159 |
| isoform_proxy_gap | 52 |
| **mtDNA mapped** | 317/388 (81.7%) |
| **nucDNA mapped** | 5688/16779 (33.9%) |

From `structural_phase_plan.md` (full panel stats):
- 33779 classified rows total
- 11049 mapped model rows
- 5562 unique variants mapped
- 113442 contacts
- 2981 cross-genome contacts
- 400 mito-nuclear interface variants
- 178 AA-level cDAVs at interface

---

## Recommended Next Steps (by priority)

1. **Fix proxy status tag bug** (P1) — low effort, correctness fix, affects 412 rows in output counts
2. **Fix `classify_anchor_failure` extended_window param** (P2) — one-line fix, improves diagnostic accuracy
3. **Add `lineages_with_disease_allele` type guard** (P3) — defensive fix, prevents silent empty-set bug
4. **Implement UniProt transit peptide annotations** (P5/P6) — medium effort, unblocks ~526 `mature_offset_candidate` rows
5. **Add contact deduplication** (P4) — medium effort, improves statistical validity of partner analysis
6. **Fix Fisher FDR scope** (P7) — low effort, tightens statistical rigor in fallback path
7. **Add structure-panel eligibility registry** (P9) — improves audit interpretability

---

## Key Files Referenced

| File | Role |
|---|---|
| [src/structural/00_map_davs_to_structure.py](src/structural/00_map_davs_to_structure.py) | Main structural mapping (1509 lines) |
| [src/structural/01_find_compensating_partners.py](src/structural/01_find_compensating_partners.py) | Co-evolution partner analysis (862 lines) |
| [notes/structural_phase_plan.md](notes/structural_phase_plan.md) | Authoritative structural-stage plan |
| [notes/structural_mapping_major_problems_report.md](notes/structural_mapping_major_problems_report.md) | Current-state failure report |
| [notes/structural_mapping_diagnosis_and_solutions.md](notes/structural_mapping_diagnosis_and_solutions.md) | Diagnosis + 4 solutions (A implemented, B/C pending) |
| [notes/structural_anchoring_refactor_notes.md](notes/structural_anchoring_refactor_notes.md) | Anchoring method provenance |
| [data/reference/structure_model_manifest.tsv](data/reference/structure_model_manifest.tsv) | Active PDB panel + priorities |
| [data/reference/structural_anchor_exception_registry.tsv](data/reference/structural_anchor_exception_registry.tsv) | Extended-anchor registry (SDH only) |
