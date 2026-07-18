# Structural Mapping: Why Things Aren't Working and Potential Solutions

## Context

The structural mapping stage is the current pipeline bottleneck. Multiple rounds of
fixes have been applied (anchoring refactor, SDH anchor registry, RCSB+sequence chain
assignment, isoform proxy table, diagnostic failure layers) yet a substantial fraction
of variants remain unmapped. This note diagnoses *why* each approach is insufficient
and what targeted changes would actually move the needle.

---

## Failure Inventory (pre-fix state, Apr 20 2026)

| Failure class | Count | % of failures |
|---|---|---|
| `no_chain_in_pdb` (model gap) | 721 | 21% |
| `large_offset_candidate_*` (diagnostic only) | 1121 | 32% |
| `unmapped_reference_position` | 214 | 6% |
| `possible_large_offset_or_gap` | 154 | 4% |
| `isoform_proxy_gap` | 49 | 1% |

---

## Why Each Tried Approach Is Not Working

### 1. Chain assignment (RCSB API + sequence fallback) → still 721 `no_chain_in_pdb`

**Root cause is not the assignment logic — it is the model choice.**

The manifest routes CIII/CIV accessory subunits through the 5XTH respirasome model.
5XTH is a lower-resolution mixed-organism supercomplex. Many accessory subunits
(`CYC1`, `UQCRC1`, `UQCRC2`, `COX4I1`, `COX5A`, `COX6A1`, `COX6B1`) simply do not
have clean chains in 5XTH regardless of how chain assignment is done. No assignment
code can find a chain that is absent from the model.

### 2. Anchor exception registry (SDH-only extended rescue) → 1121 large-offset candidates remain diagnostic

**Root cause: the non-SDH audit correctly identified that non-SDH offsets are too
diffuse/multimodal for automatic rescue, so the code is working as designed.**

The SDH registry works because SDHA/B/C/D have a clean biological reason for a large
offset (mitochondrial targeting sequence, uniform ~25–40 aa shift). Non-SDH candidates
show scattered or contradictory offsets (e.g., NDUFA9 has both +34–51 and –11), meaning
the offset is not purely a mature-protein processing artifact — it may be local alignment
drift, model gap, or a mix of causes. Expanding the registry without understanding each
gene's offset biology would silently miscall positions.

### 3. Residue anchoring (3-layer diagnostic + transcript reconciliation) → 205 `unmapped_reference_position` + 154 `possible_large_offset_or_gap`

**Root cause: NDUFV3 and several CI subunits have positions that don't exist in the
chosen structure's chain — not because of offset, but because those segments are
genuinely unresolved loops in the PDB.**

The transcript reconciliation and diagnostic layers correctly detect this. But detection
is not solution: a position absent from the PDB's residue list cannot be rescued by
widening the search window.

### 4. Isoform proxy system → 49 `isoform_proxy_gap`

**Root cause: the proxy genes themselves (`COX4I1` → proxy for `COX4I2`,
`COX6A1` → proxy for `COX6A2`) are also affected by the chain/model gap problem
in the current manifest.**

The proxy logic is correct in design, but it resolves to the proxy's chain, and that
chain is absent in the current model for some subunits. So the proxy system propagates
the model-selection problem rather than solving it.

---

## Potential Solutions

### Solution A: Revise the structure manifest to use standalone complex models for CIII and CIV subunits

**Impact: high | Effort: low**

Update `data/reference/structure_model_manifest.tsv` to route:

- CIII subunits (`CYC1`, `UQCRC1`, `UQCRC2`, `UQCRFS1`, `UQCRB`, `UQCRQ`, `UQCRH`,
  `CYB5R3`) → `9HZL` as primary
- CIV subunits (`COX4I1`, `COX5A`, `COX6A1`, `COX6B1`, `COX7A2`, `COX8A`,
  `NDUFA4`) → `9I7U` as primary
- Retain 5XTH/5XTE as secondary/validation entries only

Why this helps: `9HZL` and `9I7U` are high-resolution standalone human CIII/CIV
structures with complete chain annotations. The structural phase plan already identifies
`9I7U` as the correct primary CIV source (mature NDUFA4-bound state), but the manifest
has not been updated to reflect this.

No code changes required — only manifest changes.

**Expected impact:** Eliminates most of the 721 `no_chain_in_pdb` rows for CIII/CIV
subunits. Also resolves the isoform proxy gap for `COX4I2`/`COX6A2` since their proxies
will now have chains.

---

### Solution B: Use UniProt transit peptide annotations to anchor mature-protein offsets for non-SDH CI/CV genes

**Impact: medium | Effort: medium**

Query the UniProt API for "Transit peptide" and "Chain" feature annotations for the
high-burden non-SDH genes. Use the annotated mature chain start position to compute
an exact fixed offset per gene. Add genes with a clearly annotated single mature-chain
start to `structural_anchor_exception_registry.tsv`.

Why this is different from the failed audit approach: the non-SDH audit searched for
the offset from the *observed* distribution in failure data, which was noisy. UniProt
annotations are derived from experimental evidence (Edman sequencing, mass spec) and
give an authoritative mature-protein start residue — a biologically grounded number,
not an inference from the structure.

**Candidate genes to check first:**
- `NDUFAB1`, `NDUFS3`, `NDUFS4`, `NDUFS6`, `NDUFS7` (CI matrix-arm)
- `ATP5F1A`, `ATP5F1B`, `ATP5F1D` (CV)
- `NDUFA9`, `NDUFA10` (CI peripheral arm)

**Files to change:**
- New utility: `src/utils/uniprot_mature_chain.py`
- `data/reference/structural_anchor_exception_registry.tsv`

**Expected impact:** Rescues a subset of the 1121 large-offset candidates for CI/CV
genes that have clean UniProt-annotated transit peptide cleavage sites.

---

### Solution C: Add AlphaFold structures for unresolved-loop genes and isoform-only genes

**Impact: medium | Effort: medium**

For genes where the PDB structure has genuine unresolved segments (`NDUFV3`) and for
isoform-specific genes without direct structural representation (`COX4I2`, `COX6A2`,
`COX7A1`), download AlphaFold AFDB models and add them to the structure panel as
supplemental sources.

AlphaFold models cover the full primary sequence including loop regions and
low-complexity segments that crystallographic models leave unresolved. They are not
experimental structures but extend coverage to regions genuinely absent from PDB.

AF2-derived contacts should be labeled as `af2_model` and treated as lower-confidence
structural evidence relative to experimental PDB contacts in all downstream summaries.

**Files to change:**
- New script: `src/data_download/00k_download_alphafold_structures.py`
- `data/reference/structure_model_manifest.tsv` — add AF2 entries with source tag
- `src/structural/00_map_davs_to_structure.py` — propagate structure-source column to output

**Expected impact:** Covers the `NDUFV3` `unmapped_reference_position` burden (205 rows)
and the 49 `isoform_proxy_gap` rows for tissue-specific isoforms not rescued by Solution A.

---

### Solution D (prerequisite): Fix the two rerun blockers before anything else

**Impact: critical | Effort: low**

Even if A–C are implemented, a clean rerun will fail due to two known blockers:

**D1: Overlapping mtDNA loci collapse**

`MT-ATP6/MT-ATP8` is still resolved to `MT-ATP6` everywhere via `.split("/")[0]`.

Fix in:
- `src/structural/00_map_davs_to_structure.py`
- `src/classify/00_classify_DAV.py`
- `src/utils/variant_record.py`

Approach: add a structured `loci` field and require explicit resolution to a primary
locus based on amino-acid change context rather than positional defaulting.

**D2: Hard-coded Docker path `/app`**

`src/structural/00_map_davs_to_structure.py` uses `Path("/app")` which only works
inside Docker at the correct mount path.

Fix: replace with `Path(__file__).resolve().parents[2]`

---

## Recommended Execution Order

1. **Fix D1 and D2** — these block a valid rerun regardless of other changes
2. **Implement Solution A** (manifest revision) — highest impact, no code changes
3. Re-run `00_map_davs_to_structure.py` and check how many `no_chain_in_pdb` rows are eliminated
4. **Implement Solution B** (UniProt offset annotations) for any CI/CV genes with clean transit peptide cleavage sites
5. **Implement Solution C** (AlphaFold supplemental sources) as a last resort for genuinely unresolved regions

---

## Verification Checkpoints

- After A: `no_chain_in_pdb` count should drop substantially for `CYC1`, `UQCRC1`,
  `UQCRC2`, `COX4I1`, `COX5A`, `COX6A1`, `COX6B1`
- After B: check `large_offset_candidate_*` rows for registry-added genes — should
  convert to `ok` or `isoform_offset_*`
- After C: `NDUFV3` `unmapped_reference_position` count should drop; isoform proxy
  gap rows should convert to `af2_mapped`
- End-to-end: rerun `01_find_compensating_partners.py` and verify that
  `compensatory_partners.csv` row count changes are interpretable

## Critical files

- `data/reference/structure_model_manifest.tsv`
- `data/reference/structural_anchor_exception_registry.tsv`
- `src/structural/00_map_davs_to_structure.py`
- `src/utils/variant_record.py`
- `src/classify/00_classify_DAV.py`
- `results/structural/structure_mapping_failure_audit.csv` (use to validate fix impact)

---

## Solution A Implementation Note (May 8 2026)

The manifest was already structurally corrected before this rerun:

- `CIII` primary = `9HZL`
- `CIV` primary = `9I7U`
- `5XTE`, `5Z62`, and `5XTH` retained as non-primary context models

The remaining problem was execution policy inside
`src/structural/00_map_davs_to_structure.py`.

The old mapper still attempted non-primary models for every structurally
eligible variant. That meant `5XTH` was still being treated like a routine
first-pass mapping target for `CIII`/`CIV` rows, even though the manifest had
already been corrected conceptually.

### What changed

The mapper now applies a primary-first structure policy:

- primary models are attempted for all structurally eligible rows
- validation/reference models are attempted only if at least one primary model
  mapped successfully for that variant
- secondary models skipped for lack of primary support are written explicitly as
  `secondary_not_attempted_no_primary_support`

This keeps validation/context models in the audit trail without allowing them
to inflate the primary structural failure burden.

### Post-fix structural state (May 8 2026 rerun)

Headline counts from the rerun:

- `mapped_direct`: `10621`
- `mapped_with_extended_offset_rescue`: `443`
- `mapped_with_isoform_offset`: `412`
- `secondary_not_attempted`: `1354`
- `mature_offset_candidate`: `526`
- `chain_assignment_or_model_gap`: `509`
- `residue_anchoring_failure`: `159`
- `isoform_proxy_gap`: `52`

Variant-summary support classes:

- `primary_and_validation`: `2576`
- `primary_only`: `3429`
- `unmapped`: `11162`

Per-genome mapped variant summaries:

- `mtDNA`: `317 / 388` mapped (`81.70%`)
- `nucDNA`: `5688 / 16779` mapped (`33.90%`)

### What improved

The key improvement is interpretability:

- non-primary models are no longer counted as routine first-pass failures when a
  variant has no primary structural support
- skipped secondary attempts are now visible and auditable rather than mixed
  into true model-gap or anchoring failures
- the residual `5XTH` burden now reflects validation-stage failure for rows that
  already had primary support, not a manifest-level routing mistake

Examples of explicit secondary skips produced by the rerun:

- `9TI4` validation `CI`: `342`
- `5XTH` reference `CI`: `342`
- `8H9T` validation `CV`: `118`
- `8H9U` validation `CV`: `118`
- `9I6F` validation `CIV`: `114`
- `5Z62` reference `CIV`: `114`
- `5XTH` validation `CIV`: `114`
- `5XTE` reference `CIII`: `46`
- `5XTH` validation `CIII`: `46`

### What did not improve

The fix did not eliminate all `CIII`/`CIV` chain-gap rows, because some rows do
map in the primary model and then still fail when a validation model is tested.
That remaining burden is now narrower and more honest:

- `5XTH` validation `CIII` `no_chain_in_pdb`: `356`
- `5XTH` validation `CIV` `no_chain_in_pdb`: `205`

### Isoform proxy nuance after the primary-first fix

The slight change in `isoform_proxy_gap` from `49` to `52` is not because the
primary-first policy blocked access to `9I6F`.

Audit of the current mapping rows shows:

- `9I7U` contains the relevant proxy chains:
  - `COX4I1` chain `D`
  - `COX6A1` chain `G`
  - `COX7A2` chain `J`
- `9I6F` also contains the same proxy chains

So the proxy system still reaches the primary `CIV` model correctly.

The current `isoform_proxy_gap` rows are mostly validation-context failures in
`5XTH`, not primary-model access failures:

- `COX4I2`: `34` `5XTH` validation `no_chain_in_pdb`
- `COX6A2`: `16` `5XTH` validation `no_chain_in_pdb`
- `COX7A1`: `2` `5XTH` validation `no_chain_in_pdb`

This means:

- the proxy comments that previously named only `9I6F` were outdated and were
  corrected in code
- the remaining proxy-gap burden is now a respirasome validation-model coverage
  limitation, not a broken proxy-to-primary routing problem

So the current conclusion is:

- the manifest/source map is now aligned with the intended biology
- the mapper now respects that source hierarchy
- the remaining structural bottlenecks are no longer “wrong primary model”
  problems
- the next unresolved issues are dominated by:
  - residual validation-model chain gaps
  - non-`SDH` mature-offset candidates
  - unresolved residue anchoring, especially `NDUFV3`

## Active mapping panel restriction (May 9 2026)

The active structural denominator was then restricted further to use only the
high-resolution single-complex panel for structural mapping.

Removed from the active mapping panel:

- `5XTH` respirasome-context models
- `5XTE` legacy lower-resolution `CIII` reference
- `5Z62` legacy lower-resolution `CIV` reference

Retained in the active mapping panel:

- `CI`: `9I4I`, `9TI4`
- `CII`: `8GS8`
- `CIII`: `9HZL`
- `CIV`: `9I7U`, `9I6F`
- `CV`: `8H9S`, `8H9T`, `8H9U`

### Reduced-panel rerun outcome

Headline results from the reduced active panel:

- `dar_structure_map.csv`: `37532` rows
- `mapped model rows`: `9058`
- `dar_contacts_cbcb8A.csv`: `92894` rows
- `dar_mito_nuc_contacts.csv`: `2101` rows
- mito-nuclear interface variants: `378`
- AA-level cDAVs at the mito-nuclear interface: `165`

This panel restriction removes the respirasome-driven `5XTH` proxy/model-gap
artifacts from the active denominator and makes the structural analysis easier
to interpret as a single-complex mapping framework.
