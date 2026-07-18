## Structural Phase Plan

This note defines the next structural-analysis phase for the OXPHOS DAV study.

The structural stage is a derived mechanistic layer built on top of the
classified master table. It should not be the first stage where variants
disappear implicitly.

## Structural source validation

The local structural source map was validated against the current official RCSB
entry pages for the key human OXPHOS structures used in this project:

- `9I4I`: https://www.rcsb.org/structure/9I4I
- `9TI4`: https://www.rcsb.org/structure/9TI4
- `8GS8`: https://www.rcsb.org/structure/8GS8
- `9HZL`: https://www.rcsb.org/structure/9HZL
- `9I7U`: https://www.rcsb.org/structure/9I7U
- `9I6F`: https://www.rcsb.org/structure/9I6F
- `5Z62`: https://www.rcsb.org/structure/5Z62
- `8H9S`: https://www.rcsb.org/structure/8H9S
- `8H9T`: https://www.rcsb.org/structure/8H9T
- `8H9U`: https://www.rcsb.org/structure/8H9U
- `5XTH`: https://www.rcsb.org/structure/5XTH
- `5XTE`: https://www.rcsb.org/structure/5XTE

### Similarities between the audit map and the current local panel

The current local structural panel was already broadly aligned with the correct
human OXPHOS structure set:

- `CI`: `9I4I` and `9TI4`
- `CII`: `8GS8`
- `CIII`: `9HZL`, with `5XTE` and `5XTH` as supporting context
- `CV`: `8H9S` and `8H9T`

These were directionally correct choices for a modern human-focused structural
panel.

### Discrepancies identified during source validation

The audit exposed several important discrepancies in the previous local source
map.

#### 1. Complex IV primary/validation roles were reversed

The earlier local manifest treated:

- `9I6F` as primary
- `9I7U` as validation

This was biologically backwards.

The more defensible interpretation is:

- `9I7U` = mature `NDUFA4`-bound `CIV` state
- `9I6F` = `HIGD2A`-bound late-assembly intermediate

Therefore:

- `9I7U` should be the primary `CIV` source for mature-state mapping
- `9I6F` should be retained as a validation / assembly-intermediate source

#### 2. `5Z62` was missing from the local source map

`5Z62` is an important legacy human `CIV` reference and should remain part of
the structural provenance map, even if it is not the main primary mapping
source.

Its role is best defined as:

- historical / reference `CIV` source
- not the default primary model for the present pipeline

#### 3. `8H9U` was missing from the ATP synthase panel

The current local panel included:

- `8H9S`
- `8H9T`

but omitted:

- `8H9U`

Given that `8H9S/T/U` represent a rotational state series for human ATP
synthase, `8H9U` should be retained as an additional validation state.

#### 4. `9TI4` should not be treated as the primary respirasome source

The provided audit map grouped `9TI4` as:

- `Respirasome primary`

That is not the best representation for the current pipeline contract.

For this project, `9TI4` is better treated as:

- a `CI` validation model in mitochondrial context

while respirasome-context support should remain explicit and secondary, using:

- `5XTH` as a legacy respirasome context source
- `9I7U` and `9I6F` for `CIV` state-aware interpretation

This keeps the source map biologically clearer and avoids conflating:

- standalone / complex-focused models
- respirasome-context models
- assembly intermediates

#### 5. Resolution metadata in the local manifest was stale

The previous local manifest contained several resolution values that were too
coarse or outdated relative to the current RCSB entries.

The manifest was corrected to the current values used in the source map:

- `9I4I`: `2.63 Å`
- `9TI4`: `2.66 Å`
- `8GS8`: `2.86 Å`
- `9HZL`: `2.52 Å`
- `9I7U`: `3.15 Å`
- `9I6F`: `2.95 Å`
- `5Z62`: `3.60 Å`
- `8H9S`: `2.53 Å`
- `8H9T`: `2.77 Å`
- `8H9U`: `2.61 Å`
- `5XTH`: `3.90 Å`
- `5XTE`: `3.60 Å`

### Scientifically robust source map for this pipeline

The corrected structural source map for the current analysis should be:

| Complex | Role in pipeline | PDB | Interpretation |
| --- | --- | --- | --- |
| CI | Primary | `9I4I` | best current primary human active-state `CI` model |
| CI | Validation | `9TI4` | high-resolution alternative human `CI` model in mitochondrial context |
| CI / Respirasome context | Reference | `5XTH` | legacy lower-resolution respirasome-context source only |
| CII | Primary | `8GS8` | best current human `CII` model with bound ubiquinone |
| CIII | Primary | `9HZL` | best current human `CIII2` primary source |
| CIII | Reference | `5XTE` | legacy standalone human `CIII` reference |
| CIII / Respirasome context | Validation | `5XTH` | context validation only |
| CIV | Primary | `9I7U` | mature `NDUFA4`-bound `CIV` state |
| CIV | Validation | `9I6F` | `HIGD2A`-bound late-assembly intermediate |
| CIV | Reference | `5Z62` | legacy intact human `CIV` reference |
| CIV / Respirasome context | Validation | `5XTH` | context validation only |
| CV | Primary | `8H9S` | ATP synthase state 1 |
| CV | Validation | `8H9T` | ATP synthase state 2 |
| CV | Validation | `8H9U` | ATP synthase state 3 |

### Source-map rules for downstream use

The structural pipeline should follow these rules:

- primary models are used for default residue-to-structure mapping and default
  contact interpretation
- validation models are used to test contact and mapping robustness across
  state changes
- reference models are retained for provenance and sensitivity analysis, but
  should not override better modern primaries
- respirasome-context models should not be treated as interchangeable with
  standalone complex models
- assembly-intermediate models should not be treated as mature-state defaults
  unless the biological question explicitly concerns assembly state

## Active mapping panel policy

Current active policy for the structural mapping stage:

- use only the high-resolution single-complex structure panel for active mapping
- do not rely on legacy respirasome-context structures in `5XTH`
- do not use lower-resolution legacy reference structures (`5XTE`, `5Z62`) for
  active residue mapping

This means the active mapping panel is now:

- `CI`: `9I4I` primary, `9TI4` validation
- `CII`: `8GS8` primary
- `CIII`: `9HZL` primary
- `CIV`: `9I7U` primary, `9I6F` validation
- `CV`: `8H9S` primary, `8H9T` validation, `8H9U` validation

The older structures remain documented for provenance, but they are no longer
part of the active structural denominator.

## Study-aligned goals

The structural analysis should support the following biological questions:

- do cDAVs and uDAVs differ in structural location, interface exposure, or
  contact neighborhood
- are putative compensatory partner residues concentrated among direct
  structural neighbors of cDAVs
- are these candidate partners intra-protein, intra-genomic, or inter-genomic
- do the strongest candidate partner relationships remain stable across
  alternative structural models or conformational states

## Conceptual strengths of the current implementation

- residue-level structural mapping is explicit rather than inferred indirectly
- contact classes are annotated, not just Euclidean distance thresholds
- chain-to-gene assignment includes sequence-based fallback logic
- residue anchoring already includes a local offset-correction window
- isoform proxy usage is explicit instead of silent
- cross-genome contacts are separated from within-genome contacts
- compensatory partner testing already recognizes that phylogenetically naive
  Fisher tests are not valid as a primary inference layer

## Main limitations to fix

### 1. Old input contract

The current structural scripts still read legacy compatibility JSON outputs
under `data/annotations/curated/`. They should consume the canonical classified
Parquet layer instead:

- `data/derived/classified/variants_master_classified.parquet`
- optional downstream sensitivity subsets:
  - `classified_all.parquet`
  - `classified_clean.parquet`
  - `classified_warning.parquet`

### 2. One preferred structure per complex is too narrow

Using one static structure per complex is the main conceptual weakness in the
current structural layer.

Why this matters:

- some residues are unresolved in one model but resolved in another
- contact calls can change between states or reconstruction contexts
- interface classification can shift across structures
- isoform proxies are harder to evaluate if only one model is considered

The refactor should use a small curated structure panel per complex rather than
a single hardcoded PDB.

### 3. Protomers and repeated chain copies

The current cache stores one `gene -> chain` mapping and can overwrite repeated
copies. That is not robust for oligomeric assemblies. The new logic should
support:

- `gene -> [chain1, chain2, ...]`
- per-chain positional maps
- chain selection based on the best anchored residue match

### 4. Structural output should preserve failures

The mapping stage should emit explicit failure states instead of only retaining
successfully mapped rows. This is required for denominator-aware comparisons.

### 5. Direct contacts are one mechanism layer, not the whole mechanism

Static structure contacts are useful, but they should be interpreted as one
candidate mechanistic layer. Compensation can also involve:

- second-shell packing effects
- interface remodeling
- long-range conformational coupling
- cofactor- or membrane-mediated effects

## Curated alternative structure panel

The current structure panel should be managed through a manifest rather than
hardcoded PDB IDs in the script.

Proposed primary and validation models:

| Complex | Primary | Validation / Alternative | Notes |
| --- | --- | --- | --- |
| CI | `9I4I` | `9TI4`, `5XTH` | high-resolution human CI plus supercomplex validation |
| CII | `8GS8` | none-human-default, optional non-human comparative validation only | no strong alternative human standalone model currently selected |
| CIII | `9HZL` | `5XTE`, `5XTH` | standalone human CIII plus respirasome context |
| CIV | `9I6F` | `9I7U`, `5XTH` | HIGD2A-bound and NDUFA4-bound states plus supercomplex context |
| CV | `8H9S` | `8H9T` | human ATP synthase rotational states 1 and 2 |

Official source pages used for this panel:

- `9I4I`: https://www.rcsb.org/structure/9I4I
- `9TI4`: https://www.rcsb.org/structure/9TI4
- `8GS8`: https://www.rcsb.org/structure/8GS8
- `9HZL`: https://www.rcsb.org/structure/9HZL
- `5XTE`: https://www.rcsb.org/structure/5XTE
- `5XTH`: https://www.rcsb.org/structure/5XTH
- `9I6F`: https://www.rcsb.org/structure/9I6F
- `9I7U`: https://www.rcsb.org/structure/9I7U
- `8H9S`: https://www.rcsb.org/structure/8H9S
- `8H9T`: https://www.rcsb.org/structure/8H9T

### Primary versus validation structure usage

The structural analysis should use the structure panel asymmetrically rather
than treating all models as equivalent.

Primary structures are used for:

- the default residue-to-structure mapping when a variant maps successfully
- the default contact extraction used in downstream summaries
- the headline structural figures unless model disagreement is central to the
  point being illustrated

Validation structures are used for:

- checking whether a residue is resolved in an alternative human model or state
- checking whether a contact found in the primary model is preserved across
  another model or a supercomplex context
- identifying state-specific or model-specific contacts
- refining confidence labels for candidate compensatory partners

Interpretation rule:

- a contact present only in a validation model is hypothesis-generating
- a contact present only in the primary model is usable but should be labeled as
  single-model support
- a contact present across primary and one or more validation models has the
  strongest structural support

## Refactor targets

### Step 1. Structure manifest

Add a canonical manifest file with at least:

- `complex_id`
- `pdb_id`
- `role`
- `state_label`
- `organism`
- `method`
- `resolution_angstrom`
- `source_url`
- `active`
- `priority`

### Step 2. Structure mapping

Refactor `src/structural/00_map_davs_to_structure.py` to:

- read the canonical classified Parquet input
- use `interpreted_gene`, not `locus.split("/")`
- load the structure panel from the manifest
- support multiple candidate chains per gene
- write one mapping row per `(variant, structure model)`
- write contact rows tagged by `pdb_id` and `model_role`

Recommended mapping outputs:

- `results/structural/dar_structure_map.csv`
  - retained for compatibility, but now one row per `(variant_id, pdb_id)`
- `results/structural/dar_contacts_cbcb8A.csv`
  - all contact rows across all structure models
- `results/structural/dar_mito_nuc_contacts.csv`
  - filtered cross-genome contacts
- `results/structural/structure_model_summary.csv`
  - per-variant summary across the structure panel

### Step 3. Partner analysis

Refactor `src/structural/01_find_compensating_partners.py` to:

- read the classified Parquet master instead of compatibility JSONs
- identify cDAV rows from the canonical classified fields
- continue writing all tested pairs as the primary output
- treat significant partner views as derived outputs only

### Step 4. Confidence and consensus

Per variant or per tested pair, compute:

- `n_models_attempted`
- `n_models_mapped`
- `n_models_with_contact`
- `contact_consensus_fraction`
- `state_specific_contact`
- `proxy_mapping_used`

## Recommended next visualizations

Useful figures from the current classification stage before structural
integration:

- pipeline flow diagram: downloaded → curated → eligible → classified →
  clean/warning subsets
- grouped bars: mtDNA vs nucDNA, AA-level and NT-level cDAV proportions
- complex-level heatmap of cDAV fractions
- gene-level dot plot of AA-level vs NT-level cDAV proportions
- mtDNA source-provenance overlap summary (`MITOMAP` vs mitochondrial `ClinVar`)
- overlap-aware schematic for `MT-ATP6/MT-ATP8`

Useful figures once the structural stage is refactored:

- per-complex structure panels with mapped cDAVs and uDAVs
- contact-consensus plot across models for high-priority variants
- interface enrichment plot comparing cDAVs and uDAVs
- network view of candidate compensatory partners colored by genome of origin

## Immediate implementation boundary

This refactor should focus on:

- manifest-backed structure selection
- classified Parquet inputs
- multi-chain handling
- multi-model outputs
- compatibility-preserving CSV outputs

It should not yet attempt:

- FoldX integration
- mutagenesis prioritization changes
- new phylogenetic tests beyond making the structural output compatible with
  existing downstream scripts

## Current implementation status

Implemented in the current codebase:

- manifest-backed structure panel defined in:
  - `data/reference/structure_model_manifest.tsv`
- structural mapping now reads:
  - `data/derived/classified/variants_master_classified.parquet`
- structure mapping now writes:
  - `results/structural/dar_structure_map.csv`
  - `results/structural/dar_contacts_cbcb8A.csv`
  - `results/structural/dar_mito_nuc_contacts.csv`
  - `results/structural/structure_model_summary.csv`
- structural partner analysis was updated to load the classified Parquet layer
  and contact rows keyed by `variant_id`

Current practical state:

- the active structure panel is now present locally under `data/structures/`
- the current main limitations are no longer missing CIFs, but:
  - residue anchoring failures in specific genes and complexes
  - incomplete chain coverage for some accessory subunits in some models
  - isoform-only genes whose structural proxy is absent from the selected model

Current verified mapping summary from the refactored script:

- classified rows loaded: `33779`
- mapping rows written: `41208`
- successfully mapped model rows: `11049`
- unique variants with at least one successful structural map: `5562`
- contact rows written: `113442`
- cross-genome contact rows written: `2981`
- unique mito-nuclear interface variants: `400`
- AA-level cDAVs at the mito-nuclear interface: `178`

The current local structure cache now also includes the previously missing:

- `5Z62`
- `8H9U`

Current audit-driven improvements now implemented in
`src/structural/00_map_davs_to_structure.py`:

- nucDNA amino-acid coordinates are remapped into structure space using the
  curated transcript-position maps before residue anchoring
- partial RCSB API chain annotations are now supplemented with local
  sequence-based chain assignment instead of being treated as complete
- mapping rows now carry explicit `status_category`, `aa_coord_method`, and
  `aa_coord_status` fields
- per-variant summaries again record `mapping_support_class`
- a focused structural failure audit is written to:
  - `results/structural/structure_mapping_failure_audit.csv`

## What the audit changed

The structural audit clarified that the earlier mtDNA-dominant mapping result
was not a reliable biological conclusion. It was largely driven by mapper
contract failures on the nucDNA side.

What the audit showed:

- nucDNA mapping was being undercounted because amino-acid coordinates were not
  consistently remapped into the structure-space transcript frame before
  anchoring
- fallback handling for missing frame-specific amino-acid fields was not robust
- partial RCSB API chain annotations were incorrectly treated as complete,
  which hid valid chains in mixed-quality models

What was fixed:

- nucDNA amino-acid coordinates are now translated through the curated
  transcript-position maps before structural anchoring
- chain assignment now supplements partial API annotations with local
  sequence-based assignment rather than stopping early
- failure states are now exported explicitly and can be summarized by gene,
  model, complex, and failure class

What we learned biologically and technically:

- a large part of the previous structural absence for nucDNA was methodological
  rather than biological
- after the audit, nucDNA coverage became substantial in `CI`, `CII`, `CIII`,
  `CIV`, and `CV`, though still uneven
- the remaining structural gaps are now concentrated in a much smaller set of
  interpretable failure classes
- the structural stage is now suitable for denominator-aware analysis, but not
  yet for treating the entire classified set as equally mappable

## Current failure profile

Current residual failure burden from
`results/structural/structure_mapping_failure_audit.csv`:

- `residue_anchoring_failure`: `2315` rows
- `chain_assignment_or_model_gap`: `629` rows
- `isoform_proxy_gap`: `92` rows

Top residual genes by audited failure burden:

- `NDUFV3`: `253` residue-anchoring failures
- `SDHD`: `150` residue-anchoring failures
- `SDHA`: `125` residue-anchoring failures
- `SDHB`: `102` residue-anchoring failures
- `CYC1`: `99` chain/model-gap failures plus `42` residue-anchoring failures
- `SDHC`: `97` residue-anchoring failures
- `UQCRC1`: `94` chain/model-gap failures
- `UQCRC2`: `91` chain/model-gap failures
- `NDUFS3`: `72` residue-anchoring failures
- `NDUFS7`: `70` residue-anchoring failures
- `NDUFAB1`: `69` residue-anchoring failures
- `NDUFS6`: `66` residue-anchoring failures
- `NDUFS4`: `61` residue-anchoring failures
- `COX6A1`: split across `53` chain/model-gap failures and `75`
  residue-anchoring failures
- `COX4I2`: `49` isoform-proxy-gap failures
- `NDUFA9`: `48` residue-anchoring failures
- `NDUFA10`: `48` residue-anchoring failures

Current lowest-coverage classified genes include:

- unassigned accessory genes with no current structural representation in the
  selected human panel, such as `ATP5IF1`, `ATP5ME`, `ATP5MF`, `ATP5MG`,
  `ATP5MJ`, `ATP5MK`, `COX6B2`, `COX7B2`, `COX8C`, `COXFA4`, `COXFA4L2`, and
  `COXFA4L3`
- poorly covered mapped genes such as:
  - `NDUFV3` (`5.62%`)
  - `COX7A1` (`12.5%`)
  - `ATP5MC2` (`22.22%`)
  - `NDUFAB1` (`30.3%`)
  - `COX7A2` (`33.33%`)
  - `ATP5MC3` (`39.13%`)
  - `COX5A` (`50.0%`)
  - `NDUFS6` (`52.17%`)

## Fix plan for the remaining structural failures

The remaining work should be handled in ordered passes, because the failure
classes are now distinct.

### Pass 1. Improve residue anchoring for high-value genes

Priority genes:

- `NDUFV3`
- `SDHA`
- `SDHB`
- `SDHC`
- `SDHD`
- `NDUFS3`
- `NDUFS4`
- `NDUFS6`
- `NDUFS7`
- `NDUFA9`
- `NDUFA10`
- `NDUFAB1`

Planned changes:

- inspect whether failures are concentrated at N-termini, transit-peptide
  boundaries, or short unresolved internal segments
- replace the current uniform `±10` residue anchor window with an adaptive
  strategy:
  - direct match first
  - transcript-remapped coordinate
  - wider bounded search for genes known to show mature-protein offsets
- record the anchoring rule actually used, not only the final success/failure
- distinguish:
  - mature-protein offset rescue
  - unresolved-structure segment
  - genuine sequence-model conflict

Expected outcome:

- reduce the large `CII` burden
- improve the low-coverage `CI` subunits that are currently sequence-compatible
  but still structurally hard to anchor

### Pass 2. Separate model absence from chain-assignment failure

Priority genes:

- `CYC1`
- `UQCRC1`
- `UQCRC2`
- `COX6A1`
- `COX4I1`
- `COX5A`
- `COX6B1`
- `UQCRFS1`
- `UQCRB`
- `UQCRQ`
- `UQCRH`

Planned changes:

- explicitly annotate whether a gene is absent from a model versus present but
  unassigned
- build a per-model expected-subunit table for each selected PDB
- avoid using `5XTH` as if it were an equivalent chain-complete source for all
  CIII/CIV accessory subunits
- prefer per-complex standalone models for chain-complete mapping where
  available, and use respirasome models mainly as validation layers

Expected outcome:

- many current `chain_assignment_or_model_gap` rows will be reclassified from
  generic failure into explicit model-composition limits
- this improves scientific interpretation even when mapping counts do not rise

### Pass 3. Treat isoform-only genes as a formal sensitivity class

Priority genes:

- `COX4I2`
- `COX6A2`
- `COX7A1`

Planned changes:

- keep proxy mapping explicit and non-default
- classify proxy-derived mappings as a separate support tier
- report proxy-eligible versus proxy-unrepresented genes separately in
  structural summaries

Expected outcome:

- prevents isoform-proxy rows from being mixed silently with direct mappings
- makes structural sensitivity analyses cleaner for manuscript figures

### Pass 4. Build a structure-panel eligibility registry

Planned changes:

- add a small reference artifact describing structural eligibility by gene:
  - directly represented in at least one primary model
  - only represented in a validation model
  - only proxy-represented
  - absent from the current human panel
- use this registry to define:
  - structurally_testable_all
  - structurally_testable_direct_only
  - structurally_proxy_only
  - structurally_unrepresented

Expected outcome:

- structural enrichments and partner analyses can use explicit denominators
- unassigned/accessory genes are handled as a documented design limitation

### Pass 5. Rebuild partner analysis on top of audited mapping subsets

Planned changes:

- rerun partner extraction using:
  - all mapped rows
  - direct-only mapped rows
  - multi-model-supported rows
- compare whether partner calls are stable across these subsets

Expected outcome:

- prevents partner results from being dominated by structurally fragile or
  proxy-only mappings

## Biggest unresolved issues now

These are the main unresolved issues for the structural pipeline at this point.

### 1. Residue anchoring is still the dominant technical bottleneck

This is the largest remaining failure class by a wide margin. The current
mapper is much better than before, but a uniform anchoring heuristic is still
too simple for:

- mature-protein offsets
- unresolved loops or termini
- subunits with model-specific residue numbering drift

This is the highest-priority technical issue.

### 2. Structure-panel composition is still uneven across genes

Some genes are genuinely absent from the selected human models, especially in
the `Unassigned` set and among some accessory subunits. Without a formal
eligibility registry, these absences can still be misread as mapping failure
rather than lack of structural representation.

This is the highest-priority interpretive issue.

### 3. Respirasome models are not interchangeable with standalone models

`5XTH` is useful, but it is not a chain-complete or residue-complete surrogate
for all accessory subunits of `CIII` and `CIV`. The pipeline now shows this
more clearly, but the stage contract still needs to encode this explicitly.

This is the highest-priority model-selection issue.

### 4. Isoform proxy use is still biologically weaker than direct representation

Proxy mappings are informative, but they are not equivalent to direct mapping
of the same gene product. This matters especially for `COX4I2`, `COX6A2`, and
`COX7A1`.

This is the highest-priority sensitivity-analysis issue.

### 5. Structural coverage is still not yet a neutral denominator

The structural stage is now far more accurate, but structural mappability is
still strongly gene- and complex-dependent. Any downstream cDAV versus uDAV
claims must be conditioned on structural eligibility and support class.

This is the highest-priority downstream-analysis issue.

## Non-SDH large-offset audit outcome

A focused follow-up audit of the non-`SDH` large-offset candidates showed that
the current evidence does not support broad expansion of the rescue registry
beyond the `SDH` genes.

Main findings:

- most non-`SDH` large-offset candidates occur under `transcript_identity`,
  not transcript-map failure
- many genes show broad and multimodal offset distributions rather than a
  stable single shift
- several candidate genes also have other known caveats:
  - `COX6A1`: chain/model-gap burden in addition to offset candidates
  - `NDUFS6`, `NDUFS7`: prior transcript-model caveats from classification
  - `NDUFA10`: region-specific cross-stage discordance

Current policy decision:

- the anchor exception registry remains restricted to:
  - `SDHA`
  - `SDHB`
  - `SDHC`
  - `SDHD`
- non-`SDH` large-offset signals remain diagnostic only and are not converted
  into automatic mappings

This keeps the structural stage conservative while still preserving the audit
evidence for future manual review.

Current dominant residual failure modes:

- `residue_anchoring_failure` remains concentrated in `CII` (`SDHA`, `SDHB`,
  `SDHC`, `SDHD`) and a subset of `CI` genes such as `NDUFV3`, `NDUFS3`,
  `NDUFS4`, `NDUFS6`, `NDUFS7`, `NDUFA9`, and `NDUFA10`
- `chain_assignment_or_model_gap` remains concentrated in `5XTH` for several
  `CIII` and `CIV` accessory subunits such as `CYC1`, `UQCRC1`, `UQCRC2`,
  `COX6A1`, `COX4I1`, `COX5A`, and `COX6B1`
- `isoform_proxy_gap` remains for isoform-only genes such as `COX4I2`,
  `COX6A2`, and `COX7A1` when the proxy chain is absent from the selected
  model

## Coverage caveat and framework nuance

The structural stage is now multi-model in architecture, but that does not mean
all classified variants are equally well covered structurally.

In practice, a variant is only counted as structurally mapped in a given model
if the current code can:

- identify the relevant chain for the interpreted gene
- align the human reference protein to that chain
- anchor the variant residue position confidently
- verify that the expected human reference amino acid matches the structure

This means structural support is limited by both biology and model tractability.

Important consequence:

- successful mappings remain concentrated in the subset of models and subunits
  that support residue anchoring robustly
- this is especially true for mtDNA core membrane subunits and some well-resolved
  complexes
- many nucDNA rows remain structurally unmapped not because they lack biological
  relevance, but because the available structures do not yet support confident
  residue-level mapping for those specific subunits or positions

Therefore, structural results must be interpreted as a constrained evidence
layer:

- absence of structural support is not equivalent to absence of a structural
  mechanism
- structural enrichment claims must be conditioned on structural mappability
- cDAV versus uDAV structural comparisons should use explicit structure-coverage
  denominators

This limitation is now visible and quantifiable in the output tables, which is a
substantial improvement over the previous single-model implementation where this
coverage bias was largely hidden.
