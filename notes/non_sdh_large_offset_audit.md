# Non-SDH Large-Offset Audit

This note records the focused audit of non-`SDH` large-offset candidates after
the structural anchoring refactor and the introduction of the conservative
anchor exception registry.

## Why this audit was needed

After the staged anchoring refactor, the structural mapper began separating:

- true unresolved anchoring failures
- large-offset candidates
- explicit extended-offset rescues

An initial permissive run showed that many genes outside the `SDH` set could be
rescued by a very wide anchor window. That increased coverage, but it was not
clear that those rescues were biologically justified.

This audit was therefore used to decide whether the anchor exception registry
should be expanded beyond:

- `SDHA`
- `SDHB`
- `SDHC`
- `SDHD`

## Main result

The current evidence does **not** support broad expansion of the anchor rescue
registry beyond the `SDH` genes.

The non-`SDH` candidates are mostly characterized by:

- diffuse and multimodal offset distributions
- large positive shifts that vary widely within the same gene
- mixtures of structural-accessory and isoform-sensitive subunits
- in some cases, prior curation/classification caveats already identified

This means most current non-`SDH` large-offset signals are better treated as:

- diagnostic flags
- gene-specific audit targets

rather than as approved rescue rules.

## Audit observations

### 1. Most non-SDH large-offset candidates are not transcript-map failures

For most of the highest-burden genes, the large-offset candidates occur under:

- `transcript_identity`

rather than:

- `transcript_remapped`
- `position_not_in_enst`

This means the dominant problem is not simply ClinVar-to-TOGA transcript
reconciliation failure. It is more likely to reflect structure-space issues
such as:

- mature-protein processing
- chain/model numbering drift
- unresolved structure segments
- local structure/reference incompatibility

Examples:

- `COX6A1`: all current large-offset candidates are `transcript_identity`
- `NDUFAB1`: all current large-offset candidates are `transcript_identity`
- `NDUFS3`: all current large-offset candidates are `transcript_identity`
- `NDUFS6`: all current large-offset candidates are `transcript_identity`
- `NDUFA9`: all current large-offset candidates are `transcript_identity`

### 2. Many non-SDH genes show diffuse offset patterns rather than one stable shift

This is the main reason they were **not** added to the registry.

Examples:

#### `COX6A1`

- top offsets include `+27`, `+33`, `+19`, `+21`, `+36`, `+41`
- these are all positive, but not concentrated tightly enough to justify a
  single gene-level mature-offset rule yet

#### `NDUFAB1`

- top offsets include `+38`, `+58`, `+16`, `+41`, `+51`
- pattern is broad rather than sharply clustered

#### `NDUFS3`

- offsets include `+36`, `+49`, `+22`, but also negative values like `-12`
  and `-14`
- this argues against a simple mature-protein offset explanation

#### `NDUFA9`

- offsets include `+34`, `+31`, `+36`, `+42`, `+51`, but also `-11`
- again, not a clean single-offset case

#### `NDUFA10`

- offsets include both `+13` and `+39` as strong modes
- this fits the earlier interpretation that `NDUFA10` has region-specific
  rather than global coordinate discordance

### 3. Some genes are still plausible future rescue candidates, but not yet

These genes may warrant a second-round focused audit:

- `COX6A1`
- `NDUFAB1`
- `NDUFS6`
- `NDUFS7`
- `NDUFA10`

However, the current evidence is still insufficient for automatic rescue.

Reasons:

- `COX6A1` also has substantial `chain_assignment_or_model_gap` burden, so its
  structural issue is not purely an offset problem
- `NDUFS6` and `NDUFS7` already had transcript-model caveats in the
  curation/classification phase
- `NDUFA10` is already a known cross-stage caution gene with regional
  discordance
- `NDUFAB1` has broad, poorly concentrated offsets

### 4. Some candidate groups are likely poor rescue targets

These should remain unresolved rather than being force-rescued:

- `ATP5F1A`, `ATP5F1B`, `ATP5F1D`
- `ATP5MC2`, `ATP5MC3`
- `COX4I2`, `COX6A2`, `COX7A1`
- `COX5A`

Reasons:

- the offset patterns are broad and inconsistent
- several belong to isoform-sensitive or proxy-sensitive subunits
- some are structurally accessory and already affected by model-composition
  issues

## Registry decision

Current decision:

- keep the registry restricted to `SDHA`, `SDHB`, `SDHC`, `SDHD`
- do **not** add non-`SDH` genes automatically at this stage

This is a deliberate conservative choice.

## Implication for the current pipeline

The structural outputs should now be interpreted as three layers:

1. direct / default-window mappings
2. explicit `SDH` extended-offset rescues
3. unresolved large-offset candidates for future audit

That is a scientifically safer state than broad rescue across all high-offset
genes.

## Recommended next step

If the registry is to be expanded later, the next additions should only occur
after gene-specific manual review with local sequence windows and model
inspection, starting with:

- `COX6A1`
- `NDUFAB1`
- `NDUFS6`
- `NDUFS7`
- `NDUFA10`
