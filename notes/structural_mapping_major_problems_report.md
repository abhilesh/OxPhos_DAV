# Structural Mapping Stage: Major Problems Report

This note summarizes the major problems encountered in the structural mapping
stage, which genes are affected, how they are affected, and what the current
policy is.

It reflects the current conservative structural mapping state after:

- structure-source validation
- transcript-reconciliation audit
- anchoring refactor
- conservative anchor-exception registry gating

## 1. Current structural bottleneck

The structural stage is the current bottleneck in the analytical pipeline.

The problem is no longer:

- raw data download
- curated master-table generation
- cDAV/uDAV classification

The problem is now specifically:

- deciding which difficult residue-to-structure mappings are scientifically
  defensible
- distinguishing direct mappings from rescued mappings
- avoiding overcalling structure support for difficult nucDNA genes

## 2. Main problem classes in the structural stage

The structural failures now fall into a small number of interpretable classes.

### 2.1 `chain_assignment_or_model_gap`

Meaning:

- the structural model does not provide a confidently assignable chain for the
  interpreted gene
- or the subunit is absent / poorly represented in the selected model

Typical causes:

- accessory subunits absent from some respirasome-context models
- incomplete chain annotation in the chosen structure
- structural panel not equally suitable for all subunits in a complex

Main affected genes:

- `CYC1`
- `UQCRC1`
- `UQCRC2`
- `COX6A1`
- `COX4I1`
- `COX5A`
- `COX6B1`
- related `CIII`/`CIV` accessory genes

How they are affected:

- variants in these genes may be biologically relevant
- but the current model set does not consistently give a directly usable chain
- these rows remain structurally uninformative, not biologically negative

### 2.2 `isoform_proxy_gap`

Meaning:

- the direct gene product is not represented in the structure
- the mapper can only reason through a proxy isoform, and that proxy is itself
  absent or not confidently usable in the current model

Main affected genes:

- `COX4I2`
- `COX6A2`
- `COX7A1`

How they are affected:

- structural support for these genes is weaker than for directly represented
  genes
- they should not be mixed with direct structural evidence without explicit
  labeling

### 2.3 `residue_anchoring_failure`

Meaning:

- the chain exists
- the transcript reconciliation may succeed
- but the residue still cannot be placed confidently in structure space within
  the current conservative mapping policy

Typical causes:

- unresolved structure segments
- local numbering drift
- processed-protein offsets that are not yet approved for rescue
- local sequence-model conflicts

### 2.4 `mature_offset_candidate`

Meaning:

- a plausible larger offset exists beyond the conservative default anchor
  window
- but the gene is not currently allowed to use extended-window rescue

This is now a diagnostic class, not a mapped class.

Meaning for interpretation:

- these rows are not considered structurally mapped
- but they are important evidence for future gene-specific manual review

### 2.5 Cross-stage coordinate-system mismatch

Meaning:

for nucDNA rows, three spaces must be reconciled:

1. ClinVar / curated transcript consequence space
2. TOGA ENST human alignment space
3. structural chain residue space

The current pipeline now records this explicitly, but this remains a conceptual
source of difficulty for some genes.

## 3. Genes affected and how they are affected

### 3.1 Genes mainly affected by chain / model gaps

- `CYC1`
- `UQCRC1`
- `UQCRC2`
- `COX6A1`
- `COX4I1`
- `COX5A`
- `COX6B1`

Effect:

- these genes are often limited by model composition rather than transcript
  mapping
- the current structures do not consistently provide a chain that can be used
  as a clean residue-level mapping target

### 3.2 Genes mainly affected by isoform / proxy limitations

- `COX4I2`
- `COX6A2`
- `COX7A1`

Effect:

- the structure panel mainly represents related isoforms
- direct mapping is weak or absent
- these should remain secondary / sensitivity-class structural rows

### 3.3 Genes with the strongest remaining residue-anchoring burden

High-burden non-`SDH` structural genes:

- `COX6A1`
- `NDUFAB1`
- `NDUFS3`
- `NDUFS4`
- `NDUFS6`
- `NDUFS7`
- `NDUFA9`
- `NDUFA10`
- `NDUFV3`
- `ATP5F1A`
- `ATP5F1B`
- `ATP5F1D`
- `ATP5MC2`
- `ATP5MC3`

Effect:

- these genes repeatedly show large-offset or local anchoring problems
- most are *not* currently approved for rescue
- their structural evidence remains conservative by design

### 3.4 Genes currently approved for conservative offset rescue

- `SDHA`
- `SDHB`
- `SDHC`
- `SDHD`

Effect:

- these are the only genes currently allowed to convert audited extended-window
  offset candidates into structural mappings
- this is done through the anchor exception registry, not through silent global
  rescue logic

## 4. What the non-SDH large-offset audit found

The most important negative result was:

- the non-`SDH` large-offset candidates are not clean enough yet for automatic
  rescue

Why:

- most occur under `transcript_identity`, so they are not simple transcript-map
  failures
- many show diffuse or multimodal offsets rather than one tight, stable shift
- some overlap with already known cross-stage caution genes

Examples:

- `COX6A1`: broad positive offsets, also affected by chain/model gaps
- `NDUFAB1`: broad positive offsets, not tightly clustered
- `NDUFS3`: mixed positive and negative offsets
- `NDUFA9`: mixed positive and negative offsets
- `NDUFA10`: multiple offset modes consistent with regional discordance

Conclusion:

- do not expand the rescue registry beyond `SDH` at this stage

## 5. Current policy state

The current structural mapper uses a conservative three-layer interpretation:

1. direct/default-window mappings
2. explicit `SDH` rescued-offset mappings
3. unresolved large-offset candidates retained for future review

This is the current scientifically defensible state.

## 6. Current take-home message

The structural stage is not failing globally anymore.

It is now limited in three targeted ways:

- some genes are not well represented by the current structure panel
- some genes are isoform/proxy-limited
- some genes show unresolved large-offset patterns that still need manual
  biological review before rescue

That is a much narrower and more interpretable problem than the earlier state,
but it is still the main place where the pipeline is not yet final.
