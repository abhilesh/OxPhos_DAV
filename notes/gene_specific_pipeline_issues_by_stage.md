# Gene-Specific Pipeline Issues by Stage

This note records the specific genes that have required special attention across
the download/reference, curation/classification, and structural mapping phases.

The goal is to distinguish:

- genes that had genuine cross-stage model-mapping problems
- genes whose earlier curation/classification issues were solved, but that still
  remain structurally difficult
- genes whose current structural failures are mostly due to structure-panel
  composition rather than upstream curation errors

## 1. Short answer

Some overlap is real, but it is not one-to-one.

There are three broad patterns:

1. Certain genes had transcript-model or coordinate-model issues in
   curation/classification and still show structural mapping difficulty.
2. Some genes were repaired successfully at the curation/classification stage
   and are now only structurally difficult because of anchoring or model
   coverage.
3. Some genes were not major curation problems at all, but remain difficult in
   structure mapping because the selected human PDB panel lacks their chain or
   only represents an isoform proxy.

## 2. Genes with cross-stage overlap

These genes have evidence of problems in more than one phase of the pipeline.

### 2.1 `NDUFA10`

Stages affected:

- classification
- structural mapping

Problems observed:

- residual classification unresolved rows were split between:
  - `POSITION_NOT_IN_ENST` around residue `217`
  - `TRANSCRIPT_MISMATCH` around residues `337-350`
- one classified row remained as `REF_ALLELE_MISMATCH`:
  - `ClinVar:897590:NDUFA10` (`R337H`)
- structural mapping still shows substantial `residue_anchoring_failure`

Interpretation:

- this is not a single global failure
- `NDUFA10` appears to contain region-specific transcript/model discordance in
  classification, plus downstream residue-anchoring difficulty in structure
  space

Current status:

- curation/classification is usable, but `NDUFA10` retains warning-bearing and
  unresolved edge cases
- structural mapping is partially recovered but still incomplete

### 2.2 `NDUFS6`

Stages affected:

- classification
- structural mapping

Problems observed:

- unresolved classification rows were concentrated in
  `GENOMIC_POS_NOT_IN_ENST`
- transcript-map quality was below the threshold for transcript-first handling:
  - identity `0.871`
  - coverage `1.0`
- exception registry marks it as a transcript-model-incompatible,
  low-identity genomic-rescue gene
- structural mapping still shows `residue_anchoring_failure`

Interpretation:

- the upstream problem is transcript/ENST concordance
- the structural problem is now residue anchoring after coordinate rescue

Current status:

- still a genuine cross-stage difficult gene

### 2.3 `NDUFS7`

Stages affected:

- classification
- structural mapping

Problems observed:

- unresolved classification rows were concentrated in
  `GENOMIC_POS_NOT_IN_ENST`
- transcript-map quality was modest:
  - identity `0.887`
  - coverage `0.953`
- exception registry marks it as transcript-model-incompatible
- structural mapping still shows `residue_anchoring_failure`

Interpretation:

- same general pattern as `NDUFS6`, though less severe

Current status:

- still a real cross-stage problem gene

### 2.4 `COX5A`

Stages affected:

- classification
- structural mapping

Problems observed:

- residual classification unresolved set included `1`
  `GENOMIC_POS_NOT_IN_ENST` row
- transcript-map quality was relatively low:
  - identity `0.787`
  - coverage `0.960`
- structural mapping shows both:
  - `chain_assignment_or_model_gap`
  - `residue_anchoring_failure`

Interpretation:

- `COX5A` is a mixed case:
  - upstream transcript-model incompatibility exists
  - structural coverage is also limited in the current PDB panel

Current status:

- still a cross-stage caution gene

### 2.5 `COX4I2`

Stages affected:

- transcript-map build / earlier curation review
- structural mapping

Problems observed:

- it appeared among previously problematic genes during transcript-map repair
- after rebuild, it was considered to have a biologically coherent map
- structural mapping still shows:
  - `isoform_proxy_gap`
  - some `residue_anchoring_failure`
- current structural mapper uses a proxy:
  - `COX4I2 -> COX4I1`

Interpretation:

- the earlier sequence/transcript issue was largely solved
- the remaining problem is primarily structural-model representation, not
  curation correctness

Current status:

- no longer a major classification problem
- still a structural proxy problem

## 3. Genes that were largely fixed upstream but remain structurally difficult

These genes overlap conceptually with earlier mapping concerns, but their
current failures are mainly structural rather than curation/classification
failures.

### 3.1 `SDHA`, `SDHC`, `SDHD`

Stages affected:

- transcript-map build review
- structural mapping

Problems observed:

- these genes were explicitly noted as previously problematic before the
  transcript-map rebuild
- after rebuilding transcript maps, they were reported as having biologically
  coherent maps
- structurally, they remain prominent in `residue_anchoring_failure`

Interpretation:

- the curation/classification layer is no longer the main issue
- the present limitation is structural anchoring, likely reflecting mature
  protein numbering, model trimming, or local sequence/model mismatch

Current status:

- upstream issue mostly solved
- structural residue anchoring still needs targeted rescue

### 3.2 `NDUFA9`

Stages affected:

- earlier mapping review
- structural mapping

Problems observed:

- it was part of the broader nuclear mapping/correction effort
- no strong residual unresolved classification burden is currently recorded
- structural mapping still shows substantial `residue_anchoring_failure`

Interpretation:

- current failure is mainly structural, not classificatory

### 3.3 `NDUFV1`, `NDUFS1`, `NDUFA6`

Stages affected:

- transcript-map build review
- structural mapping for some related subunits

Problems observed:

- these were cited among genes whose transcript maps became biologically
  coherent after the rebuild
- they are not prominent residual unresolved classification genes now
- some related Complex I anchoring issues persist structurally, though these
  exact genes are not the worst current structural offenders

Interpretation:

- they demonstrate that early mapping problems and current structural failures
  are not always the same thing

## 4. Genes that are mainly structural-problem genes now

These genes are not prominent unresolved curation/classification problems in the
current pipeline state, but remain structurally limited.

### 4.1 `SDHB`

Stages affected:

- structural mapping

Problems observed:

- prominent `residue_anchoring_failure`
- not highlighted as a major residual classification problem

Interpretation:

- mostly a structural anchoring problem in the current state

### 4.2 `NDUFV3`

Stages affected:

- structural mapping

Problems observed:

- one of the strongest remaining `residue_anchoring_failure` signals
- not part of the residual unresolved classification set summarized in the
  classification notes

Interpretation:

- this is mainly a structure-space anchoring problem, not an unresolved
  curation issue

### 4.3 `NDUFS3`, `NDUFS4`

Stages affected:

- structural mapping

Problems observed:

- substantial `residue_anchoring_failure`
- not prominent residual unresolved classification genes in the current notes

Interpretation:

- current issue is mainly structural anchoring

### 4.4 `CYC1`, `UQCRC1`, `UQCRC2`, `COX6A1`, `COX4I1`, `COX6B1`

Stages affected:

- structural mapping

Problems observed:

- major `chain_assignment_or_model_gap` burden, especially in `5XTH`
- many of these genes are accessory/supercomplex-context subunits where the
  selected structure panel does not consistently expose a directly assignable
  chain in every model

Interpretation:

- these are mostly structure-panel/chain-coverage problems
- they should not be interpreted as upstream curation failures

### 4.5 `COX6A2`, `COX7A1`

Stages affected:

- structural mapping

Problems observed:

- prominent `isoform_proxy_gap`
- current structural mapper relies on proxies:
  - `COX6A2 -> COX6A1`
  - `COX7A1 -> COX7A2`

Interpretation:

- these are primarily isoform-representation problems in the structure panel

## 4A. Non-SDH large-offset candidates audited but not rescued

The non-`SDH` large-offset audit found several genes with repeated
large-offset-candidate patterns, but these were not added to the rescue
registry because the offset distributions were too broad or multimodal to be
treated as a clean mature-protein rule.

Most relevant genes:

- `COX6A1`
- `NDUFAB1`
- `NDUFS3`
- `NDUFS6`
- `NDUFS7`
- `NDUFA9`
- `NDUFA10`
- `ATP5F1A`
- `ATP5F1D`
- `ATP5MC2`
- `ATP5MC3`

Interpretation:

- these genes remain structurally interesting and should stay on the audit list
- but they are not yet strong enough candidates for automatic extended-offset
  rescue

## 5. Practical summary by gene set

### 5.1 True cross-stage caution genes

- `NDUFA10`
- `NDUFS6`
- `NDUFS7`
- `COX5A`
- `COX4I2`

These are the strongest candidates for explicit sensitivity analyses in later
mechanistic work.

### 5.2 Upstream-fixed but structurally difficult genes

- `SDHA`
- `SDHC`
- `SDHD`
- `NDUFA9`

These should not be treated as unresolved curation failures anymore.

### 5.3 Predominantly structural-panel genes

- `SDHB`
- `NDUFV3`
- `NDUFS3`
- `NDUFS4`
- `CYC1`
- `UQCRC1`
- `UQCRC2`
- `COX6A1`
- `COX4I1`
- `COX6B1`
- `COX6A2`
- `COX7A1`

These are best handled by improving structural anchoring rules, expanding the
structure panel, or documenting model-limited coverage.

## 6. Implication for downstream analysis

The structural residual set should not be interpreted as a replay of the old
curation/classification problem set.

The current state is more refined:

- curation/classification problems are now mostly narrowed to a small number of
  transcript-model discordance genes
- structural failures are now dominated by:
  - residue anchoring in otherwise valid genes
  - incomplete chain representation in certain models
  - isoform proxy limitations

This means downstream structural enrichment and partner analyses should:

- retain the gene-specific caveats above
- use structural coverage denominators explicitly
- consider sensitivity analyses excluding the true cross-stage caution genes
  from high-confidence structural claims
