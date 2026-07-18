# OxPhos cDAV Pipeline — Development Notes & Issue Log

A chronological record of every bug, mismatch, and design issue encountered during development, with root cause and resolution. Intended as a reference for debugging and as institutional memory for future contributors.

---

## Table of Contents

1. [Pagel HPC Pipeline Issues](#1-pagel-hpc-pipeline-issues)
2. [Phylogenetic Timing Issues](#2-phylogenetic-timing-issues)
3. [Structural / Compensatory Partner Issues](#3-structural--compensatory-partner-issues)
4. [Mutagenesis Pipeline Issues (New)](#4-mutagenesis-pipeline-issues-new)

---

## 1. Pagel HPC Pipeline Issues

---

### Issue 1.1 — `ROOT` undefined in `prepare_pagel_hpc.py`

**File:** `src/phylo/prepare_pagel_hpc.py`

**Symptom:** `NameError: name 'ROOT' is not defined` when running the script. All path constants (e.g. `TREE_FILE`, `OUTPUT_DIR`) were defined using `ROOT`, but `ROOT` itself was never assigned.

**Root cause:** The line `ROOT = Path(__file__).resolve().parents[2]` was missing. Every other script in the codebase has this at the top of the paths section; it was accidentally omitted during initial writing of this file.

**Fix:**
```python
# Added at top of paths section:
ROOT = Path(__file__).resolve().parents[2]
```

**Lesson:** When copying path boilerplate across scripts, always verify `ROOT` is defined before the first path constant that uses it.

---

### Issue 1.2 — All Pagel results returned `NA` (First HPC run)

**File:** `src/phylo/pagel_discrete.R`

**Symptom:** Merging HPC Pagel output (`merge_pagel_results.py`) produced 0 finite p-values. Every pair had `p_value = NA`.

**Investigation:** Ran a single chunk manually inside Docker:
```bash
conda run -n pagel Rscript src/phylo/pagel_discrete.R chunks/chunk_0001.tsv out.tsv
```
R printed repeated warnings:
```
Warning: force.ultrametric was applied because the tree is not ultrametric.
```
But crucially, `fitPagel()` was returning `NA` for all pairs.

**Root cause 1 — `is.ultrametric()` returning FALSE:** The MamPhy timetree has floating-point rounding in its branch lengths that accumulates to make the root-to-tip path lengths differ by small amounts (< 1e-6). R's `ape::is.ultrametric()` uses a strict threshold and returns `FALSE`. `phytools::fitPagel()` internally checks for ultrametricity and silently returns `NA` when the check fails.

**Root cause 2 — timeout too short:** `PAGEL_TIMEOUT_SEC = 5` (initial value) was far too low. Even well-behaved `fitPagel` calls take ~10–60 seconds on HPC nodes. The `setTimeLimit()` exception was being triggered for essentially every pair.

**Fix 1 — force ultrametric on load:**
```r
get_tree <- function(tree_file) {
  tree <- read.tree(tree_file)
  if (!is.ultrametric(tree)) {
    tree <- suppressMessages(force.ultrametric(tree, method = "extend"))
  }
  tree
}
```
The `suppressMessages()` wrapper silences the repeated "force.ultrametric" notes in logs, which had been flooding the SLURM output and obscuring the real error.

**Fix 2 — raise timeout:** Set `PAGEL_TIMEOUT_SEC <- 60L`.

---

### Issue 1.3 — All Pagel results still `NA` after timeout fix (Second HPC run)

**File:** `src/phylo/pagel_discrete.R`

**Symptom:** After raising the timeout to 60 seconds, a second full HPC run still returned all `NA`. The SLURM logs showed each chunk completing in ~2 hours — which matched the runtime limit, not the number of pairs.

**Root cause:** HPC nodes run at different CPU speeds than the local development machine. The 60s timeout was calibrated on local hardware where `fitPagel` takes ~54s per pair. On the HPC, pairs took ~64s on average. Every single pair was timing out.

**Evidence from a test chunk:**
```
Chunk: chunks/chunk_0001.tsv
...Average: 64.2 sec/pair
...Total: 3 hrs (CANCELLED AT TIME LIMIT)
```

**Fix:** Set `PAGEL_TIMEOUT_SEC <- 300L` (5× the observed HPC time per pair). This allows normal pairs to complete comfortably while still catching genuinely hung `fitPagel` calls that never converge.

**Key tension:** The timeout must simultaneously be (a) long enough that normal pairs complete and (b) short enough that pathologically hung calls don't block the batch. The 300s value works because real non-convergence hangs indefinitely, not at 65s.

---

### Issue 1.4 — `chunk_0001` running for 3+ hours without finishing

**File:** `src/phylo/pagel_discrete.R`

**Symptom:** After removing the timeout entirely as a test, chunk_0001 did not finish even after 3 hours (wall time). The job was eventually cancelled by SLURM due to the 3-hour time limit.

**Root cause:** Without any timeout, certain `fitPagel` calls enter an infinite loop during likelihood optimisation. This is a known issue with discrete character correlation models on certain degenerate binary character combinations (e.g., all species in state 0/0 or when one character has no variation). `fitPagel` uses numerical optimisation (`nlm`) that can wander without converging.

**Fix:** Keep the timeout at 300s but implement checkpointing so that cancelled jobs can resume rather than restart from the beginning.

**Checkpointing design:**
- Output file opened in append mode (`open = "a"`) so completed results survive cancellation
- `read_done()` function reads the output file at startup and returns completed `pair_id`s
- `run_batch()` skips pairs already in `done_ids` before calling `run_pagel()`
- Flush after every pair (`flush(output_con)`) to ensure disk writes are not buffered

```r
done_ids <- read_done(args[2])
out_con  <- file(args[2], open = "a")   # append; prior results intact
run_batch(manifest, out_con, done_ids = done_ids)
```

---

### Issue 1.5 — `SyntaxError: invalid syntax` at line 615 of `02_phylogenetic_timing.py`

**File:** `src/phylo/02_phylogenetic_timing.py`

**Symptom:**
```
File "/app/src/phylo/02_phylogenetic_timing.py", line 615
    for cat
           ^
SyntaxError: invalid syntax
```

**Investigation:** Reading the file showed the syntax was correct on disk. The error appeared intermittently.

**Root cause:** Docker image caching. The container had loaded a stale version of the file from an intermediate image layer. The file had been edited mid-run but the container was using a cached copy from before the edit.

**Fix:** Simply re-running the container (not `--build`, just a fresh `docker run`) picked up the current file. No actual syntax error existed.

**Lesson:** When a syntax error is reported at a specific line but the code looks correct in the editor, always verify the file is being mounted correctly (`-v $(pwd):/app`) and that Docker isn't serving a cached layer.

---

### Issue 1.6 — `force.ultrametric` repeated warning flooding SLURM logs

**File:** `src/phylo/pagel_discrete.R`

**Symptom:** Every SLURM log file contained hundreds of lines like:
```
Note: your tree is not ultrametric. force.ultrametric was used...
```
These obscured real errors and made it impossible to spot actual failures in the log.

**Root cause:** `force.ultrametric()` prints a note via `message()` by default. Because the tree is loaded and cached per chunk (many pairs per chunk), but the message is suppressed only at the call site, some code paths were not suppressed.

**Fix:** Wrap the call in `suppressMessages()` at the single load point in `get_tree()`:
```r
tree <- suppressMessages(force.ultrametric(tree, method = "extend"))
```
This is safe because the note is purely informational — we know the tree is not perfectly ultrametric and we are deliberately fixing it.

---

## 2. Phylogenetic Timing Issues

---

### Issue 2.1 — Only 458 of 6,738 pairs had age estimates

**File:** `src/phylo/02_phylogenetic_timing.py` — `node_age_mya()`

**Symptom:** After running timing annotations, `dar_origin_age_mya` was `None` for 6,280 out of 6,738 pairs (93.2%). Only internal-node origin events had age values.

**Root cause:** `node_age_mya()` required the `species_set` argument to contain ≥2 species to find an MRCA:
```python
# Original code — silently returned None for single-species sets
if len(present) < 2:
    return None
mrca = tree.common_ancestor(...)
```

Leaf-origin events (DARs arising on a terminal branch leading to a single extant species) produce a `dar_spp` set with exactly one species. These represent ~90% of all cDAR origins. The function correctly identified them as valid events but returned `None` for their age because it couldn't compute an MRCA with a single tip.

**The correct age for leaf-origin events is `0.0 Mya` (present day):** A DAR arising on the terminal branch leading to species X means it arose some time after the divergence of that species' lineage — we estimate it as occurring in "the present" relative to the phylogenetic timescale of 165 Mya.

**Fix:**
```python
def node_age_mya(tree, species_set: set[str]) -> float | None:
    if tree is None or not species_set:
        return None
    tip_names = {c.name for c in tree.get_terminals()}
    present = [s for s in species_set if s in tip_names]
    if not present:
        return None
    if len(present) == 1:
        return 0.0   # leaf-origin event: occurred on terminal branch = present day
    try:
        mrca = tree.common_ancestor(*[{"name": s} for s in present])
        root_to_mrca = tree.distance(tree.root, mrca)
        root_to_leaf = tree.distance(tree.root, tree.get_terminals()[0])
        return round(root_to_leaf - root_to_mrca, 2)
    except Exception:
        return None
```

**After fix:** Age coverage went from 458 rows to 6,544 rows (6,086 at 0.0 Mya, 458 at true ancestral ages).

---

### Issue 2.2 — `contact_first` timing was questioned as an artifact

**Context:** After reporting that 18.1% of pairs have `contact_first` timing (contact partner evolved before the DAR), a concern was raised: for leaf-origin events (single-species DARs), does the `contact_first` logic correctly identify the immediate parent branch, or does it spuriously flag all ancestor branches?

**Investigation:** The logic in `timing_for_origin()` checks:
```python
# contact_first if: contact partner changed on a branch where:
# dar_spp ⊆ parent_spp AND dar_spp ⊄ child_spp
```

For a single-species set `{X}`, this translates to:
- `{X} ⊆ parent_spp` — satisfied only if X is a descendant of the parent node
- `NOT {X} ⊆ child_spp` — violated if X is also a descendant of the child node

For ancestor branches above the immediate parent (e.g., the grandparent → parent branch), `X` is still a descendant of the child node (= the immediate parent), so `{X} ⊆ child_spp` is TRUE, meaning `NOT {X} ⊆ child_spp` is FALSE. The condition fails.

Only for the branch immediately above the terminal (immediate parent → X), `{X} ⊄ child_spp` is TRUE (because `child_spp` is `{X}` itself and the terminal branch has no further subtree).

**Conclusion:** The logic is correct. The `contact_first` result is not an artifact. The subset math uniquely identifies the immediate parent branch even for single-species origin events.

---

### Issue 2.3 — `physicochemical_type` all `"unknown"` in `compensatory_partners.csv`

**Symptom:** After regenerating `timing_annotations.csv` (following the leaf-origin age fix), all 693 pairs in `compensatory_partners.csv` had `physicochemical_type = unknown`.

**Investigation:** The `physicochemical_type` column is computed in `02_phylogenetic_timing.py` but is propagated to `compensatory_partners.csv` by a downstream join. The regenerated `timing_annotations.csv` had the correct values, but `compensatory_partners.csv` had not been re-generated since the previous (pre-fix) run of the timing script.

**Root cause:** The join between `timing_annotations.csv` and `compensatory_partners.csv` had not been re-run after updating timing annotations. Scripts 01 (find partners) and 02 (timing) needed to be run in sequence, but only script 02 was re-run.

**Fix:** Re-ran the full join script, which merged the updated timing fields back into `compensatory_partners.csv`. After the re-join, all 693 pairs had correct `physicochemical_type` values.

**Lesson:** When re-running a script that feeds into a downstream join, always re-run the join step too. Stale intermediate files are a common source of "wrong values that look right."

---

## 3. Structural / Compensatory Partner Issues

---

### Issue 3.1 — Fisher's test treats species as independent (phylogenetic pseudoreplication)

**Context:** The Fisher's exact test in `01_find_compensating_partners.py` counts species with and without the alternative amino acid at the contact position, in cDAR vs background species.

**Problem:** Species are not independent observations. A clade of 20 Muridae that all share a cDAR via common descent contributes 20 "observations" when biologically it is one substitution event. This inflates test statistics and produces false positives.

**Intended fix (Pagel's discrete):** `fitPagel()` in R accounts for phylogenetic non-independence by modelling character evolution on the tree. It is the phylogenetically valid test.

**Resolution in analysis:** Pagel's test and the branch co-occurrence test (both phylogenetically valid) were implemented as the primary significance calls. Fisher's test is retained as a fallback for low-power cases (< 20 species) where phylogenetic tests cannot run, and is explicitly labelled as not controlling for phylogenetic non-independence.

**Significance threshold:** `(pagel_fdr ≤ 0.10 OR branch_cooccur_fdr ≤ 0.10) AND low_power == False`. Fisher's FDR is used only when both phylogenetic tests are unavailable.

---

### Issue 3.2 — Isoform proxy residue position offset

**Context:** Some tissue-specific OxPhos subunits have no TOGA alignment (e.g. COX4I2, the testis-specific isoform). These are mapped to their ubiquitous paralogues (COX4I1) via a `ISOFORM_PROXY` dictionary in `00_map_davs_to_structure.py`.

**Problem:** The proxy isoform may differ in its N-terminal sequence and exact residue register from the disease-isoform. A DAR at position 59 in COX4I2 may not align to position 59 in COX4I1.

**Resolution:** A `±10` AA anchor window local alignment was implemented to find the matching residue in the proxy. The `status` field in `dar_structure_map.csv` records `isoform_offset_+N` or `isoform_offset_-N` to flag proxy-mapped residues. These should be treated with caution in downstream FoldX and MI analyses.

**Isoform proxies used:**
- COX4I2 → COX4I1
- COX6A2 → COX6A1
- COX7A1 → COX7A2
- ATP5MC2/3 → ATP5MC1

---

## 4. Mutagenesis Pipeline Issues (New)

---

### Issue 4.1 — Zero PDB IDs joined: `ann_id` format mismatch

**File:** `src/mutagenesis/00_prioritize_pairs.py`

**Symptom:** After implementing the join between `compensatory_partners.csv` and `dar_contacts_cbcb8A.csv`, zero rows had `pdb_id` populated. The script reported "Loading 0 CIF structures" and all Cβ-Cβ distances were `NaN`.

**Initial join key used:**
```python
join_keys = ["ann_id", "dar_gene", "dar_aa_coord", "contact_gene", "contact_refseq_pos"]
```

**Investigation:** Printed the `ann_id` formats from each file:

```
compensatory_partners.csv ann_id:
  1029699, 1032977, 1037130, ...    ← ClinVar numeric IDs

dar_contacts_cbcb8A.csv ann_id:
  m.3688G>A, m.13276A>G, ...        ← HGVS mtDNA notation
  729, 1029699, ...                  ← partial overlap for nucDNA
```

**Root cause:** The two files use incompatible `ann_id` formats:
- `compensatory_partners.csv` uses raw ClinVar IDs (integers) for nucDNA and MITOMAP-format strings for mtDNA
- `dar_contacts_cbcb8A.csv` uses HGVS notation for mtDNA (`m.XXXXN>N`) and a different numeric scheme for nucDNA

Checking the overlap on 5,000 rows of the contacts file showed only 49 matching IDs — far fewer than the 693 expected.

**Fix:** Changed the join key to physical residue positions, which are unambiguous and consistent across both files:
```python
join_keys = ["dar_gene", "dar_aa_coord", "dar_alt_aa", "contact_gene", "contact_refseq_pos"]
```
The `ann_id` column is no longer used for the structural coordinate join.

**After fix:** 696 rows matched (3 pairs have two structural contexts in different PDB structures), all 5 CIF structures loaded successfully, 693/693 distances computed.

---

### Issue 4.2 — Column name mismatch: `dar_locus` vs `dar_gene`

**Files:** `dar_contacts_cbcb8A.csv` vs `compensatory_partners.csv`

**Symptom:** After fixing the `ann_id` join (Issue 4.1), the join still returned zero rows because `dar_gene` was not found in `dar_contacts_cbcb8A.csv`.

**Root cause:** The contacts file uses `dar_locus` as the column name for the disease gene; the partners file uses `dar_gene`. Both contain identical values (e.g., `"MT-ATP6"`, `"NDUFS1"`), just under different names.

**Fix:** Rename before join:
```python
if "dar_locus" in contacts.columns and "dar_gene" not in contacts.columns:
    contacts = contacts.rename(columns={"dar_locus": "dar_gene"})
```

---

### Issue 4.3 — Float-formatted integer strings causing join mismatches

**Files:** `dar_contacts_cbcb8A.csv` and `compensatory_partners.csv`

**Symptom:** Even after fixing the column name, some pairs failed to join. Inspecting the values:

```
compensatory_partners.csv:   dar_aa_coord = "45"     (integer string)
dar_contacts_cbcb8A.csv:     dar_aa_coord = "45.0"   (float string from CSV parsing)
```

**Root cause:** When pandas reads a column that contains some NaN values alongside integers, it silently upcast to `float64`, storing `45` as `45.0`. This is standard pandas behaviour but produces string mismatches when both columns are cast to `str` with `.astype(str)`.

**Fix:** Normalise both columns before the join using:
```python
for df in (partners, contacts_sub):
    for col in ("dar_aa_coord", "contact_refseq_pos"):
        if col in df.columns:
            df[col] = df[col].apply(
                lambda v: str(int(float(v))) if pd.notna(v) else str(v)
            )
```
This converts `"45.0"` → `"45"` on both sides consistently.

---

### Issue 4.4 — Compiled final report had 1,017 rows instead of 693

**File:** `src/mutagenesis/03_compile_targets.py`

**Symptom:** The final `final_targets.csv` had 1,017 rows after merging `prioritized_pairs.csv` (693 rows) with `mi_scores.csv` (693 rows). Expected 693.

**Root cause:** Multiple distinct clinical variant records (`ann_id`) can produce the same physical residue change. For example:
- ClinVar variant 123456: c.576A>T → p.Ile192Thr in MT-ATP6
- ClinVar variant 789012: c.576A>C → p.Ile192Thr in MT-ATP6

Both map to the same amino acid substitution (Ile→Thr at position 192), so both appear as the same physical residue pair in the compensatory analysis. The prioritization score is identical for both, but both appear in `prioritized_pairs.csv` because they have different `ann_id` values. The merge multiplied them.

**Fix:** Deduplicate on the physical residue pair key after merging, and count the number of clinical variants per physical pair:
```python
phys_key = ["dar_gene", "dar_aa_coord", "dar_alt_aa",
            "contact_gene", "contact_refseq_pos", "contact_alt_aa"]

n_var = final.groupby(phys_key).size().rename("n_clinical_variants").reset_index()
final = final.drop_duplicates(subset=phys_key, keep="first")
final = final.merge(n_var, on=phys_key, how="left")
```

The `n_clinical_variants` column is included in the final report as a secondary signal: a physical residue pair supported by multiple independent clinical variant records has stronger disease evidence.

**After fix:** 693 unique physical pairs; the pair MT-ATP6:192T ↔ MT-ATP6:188T had `n_clinical_variants = 3`.

---

### Issue 4.5 — Wrong FoldX command for intraprotein pairs (`AnalyseComplex` vs `BuildModel`)

**Context:** The initial plan specified FoldX `AnalyseComplex` for all compensatory pairs to measure binding affinity changes.

**Problem identified:** `AnalyseComplex` computes the interaction energy between two separate polypeptide chains. It is the correct metric when a disease mutation in chain A disrupts binding to chain B. However, 594/693 pairs (85.7%) are intraprotein — both the disease residue and the compensatory contact residue are within the same polypeptide chain. There is no binding interface to measure for these.

**Consequence of using AnalyseComplex on intraprotein pairs:**
- For intraprotein pairs, `AnalyseComplex` on the same chain gives undefined or meaningless results (it measures the chain's interaction with itself)
- The relevant metric is `ΔΔG_fold` — the change in the protein's folding stability — captured by `BuildModel`

**Contact breakdown used to inform the decision:**
| Category | N | % |
|----------|---|---|
| Intraprotein | 594 | 85.7% |
| Interprotein mt-nuc | 38 | 5.5% |
| Interprotein nuc-nuc | 40 | 5.8% |
| Interprotein mt-mt | 21 | 3.0% |

**Fix:** Split the FoldX protocol by contact category:
- **Intraprotein:** `BuildModel` only → `ΔΔG_stab_DAR`, `ΔΔG_rescue_stab`, `ΔΔG_epistasis`
- **Interprotein:** `BuildModel` + `AnalyseComplex` → stability AND `ΔΔG_bind_DAR`, `ΔΔG_rescue_bind`

The `AnalyseComplex` call only triggers when `dar_chain != contact_chain`, ensuring it is never applied to intraprotein pairs that happen to share a chain letter.

---

### Issue 4.6 — EvoEF2 considered but rejected as open-source ΔΔG tool

**Context:** EvoEF2 was identified as a fully open-source, Docker-friendly alternative to FoldX (which requires academic registration).

**Problem:** EvoEF2 was designed for protein design (generating sequences with desired properties), not for predicting the effect of specific point mutations on stability. Its parameters are optimised for the design task, not for reproducing experimental ΔΔG measurements.

**Benchmarking evidence (2023–2024, Briefings in Bioinformatics):**
- FoldX: R = 0.69–0.71 with experimental ΔΔG, RMSE = 1.25 kcal/mol
- mCSM-PPI2: R = 0.82, RMSE = 1.18 kcal/mol (but web-server only)
- EvoEF2: not formally benchmarked for mutation effect prediction (ΔΔG); no validation data for interface mutations available

**Decision:** EvoEF2 rejected as primary tool. FoldX (`BuildModel` for stability) retained as primary. mCSM-PPI2 documented for manual validation of the top 10 interprotein pairs (highest accuracy, R = 0.82).

---

### Issue 4.7 — MI matrix full computation would be O(L²) per gene

**File:** `src/mutagenesis/02_mi_analysis.py`

**Context:** APC (Average Product Correction) for mutual information requires the row means and global mean of the full MI matrix for each gene — i.e., MI computed for all pairs of alignment columns within the gene, not just the specific (DAR, contact) pair of interest.

**Problem:** For a protein with L = 500 residues, the full MI matrix has L(L-1)/2 = 124,750 pairs. Computing all of them for every gene before extracting the pair of interest would be prohibitively slow (several hours for 106 OxPhos genes).

**Fix:** Sample up to 5,000 random column pairs per gene as background for APC correction:
```python
def compute_gene_mi_stats(gene, genome, sample_pairs=5000):
    if L <= 100:
        pairs = [(i, j) for i in range(L) for j in range(i+1, L)]  # full matrix
    else:
        rng = np.random.default_rng(42)
        idx = rng.integers(0, L, size=(sample_pairs, 2))
        idx = idx[idx[:, 0] != idx[:, 1]]
        pairs = [(int(a), int(b)) for a, b in idx]
```

For small proteins (L ≤ 100), the full matrix is computed exactly. For larger proteins, the sample of 5,000 pairs gives a reliable estimate of row means and global mean (typical OxPhos subunits have L = 100–700 residues, so the sample covers 1–50% of all pairs).

Results are cached per gene so that multiple pairs from the same gene share the background computation.

---

### Issue 4.8 — Inter-protein MI requires concatenated MSA; contact genome inference

**File:** `src/mutagenesis/02_mi_analysis.py`

**Context:** For the 99 interprotein pairs, the DAR and its contact reside in different proteins with different gene MSAs. MI between positions in two different MSAs cannot be computed from separate alignments — you need aligned sequences from the same set of species in the same order.

**Problem:** Simply loading two MSA files and pairing up columns would give incorrect results if the species order differs between files or if some species are present in one alignment but not the other.

**Fix:** Explicitly compute the intersection of species present in both alignments, sort them consistently, and concatenate the sequences:
```python
common_spp = sorted(set(dar_msa) & set(contact_msa))
concat_msa = {s: dar_msa[s] + contact_msa[s] for s in common_spp}
col_dar_in_concat     = dar_aa_coord - 1
col_contact_in_concat = len(dar_msa[common_spp[0]]) + contact_refseq_pos - 1
```

The column index for the contact residue in the concatenated sequence is shifted by the full length of the DAR protein sequence.

**Additional problem — contact genome inference:** For `contact_type = "mt-nuc"`, the DAR is in an mtDNA gene but the contact is in a nucDNA gene. The genome directory (`mtdna_aa` vs `toga_hg38_aa`) must be inferred from the contact_type field, not just from `dar_genome`:
```python
if contact_type == "mt-nuc":
    contact_genome = "nucDNA"   # DAR is mt, contact is nuc
elif contact_type == "mt-mt":
    contact_genome = "mtDNA"
elif contact_type == "nuc-nuc":
    contact_genome = "nucDNA"
```

---

### Issue 4.9 — 78 pairs (11.3%) had no MI value computed

**File:** `src/mutagenesis/02_mi_analysis.py`

**Symptom:** After running MI analysis, 615/693 pairs had `mi_apc` values; 78 pairs returned `NaN`.

**Root causes identified:**
1. **Alignment file not found:** Some gene names in `compensatory_partners.csv` (e.g. `COX4I2`) have isoform-proxy mappings but their primary FASTA file does not exist in `toga_hg38_aa/` — only the proxy (`COX4I1`) file exists. The alignment lookup fails silently with `return None`.

2. **Position out of range:** Some `contact_refseq_pos` values exceed the alignment length for the contact gene (possible for isoform proxies with different sequence lengths). The code returns `NaN` safely without crashing.

3. **Insufficient common species for inter-protein pairs:** A few interprotein pairs have fewer than 10 species in the intersection of their two alignments (some TOGA alignments cover only the well-assembled mammals). Threshold `< 10 species → NaN` to avoid unreliable MI estimates.

**Status:** Not fixed — acceptable at 88.7% coverage. The failing cases are predominantly isoform proxies where alignment ambiguity exists anyway. FoldX results will cover these independently.

---

## 5. Summary Table

| # | Issue | File | Symptom | Root cause | Fix |
|---|-------|------|---------|------------|-----|
| 1.1 | ROOT undefined | prepare_pagel_hpc.py | NameError | Missing ROOT assignment | Add ROOT at top of paths |
| 1.2 | All NA Pagel (HPC run 1) | pagel_discrete.R | 0 finite p-values | force.ultrametric + 5s timeout | force.ultrametric in get_tree(); timeout → 60s |
| 1.3 | All NA Pagel (HPC run 2) | pagel_discrete.R | 0 finite p-values | 60s < HPC pair time (~64s) | timeout → 300s |
| 1.4 | chunk_0001 hangs 3+ hours | pagel_discrete.R | Infinite loop in fitPagel | No timeout; non-convergence | Keep 300s timeout + checkpointing |
| 1.5 | SyntaxError line 615 | 02_phylogenetic_timing.py | SyntaxError on valid code | Docker cached stale file | Fresh docker run |
| 1.6 | force.ultrametric flooding logs | pagel_discrete.R | Hundreds of note lines | message() not suppressed | suppressMessages() wrapper |
| 2.1 | 458/6738 age estimates | 02_phylogenetic_timing.py | dar_origin_age_mya = None for leaves | Single-species sets returned None | Return 0.0 for leaf-origin events |
| 2.2 | contact_first artifact concern | 02_phylogenetic_timing.py | Believed to be wrong | Subset math misread | Algebra confirmed correct; not an artifact |
| 2.3 | physicochemical_type = unknown | compensatory_partners.csv | All "unknown" after re-run | Stale join not re-run | Re-run join step after timing update |
| 3.1 | Fisher test pseudoreplication | 01_find_compensating_partners.py | Inflated p-values | Species not independent | Pagel + branch tests as primary; Fisher as fallback |
| 3.2 | Isoform proxy position offset | 00_map_davs_to_structure.py | Residue mapped to wrong position | Isoform sequence register shift | ±10 AA anchor window; status field |
| 4.1 | 0 PDB IDs joined | 00_prioritize_pairs.py | All cbcb_dist = NaN | ann_id format mismatch between files | Join on physical positions, not ann_id |
| 4.2 | dar_locus vs dar_gene | 00_prioritize_pairs.py | KeyError on join | Different column names same data | Rename dar_locus → dar_gene |
| 4.3 | "45.0" != "45" string mismatch | 00_prioritize_pairs.py | Join misses some rows | Float-format integers from pandas NaN coercion | Normalise via int(float(v)) before str |
| 4.4 | 1017 rows instead of 693 | 03_compile_targets.py | Duplicated physical pairs | Multiple ann_ids per physical AA change | Dedup on physical key; count n_clinical_variants |
| 4.5 | AnalyseComplex wrong for intraprotein | 01_foldx_ddg.py | Protocol mismatch | 85.7% are intraprotein; no binding interface | BuildModel for all; AnalyseComplex only for interprotein |
| 4.6 | EvoEF2 rejected | design | Would give unreliable ΔΔG | Not validated for mutation effect prediction | FoldX as primary; mCSM-PPI2 for validation |
| 4.7 | O(L²) MI matrix | 02_mi_analysis.py | Would take hours | Full MI matrix per gene | Sample 5000 pairs; full matrix for L ≤ 100 |
| 4.8 | MI cross-gene species order | 02_mi_analysis.py | Wrong MI values | Separate MSAs have different species | Intersect species, sort, concatenate |
| 4.9 | 78 pairs no MI value | 02_mi_analysis.py | NaN in mi_scores.csv | Missing alignment, out-of-range pos, <10 spp | Acceptable; FoldX covers independently |
