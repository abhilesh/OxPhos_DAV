# Phylogenetic Analysis Overhaul (phylo_v2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace marginal-ASR point-estimate phylogenetics with stochastic character mapping (SIMMAP) and joint corHMM models to answer three questions: independent origin counts, mtDNA/nucDNA age comparison, and contact-vs-cDAV temporal ordering.

**Architecture:** Pure file-based Python↔R exchange (subprocess + TSV/CSV); no rpy2 dependency at runtime. Each cDAV/pair gets its own character-matrix TSV; R scripts accept file paths as CLI args; Python orchestrators parallelize via `concurrent.futures`. Old `src/phylo/` is untouched.

**Tech Stack:** Python 3.11, R 4.x (phytools, corHMM, ape, geiger, data.table), BioPython, pandas, scipy, statsmodels; Docker `oxphos_dav_analysis`

---

## Repo path corrections vs. spec

The spec names several paths that do not exist in the repo. Use these instead:

| Spec path | Actual path |
|---|---|
| `data/phylo/MamPhy_calibrated.nwk` | `data/phylo/species_tree/mammaltree_crossgenome.nwk` (289 tips, ~87 Mya root-to-tip) |
| `data/curated/cdav_species_table.csv` | Does not exist — build from `lineages_with_disease_allele` in `data/annotations/curated/cdav_classifications_{mt,nuc}DNA.json` |
| `data/curated/contact_species_table.csv` | Does not exist — build from MSAs in `data/phylo/iqtree_jobs/{GENE}/{GENE}.fasta` |
| IQ-TREE rate info | `data/phylo/ancestral_states/{GENE}/{GENE}.iqtree` (tree length, alpha, best-fit model) |

Docker already has `r-base`, `r-ape`, `r-phytools`. Still needs: `r-corhmm`, `r-geiger`, `r-data.table`.

---

## File map

| File | Responsibility |
|---|---|
| `Dockerfile` | Add missing R packages |
| `src/phylo_v2/utils/tree_io.py` | Load/prune calibrated tree, map FASTA species → tree tips |
| `src/phylo_v2/utils/msa_utils.py` | Extract per-species AA at human reference position from MSA |
| `src/phylo_v2/01_prepare_character_matrices.py` | Build per-cDAV TSVs + manifest |
| `src/phylo_v2/02_simmap_origins.R` | SIMMAP on one cDAV character vector; outputs posterior samples |
| `src/phylo_v2/02_simmap_origins_runner.py` | Parallel subprocess caller for 02 |
| `src/phylo_v2/03_per_gene_rate_calibration.py` | Parse IQ-TREE logs → per-gene subst rate |
| `src/phylo_v2/04_age_comparison.py` | Mixed-effects regression mt vs nuc ages |
| `src/phylo_v2/05_resolvability_filter.py` | Flag pairs testable for ordering |
| `src/phylo_v2/06_corhmm_joint_models.R` | Dependent vs independent corHMM on one pair |
| `src/phylo_v2/06_corhmm_runner.py` | Parallel subprocess caller for 06 |
| `src/phylo_v2/07_independence_null.py` | Tree-aware null for contact_first proportion |
| `src/phylo_v2/08_compile_phylo_results.py` | Merge all outputs into annotated pair table |
| `tests/test_phylo_v2_tree_io.py` | Tests for tree_io + msa_utils |
| `tests/test_phylo_v2_charmat.py` | Tests for character matrix building |
| `tests/test_phylo_v2_rates.py` | Tests for rate calibration parsing |
| `tests/test_phylo_v2_resolvability.py` | Tests for resolvability filter |
| `tests/test_phylo_v2_null.py` | Tests for independence null |

---

## Task 1: Dockerfile — add missing R packages

**Files:**
- Modify: `Dockerfile`

- [ ] **Step 1: Write the smoke test**

Create `tests/test_phylo_v2_r_packages.py`:

```python
"""Smoke-test that required R packages load inside the Docker env."""
import subprocess

def test_r_packages_available():
    r_code = "suppressMessages(library(corHMM)); suppressMessages(library(geiger)); suppressMessages(library(data.table)); cat('OK')"
    result = subprocess.run(
        ["conda", "run", "-n", "oxphos_dav", "Rscript", "--vanilla", "-e", r_code],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, f"R packages missing:\n{result.stderr}"
    assert "OK" in result.stdout
```

- [ ] **Step 2: Add R packages to Dockerfile**

In `Dockerfile`, after the existing `r-phytools` conda install line, add a new RUN step:

```dockerfile
# 5. Additional R packages for phylo_v2 pipeline
RUN conda run -n oxphos_dav R -e \
    "install.packages(c('corHMM','geiger','data.table'), repos='https://cloud.r-project.org', Ncpus=4)"
```

- [ ] **Step 3: Rebuild and run test**

```bash
docker build -t oxphos_dav_analysis .
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav pytest tests/test_phylo_v2_r_packages.py -v
```
Expected: PASS.

---

## Task 2: `utils/tree_io.py` — tree loading and pruning

**Files:**
- Create: `src/phylo_v2/__init__.py` (empty)
- Create: `src/phylo_v2/utils/__init__.py` (empty)
- Create: `src/phylo_v2/utils/tree_io.py`
- Test: `tests/test_phylo_v2_tree_io.py`

- [ ] **Step 1: Write tests**

```python
# tests/test_phylo_v2_tree_io.py
from pathlib import Path
import pytest
from src.phylo_v2.utils.tree_io import (
    load_calibrated_tree,
    get_tree_tip_names,
    prune_tree_to_species,
    compute_root_to_tip_dist,
)

TREE_PATH = Path("data/phylo/species_tree/mammaltree_crossgenome.nwk")

@pytest.fixture(scope="module")
def full_tree():
    return load_calibrated_tree(TREE_PATH)

def test_load_tree_tip_count(full_tree):
    tips = get_tree_tip_names(full_tree)
    assert len(tips) == 289

def test_tip_names_are_binomial(full_tree):
    tips = get_tree_tip_names(full_tree)
    for t in tips:
        assert "_" in t, f"Tip {t!r} not in Genus_species format"

def test_prune_reduces_tips(full_tree):
    keep = ["Mus_musculus", "Homo_sapiens", "Pan_troglodytes"]
    pruned = prune_tree_to_species(full_tree, keep)
    tips = get_tree_tip_names(pruned)
    assert set(tips) == set(keep)

def test_prune_unknown_species_raises(full_tree):
    with pytest.raises(ValueError, match="not in tree"):
        prune_tree_to_species(full_tree, ["Mus_musculus", "Draco_fictus"])

def test_root_to_tip_dist_positive(full_tree):
    d = compute_root_to_tip_dist(full_tree, "Mus_musculus")
    assert d > 0

def test_prune_to_one_species_raises(full_tree):
    with pytest.raises(ValueError, match="at least 4"):
        prune_tree_to_species(full_tree, ["Mus_musculus"])
```

- [ ] **Step 2: Implement `tree_io.py`**

```python
# src/phylo_v2/utils/tree_io.py
"""
Calibrated tree utilities for phylo_v2.
All functions use Bio.Phylo with the Newick format.
Branch lengths in mammaltree_crossgenome.nwk are in Mya.
"""
from __future__ import annotations
from pathlib import Path
from copy import deepcopy
from Bio import Phylo

CALIBRATED_TREE = Path("data/phylo/species_tree/mammaltree_crossgenome.nwk")
MIN_TIPS_FOR_ANALYSIS = 4


def load_calibrated_tree(path: Path | str = CALIBRATED_TREE):
    """Load the calibrated mammalian timetree (Mya branch lengths)."""
    return Phylo.read(str(path), "newick")


def get_tree_tip_names(tree) -> list[str]:
    """Return sorted list of tip (terminal) names."""
    return sorted(c.name for c in tree.get_terminals())


def prune_tree_to_species(tree, species: list[str]):
    """
    Return a copy of *tree* pruned to *species*.

    Raises ValueError if any species is absent from the tree, or if
    fewer than MIN_TIPS_FOR_ANALYSIS tips would remain.
    """
    tip_set = set(get_tree_tip_names(tree))
    missing = set(species) - tip_set
    if missing:
        raise ValueError(f"Species not in tree: {sorted(missing)}")
    if len(species) < MIN_TIPS_FOR_ANALYSIS:
        raise ValueError(
            f"Need at least {MIN_TIPS_FOR_ANALYSIS} species for phylogenetic analysis, got {len(species)}"
        )
    t = deepcopy(tree)
    to_prune = tip_set - set(species)
    for sp in to_prune:
        t.prune(sp)
    return t


def compute_root_to_tip_dist(tree, tip_name: str) -> float:
    """Distance (Mya) from root to the named tip."""
    return tree.distance(tip_name)


def write_newick(tree, path: Path | str) -> None:
    """Write *tree* to Newick file."""
    Phylo.write(tree, str(path), "newick")
```

- [ ] **Step 3: Run tests**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav pytest tests/test_phylo_v2_tree_io.py -v
```
Expected: all PASS.

---

## Task 3: `utils/msa_utils.py` — extract per-species AA at reference position

**Files:**
- Create: `src/phylo_v2/utils/msa_utils.py`
- Test: append to `tests/test_phylo_v2_tree_io.py`

The FASTA files in `data/phylo/iqtree_jobs/{GENE}/{GENE}.fasta` are protein MSAs. The human sequence ID is `Homo_sapiens`. We map a 1-based ungapped human reference position to an alignment column, then extract each species' AA.

- [ ] **Step 1: Append tests to `tests/test_phylo_v2_tree_io.py`**

```python
# append to tests/test_phylo_v2_tree_io.py
from src.phylo_v2.utils.msa_utils import (
    load_msa,
    human_refpos_to_aln_col,
    extract_states_at_position,
)

MT_CO1_FASTA = Path("data/phylo/iqtree_jobs/MT-CO1/MT-CO1.fasta")

def test_load_msa_returns_dict():
    msa = load_msa(MT_CO1_FASTA)
    assert "Homo_sapiens" in msa
    assert len(msa) >= 100

def test_refpos_to_aln_col_positive():
    msa = load_msa(MT_CO1_FASTA)
    col = human_refpos_to_aln_col(msa, 100)
    assert col >= 99  # alignment col >= ref pos - 1 (0-based)

def test_extract_states_returns_all_species():
    msa = load_msa(MT_CO1_FASTA)
    states = extract_states_at_position(msa, refpos=100)
    assert "Homo_sapiens" in states
    assert all(len(v) == 1 or v == "-" for v in states.values())

def test_extract_states_human_is_not_gap():
    msa = load_msa(MT_CO1_FASTA)
    states = extract_states_at_position(msa, refpos=100)
    assert states["Homo_sapiens"] != "-"
```

- [ ] **Step 2: Implement `msa_utils.py`**

```python
# src/phylo_v2/utils/msa_utils.py
"""
MSA utilities for phylo_v2.
Loads protein MSAs (gap char = '-') and maps human reference positions
to alignment columns so we can extract per-species amino acids.
"""
from __future__ import annotations
from pathlib import Path
from Bio import SeqIO

HUMAN_ID = "Homo_sapiens"


def load_msa(fasta_path: Path | str) -> dict[str, str]:
    """Return {species_id: aligned_sequence} from a protein FASTA MSA."""
    return {r.id: str(r.seq) for r in SeqIO.parse(str(fasta_path), "fasta")}


def human_refpos_to_aln_col(msa: dict[str, str], refpos: int) -> int:
    """
    Convert a 1-based ungapped human reference position to a 0-based
    alignment column index.

    Raises KeyError if HUMAN_ID is not in msa.
    Raises IndexError if refpos exceeds the ungapped human sequence length.
    """
    human_seq = msa[HUMAN_ID]
    ungapped_count = 0
    for col, aa in enumerate(human_seq):
        if aa != "-":
            ungapped_count += 1
            if ungapped_count == refpos:
                return col
    raise IndexError(
        f"Reference position {refpos} exceeds ungapped human sequence length {ungapped_count}"
    )


def extract_states_at_position(
    msa: dict[str, str],
    refpos: int,
) -> dict[str, str]:
    """
    Return {species: aa_char} for every species in the MSA at the
    alignment column corresponding to human reference position *refpos*.

    Gap characters ('-') are returned as-is; callers decide whether to
    treat them as missing data.
    """
    col = human_refpos_to_aln_col(msa, refpos)
    return {sp: seq[col] for sp, seq in msa.items()}
```

- [ ] **Step 3: Run tests**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav pytest tests/test_phylo_v2_tree_io.py -v
```
Expected: all PASS.

---

## Task 4: `01_prepare_character_matrices.py`

**Files:**
- Create: `src/phylo_v2/01_prepare_character_matrices.py`
- Test: `tests/test_phylo_v2_charmat.py`

**Logic:**
1. Load both JSON classification files; collect every record where `is_cdav_amino_acid == True`. The cDAV id is `ann_id`. The species carrying the cDAV is `lineages_with_disease_allele`.
2. For each cDAV, find its gene (`locus` field) and load the MSA at `data/phylo/iqtree_jobs/{locus}/{locus}.fasta`. The full species set for the cDAV column is every sequence ID in that FASTA.
3. Build `cdav_state`: 1 for species in `lineages_with_disease_allele`, 0 otherwise, NA for species absent from the MSA.
4. Load `results/structural/all_tested_pairs.csv`. For each contact tested against this cDAV, extract contact states from the contact gene's MSA using `extract_states_at_position(msa, contact_refseq_pos)`. Contact state = 1 if species AA == `contact_alt_aa`, 0 if it equals `contact_human_aa`, NA for gaps or other AAs.
5. Restrict rows to species present in the calibrated tree. Log dropped counts.
6. Write TSV to `results/phylo_v2/character_matrices/{cdav_id}.tsv`.
7. Write manifest to `results/phylo_v2/character_matrices/_manifest.csv`.

- [ ] **Step 1: Write tests**

```python
# tests/test_phylo_v2_charmat.py
import csv
from pathlib import Path

MATRIX_DIR = Path("results/phylo_v2/character_matrices")
MANIFEST = MATRIX_DIR / "_manifest.csv"


def _find_any_manifest_row(gene_filter=None):
    if not MANIFEST.exists():
        return None
    with open(MANIFEST) as f:
        rows = list(csv.DictReader(f))
    if gene_filter:
        rows = [r for r in rows if r["cdav_gene"] == gene_filter]
    return rows[0] if rows else None


def test_manifest_exists():
    assert MANIFEST.exists(), "Run 01_prepare_character_matrices.py first"


def test_manifest_has_required_columns():
    with open(MANIFEST) as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames
    required = [
        "cdav_id", "cdav_gene", "cdav_genome", "cdav_position",
        "n_species_total", "n_species_with_cdav", "n_contacts_tested",
        "informative_for_simmap", "informative_for_corhmm",
    ]
    for c in required:
        assert c in cols, f"Missing column: {c}"


def test_matrix_file_exists_for_each_manifest_row():
    with open(MANIFEST) as f:
        rows = list(csv.DictReader(f))
    for row in rows[:20]:  # spot-check first 20
        matrix_file = MATRIX_DIR / f"{row['cdav_id']}.tsv"
        assert matrix_file.exists(), f"Missing matrix: {matrix_file}"


def test_matrix_has_species_and_cdav_state_columns():
    row = _find_any_manifest_row()
    assert row is not None
    matrix_file = MATRIX_DIR / f"{row['cdav_id']}.tsv"
    with open(matrix_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        cols = reader.fieldnames
    assert "species" in cols
    assert "cdav_state" in cols


def test_cdav_state_values_are_0_1_or_na():
    row = _find_any_manifest_row()
    assert row is not None
    matrix_file = MATRIX_DIR / f"{row['cdav_id']}.tsv"
    with open(matrix_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for data_row in reader:
            assert data_row["cdav_state"] in {"0", "1", "NA"}, (
                f"Unexpected cdav_state: {data_row['cdav_state']}"
            )


def test_minimum_20_species_flag_correct():
    with open(MANIFEST) as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        n = int(row["n_species_total"])
        expected_flag = str(n >= 20).upper()  # "TRUE" or "FALSE"
        assert row["informative_for_simmap"].upper() == expected_flag, (
            f"cdav_id={row['cdav_id']}: n={n} but informative_for_simmap={row['informative_for_simmap']}"
        )


def test_no_species_duplicates_in_matrix():
    row = _find_any_manifest_row()
    assert row is not None
    matrix_file = MATRIX_DIR / f"{row['cdav_id']}.tsv"
    with open(matrix_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        species_list = [r["species"] for r in reader]
    assert len(species_list) == len(set(species_list)), "Duplicate species in matrix"
```

- [ ] **Step 2: Implement the script**

```python
# src/phylo_v2/01_prepare_character_matrices.py
"""
Build per-cDAV character matrix TSVs for phylo_v2 SIMMAP / corHMM analysis.

Usage (inside Docker):
    python src/phylo_v2/01_prepare_character_matrices.py

Outputs:
    results/phylo_v2/character_matrices/{cdav_id}.tsv
    results/phylo_v2/character_matrices/_manifest.csv
"""
from __future__ import annotations
import csv
import json
import logging
from collections import defaultdict
from pathlib import Path

from src.phylo_v2.utils.tree_io import load_calibrated_tree, get_tree_tip_names
from src.phylo_v2.utils.msa_utils import load_msa, extract_states_at_position

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
ANN_DIR = ROOT / "data" / "annotations" / "curated"
IQTREE_DIR = ROOT / "data" / "phylo" / "iqtree_jobs"
TREE_PATH = ROOT / "data" / "phylo" / "species_tree" / "mammaltree_crossgenome.nwk"
PAIRS_CSV = ROOT / "results" / "structural" / "all_tested_pairs.csv"
OUT_DIR = ROOT / "results" / "phylo_v2" / "character_matrices"
MIN_SPECIES_SIMMAP = 20


def load_cdav_records() -> list[dict]:
    """Load all AA-level cDAV records from both JSON classification files."""
    records = []
    for fname in ("cdav_classifications_mtDNA.json", "cdav_classifications_nucDNA.json"):
        path = ANN_DIR / fname
        with open(path) as f:
            data = json.load(f)
        for rec in data:
            if rec.get("is_cdav_amino_acid"):
                records.append(rec)
    return records


def load_pairs_by_cdav() -> dict[str, list[dict]]:
    """Return {ann_id: [pair_row, ...]} from all_tested_pairs.csv."""
    by_cdav: dict[str, list[dict]] = defaultdict(list)
    with open(PAIRS_CSV) as f:
        for row in csv.DictReader(f):
            by_cdav[row["ann_id"]].append(row)
    return by_cdav


def build_matrix(
    cdav_rec: dict,
    pairs: list[dict],
    tree_tips: set[str],
) -> tuple[list[dict], dict]:
    """
    Build the character matrix rows for one cDAV.

    Returns (rows, manifest_meta).
    """
    gene = cdav_rec["locus"]
    fasta_path = IQTREE_DIR / gene / f"{gene}.fasta"
    if not fasta_path.exists():
        log.warning("FASTA not found for gene %s, skipping cDAV %s", gene, cdav_rec["ann_id"])
        return [], {}

    msa = load_msa(fasta_path)
    msa_species = set(msa.keys())

    cdav_species_set = set(cdav_rec.get("lineages_with_disease_allele") or [])

    # Base species universe: intersection of MSA and tree tips
    universe = msa_species & tree_tips
    dropped = msa_species - tree_tips
    if dropped:
        log.debug("cDAV %s: dropped %d species absent from tree", cdav_rec["ann_id"], len(dropped))

    cdav_state: dict[str, str] = {
        sp: ("1" if sp in cdav_species_set else "0") for sp in universe
    }

    # --- per-contact states ---
    contact_msas: dict[str, dict] = {}
    contact_col_list: list[str] = []

    for pair in pairs:
        c_gene = pair["contact_gene"]
        c_pos = int(pair["contact_refseq_pos"])
        c_alt = pair["contact_alt_aa"]
        c_ref = pair["contact_human_aa"]
        col_name = f"contact_{c_gene}_{c_pos}_{c_alt}"
        if col_name not in contact_col_list:
            contact_col_list.append(col_name)

        if c_gene not in contact_msas:
            c_fasta = IQTREE_DIR / c_gene / f"{c_gene}.fasta"
            contact_msas[c_gene] = load_msa(c_fasta) if c_fasta.exists() else {}

    # Build row list
    rows = []
    for sp in sorted(universe):
        row: dict[str, str] = {"species": sp, "cdav_state": cdav_state[sp]}
        for pair in pairs:
            c_gene = pair["contact_gene"]
            c_pos = int(pair["contact_refseq_pos"])
            c_alt = pair["contact_alt_aa"]
            c_ref = pair["contact_human_aa"]
            col_name = f"contact_{c_gene}_{c_pos}_{c_alt}"
            c_msa = contact_msas.get(c_gene, {})
            if not c_msa:
                row[col_name] = "NA"
                continue
            try:
                aa_map = extract_states_at_position(c_msa, c_pos)
                aa = aa_map.get(sp, "-")
                if aa == "-" or aa not in (c_alt, c_ref):
                    row[col_name] = "NA"
                else:
                    row[col_name] = "1" if aa == c_alt else "0"
            except (IndexError, KeyError):
                row[col_name] = "NA"
        rows.append(row)

    n_total = len(rows)
    n_with_cdav = sum(1 for r in rows if r["cdav_state"] == "1")

    meta = {
        "cdav_id": cdav_rec["ann_id"],
        "cdav_gene": gene,
        "cdav_genome": cdav_rec.get("genome", ""),
        "cdav_position": cdav_rec.get("aa_change", ""),
        "n_species_total": n_total,
        "n_species_with_cdav": n_with_cdav,
        "n_contacts_tested": len(pairs),
        "min_contact_n_with_alt": min(
            (sum(1 for r in rows if r.get(c) == "1") for c in contact_col_list),
            default=0,
        ),
        "informative_for_simmap": str(n_total >= MIN_SPECIES_SIMMAP).upper(),
        "informative_for_corhmm": str(
            n_total >= MIN_SPECIES_SIMMAP and len(pairs) > 0
        ).upper(),
    }
    return rows, meta


def write_matrix(cdav_id: str, rows: list[dict]) -> None:
    if not rows:
        return
    cols = list(rows[0].keys())
    path = OUT_DIR / f"{cdav_id}.tsv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log.info("Loading calibrated tree...")
    tree = load_calibrated_tree(TREE_PATH)
    tree_tips = set(get_tree_tip_names(tree))

    log.info("Loading cDAV records...")
    cdav_records = load_cdav_records()
    log.info("Found %d AA-level cDAVs", len(cdav_records))

    log.info("Loading tested pairs...")
    pairs_by_cdav = load_pairs_by_cdav()

    manifest_rows = []
    skipped = 0
    for rec in cdav_records:
        ann_id = rec["ann_id"]
        pairs = pairs_by_cdav.get(ann_id, [])
        rows, meta = build_matrix(rec, pairs, tree_tips)
        if not rows:
            skipped += 1
            continue
        write_matrix(ann_id, rows)
        manifest_rows.append(meta)

    manifest_path = OUT_DIR / "_manifest.csv"
    if manifest_rows:
        with open(manifest_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()))
            writer.writeheader()
            writer.writerows(manifest_rows)

    log.info("Wrote %d character matrices (%d skipped). Manifest: %s",
             len(manifest_rows), skipped, manifest_path)


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run script then tests**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav python src/phylo_v2/01_prepare_character_matrices.py

docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav pytest tests/test_phylo_v2_charmat.py -v
```
Expected: all PASS. Log reports total matrices written and skipped count.

---

## Task 5: `02_simmap_origins.R` — SIMMAP for one cDAV

**Files:**
- Create: `src/phylo_v2/02_simmap_origins.R`

CLI: `Rscript 02_simmap_origins.R <matrix_tsv> <gene_tree_nwk> <out_dir> <nsim>`

Outputs in `<out_dir>/`:
- `{cdav_id}_simmap_summary.csv` — one row per simulation: n_origins, mean_origin_age_relative
- `{cdav_id}_simmap_posterior.csv` — posterior summaries (median, 2.5%, 97.5%) across sims

- [ ] **Step 1: Write the R script**

```r
#!/usr/bin/env Rscript
# src/phylo_v2/02_simmap_origins.R
#
# Run stochastic character mapping (SIMMAP) for one cDAV binary character.
# Outputs posterior distribution of independent origin counts and ages.
#
# Usage:
#   Rscript 02_simmap_origins.R <matrix_tsv> <gene_tree_nwk> <out_dir> <nsim>

suppressMessages({
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 02_simmap_origins.R <matrix_tsv> <gene_tree_nwk> <out_dir> <nsim>")
}

matrix_tsv  <- args[1]
tree_nwk    <- args[2]
out_dir     <- args[3]
nsim        <- as.integer(args[4])

# ── Load data ─────────────────────────────────────────────────────────────────
mat  <- read.delim(matrix_tsv, stringsAsFactors = FALSE, na.strings = "NA")
tree <- read.tree(tree_nwk)

cdav_id <- sub("\\.tsv$", "", basename(matrix_tsv))

# ── Prune tree to species with non-NA cdav_state ──────────────────────────────
valid  <- mat[!is.na(mat$cdav_state), ]
shared <- intersect(valid$species, tree$tip.label)

if (length(shared) < 10) {
  cat(sprintf("SKIP: %s has only %d species with non-NA cdav_state\n", cdav_id, length(shared)))
  quit(status = 0)
}

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, shared))

x <- setNames(
  as.integer(valid$cdav_state[match(shared, valid$species)]),
  shared
)[tree_pruned$tip.label]

# ── Fit ER model, then run SIMMAP ─────────────────────────────────────────────
fit_er <- fitMk(tree_pruned, x, model = "ER")

set.seed(42)
simmaps <- make.simmap(tree_pruned, x, model = fit_er$rates[1, 2],
                        Q = "mcmc", nsim = nsim, pi = "fitzjohn",
                        message = FALSE, verbose = FALSE)

# ── Count 0→1 origins per simulation ─────────────────────────────────────────
count_origins <- function(simmap_tree) {
  total <- 0L
  for (m in simmap_tree$maps) {
    states <- names(m)
    for (i in seq_along(states)[-1]) {
      if (states[i - 1] == "0" && states[i] == "1") total <- total + 1L
    }
  }
  total
}

# ── Mean age of origin branches (relative to gene tree root) ──────────────────
root_depth <- max(node.depth.edgelength(tree_pruned))

origin_age_relative <- function(simmap_tree) {
  ages      <- numeric(0)
  nd        <- node.depth.edgelength(simmap_tree)
  edge_mat  <- simmap_tree$edge
  for (e in seq_len(nrow(edge_mat))) {
    m      <- simmap_tree$maps[[e]]
    states <- names(m)
    child_depth <- nd[edge_mat[e, 2]]
    cumlen <- cumsum(m)
    for (i in seq_along(states)[-1]) {
      if (states[i - 1] == "0" && states[i] == "1") {
        pos_from_child <- cumlen[i] - m[i]
        ages <- c(ages, root_depth - (child_depth - pos_from_child))
      }
    }
  }
  if (length(ages) == 0) NA_real_ else mean(ages)
}

sim_results <- lapply(simmaps, function(s) list(
  n_origins                = count_origins(s),
  mean_origin_age_relative = origin_age_relative(s)
))

# ── Write outputs ──────────────────────────────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sim_df <- data.frame(
  cdav_id                  = cdav_id,
  sim_index                = seq_along(sim_results),
  n_origins                = sapply(sim_results, `[[`, "n_origins"),
  mean_origin_age_relative = sapply(sim_results, `[[`, "mean_origin_age_relative")
)
write.csv(sim_df, file.path(out_dir, paste0(cdav_id, "_simmap_summary.csv")), row.names = FALSE)

n_orig_vec <- sim_df$n_origins
age_vec    <- sim_df$mean_origin_age_relative[!is.na(sim_df$mean_origin_age_relative)]

posterior_df <- data.frame(
  cdav_id               = cdav_id,
  n_species_used        = length(shared),
  n_simulations         = nsim,
  n_origins_median      = median(n_orig_vec),
  n_origins_mean        = mean(n_orig_vec),
  n_origins_lower95     = quantile(n_orig_vec, 0.025),
  n_origins_upper95     = quantile(n_orig_vec, 0.975),
  n_origins_prob_multi  = mean(n_orig_vec > 1),
  age_relative_median   = ifelse(length(age_vec) > 0, median(age_vec),  NA),
  age_relative_lower95  = ifelse(length(age_vec) > 0, quantile(age_vec, 0.025), NA),
  age_relative_upper95  = ifelse(length(age_vec) > 0, quantile(age_vec, 0.975), NA)
)
write.csv(posterior_df, file.path(out_dir, paste0(cdav_id, "_simmap_posterior.csv")), row.names = FALSE)

cat(sprintf("Done: %s | n_origins median=%.1f (%.1f–%.1f) | prob_multi=%.3f\n",
            cdav_id,
            posterior_df$n_origins_median,
            posterior_df$n_origins_lower95,
            posterior_df$n_origins_upper95,
            posterior_df$n_origins_prob_multi))
```

- [ ] **Step 2: Smoke-test with a single cDAV**

Pick a cDAV from MT-CO1 with `informative_for_simmap=TRUE` from the manifest:

```bash
CDAV_ID=$(awk -F',' 'NR>1 && $2=="MT-CO1" && $8=="TRUE" {print $1; exit}' \
  results/phylo_v2/character_matrices/_manifest.csv)

docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav Rscript src/phylo_v2/02_simmap_origins.R \
    results/phylo_v2/character_matrices/${CDAV_ID}.tsv \
    data/phylo/iqtree_jobs/MT-CO1/MT-CO1_tree.nwk \
    results/phylo_v2/simmap/ \
    100
```
Expected: prints "Done: ... | n_origins median=..." and writes two CSV files.

---

## Task 6: `02_simmap_origins_runner.py` — parallel Python orchestrator

**Files:**
- Create: `src/phylo_v2/02_simmap_origins_runner.py`
- Test: append to `tests/test_phylo_v2_charmat.py`

Reads `_manifest.csv`, filters `informative_for_simmap=TRUE`, runs `02_simmap_origins.R` in parallel. After each R run, converts relative ages to absolute Mya:
`age_mya = age_relative × (species_tree_root_depth / gene_tree_root_depth)`

- [ ] **Step 1: Append tests**

```python
# append to tests/test_phylo_v2_charmat.py
import subprocess
from pathlib import Path

SIMMAP_DIR = Path("results/phylo_v2/simmap")

def test_simmap_runner_produces_posteriors():
    result = subprocess.run(
        ["conda", "run", "-n", "oxphos_dav", "python",
         "src/phylo_v2/02_simmap_origins_runner.py",
         "--max-cdavs", "5", "--nsim", "50"],
        capture_output=True, text=True, cwd="."
    )
    assert result.returncode == 0, f"Runner failed:\n{result.stderr}"
    posteriors = list(SIMMAP_DIR.glob("*_simmap_posterior.csv"))
    assert len(posteriors) >= 1, "No posterior files written"


def test_simmap_posterior_has_mya_column():
    posteriors = list(SIMMAP_DIR.glob("*_simmap_posterior.csv"))
    assert posteriors, "Run 02_simmap_origins_runner.py first"
    with open(posteriors[0]) as f:
        cols = csv.DictReader(f).fieldnames
    assert "age_mya_median" in cols
```

- [ ] **Step 2: Implement**

```python
# src/phylo_v2/02_simmap_origins_runner.py
"""
Parallel runner for 02_simmap_origins.R.

Usage:
    python src/phylo_v2/02_simmap_origins_runner.py [--max-cdavs N] [--nsim N] [--workers N]
"""
from __future__ import annotations
import argparse
import csv
import logging
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import Phylo

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
MANIFEST   = ROOT / "results" / "phylo_v2" / "character_matrices" / "_manifest.csv"
MATRIX_DIR = ROOT / "results" / "phylo_v2" / "character_matrices"
IQTREE_DIR = ROOT / "data" / "phylo" / "iqtree_jobs"
SPECIES_TREE = ROOT / "data" / "phylo" / "species_tree" / "mammaltree_crossgenome.nwk"
SIMMAP_DIR = ROOT / "results" / "phylo_v2" / "simmap"
R_SCRIPT   = ROOT / "src" / "phylo_v2" / "02_simmap_origins.R"


def _tree_root_depth(path: Path) -> float | None:
    if not path.exists():
        return None
    tree = Phylo.read(str(path), "newick")
    tips = tree.get_terminals()
    return tree.distance(tips[0].name) if tips else None


def run_one_cdav(cdav_id: str, gene: str, nsim: int, spp_root: float) -> dict | None:
    matrix_tsv = MATRIX_DIR / f"{cdav_id}.tsv"
    gene_tree  = IQTREE_DIR / gene / f"{gene}_tree.nwk"
    if not matrix_tsv.exists() or not gene_tree.exists():
        return None

    result = subprocess.run(
        ["conda", "run", "-n", "oxphos_dav",
         "Rscript", "--vanilla", str(R_SCRIPT),
         str(matrix_tsv), str(gene_tree), str(SIMMAP_DIR), str(nsim)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        log.error("R failed for %s:\n%s", cdav_id, result.stderr[-500:])
        return None

    posterior_path = SIMMAP_DIR / f"{cdav_id}_simmap_posterior.csv"
    if not posterior_path.exists():
        return None

    with open(posterior_path) as f:
        row = next(csv.DictReader(f))

    g_root = _tree_root_depth(gene_tree)
    if g_root and g_root > 0:
        scale = spp_root / g_root
        for suffix in ("median", "lower95", "upper95"):
            rel_key = f"age_relative_{suffix}"
            mya_key = f"age_mya_{suffix}"
            try:
                row[mya_key] = str(float(row[rel_key]) * scale)
            except (ValueError, KeyError):
                row[mya_key] = "NA"
    else:
        for suffix in ("median", "lower95", "upper95"):
            row[f"age_mya_{suffix}"] = "NA"

    with open(posterior_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)

    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-cdavs", type=int, default=None)
    parser.add_argument("--nsim", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    SIMMAP_DIR.mkdir(parents=True, exist_ok=True)

    with open(MANIFEST) as f:
        manifest = [r for r in csv.DictReader(f) if r["informative_for_simmap"] == "TRUE"]
    if args.max_cdavs:
        manifest = manifest[: args.max_cdavs]
    log.info("Running SIMMAP for %d cDAVs", len(manifest))

    spp_root = _tree_root_depth(SPECIES_TREE)
    log.info("Species tree root depth: %.2f Mya", spp_root)

    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one_cdav, r["cdav_id"], r["cdav_gene"], args.nsim, spp_root): r["cdav_id"]
                   for r in manifest}
        for fut in as_completed(futures):
            try:
                res = fut.result()
                if res:
                    results.append(res)
            except Exception as exc:
                log.error("Exception for %s: %s", futures[fut], exc)

    log.info("Completed SIMMAP for %d / %d cDAVs", len(results), len(manifest))


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run tests**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav pytest tests/test_phylo_v2_charmat.py -v
```
Expected: all PASS.

---

## Task 7: `03_per_gene_rate_calibration.py`

**Files:**
- Create: `src/phylo_v2/03_per_gene_rate_calibration.py`
- Test: `tests/test_phylo_v2_rates.py`

Parses `data/phylo/ancestral_states/{GENE}/{GENE}.iqtree` for total tree length, best-fit model, gamma alpha. Computes `subst_rate_per_mya = total_tree_length / species_tree_total_length_mya` (where the species tree is pruned to the gene's species before summing branch lengths).

- [ ] **Step 1: Write tests**

```python
# tests/test_phylo_v2_rates.py
import csv
from pathlib import Path

RATES_CSV = Path("results/phylo_v2/per_gene_rates.csv")


def test_rates_csv_exists():
    assert RATES_CSV.exists(), "Run 03_per_gene_rate_calibration.py first"


def test_rates_has_required_columns():
    with open(RATES_CSV) as f:
        cols = csv.DictReader(f).fieldnames
    for c in ["gene", "genome", "total_tree_length", "best_fit_model",
              "gamma_alpha", "subst_rate_per_mya", "n_species"]:
        assert c in cols, f"Missing column: {c}"


def test_rates_all_positive():
    with open(RATES_CSV) as f:
        for row in csv.DictReader(f):
            assert float(row["total_tree_length"]) > 0
            if row["subst_rate_per_mya"] not in ("NA", ""):
                assert float(row["subst_rate_per_mya"]) > 0


def test_mt_genes_present():
    with open(RATES_CSV) as f:
        genes = {row["gene"] for row in csv.DictReader(f)}
    assert "MT-CO1" in genes
    assert "MT-ND1" in genes


def test_nuc_genes_present():
    with open(RATES_CSV) as f:
        genes = {row["gene"] for row in csv.DictReader(f)}
    assert "ATP5F1A" in genes
```

- [ ] **Step 2: Implement**

```python
# src/phylo_v2/03_per_gene_rate_calibration.py
"""
Parse IQ-TREE logs to compute per-gene amino-acid substitution rates.

Output: results/phylo_v2/per_gene_rates.csv

Run:
    python src/phylo_v2/03_per_gene_rate_calibration.py
"""
from __future__ import annotations
import csv
import re
from copy import deepcopy
from pathlib import Path

from Bio import Phylo, SeqIO

ROOT = Path(__file__).resolve().parents[2]
ANC_DIR      = ROOT / "data" / "phylo" / "ancestral_states"
IQTREE_DIR   = ROOT / "data" / "phylo" / "iqtree_jobs"
SPECIES_TREE = ROOT / "data" / "phylo" / "species_tree" / "mammaltree_crossgenome.nwk"
OUT_CSV      = ROOT / "results" / "phylo_v2" / "per_gene_rates.csv"

MT_GENES = {
    "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3",
    "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4",
    "MT-ND4L", "MT-ND5", "MT-ND6",
}

_RE_TREE_LEN   = re.compile(r"Total tree length \(sum of branch lengths\):\s+([\d.]+)")
_RE_BEST_MODEL = re.compile(r"Best-fit model according to BIC:\s+(\S+)")
_RE_ALPHA      = re.compile(r"Gamma shape alpha:\s+([\d.]+)")


def parse_iqtree_log(log_path: Path) -> dict:
    text = log_path.read_text(errors="replace")
    tree_len_m = _RE_TREE_LEN.search(text)
    model_m    = _RE_BEST_MODEL.search(text)
    alpha_m    = _RE_ALPHA.search(text)
    return {
        "total_tree_length": float(tree_len_m.group(1)) if tree_len_m else None,
        "best_fit_model":    model_m.group(1) if model_m else "unknown",
        "gamma_alpha":       float(alpha_m.group(1)) if alpha_m else None,
    }


def fasta_species(gene: str) -> set[str]:
    path = IQTREE_DIR / gene / f"{gene}.fasta"
    if not path.exists():
        return set()
    return {r.id for r in SeqIO.parse(str(path), "fasta")}


def species_tree_length_for_gene(tree, gene_species: set[str]) -> float | None:
    all_tips = {c.name for c in tree.get_terminals()}
    keep = list(gene_species & all_tips)
    if len(keep) < 4:
        return None
    t = deepcopy(tree)
    for tip in all_tips - set(keep):
        t.prune(tip)
    return sum(c.branch_length for c in t.find_clades() if c.branch_length)


def main():
    spp_tree = Phylo.read(str(SPECIES_TREE), "newick")
    rows = []

    for log_path in sorted(ANC_DIR.glob("*/*.iqtree")):
        gene = log_path.stem
        parsed = parse_iqtree_log(log_path)
        if parsed["total_tree_length"] is None:
            continue
        gene_spp = fasta_species(gene)
        spp_tree_len = species_tree_length_for_gene(spp_tree, gene_spp)
        rate = (parsed["total_tree_length"] / spp_tree_len
                if spp_tree_len and spp_tree_len > 0 else None)
        rows.append({
            "gene":               gene,
            "genome":             "mtDNA" if gene in MT_GENES else "nucDNA",
            "total_tree_length":  parsed["total_tree_length"],
            "best_fit_model":     parsed["best_fit_model"],
            "gamma_alpha":        parsed["gamma_alpha"] if parsed["gamma_alpha"] is not None else "NA",
            "subst_rate_per_mya": rate if rate is not None else "NA",
            "n_species":          len(gene_spp),
        })

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} gene rate records to {OUT_CSV}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run script then tests**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav python src/phylo_v2/03_per_gene_rate_calibration.py

docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav pytest tests/test_phylo_v2_rates.py -v
```
Expected: all PASS.

---

## Task 8: `04_age_comparison.py` — mt vs nuc mixed-effects regression

**Files:**
- Create: `src/phylo_v2/04_age_comparison.py`
- Output: `results/phylo_v2/age_comparison.csv`, `results/phylo_v2/age_comparison_summary.txt`

Model: `log(age_mya) ~ is_mt + log(subst_rate) + (1|gene)` using `statsmodels` MixedLM.

- [ ] **Step 1: Implement**

```python
# src/phylo_v2/04_age_comparison.py
"""
Mixed-effects age comparison: mtDNA vs nucDNA cDAVs.

Inputs:
    results/phylo_v2/simmap/*_simmap_posterior.csv  (age_mya_median)
    results/phylo_v2/per_gene_rates.csv             (subst_rate_per_mya)
    results/phylo_v2/character_matrices/_manifest.csv (cdav_gene, cdav_genome)

Output:
    results/phylo_v2/age_comparison.csv
    results/phylo_v2/age_comparison_summary.txt

Run:
    python src/phylo_v2/04_age_comparison.py
"""
from __future__ import annotations
import csv
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[2]
PHY2 = ROOT / "results" / "phylo_v2"
SIMMAP_DIR    = PHY2 / "simmap"
RATES_CSV     = PHY2 / "per_gene_rates.csv"
MANIFEST_CSV  = PHY2 / "character_matrices" / "_manifest.csv"
OUT_DATA      = PHY2 / "age_comparison.csv"
OUT_SUMMARY   = PHY2 / "age_comparison_summary.txt"


def load_simmap_posteriors() -> pd.DataFrame:
    rows = []
    for p in sorted(SIMMAP_DIR.glob("*_simmap_posterior.csv")):
        with open(p) as f:
            row = next(csv.DictReader(f), None)
            if row:
                rows.append(row)
    return pd.DataFrame(rows)


def main():
    posteriors = load_simmap_posteriors()
    if posteriors.empty:
        print("No SIMMAP posteriors found. Run 02_simmap_origins_runner.py first.")
        return

    manifest = pd.read_csv(MANIFEST_CSV)
    rates    = pd.read_csv(RATES_CSV)

    df = (posteriors
          .merge(manifest[["cdav_id", "cdav_gene", "cdav_genome"]], on="cdav_id", how="left")
          .merge(rates[["gene", "subst_rate_per_mya"]].rename(columns={"gene": "cdav_gene"}),
                 on="cdav_gene", how="left"))

    for col in ["age_mya_median", "subst_rate_per_mya"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    n_before = len(df)
    df = df.dropna(subset=["age_mya_median", "subst_rate_per_mya"])
    df = df[df["age_mya_median"] > 0]
    print(f"Dropped {n_before - len(df)} rows with missing/zero age or rate; {len(df)} remain")

    df["log_age"]  = np.log(df["age_mya_median"])
    df["log_rate"] = np.log(df["subst_rate_per_mya"])
    df["is_mt"]    = (df["cdav_genome"] == "mtDNA").astype(int)
    df.to_csv(OUT_DATA, index=False)

    try:
        result = smf.mixedlm("log_age ~ is_mt + log_rate", df, groups=df["cdav_gene"]).fit(reml=True)
        summary_text = str(result.summary())
        coef_mt = result.params.get("is_mt", float("nan"))
        pval_mt = result.pvalues.get("is_mt", float("nan"))
        interp = (f"mtDNA cDAVs are {'older' if coef_mt > 0 else 'younger'} than nucDNA "
                  f"(β_mt={coef_mt:.3f}, p={pval_mt:.4f}) after controlling for substitution rate.")
    except Exception as exc:
        summary_text = f"Model failed: {exc}"
        interp = "Could not fit mixed-effects model."

    with open(OUT_SUMMARY, "w") as f:
        f.write("mtDNA vs nucDNA cDAV Age Comparison\n" + "=" * 60 + "\n\n")
        f.write(f"N cDAVs in model: {len(df)}\n")
        f.write(f"mtDNA: {(df['cdav_genome']=='mtDNA').sum()}  nucDNA: {(df['cdav_genome']=='nucDNA').sum()}\n\n")
        f.write("Model: log(age_mya) ~ is_mt + log(subst_rate) + (1|gene)\n\n")
        f.write(summary_text + "\n\nInterpretation:\n" + interp + "\n")

    print(f"Written: {OUT_SUMMARY}")
    print(interp)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav python src/phylo_v2/04_age_comparison.py
```
Expected: writes `age_comparison_summary.txt` with model results.

---

## Task 9: `05_resolvability_filter.py`

**Files:**
- Create: `src/phylo_v2/05_resolvability_filter.py`
- Test: `tests/test_phylo_v2_resolvability.py`
- Output: `results/phylo_v2/resolvable_pairs.csv`, `results/phylo_v2/resolvability_filter_log.txt`

A pair is testable if: ≥20 species have non-NA state for both columns AND ≥3 species have state=1 for each column. Do not loosen these criteria. Report exact counts at each filter step.

- [ ] **Step 1: Write tests**

```python
# tests/test_phylo_v2_resolvability.py
import csv
from pathlib import Path

RESOLVABLE = Path("results/phylo_v2/resolvable_pairs.csv")
FILTER_LOG = Path("results/phylo_v2/resolvability_filter_log.txt")


def test_resolvable_csv_exists():
    assert RESOLVABLE.exists(), "Run 05_resolvability_filter.py first"


def test_filter_log_exists():
    assert FILTER_LOG.exists()


def test_resolvable_has_required_columns():
    with open(RESOLVABLE) as f:
        cols = csv.DictReader(f).fieldnames
    required = [
        "ann_id", "dar_gene", "dar_genome", "dar_aa_coord", "dar_alt_aa",
        "contact_gene", "contact_refseq_pos", "contact_alt_aa",
        "n_species_both_known", "n_cdav_state1", "n_contact_state1",
        "resolvable",
    ]
    for c in required:
        assert c in cols, f"Missing column: {c}"


def test_all_resolvable_rows_pass_criteria():
    with open(RESOLVABLE) as f:
        for row in csv.DictReader(f):
            if row["resolvable"] != "TRUE":
                continue
            assert int(row["n_species_both_known"]) >= 20
            assert int(row["n_cdav_state1"]) >= 3
            assert int(row["n_contact_state1"]) >= 3


def test_filter_log_reports_counts():
    text = FILTER_LOG.read_text()
    assert "resolvable" in text.lower()
```

- [ ] **Step 2: Implement**

```python
# src/phylo_v2/05_resolvability_filter.py
"""
Filter tested pairs to those where temporal ordering is resolvable by corHMM.

Criteria (all must be met):
  1. ≥3 species with cdav_state=1
  2. ≥3 species with contact_state=1
  3. ≥20 species with non-NA state for BOTH columns

Outputs:
    results/phylo_v2/resolvable_pairs.csv
    results/phylo_v2/resolvability_filter_log.txt

Run:
    python src/phylo_v2/05_resolvability_filter.py
"""
from __future__ import annotations
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PHY2 = ROOT / "results" / "phylo_v2"
MATRIX_DIR = PHY2 / "character_matrices"
PAIRS_CSV  = ROOT / "results" / "structural" / "all_tested_pairs.csv"
OUT_CSV    = PHY2 / "resolvable_pairs.csv"
OUT_LOG    = PHY2 / "resolvability_filter_log.txt"

MIN_SPECIES = 20
MIN_STATE1  = 3


def contact_col_name(pair: dict) -> str:
    return f"contact_{pair['contact_gene']}_{pair['contact_refseq_pos']}_{pair['contact_alt_aa']}"


def load_matrix(cdav_id: str) -> list[dict] | None:
    path = MATRIX_DIR / f"{cdav_id}.tsv"
    if not path.exists():
        return None
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def assess_pair(pair: dict, matrix_rows: list[dict]) -> dict:
    col = contact_col_name(pair)
    both_known = [
        r for r in matrix_rows
        if r.get("cdav_state") in ("0", "1") and r.get(col) in ("0", "1")
    ]
    n_both  = len(both_known)
    n_cdav1 = sum(1 for r in both_known if r["cdav_state"] == "1")
    n_con1  = sum(1 for r in both_known if r.get(col) == "1")
    return {
        **{k: pair[k] for k in [
            "ann_id", "dar_gene", "dar_genome", "dar_aa_coord",
            "dar_ref_aa", "dar_alt_aa", "contact_gene",
            "contact_refseq_pos", "contact_human_aa", "contact_alt_aa",
            "contact_class", "contact_type", "tier",
        ]},
        "contact_col":            col,
        "n_species_both_known":   n_both,
        "n_cdav_state1":          n_cdav1,
        "n_contact_state1":       n_con1,
        "resolvable":             str(
            n_both >= MIN_SPECIES and n_cdav1 >= MIN_STATE1 and n_con1 >= MIN_STATE1
        ).upper(),
    }


def main():
    with open(PAIRS_CSV) as f:
        all_pairs = list(csv.DictReader(f))

    results = []
    n_no_matrix = 0
    n_no_col    = 0
    matrix_cache: dict[str, list[dict] | None] = {}

    for pair in all_pairs:
        cid = pair["ann_id"]
        if cid not in matrix_cache:
            matrix_cache[cid] = load_matrix(cid)
        mat = matrix_cache[cid]
        if mat is None:
            n_no_matrix += 1
            continue
        col = contact_col_name(pair)
        if mat and col not in mat[0]:
            n_no_col += 1
            continue
        results.append(assess_pair(pair, mat))

    n_resolvable = sum(1 for r in results if r["resolvable"] == "TRUE")

    PHY2.mkdir(parents=True, exist_ok=True)
    if results:
        with open(OUT_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)

    log_lines = [
        "Resolvability filter report",
        f"Total pairs tested:              {len(all_pairs)}",
        f"  Missing matrix:                {n_no_matrix}",
        f"  Contact col absent in matrix:  {n_no_col}",
        f"  Assessed:                      {len(results)}",
        f"  Resolvable (all criteria met): {n_resolvable}",
        "",
        f"Criteria: n_both>={MIN_SPECIES}, n_cdav1>={MIN_STATE1}, n_con1>={MIN_STATE1}",
        "These criteria were not loosened to increase sample size.",
    ]
    with open(OUT_LOG, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    for line in log_lines:
        print(line)


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run script then tests**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav python src/phylo_v2/05_resolvability_filter.py

docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav pytest tests/test_phylo_v2_resolvability.py -v
```
Expected: all PASS. Report the resolvable count honestly — do not go back and loosen criteria.

---

## Task 10: `06_corhmm_joint_models.R` — dependent vs independent ordering

**Files:**
- Create: `src/phylo_v2/06_corhmm_joint_models.R`

CLI: `Rscript 06_corhmm_joint_models.R <matrix_tsv> <gene_tree_nwk> <cdav_col> <contact_col> <out_dir>`

States: 00 (neither), 01 (contact only), 10 (cDAV only), 11 (both).
Fits independent (product of two 2-state ER chains) vs dependent (full 4-state ARD with only Hamming-1 transitions allowed) corHMM models.
Reports AIC, ΔAIC, and marginal-reconstruction ordering counts for all four categories.

- [ ] **Step 1: Write the R script**

```r
#!/usr/bin/env Rscript
# src/phylo_v2/06_corhmm_joint_models.R
#
# Fit independent and dependent evolution models for one cDAV-contact pair.
#
# Usage:
#   Rscript 06_corhmm_joint_models.R \
#     <matrix_tsv> <gene_tree_nwk> <cdav_col> <contact_col> <out_dir>

suppressMessages({
  library(ape)
  library(corHMM)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript 06_corhmm_joint_models.R <matrix_tsv> <gene_tree_nwk> <cdav_col> <contact_col> <out_dir>")
}

matrix_tsv  <- args[1]
tree_nwk    <- args[2]
cdav_col    <- args[3]
contact_col <- args[4]
out_dir     <- args[5]

mat  <- read.delim(matrix_tsv, stringsAsFactors = FALSE, na.strings = "NA")
tree <- read.tree(tree_nwk)

pair_id <- paste0(sub("\\.tsv$", "", basename(matrix_tsv)),
                  "_vs_", gsub("[^A-Za-z0-9]", "_", contact_col))

# ── Build 4-state character vector ────────────────────────────────────────────
valid  <- mat[!is.na(mat[[cdav_col]]) & !is.na(mat[[contact_col]]), ]
shared <- intersect(valid$species, tree$tip.label)

if (length(shared) < 10) {
  cat(sprintf("SKIP %s: only %d shared species\n", pair_id, length(shared)))
  quit(status = 0)
}

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, shared))

idx      <- match(shared, valid$species)
state_int <- as.integer(paste0(valid[[cdav_col]][idx], valid[[contact_col]][idx]))
# Recode compound integer to 1-based index: 00→1, 01→2, 10→3, 11→4
state_int <- ifelse(state_int == 0, 1,
             ifelse(state_int == 1, 2,
             ifelse(state_int == 10, 3, 4)))
names(state_int) <- shared

state_df <- data.frame(Genus_sp = shared, state = state_int[tree_pruned$tip.label],
                        stringsAsFactors = FALSE)

# ── Rate matrices ─────────────────────────────────────────────────────────────
# State order: 1=00, 2=01, 3=10, 4=11
# Independent: contact transitions (00↔01, 10↔11) share rate 1;
#              cDAV transitions    (00↔10, 01↔11) share rate 2
rate_mat_indep <- matrix(0, 4, 4)
rate_mat_indep[1,2] <- 1; rate_mat_indep[2,1] <- 1
rate_mat_indep[3,4] <- 1; rate_mat_indep[4,3] <- 1
rate_mat_indep[1,3] <- 2; rate_mat_indep[3,1] <- 2
rate_mat_indep[2,4] <- 2; rate_mat_indep[4,2] <- 2

# Dependent: all Hamming-1 transitions are free (ARD)
state_names <- c("00","01","10","11")
rate_mat_dep <- matrix(0, 4, 4)
rate_idx <- 1L
for (i in 1:4) for (j in 1:4) {
  if (i != j) {
    si <- strsplit(state_names[i], "")[[1]]
    sj <- strsplit(state_names[j], "")[[1]]
    if (sum(si != sj) == 1) {
      rate_mat_dep[i, j] <- rate_idx
      rate_idx <- rate_idx + 1L
    }
  }
}

fit_model <- function(rate_mat, label) {
  tryCatch(
    corHMM(phy = tree_pruned, data = state_df, rate.cat = 1,
           rate.mat = rate_mat, node.states = "marginal",
           get.tip.states = FALSE, quiet = TRUE),
    error = function(e) { message(sprintf("corHMM %s failed: %s", label, e$message)); NULL }
  )
}

fit_indep <- fit_model(rate_mat_indep, "independent")
fit_dep   <- fit_model(rate_mat_dep,   "dependent")

aic_indep  <- if (is.null(fit_indep)) NA_real_ else fit_indep$AIC
aic_dep    <- if (is.null(fit_dep))   NA_real_ else fit_dep$AIC
delta_aic  <- aic_dep - aic_indep
evidence   <- if (is.na(delta_aic))   "insufficient_data"
              else if (delta_aic < -10) "strong_dependent"
              else if (delta_aic < -2)  "weak_dependent"
              else if (delta_aic >  2)  "independent"
              else                      "equivocal"

# ── Infer ordering from marginal ASR of the dependent model ───────────────────
infer_ordering <- function(fit) {
  if (is.null(fit)) return(list(contact_first=NA, cdav_first=NA, co_occurring=NA, indeterminate=NA))
  node_states <- fit$states
  map_state   <- apply(node_states, 1, which.max)
  edge_mat    <- fit$phy$edge
  n_tips      <- length(fit$phy$tip.label)
  n_cf <- 0L; n_df <- 0L; n_co <- 0L; n_in <- 0L
  for (e in seq_len(nrow(edge_mat))) {
    par <- edge_mat[e, 1]; chi <- edge_mat[e, 2]
    par_s <- if (par > n_tips) state_names[map_state[par - n_tips]] else NA
    chi_s <- if (chi > n_tips) state_names[map_state[chi - n_tips]] else NA
    if (is.na(par_s) || is.na(chi_s)) next
    if (par_s == "00" && chi_s == "11") { n_co <- n_co + 1L; next }
    if ((par_s == "00" && chi_s == "01") || (par_s == "10" && chi_s == "11")) {
      n_cf <- n_cf + 1L; next
    }
    if ((par_s == "00" && chi_s == "10") || (par_s == "01" && chi_s == "11")) {
      n_df <- n_df + 1L; next
    }
    n_in <- n_in + 1L
  }
  list(contact_first=n_cf, cdav_first=n_df, co_occurring=n_co, indeterminate=n_in)
}

ordering <- infer_ordering(fit_dep)

result_df <- data.frame(
  pair_id          = pair_id,
  n_species        = length(shared),
  aic_independent  = aic_indep,
  aic_dependent    = aic_dep,
  delta_aic        = delta_aic,
  evidence_dep     = evidence,
  n_contact_first  = ordering$contact_first,
  n_cdav_first     = ordering$cdav_first,
  n_co_occurring   = ordering$co_occurring,
  n_indeterminate  = ordering$indeterminate,
  dominant_ordering = {
    vals <- c(contact_first=ordering$contact_first,
              cdav_first=ordering$cdav_first,
              co_occurring=ordering$co_occurring)
    if (all(is.na(vals))) "insufficient_data" else names(which.max(vals))
  }
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(result_df, file.path(out_dir, paste0(pair_id, "_corhmm.csv")), row.names = FALSE)

cat(sprintf("Done: %s | ΔAIC=%.2f (%s) | ordering=%s\n",
            pair_id, ifelse(is.na(delta_aic), 0, delta_aic),
            evidence, result_df$dominant_ordering))
```

- [ ] **Step 2: Smoke-test with a single resolvable pair**

```bash
# Get first resolvable pair from the CSV
read ANN_ID GENE CONTACT_COL < <(awk -F',' '
  NR>1 && $NF=="TRUE" { print $1, $2, $10; exit }
' results/phylo_v2/resolvable_pairs.csv)

docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav Rscript src/phylo_v2/06_corhmm_joint_models.R \
    results/phylo_v2/character_matrices/${ANN_ID}.tsv \
    data/phylo/iqtree_jobs/${GENE}/${GENE}_tree.nwk \
    cdav_state \
    "${CONTACT_COL}" \
    results/phylo_v2/corhmm/
```
Expected: prints "Done: ... | ΔAIC=... | ordering=..."

---

## Task 11: `06_corhmm_runner.py` — parallel corHMM orchestrator

**Files:**
- Create: `src/phylo_v2/06_corhmm_runner.py`

- [ ] **Step 1: Implement**

```python
# src/phylo_v2/06_corhmm_runner.py
"""
Parallel runner for 06_corhmm_joint_models.R.

Usage:
    python src/phylo_v2/06_corhmm_runner.py [--max-pairs N] [--workers N]
"""
from __future__ import annotations
import argparse
import csv
import logging
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
PHY2 = ROOT / "results" / "phylo_v2"
RESOLVABLE = PHY2 / "resolvable_pairs.csv"
MATRIX_DIR = PHY2 / "character_matrices"
IQTREE_DIR = ROOT / "data" / "phylo" / "iqtree_jobs"
CORHMM_DIR = PHY2 / "corhmm"
R_SCRIPT   = ROOT / "src" / "phylo_v2" / "06_corhmm_joint_models.R"


def run_one_pair(row: dict) -> str | None:
    ann_id      = row["ann_id"]
    gene        = row["dar_gene"]
    contact_col = row["contact_col"]
    matrix_tsv  = MATRIX_DIR / f"{ann_id}.tsv"
    gene_tree   = IQTREE_DIR / gene / f"{gene}_tree.nwk"
    if not matrix_tsv.exists() or not gene_tree.exists():
        return None
    result = subprocess.run(
        ["conda", "run", "-n", "oxphos_dav",
         "Rscript", "--vanilla", str(R_SCRIPT),
         str(matrix_tsv), str(gene_tree),
         "cdav_state", contact_col, str(CORHMM_DIR)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        log.error("R failed for %s / %s:\n%s", ann_id, contact_col, result.stderr[-400:])
        return None
    return result.stdout.strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-pairs", type=int, default=None)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    CORHMM_DIR.mkdir(parents=True, exist_ok=True)

    with open(RESOLVABLE) as f:
        pairs = [r for r in csv.DictReader(f) if r.get("resolvable") == "TRUE"]
    if args.max_pairs:
        pairs = pairs[: args.max_pairs]
    log.info("Running corHMM for %d resolvable pairs", len(pairs))

    completed = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one_pair, row): row for row in pairs}
        for fut in as_completed(futures):
            try:
                msg = fut.result()
                if msg:
                    completed += 1
                    log.info(msg)
            except Exception as exc:
                log.error("Exception: %s", exc)

    log.info("Completed %d / %d corHMM runs", completed, len(pairs))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run on small subset**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav python src/phylo_v2/06_corhmm_runner.py \
    --max-pairs 10 --workers 4
```
Expected: 10 `*_corhmm.csv` files in `results/phylo_v2/corhmm/`.

---

## Task 12: `07_independence_null.py` — tree-aware null for contact_first proportion

**Files:**
- Create: `src/phylo_v2/07_independence_null.py`
- Test: `tests/test_phylo_v2_null.py`
- Output: `results/phylo_v2/independence_null.csv`

Permutes contact_state tip labels (preserving marginal frequency) to generate a null distribution of the `contact_first / (contact_first + cdav_first + co_occurring)` fraction. The observed fraction comes from the corHMM result. N_PERM=200 default.

- [ ] **Step 1: Write tests**

```python
# tests/test_phylo_v2_null.py
import csv
from pathlib import Path

NULL_CSV = Path("results/phylo_v2/independence_null.csv")


def test_null_csv_exists():
    assert NULL_CSV.exists(), "Run 07_independence_null.py first"


def test_null_has_required_columns():
    with open(NULL_CSV) as f:
        cols = csv.DictReader(f).fieldnames
    required = [
        "pair_id", "n_perms", "obs_contact_first_frac",
        "null_median", "null_lower95", "null_upper95", "pval_contact_first",
    ]
    for c in required:
        assert c in cols, f"Missing: {c}"


def test_pval_bounded():
    with open(NULL_CSV) as f:
        for row in csv.DictReader(f):
            if row["pval_contact_first"] not in ("NA", ""):
                p = float(row["pval_contact_first"])
                assert 0 <= p <= 1


def test_na_obs_implies_na_pval():
    with open(NULL_CSV) as f:
        for row in csv.DictReader(f):
            if row["obs_contact_first_frac"] == "NA":
                assert row["pval_contact_first"] == "NA"
```

- [ ] **Step 2: Implement**

```python
# src/phylo_v2/07_independence_null.py
"""
Tree-aware independence null for contact_first temporal ordering.

For each resolvable pair with a corHMM result, permutes contact_state
tip labels N_PERM times, refits corHMM, and compares the observed
contact_first fraction to the null distribution.

Output: results/phylo_v2/independence_null.csv

Run:
    python src/phylo_v2/07_independence_null.py [--n-perms N] [--max-pairs N]
"""
from __future__ import annotations
import argparse
import copy
import csv
import logging
import random
import subprocess
import tempfile
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
PHY2 = ROOT / "results" / "phylo_v2"
CORHMM_DIR = PHY2 / "corhmm"
MATRIX_DIR = PHY2 / "character_matrices"
IQTREE_DIR = ROOT / "data" / "phylo" / "iqtree_jobs"
RESOLVABLE = PHY2 / "resolvable_pairs.csv"
R_SCRIPT   = ROOT / "src" / "phylo_v2" / "06_corhmm_joint_models.R"
OUT_CSV    = PHY2 / "independence_null.csv"

N_PERM_DEFAULT = 200


def contact_first_frac(row: dict) -> float | None:
    try:
        n_cf  = int(row["n_contact_first"])
        n_df  = int(row["n_cdav_first"])
        n_co  = int(row["n_co_occurring"])
        total = n_cf + n_df + n_co
        return n_cf / total if total > 0 else None
    except (TypeError, ValueError, KeyError):
        return None


def load_corhmm_result(ann_id: str, contact_col: str) -> dict | None:
    """Find the corHMM output file for this pair."""
    safe_col = contact_col.replace("/", "_").replace("-", "_")
    for p in CORHMM_DIR.glob(f"{ann_id}_vs_*_corhmm.csv"):
        if safe_col in p.name or contact_col in p.name:
            with open(p) as f:
                return next(csv.DictReader(f), None)
    # fallback: any file starting with ann_id
    candidates = list(CORHMM_DIR.glob(f"{ann_id}*_corhmm.csv"))
    if candidates:
        with open(candidates[0]) as f:
            return next(csv.DictReader(f), None)
    return None


def permute_contact_state(matrix_rows: list[dict], contact_col: str) -> list[dict]:
    rows = copy.deepcopy(matrix_rows)
    non_na_idx = [i for i, r in enumerate(rows) if r.get(contact_col) in ("0", "1")]
    values = [rows[i][contact_col] for i in non_na_idx]
    random.shuffle(values)
    for idx, val in zip(non_na_idx, values):
        rows[idx][contact_col] = val
    return rows


def run_corhmm_on_rows(matrix_rows: list[dict], gene: str,
                       ann_id: str, contact_col: str, tmp_dir: Path) -> dict | None:
    gene_tree = IQTREE_DIR / gene / f"{gene}_tree.nwk"
    if not gene_tree.exists():
        return None
    tmp_tsv = tmp_dir / f"{ann_id}_perm.tsv"
    cols = list(matrix_rows[0].keys())
    with open(tmp_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(matrix_rows)
    perm_out = tmp_dir / "out"
    perm_out.mkdir(exist_ok=True)
    result = subprocess.run(
        ["conda", "run", "-n", "oxphos_dav",
         "Rscript", "--vanilla", str(R_SCRIPT),
         str(tmp_tsv), str(gene_tree), "cdav_state", contact_col, str(perm_out)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return None
    out_files = list(perm_out.glob("*_corhmm.csv"))
    if not out_files:
        return None
    with open(out_files[0]) as f:
        return next(csv.DictReader(f), None)


def run_null_for_pair(pair_row: dict, n_perms: int) -> dict:
    ann_id      = pair_row["ann_id"]
    gene        = pair_row["dar_gene"]
    contact_col = pair_row["contact_col"]

    matrix_tsv = MATRIX_DIR / f"{ann_id}.tsv"
    if not matrix_tsv.exists():
        return {"pair_id": f"{ann_id}__{contact_col}", "n_perms": n_perms,
                "n_null_successful": 0, "obs_contact_first_frac": "NA",
                "null_median": "NA", "null_lower95": "NA",
                "null_upper95": "NA", "pval_contact_first": "NA"}

    with open(matrix_tsv) as f:
        matrix_rows = list(csv.DictReader(f, delimiter="\t"))

    obs_row  = load_corhmm_result(ann_id, contact_col)
    obs_frac = contact_first_frac(obs_row) if obs_row else None

    null_fracs = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        for _ in range(n_perms):
            perm_rows   = permute_contact_state(matrix_rows, contact_col)
            perm_result = run_corhmm_on_rows(perm_rows, gene, ann_id, contact_col, tmp_path)
            frac        = contact_first_frac(perm_result) if perm_result else None
            if frac is not None:
                null_fracs.append(frac)

    if not null_fracs or obs_frac is None:
        return {"pair_id": f"{ann_id}__{contact_col}", "n_perms": n_perms,
                "n_null_successful": len(null_fracs),
                "obs_contact_first_frac": "NA", "null_median": "NA",
                "null_lower95": "NA", "null_upper95": "NA",
                "pval_contact_first": "NA"}

    null_fracs.sort()
    n = len(null_fracs)
    return {
        "pair_id":                f"{ann_id}__{contact_col}",
        "n_perms":                n_perms,
        "n_null_successful":      n,
        "obs_contact_first_frac": obs_frac,
        "null_median":            null_fracs[n // 2],
        "null_lower95":           null_fracs[int(0.025 * n)],
        "null_upper95":           null_fracs[int(0.975 * n)],
        "pval_contact_first":     sum(1 for f in null_fracs if f >= obs_frac) / n,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-perms", type=int, default=N_PERM_DEFAULT)
    parser.add_argument("--max-pairs", type=int, default=None)
    args = parser.parse_args()

    with open(RESOLVABLE) as f:
        pairs = [r for r in csv.DictReader(f) if r.get("resolvable") == "TRUE"]
    if args.max_pairs:
        pairs = pairs[: args.max_pairs]

    log.info("Running null for %d pairs (%d perms each)", len(pairs), args.n_perms)
    results = []
    for i, row in enumerate(pairs):
        log.info("[%d/%d] %s", i + 1, len(pairs), row["ann_id"])
        results.append(run_null_for_pair(row, args.n_perms))

    if results:
        with open(OUT_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)
    log.info("Written: %s", OUT_CSV)


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run on 3 pairs then test**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav python src/phylo_v2/07_independence_null.py \
    --n-perms 50 --max-pairs 3

docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav pytest tests/test_phylo_v2_null.py -v
```
Expected: all PASS.

---

## Task 13: `08_compile_phylo_results.py` — final results table

**Files:**
- Create: `src/phylo_v2/08_compile_phylo_results.py`
- Output: `results/phylo_v2/phylo_v2_results.csv`

Joins SIMMAP posteriors, corHMM results, and independence null p-values onto `all_tested_pairs.csv`.

- [ ] **Step 1: Implement**

```python
# src/phylo_v2/08_compile_phylo_results.py
"""
Compile all phylo_v2 results into a single annotated pair table.

Output: results/phylo_v2/phylo_v2_results.csv

Run:
    python src/phylo_v2/08_compile_phylo_results.py
"""
from __future__ import annotations
import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PHY2 = ROOT / "results" / "phylo_v2"

PAIRS_CSV   = ROOT / "results" / "structural" / "all_tested_pairs.csv"
SIMMAP_DIR  = PHY2 / "simmap"
CORHMM_DIR  = PHY2 / "corhmm"
NULL_CSV    = PHY2 / "independence_null.csv"
RESOLVABLE  = PHY2 / "resolvable_pairs.csv"
OUT_CSV     = PHY2 / "phylo_v2_results.csv"


def load_simmap_posteriors() -> pd.DataFrame:
    rows = []
    for p in SIMMAP_DIR.glob("*_simmap_posterior.csv"):
        with open(p) as f:
            row = next(csv.DictReader(f), None)
            if row:
                rows.append(row)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    return df.rename(columns={
        "n_origins_median":    "simmap_n_origins_median",
        "n_origins_lower95":   "simmap_n_origins_lower95",
        "n_origins_upper95":   "simmap_n_origins_upper95",
        "n_origins_prob_multi": "simmap_prob_multi_origin",
        "age_mya_median":      "simmap_age_mya_median",
        "age_mya_lower95":     "simmap_age_mya_lower95",
        "age_mya_upper95":     "simmap_age_mya_upper95",
    })


def load_corhmm_results() -> pd.DataFrame:
    rows = []
    for p in CORHMM_DIR.glob("*_corhmm.csv"):
        with open(p) as f:
            row = next(csv.DictReader(f), None)
            if row:
                rows.append(row)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["ann_id"]          = df["pair_id"].str.split("_vs_").str[0]
    df["contact_col_key"] = df["pair_id"].str.split("_vs_").str[1]
    return df.rename(columns={
        "delta_aic":         "corhmm_delta_aic",
        "evidence_dep":      "corhmm_evidence_dep",
        "dominant_ordering": "corhmm_dominant_ordering",
        "n_contact_first":   "corhmm_n_contact_first",
        "n_cdav_first":      "corhmm_n_cdav_first",
        "n_co_occurring":    "corhmm_n_co_occurring",
    })


def main():
    pairs = pd.read_csv(PAIRS_CSV)

    simmap = load_simmap_posteriors()
    if not simmap.empty:
        simmap_cols = ["cdav_id"] + [c for c in simmap.columns if c.startswith("simmap_")]
        pairs = pairs.merge(simmap[simmap_cols].rename(columns={"cdav_id": "ann_id"}),
                            on="ann_id", how="left")

    if RESOLVABLE.exists():
        resolvable = pd.read_csv(RESOLVABLE)
        corhmm = load_corhmm_results()
        if not corhmm.empty:
            resolvable["contact_col_key"] = (
                resolvable["contact_col"].str.replace("/", "_").str.replace("-", "_")
            )
            corhmm_cols = [c for c in corhmm.columns if c.startswith("corhmm_")]
            merged_r = resolvable.merge(
                corhmm[["ann_id", "contact_col_key"] + corhmm_cols],
                on=["ann_id", "contact_col_key"], how="left"
            )
            pairs = pairs.merge(
                merged_r[["ann_id", "contact_gene", "contact_refseq_pos",
                           "contact_alt_aa"] + corhmm_cols],
                on=["ann_id", "contact_gene", "contact_refseq_pos", "contact_alt_aa"],
                how="left",
            )

    if NULL_CSV.exists():
        null_df = pd.read_csv(NULL_CSV)
        null_df[["ann_id", "_"]] = null_df["pair_id"].str.split("__", n=1, expand=True)
        null_summary = (null_df.groupby("ann_id")
                        .agg(null_pval_contact_first_min=("pval_contact_first", "min"))
                        .reset_index())
        pairs = pairs.merge(null_summary, on="ann_id", how="left")

    PHY2.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(OUT_CSV, index=False)
    new_cols = [c for c in pairs.columns if c.startswith(("simmap_", "corhmm_", "null_"))]
    print(f"Wrote {len(pairs)} rows to {OUT_CSV}")
    print(f"New phylo_v2 columns ({len(new_cols)}): {new_cols}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run**

```bash
docker run --rm -v $(pwd):/app oxphos_dav_analysis \
  conda run -n oxphos_dav python src/phylo_v2/08_compile_phylo_results.py
```
Expected: `results/phylo_v2/phylo_v2_results.csv` with all `simmap_*` and `corhmm_*` columns present.

---

## Spec coverage self-review

| Requirement | Task(s) |
|---|---|
| SIMMAP with posterior distributions (not point estimates) | Tasks 5, 6 |
| Independent-evolution null for ordering proportions | Task 12 |
| Origin counts per cDAV | Tasks 5, 6 |
| Origin ages, mt vs nuc comparison controlling for rate | Tasks 7, 8 |
| Temporal ordering — all four categories preserved | Task 10 |
| corHMM dependent vs independent joint models | Tasks 10, 11 |
| Species absent from tree propagated as NA | Task 4 |
| ≥20 species threshold flagged but not deleted | Task 4 |
| Resolvability criteria not loosened | Task 9 |
| Per-simulation files preserved | Tasks 5, 10 — individual CSVs per cDAV/pair |
| Dockerfile additions | Task 1 |
| Old `src/phylo/` untouched | Not modified anywhere |
