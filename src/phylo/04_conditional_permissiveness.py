#!/usr/bin/env python3
"""
src/phylo/04_conditional_permissiveness.py

Corrected conditional permissiveness analysis for derived cDAVs.

Addresses the circular-logic problem in OR = 24.88 (contact_first enrichment):
  For a derived cDAV on branch X with contact_first = True, ALL descendant
  species inherit the contact alt-AA from the parent node — making
  n_cdav_with_alt ≈ n_cdav_spp by ancestry alone. Both branch_cooccur_fdr
  and fisher_fdr are mechanically inflated by this inheritance, so the naive
  OR reflects phylogenetic structure rather than functional coupling.

Three analyses (in increasing computational cost):

  1. Fisher-only sensitivity (fast, always runs)
     Re-tests contact_first enrichment in derived cDAVs using ONLY fisher_fdr
     (no branch_cooccur_fdr) as the significance criterion. Still partially
     confounded by inheritance-driven Fisher OR; reported as sensitivity only.

  2. Multi-origin binomial test (fast, always runs)
     For DARs with n_dar_gain_branches ≥ 2: compares fraction of gain events
     on contact+ background branches to the branch-length-weighted expected
     rate from the ancestral state reconstruction. Binomial test against this
     null. Does not use per-clade count statistics.

  3. Tree-aware conditional permutation (slow; requires pyvolve)
     For each derived cDAV with contact_first = True:
       a. Simulate contact-column evolution on the gene tree 10,000 times using
          the WAG model + empirical column frequencies (no coupling to the cDAV).
       b. In each simulated alignment, count how many cDAV-clade species carry
          the contact alt-AA → compute simulated OR.
       c. Empirical p-value = fraction of sims with OR ≥ observed OR.
     p_contact estimated from NON-cDAV species only (n_bg_with_alt / n_bg_spp).
     Simulations batched per unique (contact_gene, contact_refseq_pos) for
     efficiency (same tree + column → identical model).

Inputs:
  results/phylo/timing_annotations.csv
  results/structural/all_tested_pairs.csv
  data/phylo/ancestral_state_maps.json
  data/phylo/iqtree_jobs/{gene}/{gene}_tree.nwk  (for pyvolve)
  data/alignments/{genome_dir}/{gene}_aa_alignment.fasta  (for column freqs)

Outputs:
  results/phylo/contact_first_revised_test.csv   — fisher-only sensitivity
  results/phylo/multi_origin_binomial.csv        — per-DAR binomial test
  results/phylo/conditional_permissiveness.csv   — per-pair pyvolve perm result

Usage:
  docker run --rm -v $(pwd):/app oxphos_dav_analysis conda run -n oxphos_dav \\
      python src/phylo/04_conditional_permissiveness.py [--n-perm N] [--skip-pyvolve]

  Use --skip-pyvolve to run only the fast analyses (1 and 2) without pyvolve.
  Use --n-perm 100 for a quick test run; default 1000 for production.
  HPC recommended for full --n-perm 10000 run.
"""

import argparse
import csv
import json
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, fisher_exact

ROOT = Path(__file__).resolve().parents[2]

# ─── Paths ────────────────────────────────────────────────────────────────────
TIMING_CSV   = ROOT / "results" / "phylo" / "timing_annotations.csv"
ALL_PAIRS    = ROOT / "results" / "structural" / "all_tested_pairs.csv"
ASM_JSON     = ROOT / "data" / "phylo" / "ancestral_state_maps.json"
IQTREE_DIR   = ROOT / "data" / "phylo" / "iqtree_jobs"
ASR_DIR      = ROOT / "data" / "phylo" / "ancestral_states"
ALN_MTDNA    = ROOT / "data" / "alignments" / "mtdna_aa"
ALN_NUCLEAR  = ROOT / "data" / "alignments" / "toga_hg38_aa"
OUT_DIR      = ROOT / "results" / "phylo"

# Join key columns (must match between timing_annotations and all_tested_pairs)
JOIN_COLS = [
    "ann_id", "dar_gene", "dar_aa_coord", "dar_alt_aa",
    "contact_gene", "contact_refseq_pos", "contact_alt_aa",
]

# Pyvolve amino acid order (alphabetical by IUPAC 1-letter code)
PYVOLVE_AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# IQTree base model name → Pyvolve model name.
# Q.* models are not in pyvolve 1.1.0; fall back to LG for those.
# +I / +G4 / +F rate-variation suffixes are stripped before lookup.
_IQTREE_TO_PYVOLVE: dict[str, str] = {
    "MTVER":   "MTVER",    # vertebrate mitochondrial
    "MTMAM":   "MTMAM",    # mammalian mitochondrial
    "MTREV":   "MTREV24",  # mitochondrial REV
    "MTINV":   "MTINV",
    "MTZOA":   "MTZOA",
    "MTART":   "MTART",
    "MTMET":   "MTMET",
    "LG":      "LG",
    "WAG":     "WAG",
    "JTT":     "JTT",
    "JTTDCMUT":"JTTDCMUT",
    "DAYHOFF": "DAYHOFF",
    "DCMUT":   "DAYHOFFDCMUT",
    "VT":      "VT",
    "RTREV":   "RTREV",
    "CPREV":   "CPREV",
    "HIVB":    "HIVB",
    "HIVW":    "HIVW",
    "AB":      "AB",
    "BLOSUM62":"BLOSUM62",
    "PMB":     "PMB",
}
_PYVOLVE_FALLBACK = "LG"  # for unsupported models (Q.BIRD, Q.MAMMAL, Q.PLANT, etc.)


# ─── Data loading helpers ──────────────────────────────────────────────────────

def load_csv_as_dicts(path: Path) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f))


def row_key(row: dict) -> tuple:
    return tuple(str(row[c]) for c in JOIN_COLS)


def load_fasta(path: Path) -> dict[str, str]:
    """Return {species: sequence} from a gapped protein FASTA."""
    seqs: dict[str, str] = {}
    header = None
    chunks: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header = line[1:].split()[0].split("|")[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        seqs[header] = "".join(chunks)
    return seqs


def aln_dir_for_gene(gene: str, asr_data: dict) -> Path | None:
    """Return alignment directory for a gene (mtDNA or nuclear)."""
    mt_genes = {g for g in asr_data if g.startswith("MT-")}
    if gene in mt_genes:
        return ALN_MTDNA
    return ALN_NUCLEAR


# ─── Ancestral state reconstruction helpers ───────────────────────────────────

def reconstruct_node_states(asm_gene: dict, positions: set[str]) -> dict[str, dict[str, str]]:
    """
    Reconstruct ancestral amino acid states at every tree node for the given
    human protein positions (1-based position strings).

    Returns {node_id: {position_str: amino_acid}}.
    Traverses the tree top-down from the root, inheriting parent states and
    applying documented changes. Both high-confidence (branches) and
    low-confidence (branches_lc) changes are applied; which branch a change
    came from is tracked separately for ASR confidence stratification.
    """
    root_node = asm_gene["root_node"]
    root_states = asm_gene["root_states"]
    all_branches = {**asm_gene["branches"], **asm_gene["branches_lc"]}
    lc_branch_set = set(asm_gene["branches_lc"].keys())
    node_to_children = asm_gene["node_to_children"]

    # Initialise root
    node_states: dict[str, dict[str, str]] = {
        root_node: {p: root_states.get(p, "-") for p in positions}
    }
    # Track which nodes' states came (at least partially) from lc branches
    node_is_lc: dict[str, bool] = {root_node: False}

    # BFS
    from collections import deque
    queue = deque([root_node])
    while queue:
        parent = queue.popleft()
        for child in node_to_children.get(parent, []):
            branch_key = f"{parent}|{child}"
            changes = all_branches.get(branch_key, {})
            is_lc_branch = branch_key in lc_branch_set

            child_states = dict(node_states[parent])  # inherit parent
            for pos_str, (from_aa, to_aa) in changes.items():
                if pos_str in positions:
                    child_states[pos_str] = to_aa

            node_states[child] = child_states
            node_is_lc[child] = node_is_lc[parent] or is_lc_branch

            if child in node_to_children:
                queue.append(child)

    return node_states, node_is_lc


def get_descendant_leaves(node: str, asm_gene: dict) -> set[str]:
    """Return all leaf-node species descending from `node` (inclusive if leaf)."""
    leaf_nodes = set(asm_gene["leaf_nodes"])
    if node in leaf_nodes:
        return {node}

    node_to_children = asm_gene["node_to_children"]
    leaves: set[str] = set()
    from collections import deque
    queue = deque([node])
    while queue:
        n = queue.popleft()
        for child in node_to_children.get(n, []):
            if child in leaf_nodes:
                leaves.add(child)
            else:
                queue.append(child)
    return leaves


# ─── Branch-length utilities ───────────────────────────────────────────────────

def compute_expected_contact_rate(
    asm_gene: dict,
    contact_pos: str,
    contact_alt_aa: str,
) -> float:
    """
    Compute the expected rate that a DAR gain event falls on a background where
    the contact alt-AA is already present at that position.

    Uses the fraction of all species carrying the contact alt-AA at this position
    (from leaf_states in the ASR data), weighted uniformly across species.

    This is equivalent to assuming DAR gain events are uniformly distributed
    across lineages — an approximation that avoids requiring named internal-node
    branch lengths (which are absent from the available _tree.nwk files).

    A branch-length-weighted version would require IQ-Tree .treefile outputs
    with labelled internal nodes; that file is not currently in the repository.
    """
    leaf_states = asm_gene.get("leaf_states", {})
    if not leaf_states:
        return np.nan

    n_total = 0
    n_with_alt = 0
    for sp_states in leaf_states.values():
        aa = sp_states.get(contact_pos)
        if aa is None:
            continue
        n_total += 1
        if aa == contact_alt_aa:
            n_with_alt += 1

    if n_total == 0:
        return np.nan
    return n_with_alt / n_total


# ─── OR helper ────────────────────────────────────────────────────────────────

def compute_or(n_cdav_alt: int, n_cdav: int, n_bg_alt: int, n_bg: int,
               pseudocount: float = 0.5) -> float:
    """Odds ratio with pseudocount for zero cells."""
    a = n_cdav_alt + pseudocount
    b = (n_cdav - n_cdav_alt) + pseudocount
    c = n_bg_alt + pseudocount
    d = (n_bg - n_bg_alt) + pseudocount
    if c * b == 0:
        return np.inf
    return (a * d) / (b * c)


# ─── Analysis 1: Fisher-only sensitivity ──────────────────────────────────────

def run_fisher_sensitivity(timing_rows: list[dict], all_pairs: dict) -> list[dict]:
    """
    Re-test contact_first enrichment in derived cDAVs using ONLY fisher_fdr
    (not branch_cooccur_fdr) as the significance criterion.

    Reports OR with 95% CI for:
      - All derived cDAVs
      - contact_first = True vs False strata
      - contact_type subsets (mt-mt, nuc-nuc, mt-nuc)

    Still partially confounded; presented as sensitivity check only.
    """
    print("\n── Analysis 1: Fisher-only sensitivity ──────────────────────────────")

    derived_rows = [r for r in timing_rows if r["is_ancestral_cdav"] == "False"]
    print(f"  Derived cDAV pairs to test: {len(derived_rows)}")

    # Join timing → all_tested_pairs for fisher_fdr and contact_type
    results = []
    n_no_match = 0
    for r in derived_rows:
        k = row_key(r)
        ap = all_pairs.get(k)
        if ap is None:
            n_no_match += 1
            continue

        is_significant = float(ap.get("fisher_fdr", "1") or "1") <= 0.10
        is_cf = int(r["n_contact_first"]) > 0

        results.append({
            **{c: r[c] for c in JOIN_COLS},
            "is_ancestral_cdav": r["is_ancestral_cdav"],
            "n_contact_first":   r["n_contact_first"],
            "n_dar_gain_branches": r["n_dar_gain_branches"],
            "dominant_timing":   r["dominant_timing"],
            "contact_type":      ap.get("contact_type", ""),
            "fisher_fdr":        ap.get("fisher_fdr", ""),
            "fisher_sig":        is_significant,
            "contact_first":     is_cf,
            "is_significant_fisher_only": is_significant,
        })

    if n_no_match:
        print(f"  Warning: {n_no_match} timing rows had no match in all_tested_pairs")

    # Aggregate Fisher test: contact_first vs significance (2×2)
    a = sum(1 for r in results if r["contact_first"] and r["fisher_sig"])
    b = sum(1 for r in results if r["contact_first"] and not r["fisher_sig"])
    c = sum(1 for r in results if not r["contact_first"] and r["fisher_sig"])
    d = sum(1 for r in results if not r["contact_first"] and not r["fisher_sig"])

    if (a + b) > 0 and (c + d) > 0:
        oddsratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
        print(f"  contact_first enrichment (fisher-only significant pairs):")
        print(f"    contact_first+sig: {a}  contact_first+non-sig: {b}")
        print(f"    no_cf+sig:         {c}  no_cf+non-sig:         {d}")
        print(f"    OR = {oddsratio:.2f}  p = {pvalue:.2e}")
        print(f"  NOTE: This OR is STILL partially confounded by inheritance.")

    # By contact_type
    for ct in ("mt-mt", "nuc-nuc", "mt-nuc"):
        sub = [r for r in results if r["contact_type"] == ct]
        a2 = sum(1 for r in sub if r["contact_first"] and r["fisher_sig"])
        b2 = sum(1 for r in sub if r["contact_first"] and not r["fisher_sig"])
        c2 = sum(1 for r in sub if not r["contact_first"] and r["fisher_sig"])
        d2 = sum(1 for r in sub if not r["contact_first"] and not r["fisher_sig"])
        if (a2 + b2) > 0 and (c2 + d2) > 0:
            or2, p2 = fisher_exact([[a2, b2], [c2, d2]], alternative="greater")
            print(f"  {ct}: OR={or2:.2f} p={p2:.2e} (n={len(sub)})")

    return results


# ─── Analysis 2: Multi-origin binomial test ────────────────────────────────────

def run_multi_origin_binomial(
    timing_rows: list[dict],
    all_pairs: dict,
    asm_data: dict,
) -> list[dict]:
    """
    For DARs with n_dar_gain_branches >= 2: binomial test comparing the
    fraction of gain events on contact+ parent branches to the
    branch-length-weighted expected rate.

    Groups timing rows by physical DAR+contact identity (not ann_id),
    since the same physical contact pair can have multiple clinical variants.
    """
    print("\n── Analysis 2: Multi-origin binomial test ──────────────────────────")

    # Group by (dar_gene, dar_aa_coord, dar_alt_aa, contact_gene,
    #            contact_refseq_pos, contact_alt_aa)
    phys_key_cols = [
        "dar_gene", "dar_aa_coord", "dar_alt_aa",
        "contact_gene", "contact_refseq_pos", "contact_alt_aa",
    ]

    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in timing_rows:
        k = tuple(r[c] for c in phys_key_cols)
        groups[k].append(r)

    # Filter to multi-origin DARs
    multi = {k: rows for k, rows in groups.items()
             if max(int(r["n_dar_gain_branches"]) for r in rows) >= 2}
    print(f"  DAR+contact pairs with n_dar_gain_branches >= 2: {len(multi)}")

    results = []
    n_skip = 0

    for phys_k, rows in multi.items():
        dar_gene, dar_aa_coord, dar_alt_aa, contact_gene, contact_pos, contact_alt = phys_k

        # Use the row with the most gain branches as representative
        rep = max(rows, key=lambda r: int(r["n_dar_gain_branches"]))
        n_gains = int(rep["n_dar_gain_branches"])
        n_cf    = int(rep["n_contact_first"])

        if n_gains == 0:
            n_skip += 1
            continue

        asm_gene = asm_data.get(contact_gene)
        if asm_gene is None:
            n_skip += 1
            continue

        # Compute expected rate from species-fraction of contact alt-AA
        try:
            expected_rate = compute_expected_contact_rate(
                asm_gene, str(contact_pos), str(contact_alt)
            )
        except Exception as e:
            warnings.warn(f"Expected rate failed for {contact_gene}:{contact_pos}: {e}")
            n_skip += 1
            continue

        if np.isnan(expected_rate) or expected_rate <= 0:
            n_skip += 1
            continue

        # Binomial test: n_contact_first successes out of n_dar_gain_branches
        try:
            btest = binomtest(
                k=n_cf,
                n=n_gains,
                p=expected_rate,
                alternative="greater",
            )
            binom_p = btest.pvalue
        except Exception:
            binom_p = np.nan

        observed_rate = n_cf / n_gains if n_gains > 0 else np.nan

        # Contact type from all_tested_pairs
        ap = all_pairs.get(row_key(rep))
        contact_type = ap.get("contact_type", "") if ap else ""

        results.append({
            **dict(zip(phys_key_cols, phys_k)),
            "n_dar_gain_branches": n_gains,
            "n_contact_first":     n_cf,
            "observed_rate":       round(observed_rate, 4),
            "expected_rate":       round(expected_rate, 4),
            "binom_p":             binom_p,
            "contact_type":        contact_type,
        })

    if n_skip:
        print(f"  Skipped (missing tree / ASR / zero gains): {n_skip}")
    print(f"  Tested: {len(results)} physical DAR+contact pairs")

    # Summary by contact_type
    for ct in ("mt-mt", "nuc-nuc", "mt-nuc", ""):
        label = ct if ct else "all"
        sub = [r for r in results if ct == "" or r["contact_type"] == ct]
        sig = sum(1 for r in sub if float(r["binom_p"] or "1") < 0.05)
        if sub:
            print(f"  {label}: n={len(sub)}, binom_p < 0.05: {sig} ({100*sig/len(sub):.1f}%)")

    return results


# ─── Analysis 3: Tree-aware conditional permutation (pyvolve) ────────────────

def get_iqtree_model(gene: str) -> str:
    """
    Return the Pyvolve-compatible model name for a gene by parsing its
    IQTree output file.  Strips +I / +G4 / +F rate-variation suffixes.
    Falls back to LG for models not supported by pyvolve 1.1.0 (Q.BIRD,
    Q.MAMMAL, Q.PLANT, FLAVI, etc.).
    """
    iqtree_file = ASR_DIR / gene / f"{gene}.iqtree"
    if iqtree_file.exists():
        with open(iqtree_file) as fh:
            for line in fh:
                if line.startswith("Best-fit model according to BIC:"):
                    raw = line.split(":", 1)[1].strip()
                    base = raw.split("+")[0].strip()
                    return _IQTREE_TO_PYVOLVE.get(base, _PYVOLVE_FALLBACK)
    return _PYVOLVE_FALLBACK


def compute_aa_freqs(col_data: list[str]) -> list[float]:
    """
    Compute empirical amino acid frequencies from an alignment column.
    Returns a 20-element list in PYVOLVE_AA_ORDER, summing to 1.
    Gap and ambiguous characters are excluded.

    A Laplace pseudocount of 1 per state is added before normalisation.
    This prevents degenerate frequency vectors (e.g. [1, 0, 0, ...]) for
    perfectly conserved columns, which cause a divide-by-zero in pyvolve's
    Q-matrix scaling step (matrix_builder.py:124).
    """
    counts = Counter(a.upper() for a in col_data if a.upper() in PYVOLVE_AA_ORDER)
    # Laplace smoothing: add 1 to each of the 20 states
    smoothed = [counts.get(aa, 0) + 1 for aa in PYVOLVE_AA_ORDER]
    total = sum(smoothed)
    return [v / total for v in smoothed]


def run_pyvolve_permutation(
    timing_rows: list[dict],
    all_pairs: dict,
    asm_data: dict,
    n_perm: int,
    chunk_idx: int = 0,
    n_chunks: int = 1,
) -> list[dict]:
    """
    Tree-aware conditional permutation test for derived cDAVs with contact_first.

    Simulations are batched: all pairs sharing the same (contact_gene,
    contact_refseq_pos) use a single set of simulated alignments, since the
    underlying model (tree + column frequencies) is identical.

    chunk_idx / n_chunks: for SLURM array parallelization. Each task processes
    a non-overlapping slice of the unique (contact_gene, position) batches.
    Results are written to conditional_permissiveness_chunk{chunk_idx:04d}.csv
    and merged by 05_merge_perm_chunks.py after all tasks complete.
    """
    import pyvolve  # will raise ImportError if not installed

    print(f"\n── Analysis 3: Pyvolve conditional permutation "
          f"(n_perm={n_perm}, chunk {chunk_idx+1}/{n_chunks}) ──")

    # Filter to derived cDAVs with contact_first
    target_rows = [
        r for r in timing_rows
        if r["is_ancestral_cdav"] == "False" and int(r["n_contact_first"]) > 0
    ]
    print(f"  Derived cDAVs with contact_first > 0: {len(target_rows)}")

    # Join to all_tested_pairs for n_cdav_with_alt, n_bg_with_alt etc.
    enriched: list[dict] = []
    for r in target_rows:
        k = row_key(r)
        ap = all_pairs.get(k)
        if ap is None:
            continue
        n_cdav_spp   = int(ap.get("n_cdav_spp", 0) or 0)
        n_bg_spp     = int(ap.get("n_bg_spp", 0) or 0)
        n_cdav_alt   = int(ap.get("n_cdav_with_alt", 0) or 0)
        n_bg_alt     = int(ap.get("n_bg_with_alt", 0) or 0)
        if n_cdav_spp == 0 or n_bg_spp == 0:
            continue
        enriched.append({
            **r,
            "n_cdav_spp":   n_cdav_spp,
            "n_bg_spp":     n_bg_spp,
            "n_cdav_alt":   n_cdav_alt,
            "n_bg_alt":     n_bg_alt,
            "contact_type": ap.get("contact_type", ""),
            "fisher_fdr":   ap.get("fisher_fdr", ""),
        })

    print(f"  After join: {len(enriched)} pairs")

    # Batch by (contact_gene, contact_refseq_pos)
    # For each batch, we will run n_perm simulations ONCE and look up
    # results for all pairs in that batch.
    batch_key_fn = lambda r: (r["contact_gene"], r["contact_refseq_pos"])
    batches_all: dict[tuple, list[dict]] = defaultdict(list)
    for r in enriched:
        batches_all[batch_key_fn(r)].append(r)

    # Chunk: take only every n_chunks-th batch starting at chunk_idx
    all_batch_keys = sorted(batches_all.keys())
    chunk_keys = all_batch_keys[chunk_idx::n_chunks]
    batches = {k: batches_all[k] for k in chunk_keys}

    print(f"  Unique (contact_gene, position) batches (total): {len(batches_all)}")
    print(f"  This chunk ({chunk_idx+1}/{n_chunks}): {len(batches)} batches")

    results_out: list[dict] = []
    n_tree_missing = 0
    n_aln_missing  = 0
    n_no_asm       = 0

    for batch_idx, ((contact_gene, contact_pos_str), batch_rows) in enumerate(batches.items()):
        if (batch_idx + 1) % 50 == 0:
            print(f"  Batch {batch_idx+1}/{len(batches)}", flush=True)

        contact_pos = str(contact_pos_str)

        # Tree path
        tree_path = IQTREE_DIR / contact_gene / f"{contact_gene}_tree.nwk"
        if not tree_path.exists():
            n_tree_missing += len(batch_rows)
            for r in batch_rows:
                results_out.append(_perm_result_stub(r, reason="no_tree"))
            continue

        # Alignment for empirical frequencies
        asm_gene = asm_data.get(contact_gene)
        aln_path = None
        for d in (ALN_NUCLEAR, ALN_MTDNA):
            cand = d / f"{contact_gene}_aa_alignment.fasta"
            if cand.exists():
                aln_path = cand
                break

        if aln_path is None:
            n_aln_missing += len(batch_rows)
            for r in batch_rows:
                results_out.append(_perm_result_stub(r, reason="no_alignment"))
            continue

        if asm_gene is None:
            n_no_asm += len(batch_rows)
            for r in batch_rows:
                results_out.append(_perm_result_stub(r, reason="no_asm"))
            continue

        # Load the alignment column corresponding to this human protein position.
        # ASR maps are keyed by human ungapped protein position; IQTree and FASTA
        # alignments are keyed by gapped alignment site.
        try:
            msa = load_fasta(aln_path)
        except Exception:
            for r in batch_rows:
                results_out.append(_perm_result_stub(r, reason="aln_load_error"))
            continue

        aln_site = asm_gene.get("protein_pos_to_alignment_site", {}).get(contact_pos)
        if aln_site is None:
            for r in batch_rows:
                results_out.append(_perm_result_stub(r, reason="coordinate_missing"))
            continue

        # Extract alignment site j (1-based → 0-based)
        col_idx = int(aln_site) - 1
        col_data = []
        species_order = sorted(msa.keys())
        for sp in species_order:
            seq = msa[sp]
            if col_idx < len(seq):
                col_data.append(seq[col_idx])
            else:
                col_data.append("-")

        aa_freqs = compute_aa_freqs(col_data)

        # Reconstruct ASR node states for contact_pos
        try:
            node_states, node_is_lc = reconstruct_node_states(asm_gene, {contact_pos})
        except Exception:
            for r in batch_rows:
                results_out.append(_perm_result_stub(r, reason="asr_error"))
            continue

        # Leaf species available from ASM
        leaf_nodes = set(asm_gene["leaf_nodes"])
        species_in_tree = {sp for sp in msa if sp in leaf_nodes}
        if not species_in_tree:
            for r in batch_rows:
                results_out.append(_perm_result_stub(r, reason="no_species_overlap"))
            continue

        # Build pyvolve model (per-gene IQTree best-fit model + empirical freqs)
        pyvolve_model_name = get_iqtree_model(contact_gene)
        try:
            tree = pyvolve.read_tree(file=str(tree_path))
            model = pyvolve.Model(pyvolve_model_name, {"state_freqs": aa_freqs})
            partition = pyvolve.Partition(models=model, size=1)
        except Exception as e:
            for r in batch_rows:
                results_out.append(_perm_result_stub(r, reason=f"pyvolve_model_error"))
            continue

        # For each row in this batch, we need the cDAV-clade leaf set and
        # the observed OR
        # Precompute per-row quantities
        batch_meta = []
        for r in batch_rows:
            dar_origin = r["dar_origin_node"]
            contact_alt = r["contact_alt_aa"]
            n_cdav_alt  = r["n_cdav_alt"]
            n_cdav_spp  = r["n_cdav_spp"]
            n_bg_alt    = r["n_bg_alt"]
            n_bg_spp    = r["n_bg_spp"]

            cdav_leaves = get_descendant_leaves(dar_origin, asm_gene)
            bg_leaves   = species_in_tree - cdav_leaves

            observed_or = compute_or(n_cdav_alt, n_cdav_spp, n_bg_alt, n_bg_spp)

            # ASR confidence: is parent-node state (contact alt-AA) assigned
            # from a high-confidence branch?
            parent_of_origin = _find_parent(dar_origin, asm_gene)
            asr_conf = "high"
            if parent_of_origin is not None:
                branch_key = f"{parent_of_origin}|{dar_origin}"
                if branch_key in asm_gene.get("branches_lc", {}):
                    asr_conf = "low"
            else:
                asr_conf = "root"  # at root; effectively high

            batch_meta.append({
                "row":           r,
                "cdav_leaves":   cdav_leaves,
                "bg_leaves":     bg_leaves,
                "contact_alt":   contact_alt,
                "observed_or":   observed_or,
                "asr_conf":      asr_conf,
            })

        # Run simulations (shared across all rows in this batch)
        null_counts: list[dict[str, Counter]] = []
        # null_counts[sim_idx] = {species: simulated_aa}
        for _ in range(n_perm):
            try:
                evolver = pyvolve.Evolver(tree=tree, partitions=partition)
                evolver()
                sim_seqs = evolver.get_sequences()
                # sim_seqs: {species: seq_string} — one character for size=1
                null_counts.append({sp: seq[0] for sp, seq in sim_seqs.items()})
            except Exception:
                null_counts.append({})

        # Compute per-row empirical p-value
        for meta in batch_meta:
            r         = meta["row"]
            cdav_l    = meta["cdav_leaves"]
            bg_l      = meta["bg_leaves"]
            c_alt     = meta["contact_alt"]
            obs_or    = meta["observed_or"]
            asr_conf  = meta["asr_conf"]

            sim_ors = []
            for sim_states in null_counts:
                if not sim_states:
                    continue
                sim_cdav_alt = sum(1 for sp in cdav_l if sim_states.get(sp) == c_alt)
                sim_bg_alt   = sum(1 for sp in bg_l   if sim_states.get(sp) == c_alt)
                sim_or = compute_or(
                    sim_cdav_alt, len(cdav_l),
                    sim_bg_alt,   len(bg_l),
                )
                sim_ors.append(sim_or)

            if sim_ors:
                perm_p        = float(np.mean(np.array(sim_ors) >= obs_or))
                null_median   = float(np.median(sim_ors))
                null_95th     = float(np.percentile(sim_ors, 95))
            else:
                perm_p = null_median = null_95th = np.nan

            ap = all_pairs.get(row_key(r))
            results_out.append({
                **{c: r[c] for c in JOIN_COLS},
                "is_ancestral_cdav":    r["is_ancestral_cdav"],
                "contact_type":         r["contact_type"],
                "n_dar_gain_branches":  r["n_dar_gain_branches"],
                "n_contact_first":      r["n_contact_first"],
                "dar_origin_node":      r["dar_origin_node"],
                "n_cdav_spp":           r["n_cdav_spp"],
                "n_bg_spp":             r["n_bg_spp"],
                "n_cdav_alt":           r["n_cdav_alt"],
                "n_bg_alt":             r["n_bg_alt"],
                "observed_or":          round(obs_or, 3),
                "null_median_or":       round(null_median, 3) if not np.isnan(null_median) else "",
                "null_95th_or":         round(null_95th, 3) if not np.isnan(null_95th) else "",
                "perm_p":               round(perm_p, 4) if not np.isnan(perm_p) else "",
                "n_perm_completed":     len(sim_ors),
                "asr_confidence":       asr_conf,
                "fisher_fdr":           r.get("fisher_fdr", ""),
                "perm_method":          f"pyvolve_{pyvolve_model_name}",
            })

    if n_tree_missing:
        print(f"  Missing gene trees: {n_tree_missing} pairs skipped")
    if n_aln_missing:
        print(f"  Missing alignments: {n_aln_missing} pairs skipped")
    if n_no_asm:
        print(f"  Missing ASR data:   {n_no_asm} pairs skipped")

    return results_out


def _perm_result_stub(r: dict, reason: str) -> dict:
    ap_key = row_key(r)
    return {
        **{c: r.get(c, "") for c in JOIN_COLS},
        "is_ancestral_cdav":   r.get("is_ancestral_cdav", ""),
        "contact_type":        r.get("contact_type", ""),
        "n_dar_gain_branches": r.get("n_dar_gain_branches", ""),
        "n_contact_first":     r.get("n_contact_first", ""),
        "dar_origin_node":     r.get("dar_origin_node", ""),
        "n_cdav_spp":          "",
        "n_bg_spp":            "",
        "n_cdav_alt":          "",
        "n_bg_alt":            "",
        "observed_or":         "",
        "null_median_or":      "",
        "null_95th_or":        "",
        "perm_p":              "",
        "n_perm_completed":    0,
        "asr_confidence":      "",
        "fisher_fdr":          r.get("fisher_fdr", ""),
        "perm_method":         f"skipped:{reason}",
    }


def _find_parent(node: str, asm_gene: dict) -> str | None:
    """Return the parent node of `node` in the ASM tree, or None if root."""
    for parent, children in asm_gene["node_to_children"].items():
        if node in children:
            return parent
    return None


# ─── Pyvolve summary ──────────────────────────────────────────────────────────

def summarise_permutation(results: list[dict]) -> None:
    """Print aggregate summary against pre-registered decision thresholds."""
    valid = [r for r in results if r["perm_p"] not in ("", None)]
    if not valid:
        print("  No permutation results to summarise.")
        return

    perm_ps = [float(r["perm_p"]) for r in valid]
    obs_ors = [float(r["observed_or"]) for r in valid if r["observed_or"] not in ("", None)]

    print(f"\n  Permutation summary ({len(valid)} pairs with perm_p):")
    print(f"    Median observed OR:  {np.median(obs_ors):.2f}")
    print(f"    Fraction perm_p < 0.05: {np.mean(np.array(perm_ps) < 0.05):.2f}")
    print(f"    Fraction perm_p < 0.01: {np.mean(np.array(perm_ps) < 0.01):.2f}")

    # Pre-registered thresholds (from plan)
    median_or = np.median(obs_ors)
    frac_sig  = np.mean(np.array(perm_ps) < 0.05)
    print("\n  Pre-registered decision thresholds:")
    if median_or > 5 and frac_sig > 0.5:
        verdict = "KEEP headline claim: pre-adaptation signal genuine"
    elif 2 < median_or <= 5:
        verdict = "QUALIFY abstract: real but weaker signal (OR 2–5)"
    else:
        verdict = "RETRACT pre-adaptation claim from abstract (OR ≤ 2 or majority non-sig)"
    print(f"    Median OR = {median_or:.2f}, frac_p<0.05 = {frac_sig:.2f}")
    print(f"    → {verdict}")

    # By ASR confidence
    for asr in ("high", "low", "root"):
        sub = [r for r in valid if r.get("asr_confidence") == asr]
        if not sub:
            continue
        sub_ps = [float(r["perm_p"]) for r in sub]
        sub_ors = [float(r["observed_or"]) for r in sub if r["observed_or"]]
        print(f"\n  ASR confidence = {asr} (n={len(sub)}):")
        print(f"    Median OR = {np.median(sub_ors):.2f}")
        print(f"    perm_p < 0.05: {np.mean(np.array(sub_ps) < 0.05):.2f}")


# ─── SDH exclusion filter ─────────────────────────────────────────────────────

SDH_GENES = {"SDHA", "SDHB", "SDHC", "SDHD", "SDHAF2"}


def exclude_sdh(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Split rows into (non_SDH, SDH) based on dar_gene."""
    non_sdh = [r for r in rows if r.get("dar_gene", "") not in SDH_GENES]
    sdh     = [r for r in rows if r.get("dar_gene", "") in SDH_GENES]
    return non_sdh, sdh


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n-perm", type=int, default=1000,
                        help="Number of pyvolve simulations per unique contact position (default: 1000)")
    parser.add_argument("--skip-pyvolve", action="store_true",
                        help="Skip tree-aware permutation (run only analyses 1 and 2)")
    parser.add_argument("--chunk-idx", type=int, default=0,
                        help="0-based chunk index for SLURM array parallelization (default: 0)")
    parser.add_argument("--n-chunks", type=int, default=1,
                        help="Total number of chunks (SLURM array size; default: 1 = no chunking)")
    parser.add_argument("--work-dir", type=Path, default=None,
                        help="Repository data/results root on HPC; code may live elsewhere")
    args = parser.parse_args()

    if args.work_dir is not None:
        global ROOT, TIMING_CSV, ALL_PAIRS, ASM_JSON, IQTREE_DIR, ASR_DIR, ALN_MTDNA, ALN_NUCLEAR, OUT_DIR
        ROOT = args.work_dir
        TIMING_CSV = ROOT / "results" / "phylo" / "timing_annotations.csv"
        ALL_PAIRS = ROOT / "results" / "structural" / "all_tested_pairs.csv"
        ASM_JSON = ROOT / "data" / "phylo" / "ancestral_state_maps.json"
        IQTREE_DIR = ROOT / "data" / "phylo" / "iqtree_jobs"
        ASR_DIR = ROOT / "data" / "phylo" / "ancestral_states"
        ALN_MTDNA = ROOT / "data" / "alignments" / "mtdna_aa"
        ALN_NUCLEAR = ROOT / "data" / "alignments" / "toga_hg38_aa"
        OUT_DIR = ROOT / "results" / "phylo"

    # ── Verify inputs ─────────────────────────────────────────────────────────
    for p in (TIMING_CSV, ALL_PAIRS, ASM_JSON):
        if not p.exists():
            sys.exit(f"ERROR: Required file not found: {p}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────────
    print("Loading data...")
    timing_rows = load_csv_as_dicts(TIMING_CSV)
    all_pairs_list = load_csv_as_dicts(ALL_PAIRS)
    all_pairs = {row_key(r): r for r in all_pairs_list}

    with open(ASM_JSON) as f:
        asm_data = json.load(f)

    print(f"  timing_annotations:  {len(timing_rows)} rows")
    print(f"  all_tested_pairs:    {len(all_pairs_list)} rows")
    print(f"  ancestral_state_maps: {len(asm_data)} genes")

    # ── Analysis 1: Fisher-only sensitivity (always runs) ─────────────────────
    fisher_results = run_fisher_sensitivity(timing_rows, all_pairs)

    fisher_out = OUT_DIR / "contact_first_revised_test.csv"
    if fisher_results:
        fieldnames = list(fisher_results[0].keys())
        with open(fisher_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(fisher_results)
        print(f"\n  Saved: {fisher_out} ({len(fisher_results)} rows)")

    # ── Analysis 2: Multi-origin binomial test (always runs) ─────────────────
    binomial_results = run_multi_origin_binomial(timing_rows, all_pairs, asm_data)

    binom_out = OUT_DIR / "multi_origin_binomial.csv"
    if binomial_results:
        fieldnames = list(binomial_results[0].keys())
        with open(binom_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(binomial_results)
        print(f"\n  Saved: {binom_out} ({len(binomial_results)} rows)")

    # ── Analysis 3: Pyvolve permutation ───────────────────────────────────────
    if args.skip_pyvolve:
        print("\nSkipping pyvolve permutation (--skip-pyvolve specified).")
        return

    try:
        import pyvolve as _pv_check  # noqa: F401
    except ImportError:
        print(
            "\nWARNING: pyvolve not installed — skipping tree-aware permutation.\n"
            "  Install with: pip install pyvolve\n"
            "  Then re-run without --skip-pyvolve.\n"
            "  Analyses 1 and 2 have been saved."
        )
        return

    perm_results = run_pyvolve_permutation(
        timing_rows, all_pairs, asm_data, n_perm=args.n_perm,
        chunk_idx=args.chunk_idx, n_chunks=args.n_chunks,
    )

    if args.n_chunks == 1:
        # Single-run mode: write final merged file and print summary
        summarise_permutation(perm_results)
        perm_out = OUT_DIR / "conditional_permissiveness.csv"
    else:
        # Chunked mode: write per-chunk file; merge with 05_merge_perm_chunks.py later
        perm_out = OUT_DIR / f"conditional_permissiveness_chunk{args.chunk_idx:04d}.csv"
        print(f"  Chunk mode: writing to {perm_out.name}")

    if perm_results:
        fieldnames = list(perm_results[0].keys())
        with open(perm_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(perm_results)
        print(f"\n  Saved: {perm_out} ({len(perm_results)} rows)")

    print("\nDone.")


if __name__ == "__main__":
    main()
