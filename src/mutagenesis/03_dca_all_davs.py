#!/usr/bin/env python3
"""
src/mutagenesis/03_dca_all_davs.py

Compute plmDCA direct coupling scores for ALL DAV contact pairs (cDAVs + uDAVs)
and compare DI score distributions between the two groups.

Scientific question
-------------------
Do contact positions of compensated disease variants (cDAVs) show higher
co-evolutionary coupling (plmDCA DI) than contact positions of uncompensated
DAVs (uDAVs)? Higher DI indicates the DAR and contact positions co-evolve
across mammals, consistent with structural compensation being selectively
maintained rather than accidental.

DI is a co-evolutionary correlation, not proof of causation. It complements
the phylogenetic evidence (Pagel, branch co-occurrence) from Stage 4.

Comparison design
-----------------
Primary:   within-gene Mann-Whitney U (controls for gene-identity confounding)
Secondary: global cDAV vs uDAV (subject to gene-identity confounding; labelled)

Inputs:
  results/structural/dar_contacts_cbcb8A.csv  -- all mapped DAV contacts
  data/alignments/toga_hg38_aa/              -- nucDNA AA alignments
  data/alignments/mtdna_aa/                  -- mtDNA AA alignments

Outputs:
  results/mutagenesis/dca_all_davs.csv
  results/mutagenesis/dca_cdav_vs_udav_comparison.csv

Usage (from project root, inside Docker container):
  conda run -n oxphos_dav python src/mutagenesis/03_dca_all_davs.py \\
      --aggregate results/mutagenesis/dca_gene_results

  --aggregate DIR   Required. Per-gene .npz files produced by 03b_dca_gene_worker
                    on HPC. plmDCA is always run on HPC; this script aggregates.
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import entropy as scipy_entropy, mannwhitneyu

ROOT = Path(__file__).resolve().parents[2]

CONTACTS_CSV = ROOT / "results" / "structural" / "dar_contacts_cbcb8A.csv"
ALN_NUC_AA   = ROOT / "data" / "alignments" / "toga_hg38_aa"
ALN_MT_AA    = ROOT / "data" / "alignments" / "mtdna_aa"
OUT_DIR      = ROOT / "results" / "mutagenesis"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# COXFA4 ↔ NDUFA4 alias bridge (contacts CSV may use either name)
_ALIAS_EQUIV = {"COXFA4": "NDUFA4", "NDUFA4": "COXFA4"}
_MASK        = {"-", "X", "!", "*"}

ALPHABET  = list("ACDEFGHIKLMNPQRSTVWY-")
AA_TO_IDX = {aa: i for i, aa in enumerate(ALPHABET)}
N_STATES  = len(ALPHABET)


# ── Alignment loading ──────────────────────────────────────────────────────────

_aln_cache: dict[tuple[str, str], dict[str, str] | None] = {}


def _aln_dir(genome: str) -> Path:
    return ALN_MT_AA if genome == "mtDNA" else ALN_NUC_AA


def load_aln(gene: str, genome: str) -> dict[str, str] | None:
    """Load {rec_id: gapped_seq} from FASTA. Caches per (gene, genome).
    Tries alias if primary name not found."""
    key = (gene, genome)
    if key in _aln_cache:
        return _aln_cache[key]
    path = _aln_dir(genome) / f"{gene}_aa_alignment.fasta"
    if not path.exists() and gene in _ALIAS_EQUIV:
        path = _aln_dir(genome) / f"{_ALIAS_EQUIV[gene]}_aa_alignment.fasta"
    if not path.exists():
        _aln_cache[key] = None
        return None
    aln = {rec.id: str(rec.seq).upper()
           for rec in SeqIO.parse(str(path), "fasta")}
    _aln_cache[key] = aln or None
    return _aln_cache[key]


# ── Position mapping ───────────────────────────────────────────────────────────

def build_pos_to_col(ref_seq: str) -> dict[int, int]:
    """1-based ungapped biological position → 0-based alignment column."""
    pos_map: dict[int, int] = {}
    pos = 0
    for col, ch in enumerate(ref_seq):
        if ch not in _MASK:
            pos += 1
            pos_map[pos] = col
    return pos_map


def get_ref_seq(aln: dict[str, str]) -> str | None:
    """Return Homo_sapiens sequence; handles 'Homo_sapiens|...' keys."""
    if "Homo_sapiens" in aln:
        return aln["Homo_sapiens"]
    for key, seq in aln.items():
        if key.startswith("Homo_sapiens"):
            return seq
    return None


# ── Column extraction and entropy ─────────────────────────────────────────────

def extract_col(seqs: list[str], c: int) -> list[str]:
    return [s[c] if c < len(s) else "-" for s in seqs]


def col_entropy(col: list[str]) -> float:
    counts: dict[str, int] = {}
    for aa in col:
        if aa not in _MASK:
            counts[aa] = counts.get(aa, 0) + 1
    if not counts:
        return 0.0
    vals = np.array(list(counts.values()), dtype=float)
    return float(scipy_entropy(vals / vals.sum(), base=20))


# ── APC-MI (secondary, phylogenetically confounded) ──────────────────────────

def _joint(ci: list[str], cj: list[str]) -> np.ndarray:
    mat = np.zeros((N_STATES, N_STATES), dtype=float)
    for a, b in zip(ci, cj):
        mat[AA_TO_IDX.get(a, N_STATES - 1), AA_TO_IDX.get(b, N_STATES - 1)] += 1.0
    n = mat.sum()
    if n > 0:
        mat /= n
    return mat


def mutual_information(ci: list[str], cj: list[str]) -> float:
    joint = _joint(ci, cj)
    mi_i  = joint.sum(axis=1)
    mi_j  = joint.sum(axis=0)
    return max(0.0, float(scipy_entropy(mi_i) + scipy_entropy(mi_j) - scipy_entropy(joint.ravel())))


_mi_stats_cache: dict = {}


def _gene_mi_stats(gene: str, genome: str, sample: int = 5000):
    key = (gene, genome)
    if key in _mi_stats_cache:
        return _mi_stats_cache[key]
    aln = load_aln(gene, genome)
    if aln is None:
        _mi_stats_cache[key] = None
        return None
    seqs = list(aln.values())
    L = len(seqs[0]) if seqs else 0
    if L < 2:
        _mi_stats_cache[key] = None
        return None
    rng = np.random.default_rng(42)
    if L <= 100:
        pairs = [(i, j) for i in range(L) for j in range(i + 1, L)]
    else:
        idx = rng.integers(0, L, size=(sample * 2, 2))
        idx = idx[idx[:, 0] != idx[:, 1]][:sample]
        pairs = [(int(a), int(b)) for a, b in idx]
    mi_vals: dict[tuple[int, int], float] = {}
    for i, j in pairs:
        v = mutual_information(extract_col(seqs, i), extract_col(seqs, j))
        mi_vals[(i, j)] = mi_vals[(j, i)] = v
    row_sum: dict[int, float] = {}
    row_cnt: dict[int, int]   = {}
    for (i, j), v in mi_vals.items():
        row_sum[i] = row_sum.get(i, 0.0) + v
        row_cnt[i] = row_cnt.get(i, 0) + 1
    row_means = {i: row_sum[i] / row_cnt[i] for i in row_sum}
    gm = float(np.mean(list(mi_vals.values()))) if mi_vals else 0.0
    result = (row_means, gm)
    _mi_stats_cache[key] = result
    return result


def apc_correct(mi_raw: float, i: int, j: int, stats) -> float:
    if stats is None:
        return mi_raw
    row_means, gm = stats
    if gm == 0:
        return mi_raw
    return float(mi_raw - (row_means.get(i, gm) * row_means.get(j, gm)) / gm)


# ── plmDCA HPC results (populated by --aggregate before pair scoring) ─────────
_hpc_di:   dict[str, dict[tuple[int, int], float]] = {}
_hpc_meff: dict[str, float]     = {}
_hpc_pct:  dict[str, np.ndarray] = {}


def load_hpc_results(npz_dir: Path) -> None:
    """Load per-gene .npz files produced by 03b_dca_gene_worker on HPC."""
    npz_files = sorted(npz_dir.glob("*.npz"))
    if not npz_files:
        sys.exit(f"ERROR: no .npz files found in {npz_dir}")
    print(f"Loading HPC DCA results from {npz_dir}  ({len(npz_files)} files)...")
    for npz_path in npz_files:
        gene = npz_path.stem
        try:
            data = np.load(str(npz_path))
        except Exception as e:
            print(f"  [WARN] {gene}: could not load {npz_path.name}: {e}")
            continue
        col_is    = data["col_i"].astype(int)
        col_js    = data["col_j"].astype(int)
        di_scores = data["di_score"].astype(float)
        meff      = float(data["meff"][0]) if "meff" in data else np.nan
        di_dict: dict[tuple[int, int], float] = {}
        for ci, cj, sc in zip(col_is, col_js, di_scores):
            di_dict[(int(ci), int(cj))] = float(sc)
        _hpc_di[gene]   = di_dict
        _hpc_meff[gene] = meff
        _hpc_pct[gene]  = di_scores
        print(f"  {gene}: {len(di_dict):,} pairs  Meff={meff:.1f}")
    print(f"Loaded {len(_hpc_di)} genes from HPC results.")


def get_gene_dca(gene: str, genome: str
                 ) -> tuple[dict[tuple[int, int], float], float, np.ndarray]:
    """Look up HPC-computed plmDCA results for a gene. Tries COXFA4/NDUFA4 alias."""
    for name in (gene, _ALIAS_EQUIV.get(gene, "")):
        if name and name in _hpc_di:
            return _hpc_di[name], _hpc_meff.get(name, np.nan), _hpc_pct.get(name, np.array([]))
    return {}, np.nan, np.array([])


def dca_percentile(di_val: float, all_di: np.ndarray) -> float | None:
    if all_di is None or len(all_di) == 0 or np.isnan(di_val):
        return None
    return float(np.mean(all_di < di_val) * 100)


def infer_contact_genome(dar_genome: str, contact_type: str) -> str:
    ct = str(contact_type)
    if ct == "mt-mt":
        return "mtDNA"
    if ct == "nuc-nuc":
        return "nucDNA"
    if ct == "mt-nuc":
        return "nucDNA"
    return dar_genome


# ── Per-pair scoring ──────────────────────────────────────────────────────────

def compute_pair(
    dar_gene: str, dar_coord: int,
    contact_gene: str, contact_coord: int,
    dar_genome: str, contact_genome: str,
) -> dict:
    out = dict(
        dca_di=np.nan, dca_di_percentile=np.nan, dca_meff=np.nan,
        dca_note="", dar_entropy=np.nan, contact_entropy=np.nan,
        mi_raw=np.nan, mi_apc=np.nan, inter_protein=False,
    )
    is_intra = (dar_gene == contact_gene)
    out["inter_protein"] = not is_intra

    dar_aln = load_aln(dar_gene, dar_genome)
    if dar_aln is None:
        out["dca_note"] = "no_alignment"
        return out

    contact_aln = dar_aln if is_intra else load_aln(contact_gene, contact_genome)
    if contact_aln is None:
        out["dca_note"] = "no_alignment"
        return out

    dar_ref = get_ref_seq(dar_aln)
    if dar_ref is None:
        out["dca_note"] = "no_human_ref"
        return out

    dar_ptc = build_pos_to_col(dar_ref)
    col_i   = dar_ptc.get(dar_coord)

    if is_intra:
        col_j = dar_ptc.get(contact_coord)
    else:
        contact_ref = get_ref_seq(contact_aln)
        if contact_ref is None:
            out["dca_note"] = "no_human_ref"
            return out
        col_j = build_pos_to_col(contact_ref).get(contact_coord)

    if col_i is None or col_j is None:
        out["dca_note"] = "pos_not_in_aln"
        return out

    dar_seqs     = list(dar_aln.values())
    contact_seqs = list(contact_aln.values())
    ci_col = extract_col(dar_seqs, col_i)
    cj_col = extract_col(contact_seqs, col_j)
    out["dar_entropy"]     = col_entropy(ci_col)
    out["contact_entropy"] = col_entropy(cj_col)

    # APC-MI (secondary)
    mi_raw = mutual_information(ci_col, cj_col)
    if is_intra:
        stats  = _gene_mi_stats(dar_gene, dar_genome)
        mi_apc = apc_correct(mi_raw, col_i, col_j, stats)
    else:
        mi_apc = mi_raw
    out["mi_raw"] = mi_raw
    out["mi_apc"] = mi_apc

    # plmDCA (from HPC results)
    notes = []
    if dar_genome != contact_genome:
        notes.append("mt_nuc_asymmetry")

    if is_intra:
        dca_result = get_gene_dca(dar_gene, dar_genome)
        if dca_result is None:
            out["dca_note"] = "plmdca_unavailable"
            return out
        di_dict, meff, all_di = dca_result
        if not di_dict:  # crashed or OOM in subprocess
            out["dca_note"] = ",".join(notes + ["plmdca_no_result"])
            out["dca_meff"] = meff
            return out
        key    = (min(col_i, col_j), max(col_i, col_j))
        di_val = di_dict.get(key, np.nan)
        if np.isnan(di_val):
            notes.append("pair_not_in_di_dict")
        pct = dca_percentile(di_val, all_di)
        out["dca_di"]            = di_val
        out["dca_di_percentile"] = pct if pct is not None else np.nan
        out["dca_meff"]          = meff
        if not np.isnan(meff) and meff < len(dar_aln):
            notes.append("insufficient_meff")
    else:
        # Interprotein: per-gene HPC .npz files are intraprotein only.
        notes.append("hpc_intra_only")
        out["dca_note"] = ",".join(notes)
        return out

    out["dca_note"] = ",".join(notes)
    return out


# ── cDAV vs uDAV comparison ────────────────────────────────────────────────────

def run_comparison(df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def _add(label: str, mask: pd.Series):
        sub   = df[mask & df["dca_di_percentile"].notna()]
        sub_c = sub[sub["is_cdav_amino_acid"]]["dca_di_percentile"]
        sub_u = sub[~sub["is_cdav_amino_acid"]]["dca_di_percentile"]
        if len(sub_c) < 3 or len(sub_u) < 3:
            return
        stat, p = mannwhitneyu(sub_c, sub_u, alternative="greater")
        rows.append({
            "stratum":        label,
            "n_cdav":         len(sub_c),
            "n_udav":         len(sub_u),
            "median_cdav":    float(sub_c.median()),
            "median_udav":    float(sub_u.median()),
            "p75_cdav":       float(sub_c.quantile(0.75)),
            "p75_udav":       float(sub_u.quantile(0.75)),
            "mannwhitney_u":  float(stat),
            "p_cdav_gt_udav": float(p),
        })

    intra = ~df["inter_protein"]
    _add("all_pairs",         pd.Series(True, index=df.index))
    _add("intraprotein_only", intra)
    _add("interprotein_only",  df["inter_protein"])

    # Primary: within-gene intraprotein (controls for gene identity)
    for gene in sorted(df["dar_gene"].unique()):
        _add(f"within_gene_{gene}", (df["dar_gene"] == gene) & intra)

    # By contact class
    for cls in ("hbond", "electrostatic", "hydrophobic", "vdw"):
        _add(f"contact_class_{cls}", df["contact_classes"].str.contains(cls, na=False) & intra)

    return pd.DataFrame(rows)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--aggregate", metavar="DIR", required=True,
                        help="Directory of per-gene .npz files from 03b_dca_gene_worker "
                             "(produced on HPC). plmDCA is always run on HPC.")
    args = parser.parse_args()

    if not CONTACTS_CSV.exists():
        sys.exit(f"ERROR: {CONTACTS_CSV} not found. Run structural analysis first.")

    load_hpc_results(Path(args.aggregate))

    # ── Load and deduplicate contacts ──────────────────────────────────────────
    print(f"Loading {CONTACTS_CSV.name}...")
    contacts = pd.read_csv(CONTACTS_CSV, dtype={"ann_id": str, "variant_id": str},
                           low_memory=False)
    print(f"  {len(contacts):,} rows")

    gene_col = "dar_locus" if "dar_locus" in contacts.columns else "dar_gene"

    agg = (
        contacts
        .assign(is_cdav=contacts["is_cdav_amino_acid"].fillna(False).astype(bool))
        .groupby([gene_col, "dar_aa_coord", "contact_gene", "contact_refseq_pos"],
                 dropna=False)
        .agg(
            is_cdav_amino_acid=("is_cdav",         "any"),
            dar_genome         =("dar_genome",       "first"),
            contact_type       =("contact_type",     "first"),
            contact_classes    =("contact_class",    lambda x: ",".join(sorted(set(x.dropna())))),
            n_source_rows      =("ann_id",           "count"),
        )
        .reset_index()
        .rename(columns={gene_col: "dar_gene"})
    )

    print(f"  {len(agg):,} unique (dar_gene, dar_aa_coord, contact_gene, contact_refseq_pos)")
    n_cdav = int(agg["is_cdav_amino_acid"].sum())
    print(f"  cDAV pairs: {n_cdav:,}   uDAV pairs: {len(agg) - n_cdav:,}")

    # ── Score each pair ────────────────────────────────────────────────────────
    results = []
    n_total = len(agg)
    for i, row in enumerate(agg.itertuples(index=False), 1):
        if i == 1 or i % 500 == 0:
            print(f"  [{i}/{n_total}] {row.dar_gene}:{row.dar_aa_coord} — "
                  f"{row.contact_gene}:{row.contact_refseq_pos}", flush=True)

        dar_genome     = str(row.dar_genome) if pd.notna(row.dar_genome) else "nucDNA"
        contact_genome = infer_contact_genome(dar_genome, str(row.contact_type))

        try:
            dar_coord     = int(float(row.dar_aa_coord))
            contact_coord = int(float(row.contact_refseq_pos))
        except (ValueError, TypeError):
            scores = dict(dca_di=np.nan, dca_di_percentile=np.nan, dca_meff=np.nan,
                          dca_note="bad_coord", dar_entropy=np.nan, contact_entropy=np.nan,
                          mi_raw=np.nan, mi_apc=np.nan, inter_protein=False)
        else:
            scores = compute_pair(
                dar_gene=row.dar_gene,     dar_coord=dar_coord,
                contact_gene=row.contact_gene, contact_coord=contact_coord,
                dar_genome=dar_genome,     contact_genome=contact_genome,
            )

        results.append({
            "dar_gene":           row.dar_gene,
            "dar_aa_coord":       row.dar_aa_coord,
            "contact_gene":       row.contact_gene,
            "contact_refseq_pos": row.contact_refseq_pos,
            "is_cdav_amino_acid": row.is_cdav_amino_acid,
            "contact_classes":    row.contact_classes,
            "n_source_rows":      row.n_source_rows,
            "dar_genome":         dar_genome,
            "contact_genome":     contact_genome,
            **scores,
        })

    df = pd.DataFrame(results)

    out_csv = OUT_DIR / "dca_all_davs.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {len(df):,} rows → {out_csv}")

    # ── Diagnostics ───────────────────────────────────────────────────────────
    note_counts = df["dca_note"].value_counts(dropna=False)
    print("\nDCA note breakdown:")
    for note, cnt in note_counts.items():
        print(f"  {str(note):35s}: {cnt:6,} ({100*cnt/len(df):.1f}%)")

    n_scored = int(df["dca_di"].notna().sum())
    print(f"\nPairs with DI score: {n_scored:,}/{len(df):,} ({100*n_scored/len(df):.1f}%)")

    pos_miss = int((df["dca_note"] == "pos_not_in_aln").sum())
    if pos_miss > 0.05 * len(df):
        warnings.warn(
            f"pos_not_in_aln rate is {100*pos_miss/len(df):.1f}% — "
            "check build_pos_to_col is using the correct reference sequence."
        )

    # ── Comparison ────────────────────────────────────────────────────────────
    print("\nRunning cDAV vs uDAV comparison...")
    cmp_df  = run_comparison(df)
    cmp_csv = OUT_DIR / "dca_cdav_vs_udav_comparison.csv"
    cmp_df.to_csv(cmp_csv, index=False)
    print(f"Wrote {len(cmp_df):,} strata → {cmp_csv}")

    print("\n" + "=" * 72)
    print(f"{'Stratum':<28} {'N_cDAV':>7} {'N_uDAV':>7} "
          f"{'med_cDAV':>9} {'med_uDAV':>9} {'p(cDAV>uDAV)':>13}")
    print("-" * 72)
    for _, r in cmp_df.iterrows():
        if r["stratum"] in ("all_pairs", "intraprotein_only", "interprotein_only"):
            print(
                f"  {r['stratum']:<26} {int(r['n_cdav']):>7} {int(r['n_udav']):>7} "
                f"{r['median_cdav']:>9.1f} {r['median_udav']:>9.1f} "
                f"{r['p_cdav_gt_udav']:>13.4g}"
            )
    print("=" * 72)
    print("\nWithin-gene strata (primary, gene-identity-controlled) → comparison CSV.")
    print("Intraprotein results are most reliable (no concatenated MSA required).")


if __name__ == "__main__":
    main()
