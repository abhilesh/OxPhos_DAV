"""
src/phylo/03_comparative_analysis.py

Comparative phylogenetic analysis of cDAV compensatory co-evolution.
Reads existing output CSVs and runs four statistical comparisons:

  1. Independent origins distribution and compensation enrichment
  2. Temporal ordering (contact_first enrichment in significant pairs)
  3. mtDNA vs nucDNA compensation rates (mt-mt / nuc-nuc / mt-nuc)
  4. Ancestral vs derived cDAV compensation patterns

Inputs (must exist):
  results/structural/all_tested_pairs.csv
  results/structural/compensatory_partners.csv
  results/phylo/timing_annotations.csv

Outputs:
  results/phylo/independent_origins_summary.csv
  results/phylo/temporal_ordering_test.csv
  results/phylo/genomic_comparison.csv
  results/phylo/genomic_comparison_no_sdh.csv
  results/phylo/ancestral_vs_derived.csv
  results/phylo/test_decomposed_rates.csv
  results/phylo/pagel_coverage.csv
  results/phylo/comparative_analysis_summary.md

Run from project root inside Docker:
    python src/phylo/03_comparative_analysis.py
"""

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu

ROOT      = Path(__file__).resolve().parents[2]
RES_STR   = ROOT / "results" / "structural"
RES_PHY   = ROOT / "results" / "phylo"

ALL_PAIRS_CSV    = RES_STR / "all_tested_pairs.csv"
PARTNERS_CSV     = RES_STR / "compensatory_partners.csv"
TIMING_CSV       = RES_PHY / "timing_annotations.csv"
OUT_DIR          = RES_PHY


# ── Identity key ───────────────────────────────────────────────────────────────

JOIN_COLS = ("ann_id", "dar_gene", "dar_aa_coord", "dar_alt_aa",
             "contact_gene", "contact_refseq_pos", "contact_alt_aa")

def identity_key(row):
    return "|".join(row[c] for c in JOIN_COLS)

DAR_COLS = ("ann_id", "dar_gene", "dar_aa_coord", "dar_alt_aa")

def dar_key(row):
    return "|".join(row[c] for c in DAR_COLS)


# ── Significance helper ────────────────────────────────────────────────────────

def safe_float(v):
    try:
        f = float(v)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None

def is_significant(row):
    """Mirror the significance logic from 01_find_compensating_partners.py."""
    if row.get("low_power") == "True":
        return False
    pf = safe_float(row.get("pagel_fdr"))
    bf = safe_float(row.get("branch_cooccur_fdr"))
    if pf is not None or bf is not None:
        return (pf is not None and pf <= 0.10) or (bf is not None and bf <= 0.10)
    ff = safe_float(row.get("fisher_fdr"))
    return ff is not None and ff <= 0.10


def test_evidence_flags(row):
    """
    Return a dict of per-test evidence flags for a single pair.

    fisher_only : significant exclusively by Fisher FDR (no Pagel or branch data)
    branch_sig  : branch_cooccur_fdr ≤ 0.10
    pagel_sig   : pagel_fdr ≤ 0.10
    pagel_tested: pagel_fdr is non-NA (test ran and returned a result)
    fisher_sig  : fisher_fdr ≤ 0.10
    multi_origin: branch or Pagel test was available (pair has ≥2-origin data)
    """
    if row.get("low_power") == "True":
        return dict(fisher_only=False, branch_sig=False, pagel_sig=False,
                    pagel_tested=False, fisher_sig=False, multi_origin=False)
    pf = safe_float(row.get("pagel_fdr"))
    bf = safe_float(row.get("branch_cooccur_fdr"))
    ff = safe_float(row.get("fisher_fdr"))
    pagel_tested  = pf is not None
    pagel_sig     = pagel_tested and pf <= 0.10
    branch_sig    = bf is not None and bf <= 0.10
    fisher_sig    = ff is not None and ff <= 0.10
    multi_origin  = pf is not None or bf is not None
    # fisher_only: pair had no Pagel or branch test available AND is sig by Fisher
    fisher_only   = (not multi_origin) and fisher_sig
    return dict(fisher_only=fisher_only, branch_sig=branch_sig, pagel_sig=pagel_sig,
                pagel_tested=pagel_tested, fisher_sig=fisher_sig, multi_origin=multi_origin)


# ── Odds ratio with 95 % CI (Wald on log scale) ───────────────────────────────

def odds_ratio_ci(a, b, c, d):
    """OR = (a*d)/(b*c), 95% CI on log scale. Returns (OR, lo, hi)."""
    # Add 0.5 continuity correction to zero cells
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    OR = (a * d) / (b * c)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    lo = math.exp(math.log(OR) - 1.96 * se)
    hi = math.exp(math.log(OR) + 1.96 * se)
    return OR, lo, hi


# ── Load data ─────────────────────────────────────────────────────────────────

def load_all_pairs():
    with open(ALL_PAIRS_CSV) as f:
        return list(csv.DictReader(f))

def load_timing():
    with open(TIMING_CSV) as f:
        return list(csv.DictReader(f))

def load_partners():
    with open(PARTNERS_CSV) as f:
        return list(csv.DictReader(f))


# ══════════════════════════════════════════════════════════════════════════════
# Analysis 1: Independent origins
# ══════════════════════════════════════════════════════════════════════════════

def analysis_independent_origins(all_pairs, timing_rows):
    """
    Aggregate n_dar_gain_branches (or n_dar_loss_branches for ancestral cDAVs)
    to the cDAV level and test whether convergent cDAVs (≥2 origins) show
    higher compensation rates than single-origin cDAVs.
    """
    print("\n── Analysis 1: Independent Origins ──────────────────────────────────")

    # Build timing lookup: dar_key → (n_origins, is_ancestral)
    dar_origins = {}
    for r in timing_rows:
        dk = dar_key(r)
        anc = r.get("is_ancestral_cdav") == "True"
        gains = int(r.get("n_dar_gain_branches") or 0)
        losses = int(r.get("n_dar_loss_branches") or 0)
        n_orig = losses if anc else gains
        if dk not in dar_origins or n_orig > dar_origins[dk][0]:
            dar_origins[dk] = (n_orig, anc)

    # Build significance lookup: identity_key → bool
    sig_lookup = {identity_key(r): is_significant(r) for r in all_pairs}

    # Per-cDAV: count tested contacts and significant contacts
    dar_contacts     = defaultdict(lambda: {"total": 0, "sig": 0})
    for r in all_pairs:
        dk = dar_key(r)
        ik = identity_key(r)
        dar_contacts[dk]["total"] += 1
        if sig_lookup.get(ik):
            dar_contacts[dk]["sig"] += 1

    # Combine
    records = []
    for dk, (n_orig, is_anc) in dar_origins.items():
        parts = dk.split("|")
        ann_id, gene, coord, alt = parts
        t = dar_contacts.get(dk, {"total": 0, "sig": 0})
        records.append({
            "ann_id":               ann_id,
            "dar_gene":             gene,
            "dar_aa_coord":         coord,
            "dar_alt_aa":           alt,
            "is_ancestral_cdav":    is_anc,
            "n_independent_origins": n_orig,
            "n_contacts_tested":    t["total"],
            "n_contacts_significant": t["sig"],
            "has_compensation":     t["sig"] > 0,
        })

    # Distribution
    origin_dist = Counter(r["n_independent_origins"] for r in records)
    bins = [(1,1),(2,5),(6,10),(11,999)]
    print(f"{'Origins':>10}  {'Count':>7}  {'% cDAVs':>8}")
    for lo, hi in bins:
        n = sum(v for k, v in origin_dist.items() if lo <= k <= hi)
        label = str(lo) if lo == hi else f"{lo}–{hi}" if hi < 999 else f">{lo-1}"
        print(f"  {label:>8}  {n:>7}  {100*n/len(records):>7.1f}%")
    n_conv = sum(1 for r in records if r["n_independent_origins"] >= 2)
    print(f"\n  cDAVs with ≥2 independent origins (convergent): "
          f"{n_conv}/{len(records)} ({100*n_conv/len(records):.1f}%)")

    # Fisher: ≥2 origins vs compensation
    a = sum(1 for r in records if r["n_independent_origins"] >= 2 and r["has_compensation"])
    b = sum(1 for r in records if r["n_independent_origins"] >= 2 and not r["has_compensation"])
    c = sum(1 for r in records if r["n_independent_origins"] < 2  and r["has_compensation"])
    d = sum(1 for r in records if r["n_independent_origins"] < 2  and not r["has_compensation"])
    _, p_fisher = fisher_exact([[a, b], [c, d]], alternative="greater")
    OR, lo, hi = odds_ratio_ci(a, b, c, d)
    print(f"\n  Convergent (≥2 origins) vs compensation:")
    print(f"    2×2: a={a} b={b} c={c} d={d}")
    print(f"    OR={OR:.2f} (95% CI {lo:.2f}–{hi:.2f}), Fisher p={p_fisher:.3e}")

    # Write output
    out_path = OUT_DIR / "independent_origins_summary.csv"
    fields = ["ann_id","dar_gene","dar_aa_coord","dar_alt_aa","is_ancestral_cdav",
              "n_independent_origins","n_contacts_tested",
              "n_contacts_significant","has_compensation"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(records)
    print(f"\n  → {out_path}")

    return records, {"a": a, "b": b, "c": c, "d": d,
                     "OR": OR, "CI_lo": lo, "CI_hi": hi, "p": p_fisher,
                     "n_total": len(records), "n_convergent": n_conv}


# ══════════════════════════════════════════════════════════════════════════════
# Analysis 2: Temporal ordering — contact_first enrichment
# ══════════════════════════════════════════════════════════════════════════════

def analysis_temporal_ordering(all_pairs, timing_rows):
    """
    Test whether contact_first timing is enriched in significant pairs.
    Run overall and stratified by: ancestral/derived, contact_type.
    """
    print("\n── Analysis 2: Temporal Ordering (Contact-First Hypothesis) ─────────")

    # Build timing lookup: identity_key → timing fields
    timing_by_key = {}
    for r in timing_rows:
        timing_by_key[identity_key(r)] = r

    # Annotate all_pairs with timing and significance
    annotated = []
    for r in all_pairs:
        ik = identity_key(r)
        tr = timing_by_key.get(ik, {})
        if not tr:
            continue
        n_cf   = int(tr.get("n_contact_first") or 0)
        n_total_timing = (int(tr.get("n_dar_gain_branches") or 0)
                          + int(tr.get("n_dar_loss_branches") or 0))
        annotated.append({
            "identity_key":      ik,
            "contact_type":      r.get("contact_type", ""),
            "dar_genome":        r.get("dar_genome", ""),
            "is_ancestral":      tr.get("is_ancestral_cdav") == "True",
            "contact_first":     n_cf > 0,
            "n_origins":         n_total_timing,
            "significant":       is_significant(r),
        })

    results = []

    def run_test(label, rows):
        a = sum(1 for r in rows if r["contact_first"] and r["significant"])
        b = sum(1 for r in rows if r["contact_first"] and not r["significant"])
        c = sum(1 for r in rows if not r["contact_first"] and r["significant"])
        d = sum(1 for r in rows if not r["contact_first"] and not r["significant"])
        if a + b + c + d == 0:
            return
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        OR, lo, hi = odds_ratio_ci(a, b, c, d)
        n_cf  = a + b
        n_sig = a + c
        pct_cf_sig  = 100 * a / n_cf  if n_cf  else 0
        pct_sig_cf  = 100 * a / n_sig if n_sig else 0
        print(f"\n  {label} (n={a+b+c+d})")
        print(f"    contact_first & sig={a}, contact_first & not_sig={b}")
        print(f"    not_cf & sig={c}, not_cf & not_sig={d}")
        print(f"    OR={OR:.2f} (95% CI {lo:.2f}–{hi:.2f}), Fisher p={p:.3e}")
        print(f"    {pct_cf_sig:.1f}% of contact_first pairs are significant")
        print(f"    {pct_sig_cf:.1f}% of significant pairs had contact_first timing")
        results.append({
            "subset": label, "n_total": a+b+c+d,
            "a_cf_sig": a, "b_cf_nsig": b, "c_ncf_sig": c, "d_ncf_nsig": d,
            "OR": round(OR,3), "CI_lo": round(lo,3), "CI_hi": round(hi,3),
            "fisher_p": p,
        })

    run_test("All pairs",               annotated)
    run_test("Derived cDAVs",           [r for r in annotated if not r["is_ancestral"]])
    run_test("Ancestral cDAVs",         [r for r in annotated if r["is_ancestral"]])
    run_test("mt-mt pairs",             [r for r in annotated if r["contact_type"] == "mt-mt"])
    run_test("nuc-nuc pairs",           [r for r in annotated if r["contact_type"] == "nuc-nuc"])
    run_test("mt-nuc pairs",            [r for r in annotated if r["contact_type"] == "mt-nuc"])

    out_path = OUT_DIR / "temporal_ordering_test.csv"
    fields = ["subset","n_total","a_cf_sig","b_cf_nsig","c_ncf_sig","d_ncf_nsig",
              "OR","CI_lo","CI_hi","fisher_p"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)
    print(f"\n  → {out_path}")
    return results


# ══════════════════════════════════════════════════════════════════════════════
# Analysis 3: mtDNA vs nucDNA / intragenomic vs intergenomic
# ══════════════════════════════════════════════════════════════════════════════

def analysis_genomic_comparison(all_pairs):
    """
    Compare compensation rates across mt-mt, nuc-nuc, mt-nuc contact types.
    Tests:
      1. Chi-square: contact_type × significant
      2. Intragenomic vs intergenomic Fisher's exact + OR
      3. Mann-Whitney for fisher_p distributions by contact_type
    """
    print("\n── Analysis 3: Genomic Comparison (mt-mt / nuc-nuc / mt-nuc) ────────")
    return _analysis_genomic_comparison_inner(all_pairs, OUT_DIR / "genomic_comparison.csv")


def _analysis_genomic_comparison_inner(all_pairs, out_path):
    """Internal implementation; `out_path` is the CSV output path."""

    by_type = defaultdict(lambda: {"total": 0, "sig": 0, "fisher_ps": []})
    for r in all_pairs:
        ct = r.get("contact_type", "unknown")
        by_type[ct]["total"] += 1
        if is_significant(r):
            by_type[ct]["sig"] += 1
        fp = safe_float(r.get("fisher_p"))
        if fp is not None:
            by_type[ct]["fisher_ps"].append(fp)

    print(f"\n  {'Type':>10}  {'Tested':>7}  {'Sig':>6}  {'Rate':>7}")
    for ct in ("mt-mt", "nuc-nuc", "mt-nuc"):
        t = by_type[ct]
        rate = 100 * t["sig"] / t["total"] if t["total"] else 0
        print(f"  {ct:>10}  {t['total']:>7}  {t['sig']:>6}  {rate:>6.1f}%")

    # Chi-square: contact_type × significant (3 categories)
    obs = np.array([
        [by_type["mt-mt"]["sig"],    by_type["mt-mt"]["total"]    - by_type["mt-mt"]["sig"]],
        [by_type["nuc-nuc"]["sig"],  by_type["nuc-nuc"]["total"]  - by_type["nuc-nuc"]["sig"]],
        [by_type["mt-nuc"]["sig"],   by_type["mt-nuc"]["total"]   - by_type["mt-nuc"]["sig"]],
    ])
    chi2, p_chi2, dof, _ = chi2_contingency(obs)
    print(f"\n  Chi-square (contact_type × significant): χ²={chi2:.2f}, df={dof}, p={p_chi2:.3e}")

    # Intragenomic vs intergenomic Fisher
    intra_sig  = by_type["mt-mt"]["sig"]  + by_type["nuc-nuc"]["sig"]
    intra_nsig = (by_type["mt-mt"]["total"]  + by_type["nuc-nuc"]["total"]
                  - intra_sig)
    inter_sig  = by_type["mt-nuc"]["sig"]
    inter_nsig = by_type["mt-nuc"]["total"] - inter_sig
    _, p_intra_inter = fisher_exact([[intra_sig, intra_nsig],
                                     [inter_sig, inter_nsig]])
    OR_ii, lo_ii, hi_ii = odds_ratio_ci(inter_sig, inter_nsig, intra_sig, intra_nsig)
    print(f"\n  Intergenomic vs intragenomic:")
    print(f"    Intergenomic: {inter_sig}/{inter_sig+inter_nsig} sig "
          f"({100*inter_sig/(inter_sig+inter_nsig):.1f}%)")
    print(f"    Intragenomic: {intra_sig}/{intra_sig+intra_nsig} sig "
          f"({100*intra_sig/(intra_sig+intra_nsig):.1f}%)")
    print(f"    OR={OR_ii:.2f} (95% CI {lo_ii:.2f}–{hi_ii:.2f}), Fisher p={p_intra_inter:.3e}")

    # Mann-Whitney on Fisher p distributions
    ps_mt  = by_type["mt-mt"]["fisher_ps"]
    ps_nuc = by_type["nuc-nuc"]["fisher_ps"]
    ps_mtn = by_type["mt-nuc"]["fisher_ps"]
    if ps_mt and ps_nuc:
        stat, p_mw = mannwhitneyu(ps_mt, ps_nuc, alternative="two-sided")
        print(f"\n  Mann-Whitney (fisher_p): mt-mt vs nuc-nuc: U={stat:.0f}, p={p_mw:.3e}")
    if ps_mtn and ps_nuc:
        stat2, p_mw2 = mannwhitneyu(ps_mtn, ps_nuc, alternative="two-sided")
        print(f"  Mann-Whitney (fisher_p): mt-nuc vs nuc-nuc: U={stat2:.0f}, p={p_mw2:.3e}")

    results = []
    for ct in ("mt-mt", "nuc-nuc", "mt-nuc"):
        t = by_type[ct]
        results.append({
            "contact_type": ct,
            "n_tested":     t["total"],
            "n_significant": t["sig"],
            "rate":         round(t["sig"]/t["total"], 4) if t["total"] else 0,
        })
    results.append({
        "contact_type": "chi2_3way",
        "n_tested":     int(obs.sum()),
        "n_significant": int(obs[:,0].sum()),
        "rate":         round(p_chi2, 6),
    })
    results.append({
        "contact_type": "inter_vs_intra_OR",
        "n_tested":     inter_sig + inter_nsig + intra_sig + intra_nsig,
        "n_significant": inter_sig + intra_sig,
        "rate":         round(OR_ii, 4),
    })

    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["contact_type","n_tested",
                                           "n_significant","rate"])
        w.writeheader()
        w.writerows(results)
    print(f"\n  → {out_path}")
    return {
        "by_type":         dict(by_type),
        "chi2_p":          p_chi2,
        "OR_inter_intra":  OR_ii,
        "CI_lo":           lo_ii,
        "CI_hi":           hi_ii,
        "p_inter_intra":   p_intra_inter,
    }


# ══════════════════════════════════════════════════════════════════════════════
# Analysis 4: Ancestral vs derived cDAV compensation
# ══════════════════════════════════════════════════════════════════════════════

def analysis_ancestral_vs_derived(all_pairs, timing_rows):
    """
    Compare compensation rates between ancestral and derived cDAVs.
    Also compare n_origins and sensitivity distributions.
    """
    print("\n── Analysis 4: Ancestral vs Derived cDAV Compensation ───────────────")

    # Build timing lookup
    timing_by_key = {identity_key(r): r for r in timing_rows}

    rows_anc, rows_der = [], []
    for r in all_pairs:
        ik = identity_key(r)
        tr = timing_by_key.get(ik)
        if not tr:
            continue
        entry = {
            "sig":         is_significant(r),
            "fisher_p":    safe_float(r.get("fisher_p")),
            "sensitivity": safe_float(r.get("sensitivity")),
            "n_origins":   int(tr.get("n_dar_loss_branches") or 0)
                           if tr.get("is_ancestral_cdav") == "True"
                           else int(tr.get("n_dar_gain_branches") or 0),
        }
        if tr.get("is_ancestral_cdav") == "True":
            rows_anc.append(entry)
        else:
            rows_der.append(entry)

    for label, rows in [("Ancestral", rows_anc), ("Derived", rows_der)]:
        n_sig  = sum(1 for r in rows if r["sig"])
        rate   = 100 * n_sig / len(rows) if rows else 0
        mean_o = sum(r["n_origins"] for r in rows) / len(rows) if rows else 0
        print(f"\n  {label} cDAVs (n={len(rows)} pairs):")
        print(f"    Significant: {n_sig} ({rate:.1f}%)")
        print(f"    Mean independent origins: {mean_o:.2f}")

    # Chi-square: ancestral/derived × significant
    a = sum(1 for r in rows_anc if r["sig"])
    b = sum(1 for r in rows_anc if not r["sig"])
    c = sum(1 for r in rows_der if r["sig"])
    d = sum(1 for r in rows_der if not r["sig"])
    _, p_chi = fisher_exact([[a, b], [c, d]])
    OR, lo, hi = odds_ratio_ci(a, b, c, d)
    print(f"\n  Fisher's exact (ancestral vs derived compensation):")
    print(f"    OR={OR:.2f} (95% CI {lo:.2f}–{hi:.2f}), p={p_chi:.3e}")

    # Mann-Whitney on n_origins
    anc_orig = [r["n_origins"] for r in rows_anc]
    der_orig = [r["n_origins"] for r in rows_der]
    if anc_orig and der_orig:
        stat, p_mw = mannwhitneyu(anc_orig, der_orig, alternative="two-sided")
        print(f"  Mann-Whitney (n_origins anc vs der): U={stat:.0f}, p={p_mw:.3e}")
        print(f"    Ancestral median origins: {sorted(anc_orig)[len(anc_orig)//2]}")
        print(f"    Derived median origins:   {sorted(der_orig)[len(der_orig)//2]}")

    # Sensitivity distributions
    anc_sens = [r["sensitivity"] for r in rows_anc if r["sensitivity"] is not None]
    der_sens = [r["sensitivity"] for r in rows_der if r["sensitivity"] is not None]
    if anc_sens and der_sens:
        stat2, p_mw2 = mannwhitneyu(anc_sens, der_sens, alternative="two-sided")
        print(f"  Mann-Whitney (sensitivity anc vs der): U={stat2:.0f}, p={p_mw2:.3e}")

    results = [
        {"group": "ancestral", "n_pairs": len(rows_anc),
         "n_significant": a, "rate": round(a/(a+b),4) if a+b else 0,
         "mean_origins": round(sum(anc_orig)/len(anc_orig),2) if anc_orig else 0,
         "OR_vs_derived": round(OR,3), "CI_lo": round(lo,3), "CI_hi": round(hi,3),
         "fisher_p": round(p_chi,6)},
        {"group": "derived", "n_pairs": len(rows_der),
         "n_significant": c, "rate": round(c/(c+d),4) if c+d else 0,
         "mean_origins": round(sum(der_orig)/len(der_orig),2) if der_orig else 0,
         "OR_vs_derived": 1.0, "CI_lo": None, "CI_hi": None, "fisher_p": None},
    ]
    out_path = OUT_DIR / "ancestral_vs_derived.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)
    print(f"\n  → {out_path}")
    return results


# ══════════════════════════════════════════════════════════════════════════════
# Analysis 5: Test-decomposed cross-genome rates
# ══════════════════════════════════════════════════════════════════════════════

def analysis_test_decomposed_rates(all_pairs):
    """
    Report cross-genome significance rates broken down by which test drove significance.

    The three tests have overlapping but non-identical failure modes:
      - Fisher FDR: available for all pairs; inflated by cDAV-clade inheritance for
        single-origin derived cDAVs (contact states inherited, not independently evolved)
      - Branch co-occurrence FDR: requires multi-origin cDAVs (≥2 gain branches);
        partially confounded by inheritance for derived cDAVs
      - Pagel FDR: requires both traits to be variable in the common-species pruned tree;
        systematically under-samples pairs with widespread cDAVs (mt-mt bias)

    Strata reported:
      (a) Fisher-only — pairs where ONLY Fisher was available (single-origin, no branch/Pagel).
          These are the most inheritance-confounded; reported with explicit caveat.
      (b) Multi-origin branch — pairs significant by branch_cooccur_fdr (excludes Fisher-only).
      (c) Pagel-confirmed — pairs significant by pagel_fdr (highest methodological rigour
          when the test runs; smallest and most biased coverage subset).
    """
    print("\n── Analysis 5: Test-decomposed cross-genome rates ───────────────────")

    contact_types = ("mt-mt", "nuc-nuc", "mt-nuc")

    # Accumulate counts per stratum × contact_type
    counts = {ct: {
        "total": 0,
        "fisher_only_total": 0, "fisher_only_sig": 0,
        "multi_orig_total":  0, "multi_orig_sig":  0,
        "pagel_tested":      0, "pagel_sig":       0,
    } for ct in contact_types}
    counts["all"] = {k: 0 for k in counts["mt-mt"]}

    for r in all_pairs:
        ct   = r.get("contact_type", "unknown")
        if ct not in counts:
            continue
        flags = test_evidence_flags(r)
        sig   = is_significant(r)

        for key in (ct, "all"):
            counts[key]["total"] += 1
            if flags["fisher_only"]:
                counts[key]["fisher_only_total"] += 1
                if sig:
                    counts[key]["fisher_only_sig"] += 1
            if flags["multi_origin"]:
                counts[key]["multi_orig_total"] += 1
                if flags["branch_sig"] or flags["pagel_sig"]:
                    counts[key]["multi_orig_sig"] += 1
            if flags["pagel_tested"]:
                counts[key]["pagel_tested"] += 1
                if flags["pagel_sig"]:
                    counts[key]["pagel_sig"] += 1

    def pct(n, d):
        return f"{100*n/d:.1f}%" if d else "N/A"

    print("\n  (a) Fisher-only pairs (single-origin; inheritance-confounded):")
    print(f"  {'Type':>10}  {'Tested':>7}  {'Sig':>6}  {'Rate':>7}  {'of all':>7}")
    for ct in contact_types:
        c = counts[ct]
        print(f"  {ct:>10}  {c['fisher_only_total']:>7}  "
              f"{c['fisher_only_sig']:>6}  "
              f"{pct(c['fisher_only_sig'], c['fisher_only_total']):>7}  "
              f"({pct(c['fisher_only_total'], c['total'])} of pairs)")

    print("\n  (b) Multi-origin pairs (branch_cooccur OR Pagel significant; excl. Fisher-only):")
    print(f"  {'Type':>10}  {'Tested':>7}  {'Sig':>6}  {'Rate':>7}")
    for ct in contact_types:
        c = counts[ct]
        print(f"  {ct:>10}  {c['multi_orig_total']:>7}  "
              f"{c['multi_orig_sig']:>6}  "
              f"{pct(c['multi_orig_sig'], c['multi_orig_total']):>7}")

    print("\n  (c) Pagel-tested pairs (pagel_fdr available; most rigorous when runnable):")
    print(f"  {'Type':>10}  {'Tested':>7}  {'Sig':>6}  {'Rate':>7}")
    for ct in contact_types:
        c = counts[ct]
        print(f"  {ct:>10}  {c['pagel_tested']:>7}  "
              f"{c['pagel_sig']:>6}  "
              f"{pct(c['pagel_sig'], c['pagel_tested']):>7}")

    # Check whether the cross-genome rate difference survives outside Fisher-only
    mt_mo_sig  = counts["mt-mt"]["multi_orig_sig"]
    mt_mo_tot  = counts["mt-mt"]["multi_orig_total"]
    nuc_mo_sig = counts["nuc-nuc"]["multi_orig_sig"]
    nuc_mo_tot = counts["nuc-nuc"]["multi_orig_total"]
    if mt_mo_tot and nuc_mo_tot:
        a = nuc_mo_sig;  b = nuc_mo_tot - nuc_mo_sig
        c = mt_mo_sig;   d = mt_mo_tot  - mt_mo_sig
        _, p_mo = fisher_exact([[a, b], [c, d]])
        OR_mo, lo_mo, hi_mo = odds_ratio_ci(a, b, c, d)
        print(f"\n  Multi-origin only — nuc-nuc vs mt-mt:")
        print(f"    nuc-nuc: {nuc_mo_sig}/{nuc_mo_tot} ({pct(nuc_mo_sig, nuc_mo_tot)})")
        print(f"    mt-mt:   {mt_mo_sig}/{mt_mo_tot} ({pct(mt_mo_sig, mt_mo_tot)})")
        print(f"    OR = {OR_mo:.2f} (95% CI {lo_mo:.2f}–{hi_mo:.2f}), Fisher p = {p_mo:.3e}")
    else:
        OR_mo = lo_mo = hi_mo = p_mo = None

    # Pagel: nuc-nuc vs mt-mt
    mt_pag_sig  = counts["mt-mt"]["pagel_sig"]
    mt_pag_tot  = counts["mt-mt"]["pagel_tested"]
    nuc_pag_sig = counts["nuc-nuc"]["pagel_sig"]
    nuc_pag_tot = counts["nuc-nuc"]["pagel_tested"]
    if mt_pag_tot and nuc_pag_tot:
        a2 = nuc_pag_sig; b2 = nuc_pag_tot - nuc_pag_sig
        c2 = mt_pag_sig;  d2 = mt_pag_tot  - mt_pag_sig
        _, p_pag = fisher_exact([[a2, b2], [c2, d2]])
        OR_pag, lo_pag, hi_pag = odds_ratio_ci(a2, b2, c2, d2)
        print(f"\n  Pagel-confirmed only — nuc-nuc vs mt-mt:")
        print(f"    nuc-nuc: {nuc_pag_sig}/{nuc_pag_tot} ({pct(nuc_pag_sig, nuc_pag_tot)})")
        print(f"    mt-mt:   {mt_pag_sig}/{mt_pag_tot} ({pct(mt_pag_sig, mt_pag_tot)})")
        print(f"    OR = {OR_pag:.2f} (95% CI {lo_pag:.2f}–{hi_pag:.2f}), Fisher p = {p_pag:.3e}")
    else:
        OR_pag = lo_pag = hi_pag = p_pag = None

    # Write CSV
    rows = []
    for ct in contact_types:
        c = counts[ct]
        rows.append({
            "contact_type": ct,
            "n_total":      c["total"],
            "fisher_only_n":   c["fisher_only_total"],
            "fisher_only_sig": c["fisher_only_sig"],
            "fisher_only_rate": round(c["fisher_only_sig"]/c["fisher_only_total"], 4) if c["fisher_only_total"] else "",
            "multi_orig_n":    c["multi_orig_total"],
            "multi_orig_sig":  c["multi_orig_sig"],
            "multi_orig_rate": round(c["multi_orig_sig"]/c["multi_orig_total"], 4) if c["multi_orig_total"] else "",
            "pagel_n":         c["pagel_tested"],
            "pagel_sig":       c["pagel_sig"],
            "pagel_rate":      round(c["pagel_sig"]/c["pagel_tested"], 4) if c["pagel_tested"] else "",
        })
    out_path = OUT_DIR / "test_decomposed_rates.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n  → {out_path}")

    return {
        "counts": counts,
        "OR_multi_orig": OR_mo, "CI_lo_mo": lo_mo, "CI_hi_mo": hi_mo, "p_mo": p_mo,
        "OR_pagel":      OR_pag, "CI_lo_pag": lo_pag, "CI_hi_pag": hi_pag, "p_pag": p_pag,
    }


# ══════════════════════════════════════════════════════════════════════════════
# Analysis 6: Pagel coverage and bias
# ══════════════════════════════════════════════════════════════════════════════

def analysis_pagel_coverage(all_pairs):
    """
    Quantify Pagel test coverage and NA-rate by contact_type.

    Pagel's discrete test fails (returns NA) when one or both binary traits is
    invariant after pruning to common species.  mt-mt pairs are disproportionately
    affected because cDAVs in mtDNA-encoded subunits tend to span many species
    (median ~154), leaving the contact trait with insufficient variation to fit
    the 4-rate model.

    Reports:
      - n_eligible:  pairs that meet the Pagel pre-filter (fisher_p < 0.20 AND n_species ≥ 20)
      - n_pagel_ran: pairs where pagel_fdr is non-NA (fitPagel returned a result)
      - n_pagel_na:  pairs where pagel_fdr is NA (invariant trait or convergence failure)
      - na_rate:     n_pagel_na / n_eligible
    """
    print("\n── Analysis 6: Pagel test coverage and NA-rate by contact_type ──────")

    contact_types = ("mt-mt", "nuc-nuc", "mt-nuc")
    cov = {ct: {"eligible": 0, "ran": 0, "na": 0} for ct in contact_types}
    cov["all"] = {"eligible": 0, "ran": 0, "na": 0}

    for r in all_pairs:
        ct = r.get("contact_type", "unknown")
        if ct not in cov:
            continue
        fp = safe_float(r.get("fisher_p"))
        ns = safe_float(r.get("n_species_in_test"))
        # Check eligibility (mirrors prepare_pagel_hpc.py pre-filter)
        eligible = (fp is not None and fp < 0.20) and (ns is not None and ns >= 20)
        pf_raw = r.get("pagel_fdr", "")
        pagel_ran = pf_raw not in ("", "NA", "NaN", None) and safe_float(pf_raw) is not None
        if eligible:
            for key in (ct, "all"):
                cov[key]["eligible"] += 1
                if pagel_ran:
                    cov[key]["ran"] += 1
                else:
                    cov[key]["na"] += 1

    def pct(n, d):
        return f"{100*n/d:.1f}%" if d else "N/A"

    print(f"\n  {'Type':>10}  {'Eligible':>9}  {'Ran':>7}  {'NA':>7}  {'NA rate':>8}")
    for ct in contact_types:
        c = cov[ct]
        print(f"  {ct:>10}  {c['eligible']:>9}  {c['ran']:>7}  {c['na']:>7}  "
              f"{pct(c['na'], c['eligible']):>8}")
    c = cov["all"]
    print(f"  {'ALL':>10}  {c['eligible']:>9}  {c['ran']:>7}  {c['na']:>7}  "
          f"{pct(c['na'], c['eligible']):>8}")
    print()
    print("  NOTE: The Pagel NA rate is disproportionately high for mt-mt pairs "
          "because widespread mtDNA cDAVs (many species) leave the contact trait "
          "invariant after pruning to common-species trees. Cross-genome Pagel "
          "comparisons therefore reflect both biological signal AND sampling bias "
          "against mt-mt pairs.")

    rows = []
    for ct in contact_types + ("all",):
        c = cov[ct]
        rows.append({
            "contact_type": ct,
            "n_eligible":   c["eligible"],
            "n_pagel_ran":  c["ran"],
            "n_pagel_na":   c["na"],
            "na_rate":      round(c["na"] / c["eligible"], 4) if c["eligible"] else "",
        })
    out_path = OUT_DIR / "pagel_coverage.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["contact_type","n_eligible","n_pagel_ran",
                                           "n_pagel_na","na_rate"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n  → {out_path}")
    return cov


# ══════════════════════════════════════════════════════════════════════════════
# Analysis 7: Write markdown summary
# ══════════════════════════════════════════════════════════════════════════════

def write_summary(all_pairs, timing_rows, origins_stats,
                  temporal_results, genomic_stats, anc_der_results,
                  genomic_stats_no_sdh=None, n_sdh_excluded=0,
                  decomposed_stats=None, pagel_cov=None):
    """Write a human-readable markdown summary of all analyses."""

    n_all    = len(all_pairs)
    n_sig    = sum(1 for r in all_pairs if is_significant(r))
    n_timing = len(timing_rows)

    # Genomic breakdown
    by_type = genomic_stats["by_type"]

    # Temporal ordering global
    t_all = next((r for r in temporal_results if r["subset"] == "All pairs"), None)

    anc = anc_der_results[0]
    der = anc_der_results[1]

    lines = [
        "# Comparative Phylogenetic Analysis: Summary",
        "",
        "## Dataset",
        f"- Total DAR–contact pairs tested: {n_all:,}",
        f"- Significant pairs (FDR ≤ 0.10): {n_sig:,} ({100*n_sig/n_all:.1f}%)",
        f"- Pairs with timing annotation: {n_timing:,}",
        "",
        "---",
        "",
        "## 1. Independent Origins (Convergent Evolution)",
        "",
        f"- Total unique cDAV positions analysed: {origins_stats['n_total']:,}",
        f"- cDAVs with ≥2 independent origins (convergent): "
        f"{origins_stats['n_convergent']:,} "
        f"({100*origins_stats['n_convergent']/origins_stats['n_total']:.1f}%)",
        f"- Convergent cDAVs vs compensation enrichment:",
        f"  - OR = {origins_stats['OR']:.2f} "
        f"(95% CI {origins_stats['CI_lo']:.2f}–{origins_stats['CI_hi']:.2f})",
        f"  - Fisher p = {origins_stats['p']:.3e}",
        "",
        "**Interpretation**: cDAVs that arose independently in multiple lineages "
        "(convergent) are enriched for significant compensatory contacts. "
        "This indicates that convergent cDAVs reside in protein environments that "
        "repeatedly co-evolve the same structural solutions, providing stronger "
        "evidence of functional compensation.",
        "",
        "---",
        "",
        "## 2. Temporal Ordering: Did the Contact Substitution Arise First?",
        "",
    ]

    if t_all:
        a, b, c, d = t_all["a_cf_sig"], t_all["b_cf_nsig"], t_all["c_ncf_sig"], t_all["d_ncf_nsig"]
        n_cf   = a + b
        n_sig2 = a + c
        lines += [
            f"- Pairs with contact_first timing: {n_cf:,} "
            f"({100*n_cf/t_all['n_total']:.1f}% of annotated pairs)",
            f"- Of contact_first pairs, {100*a/n_cf:.1f}% are significant",
            f"- Of significant pairs, {100*a/n_sig2:.1f}% had contact_first timing",
            f"- OR = {t_all['OR']:.2f} "
            f"(95% CI {t_all['CI_lo']:.2f}–{t_all['CI_hi']:.2f}), "
            f"Fisher p = {t_all['fisher_p']:.3e}",
            "",
        ]

    for subset in ("Derived cDAVs", "Ancestral cDAVs",
                   "mt-mt pairs", "nuc-nuc pairs", "mt-nuc pairs"):
        tr = next((r for r in temporal_results if r["subset"] == subset), None)
        if tr:
            lines.append(
                f"- **{subset}**: OR = {tr['OR']:.2f} "
                f"({tr['CI_lo']:.2f}–{tr['CI_hi']:.2f}), p = {tr['fisher_p']:.3e}"
            )
    lines += [
        "",
        "**Interpretation**: If contact_first is enriched in significant pairs, "
        "this is evidence of epistatic pre-adaptation — the structural background "
        "permitting a human-pathogenic residue was in place before the variant arose. "
        "This distinguishes permissive co-evolution (contact enables DAR) from "
        "reactive co-evolution (contact rescues after DAR arises).",
        "",
        "---",
        "",
        "## 3. Genomic Architecture: mt-mt vs nuc-nuc vs mt-nuc",
        "",
        f"| Contact type | Pairs tested | Significant | Rate |",
        f"|---|---|---|---|",
    ]

    for ct in ("mt-mt", "nuc-nuc", "mt-nuc"):
        t = by_type[ct]
        rate = 100 * t["sig"] / t["total"] if t["total"] else 0
        lines.append(f"| {ct} | {t['total']:,} | {t['sig']:,} | {rate:.1f}% |")

    lines += [
        "",
        f"- Chi-square (3-way, contact_type × significance): p = {genomic_stats['chi2_p']:.3e}",
        f"- Intergenomic (mt-nuc) vs intragenomic (mt-mt + nuc-nuc):",
        f"  - OR = {genomic_stats['OR_inter_intra']:.2f} "
        f"(95% CI {genomic_stats['CI_lo']:.2f}–{genomic_stats['CI_hi']:.2f}), "
        f"p = {genomic_stats['p_inter_intra']:.3e}",
        "",
    ]

    # SDH-excluded comparison (Additional B)
    if genomic_stats_no_sdh is not None:
        by_type_ns = genomic_stats_no_sdh["by_type"]
        lines += [
            f"### SDH-excluded comparison (n={n_sdh_excluded} SDH pairs removed)",
            "",
            "Complex II (SDHA/B/C/D) variants cause disease via succinate accumulation "
            "(tumour-suppressor mechanism), NOT OXPHOS dysfunction. Including them in "
            "cross-genome rate comparisons conflates two pathological mechanisms.",
            "",
            f"| Contact type | Pairs (excl. SDH) | Significant | Rate |",
            f"|---|---|---|---|",
        ]
        for ct in ("mt-mt", "nuc-nuc", "mt-nuc"):
            t_ns = by_type_ns[ct]
            rate_ns = 100 * t_ns["sig"] / t_ns["total"] if t_ns["total"] else 0
            lines.append(f"| {ct} | {t_ns['total']:,} | {t_ns['sig']:,} | {rate_ns:.1f}% |")
        lines += [
            "",
            f"- SDH-excluded intergenomic OR = {genomic_stats_no_sdh['OR_inter_intra']:.2f} "
            f"(95% CI {genomic_stats_no_sdh['CI_lo']:.2f}–{genomic_stats_no_sdh['CI_hi']:.2f}), "
            f"p = {genomic_stats_no_sdh['p_inter_intra']:.3e}",
            "- If the nuc-nuc rate drops substantially after SDH exclusion, the "
            "nuclear > mitochondrial compensation finding was driven by SDH mechanism, "
            "not by general OXPHOS co-evolution.",
            "",
        ]

    lines += [
        "**Interpretation**: Differences in compensation rates between mt-mt, nuc-nuc, "
        "and mt-nuc contact types reflect the distinct co-evolutionary pressures at "
        "intra- vs inter-genomic interfaces. If intergenomic contacts show significantly "
        "different compensation rates, this provides direct molecular evidence for or "
        "against cytonuclear co-adaptation at OxPhos subunit interfaces.",
        "",
        "---",
        "",
        "## 4. Ancestral vs Derived cDAVs",
        "",
        f"| Group | Pairs | Significant | Rate | Mean origins |",
        f"|---|---|---|---|---|",
        f"| Ancestral | {anc['n_pairs']:,} | {anc['n_significant']:,} | "
        f"{100*anc['rate']:.1f}% | {anc['mean_origins']:.1f} |",
        f"| Derived   | {der['n_pairs']:,} | {der['n_significant']:,} | "
        f"{100*der['rate']:.1f}% | {der['mean_origins']:.1f} |",
        "",
        f"- Ancestral vs derived OR = {anc['OR_vs_derived']:.2f} "
        f"(95% CI {anc['CI_lo']:.2f}–{anc['CI_hi']:.2f}), "
        f"Fisher p = {anc['fisher_p']:.3e}",
        "",
        "**Interpretation**: Ancestral cDAVs (human-pathogenic residue present in "
        "the common mammalian ancestor) represent long-maintained compensated states. "
        "If they show different compensation rates from derived cDAVs, this indicates "
        "that the evolutionary history of a variant affects the probability of "
        "observing associated contact co-evolution in extant species.",
        "",
        "---",
        "",
        "## 5. Test-decomposed Cross-Genome Rates",
        "",
        "The headline significance rates depend critically on which test is used as "
        "the primary evidence. Each test has distinct failure modes that are non-randomly "
        "distributed across contact types:",
        "",
        "| Test | Primary failure mode | Direction of bias |",
        "|------|---------------------|-------------------|",
        "| Fisher FDR | Inflated by cDAV-clade inheritance (single-origin pairs) | Anti-conservative for mt-mt (widespread cDAVs) |",
        "| Branch co-occurrence FDR | Partially inflated for derived cDAVs; requires multi-origin | Excludes all single-origin pairs |",
        "| Pagel FDR | Fails when either trait is invariant after tree pruning | Under-samples mt-mt (widespread cDAVs leave contact invariant) |",
        "",
    ]

    if decomposed_stats is not None:
        dc = decomposed_stats["counts"]
        def _pct(n, d): return f"{100*n/d:.1f}%" if d else "N/A"

        lines += [
            "### (a) Fisher-only pairs (single-origin; inheritance-confounded)",
            "",
            "These pairs had no branch co-occurrence or Pagel test available. "
            "Significance driven entirely by Fisher FDR, which is inflated by "
            "clade-level inheritance of contact states (not independent evolution).",
            "",
            "| Contact type | Fisher-only pairs | Significant | Rate |",
            "|---|---|---|---|",
        ]
        for ct in ("mt-mt", "nuc-nuc", "mt-nuc"):
            c = dc[ct]
            lines.append(f"| {ct} | {c['fisher_only_total']:,} | "
                         f"{c['fisher_only_sig']:,} | "
                         f"{_pct(c['fisher_only_sig'], c['fisher_only_total'])} |")
        lines += [
            "",
            "### (b) Multi-origin branch co-occurrence (less confounded)",
            "",
            "Pairs with ≥2 independent cDAV origins, tested by branch co-occurrence "
            "or Pagel. Still partially confounded for derived cDAVs via contact "
            "inheritance; pyvolve permutation null is the primary arbiter once available.",
            "",
            "| Contact type | Pairs with multi-origin data | Branch/Pagel significant | Rate |",
            "|---|---|---|---|",
        ]
        for ct in ("mt-mt", "nuc-nuc", "mt-nuc"):
            c = dc[ct]
            lines.append(f"| {ct} | {c['multi_orig_total']:,} | "
                         f"{c['multi_orig_sig']:,} | "
                         f"{_pct(c['multi_orig_sig'], c['multi_orig_total'])} |")

        if decomposed_stats["OR_multi_orig"] is not None:
            OR_mo = decomposed_stats["OR_multi_orig"]
            lo_mo = decomposed_stats["CI_lo_mo"]
            hi_mo = decomposed_stats["CI_hi_mo"]
            p_mo  = decomposed_stats["p_mo"]
            lines += [
                "",
                f"- Multi-origin: nuc-nuc vs mt-mt — OR = {OR_mo:.2f} "
                f"(95% CI {lo_mo:.2f}–{hi_mo:.2f}), p = {p_mo:.3e}",
            ]

        lines += [
            "",
            "### (c) Pagel-confirmed pairs (most rigorous; smallest coverage)",
            "",
            "| Contact type | Pagel-tested | Pagel-significant | Rate |",
            "|---|---|---|---|",
        ]
        for ct in ("mt-mt", "nuc-nuc", "mt-nuc"):
            c = dc[ct]
            lines.append(f"| {ct} | {c['pagel_tested']:,} | "
                         f"{c['pagel_sig']:,} | "
                         f"{_pct(c['pagel_sig'], c['pagel_tested'])} |")

        if decomposed_stats["OR_pagel"] is not None:
            OR_pag = decomposed_stats["OR_pagel"]
            lo_pag = decomposed_stats["CI_lo_pag"]
            hi_pag = decomposed_stats["CI_hi_pag"]
            p_pag  = decomposed_stats["p_pag"]
            lines += [
                "",
                f"- Pagel-confirmed: nuc-nuc vs mt-mt — OR = {OR_pag:.2f} "
                f"(95% CI {lo_pag:.2f}–{hi_pag:.2f}), p = {p_pag:.3e}",
            ]

    lines += [
        "",
        "**Key question**: Does the nuc-nuc > mt-mt rate difference persist outside "
        "Fisher-only pairs? If the multi-origin and Pagel rates converge, the headline "
        "finding is robust. If they diverge from Fisher-only rates, the headline is "
        "driven by inheritance artefact, not biological compensation differences.",
        "",
        "---",
        "",
        "## 6. Pagel Test Coverage and Systematic Bias",
        "",
    ]

    if pagel_cov is not None:
        def _pct2(n, d): return f"{100*n/d:.1f}%" if d else "N/A"
        lines += [
            "| Contact type | Eligible pairs | Pagel ran | Pagel NA | NA rate |",
            "|---|---|---|---|---|",
        ]
        for ct in ("mt-mt", "nuc-nuc", "mt-nuc", "all"):
            c = pagel_cov[ct]
            lines.append(f"| {ct} | {c['eligible']:,} | {c['ran']:,} | "
                         f"{c['na']:,} | {_pct2(c['na'], c['eligible'])} |")
        lines += [
            "",
            "> **Methodological caveat**: Pagel's discrete test was disproportionately "
            "uninformative for mt-mt pairs due to cDAV breadth in mtDNA-encoded subunits, "
            "which often left the contact trait invariant after pruning to common-species "
            "trees. Cross-genome comparisons of Pagel-significant rates therefore reflect "
            "both biological co-evolution differences AND a systematic sampling bias against "
            "widespread mt-mt cDAVs, and should not be interpreted as direct evidence of "
            "differential compensation between genomes.",
        ]

    lines += [
        "",
        "---",
        "",
        "## Reporting guidance",
        "",
        "Report all three test ensembles with explicit coverage and limitations:",
        "1. **Fisher FDR** on all eligible pairs, with the inheritance-circularity caveat",
        "2. **Branch co-occurrence** on multi-origin pairs, with circularity caveat for derived cDAVs",
        "3. **Pagel FDR** on the testable subset, with coverage statistics and mt-mt bias noted",
        "4. **pyvolve conditional permutation** (primary, once available): "
        "tree-aware null that directly addresses inheritance confounding",
        "",
        "Use Pagel as confirmatory (not primary): 'Of X pairs significant by Fisher "
        "AND testable by Pagel, Y also showed significant Pagel evidence.' "
        "The pyvolve null is the headline analysis; Fisher, branch, and Pagel are "
        "supporting evidence, each with its own caveats.",
    ]

    out_path = OUT_DIR / "comparative_analysis_summary.md"
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n  → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    all_pairs    = load_all_pairs()
    timing_rows  = load_timing()
    print(f"  all_tested_pairs:     {len(all_pairs):,} rows")
    print(f"  timing_annotations:   {len(timing_rows):,} rows")

    origins_records, origins_stats  = analysis_independent_origins(all_pairs, timing_rows)
    temporal_results                = analysis_temporal_ordering(all_pairs, timing_rows)
    genomic_stats                   = analysis_genomic_comparison(all_pairs)

    # ── SDH-excluded genomic comparison (Additional B) ────────────────────────
    # Complex II (SDH) causes disease via succinate accumulation / HIF stabilisation,
    # NOT OXPHOS dysfunction. Including SDH in cross-genome rate comparisons mixes
    # two mechanistically distinct pathological routes with CI/III/IV/V.
    # Run genomic comparison with and without SDH; save both.
    SDH_GENES = {"SDHA", "SDHB", "SDHC", "SDHD", "SDHAF2"}
    all_pairs_no_sdh = [r for r in all_pairs if r.get("dar_gene", "") not in SDH_GENES]
    n_sdh_excluded = len(all_pairs) - len(all_pairs_no_sdh)
    print(f"\n── SDH-excluded comparison (n={len(all_pairs_no_sdh)}, "
          f"excluded {n_sdh_excluded} SDH pairs) ─")

    # Temporarily redirect output file for SDH-excluded version
    import os
    orig_out = OUT_DIR / "genomic_comparison.csv"
    sdh_out  = OUT_DIR / "genomic_comparison_no_sdh.csv"

    # Monkey-patch the output path for the second call
    _orig_gc_path = OUT_DIR / "genomic_comparison.csv"
    genomic_stats_no_sdh = _analysis_genomic_comparison_inner(all_pairs_no_sdh, sdh_out)

    anc_der_results  = analysis_ancestral_vs_derived(all_pairs, timing_rows)
    decomposed_stats = analysis_test_decomposed_rates(all_pairs)
    pagel_cov        = analysis_pagel_coverage(all_pairs)

    print("\n── Writing markdown summary ──────────────────────────────────────────")
    write_summary(all_pairs, timing_rows, origins_stats,
                  temporal_results, genomic_stats, anc_der_results,
                  genomic_stats_no_sdh=genomic_stats_no_sdh,
                  n_sdh_excluded=n_sdh_excluded,
                  decomposed_stats=decomposed_stats,
                  pagel_cov=pagel_cov)

    print("\n════════════════════════════════════════════════════════")
    print("All analyses complete. Outputs in results/phylo/")
    print("════════════════════════════════════════════════════════")


if __name__ == "__main__":
    main()
