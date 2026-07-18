"""
src/phylo/09_pyvolve_rate_controlled.py

Step 3.2 — Rate-controlled Pyvolve OR comparison: mtDNA vs non-SDH nuclear.

Tests whether non-SDH nuclear DAR pairs have higher compensatory OR than
mtDNA DAR pairs *after controlling for n_cdav_spp* (lower n_cdav_spp inflates
OR because the comparison group is smaller).

Stratifies by n_cdav_spp bins: [1-5], [6-20], [21-50], [51-200], [>200].
Within each bin: compares median OR and sig% for mtDNA vs non-SDH nuclear.
Runs a Cochran-Mantel-Haenszel (CMH) test for the combined within-bin OR.

Output: results/phylo/pyvolve_rate_controlled.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

CP_FILE  = ROOT / "results" / "phylo" / "conditional_permissiveness.csv"
OUT_FILE = ROOT / "results" / "phylo" / "pyvolve_rate_controlled.csv"

MTDNA = {"MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
         "MT-CYB","MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8"}
SDH   = {"SDHA","SDHB","SDHC","SDHD","SDHAF2"}

BINS = [(1, 5), (6, 20), (21, 50), (51, 200), (201, 10000)]
BIN_LABELS = ["1-5", "6-20", "21-50", "51-200", ">200"]


def cochran_mantel_haenszel(tables: list[np.ndarray]) -> tuple[float, float]:
    """CMH test for a series of 2x2 tables [[a,b],[c,d]]."""
    numerator = 0.0
    denominator = 0.0
    for t in tables:
        a, b, c, d = t[0,0], t[0,1], t[1,0], t[1,1]
        n = a + b + c + d
        if n < 2:
            continue
        numerator   += a - (a + b) * (a + c) / n
        denominator += (a + b) * (c + d) * (a + c) * (b + d) / (n**2 * (n - 1))
    if denominator == 0:
        return float("nan"), float("nan")
    from scipy.stats import chi2
    chi2_stat = numerator**2 / denominator
    p = chi2.sf(chi2_stat, df=1)
    return float(chi2_stat), float(p)


def main() -> None:
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    cp = pd.read_csv(CP_FILE)
    cp["dar_compartment"] = cp["dar_gene"].apply(
        lambda g: "mtDNA" if g in MTDNA else "nucDNA")
    cp["is_sdh"] = cp["dar_gene"].isin(SDH)
    cp["sig"] = cp["perm_p"] < 0.05
    cp["log_or"] = np.log(cp["observed_or"].clip(lower=0.001))

    mt   = cp[cp["dar_compartment"] == "mtDNA"]
    nuc  = cp[(cp["dar_compartment"] == "nucDNA") & ~cp["is_sdh"]]

    print(f"mtDNA pairs:           {len(mt)}")
    print(f"Non-SDH nuclear pairs: {len(nuc)}")
    print()

    rows = []
    cmh_tables = []

    print(f"{'Bin':<10} {'Genome':<10} {'n':<6} {'sig%':<10} {'med_OR':<10} {'p_bin'}")
    print("-" * 60)

    for (lo, hi), label in zip(BINS, BIN_LABELS):
        for subset, tag in [(mt, "mtDNA"), (nuc, "non-SDH-nuc")]:
            grp = subset[(subset["n_cdav_spp"] >= lo) & (subset["n_cdav_spp"] <= hi)]
            if len(grp) == 0:
                rows.append({"bin": label, "genome": tag, "n": 0,
                             "n_sig": 0, "pct_sig": np.nan,
                             "median_or": np.nan, "fisher_p": np.nan})
                continue
            n_sig = grp["sig"].sum()
            pct   = 100 * n_sig / len(grp)
            med   = grp["observed_or"].median()
            rows.append({"bin": label, "genome": tag, "n": len(grp),
                         "n_sig": int(n_sig), "pct_sig": pct,
                         "median_or": med, "fisher_p": np.nan})
            print(f"  {label:<8} {tag:<12} {len(grp):<6} {pct:<10.1f} {med:<10.2f}")

        # Within-bin CMH: sig vs non-sig, mtDNA vs nuclear
        mt_bin  = mt[(mt["n_cdav_spp"] >= lo)  & (mt["n_cdav_spp"] <= hi)]
        nuc_bin = nuc[(nuc["n_cdav_spp"] >= lo) & (nuc["n_cdav_spp"] <= hi)]
        a = mt_bin["sig"].sum()
        b = nuc_bin["sig"].sum()
        c = len(mt_bin) - a
        d = len(nuc_bin) - b
        if (a + b + c + d) > 0:
            cmh_tables.append(np.array([[b, a],[d, c]]))  # nuc vs mt (nuc=row1)

    print()
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))

    # Overall OR difference (Mann-Whitney on log-OR)
    from scipy.stats import mannwhitneyu
    stat, p_mw = mannwhitneyu(nuc["log_or"].dropna(),
                               mt["log_or"].dropna(), alternative="greater")
    print(f"\nMann-Whitney (nuc log-OR > mtDNA log-OR): U={stat:.0f}, p={p_mw:.3e}")
    print(f"  Median OR — nuclear: {nuc['observed_or'].median():.3f}  "
          f"mtDNA: {mt['observed_or'].median():.3f}")

    # CMH test
    if cmh_tables:
        chi2_cmh, p_cmh = cochran_mantel_haenszel(cmh_tables)
        print(f"\nCMH test (nuc sig% > mtDNA, stratified by n_cdav_spp bins):")
        print(f"  χ²={chi2_cmh:.3f}, p={p_cmh:.3e}")

    # Correlation n_cdav_spp vs log_or
    from scipy.stats import spearmanr
    r, p_r = spearmanr(cp["n_cdav_spp"], cp["log_or"])
    print(f"\nSpearman r (n_cdav_spp vs log_OR, all pairs): r={r:.3f}, p={p_r:.3e}")

    df.to_csv(OUT_FILE, index=False)
    print(f"\nSaved {len(df)} rows → {OUT_FILE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
