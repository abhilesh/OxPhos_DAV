#!/usr/bin/env python3
"""
src/mutagenesis/01_foldx_ddg.py

Run FoldX ΔΔG calculations on the top-ranked compensatory pairs.

Protocol by contact category:
  Intraprotein (dar_gene == contact_gene, 85.7%):
    FoldX BuildModel → ΔΔG_stability for single and double mutants.
    Reports ΔΔG_rescue_stab = ΔΔG_stab_DAR − ΔΔG_stab_double
            ΔΔG_epistasis   = ΔΔG_stab_double − (ΔΔG_stab_DAR + ΔΔG_stab_contact)

  Interprotein (dar_gene != contact_gene, 14.3%):
    Same BuildModel runs PLUS AnalyseComplex on each mutant vs. WT
    to capture change in binding affinity between chains.
    Reports ΔΔG_rescue_bind = ΔΔG_bind_DAR − ΔΔG_bind_double

Prerequisites:
  FoldX binary must be present at $FOLDX_PATH or tools/foldx/foldx.
  Download from https://foldxsuite.crg.eu (free academic registration).

Accuracy note:
  FoldX AnalyseComplex: R ≈ 0.70 with experimental ΔΔG_binding (RMSE ≈ 1.25 kcal/mol).
  Structures at 4 Å resolution (9I4I, Complex I) have higher prediction noise.
  For top interprotein pairs, validate with mCSM-PPI2 (R ≈ 0.82):
    https://biosig.unimelb.edu.au/mcsm_ppi2/

Outputs:
  results/mutagenesis/foldx_ddg.csv

Usage:
  FOLDX_PATH=/path/to/foldx docker run --rm -v $(pwd):/app \\
      -e FOLDX_PATH=$FOLDX_PATH oxphos_dav_analysis \\
      conda run -n oxphos_dav python src/mutagenesis/01_foldx_ddg.py
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

# ─── Paths ────────────────────────────────────────────────────────────────────
# NOTE: these are module-level defaults for local/interactive use. For SLURM
# array / HPC scale-up runs, main() rebinds all of these (via --work-dir,
# --input, --workdir-name) *before* any path is touched, so downstream code
# that references the module globals still works unmodified.
TOP_TARGETS_CSV = ROOT / "results" / "mutagenesis" / "top_targets.csv"
CONTACTS_CSV    = ROOT / "results" / "structural" / "dar_contacts_cbcb8A.csv"
STRUCTURES_DIR  = ROOT / "data" / "structures"
WORKDIR         = ROOT / "results" / "mutagenesis" / "foldx_work"
OUT_DIR         = ROOT / "results" / "mutagenesis"

OUT_DIR.mkdir(parents=True, exist_ok=True)
WORKDIR.mkdir(parents=True, exist_ok=True)

# FoldX binary resolution — try explicit env var, then versioned name, then plain name
def _find_foldx() -> str:
    env = os.environ.get("FOLDX_PATH")
    if env and Path(env).exists():
        return env
    foldx_dir = ROOT / "tools" / "foldx"
    # Prefer versioned binary (e.g. foldx_20270131) over plain 'foldx'
    for candidate in sorted(foldx_dir.glob("foldx_*"), reverse=True):
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    plain = foldx_dir / "foldx"
    if plain.exists():
        return str(plain)
    return str(foldx_dir / "foldx")  # will trigger error in check_foldx()

_FOLDX_PATH = _find_foldx()
# FoldX 5 uses a molecules/ directory (co-located with binary) instead of rotabase.txt
_FOLDX_DIR  = Path(_FOLDX_PATH).parent

# ─── Verification ─────────────────────────────────────────────────────────────
def check_foldx() -> None:
    if not Path(_FOLDX_PATH).exists():
        sys.exit(
            f"\nERROR: FoldX binary not found at {_FOLDX_PATH}\n"
            "  Download FoldX from https://foldxsuite.crg.eu (free academic)\n"
            "  Place the binary at tools/foldx/foldx  OR\n"
            "  Set FOLDX_PATH=/path/to/foldx in the environment"
        )
    print(f"FoldX binary: {_FOLDX_PATH}")


# ─── CIF → PDB conversion ─────────────────────────────────────────────────────
_pdb_converted: dict[str, Path] = {}


def cif_to_pdb(pdb_id: str, chains: list[str]) -> Path | None:
    """
    Convert CIF to PDB extracting only the requested chains.
    Returns path to the output PDB file.
    """
    key = (pdb_id, tuple(sorted(set(chains))))
    if key in _pdb_converted:
        return _pdb_converted[key]

    cif_path = STRUCTURES_DIR / f"{pdb_id}.cif"
    if not cif_path.exists():
        _pdb_converted[key] = None
        return None

    chain_tag = "_".join(sorted(set(chains)))
    out_pdb = WORKDIR / f"{pdb_id}_{chain_tag}.pdb"

    if not out_pdb.exists():
        try:
            from Bio.PDB import MMCIFParser, PDBIO, Select

            class ChainSelect(Select):
                def accept_chain(self, chain):
                    return chain.id in chains

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                parser = MMCIFParser(QUIET=True)
                structure = parser.get_structure(pdb_id, str(cif_path))

            io = PDBIO()
            io.set_structure(structure)
            io.save(str(out_pdb), ChainSelect())
        except Exception as e:
            print(f"  WARNING: CIF→PDB conversion failed for {pdb_id}: {e}")
            _pdb_converted[key] = None
            return None

    _pdb_converted[key] = out_pdb
    return out_pdb


# ─── FoldX helpers ────────────────────────────────────────────────────────────
def _ensure_molecules_dir(work_dir: Path) -> None:
    """
    FoldX 5 looks for molecules/ relative to cwd. Copy it into the work
    directory (428 KB; symlinks fail on macOS Docker bind mounts).
    """
    dest   = work_dir / "molecules"
    source = _FOLDX_DIR / "molecules"
    if not dest.exists() and source.exists():
        shutil.copytree(str(source), str(dest))


def run_foldx_cmd(args: list[str], cwd: Path) -> subprocess.CompletedProcess:
    _ensure_molecules_dir(cwd)
    cmd = [_FOLDX_PATH] + args
    return subprocess.run(
        cmd, cwd=str(cwd), capture_output=True, text=True, timeout=600
    )


def repair_pdb(pdb_path: Path) -> Path | None:
    """Run FoldX RepairPDB. Returns path to repaired PDB."""
    repaired = pdb_path.parent / f"{pdb_path.stem}_Repair.pdb"
    if repaired.exists():
        return repaired

    run_foldx_cmd(["--command=RepairPDB", f"--pdb={pdb_path.name}"], cwd=pdb_path.parent)
    if repaired.exists():
        return repaired
    return None


def parse_buildmodel_ddg(fxout_path: Path) -> tuple[float, float] | None:
    """
    Parse (ΔΔG_mean, ΔΔG_sd) from Average_*.fxout (FoldX 5 tab-delimited format).

    FoldX 5 header: Pdb  SD  total energy  Backbone Hbond  ...
    col 0 = Pdb name, col 1 = SD across replicates, col 2 = mean ΔΔG.
    Returns (ddg_mean, ddg_sd) or None on failure.
    """
    try:
        with open(fxout_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                if parts[0].startswith("Pdb") or not parts[0].strip():
                    continue
                try:
                    ddg = float(parts[2])   # col 2 = mean total energy = ΔΔG
                    sd  = float(parts[1])   # col 1 = SD across replicates
                    return (ddg, sd)
                except ValueError:
                    continue
    except Exception:
        pass
    return None


def parse_buildmodel_dif_ddg(fxout_path: Path) -> list[float]:
    """
    Parse individual replicate ΔΔG values from Dif_*.fxout.

    FoldX 5 Dif format: Pdb  total energy  ...  (no SD column)
    col 0 = Pdb name (e.g. 8H9S_N_Repair_1_0.pdb), col 1 = total energy.
    Returns list of up to 3 individual values; empty list on failure.
    """
    vals = []
    try:
        with open(fxout_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                if parts[0].startswith("Pdb") or not parts[0].strip():
                    continue
                try:
                    vals.append(float(parts[1]))
                except ValueError:
                    continue
    except Exception:
        pass
    return vals


def parse_analysecomplex_dg(fxout_path: Path) -> float | None:
    """Parse ΔG_binding from Summary_*.fxout (column 'Interaction Energy')."""
    try:
        with open(fxout_path) as f:
            header = None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if header is None:
                    header = line.split("\t")
                    try:
                        ie_col = header.index("Interaction Energy")
                    except ValueError:
                        ie_col = 5  # fallback position
                    continue
                parts = line.split("\t")
                if len(parts) > ie_col:
                    try:
                        return float(parts[ie_col])
                    except ValueError:
                        continue
    except Exception:
        pass
    return None


def run_buildmodel(
    repaired_pdb: Path, mutation_str: str, run_id: str, n_runs: int = 3
) -> tuple[float | None, float | None, list[float]]:
    """
    Run FoldX BuildModel for one mutation set.
    mutation_str: FoldX notation, e.g. "SA3P;" or "SA3P,VB12F;"
    Returns (ddg_mean, ddg_sd, individual_vals) — any element may be None/[].
    """
    work = WORKDIR / run_id
    work.mkdir(exist_ok=True)

    # Copy repaired PDB into work dir
    pdb_copy = work / repaired_pdb.name
    shutil.copy(repaired_pdb, pdb_copy)

    # Write individual_list.txt
    mut_file = work / "individual_list.txt"
    mut_file.write_text(mutation_str + "\n")

    run_foldx_cmd(
        ["--command=BuildModel",
         f"--pdb={repaired_pdb.name}",
         "--mutant-file=individual_list.txt",
         f"--numberOfRuns={n_runs}"],
        cwd=work,
    )

    avg_files = list(work.glob(f"Average_{repaired_pdb.stem}*.fxout"))
    dif_files = list(work.glob(f"Dif_{repaired_pdb.stem}*.fxout"))

    result_tuple = parse_buildmodel_ddg(avg_files[0]) if avg_files else None
    if result_tuple is None:
        return (None, None, [])
    ddg_mean, ddg_sd = result_tuple
    individual = parse_buildmodel_dif_ddg(dif_files[0]) if dif_files else []
    return (ddg_mean, ddg_sd, individual)


def run_analysecomplex(pdb_path: Path, chain1: str, chain2: str, run_id: str) -> float | None:
    """
    Run FoldX AnalyseComplex on a (mutant) PDB for the specified chain pair.
    Returns binding energy (ΔG, kcal/mol).
    """
    work = WORKDIR / run_id
    work.mkdir(exist_ok=True)

    pdb_copy = work / pdb_path.name
    shutil.copy(pdb_path, pdb_copy)

    res = run_foldx_cmd(
        ["--command=AnalyseComplex",
         f"--pdb={pdb_path.name}",
         f"--analyseComplexChains={chain1},{chain2}"],
        cwd=work,
    )

    sum_files = list(work.glob("Summary_*.fxout"))
    if not sum_files:
        return None
    return parse_analysecomplex_dg(sum_files[0])


# ─── Per-pair ΔΔG calculation ──────────────────────────────────────────────────
def foldx_notation(ref_aa: str, chain: str, resnum, alt_aa: str) -> str:
    """FoldX single mutation: {ref}{chain}{resnum}{alt}"""
    return f"{ref_aa}{chain}{int(float(resnum))}{alt_aa}"


def compute_ddg_for_pair(row: pd.Series, pair_idx: int) -> dict:
    """
    Run FoldX for one compensatory pair.
    Returns a dict of ΔΔG metrics.
    """
    result = {
        "ann_id":             row["ann_id"],
        "dar_gene":           row["dar_gene"],
        "dar_aa_coord":       row["dar_aa_coord"],
        "dar_ref_aa":         row.get("dar_ref_aa"),
        "dar_alt_aa":         row["dar_alt_aa"],
        "contact_gene":       row["contact_gene"],
        "contact_refseq_pos": row["contact_refseq_pos"],
        "contact_human_aa":   row.get("contact_human_aa"),
        "contact_alt_aa":     row["contact_alt_aa"],
        "pdb_id":             row.get("pdb_id"),
        "dar_chain":          row.get("dar_chain"),
        "contact_chain":      row.get("contact_chain"),
        "is_intraprotein":    row.get("is_intraprotein", row["dar_gene"] == row["contact_gene"]),
        "ddg_stab_dar":         np.nan,
        "ddg_stab_dar_sd":      np.nan,
        "ddg_stab_contact":     np.nan,
        "ddg_stab_contact_sd":  np.nan,
        "ddg_stab_double":      np.nan,
        "ddg_stab_double_sd":   np.nan,
        "ddg_rescue_stab":      np.nan,
        "ddg_epistasis":        np.nan,
        "ddg_epistasis_sd":     np.nan,
        "ddg_epistasis_sig":    False,   # |epistasis| > 2 * epistasis_sd
        "ddg_bind_dar":         np.nan,
        "ddg_bind_double":      np.nan,
        "ddg_rescue_bind":      np.nan,
        "foldx_status":         "skipped",
    }

    pdb_id    = str(row.get("pdb_id", ""))
    dar_chain = str(row.get("dar_chain", "")).strip()
    dar_res   = row.get("dar_struct_res")
    con_chain = str(row.get("contact_chain", "")).strip()
    con_res   = row.get("contact_resnum")

    if not pdb_id or pd.isna(row.get("pdb_id")):
        result["foldx_status"] = "no_pdb"
        return result
    if not dar_chain or pd.isna(dar_res) or not con_chain or pd.isna(con_res):
        result["foldx_status"] = "missing_coords"
        return result

    # Amino acid sanity: need single-letter codes
    dar_ref = str(row.get("dar_ref_aa", "")).strip()
    dar_alt = str(row.get("dar_alt_aa", "")).strip()
    con_ref = str(row.get("contact_human_aa", "")).strip()
    con_alt = str(row.get("contact_alt_aa", "")).strip()

    if not all(len(x) == 1 for x in [dar_ref, dar_alt, con_ref, con_alt]):
        result["foldx_status"] = "invalid_aa"
        return result

    # Build PDB (chain extraction + RepairPDB)
    chains = list({dar_chain, con_chain})
    raw_pdb = cif_to_pdb(pdb_id, chains)
    if raw_pdb is None:
        result["foldx_status"] = "cif_conversion_failed"
        return result

    repaired = repair_pdb(raw_pdb)
    if repaired is None:
        result["foldx_status"] = "repair_failed"
        return result

    prefix = f"pair{pair_idx:04d}"

    # ─── BuildModel: 3 mutation conditions ────────────────────────────────
    mut_dar     = foldx_notation(dar_ref, dar_chain, dar_res, dar_alt)
    mut_contact = foldx_notation(con_ref, con_chain, con_res, con_alt)
    mut_double  = f"{mut_dar},{mut_contact}"

    ddg_dar,     sd_dar,     _ = run_buildmodel(repaired, f"{mut_dar};",     f"{prefix}_dar")
    ddg_contact, sd_contact, _ = run_buildmodel(repaired, f"{mut_contact};", f"{prefix}_contact")
    ddg_double,  sd_double,  _ = run_buildmodel(repaired, f"{mut_double};",  f"{prefix}_double")

    result["ddg_stab_dar"]        = ddg_dar
    result["ddg_stab_dar_sd"]     = sd_dar
    result["ddg_stab_contact"]    = ddg_contact
    result["ddg_stab_contact_sd"] = sd_contact
    result["ddg_stab_double"]     = ddg_double
    result["ddg_stab_double_sd"]  = sd_double

    if all(v is not None and pd.notna(v) for v in [ddg_dar, ddg_double]):
        result["ddg_rescue_stab"] = float(ddg_dar) - float(ddg_double)

    if all(v is not None and pd.notna(v) for v in [ddg_dar, ddg_contact, ddg_double]):
        epi = float(ddg_double) - (float(ddg_dar) + float(ddg_contact))
        result["ddg_epistasis"] = epi
        # Propagated SD: sqrt(SD_DAR² + SD_contact² + SD_double²)
        if all(v is not None and pd.notna(v) for v in [sd_dar, sd_contact, sd_double]):
            epi_sd = float(np.sqrt(sd_dar**2 + sd_contact**2 + sd_double**2))
            result["ddg_epistasis_sd"]  = epi_sd
            result["ddg_epistasis_sig"] = (epi_sd > 0 and abs(epi) > 2 * epi_sd)

    result["foldx_status"] = "stability_done"

    # ─── AnalyseComplex: interprotein pairs only ───────────────────────────
    is_intra = (row["dar_gene"] == row["contact_gene"])
    if not is_intra and dar_chain != con_chain:
        # WT binding energy
        dg_wt = run_analysecomplex(repaired, dar_chain, con_chain, f"{prefix}_ac_wt")

        # DAR mutant
        dar_mut_pdb = WORKDIR / f"{prefix}_dar" / f"WT_{repaired.stem}_1.pdb"
        double_mut_pdb = WORKDIR / f"{prefix}_double" / f"WT_{repaired.stem}_1.pdb"

        dg_dar_mut    = run_analysecomplex(dar_mut_pdb, dar_chain, con_chain, f"{prefix}_ac_dar")    if dar_mut_pdb.exists() else None
        dg_double_mut = run_analysecomplex(double_mut_pdb, dar_chain, con_chain, f"{prefix}_ac_double") if double_mut_pdb.exists() else None

        if all(v is not None and pd.notna(v) for v in [dg_wt, dg_dar_mut]):
            result["ddg_bind_dar"]    = float(dg_dar_mut) - float(dg_wt)
        if all(v is not None and pd.notna(v) for v in [dg_wt, dg_double_mut]):
            result["ddg_bind_double"] = float(dg_double_mut) - float(dg_wt)
        if pd.notna(result["ddg_bind_dar"]) and pd.notna(result["ddg_bind_double"]):
            result["ddg_rescue_bind"] = float(result["ddg_bind_dar"]) - float(result["ddg_bind_double"])

        result["foldx_status"] = "full"

    return result


# ─── Reparse mode ─────────────────────────────────────────────────────────────
def reparse_existing() -> None:
    """
    Re-read Average_*.fxout and Dif_*.fxout from existing foldx_work/pair*
    directories to extract SD values without re-running FoldX.
    Updates foldx_ddg.csv in-place with new SD and epistasis_sd columns.
    """
    existing_csv = OUT_DIR / "foldx_ddg.csv"
    if not existing_csv.exists():
        sys.exit(f"\nERROR: {existing_csv} not found. Run full FoldX pipeline first.")

    print(f"Reparsing existing FoldX outputs from {WORKDIR}")
    df = pd.read_csv(existing_csv, dtype={"ann_id": str})
    print(f"  {len(df)} rows in foldx_ddg.csv")

    # Initialise new columns
    for col in ["ddg_stab_dar_sd", "ddg_stab_contact_sd", "ddg_stab_double_sd",
                "ddg_epistasis_sd", "ddg_epistasis_sig"]:
        if col not in df.columns:
            df[col] = np.nan if col != "ddg_epistasis_sig" else False

    for i, row in df.iterrows():
        prefix = f"pair{i:04d}"
        updated = {}

        for cond in ("dar", "contact", "double"):
            work = WORKDIR / f"{prefix}_{cond}"
            avg_files = list(work.glob("Average_*.fxout")) if work.exists() else []
            dif_files = list(work.glob("Dif_*.fxout"))     if work.exists() else []
            if avg_files:
                parsed = parse_buildmodel_ddg(avg_files[0])
                if parsed is not None:
                    ddg, sd = parsed
                    updated[f"ddg_stab_{cond}"]    = ddg
                    updated[f"ddg_stab_{cond}_sd"] = sd

        # Propagate epistasis SD
        sd_dar     = updated.get("ddg_stab_dar_sd",     df.at[i, "ddg_stab_dar_sd"]     if "ddg_stab_dar_sd"     in df.columns else np.nan)
        sd_contact = updated.get("ddg_stab_contact_sd", df.at[i, "ddg_stab_contact_sd"] if "ddg_stab_contact_sd" in df.columns else np.nan)
        sd_double  = updated.get("ddg_stab_double_sd",  df.at[i, "ddg_stab_double_sd"]  if "ddg_stab_double_sd"  in df.columns else np.nan)
        epi        = df.at[i, "ddg_epistasis"] if "ddg_epistasis" in df.columns else np.nan

        if all(pd.notna(v) for v in [sd_dar, sd_contact, sd_double]):
            epi_sd = float(np.sqrt(float(sd_dar)**2 + float(sd_contact)**2 + float(sd_double)**2))
            updated["ddg_epistasis_sd"]  = epi_sd
            updated["ddg_epistasis_sig"] = (epi_sd > 0 and pd.notna(epi) and abs(float(epi)) > 2 * epi_sd)

        for col, val in updated.items():
            df.at[i, col] = val

    df.to_csv(existing_csv, index=False)
    print(f"\nUpdated {existing_csv}")
    _print_summary(df)


def _print_summary(df: pd.DataFrame) -> None:
    """Print threshold summary at both 0.5 and 1.0 kcal/mol."""
    n_rescue_05 = (df["ddg_rescue_stab"] > 0.5).sum() if "ddg_rescue_stab" in df.columns else 0
    n_rescue_10 = (df["ddg_rescue_stab"] > 1.0).sum() if "ddg_rescue_stab" in df.columns else 0
    n_dar_05    = (df["ddg_stab_dar"] > 0.5).sum()    if "ddg_stab_dar" in df.columns else 0
    n_dar_10    = (df["ddg_stab_dar"] > 1.0).sum()    if "ddg_stab_dar" in df.columns else 0
    n_both_05   = ((df["ddg_rescue_stab"] > 0.5) & (df["ddg_stab_dar"] > 0.5)).sum() \
                  if all(c in df.columns for c in ["ddg_rescue_stab", "ddg_stab_dar"]) else 0
    n_both_10   = ((df["ddg_rescue_stab"] > 1.0) & (df["ddg_stab_dar"] > 1.0)).sum() \
                  if all(c in df.columns for c in ["ddg_rescue_stab", "ddg_stab_dar"]) else 0
    n_epi_sig   = df["ddg_epistasis_sig"].sum() if "ddg_epistasis_sig" in df.columns else 0

    print("\n── FoldX threshold summary ──────────────────────────────────")
    print(f"  ΔΔG_DAR > 0.5 kcal/mol (destabilising):       {n_dar_05}")
    print(f"  ΔΔG_DAR > 1.0 kcal/mol (strong destabilise):  {n_dar_10}")
    print(f"  ΔΔG_rescue > 0.5 kcal/mol:                    {n_rescue_05}  (within FoldX error, not threshold)")
    print(f"  ΔΔG_rescue > 1.0 kcal/mol (reported threshold): {n_rescue_10}")
    print(f"  Both rescue+DAR > 0.5:  {n_both_05}  |  Both > 1.0: {n_both_10}")
    print(f"  |ΔΔG_epistasis| > 2×SD (significant non-additive): {n_epi_sig}")

    if "ddg_epistasis_sd" in df.columns:
        epi_check = df[["ddg_stab_dar", "dar_gene", "dar_aa_coord", "dar_alt_aa",
                         "contact_gene", "contact_refseq_pos", "contact_alt_aa",
                         "ddg_epistasis", "ddg_epistasis_sd", "ddg_epistasis_sig"]].copy()
        epi_check = epi_check[epi_check["ddg_epistasis_sig"] == True].sort_values(
            "ddg_epistasis", key=abs, ascending=False)
        if len(epi_check):
            print("\n  Pairs with significant non-additive epistasis (|epi| > 2×SD):")
            for _, r in epi_check.iterrows():
                print(f"    {r.get('dar_gene','?')}:{r.get('dar_aa_coord','?')}:{r.get('dar_alt_aa','?')} <-> "
                      f"{r.get('contact_gene','?')}:{r.get('contact_refseq_pos','?')}:{r.get('contact_alt_aa','?')}"
                      f"  epistasis={r['ddg_epistasis']:.3f}  SD={r['ddg_epistasis_sd']:.3f}")
        else:
            print("\n  No pairs pass |epistasis| > 2×SD threshold.")


# ─── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(description="FoldX ΔΔG for compensatory pairs")
    parser.add_argument(
        "--reparse", action="store_true",
        help="Re-parse existing Average_*.fxout files to extract SD; no FoldX re-run"
    )
    parser.add_argument(
        "--input", type=str, default=None,
        help="Path to input pairs CSV (default: results/mutagenesis/top_targets.csv). "
             "Must be structurally anchored (have a pdb_id column) — e.g. "
             "results/mutagenesis/prioritized_pairs.csv for the full 694-pair scale-up.",
    )
    parser.add_argument(
        "--chunk-idx", type=int, default=0,
        help="0-based chunk index for SLURM array parallelization (default: 0)",
    )
    parser.add_argument(
        "--n-chunks", type=int, default=1,
        help="Total number of chunks (SLURM array size; default: 1 = no chunking)",
    )
    parser.add_argument(
        "--work-dir", type=Path, default=None,
        help="Repository root on HPC (data/ + results/ live here; code may live elsewhere)",
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip pairs (by physical dar/contact/alt-aa key) already present in "
             "results/mutagenesis/foldx_ddg.csv — avoids re-running pairs a prior "
             "local/array run already completed. No effect if that file doesn't exist.",
    )
    args = parser.parse_args()

    global ROOT, TOP_TARGETS_CSV, CONTACTS_CSV, STRUCTURES_DIR, WORKDIR, OUT_DIR, \
        _FOLDX_PATH, _FOLDX_DIR

    if args.work_dir is not None:
        ROOT = args.work_dir
        STRUCTURES_DIR = ROOT / "data" / "structures"
        OUT_DIR = ROOT / "results" / "mutagenesis"
        # keep foldx binary resolution relative to the *code* location (unchanged
        # via FOLDX_PATH / tools/foldx under the code checkout), not the work-dir
        OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.input is not None:
        TOP_TARGETS_CSV = Path(args.input)
    elif args.work_dir is not None:
        TOP_TARGETS_CSV = ROOT / "results" / "mutagenesis" / "top_targets.csv"

    # Per-chunk work directory keeps concurrent SLURM array tasks from
    # colliding on foldx_work/pairNNNN_* subdirectories.
    if args.n_chunks > 1:
        WORKDIR = OUT_DIR / "foldx_work" / f"chunk{args.chunk_idx:04d}"
    else:
        WORKDIR = OUT_DIR / "foldx_work"
    WORKDIR.mkdir(parents=True, exist_ok=True)

    if args.reparse:
        reparse_existing()
        return

    check_foldx()

    if not TOP_TARGETS_CSV.exists():
        sys.exit(
            f"\nERROR: {TOP_TARGETS_CSV} not found.\n"
            "  Run 00_prioritize_pairs.py first, or pass --input <path>."
        )

    print(f"Loading pairs from {TOP_TARGETS_CSV} ...")
    targets = pd.read_csv(TOP_TARGETS_CSV, dtype={"ann_id": str})
    print(f"  {len(targets)} pairs total")

    has_coords = targets["pdb_id"].notna().sum()
    print(f"  {has_coords} pairs have PDB structural coordinates")

    if has_coords == 0:
        print(
            "\nWARNING: No structural coordinates in input CSV.\n"
            "  Re-run 00_prioritize_pairs.py to ensure the contacts join was successful."
        )

    # Preserve original row index (pair identity + foldx_work/pairNNNN_* naming)
    # before chunking, so array tasks never reuse another task's directory names.
    targets = targets.reset_index(drop=False).rename(columns={"index": "orig_idx"})

    if args.skip_existing:
        existing_path = OUT_DIR / "foldx_ddg.csv"
        if existing_path.exists():
            existing = pd.read_csv(existing_path, dtype={"ann_id": str})
            phys_key = ["dar_gene", "dar_aa_coord", "dar_alt_aa",
                        "contact_gene", "contact_refseq_pos", "contact_alt_aa"]
            for col in phys_key:
                existing[col] = existing[col].astype(str)
                targets[col] = targets[col].astype(str)
            done_keys = set(existing[phys_key].agg("|".join, axis=1))
            before = len(targets)
            targets = targets[~targets[phys_key].agg("|".join, axis=1).isin(done_keys)].copy()
            print(f"  --skip-existing: {before - len(targets)} pairs already in foldx_ddg.csv, "
                  f"{len(targets)} remaining")
        else:
            print("  --skip-existing: no existing foldx_ddg.csv found, running all input pairs")

    if args.n_chunks > 1:
        targets = targets.iloc[args.chunk_idx::args.n_chunks].copy()
        print(f"  Chunk {args.chunk_idx}/{args.n_chunks}: {len(targets)} pairs assigned to this task")

    results = []
    for _, row in targets.iterrows():
        i = int(row["orig_idx"])
        print(
            f"  [orig_idx={i:4d}] "
            f"{row['dar_gene']}:{row['dar_aa_coord']}:{row['dar_alt_aa']} <-> "
            f"{row['contact_gene']}:{row['contact_refseq_pos']}:{row['contact_alt_aa']}",
            flush=True,
        )
        r = compute_ddg_for_pair(row, i)
        r["orig_idx"] = i
        print(
            f"    ΔΔG_stab_DAR={r['ddg_stab_dar']:.2f} "
            f"rescue_stab={r['ddg_rescue_stab']:.2f} "
            f"epistasis={r['ddg_epistasis']:.2f}"
            if all(pd.notna(r.get(k)) for k in ["ddg_stab_dar", "ddg_rescue_stab", "ddg_epistasis"])
            else f"    status={r['foldx_status']}"
        )
        results.append(r)

    out_df = pd.DataFrame(results)

    if args.n_chunks > 1:
        out_path = OUT_DIR / f"foldx_ddg_chunk{args.chunk_idx:04d}.csv"
        out_df.to_csv(out_path, index=False)
        print(f"\nSaved chunk → {out_path}")
        print("  Merge all chunks with: python src/mutagenesis/04_merge_foldx_chunks.py")
        return

    out_path = OUT_DIR / "foldx_ddg.csv"
    out_df.to_csv(out_path, index=False)
    print(f"\nSaved → {out_path}")
    _print_summary(out_df)

    print("\nStatus breakdown:")
    print(out_df["foldx_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
