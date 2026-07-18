"""
src/phylo/00a_download_vertlife_tree.py

Obtains a time-calibrated ultrametric species tree for the 312-species
cross-genome overlap (species with data in BOTH nucDNA TOGA and mtDNA
alignments) and prunes/saves it for use in IQTree jobs and Pagel's test.

Tree source strategy (in order):
  1. GitHub direct download — MamPhy v1 MCC tree (Upham et al. 2019), 5,911
     mammal species, Nexus format. Downloaded from the MamPhy_v1 GitHub repo
     and converted to Newick. No authentication required.
  2. TimeTree API — queries timetree.org/api with the species list to get a
     dated Newick tree directly for our cross-genome species. No intermediate
     full-tree download needed.
  3. VertLife phylosubsets (manual) — if both automated downloads fail,
     instructions are printed for manually downloading a pruned tree from
     vertlife.org/phylosubsets/.
  4. Local file — if mammaltree_full.nwk already exists (e.g. manually placed),
     prune it to the cross-genome species.

Outputs:
  data/phylo/species_tree/mammaltree_crossgenome.nwk  -- tree for cross-genome species
  data/phylo/cross_genome_species.txt                 -- species list retained in tree

Run from project root inside the Docker container:
    python src/phylo/00a_download_vertlife_tree.py
"""

import csv
import json
import ssl
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path

from Bio import Phylo

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT       = Path(__file__).resolve().parents[2]
DATA_DIR   = ROOT / "data"
REF_DIR    = DATA_DIR / "reference"
PHYLO_DIR  = DATA_DIR / "phylo" / "species_tree"
TAX_MAP    = REF_DIR / "taxid_species_mapping.csv"

FULL_TREE       = PHYLO_DIR / "mammaltree_full.nwk"      # optional: manually placed VertLife tree
FULL_TREE_NEXUS = PHYLO_DIR / "mammaltree_full.tre"      # MamPhy Nexus download
CG_TREE         = PHYLO_DIR / "mammaltree_crossgenome.nwk"
CROSS_SPP       = DATA_DIR / "phylo" / "cross_genome_species.txt"

MAMPHY_URL   = ("https://raw.githubusercontent.com/n8upham/MamPhy_v1/master/_DATA/"
                "MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")
TIMETREE_API = "https://timetree.org/api/tree"

VERTLIFE_INSTRUCTIONS = """
─────────────────────────────────────────────────────
VertLife manual download instructions (alternative):
  1. Go to https://vertlife.org/phylosubsets/
  2. Upload the species list: data/phylo/cross_genome_species.txt
     (written by this script before the download attempt)
  3. Download 1 tree from the DNA-only distribution (NDexp model)
  4. Save it as: data/phylo/species_tree/mammaltree_full.nwk
  5. Re-run this script — it will detect the file and prune it.
─────────────────────────────────────────────────────
"""

_ctx = ssl.create_default_context()
_ctx.check_hostname = False
_ctx.verify_mode = ssl.CERT_NONE


# ── Cross-genome species ───────────────────────────────────────────────────────

def load_cross_genome_species() -> list[str]:
    """
    Load cross-genome overlap from taxid_species_mapping.csv.
    Uses toga_species name (canonical TOGA identifier) from Exact_TaxID_Match rows.
    TaxID-based matching ensures correctness; species names may differ slightly
    between TOGA and mtDNA headers for the same organism.
    """
    if not TAX_MAP.exists():
        print(f"ERROR: {TAX_MAP} not found — run src/data_prep/00d_taxid_species_map.py first.")
        sys.exit(1)
    overlap = sorted({
        row["toga_species"]
        for row in csv.DictReader(open(TAX_MAP))
        if row.get("match_type") == "Exact_TaxID_Match" and row.get("toga_species")
    })
    print(f"Cross-genome species (TaxID-matched, TOGA names): {len(overlap)}")
    return overlap


# ── TimeTree API ───────────────────────────────────────────────────────────────

def fetch_timetree(species: list[str], out_path: Path) -> bool:
    """
    Query TimeTree API with a species list, get a dated Newick tree.
    Returns True on success.

    TimeTree accepts POST to /api/tree with JSON body:
      {"taxa": ["Homo sapiens", "Mus musculus", ...]}
    Returns JSON with "newick" field.
    """
    print(f"Querying TimeTree API for {len(species)} species ...")
    # TimeTree uses space-separated binomial names (underscores → spaces)
    taxa = [sp.replace("_", " ") for sp in species]
    payload = json.dumps({"taxa": taxa}).encode("utf-8")
    req = urllib.request.Request(
        TIMETREE_API,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "OxPhos_DAV/1.0",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, context=_ctx, timeout=300) as r:
            data = json.loads(r.read())
        newick = data.get("newick") or data.get("tree")
        if not newick:
            print(f"  TimeTree response missing 'newick' field. Keys: {list(data.keys())}")
            return False
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(newick)
        n_tips = newick.count(",") + 1
        print(f"  TimeTree tree (~{n_tips} tips) → {out_path}")
        return True
    except Exception as e:
        print(f"  TimeTree API failed: {e}")
        return False


# ── MamPhy GitHub download ────────────────────────────────────────────────────

def download_mamphy_github(out_path: Path) -> bool:
    """
    Download the MamPhy v1 MCC tree (Nexus) directly from GitHub.
    Returns True on success, False on any error.
    """
    print(f"Downloading MamPhy MCC tree from GitHub ...")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        urllib.request.urlretrieve(MAMPHY_URL, str(out_path))
        size_mb = out_path.stat().st_size / 1e6
        print(f"  Downloaded {size_mb:.1f} MB → {out_path}")
        return True
    except Exception as e:
        print(f"  GitHub download failed: {e}")
        return False


def nexus_to_newick(nexus_path: Path, newick_path: Path) -> None:
    """Convert a Nexus tree file to Newick format using BioPython Phylo.

    MamPhy trees contain BEAST-style inline annotations on every branch
    (e.g. [&height_95%_HPD={...}]). These are stored as clade.comment by
    BioPython and written into the Newick output, causing IQ-TREE to reject
    the file. Strip all comments and confidence values before writing.
    """
    print(f"Converting Nexus → Newick ...")
    tree = Phylo.read(str(nexus_path), "nexus")
    # Strip BEAST annotations from every clade
    for clade in tree.find_clades():
        clade.comment = None
        if not clade.is_terminal():
            clade.confidence = None   # remove internal node labels/posteriors
    newick_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(newick_path), "newick")
    print(f"  Newick tree → {newick_path}")


# ── VertLife pruning (if full tree provided manually) ─────────────────────────

def _tip_to_binomial(tip_name: str) -> str:
    """
    Strip trailing taxonomy suffixes from MamPhy tip names.
    MamPhy format: 'Genus_species_FAMILY_ORDER' → 'Genus_species'
    Plain format:  'Genus_species' → 'Genus_species' (unchanged)
    """
    parts = tip_name.split("_")
    # MamPhy tips have ≥4 underscore-delimited parts; last two are FAMILY and ORDER
    # (both are all-uppercase).  Strip trailing all-caps tokens.
    while len(parts) > 2 and parts[-1].isupper():
        parts.pop()
    return "_".join(parts)


def prune_vertlife(full_tree_path: Path, keep: set[str], out_path: Path) -> set[str]:
    """Prune the VertLife full tree to `keep` species. Returns retained set.

    MamPhy tip names include '_FAMILY_ORDER' suffixes (e.g.
    'Acinonyx_jubatus_FELIDAE_CARNIVORA'); these are stripped before matching
    against the cross-genome species set (which uses plain 'Genus_species' names).
    Homo_sapiens is always retained as the outgroup/reference regardless of
    whether it is in the cross-genome overlap list.
    The pruned tree retains original full tip names for downstream compatibility.
    """
    print(f"\nLoading VertLife tree from {full_tree_path} ...")
    tree = Phylo.read(str(full_tree_path), "newick")
    terminals = tree.get_terminals()
    print(f"  {len(terminals)} tips in full tree.")

    # Build map: binomial → full tip name (original)
    binomial_to_full = {_tip_to_binomial(c.name): c.name for c in terminals}

    # Always keep Homo_sapiens as reference
    keep_with_human = keep | {"Homo_sapiens"}

    retain_full = set()   # full tip names to keep
    missing     = set()   # keep-species not found in tree
    for sp in keep_with_human:
        if sp in binomial_to_full:
            retain_full.add(binomial_to_full[sp])
        else:
            missing.add(sp)

    print(f"  Cross-genome species in tree : {len(retain_full)}")
    if missing:
        print(f"  Not in VertLife tree (dropped): {len(missing)}")
        for sp in sorted(missing)[:10]:
            print(f"    {sp}")
        if len(missing) > 10:
            print(f"    ... and {len(missing) - 10} more")

    all_full = {c.name for c in tree.get_terminals()}
    print(f"Pruning to {len(retain_full)} species (including Homo_sapiens) ...")
    for name in all_full - retain_full:
        tree.prune(name)

    # Rename tips back to plain binomial names for downstream use
    for clade in tree.get_terminals():
        clade.name = _tip_to_binomial(clade.name)

    remaining = {c.name for c in tree.get_terminals()}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(out_path), "newick")
    print(f"  Pruned tree → {out_path}  ({len(remaining)} tips)")
    return remaining


# ── Validate tree ──────────────────────────────────────────────────────────────

def validate_and_report(tree_path: Path, expected_spp: set[str]) -> set[str]:
    """Read the written tree, report tip coverage vs expected species."""
    tree = Phylo.read(str(tree_path), "newick")
    tips = {c.name for c in tree.get_terminals()}
    # TimeTree may return species with spaces — normalize
    tips_norm = {t.replace(" ", "_") for t in tips}
    matched = expected_spp & tips_norm
    missing = expected_spp - tips_norm
    print(f"\nTree validation:")
    print(f"  Expected cross-genome species : {len(expected_spp)}")
    print(f"  Found in tree                 : {len(matched)}")
    if missing:
        print(f"  Missing from tree             : {len(missing)}")
        for sp in sorted(missing)[:5]:
            print(f"    {sp}")
    return matched


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    PHYLO_DIR.mkdir(parents=True, exist_ok=True)
    (DATA_DIR / "phylo").mkdir(parents=True, exist_ok=True)

    # ── Step 1: load cross-genome species and write list ──────────────────────
    cross_spp = load_cross_genome_species()
    CROSS_SPP.parent.mkdir(parents=True, exist_ok=True)
    CROSS_SPP.write_text("\n".join(cross_spp) + "\n")
    print(f"Species list → {CROSS_SPP}")

    cross_spp_set = set(cross_spp)

    # ── Step 2: obtain tree ────────────────────────────────────────────────────
    if CG_TREE.exists():
        print(f"\nCross-genome tree already present: {CG_TREE}")
        retained = validate_and_report(CG_TREE, cross_spp_set)

    elif FULL_TREE.exists():
        # VertLife full Newick tree placed manually
        print(f"\nVertLife full tree found: {FULL_TREE}")
        retained = prune_vertlife(FULL_TREE, cross_spp_set, CG_TREE)

    elif FULL_TREE_NEXUS.exists():
        # MamPhy Nexus already downloaded — convert and prune
        print(f"\nMamPhy Nexus tree found: {FULL_TREE_NEXUS}")
        nexus_to_newick(FULL_TREE_NEXUS, FULL_TREE)
        retained = prune_vertlife(FULL_TREE, cross_spp_set, CG_TREE)

    else:
        # 1. Try GitHub direct download (MamPhy v1)
        ok = download_mamphy_github(FULL_TREE_NEXUS)
        if ok:
            nexus_to_newick(FULL_TREE_NEXUS, FULL_TREE)
            retained = prune_vertlife(FULL_TREE, cross_spp_set, CG_TREE)
        else:
            # 2. Fall back to TimeTree API
            ok = fetch_timetree(cross_spp, CG_TREE)
            if ok:
                retained = validate_and_report(CG_TREE, cross_spp_set)
            else:
                # 3. Print manual download instructions
                print(VERTLIFE_INSTRUCTIONS)
                print(f"Species list for upload → {CROSS_SPP}")
                sys.exit(1)

    # ── Step 3: update species list to exactly those in the tree ─────────────
    # Normalize tip names (TimeTree may use spaces)
    tree = Phylo.read(str(CG_TREE), "newick")
    final_tips = {c.name.replace(" ", "_") for c in tree.get_terminals()}
    CROSS_SPP.write_text("\n".join(sorted(final_tips)) + "\n")
    print(f"\nFinal cross-genome species list → {CROSS_SPP}  ({len(final_tips)} species)")
    print("\nDone. Next step:")
    print("  python src/phylo/00_prep_iqtree_jobs.py")


if __name__ == "__main__":
    main()