import subprocess
from pathlib import Path

# ==== Configuration ====
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"

# Input
RAW_MT_DIR = DATA_DIR / "alignments" / "mtdna_raw_cds"

# Outputs
CODON_OUT_DIR = DATA_DIR / "alignments" / "mtdna_codon"
AA_OUT_DIR = DATA_DIR / "alignments" / "mtdna_aa"

CODON_OUT_DIR.mkdir(parents=True, exist_ok=True)
AA_OUT_DIR.mkdir(parents=True, exist_ok=True)


def align_gene(raw_fasta: Path):
    """Runs native MACSE with the mammalian mitochondrial genetic code."""
    gene_name = raw_fasta.name.replace("_raw_cds.fasta", "")

    out_nt = CODON_OUT_DIR / f"{gene_name}_codon_alignment.fasta"
    out_aa = AA_OUT_DIR / f"{gene_name}_aa_alignment.fasta"

    if out_nt.exists() and out_aa.exists() and out_nt.stat().st_size > 0:
        print(f"  [SKIP] {gene_name} (Already aligned)")
        return True

    print(f"  [ALIGN] Processing {gene_name} (This may take several minutes)...")

    # Construct native MACSE command
    # -gc_def 2 : Vertebrate Mitochondrial Translation Table
    cmd = [
        "macse",
        "-prog",
        "alignSequences",
        "-gc_def",
        "2",
        "-seq",
        str(raw_fasta),
        "-out_NT",
        str(out_nt),
        "-out_AA",
        str(out_aa),
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  [FAIL] MACSE failed on {gene_name}:")
        print(e.stderr[:500] + "...\n")
        return False


def main():
    print("Initializing Native mtDNA Codon Alignment Pipeline...\n")

    raw_files = sorted(list(RAW_MT_DIR.glob("*_raw_cds.fasta")))
    if not raw_files:
        print(f"No raw CDS files found in {RAW_MT_DIR.relative_to(ROOT)}.")
        return

    print(f"Found {len(raw_files)} mitochondrial genes to align.\n")

    success_count = 0
    for raw_fasta in raw_files:
        if align_gene(raw_fasta):
            success_count += 1

    print(f"\n{'='*50}")
    print("ALIGNMENT SUMMARY")
    print(f"{'='*50}")
    print(f"Successfully aligned: {success_count} / {len(raw_files)} genes.")


if __name__ == "__main__":
    main()
