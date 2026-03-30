"""
Alignment File Integrity Check — mtDNA and nucDNA
==================================================
Verifies structural integrity of all FASTA alignment files across all four
alignment directories:

  toga_hg38_aa     — nucDNA amino acid alignments
  toga_hg38_codon  — nucDNA codon (nucleotide) alignments
  mtdna_aa         — mtDNA amino acid alignments
  mtdna_codon      — mtDNA codon (nucleotide) alignments

Checks per file:
  1. Human reference sequence present (header starts with "Homo_sapiens")
  2. Species count >= MIN_SPECIES
  3. Codon files: sequence length divisible by 3
  4. Nucleotide alphabet: only {A,C,T,G,N,X,-,!,*} — no amino acid letters
  5. AA alphabet: only standard amino acids + mask chars {-,X,!,*}
  6. Paired codon/AA files exist for the same gene and have matching species sets

Run from project root:
    python tests/test_check_alignment.py
"""

import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DATA_DIR    = ROOT / "data" / "alignments"
TOGA_AA_DIR = DATA_DIR / "toga_hg38_aa"
TOGA_NT_DIR = DATA_DIR / "toga_hg38_codon"
MT_AA_DIR   = DATA_DIR / "mtdna_aa"
MT_NT_DIR   = DATA_DIR / "mtdna_codon"

MIN_SPECIES = 100  # minimum non-human species expected in any alignment

VALID_NT = set("ACTGNX-!*RYSWKMBDHV")  # includes IUPAC ambiguity codes
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-X!*")


# ==== FASTA parser (no BioPython dependency) ====

def parse_fasta(path: Path) -> dict:
    """Returns {header_id: sequence} — sequence uppercased, gaps preserved."""
    seqs = {}
    current = None
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]  # id = everything before first space
                seqs[current] = []
            elif current is not None:
                seqs[current].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def check_file(path: Path, is_codon: bool) -> list:
    """Returns list of error strings (empty = all checks passed)."""
    errors = []
    seqs = parse_fasta(path)

    if not seqs:
        return ["EMPTY: no sequences parsed"]

    # 1. Human reference present — prefer Homo_sapiens|9606 (modern human, taxon ID)
    # to avoid selecting archaic humans (neanderthalensis, Denisova) that also start
    # with "Homo_sapiens" but have different sequences from the reference genome.
    human_keys = [k for k in seqs if k.startswith("Homo_sapiens|9606")]
    if not human_keys:
        human_keys = [k for k in seqs if k.startswith("Homo_sapiens")]
    if not human_keys:
        errors.append("MISSING human reference (no 'Homo_sapiens' header)")
        return errors  # remaining checks need ref sequence

    ref_seq  = seqs[human_keys[0]]
    n_species = len(seqs) - 1  # exclude human

    # 2. Species count
    if n_species < MIN_SPECIES:
        errors.append(f"LOW SPECIES COUNT: {n_species} (expected >= {MIN_SPECIES})")

    # 3. Codon architecture (check ungapped human reference)
    if is_codon:
        ungapped = ref_seq.replace("-", "")
        if len(ungapped) % 3 != 0:
            errors.append(
                f"NOT CODON-AWARE: ungapped length {len(ungapped)} not divisible by 3"
            )

    # 4 / 5. Alphabet check across all sequences
    valid = VALID_NT if is_codon else VALID_AA
    rogue = set()
    for seq in seqs.values():
        rogue |= set(seq) - valid
    if rogue:
        errors.append(f"UNEXPECTED CHARACTERS: {sorted(rogue)}")

    return errors


def check_pairing(aa_dir: Path, nt_dir: Path, label: str) -> list:
    """Returns error strings for missing paired files or species set mismatches."""
    errors = []
    aa_genes = {f.stem.replace("_aa_alignment",    ""): f for f in aa_dir.glob("*_aa_alignment.fasta")}
    nt_genes = {f.stem.replace("_codon_alignment", ""): f for f in nt_dir.glob("*_codon_alignment.fasta")}

    for gene in sorted(set(aa_genes) - set(nt_genes)):
        errors.append(f"{gene}: AA file exists but codon file missing")
    for gene in sorted(set(nt_genes) - set(aa_genes)):
        errors.append(f"{gene}: codon file exists but AA file missing")

    for gene in sorted(set(aa_genes) & set(nt_genes)):
        aa_ids = set(parse_fasta(aa_genes[gene]).keys())
        nt_ids = set(parse_fasta(nt_genes[gene]).keys())
        diff   = aa_ids.symmetric_difference(nt_ids)
        if diff:
            errors.append(
                f"{gene}: AA/codon species mismatch "
                f"({len(aa_ids)} vs {len(nt_ids)} sequences, {len(diff)} differing IDs)"
            )

    return errors


def check_directory(aa_dir: Path, nt_dir: Path, label: str):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")

    aa_files = sorted(aa_dir.glob("*_aa_alignment.fasta"))
    nt_files = sorted(nt_dir.glob("*_codon_alignment.fasta"))

    if not aa_files and not nt_files:
        print("  No alignment files found.")
        return

    print(f"  AA files   : {len(aa_files)}  ({aa_dir.name})")
    print(f"  Codon files: {len(nt_files)}  ({nt_dir.name})")

    total_errors = []

    for path in aa_files:
        for e in check_file(path, is_codon=False):
            total_errors.append(f"[AA]    {path.stem}: {e}")

    for path in nt_files:
        for e in check_file(path, is_codon=True):
            total_errors.append(f"[Codon] {path.stem}: {e}")

    for e in check_pairing(aa_dir, nt_dir, label):
        total_errors.append(f"[Pair]  {e}")

    if total_errors:
        print(f"\n  {len(total_errors)} issue(s) found:")
        for e in total_errors:
            print(f"    ✗ {e}")
    else:
        print(f"\n  All {len(aa_files) + len(nt_files)} files passed.")

    # Special character inventory across all non-human sequences
    print(f"\n  Special character inventory (non-human sequences):")
    for dir_path, is_codon, tag in [(aa_dir, False, "AA"), (nt_dir, True, "Codon")]:
        chars: Counter = Counter()
        total_bases = 0
        for path in dir_path.glob("*.fasta"):
            seqs = parse_fasta(path)
            human = next((k for k in seqs if k.startswith("Homo_sapiens|9606")), None) or \
                    next((k for k in seqs if k.startswith("Homo_sapiens")), None)
            for k, seq in seqs.items():
                if k == human:
                    continue
                for c in seq:
                    chars[c] += 1
                    total_bases += 1
        special = {c: n for c, n in chars.items()
                   if c not in "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy"}
        if special:
            print(f"    {tag} ({total_bases:,} total chars):")
            for c, n in sorted(special.items(), key=lambda x: -x[1]):
                print(f"      {repr(c)}: {n:>8,}  ({100*n/total_bases:.3f}%)")
        else:
            print(f"    {tag}: no special characters found")


def main():
    print("Alignment File Integrity Check")
    print("=" * 60)

    check_directory(TOGA_AA_DIR, TOGA_NT_DIR, "nucDNA  (TOGA hg38)")
    check_directory(MT_AA_DIR,   MT_NT_DIR,   "mtDNA")

    print(f"\n{'='*60}")
    print("Done.")


if __name__ == "__main__":
    main()
