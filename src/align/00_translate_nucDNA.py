from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ==== Configuration ====
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"

# Input
NUC_CODON_DIR = DATA_DIR / "alignments" / "toga_hg38_codon"

# Output
NUC_AA_DIR = DATA_DIR / "alignments" / "toga_hg38_aa"
NUC_AA_DIR.mkdir(parents=True, exist_ok=True)


def translate_codon_sequence(nt_seq: str) -> str:
    """
    Translates a codon-aware nucleotide string into an amino acid string.
    Strictly handles MACSE/TOGA gap padding and Senckenberg missing data masks.
    """
    aa_seq = []

    # Process in chunks of 3 (codons)
    for i in range(0, len(nt_seq), 3):
        codon = nt_seq[i : i + 3]

        # Handle structural gaps
        if codon == "---":
            aa_seq.append("-")
            continue

        # Handle Senckenberg masking (N or X indicates damaged/missing sequence)
        if "N" in codon or "X" in codon or "-" in codon:
            aa_seq.append("X")
            continue

        # Standard Nuclear Translation (Table 1)
        try:
            # Biopython translate automatically handles standard codons
            aa = str(Seq(codon).translate(table=1))
            aa_seq.append(aa)
        except Exception:
            aa_seq.append("X")  # Catch-all for anomalous triplets

    return "".join(aa_seq)


def process_alignment(codon_fasta: Path) -> bool:
    """Reads a TOGA codon alignment and writes the corresponding AA alignment."""
    gene_name = codon_fasta.name.replace("_codon_alignment.fasta", "")
    out_file = NUC_AA_DIR / f"{gene_name}_aa_alignment.fasta"

    if out_file.exists() and out_file.stat().st_size > 0:
        return True

    try:
        records = list(SeqIO.parse(codon_fasta, "fasta"))
        aa_records = []

        for record in records:
            nt_string = str(record.seq).upper()

            # Validation: Must be codon-aware
            if len(nt_string) % 3 != 0:
                print(f"  [ERROR] {gene_name} sequence length is not divisible by 3.")
                return False

            aa_string = translate_codon_sequence(nt_string)

            aa_record = SeqRecord(Seq(aa_string), id=record.id, description="")
            aa_records.append(aa_record)

        SeqIO.write(aa_records, out_file, "fasta")
        print(f"  [OK] Translated {gene_name} ({len(records)} species)")
        return True

    except Exception as e:
        print(f"  [FAIL] Error translating {gene_name}: {e}")
        return False


def main():
    print("Initializing Nuclear DNA Translation Pipeline...\n")

    codon_files = sorted(list(NUC_CODON_DIR.glob("*_codon_alignment.fasta")))
    if not codon_files:
        print(f"No TOGA alignments found in {NUC_CODON_DIR.relative_to(ROOT)}.")
        return

    print(f"Found {len(codon_files)} nuclear alignments to translate.\n")

    success_count = 0
    for fasta in codon_files:
        if process_alignment(fasta):
            success_count += 1

    print(f"\n{'='*50}")
    print("TRANSLATION SUMMARY")
    print(f"{'='*50}")
    print(f"Successfully translated: {success_count} / {len(codon_files)} genes.")
    print(f"AA Alignments saved to: {NUC_AA_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
