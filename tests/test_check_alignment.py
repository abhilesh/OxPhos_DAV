from pathlib import Path
from collections import Counter


def verify_alignment():
    # Target the alignments directory
    align_dir = Path("data/alignments/toga_hg38_codon")

    # Grab the first available FASTA file (e.g., SDHA_codon_alignment.fasta)
    fasta_files = list(align_dir.glob("*_codon_alignment.fasta"))

    if not fasta_files:
        print("No FASTA files found in the directory.")
        return

    sample_file = fasta_files[0]
    print(f"{'='*50}")
    print(f"INSPECTING: {sample_file.name}")
    print(f"{'='*50}")

    num_species = 0
    ref_seq = ""
    temp_seq = []

    # Parse the FASTA file
    with open(sample_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                num_species += 1
                # Save the first sequence (the Human REFERENCE) for alphabet testing
                if num_species == 2:
                    ref_seq = "".join(temp_seq)
            elif num_species == 1:
                temp_seq.append(line.upper())

    # Catch the sequence if there was only one (unlikely in TOGA)
    if not ref_seq and temp_seq:
        ref_seq = "".join(temp_seq)

    # 1. Check Species Count
    print(f"Total Species (Sequences): {num_species}")

    # 2. Check Codon Architecture
    seq_length = len(ref_seq)
    is_codon_aware = seq_length % 3 == 0
    print(f"\nSequence Length: {seq_length} bp")
    print(f"Divisible by 3 (Codon-Aware)?: {is_codon_aware}")

    # 3. Check Alphabet (Nucleotide vs Protein)
    char_counts = Counter(ref_seq)
    print("\nCharacter Distribution in Human Reference:")
    for char, count in char_counts.most_common():
        print(f"  '{char}': {count}")

    # Define strict nucleotide characters (including Senckenberg's 'N' and 'X' masks, and '-' for gaps)
    allowed_nucleotides = set("ACTGNX-")
    found_chars = set(char_counts.keys())
    rogue_chars = found_chars - allowed_nucleotides

    print(f"\n{'='*50}")
    if rogue_chars:
        print(f"WARNING: Found unexpected characters: {rogue_chars}")
        print("This may be an amino acid sequence or contain unmasked artifacts.")
    else:
        print("VERIFIED: Strict Nucleotide (DNA) Alphabet Confirmed.")
    print(f"{'='*50}")


if __name__ == "__main__":
    verify_alignment()
