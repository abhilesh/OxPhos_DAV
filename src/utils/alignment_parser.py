from pathlib import Path
from Bio import SeqIO
from collections import defaultdict


class AlignmentParser:
    def __init__(self, aa_fasta: Path, nt_fasta: Path, genome: str):
        self.genome = genome
        self.aa_alignment = self._load_fasta(aa_fasta)
        self.nt_alignment = self._load_fasta(nt_fasta)

        self.ref_header = self._identify_reference()
        self.transcript_id = (
            self.ref_header.split("|")[2]
            if len(self.ref_header.split("|")) > 2
            else "UNKNOWN"
        )

        # Maps biological positions to alignment column indices
        self.aa_map = self._build_coordinate_map(self.aa_alignment[self.ref_header])
        self.nt_map = self._build_coordinate_map(self.nt_alignment[self.ref_header])

    def _load_fasta(self, file_path: Path) -> dict:
        seqs = {}
        for record in SeqIO.parse(file_path, "fasta"):
            # Ensure we keep the full header exactly as mapped in 03_sanitize_all_alignments.py
            seqs[record.id] = str(record.seq).upper()
        return seqs

    def _identify_reference(self) -> str:
        # Match Homo_sapiens|9606| (taxon ID 9606 = modern human) to avoid
        # accidentally selecting Homo_sapiens_neanderthalensis or other archaic humans
        # that also start with "Homo_sapiens" but have different sequences.
        for header in self.aa_alignment.keys():
            if header.startswith("Homo_sapiens|9606"):
                return header
        # Fallback: any Homo_sapiens entry (e.g. TOGA alignments without taxon ID)
        for header in self.aa_alignment.keys():
            if header.startswith("Homo_sapiens"):
                return header
        raise ValueError("Human reference not found in alignment.")

    def _build_coordinate_map(self, ref_seq: str) -> dict:
        coord_map = {}
        bio_pos = 1
        for col_idx, char in enumerate(ref_seq):
            if char != "-":
                coord_map[bio_pos] = col_idx
                bio_pos += 1
        return coord_map

    def find_sequence_anchor(self, aa_pos: int, expected_ref_aa: str) -> int:
        """
        Dynamically corrects transcript/isoform shifts by searching the
        alignment for the expected amino acid context.
        """
        ref_seq = self.aa_alignment[self.ref_header]

        # 1. Try strict coordinate first
        if aa_pos in self.aa_map:
            col = self.aa_map[aa_pos]
            if ref_seq[col] == expected_ref_aa:
                return aa_pos

        # 2. Sliding window search for isoform offsets
        search_window = 10
        start_search = max(1, aa_pos - search_window)
        end_search = aa_pos + search_window

        for i in range(start_search, end_search + 1):
            if i in self.aa_map:
                col = self.aa_map[i]
                if ref_seq[col] == expected_ref_aa:
                    return i

        return None  # Unrecoverable mismatch

    def extract_mutant_codon(self, nt_pos: int, alt_nt: str) -> str:
        """
        Constructs the expected mutant codon by injecting the alt nucleotide
        into the correct reading frame of the reference sequence.
        """
        if nt_pos not in self.nt_map:
            return None

        # Snap to the correct codon reading frame (1-indexed)
        frame_offset = (nt_pos - 1) % 3
        codon_start_pos = nt_pos - frame_offset

        # Verify all 3 bases of the codon exist in the map
        for i in range(3):
            if (codon_start_pos + i) not in self.nt_map:
                return None

        ref_seq = self.nt_alignment[self.ref_header]
        codon_bases = [
            ref_seq[self.nt_map[codon_start_pos]],
            ref_seq[self.nt_map[codon_start_pos + 1]],
            ref_seq[self.nt_map[codon_start_pos + 2]],
        ]

        # Inject the mutant base at the exact offset
        codon_bases[frame_offset] = alt_nt.upper()

        return "".join(codon_bases)

    def check_compensation(
        self, reported_aa_pos: int, wt_aa: str, mut_aa: str, nt_pos: int, mut_codon: str
    ) -> dict:
        """Evaluates c-DARs using sequence-anchored coordinates.

        wt_aa: the wild-type amino acid from the variant record (used as anchor).
               This is the authoritative source — do NOT derive it from the alignment,
               because MACSE alignments may start partway into the CDS (coordinate offset).
        """
        results = {
            "aa_cdar": False,
            "nt_cdar": False,
            "aa_species": [],
            "nt_species": [],
        }

        true_aa_pos = self.find_sequence_anchor(reported_aa_pos, wt_aa)

        if not true_aa_pos:
            return results

        aa_col = self.aa_map[true_aa_pos]

        # If the anchor shifted the AA position, shift the NT position by the same amount
        # (MACSE alignments may start partway into the CDS, creating a constant offset)
        aa_shift = true_aa_pos - reported_aa_pos
        corrected_nt_pos = nt_pos + aa_shift * 3

        # Snap nucleotide positions to the correct reading frame
        frame_offset = (corrected_nt_pos - 1) % 3
        codon_start_pos = corrected_nt_pos - frame_offset

        nt_cols = [
            self.nt_map.get(codon_start_pos),
            self.nt_map.get(codon_start_pos + 1),
            self.nt_map.get(codon_start_pos + 2),
        ]

        if None in nt_cols:
            return results

        for species, aa_seq in self.aa_alignment.items():
            if species == self.ref_header:
                continue

            # Skip gap, frameshift (!), stop (*), and uncertain (X) positions
            _MASK = {"-", "!", "*", "X"}
            if aa_seq[aa_col] in _MASK or any(
                self.nt_alignment[species][c] in _MASK for c in nt_cols
            ):
                continue

            # AA Level Check
            if aa_seq[aa_col] == mut_aa:
                results["aa_cdar"] = True
                sp_name = species.split("|")[0]
                results["aa_species"].append(sp_name)

                # NT Level Check (Strict subset of AA check)
                species_codon = "".join(
                    [self.nt_alignment[species][c] for c in nt_cols]
                )
                if species_codon == mut_codon:
                    results["nt_cdar"] = True
                    results["nt_species"].append(sp_name)

        return results
