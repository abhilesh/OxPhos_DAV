import re
from pathlib import Path
from Bio import SeqIO


class AlignmentParser:
    def __init__(self, aa_fasta: Path, nt_fasta: Path, genome: str):
        self.genome = genome
        self.transcript_id = None
        self.aa_alignment = self._load_fasta(aa_fasta)
        self.nt_alignment = self._load_fasta(nt_fasta)

        # Identify the human reference sequence
        self.ref_header = self._identify_reference()

        # Build coordinate maps (Biological 1-indexed pos -> Alignment 0-indexed col)
        if self.ref_header:
            self.aa_map = self._build_coordinate_map(self.aa_alignment)
            self.nt_map = self._build_coordinate_map(self.nt_alignment)
        else:
            self.aa_map = {}
            self.nt_map = {}

    def _load_alignment(self, fasta_path: Path) -> dict:
        alignment = {}
        with open(fasta_path, "r") as f:
            header = None
            seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        alignment[header] = "".join(seq)

                    # Parse the header. Example: >Homo_sapiens|hg38|ENST00000263623
                    full_id = line[1:].split()[0]
                    parts = full_id.split("|")

                    # Keep the species name as the primary dictionary key
                    header = parts[0]

                    # Extract transcript ID from the human reference sequence
                    if header == "Homo_sapiens" and len(parts) >= 3:
                        self.transcript_id = parts[2]

                    seq = []
                else:
                    seq.append(line)
            if header:
                alignment[header] = "".join(seq)
        return alignment

    def _load_fasta(self, fasta_path: Path) -> dict:
        if not fasta_path.exists():
            return {}
        return {
            record.id: str(record.seq).upper()
            for record in SeqIO.parse(fasta_path, "fasta")
        }

    def _identify_reference(self) -> str:
        """Finds the human anchor sequence depending on the database format."""
        if not self.aa_alignment:
            return None

        # Look for the standardized human header we built
        human_header = next(
            (
                h
                for h in self.aa_alignment.keys()
                if "Homo_sapiens" in h or "REFERENCE" in h
            ),
            None,
        )
        if human_header:
            return human_header

        # Fallback for un-standardized files
        return list(self.aa_alignment.keys())[0]

    def _build_coordinate_map(self, alignment: dict) -> dict:
        """Maps biological position (ignoring gaps) to the alignment column matrix."""
        coord_map = {}
        biological_pos = 1
        ref_seq = alignment.get(self.ref_header, "")

        for col_index, char in enumerate(ref_seq):
            if char != "-":
                coord_map[biological_pos] = col_index
                biological_pos += 1

        return coord_map

    def extract_mutant_codon(self, nt_pos: int, alt_allele: str) -> str:
        """Determines the exact 3-letter mutant codon created by the variant."""
        if nt_pos not in self.nt_map:
            return None

        col_idx = self.nt_map[nt_pos]

        # Calculate codon boundaries (0, 1, or 2 positions back from the mutation)
        # Because we built the map ignoring gaps, we have to look at the raw biological reading frame
        codon_start_bio = nt_pos - ((nt_pos - 1) % 3)

        if codon_start_bio not in self.nt_map:
            return None

        col_start = self.nt_map[codon_start_bio]

        # Extract the wildtype codon from the alignment
        wt_codon = self.nt_alignment[self.ref_header][col_start : col_start + 3]

        # Inject the mutation
        pos_in_codon = (nt_pos - 1) % 3
        mut_codon = wt_codon[:pos_in_codon] + alt_allele + wt_codon[pos_in_codon + 1 :]

        return mut_codon

    def check_compensation(
        self, aa_pos: int, mut_aa: str, nt_pos: int, mut_codon: str
    ) -> dict:
        """Scans the mammalian tree for AA and NT level compensation."""
        result = {
            "aa_cdar": False,
            "nt_cdar": False,
            "aa_species": [],
            "nt_species": [],
        }

        if aa_pos not in self.aa_map or nt_pos not in self.nt_map:
            return result

        aa_col = self.aa_map[aa_pos]
        nt_col_start = self.nt_map[nt_pos - ((nt_pos - 1) % 3)]  # Start of the codon

        for species, aa_seq in self.aa_alignment.items():
            if species == self.ref_header:
                continue

            species_aa = aa_seq[aa_col]
            species_codon = self.nt_alignment[species][nt_col_start : nt_col_start + 3]

            # Skip missing data masks (CRITICAL: Includes MACSE '!' frameshift masks)
            if (
                species_aa in ["X", "-", "?", "!"]
                or "N" in species_codon
                or "!" in species_codon
            ):
                continue

            # Check AA level
            if species_aa == mut_aa:
                result["aa_cdar"] = True
                result["aa_species"].append(species)

            # Check NT level
            if species_codon == mut_codon:
                result["nt_cdar"] = True
                result["nt_species"].append(species)

        return result
