from pathlib import Path
from Bio import SeqIO
from collections import defaultdict


class AlignmentParser:
    def __init__(self, aa_fasta: Path, nt_fasta: Path, genome: str,
                 tx_pos_map: dict | None = None):
        """
        aa_fasta, nt_fasta : TOGA alignment files for this gene.
        genome             : "mtDNA" | "nucDNA"
        tx_pos_map         : optional dict {nm_aa_pos (int): enst_aa_pos (int|None)}
                             pre-built by 00f_build_transcript_position_maps.py.
                             When provided, check_compensation uses it as a direct
                             NM_→ENST position lookup instead of find_sequence_anchor.
        """
        self.genome = genome
        self.aa_alignment = self._load_fasta(aa_fasta)
        self.nt_alignment = self._load_fasta(nt_fasta)

        self.ref_header = self._identify_reference()
        # Extract ENST from anywhere in the header (format varies between TOGA codon
        # and AA alignments, but ENST always appears as a distinct "|"-separated field
        # or as the first word in a field).  Fall back to "UNKNOWN" if absent.
        import re as _re
        _enst = _re.search(r"(ENST\d+)", self.ref_header)
        self.transcript_id = _enst.group(1) if _enst else "UNKNOWN"

        # Maps biological positions to alignment column indices
        self.aa_map = self._build_coordinate_map(self.aa_alignment[self.ref_header])
        self.nt_map = self._build_coordinate_map(self.nt_alignment[self.ref_header])

        # tx_pos_map: int keys (JSON stores them as str → convert on load)
        self.tx_pos_map: dict[int, int | None] | None = (
            {int(k): v for k, v in tx_pos_map.items()} if tx_pos_map else None
        )

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
        self, reported_aa_pos: int, wt_aa: str, mut_aa: str, nt_pos: int, alt_nt: str
    ) -> dict:
        """Evaluates cDAVs using sequence-anchored coordinates.

        Codon construction happens *after* anchor correction, so it uses the
        TOGA-aligned position rather than the raw ClinVar c. coordinate (which
        is in NM_ space; the alignment may use a different ENST isoform).

        wt_aa: wild-type AA from the variant record — used as the positional
               anchor. Do NOT derive from the alignment; TOGA alignments may
               start partway into the CDS.

        Returns a dict with:
            aa_cdav          bool        mutant AA found in ≥1 non-human species
            nt_cdav          bool        mutant codon found in ≥1 non-human species
            aa_species       list[str]   species names with AA cDAV
            nt_species       list[str]   species names with NT cDAV
            anchor_found     bool        False → wt_aa not in alignment within ±10 aa
            mut_codon        str|None    mutant codon at corrected position; None if
                                         corrected position is outside the alignment
            ref_base_found   str         actual ref base at corrected position
                                         ("ANCHOR_NOT_FOUND" or "POS_NOT_IN_MAP" on failure)
            corrected_nt_pos int         NT position after anchor shift
        """
        ref_seq = self.nt_alignment[self.ref_header]

        results = {
            "aa_cdav": False,
            "nt_cdav": False,
            "aa_species": [],
            "nt_species": [],
            "anchor_found": False,
            "position_not_in_enst": False,
            "mut_codon": None,
            "ref_base_found": "ANCHOR_NOT_FOUND",
            "corrected_nt_pos": nt_pos,
        }

        # ── Position resolution ───────────────────────────────────────────────
        # Strategy 1 (preferred): global NM_→ENST map pre-built by
        #   00f_build_transcript_position_maps.py.  Direct O(1) lookup; handles
        #   large offsets and fully different N-termini.
        # Strategy 2 (fallback): local ±10-residue sequence anchor search.
        #   Used when no map is available (mtDNA, genes not yet in the map file).
        if self.tx_pos_map is not None:
            enst_aa_pos = self.tx_pos_map.get(reported_aa_pos)
            if enst_aa_pos is None:
                # Position exists only in NM_ isoform (aligns to a gap in ENST)
                results["position_not_in_enst"] = True
                return results
            true_aa_pos = enst_aa_pos
        else:
            true_aa_pos = self.find_sequence_anchor(reported_aa_pos, wt_aa)
            if not true_aa_pos:
                return results

        results["anchor_found"] = True
        aa_col = self.aa_map.get(true_aa_pos)
        if aa_col is None:
            # Mapped ENST position is outside the alignment range
            return results

        # NT shift: same amount the AA position was corrected.
        aa_shift = true_aa_pos - reported_aa_pos
        corrected_nt_pos = nt_pos + aa_shift * 3
        results["corrected_nt_pos"] = corrected_nt_pos

        if corrected_nt_pos not in self.nt_map:
            results["ref_base_found"] = "POS_NOT_IN_MAP"
            return results

        # Build mutant codon at the corrected position
        frame_offset = (corrected_nt_pos - 1) % 3
        codon_start_pos = corrected_nt_pos - frame_offset
        nt_cols = [
            self.nt_map.get(codon_start_pos),
            self.nt_map.get(codon_start_pos + 1),
            self.nt_map.get(codon_start_pos + 2),
        ]
        if None in nt_cols:
            results["ref_base_found"] = "POS_NOT_IN_MAP"
            return results

        ref_codon_bases = [ref_seq[c] for c in nt_cols]
        results["ref_base_found"] = ref_codon_bases[frame_offset]  # base at exact nt_pos

        mut_codon_bases = ref_codon_bases[:]
        mut_codon_bases[frame_offset] = alt_nt.upper()
        mut_codon = "".join(mut_codon_bases)
        results["mut_codon"] = mut_codon

        # Species scan
        _MASK = {"-", "!", "*", "X"}
        for species, aa_seq in self.aa_alignment.items():
            if species == self.ref_header:
                continue
            if aa_seq[aa_col] in _MASK or any(
                self.nt_alignment[species][c] in _MASK for c in nt_cols
            ):
                continue

            if aa_seq[aa_col] == mut_aa:
                results["aa_cdav"] = True
                sp_name = species.split("|")[0]
                results["aa_species"].append(sp_name)

                species_codon = "".join(self.nt_alignment[species][c] for c in nt_cols)
                if species_codon == mut_codon:
                    results["nt_cdav"] = True
                    results["nt_species"].append(sp_name)

        return results
