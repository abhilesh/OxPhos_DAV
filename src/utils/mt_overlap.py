from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "*", "AGG": "*",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass
class MtGeneSequence:
    gene: str
    start: int
    end: int
    strand: str
    sequence: str

    def cds_index_for_genomic_pos(self, genomic_pos: int) -> int | None:
        if not (self.start <= genomic_pos <= self.end):
            return None
        if self.strand == "+":
            return genomic_pos - self.start
        return self.end - genomic_pos

    def codon_index_for_genomic_pos(self, genomic_pos: int) -> int | None:
        cds_index = self.cds_index_for_genomic_pos(genomic_pos)
        if cds_index is None:
            return None
        return cds_index // 3 + 1

    def consequence_for_variant(self, genomic_pos: int, alt_nt: str) -> dict:
        cds_index = self.cds_index_for_genomic_pos(genomic_pos)
        if cds_index is None:
            return {
                "coordinate_status": "position_outside_gene",
                "codon_index": None,
                "ref_aa": "",
                "alt_aa": "",
                "hgvs_p": "",
                "aa_change": "",
                "is_synonymous": None,
                "is_missense": None,
            }

        seq = self.sequence
        if cds_index >= len(seq):
            return {
                "coordinate_status": "position_outside_sequence",
                "codon_index": None,
                "ref_aa": "",
                "alt_aa": "",
                "hgvs_p": "",
                "aa_change": "",
                "is_synonymous": None,
                "is_missense": None,
            }

        codon_start = (cds_index // 3) * 3
        codon_end = codon_start + 3
        if codon_end > len(seq):
            return {
                "coordinate_status": "partial_terminal_codon",
                "codon_index": cds_index // 3 + 1,
                "ref_aa": "",
                "alt_aa": "",
                "hgvs_p": "",
                "aa_change": "",
                "is_synonymous": None,
                "is_missense": None,
            }

        ref_codon = seq[codon_start:codon_end]
        offset = cds_index - codon_start
        mutated_codon = list(ref_codon)
        mutated_codon[offset] = alt_nt.upper()
        mutated_codon = "".join(mutated_codon)

        ref_aa = CODON_TABLE.get(ref_codon, "X")
        alt_aa = CODON_TABLE.get(mutated_codon, "X")
        codon_index = cds_index // 3 + 1
        aa_change = f"{ref_aa}{codon_index}{alt_aa}" if ref_aa and alt_aa else ""

        if ref_aa == "X" or alt_aa == "X":
            status = "translation_failed"
        else:
            status = "resolved"

        return {
            "coordinate_status": status,
            "codon_index": codon_index,
            "ref_aa": ref_aa,
            "alt_aa": alt_aa,
            "hgvs_p": f"p.{aa_change}" if aa_change else "",
            "aa_change": aa_change,
            "is_synonymous": ref_aa == alt_aa if aa_change else None,
            "is_missense": (ref_aa != alt_aa and alt_aa != "*") if aa_change else None,
        }


def load_mtdna_gene_coords(path: Path) -> dict[str, tuple[int, int, str]]:
    coords: dict[str, tuple[int, int, str]] = {}
    with open(path, encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            gene, start, end, strand = line.rstrip("\n").split("\t")
            coords[gene] = (int(start), int(end), strand)
    return coords


def load_human_mt_sequences(alignment_dir: Path, coords: dict[str, tuple[int, int, str]]) -> dict[str, MtGeneSequence]:
    sequences: dict[str, MtGeneSequence] = {}
    for fasta in sorted(alignment_dir.glob("MT-*_codon_alignment.fasta")):
        gene = fasta.name.replace("_codon_alignment.fasta", "")
        if gene not in coords:
            continue
        human_seq = ""
        with open(fasta, encoding="utf-8") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline().strip()
                if header.startswith(">Homo_sapiens|9606|"):
                    human_seq = seq.replace("-", "").replace("*", "")
                    break
        if not human_seq:
            continue
        start, end, strand = coords[gene]
        sequences[gene] = MtGeneSequence(gene=gene, start=start, end=end, strand=strand, sequence=human_seq)
    return sequences


def parse_simple_aa_change(aa_change: str) -> tuple[str, int, str] | None:
    match = re.match(r"^([A-Za-z*]+)(\d+)([A-Za-z*]+)$", aa_change.strip())
    if not match:
        return None
    return match.group(1).upper(), int(match.group(2)), match.group(3).upper()
