import importlib.util
from pathlib import Path


def load_asr_parser():
    path = Path(__file__).resolve().parents[1] / "src" / "phylo" / "01_parse_ancestral_states.py"
    spec = importlib.util.spec_from_file_location("parse_ancestral_states", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_fasta(path: Path):
    path.write_text(
        ">Homo_sapiens\n"
        "MAALLLR-----------HISY\n"
        ">Mus_musculus\n"
        "MPALLLR-----------KIS-\n"
        ">Canis_lupus\n"
        "MAALLLR-----------H-SY\n"
    )


def test_human_coordinate_maps_skip_human_gap_columns(tmp_path):
    parser = load_asr_parser()
    fasta = tmp_path / "TEST.fasta"
    write_fasta(fasta)

    maps = parser.build_human_coordinate_maps(fasta)

    assert maps["coordinate_system"] == "human_protein_position"
    assert maps["protein_pos_to_alignment_site"]["8"] == "19"
    assert maps["alignment_site_to_protein_pos"]["19"] == "8"
    assert "8" not in maps["alignment_site_to_protein_pos"]
    assert maps["n_alignment_sites"] == 22
    assert maps["n_human_protein_positions"] == 11
    assert maps["n_human_gap_alignment_sites_dropped"] == 11


def test_leaf_states_use_human_alignment_columns_not_species_ungapped_positions(tmp_path):
    parser = load_asr_parser()
    fasta = tmp_path / "TEST.fasta"
    write_fasta(fasta)
    maps = parser.build_human_coordinate_maps(fasta)

    leaf_states = parser.read_leaf_states(fasta, maps["alignment_site_to_protein_pos"])

    assert leaf_states["Homo_sapiens"]["8"] == "H"
    assert leaf_states["Mus_musculus"]["8"] == "K"
    assert "9" not in leaf_states["Canis_lupus"]


def test_node_states_and_branch_changes_are_keyed_by_human_protein_position(tmp_path):
    parser = load_asr_parser()
    fasta = tmp_path / "TEST.fasta"
    write_fasta(fasta)
    maps = parser.build_human_coordinate_maps(fasta)

    node_states = {
        "Root": {19: "H", 20: "I", 8: "Q"},
        "Clade": {19: "K", 20: "I", 8: "R"},
    }
    remapped = parser.remap_node_states_to_human_positions(
        node_states,
        maps["alignment_site_to_protein_pos"],
    )

    assert remapped == {
        "Root": {8: "H", 9: "I"},
        "Clade": {8: "K", 9: "I"},
    }
    branches = parser.build_branch_changes(
        remapped,
        {"Clade": "Root"},
    )
    assert branches == {"Root|Clade": {"8": ["H", "K"]}}
