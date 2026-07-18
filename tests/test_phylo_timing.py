import importlib.util
from pathlib import Path


def load_timing_module():
    path = Path(__file__).resolve().parents[1] / "src" / "phylo" / "02_phylogenetic_timing.py"
    spec = importlib.util.spec_from_file_location("phylogenetic_timing", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def simple_subtree_map():
    return {
        "Root": frozenset({"A", "B", "C", "D"}),
        "CladeAB": frozenset({"A", "B"}),
        "CladeCD": frozenset({"C", "D"}),
        "A": frozenset({"A"}),
        "B": frozenset({"B"}),
        "C": frozenset({"C"}),
        "D": frozenset({"D"}),
    }


def test_timing_classifies_preexisting_contact_state_as_contact_first():
    timing = load_timing_module()
    subtree = simple_subtree_map()
    contact_map = {
        "root_node": "Root",
        "node_to_children": {"Root": ["CladeAB", "CladeCD"], "CladeAB": ["A", "B"], "CladeCD": ["C", "D"]},
        "leaf_nodes": ["A", "B", "C", "D"],
        "root_states": {"10": "K"},
        "branches": {},
    }

    assert timing.timing_for_origin(
        "CladeAB",
        subtree,
        contact_map,
        subtree,
        10,
        "K",
    ) == "contact_first"


def test_timing_distinguishes_cooccurring_from_contact_after():
    timing = load_timing_module()
    subtree = simple_subtree_map()

    cooccurring_map = {
        "root_node": "Root",
        "node_to_children": {"Root": ["CladeAB", "CladeCD"], "CladeAB": ["A", "B"], "CladeCD": ["C", "D"]},
        "leaf_nodes": ["A", "B", "C", "D"],
        "root_states": {"10": "R"},
        "branches": {"Root|CladeAB": {"10": ["R", "K"]}},
    }
    assert timing.timing_for_origin(
        "CladeAB", subtree, cooccurring_map, subtree, 10, "K"
    ) == "co_occurring"

    contact_after_map = {
        "root_node": "Root",
        "node_to_children": {"Root": ["CladeAB", "CladeCD"], "CladeAB": ["A", "B"], "CladeCD": ["C", "D"]},
        "leaf_nodes": ["A", "B", "C", "D"],
        "root_states": {"10": "R"},
        "branches": {"CladeAB|A": {"10": ["R", "K"]}},
    }
    assert timing.timing_for_origin(
        "CladeAB", subtree, contact_after_map, subtree, 10, "K"
    ) == "contact_after"


def test_directional_class_labels_mechanistic_interpretation():
    timing = load_timing_module()
    counter = timing.Counter

    assert timing.directional_class(
        False, counter({"contact_first": 2})
    ) == "permissive_background"
    assert timing.directional_class(
        False, counter({"contact_after": 1})
    ) == "responding_secondary"
    assert timing.directional_class(
        False, counter({"co_occurring": 1})
    ) == "co_occurring_unresolved"
    assert timing.directional_class(
        False, counter({"contact_first": 1, "contact_after": 1})
    ) == "mixed_timing"
    assert timing.directional_class(
        True, counter({"contact_first": 3})
    ) == "ancestral_cdav_not_directional"
