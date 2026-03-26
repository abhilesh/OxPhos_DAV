import sys
from pathlib import Path
import unittest
from unittest.mock import patch
import pandas as pd

# Add src to the Python path so it can find your utils module
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from utils.variant_annotators import MitoVariantAnnotator, NucVariantAnnotator


class TestVariantAnnotators(unittest.TestCase):

    @patch("utils.variant_annotators.MitoVariantAnnotator._load_mitimpact")
    @patch("utils.variant_annotators.MitoVariantAnnotator._load_phylotree")
    def test_mtdna_curation(self, mock_phylo, mock_mit):
        # Mock class initialization
        mock_mit.return_value = {}
        mock_phylo.return_value = set()
        annotator = MitoVariantAnnotator("dummy.zip", "dummy.zip")

        # Inject custom lookup data
        annotator.mitimpact_lookup = {
            (3308, "T", "C"): (0.80, "Pathogenic"),  # High Score
            (3309, "G", "A"): (0.10, "Benign"),  # Low Score
        }
        annotator.phylotree_markers = {(3310, "A")}  # Pop marker

        raw_data = pd.DataFrame(
            [
                {
                    "pos": 100,
                    "ref": "A",
                    "alt": "T",
                    "aachange": "M1T",
                    "status": "",
                },  # 1. Non-OXPHOS range
                {
                    "pos": 3308,
                    "ref": "AA",
                    "alt": "T",
                    "aachange": "M1T",
                    "status": "",
                },  # 2. Indel (len > 1)
                {
                    "pos": 3308,
                    "ref": "A",
                    "alt": "T",
                    "aachange": "frameshift",
                    "status": "",
                },  # 3. Frameshift
                {
                    "pos": 3308,
                    "ref": "A",
                    "alt": "T",
                    "aachange": "M1Ter",
                    "status": "",
                },  # 4. Nonsense 'Ter'
                {
                    "pos": 3308,
                    "ref": "A",
                    "alt": "T",
                    "aachange": "M1*",
                    "status": "",
                },  # 5. Nonsense '*'
                {
                    "pos": 3308,
                    "ref": "A",
                    "alt": "T",
                    "aachange": "M1M",
                    "status": "",
                },  # 6. Synonymous
                {
                    "pos": 3308,
                    "ref": "T",
                    "alt": "C",
                    "aachange": "M1T",
                    "status": "Cfrm [P]",
                },  # 7. Tier 1
                {
                    "pos": 3308,
                    "ref": "T",
                    "alt": "C",
                    "aachange": "M1T",
                    "status": "Reported",
                },  # 8. Tier 2 (Apogee >0.716)
                {
                    "pos": 3309,
                    "ref": "G",
                    "alt": "A",
                    "aachange": "M1T",
                    "status": "Reported",
                },  # 9. Tier 3 (Apogee <0.716)
                {
                    "pos": 3310,
                    "ref": "C",
                    "alt": "A",
                    "aachange": "M1T",
                    "status": "Reported",
                },  # 10. Discarded (PhyloTree)
                {
                    "pos": 3311,
                    "ref": "C",
                    "alt": "T",
                    "aachange": "M1T",
                    "status": "Benign",
                },  # 11. Discarded (Status)
            ]
        )

        results = annotator.curate(raw_data)

        # Exactly 5 variants should survive the strict missense filters
        self.assertEqual(len(results), 5, "Strict missense filters failed.")

        tiers = [v.tier for v in results]
        self.assertIn("Tier 1", tiers)
        self.assertIn("Tier 2", tiers)
        self.assertIn("Tier 3", tiers)
        self.assertEqual(tiers.count("Discarded"), 2)

    @patch("utils.variant_annotators.NucVariantAnnotator._load_myvariant")
    def test_nucdna_curation(self, mock_myvar):
        mock_myvar.return_value = {}
        annotator = NucVariantAnnotator("dummy.json", min_clinvar_stars=0)

        # Inject custom lookup data
        annotator.annotation_lookup = {
            "POP_ID": {"revel": 0.5, "af": 0.05}  # High AF (5%)
        }

        raw_data = pd.DataFrame(
            [
                {
                    "ReviewStatus": "practice guideline",
                    "Name": "c.1A>T (p.M1Ter)",
                    "AlleleID": "1",
                },  # 1. Stop-gain
                {
                    "ReviewStatus": "practice guideline",
                    "Name": "c.1A>T (p.Ter1M)",
                    "AlleleID": "2",
                },  # 2. Stop-loss
                {
                    "ReviewStatus": "practice guideline",
                    "Name": "c.1A>T (p.M1M)",
                    "AlleleID": "3",
                },  # 3. Synonymous
                {
                    "ReviewStatus": "no assertion provided",
                    "Name": "c.1A>T (p.M1T)",
                    "AlleleID": "4",
                },  # 4. Tier 3 (0 star)
                {
                    "ReviewStatus": "criteria provided, single submitter",
                    "Name": "c.1A>T (p.M1T)",
                    "AlleleID": "5",
                },  # 5. Tier 2 (1 star)
                {
                    "ReviewStatus": "practice guideline",
                    "Name": "c.1A>T (p.M1T)",
                    "AlleleID": "6",
                },  # 6. Tier 1 (4 star)
                {
                    "ReviewStatus": "practice guideline",
                    "Name": "c.1A>T (p.M1T)",
                    "AlleleID": "POP_ID",
                },  # 7. Discarded (AF > 0.01)
            ]
        )

        results = annotator.curate(raw_data)

        # Exactly 4 variants should survive the strict missense filters
        self.assertEqual(len(results), 4, "Nuc missense filters failed.")

        tiers = {v.ann_id: v.tier for v in results}
        self.assertEqual(tiers["4"], "Tier 3")
        self.assertEqual(tiers["5"], "Tier 2")
        self.assertEqual(tiers["6"], "Tier 1")
        self.assertEqual(tiers["POP_ID"], "Discarded")


if __name__ == "__main__":
    unittest.main(verbosity=2)
