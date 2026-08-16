import importlib.util
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))

from rh_weil.external.connes_cvs_adapter import dependency_info, TESTED_VERSION


class TestConnesCVSAdapterContract(unittest.TestCase):
    def test_dependency_metadata_is_noncanonical(self):
        info = dependency_info()
        self.assertEqual(info.tested_version, TESTED_VERSION)
        self.assertFalse(info.canonical_proof_engine)
        self.assertEqual(info.role, "external_cross_validation_oracle")

    def test_availability_matches_importability(self):
        info = dependency_info()
        available = importlib.util.find_spec("connes_cvs") is not None
        self.assertEqual(info.available, available)


if __name__ == "__main__":
    unittest.main()
