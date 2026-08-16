"""Optional Connes-CvS XC acceptance tests (skipped without the oracle)."""
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "external"))
sys.path.insert(0, str(ROOT / "src"))

HAS_CVS = importlib.util.find_spec("connes_cvs") is not None
HAS_MPMATH = importlib.util.find_spec("mpmath") is not None


@unittest.skipUnless(HAS_CVS and HAS_MPMATH, "optional connes-cvs/mpmath not installed")
class TestConnesCVSAcceptance(unittest.TestCase):
    """Acceptance gate: XC-01..03 at two precision levels."""

    def test_xc01_xc02_xc03_two_precisions(self):
        import crosschecks as xc

        for dps in xc.PRECISION_LEVELS:
            r1 = xc.xc01_archimedean_multiplier(dps)
            r2 = xc.xc02_scalar_transform(dps)
            r3 = xc.xc03_prime_power_ledger(dps)
            self.assertEqual(r1.status, "pass", r1.detail)
            self.assertEqual(r2.status, "pass", r2.detail)
            self.assertEqual(r3.status, "pass", r3.detail)


@unittest.skipUnless(HAS_MPMATH, "mpmath not installed")
class TestAtlasMpmathCoreSmoke(unittest.TestCase):
    def test_h0_at_zero(self):
        import mpmath as mp
        import mpmath_core as atlas

        L = mp.log(3)
        h0 = atlas.atlas_h0(0, L)
        self.assertEqual(h0.real, L)
        self.assertEqual(h0.imag, 0)


if __name__ == "__main__":
    unittest.main()
