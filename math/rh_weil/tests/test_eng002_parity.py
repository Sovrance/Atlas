import os
import sys
import unittest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

from interval_backend import flint_available


@unittest.skipUnless(flint_available(), "python-flint not installed")
class TestIntervalBackend(unittest.TestCase):
    def test_require_flint_and_h_plus(self):
        from archimedean import h_plus
        from interval_backend import require_flint, set_precision_bits

        require_flint()
        set_precision_bits(128)
        hp = h_plus(1)
        self.assertTrue(hp.is_finite())


@unittest.skipUnless(flint_available(), "python-flint not installed")
class TestFourierJets(unittest.TestCase):
    def test_H0_jet_matches_analytic(self):
        import mpmath as mp
        import fourier_jets as fj

        t, L = mp.mpf("0.8"), mp.mpf("1.25")
        jets = fj.H0_L_jets(t, L, order=3)
        self.assertTrue(abs(jets[1] - mp.exp(1j * t * L)) < 1e-12)
        self.assertTrue(abs(jets[2] - (1j * t) * mp.exp(1j * t * L)) < 1e-12)

    def test_Hb_second_jet(self):
        import mpmath as mp
        import fourier_jets as fj

        t, L = mp.mpf("1.3"), mp.mpf("1.1")
        jets = fj.Hb_L_jets(t, L, order=2)
        expect = L * mp.exp(1j * t * L)
        self.assertTrue(abs(jets[2] - expect) < 1e-10)


@unittest.skipUnless(flint_available(), "python-flint not installed")
class TestFiniteWeilAPI(unittest.TestCase):
    def test_even_block_keys(self):
        import math
        from finite_weil import finite_weil_even_block

        blk = finite_weil_even_block(math.log(3), T=84, precision_bits=128, n_quad=64)
        for k in ("G00", "G0b", "Gbb", "E2", "normalization", "cutoff_T"):
            self.assertIn(k, blk)
        self.assertFalse(blk["rh_proof_claim"])
        self.assertEqual(blk["cutoff_T"], 84)


class TestE3ProbeQuarantine(unittest.TestCase):
    def test_probe_flagged(self):
        import fourier

        ent = fourier.fourier_energy_probe(1.2, T=10, dps=15)
        self.assertEqual(ent["evidence_class"], "E3")
        self.assertTrue(ent.get("quarantined"))


if __name__ == "__main__":
    unittest.main()
