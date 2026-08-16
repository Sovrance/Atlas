import importlib.util
import os
import sys
import unittest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

from interval_backend import flint_available

HAS_MPMATH = importlib.util.find_spec("mpmath") is not None


@unittest.skipUnless(flint_available(), "python-flint not installed")
class TestIntervalBackend(unittest.TestCase):
    def test_require_flint_and_h_plus(self):
        from archimedean import h_plus
        from interval_backend import require_flint, set_precision_bits

        require_flint()
        set_precision_bits(128)
        hp = h_plus(1)  # must preserve caller context (128), not reset to 256
        self.assertTrue(hp.is_finite())
        _, _, _, ctx = require_flint()
        self.assertEqual(int(ctx.prec), 128)


@unittest.skipUnless(HAS_MPMATH, "mpmath not installed")
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

    def test_Hb_zero_freq_third_derivative(self):
        import mpmath as mp
        import fourier_jets as fj

        jets = fj.Hb_L_jets(0, mp.mpf("1.2"), order=4)
        self.assertEqual(jets[3], mp.mpc(1, 0))
        self.assertEqual(jets[4], mp.mpc(0, 0))


@unittest.skipUnless(flint_available(), "python-flint not installed")
class TestFiniteWeilAPI(unittest.TestCase):
    def test_even_block_keys(self):
        import math
        from finite_weil import finite_weil_even_block

        blk = finite_weil_even_block(
            math.log(3), T=84, precision_bits=128, n_quad=64, backend="flint"
        )
        for k in ("G00", "G0b", "Gbb", "E2", "normalization", "cutoff_T"):
            self.assertIn(k, blk)
        self.assertFalse(blk["rh_proof_claim"])
        self.assertEqual(blk["cutoff_T"], 84)
        self.assertEqual(blk["pole_scale"], "sqrt(3)/2")


@unittest.skipUnless(flint_available(), "python-flint not installed")
class TestEvenPoleOuterProduct(unittest.TestCase):
    def test_pole_det_zero_and_log3_regression(self):
        import math
        from finite_weil import g0_even_block, finite_weil_even_block
        from interval_backend import require_flint, set_precision_bits

        _, arb, _, _ = require_flint()
        set_precision_bits(192)
        L = arb(math.log(3))
        g00, g0b, gbb = g0_even_block(L, arb)
        det = g00 * gbb - g0b * g0b
        self.assertLess(abs(float(det.mid())), 1e-20)
        blk = finite_weil_even_block(
            math.log(3), T=84, n_quad=2048, precision_bits=192, rigorous=False
        )
        self.assertAlmostEqual(float(blk["G00"].mid()), 0.107356700414591762, places=6)
        self.assertAlmostEqual(float(blk["E2"].mid()), 3.4640947469748e-6, places=10)


@unittest.skipUnless(flint_available(), "python-flint not installed")
class TestArchimedeanTail(unittest.TestCase):
    def test_tail_positive(self):
        import math
        from finite_weil import archimedean_tail_even
        from interval_backend import require_flint, set_precision_bits

        _, arb, _, _ = require_flint()
        set_precision_bits(128)
        t00, t0b, tbb = archimedean_tail_even(arb(math.log(3)), 84.0, arb)
        self.assertGreater(float(t00.lower()), 0)
        self.assertGreater(float(tbb.lower()), 0)



@unittest.skipUnless(HAS_MPMATH, "mpmath not installed")
class TestE3ProbeQuarantine(unittest.TestCase):
    def test_probe_flagged(self):
        import fourier

        ent = fourier.fourier_energy_probe(1.2, T=10, dps=15)
        self.assertEqual(ent["evidence_class"], "E3")
        self.assertTrue(ent.get("quarantined"))


if __name__ == "__main__":
    unittest.main()
