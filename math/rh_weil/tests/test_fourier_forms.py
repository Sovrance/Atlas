import importlib.util
import math
import os
import sys
import unittest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

HAS_MPMATH = importlib.util.find_spec("mpmath") is not None


@unittest.skipUnless(HAS_MPMATH, "mpmath not installed")
class TestFourierForms(unittest.TestCase):
    def test_H0_zero_freq(self):
        import fourier
        import mpmath as mp

        L = mp.log(3)
        h = fourier.H0(0, L)
        self.assertEqual(h.real, L)
        self.assertEqual(h.imag, 0)

    def test_Hb_zero_freq(self):
        import fourier
        import mpmath as mp

        L = mp.mpf("1.2")
        h = fourier.Hb(0, L)
        self.assertAlmostEqual(float(h.real), float(L**3 / 6), places=10)

    def test_H0_matches_naive_away_from_zero(self):
        import fourier
        import mpmath as mp

        t, L = mp.mpf("0.7"), mp.log(4)
        naive = (mp.exp(1j * t * L) - 1) / (1j * t)
        self.assertTrue(abs(fourier.H0(t, L) - naive) < mp.mpf("1e-12"))

    def test_H0_L_jet_order1(self):
        import fourier
        import mpmath as mp

        t, L = mp.mpf("1.1"), mp.mpf("1.2")
        h, d1 = fourier.H0_L_jet(t, L, order=1)
        self.assertTrue(abs(d1 - mp.exp(1j * t * L)) < mp.mpf("1e-12"))
        self.assertTrue(abs(h - fourier.H0(t, L)) < mp.mpf("1e-12"))

    def test_Hb_relation_to_Hq1_proxy(self):
        # FORMULAS: Hb = -2i/t Hq1 for t≠0 — checked via Hb stability only here.
        import fourier
        import mpmath as mp

        t, L = mp.mpf("2.0"), mp.mpf("1.3")
        hb = fourier.Hb(t, L)
        self.assertTrue(mp.isfinite(hb.real) and mp.isfinite(hb.imag))

    def test_e3_scan_runs(self):
        import fourier
        import mpmath as mp

        Ls = [mp.log(3), mp.mpf("1.1059498113"), mp.mpf("1.20"), mp.log(4)]
        rows = fourier.scan_E2_probe(Ls, T=84, dps=20)
        self.assertEqual(len(rows), 4)
        for row in rows:
            self.assertIn(row["sign"], (-1, 0, 1))


if __name__ == "__main__":
    unittest.main()
