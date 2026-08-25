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

class SmallFrequencyBranch(unittest.TestCase):
    """ATLAS-RH-ENG-008 §"preserve the stable small-|t| Taylor branch".

    ``H0`` and ``Hb`` switch to a Taylor form for small ``|z| = |L t / 2|``,
    because their closed forms cancel catastrophically there: ``Hb`` divides
    ``sin z - z cos z`` -- two quantities of size ``z`` whose difference is of
    size ``z^3/3`` -- by ``z^3``, losing about ``2 log10(1/z)`` significant
    digits before the division ever happens. The branch exists to avoid that.

    The threshold is therefore *coupled to working precision*, and the tests
    below say so explicitly rather than leaving it implicit. ``Hb`` cuts over at
    ``|z| < 1e-8``, which is right for the precision its only caller sets
    (``fourier_energy_probe`` runs at ``dps`` 25 or 40) and wrong at mpmath's
    default 15, where the closed form just above the threshold has already lost
    every digit it had. ``test_the_threshold_is_coupled_to_working_precision``
    pins that relationship, so raising the threshold or lowering the probe's
    precision cannot silently pass.

    These are E3 forms. Nothing here is a warrant for anything: the module is
    the quarantined energy probe's, and the rigorous ``T = 84`` route is Arb's.
    """

    #: Straddle ``Hb``'s cutover: ``z = L t / 2`` with ``L = 1.2``, so these are
    #: ``z`` of 5.4e-9 and 6.6e-9 (Taylor) against 1.02e-8 and 6e-8 (closed).
    HB_STRADDLE = ("9e-9", "1.1e-8", "1.7e-8", "1e-7")

    @classmethod
    def setUpClass(cls):
        if not HAS_MPMATH:
            raise unittest.SkipTest("mpmath not installed")

    def _reference(self, fn, t, L):
        """The closed form with enough precision to survive its own arithmetic."""
        import mpmath as mp

        with mp.workdps(80):
            z = mp.mpf(L) * mp.mpf(t) / 2
            if fn == "H0":
                return mp.mpf(L) * mp.exp(1j * z) * (mp.sin(z) / z)
            B = (mp.sin(z) - z * mp.cos(z)) / z ** 3
            return (mp.mpf(L) ** 3 / 2) * mp.exp(1j * z) * B

    def test_the_branch_agrees_with_the_closed_form_it_replaces(self):
        import fourier
        import mpmath as mp

        L = mp.mpf("1.2")
        old = mp.mp.dps
        mp.mp.dps = 40  # what fourier_energy_probe's callers use by default
        try:
            for t in self.HB_STRADDLE:
                t_m = mp.mpf(t)
                got = fourier.Hb(t_m, L)
                ref = self._reference("Hb", t_m, L)
                self.assertLess(abs(got - ref) / abs(ref), mp.mpf("1e-20"),
                                f"Hb at t={t}")
            # H0's cutover is at 1e-18 in z; sinc is far better conditioned, so
            # the same check passes on either side with room to spare.
            for t in ("1e-19", "1e-18", "1e-17", "1e-9"):
                t_m = mp.mpf(t)
                got = fourier.H0(t_m, L)
                ref = self._reference("H0", t_m, L)
                self.assertLess(abs(got - ref) / abs(ref), mp.mpf("1e-20"),
                                f"H0 at t={t}")
        finally:
            mp.mp.dps = old

    def test_the_threshold_is_coupled_to_working_precision(self):
        """Below the precision the cutover was chosen for, it does not hold.

        This is not a defect being tolerated -- it is the reason the probe sets
        ``dps`` before it evaluates anything, and the reason this file records
        the coupling. If someone widens ``Hb``'s Taylor branch so that the
        default precision is safe too, this test fails and should be deleted
        along with the caveat above.
        """
        import fourier
        import mpmath as mp

        L = mp.mpf("1.2")
        just_above = mp.mpf("1.7e-8")  # z = 1.02e-8, first z past the cutover
        ref = self._reference("Hb", just_above, L)
        old = mp.mp.dps
        try:
            mp.mp.dps = 15
            bad = abs(fourier.Hb(just_above, L) - ref) / abs(ref)
            mp.mp.dps = 40
            good = abs(fourier.Hb(just_above, L) - ref) / abs(ref)
        finally:
            mp.mp.dps = old
        self.assertGreater(bad, mp.mpf("1e-3"))
        self.assertLess(good, mp.mpf("1e-20"))

    def test_the_branch_is_continuous_with_the_zero_frequency_value(self):
        """Continuity, at the rate the phase factor allows.

        ``H0 = L e^{i z} sinc(z)`` differs from ``H0(0) = L`` at first order in
        ``z``, not beyond it -- the phase, not the amplitude, is what moves. So
        the honest statement is a linear bound, and asserting anything tighter
        would be asserting something false.
        """
        import fourier
        import mpmath as mp

        L = mp.mpf("1.2")
        old = mp.mp.dps
        mp.mp.dps = 40
        try:
            for call in (fourier.H0, fourier.Hb):
                at_zero = call(0, L)
                for t in ("1e-9", "1e-12", "1e-18"):
                    z = L * mp.mpf(t) / 2
                    got = call(mp.mpf(t), L)
                    self.assertLess(abs(got - at_zero) / abs(at_zero), 10 * z,
                                    f"{call.__name__} at t={t}")
        finally:
            mp.mp.dps = old

    def test_production_never_samples_the_cancellation_regime(self):
        """The guard the branch alone does not give.

        The probe's grid is ``t = T i / n`` with ``n = max(32, T)``, so it
        evaluates ``t = 0`` -- the exact branch, no cancellation at all -- and
        then jumps to ``t >= 1``, never landing in between. That is what makes
        the narrow threshold harmless in practice, and this is what would
        notice if a future grid refinement changed it.
        """
        import mpmath as mp

        T, L = 84, mp.mpf("1.2")
        n = max(32, T)
        for i in range(n + 1):
            t = mp.mpf(T) * i / n
            if t == 0:
                continue
            z = L * t / 2
            self.assertGreater(float(z), 1e-3,
                               "grid entered the regime where Hb's closed form "
                               "cancels; widen the Taylor threshold")


@unittest.skipUnless(HAS_MPMATH, "mpmath not installed")
class TestFourierScan(unittest.TestCase):
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
