"""WO-RH-17 / WO-RH-18 regression tests.

Locks the adjudicated normalization in place:

* the adopted pole is reproduced by an **independent** real-space route;
* the parity lemma and the parity block-diagonality hold;
* the rejected ``sqrt(3)/2`` block is exactly ``(sqrt(3)/2)cosh(L/2)`` times the
  adopted one, hence agrees at ``L = log 3`` and nowhere else;
* the repository's *odd* pivot already equalled the adopted pole;
* every certificate that depended on the disputed block is quarantined and the
  PIR guard refuses to promote it;
* no adopted-path quantity depends on the fitted constant.
"""
import json
import math
import os
import sys
import unittest

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import normalization as N  # noqa: E402

CERT_DIR = os.path.join(ROOT, "certificates")
L_POINTS = [math.log(3.0), 1.1059498113, 1.20, math.log(4.0)]
LOG3 = math.log(3.0)


class PoleDerivationTests(unittest.TestCase):
    def test_pole_matches_independent_realspace_route(self):
        """G0_ij = int_0^L K_ij(a) 2cosh(a/2) da — the same K_ij the prime block uses."""
        try:
            import mpmath as mp
        except ImportError:  # pragma: no cover
            self.skipTest("mpmath unavailable")
        mp.mp.dps = 40
        for L in L_POINTS:
            for i in N.BASIS_NAMES:
                for j in N.BASIS_NAMES:
                    closed = N.pole_entry(i, j, L)
                    quad = float(
                        mp.quad(lambda a: N.kernel_K(i, j, float(a), L) * 2 * mp.cosh(a / 2), [0, L])
                    )
                    self.assertLess(
                        abs(closed - quad) / max(1.0, abs(closed)), 1e-12, (i, j, L)
                    )

    def test_parity_lemma(self):
        """h(L-x) = ±h(x)  =>  E^- = ± e^{-L/2} E^+."""
        for L in L_POINTS:
            for name in N.BASIS_NAMES:
                sign = 1 if N.basis_parity(name) == "even" else -1
                lhs = N.E_pm(name, L, -1)
                rhs = sign * math.exp(-L / 2) * N.E_pm(name, L, 1)
                self.assertLess(abs(lhs - rhs), 1e-12 * max(1.0, abs(rhs)), (name, L))

    def test_pole_is_parity_block_diagonal(self):
        for L in L_POINTS:
            for even in ("one", "b"):
                self.assertLess(abs(N.pole_entry(even, "q1", L)), 1e-12, (even, L))

    def test_even_block_is_rank_one(self):
        """Delta = E0+Eb- - E0-Eb+ vanishes, so both candidates are rank 1.

        This is why a determinant/rank regression could not discriminate them.
        """
        for L in L_POINTS:
            delta = N.E_pm("one", L, 1) * N.E_pm("b", L, -1) - N.E_pm("one", L, -1) * N.E_pm("b", L, 1)
            self.assertLess(abs(delta), 1e-12, L)

    def test_odd_pivot_already_matched_the_repository(self):
        """E_q1^+ = 2 e^{L/4} A(L)  =>  G0[q1,q1] = -8 A^2 (the value already shipped)."""
        for L in L_POINTS:
            A = L * math.cosh(L / 4) - 4 * math.sinh(L / 4)
            self.assertLess(abs(N.E_pm("q1", L, 1) - 2 * math.exp(L / 4) * A), 1e-12 * max(1.0, abs(A)))
            self.assertLess(abs(N.pole_entry("q1", "q1", L) - (-8 * A * A)), 1e-12 * max(1.0, A * A))


class LegacyRejectionTests(unittest.TestCase):
    def test_quotient_is_sqrt3_over_2_cosh(self):
        for L in L_POINTS:
            ratio = N.legacy_pole_entry("one", "one", L) / N.pole_entry("one", "one", L)
            self.assertLess(abs(ratio - N.legacy_over_adopted_ratio(L)), 1e-12, L)

    def test_agreement_only_at_log3(self):
        self.assertLess(abs(N.legacy_over_adopted_ratio(LOG3) - 1.0), 1e-14)
        for L in (1.1059498113, 1.20, math.log(4.0)):
            self.assertGreater(abs(N.legacy_over_adopted_ratio(L) - 1.0), 1e-4, L)

    def test_discrepancy_is_L_dependent_so_not_a_basis_change(self):
        """A change of basis is a constant congruence; it cannot produce a factor
        that varies with L. The quotient does vary, so B is not a renormalised A."""
        r1 = N.legacy_over_adopted_ratio(1.20)
        r2 = N.legacy_over_adopted_ratio(math.log(4.0))
        self.assertGreater(abs(r1 - r2), 1e-3)

    def test_adopted_path_carries_no_fitted_constant(self):
        """The adopted pole is a pure E^± product: rescaling h rescales it
        bilinearly, with no residual constant left over."""
        L = 1.20
        direct = N.E_pm("one", L, 1) * N.E_pm("b", L, -1) + N.E_pm("one", L, -1) * N.E_pm("b", L, 1)
        self.assertLess(abs(N.pole_entry("one", "b", L) - direct), 1e-15 * max(1.0, abs(direct)))

    def test_legacy_model_reproduces_the_shipped_block(self):
        """The audit model must be faithful to the code it rejects."""
        try:
            from interval_backend import require_flint

            _flint, arb, _acb, _ctx = require_flint()
            import finite_weil as FW
        except Exception:  # pragma: no cover
            self.skipTest("python-flint unavailable")
        for L in (LOG3, 1.20, math.log(4.0)):
            g00, g0b, gbb = FW.g0_even_block(arb(L), arb)
            for shipped, mine in (
                (g00, N.legacy_pole_entry("one", "one", L)),
                (g0b, N.legacy_pole_entry("one", "b", L)),
                (gbb, N.legacy_pole_entry("b", "b", L)),
            ):
                self.assertLess(abs(float(shipped.mid()) - mine), 1e-9 * max(1.0, abs(mine)), L)


class QuarantineTests(unittest.TestCase):
    AFFECTED = [
        "e1_scalar_log3_log4.json",
        "e1_degree1_log3_log4.json",
        "e1_degree2_compact_log3_log4.json",
        "e1_fourier_T84_points.json",
        "e1_fourier_T84_uniform_degree2.json",
    ]

    def _load(self, name):
        with open(os.path.join(CERT_DIR, name), encoding="utf-8") as fh:
            return json.load(fh)

    def test_affected_certificates_are_quarantined(self):
        for name in self.AFFECTED:
            cert = self._load(name)
            self.assertEqual(cert.get("promotion_state"), "QUARANTINED_NORMALIZATION_ADJUDICATION", name)
            self.assertFalse(cert.get("hard_constraints_certified"), name)
            self.assertFalse(cert.get("rh_proof_claim", False), name)

    def test_quarantine_preserves_history(self):
        """Nothing is deleted and the claimed evidence class is not rewritten."""
        for name in self.AFFECTED:
            cert = self._load(name)
            q = cert.get("quarantine", {})
            self.assertIn("prior_state", q, name)
            self.assertIn("reason", q, name)
            self.assertIsNotNone(cert.get("evidence_class"), name)

    def test_pir_guard_refuses_quarantined_certificates(self):
        import pir_bridge

        refused = {r["certificate_file"] for r in pir_bridge.refused_promotions()}
        for name in self.AFFECTED:
            self.assertIn(name, refused, name)

    def test_work_order_status_flags_quarantine(self):
        st = self._load("work_order_status.json")
        for wo in ["WO-RH-05"] + [f"WO-RH-{n:02d}" for n in range(9, 16)]:
            if wo in st["orders"]:
                self.assertEqual(st["orders"][wo], "quarantined_pending_WO-RH-17", wo)


class AdjudicationCertificateTests(unittest.TestCase):
    def _cert(self, name):
        path = os.path.join(CERT_DIR, name)
        if not os.path.exists(path):
            self.skipTest(f"{name} not generated yet")
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)

    def test_adjudication_certificate_shape(self):
        c = self._cert("normalization_adjudication.json")
        for field in (
            "certificate_version", "program", "status", "active_normalization_id",
            "fourier_convention", "tilde_convention", "pole_formula", "prime_formula",
            "archimedean_formula", "crosschecks", "rh_proof_claim",
        ):
            self.assertIn(field, c, field)
        self.assertEqual(c["program"], "RH/Weil normalization adjudication")
        self.assertEqual(c["status"], "ADJUDICATED")
        self.assertIs(c["rh_proof_claim"], False)

    def test_dispositions_are_recorded(self):
        c = self._cert("normalization_adjudication.json")
        d = c["candidate_dispositions"]
        self.assertEqual(d["candidate_A_explicit_formula"]["disposition"], "ADOPTED")
        self.assertEqual(d["candidate_B_repo_sqrt3_over_2"]["disposition"], "REJECTED")

    def test_normalization_id_is_stable_and_active(self):
        c = self._cert("normalization_adjudication.json")
        self.assertEqual(c["active_normalization_id"], N.normalization_id())
        self.assertEqual(N.normalization_id(), N.normalization_id())

    def test_crosscheck_certificate_has_no_disagreements(self):
        c = self._cert("normalization_crosscheck.json")
        self.assertEqual(c["disagreements"], [], c["disagreements"][:2])
        self.assertEqual(c["status"], "AGREE")
        self.assertGreater(c["coverage"]["cells_compared"], 0)


if __name__ == "__main__":
    unittest.main()
