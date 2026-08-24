"""ENG-005 acceptance tests — recovery of the Candidate-A E1 chain.

Covers §1 (reproducible Ginf'' derivation, without overclaiming), §2 (canonical
deterministic panel integration), §3 (the invalid tail assumption stays
rejected), §4/§5 (degree-1 and compact degree-2), §6/§9 (fresh T=84 topology and
exact jets), §7/§10 (T=84 point and uniform E1), and the standing rule that
Candidate-B code is unreachable from any derivative path.

The heavy T=84 integrals are exercised once, not per-assertion: several tests read
the committed certificates rather than recomputing them, which is also what makes
them meaningful as *regression* tests of the artifacts the program ships.
"""
from __future__ import annotations

import json
import math
import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import promotion  # noqa: E402
import rigorous_integration as RI  # noqa: E402

CERT_DIR = ROOT / "certificates"


def _load(name: str):
    p = CERT_DIR / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def _flint():
    try:
        from interval_backend import require_flint, set_precision_bits

        _, arb, acb, _ = require_flint()
        set_precision_bits(200)
        return arb, acb
    except Exception:
        raise unittest.SkipTest("python-flint unavailable")


# --------------------------------------------------------------------------- #
# §1 — the curvature derivation                                                #
# --------------------------------------------------------------------------- #
class CurvatureDerivation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import curvature_derivation as CD
        except Exception:
            raise unittest.SkipTest("curvature_derivation unavailable")
        cls.CD = CD
        cls.report = CD.derivation_report()

    def test_every_symbolic_step_reproduces(self):
        for step in self.report["steps"]:
            self.assertTrue(step["verified"], step["step"])
        self.assertEqual(self.report["status"], "VERIFIED")

    def test_independent_series_regression_agrees(self):
        self.assertTrue(self.report["all_regression_rows_agree"])
        for row in self.report["independent_series_regression"]:
            self.assertLess(float(row["abs_diff"]), 1e-12)

    def test_does_not_overclaim(self):
        """§1: the interchange is argued, not machine-checked, and must say so."""
        self.assertEqual(self.report["evidence_class"], "E0")
        self.assertTrue(self.report["analytic_hypotheses"])
        self.assertIn("NOT", self.CD.INTERCHANGE_HYPOTHESIS)

    def test_geometric_sum_is_finite_identity_plus_limit(self):
        """Asserting the closed form without the ratio condition would be unearned."""
        step = next(s for s in self.report["steps"] if s["step"] == "geometric sum over n")
        self.assertTrue(step["finite_sum_verified"])
        self.assertTrue(step["limit_verified"])

    def test_sympy_is_required_not_optional(self):
        self.assertTrue(issubclass(self.CD.SymPyRequired, RuntimeError))


# --------------------------------------------------------------------------- #
# §2 — canonical panel integration                                             #
# --------------------------------------------------------------------------- #
class PanelIntegration(unittest.TestCase):
    def test_t84_schedule_is_the_one_the_spec_names(self):
        self.assertEqual(RI.panel_schedule(84), list(RI.PANELS_T84))
        self.assertEqual(RI.PANELS_T84[0], (0.0, 1.0))
        self.assertEqual(RI.PANELS_T84[-1], (64.0, 84.0))

    def test_schedule_is_deterministic(self):
        for T in (84, 1000, 20000, 200000):
            self.assertEqual(RI.panel_schedule(T), RI.panel_schedule(T))

    def test_schedule_covers_exactly_zero_to_T(self):
        """The bug this replaced integrated past T for any T below its last edge."""
        for T in (84, 500, 20000, 50000, 200000):
            sched = RI.panel_schedule(T)
            self.assertEqual(sched[0][0], 0.0, T)
            self.assertAlmostEqual(sched[-1][1], float(T), places=9, msg=T)
            for (a1, b1), (a2, b2) in zip(sched, sched[1:]):
                self.assertEqual(b1, a2, "panels must abut")

    def test_a_schedule_that_misses_T_is_refused(self):
        arb, acb = _flint()
        with self.assertRaises(ValueError):
            RI.rigorous_panel_integral(lambda z, a: acb(1), 84, acb,
                                       panels=[(0.0, 10.0)])

    def test_trapezoid_path_may_not_emit_e1(self):
        with self.assertRaises(RI.QuadratureFailure):
            RI.assert_not_trapezoid_path({"method": "trapezoid_global_M2",
                                          "trapezoid_path_used": True})
        RI.assert_not_trapezoid_path({"method": "arb_acb_integral_panelled"})

    def test_non_finite_panel_raises_rather_than_returning_infinity(self):
        arb, acb = _flint()

        def bad(z, _a):
            return acb(arb("inf"))

        with self.assertRaises(RI.QuadratureFailure):
            RI.rigorous_panel_integral(bad, 84, acb)


# --------------------------------------------------------------------------- #
# §3 — the tail lemma                                                          #
# --------------------------------------------------------------------------- #
class TailLemma(unittest.TestCase):
    def test_invalid_assumption_is_rejected_near_t_equals_two(self):
        arb, _acb = _flint()
        import scalar_canary as SC

        row = SC.invalid_assumption_is_rejected(arb)
        self.assertEqual(row["verdict"], "REJECTED")
        self.assertTrue(row["exceeds_one"])
        self.assertGreater(float(row["t_h_plus_prime_enclosure"][0]), 1.0)

    def test_kappa_depends_on_the_tail_domain(self):
        import scalar_canary as SC

        self.assertGreater(SC.lemma_A_constant(2.0), 1.0)
        self.assertLess(SC.lemma_A_constant(200000.0), 1.0001)
        self.assertGreater(SC.lemma_A_constant(2.0), SC.lemma_A_constant(200.0))

    def test_kappa_rejects_a_nonpositive_domain(self):
        import scalar_canary as SC

        with self.assertRaises(ValueError):
            SC.lemma_A_constant(0.0)


# --------------------------------------------------------------------------- #
# §4/§5 — degree-1 and compact degree-2                                        #
# --------------------------------------------------------------------------- #
class CutoffFreeRecovery(unittest.TestCase):
    FILES = {"degree1": "e1_degree1_log3_log4.json",
             "degree2": "e1_degree2_compact_log3_log4.json"}

    def test_both_are_promoted_with_positive_bounds(self):
        for label, name in self.FILES.items():
            cert = _load(name)
            self.assertIsNotNone(cert, name)
            self.assertEqual(cert["status"], "PASS", label)
            self.assertEqual(cert["promotion_state"], promotion.PROMOTED_STATE, label)
            self.assertGreater(float(cert["certified_lower_bound"]), 0.0, label)
            self.assertIsNone(promotion.promotion_refusal(cert), label)

    def test_they_are_cutoff_free_not_truncated(self):
        for name in self.FILES.values():
            cert = _load(name)
            self.assertTrue(cert["cutoff_free"])
            self.assertIn("no frequency cutoff", cert["archimedean_route"]["method"])

    def test_parity_identities_hold(self):
        for name in self.FILES.values():
            cert = _load(name)
            self.assertTrue(cert["parity_identities"]["all_hold"], name)

    def test_uniform_over_the_closed_cell(self):
        for name in self.FILES.values():
            cert = _load(name)
            self.assertEqual(cert["domain"]["L_interval"], ["log(3)", "log(4)"])
            self.assertTrue(cert["domain"]["closed"])
            self.assertGreater(cert["subdivision_statistics"]["boxes_examined"], 0)

    def test_historical_values_are_not_a_warrant(self):
        for name in self.FILES.values():
            self.assertIn("regression evidence only", _load(name)["historical_values_note"])

    def test_real_space_and_frequency_routes_agree(self):
        """The two archimedean routes are independent; they must not diverge."""
        arb, acb = _flint()
        import archimedean_realspace as AR
        import weil_entries as WE

        L = arb("1.2")
        for i, j in (("one", "b"), ("b", "b")):
            rs, _ = AR.arch_entry_realspace(i, j, L, arb, acb)
            fs, _ = WE.arch_entry(i, j, L, 200000, arb, acb)
            # These entries decay fast enough that T=2e5 truncation is negligible.
            self.assertLess(abs(float(rs) - float(fs)), 1e-12, (i, j))


# --------------------------------------------------------------------------- #
# §6/§9 — fresh topology and exact jets                                        #
# --------------------------------------------------------------------------- #
class T84TopologyAndJets(unittest.TestCase):
    def test_scan_is_candidate_A(self):
        scan = _load("e3_fourier_T84_scan.json")
        self.assertIsNotNone(scan)
        self.assertEqual(scan["pole_candidate"], "A")
        self.assertFalse(scan["topology"]["candidate_b_topology_reused"])

    def test_superseded_candidate_b_scan_is_preserved(self):
        kept = CERT_DIR / "history" / "e3_fourier_T84_scan_candidateB_superseded.json"
        if not kept.exists():
            self.skipTest("no Candidate-B scan was present to supersede")
        body = json.loads(kept.read_text(encoding="utf-8"))
        self.assertEqual(body["pole_candidate"], "B_REJECTED")
        self.assertIn("provenance", body["retained_as"])

    def test_exact_jets_match_finite_differences(self):
        """§9 forbids finite differences in E1 — so they are only a *check* here."""
        arb, acb = _flint()
        import t84

        L, h = 1.25, 1e-4
        for key, (i, j) in (("Gbb", ("b", "b")), ("O1", ("q1", "q1"))):
            d1, _ = t84.entry_jet(i, j, 1, arb(repr(L)), arb, acb)
            vp, _ = t84.entry_jet(i, j, 0, arb(repr(L + h)), arb, acb)
            vm, _ = t84.entry_jet(i, j, 0, arb(repr(L - h)), arb, acb)
            fd = float((vp - vm) / (2 * h))
            self.assertLess(abs(float(d1) - fd), 1e-6 * max(1.0, abs(fd)), key)

    def test_second_kernel_derivative_coefficients_are_right(self):
        """A wrong L/2-for-L here threw d^2 O1 off by 0.7 and nothing else caught it."""
        arb, acb = _flint()
        import t84

        L = arb("1.25")
        # d^2/dL^2 of K_q1q1 = L^3/6 - (L^2/2)u + u^3/3 is [L, -1, 0, 0].
        coeffs = t84.kernel_coeffs_d2L_in_u("q1", "q1", L, acb)
        self.assertAlmostEqual(float(coeffs[0].real), 1.25, places=12)
        self.assertAlmostEqual(float(coeffs[1].real), -1.0, places=12)

    def test_pole_second_derivative_matches_a_stencil(self):
        import pole

        h = 1e-3
        for L in (math.log(3.0), 1.2, math.log(4.0)):
            for i in ("one", "q1", "b"):
                an = pole.pole_gram_entry_d2L(i, i, L)
                f = lambda x: pole.pole_gram_entry(i, i, x)  # noqa: E731
                fd = (-f(L - 2 * h) + 16 * f(L - h) - 30 * f(L)
                      + 16 * f(L + h) - f(L + 2 * h)) / (12 * h * h)
                self.assertLess(abs(an - fd), 1e-6 * max(1.0, abs(fd)), (i, L))


# --------------------------------------------------------------------------- #
# §7/§10 — T=84 point and uniform E1                                           #
# --------------------------------------------------------------------------- #
class T84Certificates(unittest.TestCase):
    def test_point_certificate_is_promoted_and_point_scoped(self):
        cert = _load("e1_fourier_T84_points.json")
        self.assertIsNotNone(cert)
        self.assertEqual(cert["status"], "PASS")
        self.assertEqual(cert["promotion_state"], promotion.PROMOTED_STATE)
        self.assertTrue(cert["point_scoped"])
        self.assertIsNone(promotion.promotion_refusal(cert))
        for row in cert["points"]:
            self.assertTrue(row["E2_definitely_positive"], row["label"])
            self.assertTrue(row["O1_definitely_positive"], row["label"])

    def test_points_include_the_ones_the_spec_names(self):
        cert = _load("e1_fourier_T84_points.json")
        labels = {r["label"] for r in cert["points"]}
        self.assertIn("log3", labels)
        self.assertIn("1.20", labels)
        self.assertIn("log4", labels)

    def test_uniform_certificate_is_promoted_with_a_positive_bound(self):
        cert = _load("e1_fourier_T84_uniform_degree2.json")
        self.assertIsNotNone(cert)
        self.assertEqual(cert["status"], "PASS")
        self.assertEqual(cert["promotion_state"], promotion.PROMOTED_STATE)
        self.assertGreater(float(cert["certified_lower_bound"]), 0.0)
        self.assertIsNone(promotion.promotion_refusal(cert))

    def test_uniform_certificate_describes_the_topology_it_proved(self):
        """§8: no precommitment; say what was actually established."""
        cert = _load("e1_fourier_T84_uniform_degree2.json")
        topo = cert["certified_topology"]
        self.assertIn("classification", topo)
        self.assertFalse(topo["candidate_b_topology_reused"])
        self.assertIn("assumes no topology", topo["warrant"])
        cov = cert["interval_coverage"]
        self.assertIn("no convexity or monotonicity assumed", cov["topology_proved"])

    def test_uniform_certificate_records_its_panel_schedule(self):
        cert = _load("e1_fourier_T84_uniform_degree2.json")
        self.assertEqual([tuple(p) for p in cert["quadrature"]["panel_schedule"]],
                         list(RI.PANELS_T84))


# --------------------------------------------------------------------------- #
# Standing rules                                                               #
# --------------------------------------------------------------------------- #
class StandingRules(unittest.TestCase):
    RECOVERED = ("e1_scalar_log3_log4.json", "e1_degree1_log3_log4.json",
                 "e1_degree2_compact_log3_log4.json", "e1_fourier_T84_points.json",
                 "e1_fourier_T84_uniform_degree2.json")

    def test_no_certificate_claims_rh(self):
        for p in sorted(CERT_DIR.rglob("*.json")):
            body = json.loads(p.read_text(encoding="utf-8"))
            if isinstance(body, dict):
                self.assertFalse(body.get("rh_proof_claim", False), p.name)

    def test_every_recovered_certificate_binds_the_active_normalization_id(self):
        active = promotion.active_normalization_id()
        for name in self.RECOVERED:
            cert = _load(name)
            self.assertIsNotNone(cert, name)
            self.assertEqual(cert[promotion.NORMALIZATION_ID_FIELD], active, name)

    def test_no_recovered_certificate_used_mpmath(self):
        for name in self.RECOVERED:
            self.assertFalse(_load(name).get("mpmath_used", False), name)

    def test_derivative_paths_do_not_reach_candidate_b(self):
        """Candidate-B derivative code must not be reachable (§9)."""
        import ast

        for mod in ("t84.py", "archimedean_realspace.py", "e1_t84.py",
                    "e1_cutoff_free.py", "interval_cover.py"):
            src = (ROOT / "src" / mod).read_text(encoding="utf-8")
            tree = ast.parse(src)
            names = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    names.update(a.name.split(".")[0] for a in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    names.add(node.module.split(".")[0])
            self.assertNotIn("rejected_pole", names, mod)

    def test_source_hashes_are_current(self):
        for name in self.RECOVERED:
            self.assertEqual(promotion.stale_dependencies(_load(name)), [], name)


if __name__ == "__main__":
    unittest.main()
