#!/usr/bin/env python3
"""ATLAS-RH-ENG-006 §6 — spectral moments and the Atlas B1 adapter (WO-RH-31).

Two things are being defended here. First that the moments themselves are
enclosures of traces of matrix powers, with no eigenvalue solver anywhere near
an E1 path (§14.4). Second, and more easily got wrong, that the adapter reports
the *available* conclusion and not the convenient one: a truncated localizing
matrix being PSD is necessary and not sufficient for a non-negative spectrum, so
"the moments force PSD" must come back INSUFFICIENT_INFORMATION even when the
matrix is in fact positive definite.
"""
from __future__ import annotations

import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import inertia.congruence as C  # noqa: E402
from inertia.ldl import exact_inertia  # noqa: E402
from moments.adapter import analyse  # noqa: E402
from moments.feasible_spectrum import (  # noqa: E402
    CONCLUSIVE,
    INSUFFICIENT,
    b1_available,
    eigenvalue_bounds_from_two_moments,
)
from moments.spectral_moments import (  # noqa: E402
    sanity_violations,
    spectral_moments,
    trace_of_power,
)


def _flint():
    try:
        from flint import arb, ctx

        ctx.prec = 140
        return arb
    except ImportError:  # pragma: no cover
        return None


def _q(report, name):
    return next(q for q in report["b1_queries"] if q["query"] == name)


class MomentExtraction(unittest.TestCase):
    def test_moments_match_sums_of_eigenvalue_powers(self):
        for signs in ([2, 3], [2, -3], [-2, -3], [1, 1, -1], [5, -1, 2]):
            G = C.diagonal(signs)
            ms = spectral_moments(G)
            for k in range(1, 5):
                with self.subTest(signs=tuple(signs), k=k):
                    self.assertEqual(ms[f"m{k}"], sum(F(s) ** k for s in signs))

    def test_moments_are_invariant_under_orthogonal_like_congruence(self):
        """tr(G^k) is a similarity invariant; check on a nontrivially rotated block."""
        D = C.diagonal([2, -3])
        # S^-1 A S with S unimodular is a similarity, so moments must not move.
        S = C.unimodular(2, 5)
        Sinv_A_S = C.matmul(C.matmul(_inverse(S), D), S)
        self.assertEqual(spectral_moments(Sinv_A_S), spectral_moments(D))

    def test_m2_equals_the_hilbert_schmidt_norm_squared(self):
        G = C.congruence(C.diagonal([1, -2, 3]), C.unimodular(3, 4))
        hs = sum(G[i][j] * G[i][j] for i in range(3) for j in range(3))
        self.assertEqual(trace_of_power(G, 2), hs)

    def test_even_moments_are_non_negative(self):
        for seed in range(1, 6):
            G = C.congruence(C.diagonal([1, -1, 2]), C.unimodular(3, seed))
            ms = spectral_moments(G)
            with self.subTest(seed=seed):
                self.assertGreaterEqual(ms["m2"], 0)
                self.assertGreaterEqual(ms["m4"], 0)
                self.assertEqual(sanity_violations(ms), [])

    @unittest.skipIf(_flint() is None, "python-flint not installed")
    def test_interval_moments_enclose_the_exact_ones(self):
        arb = _flint()
        G = C.congruence(C.diagonal([1, -1, 2]), C.unimodular(3, 3))
        exact = spectral_moments(G)
        ball = spectral_moments([[arb(str(x)) for x in row] for row in G])
        for k in range(1, 5):
            with self.subTest(k=k):
                lo, hi = float(ball[f"m{k}"].lower()), float(ball[f"m{k}"].upper())
                self.assertLessEqual(lo, float(exact[f"m{k}"]))
                self.assertGreaterEqual(hi, float(exact[f"m{k}"]))


def _inverse(S):
    """Exact 2x2/3x3 inverse via adjugate, for the similarity test."""
    n = len(S)
    det = C.determinant(S)
    assert det != 0
    cof = [[((-1) ** (i + j)) * C.determinant(
        [[S[r][c] for c in range(n) if c != i] for r in range(n) if r != j])
        for j in range(n)] for i in range(n)]
    return [[cof[i][j] / det for j in range(n)] for i in range(n)]


@unittest.skipIf(_flint() is None, "python-flint not installed")
class AvailableConclusionsOnly(unittest.TestCase):
    def test_a_positive_definite_block_is_not_claimed_psd_from_moments(self):
        """The load-bearing refusal: necessary is not sufficient."""
        r = analyse([[F(2), F(0)], [F(0), F(3)]], observed_inertia=[2, 0, 0])
        q = _q(r, "moments_force_psd")
        self.assertEqual(q["status"], INSUFFICIENT)
        self.assertIsNone(q["answer"])
        self.assertIn("not sufficient", q["reason"])

    def test_an_indefinite_block_is_conclusively_refuted(self):
        """The direction that *is* available: a localizing obstruction proves it."""
        r = analyse([[F(2), F(0)], [F(0), F(-3)]], observed_inertia=[1, 1, 0])
        q = _q(r, "moments_force_psd")
        self.assertEqual(q["status"], CONCLUSIVE)
        self.assertIs(q["answer"], False)
        self.assertEqual(_q(r, "minimum_negative_eigenvalue_count")["minimum"], 1)

    def test_2x2_moments_determine_the_inertia(self):
        for signs, want in ([2, 3], [2, 0, 0]), ([2, -3], [1, 1, 0]), ([-2, -3], [0, 2, 0]):
            r = analyse(C.diagonal(signs), observed_inertia=want)
            q = _q(r, "inertia_determined_by_moments")
            with self.subTest(signs=tuple(signs)):
                self.assertEqual(q["status"], CONCLUSIVE)
                self.assertTrue(q["determined"])
                self.assertEqual(q["implied_inertia"], want)
                self.assertTrue(q["matches_observed"])

    def test_3x3_moments_do_not_determine_the_inertia(self):
        G = C.congruence(C.diagonal([1, 1, -1]), C.unimodular(3, 5))
        r = analyse(G, observed_inertia=list(exact_inertia(G).signature))
        q = _q(r, "inertia_determined_by_moments")
        self.assertEqual(q["status"], INSUFFICIENT)
        self.assertFalse(q["determined"])
        self.assertIn("not injective", q["reason"])

    def test_smallest_eigenvalue_bounds_bracket_the_truth(self):
        for signs in ([2, 3], [2, -3], [-2, -3], [1, 1, -4], [5, -1, 2]):
            r = analyse(C.diagonal(signs))
            q = _q(r, "smallest_eigenvalue_bounds")
            with self.subTest(signs=tuple(signs)):
                self.assertEqual(q["status"], CONCLUSIVE)
                lo = float(q["lambda_min"]["lo"])
                hi = float(q["lambda_min"]["hi"])
                self.assertLessEqual(lo, min(signs) + 1e-9)
                self.assertGreaterEqual(hi, min(signs) - 1e-9)
                self.assertLessEqual(float(q["lambda_max"]["lo"]), max(signs) + 1e-9)

    def test_bounds_are_tight_at_n_equals_2_and_loose_above(self):
        self.assertTrue(eigenvalue_bounds_from_two_moments(2, F(5), F(13))["tight"])
        self.assertFalse(eigenvalue_bounds_from_two_moments(3, F(5), F(13))["tight"])

    def test_interval_input_reaches_the_same_conclusions(self):
        arb = _flint()
        G = [[arb("2"), arb("0.5")], [arb("0.5"), arb("-3")]]
        r = analyse(G, observed_inertia=[1, 1, 0])
        self.assertEqual(_q(r, "moments_force_psd")["status"], CONCLUSIVE)
        q = _q(r, "inertia_determined_by_moments")
        self.assertEqual(q["implied_inertia"], [1, 1, 0])

    def test_no_report_claims_rh(self):
        r = analyse(C.diagonal([2, 3]))
        self.assertIs(r["rh_proof_claim"], False)
        self.assertEqual(r["claim_scope"], "finite_dimensional_weil_compression")


class B1IsTheMomentEngine(unittest.TestCase):
    @unittest.skipUnless(b1_available(), "b1_moment_solver not importable")
    def test_exact_moments_are_routed_through_b1(self):
        r = analyse(C.diagonal([F(2), F(-3)]))
        q = _q(r, "moments_force_psd")
        self.assertIn("b1_moment_solver", q["localizing_matrix_status"]["engine"])

    @unittest.skipUnless(b1_available(), "b1_moment_solver not importable")
    def test_b1_hankel_view_reports_rank_and_flatness(self):
        r = analyse(C.diagonal([F(2), F(3)]))
        v = r["b1_hankel_view"]
        self.assertTrue(v["available"])
        self.assertEqual(v["rank_M2"], 2)
        self.assertEqual(v["distinct_eigenvalues_at_least"], 2)

    @unittest.skipUnless(b1_available(), "b1_moment_solver not importable")
    def test_hankel_rank_lower_bounds_the_distinct_eigenvalue_count(self):
        for signs in ([2, 3], [1, 1], [1, 2, 3], [1, 1, 1]):
            r = analyse(C.diagonal(signs))
            v = r["b1_hankel_view"]
            with self.subTest(signs=tuple(signs)):
                self.assertLessEqual(v["distinct_eigenvalues_at_least"],
                                     len(set(signs)))

    @unittest.skipIf(_flint() is None, "python-flint not installed")
    def test_interval_moments_do_not_pretend_to_be_exact(self):
        arb = _flint()
        r = analyse([[arb("2"), arb("0")], [arb("0"), arb("3")]])
        self.assertFalse(r["b1_hankel_view"]["available"])
        self.assertIn("exact-rational", r["b1_hankel_view"]["reason"])


if __name__ == "__main__":
    unittest.main(verbosity=1)
