#!/usr/bin/env python3
"""ATLAS-RH-ENG-009 §WO-RH-57/58 — reference metric and generalized gap, exact half.

Everything here is E0 unless a test says otherwise: exact rational arithmetic,
no floating point on any load-bearing path. The two facts that make the whole
work order sound are proved *as identities on samples* here and as theorems in
the Lean layer:

1. the reference metric is what it says it is (the L^2 Gram of the basis) and
   is positive definite for every L > 0;
2. the generalized spectrum of the pencil (G, M) is invariant under
   simultaneous congruence, which is what raw eigenvalues are not.
"""
from __future__ import annotations

import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for extra in (ROOT, ROOT / "src"):
    sys.path.insert(0, str(extra))

import basis_algebra as BA  # noqa: E402
import generalized_gap as GG  # noqa: E402
import reference_metric as RM  # noqa: E402

LS = (F(7, 6), F(11, 8), F(6, 5))

BLOCKS = (("one",), ("q1",), ("one", "b"), ("q1", "b3"), ("one", "b", "b2"))


def h_eval(name, x, L):
    """A basis element at exact rational (x, L), straight off the primitive table."""
    total = F(0)
    for xpow, lpoly in enumerate(BA.BASIS_L_POLY[name]):
        for lpow, c in lpoly.items():
            total += c * x ** xpow * L ** lpow
    return total


class TheMetricIsWhatItSaysItIs(unittest.TestCase):
    """M_ij really is the L^2 inner product, checked against direct integration."""

    def test_against_exact_quadrature_of_the_defining_integral(self):
        # h_i h_j is a polynomial of degree <= 8, so Gauss-Legendre would do,
        # but exact power-sum integration is simpler: integrate monomials.
        for L in LS:
            for i in BA.BASIS_NAMES:
                for j in BA.BASIS_NAMES:
                    # exact integral via monomial expansion at this L
                    prod: dict = {}
                    for xi, li in enumerate(BA.BASIS_L_POLY[i]):
                        for xj, lj in enumerate(BA.BASIS_L_POLY[j]):
                            ci = sum(c * L ** k for k, c in li.items())
                            cj = sum(c * L ** k for k, c in lj.items())
                            prod[xi + xj] = prod.get(xi + xj, F(0)) + ci * cj
                    integral = sum(c * L ** (k + 1) / (k + 1)
                                   for k, c in prod.items())
                    self.assertEqual(RM.metric_exact(i, j, L), integral, (i, j, L))

    def test_every_entry_is_the_predicted_monomial(self):
        for i in BA.BASIS_NAMES:
            for j in BA.BASIS_NAMES:
                coeff, power = RM.metric_monomial(i, j)
                self.assertEqual(power,
                                 RM.BASIS_DEGREE[i] + RM.BASIS_DEGREE[j] + 1)

    def test_cross_parity_entries_vanish_exactly(self):
        for i in BA.BASIS_NAMES:
            for j in BA.BASIS_NAMES:
                if BA.BASIS_PARITY[i] != BA.BASIS_PARITY[j]:
                    self.assertEqual(RM.metric_monomial(i, j)[0], F(0), (i, j))

    def test_same_parity_diagonals_are_positive(self):
        for i in BA.BASIS_NAMES:
            self.assertGreater(RM.metric_monomial(i, i)[0], 0, i)


class TheMetricIsPositiveDefinite(unittest.TestCase):
    """The E0 certificate: exact Sylvester on M(1), congruence for L > 0."""

    def test_every_block_certifies(self):
        for basis in BLOCKS:
            rec = RM.certify_positive_definite(basis)
            self.assertEqual(rec["evidence_class"], "E0")
            for minor in rec["unit_leading_minors"]:
                self.assertGreater(F(minor), 0, basis)

    def test_the_certificate_recomputes_rather_than_stores(self):
        # The record's unit matrix must equal the derived one, entry for entry.
        rec = RM.certify_positive_definite(("one", "b", "b2"))
        derived = RM.unit_matrix(("one", "b", "b2"))
        for r, row in enumerate(rec["unit_matrix"]):
            for c, val in enumerate(row):
                self.assertEqual(F(val), derived[r][c])

    def test_positivity_at_exact_sample_points_follows(self):
        # Spot-check the conclusion the congruence licenses: v^T M(L) v > 0.
        vs = ((F(1), F(-3), F(2)), (F(0), F(1), F(-1)), (F(5), F(1), F(7)))
        for L in LS:
            m = RM.metric_matrix_exact(("one", "b", "b2"), L)
            for v in vs:
                q = sum(v[a] * m[a][b] * v[b] for a in range(3) for b in range(3))
                self.assertGreater(q, 0)


class TheGeneralizedSpectrumIsInvariant(unittest.TestCase):
    """det(S^T G S - lam S^T M S) = det(S)^2 det(G - lam M), exactly.

    This is §WO-RH-57's load-bearing claim -- the reason a generalized gap can
    be compared across bases when a raw eigenvalue cannot. It is checked as a
    polynomial identity in lam: all coefficients, not a sample of them.
    """

    G = [[F(2), F(1, 3), F(0)], [F(1, 3), F(5, 7), F(-1)], [F(0), F(-1), F(3)]]
    M = [[F(1), F(1, 6), F(1, 30)], [F(1, 6), F(1, 30), F(1, 140)],
         [F(1, 30), F(1, 140), F(1, 630)]]
    S = [[F(1), F(2), F(0)], [F(0), F(1, 2), F(3)], [F(1), F(0), F(1)]]

    @staticmethod
    def _congruence(S, A):
        n = len(A)
        SA = [[sum(S[k][a] * A[k][l] for k in range(n)) for l in range(n)]
              for a in range(n)]
        return [[sum(SA[a][k] * S[k][b] for k in range(n)) for b in range(n)]
                for a in range(n)]

    @staticmethod
    def _pencil_poly(G, M):
        """Coefficients of det(G - lam M) as an exact polynomial in lam."""
        n = len(G)
        # Represent entries as polynomials in lam: (g, -m).
        entries = [[(G[a][b], -M[a][b]) for b in range(n)] for a in range(n)]

        def poly_mul(p, q):
            out = [F(0)] * (len(p) + len(q) - 1)
            for i, a in enumerate(p):
                for j, b in enumerate(q):
                    out[i + j] += a * b
            return out

        def det(rows):
            if len(rows) == 1:
                return list(rows[0][0])
            total = None
            for col in range(len(rows)):
                minor = [r[:col] + r[col + 1:] for r in rows[1:]]
                term = poly_mul(rows[0][col], det(minor))
                if col % 2:
                    term = [-c for c in term]
                if total is None:
                    total = term
                else:
                    total = [a + b for a, b in
                             zip(total + [F(0)] * len(term),
                                 term + [F(0)] * len(total))][:max(len(total), len(term))]
            return total

        return det(entries)

    @staticmethod
    def _det3(A):
        return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))

    def test_the_pencil_polynomial_transforms_by_det_squared(self):
        Gp = self._congruence(self.S, self.G)
        Mp = self._congruence(self.S, self.M)
        p = self._pencil_poly(self.G, self.M)
        q = self._pencil_poly(Gp, Mp)
        scale = self._det3(self.S) ** 2
        self.assertNotEqual(scale, 1)  # a congruence that actually moves things
        self.assertEqual(len(p), len(q))
        for a, b in zip(p, q):
            self.assertEqual(b, a * scale)

    def test_raw_eigenvalues_are_not_invariant_and_the_module_says_so(self):
        # The contrast that motivates the pencil: trace (sum of raw eigenvalues)
        # moves under the same congruence.
        Gp = self._congruence(self.S, self.G)
        tr = sum(self.G[a][a] for a in range(3))
        trp = sum(Gp[a][a] for a in range(3))
        self.assertNotEqual(tr, trp)
        self.assertIn("basis", GG.__doc__)


class ShiftedPositivityMachinery(unittest.TestCase):
    """The runtime pieces the certifier composes, on exact carriers."""

    def test_shifted_matrix_is_g_minus_lambda_m(self):
        basis = ("one", "b")
        L = F(6, 5)
        ent = {("one", "one"): F(2), ("one", "b"): F(1, 4), ("b", "b"): F(1, 9)}
        lam = F(1, 3)
        m = GG.shifted_matrix(basis, ent, lam, L)
        self.assertEqual(m[0][0], F(2) - lam * RM.metric_exact("one", "one", L))
        self.assertEqual(m[0][1], F(1, 4) - lam * RM.metric_exact("one", "b", L))
        self.assertEqual(m[1][1], F(1, 9) - lam * RM.metric_exact("b", "b", L))

    def test_leading_minors_match_the_even3_implementation(self):
        import even3 as E3
        m = [[F(4), F(1), F(0)], [F(1), F(3), F(-1)], [F(0), F(-1), F(2)]]
        self.assertEqual(GG.leading_minors(m), E3.leading_minors(m))

    def test_preconditioning_scales_minors_by_exact_powers_of_two(self):
        m = [[4.0, 1.0, 0.5], [1.0, 3.0, -1.0], [0.5, -1.0, 2.0]]
        pre = GG.precondition(m, (-1, -2, -3))
        raw = GG.leading_minors(m)
        scaled = GG.leading_minors(pre)
        for k, (r, s) in enumerate(zip(raw, scaled), start=1):
            factor = 4.0 ** (-sum(range(1, k + 1)) * 0 - sum([1, 2, 3][:k]))
            self.assertAlmostEqual(s, r * factor, places=12, msg=str(k))

    def test_rayleigh_upper_is_exact_on_rational_data(self):
        basis = ("one", "b")
        L = F(6, 5)
        ent = {("one", "one"): F(2), ("one", "b"): F(1, 4), ("b", "b"): F(1, 9)}
        v = (F(1), F(-2))
        got = GG.rayleigh_upper(basis, ent, v, L)
        num = (F(2) - 2 * F(1, 4) * 2 + 4 * F(1, 9))
        den = (RM.metric_exact("one", "one", L)
               - 4 * RM.metric_exact("one", "b", L)
               + 4 * RM.metric_exact("b", "b", L))
        self.assertEqual(got, num / den)

    def test_scout_agrees_with_the_exact_crossing_on_a_solvable_case(self):
        # 1x1: G - lam M crosses zero at exactly G / M.
        basis = ("one",)
        ent = {("one", "one"): 0.07}
        lam = GG.scout_gap_at(basis, ent, 1.25)
        self.assertAlmostEqual(lam, 0.07 / 1.25, places=9)


class NoClaims(unittest.TestCase):
    def test_the_modules_state_their_boundary(self):
        self.assertIn("No RH proof claim", RM.__doc__)
        self.assertIn("No RH proof claim", GG.__doc__)

    def test_the_gap_kind_is_not_psd_licensable(self):
        import content_kinds as CK
        for kind in (CK.KIND_GENERALIZED_GAP, CK.KIND_STRUCTURAL_DIAGNOSTIC,
                     CK.KIND_SCALING_MODEL, CK.KIND_NEXT_BLOCK_SELECTION):
            self.assertTrue(CK.is_registered(kind))
            self.assertFalse(CK.psd_licensable(kind))


if __name__ == "__main__":
    unittest.main(verbosity=2)
