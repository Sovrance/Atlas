"""Canonical Candidate-A pole primitive and its exact invariants (ENG-004 §2).

The invariants are the ones the adjudication turned on:

* the midpoint parity relation ``E^- = ± e^{-L/2} E^+``;
* vanishing even/odd pole cross terms;
* the odd pivot ``G0[q1,q1] = -8 A(L)^2``;
* agreement with the independent real-space representation using ``2cosh(a/2)``;
* Candidate B disagreeing away from ``L = log 3``.

Every one is checked at the adjudication points **and** at an out-of-cell point,
because a calibration fitted at a single ``L`` is exactly what WO-RH-17 rejected.
"""
from __future__ import annotations

import math
import os
import sys
import unittest
from pathlib import Path

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import pole  # noqa: E402
import normalization as N  # noqa: E402
import rejected_pole as RP  # noqa: E402  (archival; this is an audit)

LOG3 = math.log(3.0)
LOG4 = math.log(4.0)

#: Adjudication points plus an out-of-cell point (ENG-004 §2).
POINTS = (LOG3, 1.1059498113, 1.20, LOG4, 3.5)
OUT_OF_CELL = 3.5

BASIS = ("one", "q1", "b")


def _real_space_pole(i: str, j: str, L: float, n: int = 4000) -> float:
    """``int_0^L K_ij(a) 2cosh(a/2) da`` — the independent representation.

    Uses the *same* ``K_ij`` as the prime block, so agreement ties the pole and
    prime halves of the assembly to one kernel. Composite Simpson; the tolerance
    below is set well above its truncation error.
    """
    h = L / n
    total = 0.0
    for k in range(n + 1):
        a = k * h
        w = 1 if k in (0, n) else (4 if k % 2 else 2)
        total += w * N.kernel_K(i, j, a, L) * 2 * math.cosh(a / 2)
    return total * h / 3


class PoleFormulaTests(unittest.TestCase):
    def test_module_declares_the_adopted_candidate(self):
        self.assertEqual(pole.POLE_CANDIDATE, "A")
        self.assertEqual(pole.POLE_STATUS, "ADOPTED_WO_RH_17")

    def test_entry_is_the_explicit_formula_product(self):
        for L in POINTS:
            for i in BASIS:
                for j in BASIS:
                    want = (pole.laplace_plus(i, L) * pole.laplace_minus(j, L)
                            + pole.laplace_minus(i, L) * pole.laplace_plus(j, L))
                    got = pole.pole_gram_entry(i, j, L)
                    self.assertAlmostEqual(got, want, delta=1e-12 * max(1.0, abs(want)))

    def test_matrix_matches_entrywise(self):
        for L in POINTS:
            m = pole.pole_gram_matrix(BASIS, L)
            for a, i in enumerate(BASIS):
                for b, j in enumerate(BASIS):
                    self.assertAlmostEqual(
                        m[a][b], pole.pole_gram_entry(i, j, L),
                        delta=1e-12 * max(1.0, abs(m[a][b])))

    def test_normalization_delegates_rather_than_duplicating(self):
        for L in POINTS:
            for i in BASIS:
                for j in BASIS:
                    self.assertEqual(N.pole_entry(i, j, L), pole.pole_gram_entry(i, j, L))


class ParityInvariants(unittest.TestCase):
    def test_midpoint_parity_relation(self):
        """``h(L-x) = ±h(x)`` forces ``E^- = ± e^{-L/2} E^+``."""
        for L in POINTS:
            damp = math.exp(-L / 2)
            for name in BASIS:
                sign = -1.0 if pole.basis_parity(name) == "odd" else 1.0
                lhs = pole.laplace_minus(name, L)
                rhs = sign * damp * pole.laplace_plus(name, L)
                self.assertAlmostEqual(lhs, rhs, delta=1e-12 * max(1.0, abs(rhs)), msg=(name, L))

    def test_even_odd_cross_terms_vanish(self):
        """The pole matrix is parity block diagonal."""
        for L in POINTS:
            for even in ("one", "b"):
                scale = max(1.0, abs(pole.pole_gram_entry(even, even, L)))
                self.assertLess(abs(pole.pole_gram_entry(even, "q1", L)), 1e-12 * scale, (even, L))
                self.assertLess(abs(pole.pole_gram_entry("q1", even, L)), 1e-12 * scale, (even, L))

    def test_odd_pivot_is_minus_eight_A_squared(self):
        """``A(L) = L cosh(L/4) - 4 sinh(L/4)`` and ``G0[q1,q1] = -8A^2``."""
        for L in POINTS:
            A = L * math.cosh(L / 4) - 4 * math.sinh(L / 4)
            self.assertAlmostEqual(
                pole.pole_gram_entry("q1", "q1", L), -8 * A * A,
                delta=1e-11 * max(1.0, A * A), msg=L)


class RealSpaceAgreement(unittest.TestCase):
    def test_candidate_A_matches_the_real_space_route(self):
        """Independent representation, same ``K_ij`` as the prime block."""
        for L in (LOG3, 1.20, LOG4, OUT_OF_CELL):
            for i, j in (("one", "one"), ("one", "b"), ("b", "b"), ("q1", "q1")):
                closed = pole.pole_gram_entry(i, j, L)
                direct = _real_space_pole(i, j, L)
                self.assertAlmostEqual(closed, direct, delta=1e-7 * max(1.0, abs(closed)),
                                       msg=(i, j, L))

    def test_out_of_cell_point_is_covered(self):
        """Guards against an invariant that only holds on the research cell."""
        self.assertGreater(OUT_OF_CELL, LOG4)


class CandidateBRejection(unittest.TestCase):
    def test_ratio_is_sqrt3_over_2_cosh(self):
        for L in POINTS:
            for i, j in (("one", "one"), ("one", "b"), ("b", "b")):
                ratio = RP.legacy_pole_entry(i, j, L) / pole.pole_gram_entry(i, j, L)
                self.assertAlmostEqual(ratio, RP.legacy_over_adopted_ratio(L), delta=1e-12,
                                       msg=(i, j, L))

    def test_agrees_only_at_log3(self):
        self.assertAlmostEqual(RP.legacy_over_adopted_ratio(LOG3), 1.0, delta=1e-14)
        for L in (1.1059498113, 1.20, LOG4, OUT_OF_CELL):
            self.assertGreater(abs(RP.legacy_over_adopted_ratio(L) - 1.0), 1e-4, L)

    def test_the_factor_is_L_dependent_so_not_a_change_of_basis(self):
        """A change of basis is a *constant* congruence."""
        r1 = RP.legacy_over_adopted_ratio(1.20)
        r2 = RP.legacy_over_adopted_ratio(LOG4)
        self.assertGreater(abs(r1 - r2), 1e-3)


class PoleDerivative(unittest.TestCase):
    """``pole_gram_entry_dL`` replaced a second, rejected copy in the jet module."""

    def test_matches_central_finite_difference(self):
        h = 1e-6
        for L in (LOG3, 1.20, LOG4, OUT_OF_CELL):
            for i in BASIS:
                for j in BASIS:
                    an = pole.pole_gram_entry_dL(i, j, L)
                    fd = (pole.pole_gram_entry(i, j, L + h)
                          - pole.pole_gram_entry(i, j, L - h)) / (2 * h)
                    scale = max(1.0, abs(pole.pole_gram_entry(i, j, L)))
                    self.assertLess(abs(an - fd), 1e-6 * scale, msg=(i, j, L))

    def test_jet_module_uses_the_primitive(self):
        import weil_fourier_jets

        src = Path(weil_fourier_jets.__file__).read_text(encoding="utf-8")
        self.assertIn("pole.pole_gram_entry_dL", src)


class IntervalCarrier(unittest.TestCase):
    def test_arb_enclosure_is_outward(self):
        try:
            from flint import arb, ctx
        except ImportError:
            self.skipTest("python-flint unavailable")
        ctx.prec = 256
        box = arb("1.20", "0.05")
        enclosure = pole.pole_gram_entry("one", "one", box, backend="flint")
        for t in (1.15, 1.20, 1.25):
            v = float(pole.pole_gram_entry("one", "one", arb(repr(t))))
            self.assertLessEqual(float(enclosure.lower()), v)
            self.assertLessEqual(v, float(enclosure.upper()))

    def test_flint_backend_rejects_a_float_carrier(self):
        from interval_backend import FlintUnavailable

        with self.assertRaises(FlintUnavailable):
            pole.pole_gram_entry("one", "one", 1.2, backend="flint")


if __name__ == "__main__":
    unittest.main()
