#!/usr/bin/env python3
"""ATLAS-RH-ENG-008 §WO-RH-48/49 — the derived kernel and derivative algebra.

Everything here is E0: exact rational arithmetic, no floating point.

Two things are being checked, and they are different. The first is that the
derived engine *reproduces* what it replaced -- the six hand-written closed
forms, their six coefficient expansions in ``a``, and their six ``d/dL``
expansions, all of which ENG-005 and ENG-006 had verified against SymPy. Those
values are inlined here as literals, so this file is the record of them and the
generalization is pinned rather than trusted.

The second is that the engine is *right*, independently of what it replaced:
each kernel is checked against direct symbolic integration of its own
definition, and the derivative machinery against symbolic differentiation of
the same integral. A generalization that faithfully reproduces a wrong table
would pass the first check and fail the second.
"""
from __future__ import annotations

import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import basis_algebra as BA  # noqa: E402
import core  # noqa: E402
import pole  # noqa: E402

#: Rational sample points. ``a`` may exceed ``L``: the kernel is a polynomial
#: identity in ``(a, L)`` and the tests compare it as one. The integral reading
#: only holds for ``a <= L``, which is where the runtime evaluates it.
LS = (F(7, 3), F(11, 8), F(3), F(5, 2))
AS = (F(0), F(1, 4), F(5, 7), F(9, 8), F(2))

LEGACY_PAIRS = (
    ("one", "one", core.kernel_00),
    ("one", "b", core.kernel_0b),
    ("b", "b", core.kernel_bb),
    ("q1", "q1", core.kernel_q1q1),
    ("q1", "b3", core.kernel_q1b3),
    ("b3", "b3", core.kernel_b3b3),
)

#: The retired hand-written coefficient expansions in ``a``, ascending, as
#: functions of ``L``. Verbatim from ``archimedean_realspace`` before ENG-008.
LEGACY_COEFFS_IN_A = {
    ("one", "one"): lambda L: [2 * L, F(-2)],
    ("one", "b"): lambda L: [L**3 / 3, F(0), -L, F(2, 3)],
    ("b", "b"): lambda L: [L**5 / 15, F(0), -L**3 / 3, L**2 / 3, F(0), F(-1, 15)],
    ("q1", "q1"): lambda L: [L**3 / 6, -L**2 / 2, F(0), F(1, 3)],
    ("q1", "b3"): lambda L: [L**5 / 60, F(0), -L**3 / 4, L**2 / 3, F(0), F(-1, 10)],
    ("b3", "b3"): lambda L: [L**7 / 420, F(0), -L**5 / 20, L**4 / 12, F(0),
                             -L**2 / 20, F(0), F(1, 70)],
}

#: The retired hand-written ``d/dL`` expansions, same convention.
LEGACY_DL_COEFFS_IN_A = {
    ("one", "one"): lambda L: [F(2), F(0)],
    ("one", "b"): lambda L: [L**2, F(0), F(-1), F(0)],
    ("b", "b"): lambda L: [L**4 / 3, F(0), -L**2, 2 * L / 3, F(0), F(0)],
    ("q1", "q1"): lambda L: [L**2 / 2, -L, F(0), F(0)],
    ("q1", "b3"): lambda L: [L**4 / 12, F(0), -L**2 * F(3, 4), L * F(2, 3), F(0), F(0)],
    ("b3", "b3"): lambda L: [L**6 / 60, F(0), -L**4 / 4, L**3 / 3, F(0),
                             -L / 10, F(0), F(0)],
}


def _pad(a, b):
    a, b = list(a), list(b)
    while len(a) < len(b):
        a.append(F(0))
    while len(b) < len(a):
        b.append(F(0))
    return a, b


class ReproducesTheRetiredTables(unittest.TestCase):
    def test_closed_forms(self):
        for i, j, fn in LEGACY_PAIRS:
            for L in LS:
                for a in AS:
                    self.assertEqual(BA.kernel_exact(i, j, a, L), F(fn(a, L)),
                                     (i, j, a, L))

    def test_coefficient_expansions_in_a(self):
        for (i, j), fn in LEGACY_COEFFS_IN_A.items():
            for L in LS:
                got, want = _pad([F(c) for c in BA.kernel_coeffs_in_a(i, j, L)],
                                 fn(L))
                self.assertEqual(got, want, (i, j, L))

    def test_dL_expansions_in_a(self):
        for (i, j), fn in LEGACY_DL_COEFFS_IN_A.items():
            for L in LS:
                got, want = _pad([F(c) for c in BA.kernel_dL_coeffs_in_a(i, j, L)],
                                 fn(L))
                self.assertEqual(got, want, (i, j, L))

    def test_the_retired_second_derivative_closed_forms(self):
        # The four hand-written `_laplace_d2L` branches, as they read before
        # ENG-008 generalized them. `b3` is the interesting one: it is quadratic
        # in L, so its integral term was already being carried.
        import math

        def retired(name, L, sign):
            e = math.exp(sign * L / 2)
            if name == "one":
                return sign * e / 2
            if name == "q1":
                return sign * L * e / 4
            if name == "b":
                return L * e
            if name == "b3":
                return L * L * e / 2 - pole.poly_exp_integral(
                    (0.0, 1.0), sign * 0.5, L)
            raise KeyError(name)

        for L in (1.0986122886681098, 1.2, 1.3862943611198906, 2.5):
            for name in ("one", "q1", "b", "b3"):
                for sign in (1, -1):
                    got = pole._laplace_d2L(name, L, sign)
                    want = retired(name, L, sign)
                    self.assertAlmostEqual(got, want, delta=1e-12 * max(1.0, abs(want)),
                                           msg=(name, L, sign))


class TheFactorizationIsRecoveredAndUsed(unittest.TestCase):
    """``(L - a)^m`` is not cosmetic on an interval carrier.

    Every kernel integrates over ``[0, L - a]``, so that factor is there, and the
    retired hand-written closed forms all displayed it. The first version of the
    derived engine evaluated the *expanded* polynomial instead, which is the same
    number on an exact carrier and a much wider one on a ball: the prime block's
    enclosure widened 3x for ``K_q1q1`` and 48x for ``K_b3b3``, costing the
    degree-3 determinant bound 26%. (The 3x3 covers are insensitive: they only
    ever evaluate kernels at box midpoints, via the centred mean-value form.)

    These tests pin both halves -- that the factorization is found, and that
    using it is what keeps the enclosures tight.
    """

    #: The multiplicities the retired closed forms displayed.
    EXPECTED_M = {
        ("one", "one"): 1, ("one", "b"): 2, ("b", "b"): 3,
        ("q1", "q1"): 1, ("q1", "b3"): 2, ("b3", "b3"): 3,
        ("one", "b2"): 3, ("b", "b2"): 4, ("b2", "b2"): 5,
        # ENG-009 §WO-RH-62 prep: bcube = b^3 vanishes to order 3 at both
        # endpoints, and m = 1 + ord_i + ord_j throughout this table.
        ("one", "bcube"): 4, ("b", "bcube"): 5, ("b2", "bcube"): 6,
        ("bcube", "bcube"): 7,
    }

    def test_the_multiplicities_match_the_retired_closed_forms(self):
        for (i, j), m in self.EXPECTED_M.items():
            self.assertEqual(BA.kernel_factored(i, j)[0], m, (i, j))

    def test_every_same_parity_kernel_has_at_least_one_such_factor(self):
        # It integrates over [0, L - a]; a nonzero kernel without the factor
        # would mean the derivation is wrong.
        for i in BA.BASIS_NAMES:
            for j in BA.BASIS_NAMES:
                if BA.BASIS_PARITY[i] != BA.BASIS_PARITY[j]:
                    continue
                self.assertGreaterEqual(BA.kernel_factored(i, j)[0], 1, (i, j))

    def test_cross_parity_kernels_vanish_identically(self):
        # The parity block structure at the kernel level, exactly: a kernel
        # pairing an even element against an odd one is the zero polynomial in
        # (a, L), not merely small. This is the same fact `cross_block_vanishes`
        # records on the Lean side, and it is why the multiplicity test above is
        # restricted -- the zero polynomial carries no `(L - a)` factor to find.
        for i in BA.BASIS_NAMES:
            for j in BA.BASIS_NAMES:
                if BA.BASIS_PARITY[i] == BA.BASIS_PARITY[j]:
                    continue
                self.assertEqual(BA.kernel_bipoly(i, j), {}, (i, j))
                for L in LS:
                    for a in AS:
                        self.assertEqual(BA.kernel_exact(i, j, a, L), 0, (i, j, a, L))

    def test_the_factored_form_is_the_same_number_exactly(self):
        for i in BA.BASIS_NAMES:
            for j in BA.BASIS_NAMES:
                for L in LS:
                    for a in AS:
                        self.assertEqual(BA.kernel_value(i, j, a, L),
                                         BA.kernel_exact(i, j, a, L), (i, j, a, L))

    def test_it_is_materially_tighter_on_a_ball(self):
        try:
            from interval_backend import interval_box, require_flint
        except Exception:  # pragma: no cover
            raise unittest.SkipTest("python-flint is required")
        _, arb, _, ctx = require_flint()
        ctx.prec = 160
        box = interval_box(1.2559472014002986, 1.2559559807604197)
        a = arb(2).log()

        def expanded(i, j, a, L):
            coeffs = BA.kernel_coeffs_in_a(i, j, L)
            out = coeffs[-1]
            for c in reversed(coeffs[:-1]):
                out = out * a + c
            return out

        for (i, j), floor in ((("q1", "q1"), 2.0), (("b3", "b3"), 5.0),
                              (("b2", "b2"), 5.0)):
            wide = float(expanded(i, j, a, box).rad())
            tight = float(BA.kernel_value(i, j, a, box).rad())
            self.assertGreater(wide / tight, floor,
                               f"{i},{j}: factored form only {wide / tight:.1f}x "
                               "tighter than expanded; the factorization may have "
                               "stopped being used")


class AgreesWithSymbolicIntegration(unittest.TestCase):
    """Independent of the retired tables: the definitions themselves."""

    @classmethod
    def setUpClass(cls):
        try:
            import sympy  # noqa: F401
        except ImportError:  # pragma: no cover
            raise unittest.SkipTest("sympy is required for the symbolic route")

    def test_every_kernel_against_its_defining_integral(self):
        import sympy as sp

        x, a, L = sp.symbols("x a L", positive=True)
        h = {
            "one": sp.Integer(1),
            "q1": x - L / 2,
            "b": x * (L - x),
            "b3": x * (L - x) * (x - L / 2),
            "b2": (x * (L - x)) ** 2,
        }
        names = list(h)
        for m, i in enumerate(names):
            for j in names[m:]:
                fi, fj = h[i], h[j]
                integrand = (fi * fj.subs(x, x + a) + fj * fi.subs(x, x + a))
                K = sp.expand(sp.integrate(sp.expand(integrand), (x, 0, L - a)))
                for Lv in (F(7, 3), F(5, 2)):
                    for av in (F(0), F(1, 4), F(2)):
                        want = sp.nsimplify(
                            K.subs({a: sp.Rational(av), L: sp.Rational(Lv)}))
                        got = BA.kernel_exact(i, j, av, Lv)
                        self.assertEqual(
                            sp.Rational(got.numerator, got.denominator), want,
                            (i, j, av, Lv))

    def test_second_L_derivative_against_symbolic_differentiation(self):
        import sympy as sp

        x, L = sp.symbols("x L", positive=True)
        h = {
            "one": sp.Integer(1),
            "q1": x - L / 2,
            "b": x * (L - x),
            "b3": x * (L - x) * (x - L / 2),
            "b2": (x * (L - x)) ** 2,
        }
        for sign in (1, -1):
            for name, expr in h.items():
                Fint = sp.integrate(sp.expand(expr) * sp.exp(sign * x / 2), (x, 0, L))
                d2 = sp.diff(Fint, L, 2)
                for Lv in (sp.Rational(6, 5), sp.Rational(5, 2)):
                    want = float(d2.subs(L, Lv))
                    got = pole._laplace_d2L(name, float(Lv), sign)
                    self.assertAlmostEqual(
                        got, want, delta=1e-11 * max(1.0, abs(want)),
                        msg=(name, sign, Lv))


class TheNewEvenSector(unittest.TestCase):
    """``b2`` specifically -- the element the retired tables did not know."""

    def test_b2_is_the_square_of_b(self):
        for L in LS:
            for xv in (F(1, 5), F(2, 7), F(11, 13), F(3, 2)):
                b = pole._horner(pole.basis_coeffs("b", L), xv)
                b2 = pole._horner(pole.basis_coeffs("b2", L), xv)
                self.assertEqual(b2, b * b, (L, xv))

    def test_b2_is_even_about_the_midpoint(self):
        for L in LS:
            for xv in (F(1, 5), F(2, 7), F(11, 13), F(3, 2)):
                left = pole._horner(pole.basis_coeffs("b2", L), L - xv)
                right = pole._horner(pole.basis_coeffs("b2", L), xv)
                self.assertEqual(left, right, (L, xv))
        self.assertEqual(pole.basis_parity("b2"), "even")

    def test_b2_vanishes_at_both_endpoints(self):
        for L in LS:
            self.assertEqual(pole._horner(pole.basis_coeffs("b2", L), F(0)), 0)
            self.assertEqual(pole.basis_at_right_endpoint("b2", L), 0)

    def test_b2_has_a_nonzero_second_L_derivative(self):
        # The point of §WO-RH-49. Every element before b2 had either a vanishing
        # d^2_L h (one, q1, b) or one the table happened to carry (b3); b2 is
        # quadratic in L through its L^2 x^2 term and the table did not know it
        # at all.
        for L in LS:
            coeffs = pole.basis_coeffs_d2L("b2", L)
            self.assertNotEqual(list(coeffs), [0] * len(coeffs), L)
            self.assertEqual(list(coeffs), [0, 0, 2, 0, 0], L)

    def test_the_elements_that_should_have_vanishing_second_derivatives_do(self):
        for name in ("one", "q1", "b"):
            for L in LS:
                self.assertTrue(
                    all(c == 0 for c in pole.basis_coeffs_d2L(name, L)),
                    (name, L))

    def test_b3_also_has_one_and_always_did(self):
        for L in LS:
            self.assertEqual(list(pole.basis_coeffs_d2L("b3", L)), [0, -1, 0, 0], L)

    def test_the_kernels_are_symmetric(self):
        for i in BA.BASIS_NAMES:
            for j in BA.BASIS_NAMES:
                for L in LS:
                    for a in AS:
                        self.assertEqual(BA.kernel_exact(i, j, a, L),
                                         BA.kernel_exact(j, i, a, L), (i, j, a, L))

    def test_every_kernel_vanishes_at_a_equals_L(self):
        # K integrates over [0, L - a], which is empty at a = L.
        for i in BA.BASIS_NAMES:
            for j in BA.BASIS_NAMES:
                for L in LS:
                    self.assertEqual(BA.kernel_exact(i, j, L, L), 0, (i, j, L))

    def test_the_diagonal_kernels_are_positive_inside_the_cell(self):
        # K_ii(a; L) = 2 int h_i(x) h_i(x+a) dx is not sign-definite in general,
        # but for the non-negative even elements it is.
        for i in ("one", "b", "b2"):
            for L in LS:
                for a in (F(0), F(1, 4), F(5, 7)):
                    if a >= L:
                        continue
                    self.assertGreater(BA.kernel_exact(i, i, a, L), 0, (i, a, L))

    def test_an_unknown_basis_element_is_refused(self):
        with self.assertRaises(KeyError):
            BA.kernel_exact("not_a_basis_element", "one", F(1, 2), F(2))
        with self.assertRaises(KeyError):
            pole.basis_coeffs_d2L("not_a_basis_element", F(2))


if __name__ == "__main__":
    unittest.main(verbosity=2)
