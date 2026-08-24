#!/usr/bin/env python3
"""ATLAS-RH-ENG-006 §7 — exact degree-3 identities (WO-RH-32), evidence class E0.

The two kernels the work order supplies are re-derived from the basis functions
with SymPy on every run rather than trusted as transcribed constants::

    K_q1b3 = (L-a)^2 (L^3 + 2L^2 a - 12L a^2 - 6a^3) / 60
    K_b3b3 = (L-a)^3 (L^4 + 3L^3 a - 15L^2 a^2 - 18L a^3 - 6a^4) / 420

A transcription error in either would be invisible downstream: the assembled
block would still be a smooth positive-looking 2x2, and every certificate built
on it would be internally consistent and wrong. Re-deriving is cheap; not
re-deriving is the kind of thing that survives to publication.

The kernel convention matters and is checked explicitly. The repository's
``K_ij`` is the *symmetrized* overlap, covering both the ``+log q`` and
``-log q`` shifts the explicit formula sums over. Taking the one-sided integral
instead gives exactly half, which would be a silent factor-2 error throughout
the prime block; the test pins the convention by also re-deriving ENG-005's
``K_q1q1``, whose value is already fixed by merged work.
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

try:
    import sympy as sp
except ImportError:  # pragma: no cover
    sp = None

import core  # noqa: E402
import degree3 as D3  # noqa: E402
import pole  # noqa: E402
import weil_entries as WE  # noqa: E402
from inertia.congruence import inertia_2x2  # noqa: E402


def _flint():
    try:
        from flint import arb, ctx

        ctx.prec = 160
        return arb
    except ImportError:  # pragma: no cover
        return None


@unittest.skipIf(sp is None, "sympy is required for the exact identities")
class KernelIdentities(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x, cls.a, cls.L = sp.symbols("x a L", positive=True)
        x, L = cls.x, cls.L
        cls.basis = {
            "one": sp.Integer(1),
            "q1": x - L / 2,
            "b": x * (L - x),
            "b3": x * (L - x) * (x - L / 2),
        }

    def _Ksym(self, i, j):
        """The repository convention: both +a and -a shifts, as the formula sums."""
        x, a, L = self.x, self.a, self.L
        f, g = self.basis[i], self.basis[j]
        return sp.expand(sp.integrate(f * g.subs(x, x + a) + f.subs(x, x + a) * g,
                                      (x, 0, L - a)))

    def test_K_q1b3_matches_the_work_order(self):
        a, L = self.a, self.L
        want = (L - a) ** 2 * (L**3 + 2 * L**2 * a - 12 * L * a**2 - 6 * a**3) / 60
        self.assertEqual(sp.simplify(self._Ksym("q1", "b3") - want), 0)

    def test_K_b3b3_matches_the_work_order(self):
        a, L = self.a, self.L
        want = ((L - a) ** 3
                * (L**4 + 3 * L**3 * a - 15 * L**2 * a**2 - 18 * L * a**3 - 6 * a**4)
                / 420)
        self.assertEqual(sp.simplify(self._Ksym("b3", "b3") - want), 0)

    def test_the_symmetrized_convention_reproduces_ENG005_K_q1q1(self):
        """Pins the factor of 2: the one-sided integral would give half of this."""
        a, L = self.a, self.L
        want = (L - a) * (L**2 - 2 * L * a - 2 * a**2) / 6
        self.assertEqual(sp.simplify(self._Ksym("q1", "q1") - want), 0)
        one_sided = sp.expand(sp.integrate(
            self.basis["q1"] * self.basis["q1"].subs(self.x, self.x + a),
            (self.x, 0, L - a)))
        self.assertEqual(sp.simplify(one_sided - want / 2), 0)

    def test_the_shipped_kernel_functions_agree_with_the_derivation(self):
        a, L = self.a, self.L
        for i, j, fn in (("q1", "b3", core.kernel_q1b3), ("b3", "b3", core.kernel_b3b3),
                         ("q1", "q1", core.kernel_q1q1)):
            with self.subTest(pair=(i, j)):
                self.assertEqual(sp.simplify(sp.expand(fn(a, L)) - self._Ksym(i, j)), 0)

    def test_the_kernel_table_routes_both_orderings(self):
        a, L = self.a, self.L
        self.assertEqual(sp.simplify(WE.kernel("q1", "b3", a, L)
                                     - WE.kernel("b3", "q1", a, L)), 0)

    def test_u_expansions_match_the_closed_forms(self):
        """The interval path uses polynomial coefficient lists, not the closed form."""
        import archimedean_realspace as AR

        a, L = self.a, self.L
        for i, j in (("q1", "q1"), ("q1", "b3"), ("b3", "b3")):
            with self.subTest(pair=(i, j)):
                coeffs = AR.kernel_coeffs_in_u(i, j, L, lambda v: sp.sympify(v))
                rebuilt = sum(c * a**k for k, c in enumerate(coeffs))
                self.assertEqual(sp.simplify(sp.expand(rebuilt - self._Ksym(i, j))), 0)

    def test_dL_expansions_are_the_derivative_of_the_kernels(self):
        import archimedean_realspace as AR

        a, L = self.a, self.L
        for i, j in (("q1", "q1"), ("q1", "b3"), ("b3", "b3")):
            with self.subTest(pair=(i, j)):
                coeffs = AR.kernel_coeffs_dL_in_u(i, j, L, lambda v: sp.sympify(v))
                rebuilt = sum(c * a**k for k, c in enumerate(coeffs))
                want = sp.diff(self._Ksym(i, j), L)
                self.assertEqual(sp.simplify(sp.expand(rebuilt - want)), 0)


@unittest.skipIf(sp is None, "sympy is required for the exact identities")
class BasisAndParity(unittest.TestCase):
    def test_pole_basis_coefficients_rebuild_the_basis_functions(self):
        x, L = sp.symbols("x L", positive=True)
        want = {"one": sp.Integer(1), "q1": x - L / 2, "b": x * (L - x),
                "b3": x * (L - x) * (x - L / 2)}
        for name, h in want.items():
            with self.subTest(name):
                coeffs = pole.basis_coeffs(name, L)
                rebuilt = sum(sp.expand(c) * x**k for k, c in enumerate(coeffs))
                self.assertEqual(sp.simplify(sp.expand(rebuilt - h)), 0)

    def test_pole_dL_coefficients_are_the_L_derivative(self):
        x, L = sp.symbols("x L", positive=True)
        want = {"one": sp.Integer(1), "q1": x - L / 2, "b": x * (L - x),
                "b3": x * (L - x) * (x - L / 2)}
        for name, h in want.items():
            with self.subTest(name):
                coeffs = pole.basis_coeffs_dL(name, L)
                rebuilt = sum(sp.expand(c) * x**k for k, c in enumerate(coeffs))
                self.assertEqual(sp.simplify(sp.expand(rebuilt - sp.diff(h, L))), 0)

    def test_b3_vanishes_at_both_endpoints(self):
        x, L = sp.symbols("x L", positive=True)
        b3 = x * (L - x) * (x - L / 2)
        self.assertEqual(sp.simplify(b3.subs(x, 0)), 0)
        self.assertEqual(sp.simplify(b3.subs(x, L)), 0)
        self.assertEqual(sp.simplify(pole.basis_at_right_endpoint("b3", L)), 0)

    def test_parity_about_the_midpoint_splits_the_basis(self):
        x, L = sp.symbols("x L", positive=True)
        want = {"one": sp.Integer(1), "q1": x - L / 2, "b": x * (L - x),
                "b3": x * (L - x) * (x - L / 2)}
        for name, h in want.items():
            reflected = sp.simplify(sp.expand(h.subs(x, L - x)))
            with self.subTest(name):
                if D3.PARITY[name] == "even":
                    self.assertEqual(sp.simplify(reflected - h), 0)
                else:
                    self.assertEqual(sp.simplify(reflected + h), 0)
                self.assertEqual(pole.basis_parity(name), D3.PARITY[name])

    def test_q1_and_b3_are_the_odd_block(self):
        ids = D3.parity_identities()
        self.assertEqual(set(ids["basis_parity"]), {"one", "q1", "b", "b3"})
        odd = {n for n, p in ids["basis_parity"].items() if p == "odd"}
        self.assertEqual(odd, {"q1", "b3"})


@unittest.skipIf(_flint() is None, "python-flint not installed")
class PrimeShiftIndefiniteness(unittest.TestCase):
    """§7: preserve this. It is a fact about the decomposition, not an obstacle."""

    def test_every_active_prime_shift_block_is_indefinite_on_the_cell(self):
        arb = _flint()
        for Lv in (D3.CELL[0], 1.15, 1.20, 1.30, D3.CELL[1]):
            L = arb(repr(Lv))
            shifts = D3.active_shifts(Lv)
            self.assertTrue(shifts, f"no active shifts at L={Lv}")
            for q, _p in shifts:
                with self.subTest(L=Lv, q=q):
                    det = D3.prime_shift_determinant(q, L, arb)
                    self.assertLess(float(det.upper()), 0.0,
                                    "a prime-shift block must stay indefinite")
                    M = D3.prime_shift_block(q, L, arb)
                    tr = float(M[0][0] + M[1][1])
                    self.assertEqual(inertia_2x2(round(tr, 15),
                                                 round(float(det), 15)), (1, 1, 0))

    def test_termwise_psd_domination_is_not_available(self):
        """A guard against re-deriving a bound by dominating the assembly termwise.

        If any single shift block were PSD this test would pass vacuously, so it
        asserts the opposite: they are all indefinite, which is precisely why the
        assembled entry has to be bounded as a whole.
        """
        arb = _flint()
        L = arb(repr(1.20))
        dets = [float(D3.prime_shift_determinant(q, L, arb))
                for q, _ in D3.active_shifts(1.20)]
        self.assertTrue(dets)
        self.assertTrue(all(d < 0 for d in dets), dets)


if __name__ == "__main__":
    unittest.main(verbosity=1)
