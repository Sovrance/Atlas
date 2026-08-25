#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §15 (WO-RH-46) — exact identities for the 3x3 pilot block.

Everything here is E0: exact rational arithmetic, no floating point, no
quadrature. The three new prime-overlap kernels are checked against exact
integration of their own integrand, so the closed forms are verified rather
than asserted, and against their coefficient expansions, so the two
representations the assembly routes use cannot drift apart.

The last class checks the thing that would quietly break a certified result:
that this module has not modified the basis, the kernels or the pole primitive
the promoted E1 certificates depend on.
"""
from __future__ import annotations

import math
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import core  # noqa: E402
import pilot3 as P  # noqa: E402

#: Rational test points. `a` must stay below `L`; the kernels are integrals over
#: `[0, L - a]` and are not claimed outside that.
LS = (F(7, 3), F(11, 8), F(5, 2), F(3, 1))
AS = (F(0), F(1, 4), F(5, 7), F(9, 8), F(2))
XS = (F(1, 5), F(2, 7), F(11, 13), F(3, 2), F(17, 4))


class Parity(unittest.TestCase):
    def test_the_even_sector_is_even_about_the_midpoint(self):
        for L in LS:
            for name in P.EVEN_BASIS:
                self.assertTrue(P.is_even_about_midpoint(name, L, XS), (name, L))

    def test_the_odd_sector_is_odd_about_the_midpoint(self):
        for L in LS:
            for name in P.ODD_BASIS:
                self.assertTrue(P.is_odd_about_midpoint(name, L, XS), (name, L))

    def test_b2_is_the_square_of_b(self):
        for L in LS:
            for x in XS:
                self.assertEqual(P.evaluate("b2", x, L), P.evaluate("b", x, L) ** 2)

    def test_the_even_sector_is_three_dimensional(self):
        # In u = x - L/2 the sector is span{1, u^2, u^4}. Independence is the
        # whole reason this block is worth building: a dependent third element
        # would give a singular Gram and nothing to learn.
        for L in LS:
            rows = []
            for x in XS[:3]:
                rows.append([P.evaluate(n, x, L) for n in P.EVEN_BASIS])
            det = (rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
                   - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
                   + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0]))
            self.assertNotEqual(det, 0, L)

    def test_q1_cubed_adds_nothing_to_the_odd_sector(self):
        # (x - L/2)^3 = (L^2/4) q1 - b3, so extending the odd sector at degree 3
        # is impossible; this is why the pilot extends the even sector instead.
        for L in LS:
            for x in XS:
                u = x - L / 2
                self.assertEqual(
                    u ** 3,
                    (L * L / 4) * P.evaluate("q1", x, L) - P.evaluate("b3", x, L),
                )


class NewKernels(unittest.TestCase):
    """The closed forms, against exact integration of the integrand."""

    PAIRS = (("one", "b2"), ("b", "b2"), ("b2", "b2"))

    def test_closed_form_matches_exact_integration(self):
        for L in LS:
            for a in AS:
                if a >= L:
                    continue
                for i, j in self.PAIRS:
                    self.assertEqual(
                        F(P.kernel(i, j, a, L)),
                        P.kernel_by_quadrature(i, j, a, L),
                        (i, j, a, L),
                    )

    def test_coefficient_expansion_matches_the_closed_form(self):
        for L in LS:
            for a in AS:
                if a >= L:
                    continue
                for i, j in self.PAIRS:
                    self.assertEqual(
                        P.kernel_from_coeffs(i, j, a, L),
                        F(P.kernel(i, j, a, L)),
                        (i, j, a, L),
                    )

    def test_the_kernels_are_symmetric_in_their_arguments(self):
        for L in LS:
            for a in AS:
                if a >= L:
                    continue
                for i, j in self.PAIRS:
                    self.assertEqual(P.kernel(i, j, a, L), P.kernel(j, i, a, L))

    def test_every_kernel_vanishes_at_the_right_endpoint(self):
        # K_ij(L; L) integrates over the empty interval [0, 0].
        for L in LS:
            for i, j in self.PAIRS:
                self.assertEqual(P.kernel(i, j, L, L), 0, (i, j, L))

    def test_the_diagonal_kernel_is_positive_inside_the_cell(self):
        # K_b2b2(a; L) = 2 int_0^{L-a} b2(x) b2(x+a) dx with b2 >= 0.
        for L in LS:
            for a in AS:
                if a >= L:
                    continue
                self.assertGreater(P.kernel("b2", "b2", a, L), 0, (a, L))


class LegacyKernelsUnchanged(unittest.TestCase):
    """The pilot reuses the ENG-005/006 kernels; it does not redefine them."""

    LEGACY = (("one", "one", core.kernel_00), ("one", "b", core.kernel_0b),
              ("b", "b", core.kernel_bb), ("q1", "q1", core.kernel_q1q1),
              ("q1", "b3", core.kernel_q1b3), ("b3", "b3", core.kernel_b3b3))

    def test_pilot_kernel_dispatch_agrees_with_core(self):
        for L in LS:
            for a in AS:
                if a >= L:
                    continue
                for i, j, fn in self.LEGACY:
                    self.assertEqual(P.kernel(i, j, a, L), fn(a, L), (i, j, a, L))

    def test_legacy_kernels_also_pass_the_independent_route(self):
        # Not strictly the pilot's business, but it costs nothing and it means a
        # failure here localizes to core.py rather than to the new code.
        for L in LS:
            for a in AS:
                if a >= L:
                    continue
                for i, j, fn in self.LEGACY:
                    self.assertEqual(
                        F(fn(a, L)), P.kernel_by_quadrature(i, j, a, L), (i, j, a, L)
                    )

    def test_the_pilot_basis_agrees_with_the_production_basis(self):
        import pole

        for L in LS:
            for name in ("one", "b", "q1", "b3"):
                self.assertEqual(
                    tuple(F(c) for c in P.basis_coeffs(name, L)),
                    tuple(F(c) for c in pole.basis_coeffs(name, L)),
                    name,
                )

    def test_production_now_carries_b2_as_a_frozen_basis_element(self):
        # ENG-007 kept b2 out of production deliberately: extending
        # pole.basis_coeffs would have changed pole.py's source hash and staled
        # every certificate that binds it, which a preparatory pilot must not do.
        # ENG-008 §WO-RH-47 freezes b2 into the canonical basis on purpose, and
        # regenerated the chain to pay for it. This test records that reversal
        # rather than leaving the old prohibition standing.
        import pole

        self.assertIn("b2", pole.BASIS_NAMES)
        self.assertEqual(pole.basis_parity("b2"), "even")
        for L in LS:
            self.assertEqual(
                tuple(F(c) for c in P.basis_coeffs("b2", L)),
                tuple(F(c) for c in pole.basis_coeffs("b2", L)),
                L)

    def test_the_pilot_kernels_still_agree_with_production(self):
        # The pilot keeps its own kernel table, and that is still useful: it is
        # a second derivation of the same polynomials. What changed is only that
        # production now has them too.
        import weil_entries as WE

        for i, j in (("one", "b2"), ("b", "b2"), ("b2", "b2")):
            for L in LS:
                for a in AS:
                    if a >= L:
                        continue
                    self.assertEqual(F(P.kernel(i, j, a, L)),
                                     F(WE.kernel(i, j, a, L)), (i, j, a, L))


class ParityBlockStructure(unittest.TestCase):
    """Cross-parity entries vanish, exactly, in both the prime and pole blocks."""

    def test_cross_parity_prime_kernels_are_not_defined(self):
        # There is no even/odd kernel because there is no even/odd entry: the
        # parity argument (AtlasRH.cross_block_vanishes) says it is zero.
        for i in P.EVEN_BASIS:
            for j in P.ODD_BASIS:
                with self.assertRaises(KeyError):
                    P.kernel(i, j, F(1, 2), F(5, 2))

    def test_the_parity_map_covers_the_whole_basis(self):
        self.assertEqual(
            set(P.PARITY), set(P.EVEN_BASIS) | set(P.ODD_BASIS)
        )
        for n in P.EVEN_BASIS:
            self.assertEqual(P.PARITY[n], "even")
        for n in P.ODD_BASIS:
            self.assertEqual(P.PARITY[n], "odd")


class NoClaims(unittest.TestCase):
    def test_the_module_states_its_evidence_boundary(self):
        self.assertIn("E0", P.EVIDENCE_NOTE)
        self.assertIn("E3", P.EVIDENCE_NOTE)
        self.assertIn("no E1", P.EVIDENCE_NOTE)

    def test_the_preview_certificate_is_e3_and_claims_nothing(self):
        import json

        path = ROOT / "certificates" / "e3_pilot3_even_conditioning_log3_log4.json"
        if not path.exists():
            self.skipTest("preview not generated")
        cert = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(cert["evidence_class"], "E3")
        self.assertIs(cert["rh_proof_claim"], False)
        self.assertIs(cert["rigorous"], False)
        self.assertIs(cert["psd_claim"], False)
        self.assertIs(cert["hard_constraints_certified"], False)
        self.assertIs(cert["mpmath_used"], True)
        self.assertEqual(cert["status"], "PREVIEW")


if __name__ == "__main__":
    unittest.main(verbosity=2)
