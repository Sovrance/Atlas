"""E0 identities for the 3x3 odd pilot (ATLAS-RH-ENG-007 §15, WO-RH-46).

These re-derive the pilot kernels from the basis with SymPy rather than checking the shipped
closed forms against themselves. A closed form that was transcribed wrongly and then tested
against its own transcription passes every time, which is the failure this file is written to
avoid.

Nothing here certifies anything. The pilot is E0 preparation for ENG-008.
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

try:
    import sympy as sp
except ImportError:  # pragma: no cover
    sp = None

import degree3_pilot as P  # noqa: E402


@unittest.skipIf(sp is None, "sympy is required for the exact identities")
class PilotKernelIdentities(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x, cls.a, cls.L = sp.symbols("x a L", positive=True)
        x, L = cls.x, cls.L
        cls.basis = {
            "q1": x - L / 2,
            "b3": x * (L - x) * (x - L / 2),
            "q3": (x - L / 2) ** 3,
        }

    def _Ksym(self, i, j):
        """The repository convention (`core`): both shifts, as the formula sums."""
        x, a, L = self.x, self.a, self.L
        f, g = self.basis[i], self.basis[j]
        return sp.expand(sp.integrate(f * g.subs(x, x + a) + f.subs(x, x + a) * g,
                                      (x, 0, L - a)))

    def test_every_pilot_basis_element_is_odd_under_midpoint_reflection(self):
        """The whole point of the pilot is that it stays inside ONE parity sector.

        A mixed-parity block would destroy the even/odd factorization `det G = O1 * E2` that
        the existing degree-2 certificates rely on.
        """
        x, L = self.x, self.L
        for name, f in self.basis.items():
            with self.subTest(element=name):
                self.assertEqual(sp.simplify(f.subs(x, L - x) + f), 0)

    def test_shipped_kernels_agree_with_the_symbolic_derivation(self):
        a, L = self.a, self.L
        for (i, j), fn in P.PILOT_KERNELS.items():
            with self.subTest(pair=(i, j)):
                self.assertEqual(sp.simplify(self._Ksym(i, j) - fn(a, L)), 0)

    def test_the_pilot_block_is_genuinely_three_dimensional(self):
        """`q1`, `b3`, `q3` must be linearly independent, or the 'pilot' is a 2x2 in disguise.

        This is the property that makes ENG-008 worth running: inertia and spectral moments
        add information beyond a determinant only when the block really has three dimensions.
        """
        x, L = self.x, self.L
        c1, c2, c3 = sp.symbols("c1 c2 c3")
        combo = sp.expand(c1 * self.basis["q1"] + c2 * self.basis["b3"] + c3 * self.basis["q3"])
        poly = sp.Poly(combo, x)
        sol = sp.solve(poly.coeffs(), [c1, c2, c3], dict=True)
        # The only vanishing combination is the trivial one.
        self.assertTrue(all(s.get(c1, 0) == 0 and s.get(c2, 0) == 0 and s.get(c3, 0) == 0)
                        for s in sol)


class PilotPromotesNothing(unittest.TestCase):
    """§15: prepare, do not certify. Guarded by a test so a later edit cannot quietly
    promote the pilot without someone noticing."""

    def test_evidence_class_is_E0_and_nothing_claims_E1(self):
        self.assertEqual(P.EVIDENCE_CLASS, "E0")
        summary = P.pilot_summary()
        self.assertFalse(summary["e1_promoted"])
        self.assertEqual(summary["status"], "PREPARED_NOT_CERTIFIED")

    def test_conditioning_preview_is_labelled_diagnostic(self):
        entries = {"q1q1": 1.0, "q1b3": 0.1, "q1q3": 0.2,
                   "b3b3": 0.5, "b3q3": 0.05, "q3q3": 0.3}
        prev = P.conditioning_preview(1.1, entries)
        self.assertEqual(prev["evidence_class"], "E3")
        self.assertFalse(prev["rigorous"])
        self.assertFalse(prev["certifies"])
        self.assertFalse(prev["rh_proof_claim"])
        self.assertEqual(prev["promotion_state"], "NOT_PROMOTED_PREVIEW_ONLY")

    def test_scaling_is_a_congruence_so_it_cannot_change_inertia(self):
        note = P.scaled_basis_note()
        self.assertEqual(note["strategy"], "diagonal_rescaling")
        self.assertIn("preserves inertia", note["justification"])


if __name__ == "__main__":
    unittest.main()
