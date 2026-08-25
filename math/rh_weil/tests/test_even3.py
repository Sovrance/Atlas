#!/usr/bin/env python3
"""ATLAS-RH-ENG-008 — the 3x3 even block, its preconditioner, and its limits.

The tests that matter here are the ones that would catch a *wrong* certificate
rather than a broken one. Three groups do most of that work:

* the preconditioner is a congruence and nothing else -- exactly invertible,
  exactly applied, and inertia-preserving on matrices whose signature is known
  independently;
* the two routes to the signature, interval LDL and Sylvester's minors, agree on
  matrices where the answer is known, and both fail closed on a singular one;
* the independent assembly really is independent -- it imports none of the code
  it checks, asserted by reading its imports.

Anything that needs python-flint skips without it, so the exact half of this
file runs anywhere.
"""
from __future__ import annotations

import ast
import json
import math
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import basis_algebra as BA  # noqa: E402,F401
import even3  # noqa: E402
from inertia.congruence import charpoly_inertia  # noqa: E402
from inertia.ldl import exact_inertia, interval_inertia  # noqa: E402

CERT_DIR = ROOT / "certificates"


def _flint():
    try:
        from interval_backend import require_flint

        return require_flint()
    except Exception:  # pragma: no cover
        raise unittest.SkipTest("python-flint is required")


def _load(name):
    p = CERT_DIR / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


class BasisIsFrozen(unittest.TestCase):
    def test_the_basis_is_the_one_the_work_order_names(self):
        self.assertEqual(even3.EVEN3_BASIS, ("one", "b", "b2"))

    def test_every_element_is_even(self):
        import pole

        for name in even3.EVEN3_BASIS:
            self.assertEqual(pole.basis_parity(name), "even")

    def test_the_basis_identity_records_what_it_is(self):
        ident = even3.basis_identity()
        self.assertEqual(ident["basis_id"], even3.EVEN3_BASIS_ID)
        self.assertEqual(ident["elements"], ["one", "b", "b2"])
        self.assertIn("u^4", ident["spans_in_u"])

    def test_the_six_entries_cover_the_symmetric_block(self):
        pairs = {tuple(sorted(p)) for _, p in even3.ENTRY_KEYS}
        expected = {tuple(sorted((i, j)))
                    for i in even3.EVEN3_BASIS for j in even3.EVEN3_BASIS}
        self.assertEqual(pairs, expected)
        self.assertEqual(len(even3.ENTRY_KEYS), 6)


class Preconditioner(unittest.TestCase):
    """It has to be a congruence, and it has to be exactly one."""

    def test_the_exponents_are_frozen_for_the_cell(self):
        self.assertEqual(len(even3.PRECONDITIONER_EXPONENTS), 3)
        for e in even3.PRECONDITIONER_EXPONENTS:
            self.assertIsInstance(e, int)

    def test_the_record_states_invertibility_and_what_licenses_it(self):
        rec = even3.preconditioner_record(even3.PRECONDITIONER_EXPONENTS)
        self.assertIs(rec["invertible"], True)
        self.assertIs(rec["exactly_representable"], True)
        self.assertIn("AtlasRH.posIndexAtLeast_congruence_iff", rec["licensed_by"])
        self.assertIn("rank_congruence", " ".join(rec["licensed_by"]))

    def test_applying_it_is_exact_on_rationals(self):
        # Powers of two, so no rounding: the identity holds in exact arithmetic.
        A = [[F(3), F(1), F(-2)], [F(1), F(5), F(7)], [F(-2), F(7), F(11)]]
        exps = [-2, -6, -9]
        got = even3.apply_preconditioner(A, exps)
        for i in range(3):
            for j in range(3):
                want = A[i][j] * F(2) ** (-exps[i]) * F(2) ** (-exps[j])
                self.assertEqual(F(got[i][j]), want, (i, j))

    def test_it_preserves_inertia_on_matrices_with_a_known_signature(self):
        cases = [
            ([[F(2), F(0), F(0)], [F(0), F(3), F(0)], [F(0), F(0), F(5)]], (3, 0, 0)),
            ([[F(-2), F(0), F(0)], [F(0), F(3), F(0)], [F(0), F(0), F(5)]], (2, 1, 0)),
            ([[F(1), F(2), F(0)], [F(2), F(1), F(0)], [F(0), F(0), F(-7)]], (1, 2, 0)),
            ([[F(-1), F(0), F(0)], [F(0), F(-3), F(0)], [F(0), F(0), F(-5)]], (0, 3, 0)),
        ]
        for A, want in cases:
            self.assertEqual(exact_inertia(A).signature, want, A)
            for exps in ([-2, -6, -9], [3, -1, 4], [0, 0, 0]):
                B = even3.apply_preconditioner(A, exps)
                self.assertEqual(exact_inertia(B).signature, want, (A, exps))
                # And the independent oracle agrees, so this is not a property
                # of the elimination order.
                self.assertEqual(charpoly_inertia(B), want, (A, exps))

    def test_the_minor_scale_factors_are_exact_powers_of_two(self):
        scales = even3.minor_scale_factors([-2, -6, -9])
        self.assertEqual(scales, [16.0, 65536.0, 17179869184.0])
        for s in scales:
            self.assertEqual(s, 2 ** round(math.log2(s)))

    def test_scaling_never_changes_the_sign_of_a_minor(self):
        A = [[F(3), F(1), F(-2)], [F(1), F(5), F(7)], [F(-2), F(7), F(11)]]
        for exps in ([-2, -6, -9], [5, 0, -3]):
            B = even3.apply_preconditioner(A, exps)
            for a, b in zip(even3.leading_minors(A), even3.leading_minors(B)):
                self.assertEqual((F(a) > 0), (F(b) > 0))


class SignatureRoutesAgree(unittest.TestCase):
    """LDL and the leading minors, on matrices whose answer is known."""

    def test_sylvester_and_ldl_agree_on_definite_matrices(self):
        cases = [
            [[F(2), F(0), F(0)], [F(0), F(3), F(0)], [F(0), F(0), F(5)]],
            [[F(4), F(1), F(1)], [F(1), F(4), F(1)], [F(1), F(1), F(4)]],
            [[F(1), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(1)]],
        ]
        for A in cases:
            minors = [F(m) for m in even3.leading_minors(A)]
            self.assertTrue(all(m > 0 for m in minors), A)
            self.assertEqual(exact_inertia(A).signature, (3, 0, 0), A)

    def test_a_failing_minor_means_not_definite(self):
        cases = [
            [[F(1), F(2), F(0)], [F(2), F(1), F(0)], [F(0), F(0), F(1)]],
            [[F(-1), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(1)]],
        ]
        for A in cases:
            minors = [F(m) for m in even3.leading_minors(A)]
            self.assertFalse(all(m > 0 for m in minors), A)
            self.assertNotEqual(exact_inertia(A).signature, (3, 0, 0), A)

    def test_a_singular_matrix_is_inconclusive_on_the_interval_path(self):
        # §Stop conditions: no exact zero may be inferred from numerical
        # smallness. A ball straddling zero must not resolve.
        _, arb, _, ctx = _flint()
        ctx.prec = 128
        eps = arb(0, 1e-30)  # a ball centred on zero
        A = [[arb(1), arb(0), arb(0)],
             [arb(0), arb(1), arb(0)],
             [arb(0), arb(0), eps]]
        res = interval_inertia(A)
        self.assertEqual(res.status, "INCONCLUSIVE")
        self.assertIsNone(res.signature)

    def test_a_definite_interval_matrix_resolves(self):
        _, arb, _, ctx = _flint()
        ctx.prec = 128
        A = [[arb(4), arb(1), arb(1)],
             [arb(1), arb(4), arb(1)],
             [arb(1), arb(1), arb(4)]]
        res = interval_inertia(A)
        self.assertEqual(res.status, "PASS")
        self.assertEqual(res.signature, (3, 0, 0))


class TheIndependentAssemblyIsIndependent(unittest.TestCase):
    """A cross-check that shares code is a re-run, not a check."""

    FORBIDDEN = ("basis_algebra", "pole", "weil_entries",
                 "archimedean_realspace", "even3", "core", "interval_backend")

    def test_it_imports_none_of_the_code_it_checks(self):
        src = (ROOT / "src" / "independent_even3.py").read_text(encoding="utf-8")
        tree = ast.parse(src)
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])
        for name in self.FORBIDDEN:
            self.assertNotIn(name, imported,
                             f"the independent assembly imports {name}")

    def test_it_declares_itself_e3(self):
        import independent_even3 as IE

        self.assertEqual(IE.EVIDENCE_CLASS, "E3")

    def test_it_agrees_with_the_rigorous_assembly_at_the_cell_midpoint(self):
        _flint()
        import independent_even3 as IE

        try:
            mp = IE.require_mpmath()
            IE.require_sympy()
        except ImportError:  # pragma: no cover
            raise unittest.SkipTest("mpmath and sympy are required")
        L = (even3.CELL[0] + even3.CELL[1]) / 2
        rigorous = even3.assemble_even3_arb(L, precision_bits=160)
        M = IE.gram_matrix(L, mp, dps=40)
        indep = {"G00": M[0][0], "G01": M[0][1], "G02": M[0][2],
                 "G11": M[1][1], "G12": M[1][2], "G22": M[2][2]}
        for key, _ in even3.ENTRY_KEYS:
            a = rigorous["entries"][key]
            v = float(indep[key])
            self.assertLessEqual(float(a.lower()), v, key)
            self.assertGreaterEqual(float(a.upper()), v, key)


class CertificateSemantics(unittest.TestCase):
    def setUp(self):
        self.inertia = _load("e1_degree4_even3_inertia_log3_log4.json")
        self.positivity = _load("e1_degree4_even3_positivity_log3_log4.json")
        self.crosscheck = _load("e3_degree4_even3_crosscheck.json")
        if self.inertia is None:
            self.skipTest("run scripts/certify_even3.py first")

    def test_the_inertia_certificate_never_claims_psd(self):
        # ENG-006 §11 is categorical: an inertia artifact refuses a PSD consumer
        # whatever its signature. The positivity artifact is what claims it.
        from inertia.certificate import satisfies_psd_requirement

        self.assertIs(self.inertia["psd_claim"], False)
        self.assertFalse(satisfies_psd_requirement(self.inertia))

    def test_the_positivity_certificate_does_claim_psd_and_is_licensed(self):
        from inertia.certificate import satisfies_psd_requirement

        if self.positivity is None:
            self.skipTest("positivity was not certified")
        self.assertIs(self.positivity["psd_claim"], True)
        self.assertTrue(satisfies_psd_requirement(self.positivity))

    def test_the_crosscheck_is_e3_and_can_never_promote(self):
        from inertia.certificate import satisfies_psd_requirement

        if self.crosscheck is None:
            self.skipTest("cross-check not generated")
        self.assertEqual(self.crosscheck["evidence_class"], "E3")
        self.assertIs(self.crosscheck["rigorous"], False)
        self.assertIs(self.crosscheck["mpmath_used"], True)
        self.assertFalse(satisfies_psd_requirement(self.crosscheck))

    def test_the_two_routes_agree(self):
        if self.positivity is None:
            self.skipTest("positivity was not certified")
        self.assertEqual(self.inertia["signatures_seen"], [[3, 0, 0]])
        self.assertEqual(
            [self.positivity["n_positive"], self.positivity["n_negative"],
             self.positivity["n_zero"]], [3, 0, 0])
        for cover in self.positivity["leading_minors"]:
            self.assertGreater(float(cover["certified_lower_bound"]), 0.0,
                               cover["minor"])

    def test_no_certificate_claims_an_rh_proof(self):
        for cert in (self.inertia, self.positivity, self.crosscheck):
            if cert is not None:
                self.assertIs(cert["rh_proof_claim"], False)

    def test_no_eigenvalue_solver_on_the_rigorous_path(self):
        self.assertIs(self.inertia["mpmath_used"], False)
        self.assertIn("no eigenvalue solver", self.inertia["method"])

    def test_delta2_reconciles_with_the_certified_two_by_two_block(self):
        # §WO-RH-51: the Delta2 route must reconcile with the previously
        # certified even 2x2, whose E2 is literally this block's second leading
        # minor. Both are valid lower bounds on the same quantity over the same
        # cell, so they must not contradict -- and they should be within an
        # order of magnitude, since they bound the same minimum.
        if self.positivity is None:
            self.skipTest("positivity was not certified")
        old = _load("e1_degree2_compact_log3_log4.json")
        if old is None:
            self.skipTest("degree-2 certificate absent")
        d2 = next(c for c in self.positivity["leading_minors"]
                  if c["minor"] == "Delta2")
        mine = float(d2["implied_raw_lower_bound"])
        theirs = float(old["certified_lower_bound"])
        self.assertGreater(mine, 0.0)
        self.assertGreater(theirs, 0.0)
        self.assertLess(abs(math.log10(mine) - math.log10(theirs)), 1.0,
                        f"Delta2 bound {mine} vs certified E2 bound {theirs}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
