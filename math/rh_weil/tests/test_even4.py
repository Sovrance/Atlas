#!/usr/bin/env python3
"""ATLAS-RH-ENG-010 §15 — the 4x4 even block, exact half.

Everything here runs without python-flint unless a test says otherwise.
"""
from __future__ import annotations

import ast
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for extra in (ROOT, ROOT / "src"):
    sys.path.insert(0, str(extra))

import basis_algebra as BA  # noqa: E402
import even4  # noqa: E402
import reference_metric as RM  # noqa: E402


class BasisIsFrozen(unittest.TestCase):
    def test_the_basis_is_the_one_the_work_order_names(self):
        self.assertEqual(even4.EVEN4_BASIS, ("one", "b", "b2", "bcube"))

    def test_every_element_is_even(self):
        import pole
        for name in even4.EVEN4_BASIS:
            self.assertEqual(pole.basis_parity(name), "even", name)

    def test_the_identity_records_what_it_is(self):
        ident = even4.basis_identity()
        self.assertEqual(ident["basis_id"], even4.EVEN4_BASIS_ID)
        self.assertEqual(ident["elements"], ["one", "b", "b2", "bcube"])
        self.assertEqual(ident["reference_metric_id"], "l2_gram_on_support")

    def test_the_ten_entries_cover_the_symmetric_block(self):
        pairs = {tuple(sorted(p)) for _, p in even4.ENTRY_KEYS}
        want = {tuple(sorted((i, j)))
                for i in even4.EVEN4_BASIS for j in even4.EVEN4_BASIS}
        self.assertEqual(pairs, want)
        self.assertEqual(len(even4.ENTRY_KEYS), 10)

    def test_bcube_is_the_cube_of_b_exactly(self):
        for L in (F(7, 6), F(5, 4)):
            for x in (F(0), F(1, 3), F(6, 7)):
                b = sum(c * x ** xp * L ** lp
                        for xp, lp_ in enumerate(BA.BASIS_L_POLY["b"])
                        for lp, c in lp_.items())
                bc = sum(c * x ** xp * L ** lp
                         for xp, lp_ in enumerate(BA.BASIS_L_POLY["bcube"])
                         for lp, c in lp_.items())
                self.assertEqual(bc, b ** 3)


class Preconditioner(unittest.TestCase):
    def test_the_exponents_are_frozen_for_the_cell(self):
        self.assertEqual(even4.PRECONDITIONER_EXPONENTS, (-2, -6, -10, -10))

    def test_the_record_states_invertibility_and_the_pencil_licence(self):
        rec = even4.preconditioner_record(even4.PRECONDITIONER_EXPONENTS)
        self.assertTrue(rec["invertible"])
        self.assertTrue(rec["frozen_for_cell"])
        self.assertIn("AtlasRH.generalized_pencil_congruence", rec["licensed_by"])

    def test_applying_it_is_exact_on_dyadics(self):
        # Dyadic entries on purpose: the scaling is by powers of two, and on a
        # binary carrier (Fraction-with-power-of-two-denominator, float, or an
        # Arb ball) it is exact. A denominator of 3 would round through the
        # float product, and nothing in production ever feeds one.
        exps = (-1, 0, 2, 1)
        m = [[F(3), F(1, 2), F(0), F(1)],
             [F(1, 2), F(2), F(-1), F(0)],
             [F(0), F(-1), F(5), F(3, 8)],
             [F(1), F(0), F(3, 8), F(7)]]
        out = even4.apply_preconditioner(m, exps)
        for a in range(4):
            for b in range(4):
                scale = F(2) ** (-exps[a]) * F(2) ** (-exps[b])
                self.assertEqual(F(out[a][b]), m[a][b] * scale, (a, b))

    def test_minor_scale_factors_are_exact_powers_of_two(self):
        import math
        for s in even4.minor_scale_factors(even4.PRECONDITIONER_EXPONENTS):
            frac, _ = math.frexp(s)
            self.assertEqual(frac, 0.5)

    def test_scaling_never_changes_the_sign_of_a_minor(self):
        m = [[F(3), F(1, 2), F(0), F(1)],
             [F(1, 2), F(2), F(-1), F(0)],
             [F(0), F(-1), F(-5), F(1, 3)],
             [F(1), F(0), F(1, 3), F(7)]]
        raw = even4.leading_minors(m)
        pre = even4.leading_minors(even4.apply_preconditioner(m, (-1, 3, -2, 4)))
        for r, p in zip(raw, pre):
            self.assertEqual(F(r) > 0, F(p) > 0)


class Minors(unittest.TestCase):
    def test_leading_minors_agree_with_fraction_determinants(self):
        m = [[F(4), F(1), F(0), F(2)],
             [F(1), F(3), F(-1), F(0)],
             [F(0), F(-1), F(2), F(1)],
             [F(2), F(0), F(1), F(5)]]
        got = even4.leading_minors(m)
        self.assertEqual(got[0], F(4))
        self.assertEqual(got[1], F(11))
        # cross-check Delta3/Delta4 against an independent expansion
        import generalized_gap as GG
        self.assertEqual(got[2], GG._det([r[:3] for r in m[:3]]))
        self.assertEqual(got[3], GG._det(m))

    def test_a_known_signature_is_respected_by_the_criterion(self):
        # diag(1, 1, 1, -1) fails exactly at the fourth minor
        m = [[F(1), F(0), F(0), F(0)],
             [F(0), F(1), F(0), F(0)],
             [F(0), F(0), F(1), F(0)],
             [F(0), F(0), F(0), F(-1)]]
        minors = even4.leading_minors(m)
        self.assertTrue(all(v > 0 for v in minors[:3]))
        self.assertLess(minors[3], 0)


class ReferenceMetric(unittest.TestCase):
    def test_the_even4_block_metric_is_pd_with_the_expected_minors(self):
        rec = RM.certify_positive_definite(even4.EVEN4_BASIS)
        self.assertEqual(rec["unit_leading_minors"],
                         ["1", "1/180", "1/7938000", "1/88104560544000"])

    def test_shifted_matrix_is_g_minus_lambda_m_exactly(self):
        L = F(6, 5)
        m_ref = [[RM.metric_exact(i, j, L) for j in even4.EVEN4_BASIS]
                 for i in even4.EVEN4_BASIS]
        raw = [[F(a + 1, b + 2) for b in range(4)] for a in range(4)]
        out = even4.shifted_matrix(raw, m_ref, 3, 7)
        for a in range(4):
            for b in range(4):
                self.assertEqual(out[a][b], raw[a][b] - F(3, 7) * m_ref[a][b])


class IndependentAssemblyIsIndependent(unittest.TestCase):
    FORBIDDEN = {"basis_algebra", "pole", "weil_entries",
                 "archimedean_realspace", "even3", "even4", "core",
                 "interval_backend", "generalized_gap", "reference_metric"}

    def _imports(self, relpath):
        tree = ast.parse((ROOT / relpath).read_text(encoding="utf-8"))
        out = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                out |= {a.name.split(".")[0] for a in node.names}
            elif isinstance(node, ast.ImportFrom) and node.module:
                out.add(node.module.split(".")[0])
        return out

    def test_it_imports_none_of_the_code_it_checks(self):
        for rel in ("src/independent_even4.py", "src/independent_even3.py"):
            got = self._imports(rel)
            self.assertFalse(got & self.FORBIDDEN, (rel, got & self.FORBIDDEN))

    def test_it_declares_itself_e3(self):
        import independent_even4 as IE4
        self.assertEqual(IE4.EVIDENCE_CLASS, "E3")


class ScalingModelImmutability(unittest.TestCase):
    """§WO-RH-71: the preregistered predictions may not move under adjudication."""

    BASELINE_CONTENT_HASH = (
        "2719ade96da77279ea6350fdfb19f49c1a22434e970fafc302036480a85d6a23")

    def test_the_eng009_models_artifact_is_bitwise_the_committed_one(self):
        import json
        cert = json.loads((ROOT / "certificates" / "e3_eng009_scaling_models.json"
                           ).read_text(encoding="utf-8"))
        self.assertEqual(cert["content_hash"], self.BASELINE_CONTENT_HASH,
                         "the ENG-009 scaling models changed after "
                         "preregistration; adjudication must run against the "
                         "original artifact (§Stop conditions)")


class NoClaims(unittest.TestCase):
    def test_the_module_states_its_boundary(self):
        self.assertIn("No RH proof claim", even4.__doc__)

    def test_new_content_kinds_have_explicit_psd_answers(self):
        import content_kinds as CK
        self.assertTrue(CK.psd_licensable(CK.KIND_DEGREE6_POSITIVITY))
        self.assertFalse(CK.psd_licensable(CK.KIND_SCALING_ADJUDICATION))


if __name__ == "__main__":
    unittest.main(verbosity=2)
