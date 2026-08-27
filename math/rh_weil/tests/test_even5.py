#!/usr/bin/env python3
"""ATLAS-RH-ENG-011 §Required tests — the 5x5 even block, exact half.

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
import even5  # noqa: E402
import reference_metric as RM  # noqa: E402


class BasisIsFrozen(unittest.TestCase):
    def test_the_basis_is_the_one_the_work_order_names(self):
        self.assertEqual(even5.EVEN5_BASIS, ("one", "b", "b2", "bcube", "bquart"))

    def test_every_element_is_even(self):
        import pole
        for name in even5.EVEN5_BASIS:
            self.assertEqual(pole.basis_parity(name), "even", name)

    def test_the_identity_records_what_it_is(self):
        ident = even5.basis_identity()
        self.assertEqual(ident["basis_id"], even5.EVEN5_BASIS_ID)
        self.assertEqual(ident["elements"], ["one", "b", "b2", "bcube", "bquart"])
        self.assertEqual(ident["reference_metric_id"], "l2_gram_on_support")

    def test_the_ten_entries_cover_the_symmetric_block(self):
        pairs = {tuple(sorted(p)) for _, p in even5.ENTRY_KEYS}
        want = {tuple(sorted((i, j)))
                for i in even5.EVEN5_BASIS for j in even5.EVEN5_BASIS}
        self.assertEqual(pairs, want)
        self.assertEqual(len(even5.ENTRY_KEYS), 15)

    def test_bquart_is_the_fourth_power_of_b_exactly(self):
        for L in (F(7, 6), F(5, 4)):
            for x in (F(0), F(1, 3), F(6, 7)):
                b = sum(c * x ** xp * L ** lp
                        for xp, lp_ in enumerate(BA.BASIS_L_POLY["b"])
                        for lp, c in lp_.items())
                bq = sum(c * x ** xp * L ** lp
                         for xp, lp_ in enumerate(BA.BASIS_L_POLY["bquart"])
                         for lp, c in lp_.items())
                self.assertEqual(bq, b ** 4)


class Preconditioner(unittest.TestCase):
    def test_the_exponents_are_frozen_for_the_cell(self):
        self.assertEqual(even5.PRECONDITIONER_EXPONENTS, (-2, -6, -10, -10, -13))

    def test_the_record_states_invertibility_and_the_pencil_licence(self):
        rec = even5.preconditioner_record(even5.PRECONDITIONER_EXPONENTS)
        self.assertTrue(rec["invertible"])
        self.assertTrue(rec["frozen_for_cell"])
        self.assertIn("AtlasRH.generalized_pencil_congruence", rec["licensed_by"])

    def test_applying_it_is_exact_on_dyadics(self):
        # Dyadic entries on purpose: the scaling is by powers of two, and on a
        # binary carrier (Fraction-with-power-of-two-denominator, float, or an
        # Arb ball) it is exact. A denominator of 3 would round through the
        # float product, and nothing in production ever feeds one.
        exps = (-1, 0, 2, 1, -2)
        m = [[F(3), F(1, 2), F(0), F(1), F(1, 4)],
             [F(1, 2), F(2), F(-1), F(0), F(0)],
             [F(0), F(-1), F(5), F(3, 8), F(1)],
             [F(1), F(0), F(3, 8), F(7), F(1, 2)],
             [F(1, 4), F(0), F(1), F(1, 2), F(9)]]
        out = even5.apply_preconditioner(m, exps)
        for a in range(5):
            for b in range(5):
                scale = F(2) ** (-exps[a]) * F(2) ** (-exps[b])
                self.assertEqual(F(out[a][b]), m[a][b] * scale, (a, b))

    def test_minor_scale_factors_are_exact_powers_of_two(self):
        import math
        for s in even5.minor_scale_factors(even5.PRECONDITIONER_EXPONENTS):
            frac, _ = math.frexp(s)
            self.assertEqual(frac, 0.5)

    def test_scaling_never_changes_the_sign_of_a_minor(self):
        m = [[F(3), F(1, 2), F(0), F(1), F(1, 4)],
             [F(1, 2), F(2), F(-1), F(0), F(0)],
             [F(0), F(-1), F(-5), F(1, 3), F(1)],
             [F(1), F(0), F(1, 3), F(7), F(1, 2)],
             [F(1, 4), F(0), F(1), F(1, 2), F(-9)]]
        raw = even5.leading_minors(m)
        pre = even5.leading_minors(even5.apply_preconditioner(m, (-1, 3, -2, 4, 2)))
        for r, p in zip(raw, pre):
            self.assertEqual(F(r) > 0, F(p) > 0)


class Minors(unittest.TestCase):
    def test_leading_minors_agree_with_fraction_determinants(self):
        m = [[F(4), F(1), F(0), F(2), F(1)],
             [F(1), F(3), F(-1), F(0), F(0)],
             [F(0), F(-1), F(2), F(1), F(1, 2)],
             [F(2), F(0), F(1), F(5), F(0)],
             [F(1), F(0), F(1, 2), F(0), F(6)]]
        got = even5.leading_minors(m)
        self.assertEqual(got[0], F(4))
        self.assertEqual(got[1], F(11))
        import generalized_gap as GG
        for k in range(3, 6):
            self.assertEqual(got[k - 1], GG._det([r[:k] for r in m[:k]]))

    def test_a_known_signature_is_respected_by_the_criterion(self):
        # diag(1, 1, 1, 1, -1) fails exactly at the fifth minor
        m = [[F(1) if a == b else F(0) for b in range(5)] for a in range(5)]
        m[4][4] = F(-1)
        minors = even5.leading_minors(m)
        self.assertTrue(all(v > 0 for v in minors[:4]))
        self.assertLess(minors[4], 0)


class ReferenceMetric(unittest.TestCase):
    def test_the_even5_block_metric_is_pd_with_the_expected_minors(self):
        rec = RM.certify_positive_definite(even5.EVEN5_BASIS)
        self.assertEqual(rec["unit_leading_minors"],
                         ["1", "1/180", "1/7938000", "1/88104560544000",
                          "1/248087226834298051200000"])

    def test_shifted_matrix_is_g_minus_lambda_m_exactly(self):
        L = F(6, 5)
        m_ref = [[RM.metric_exact(i, j, L) for j in even5.EVEN5_BASIS]
                 for i in even5.EVEN5_BASIS]
        raw = [[F(a + 1, b + 2) for b in range(5)] for a in range(5)]
        out = even5.shifted_matrix(raw, m_ref, 3, 7)
        for a in range(5):
            for b in range(5):
                self.assertEqual(out[a][b], raw[a][b] - F(3, 7) * m_ref[a][b])


class IndependentAssemblyIsIndependent(unittest.TestCase):
    FORBIDDEN = {"basis_algebra", "pole", "weil_entries",
                 "archimedean_realspace", "even3", "even5", "core",
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
        for rel in ("src/independent_even5.py", "src/independent_even3.py"):
            got = self._imports(rel)
            self.assertFalse(got & self.FORBIDDEN, (rel, got & self.FORBIDDEN))

    def test_it_declares_itself_e3(self):
        import independent_even5 as IE5
        self.assertEqual(IE5.EVIDENCE_CLASS, "E3")


class ScalingModelImmutability(unittest.TestCase):
    """§WO-RH-83: the frozen refit predictions may not move under adjudication."""

    REFIT_CONTENT_HASH = (
        "110d87a6f25eb196a949e71b116429f9f29c9f46ced74938daede6804e97fe60")

    def test_the_frozen_eng010_refit_artifact_is_bitwise_the_committed_one(self):
        import json
        cert = json.loads((ROOT / "certificates" /
                           "e3_eng010_scaling_models_refit.json"
                           ).read_text(encoding="utf-8"))
        self.assertEqual(cert["content_hash"], self.REFIT_CONTENT_HASH,
                         "the frozen ENG-010 refit models changed after "
                         "preregistration; adjudication must run against the "
                         "original artifact (§Stop conditions)")


class NoClaims(unittest.TestCase):
    def test_the_module_states_its_boundary(self):
        self.assertIn("No RH proof claim", even5.__doc__)

    def test_new_content_kinds_have_explicit_psd_answers(self):
        import content_kinds as CK
        self.assertTrue(CK.psd_licensable(CK.KIND_DEGREE8_POSITIVITY))
        self.assertFalse(CK.psd_licensable(CK.KIND_SCALING_ADJUDICATION))


if __name__ == "__main__":
    unittest.main(verbosity=2)
