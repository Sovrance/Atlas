#!/usr/bin/env python3
"""ATLAS-RH-ENG-006 §5 — rank-trace theorem runtime checks (WO-RH-30).

The interesting tests here are the refusals. An inequality carried between
settings loses its hypotheses first and its validity second, so most of this
file checks that the engine declines to produce a number when it should:
unverified hypotheses, a normalization that does not match, a positive-index
bound that is violated, a theorem id it does not implement.
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
from ranktrace.theorem import (  # noqa: E402
    HYPOTHESES,
    NORMALIZATION_TAG,
    THEOREM_ID,
    hypotheses_from_matrices,
    rank_trace_lower_bound,
)

ALL_VERIFIED = {name: {"verified": True, "evidence": {"by": "test"}}
                for name, _ in HYPOTHESES}


def hs_sq(M):
    return sum(M[i][j] * M[i][j] for i in range(len(M)) for j in range(len(M)))


def trace(M):
    return sum(M[i][i] for i in range(len(M)))


def add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A))] for i in range(len(A))]


class TheoremEvaluation(unittest.TestCase):
    def test_a_fully_verified_call_produces_a_bound(self):
        cert = rank_trace_lower_bound(
            trace_P=F(10), trace_Q=F(1), hs_sq_P_plus_Q=F(3),
            positive_index_Q_bound=1, hypotheses=ALL_VERIFIED)
        self.assertEqual(cert.status, "PASS")
        # 2*10 + 4*1 - 4*1 - 3 = 17
        self.assertEqual(cert.result["certified_rank_lower_bound"], 17)
        self.assertFalse(cert.result["trivial"])

    def test_the_bound_matches_the_stated_formula_on_exact_inputs(self):
        for tp, tq, hs, b in [(F(10), F(1), F(3), 1), (F(5), F(0), F(2), 0),
                              (F(7), F(2), F(11), 3)]:
            with self.subTest(tp=tp, tq=tq, hs=hs, b=b):
                cert = rank_trace_lower_bound(
                    trace_P=tp, trace_Q=tq, hs_sq_P_plus_Q=hs,
                    positive_index_Q_bound=b, hypotheses=ALL_VERIFIED)
                want = 2 * float(tp) + 4 * float(tq) - 4 * b - float(hs)
                self.assertAlmostEqual(float(cert.result["rhs_enclosure"]["lo"]), want)

    def test_a_non_positive_right_hand_side_is_reported_as_trivial(self):
        """§10: a null result must be preserved, not dressed up as a finding."""
        cert = rank_trace_lower_bound(
            trace_P=F(1), trace_Q=F(0), hs_sq_P_plus_Q=F(50),
            positive_index_Q_bound=2, hypotheses=ALL_VERIFIED)
        self.assertEqual(cert.status, "PASS")
        self.assertTrue(cert.result["trivial"])
        self.assertEqual(cert.result["certified_rank_lower_bound"], 0)
        self.assertIn("says nothing", cert.result["interpretation"])

    def test_the_bound_is_taken_at_the_worst_case_of_an_enclosure(self):
        """Lower bounds must use the enclosure's pessimistic end, not its centre."""
        try:
            from flint import arb, ctx

            ctx.prec = 120
        except ImportError:
            self.skipTest("python-flint not installed")
        cert = rank_trace_lower_bound(
            trace_P=arb("10", "0.5"), trace_Q=arb("1", "0.25"),
            hs_sq_P_plus_Q=arb("3", "0.5"), positive_index_Q_bound=1,
            hypotheses=ALL_VERIFIED)
        lo = float(cert.result["rhs_enclosure"]["lo"])
        # Worst case is 2*9.5 + 4*0.75 - 4 - 3.5 = 14.5, well below the value 17
        # a midpoint evaluation would give.
        self.assertLess(lo, 17.0)
        # Outward rounding must push a *lower* bound down, never up: asserting
        # exact equality with 14.5 would be asserting that the enclosure is not
        # conservative, which is the opposite of what is wanted.
        self.assertLessEqual(lo, 14.5)
        self.assertAlmostEqual(lo, 14.5, places=6)
        hi = float(cert.result["rhs_enclosure"]["hi"])
        self.assertGreaterEqual(hi, 19.5)


class RefusalsAreTheProduct(unittest.TestCase):
    def test_every_single_unverified_hypothesis_blocks_the_result(self):
        for name, _ in HYPOTHESES:
            hyp = {k: dict(v) for k, v in ALL_VERIFIED.items()}
            hyp[name] = {"verified": False, "evidence": {"why": "deliberate"}}
            with self.subTest(name):
                cert = rank_trace_lower_bound(
                    trace_P=F(10), trace_Q=F(1), hs_sq_P_plus_Q=F(3),
                    positive_index_Q_bound=1, hypotheses=hyp)
                self.assertEqual(cert.status, "INCONCLUSIVE")
                self.assertIn(name, cert.blocker)
                self.assertEqual(cert.result, {})

    def test_missing_hypotheses_block_the_result(self):
        cert = rank_trace_lower_bound(
            trace_P=F(10), trace_Q=F(1), hs_sq_P_plus_Q=F(3),
            positive_index_Q_bound=1, hypotheses={})
        self.assertEqual(cert.status, "INCONCLUSIVE")
        for name, _ in HYPOTHESES:
            self.assertIn(name, cert.blocker)

    def test_a_normalization_mismatch_is_refused(self):
        cert = rank_trace_lower_bound(
            trace_P=F(10), trace_Q=F(1), hs_sq_P_plus_Q=F(3),
            positive_index_Q_bound=1, hypotheses=ALL_VERIFIED,
            normalization="some_other_papers_convention")
        self.assertEqual(cert.status, "INCONCLUSIVE")
        self.assertIn("normalization mismatch", cert.blocker)

    def test_an_unimplemented_theorem_id_is_refused(self):
        cert = rank_trace_lower_bound(
            trace_P=F(10), trace_Q=F(1), hs_sq_P_plus_Q=F(3),
            positive_index_Q_bound=1, hypotheses=ALL_VERIFIED,
            theorem_id="some_other_theorem_v9")
        self.assertEqual(cert.status, "INCONCLUSIVE")
        self.assertIn("refusing", cert.blocker)

    def test_a_negative_positive_index_bound_is_refused(self):
        cert = rank_trace_lower_bound(
            trace_P=F(10), trace_Q=F(1), hs_sq_P_plus_Q=F(3),
            positive_index_Q_bound=-1, hypotheses=ALL_VERIFIED)
        self.assertEqual(cert.status, "INCONCLUSIVE")

    def test_no_certificate_ever_claims_rh(self):
        for hyp in (ALL_VERIFIED, {}):
            cert = rank_trace_lower_bound(
                trace_P=F(10), trace_Q=F(1), hs_sq_P_plus_Q=F(3),
                positive_index_Q_bound=1, hypotheses=hyp)
            self.assertIs(cert.to_dict()["rh_proof_claim"], False)
            self.assertEqual(cert.to_dict()["theorem_id"], THEOREM_ID)


class HypothesesCheckedAgainstRealMatrices(unittest.TestCase):
    def test_a_psd_p_and_controlled_signature_q_verify(self):
        P = C.congruence(C.diagonal([1, 1, 0]), C.unimodular(3, 5))
        Q = C.congruence(C.diagonal([1, -1, -1]), C.unimodular(3, 11))
        hyp = hypotheses_from_matrices(P, Q)
        self.assertTrue(hyp["P_positive_semidefinite"]["verified"])
        self.assertTrue(hyp["Q_hermitian"]["verified"])
        self.assertEqual(hyp["_positive_index_Q"], 1)
        cert = rank_trace_lower_bound(
            trace_P=trace(P), trace_Q=trace(Q),
            hs_sq_P_plus_Q=hs_sq(add(P, Q)),
            positive_index_Q_bound=hyp["_positive_index_Q"], hypotheses=hyp)
        self.assertEqual(cert.status, "PASS")

    def test_a_non_psd_p_is_caught_and_blocks_the_bound(self):
        """The hypothesis is checked against the matrix, not taken on trust."""
        P = C.congruence(C.diagonal([1, -1]), C.unimodular(2, 3))
        Q = C.diagonal([1, 0])
        hyp = hypotheses_from_matrices(P, Q)
        self.assertFalse(hyp["P_positive_semidefinite"]["verified"])
        cert = rank_trace_lower_bound(
            trace_P=trace(P), trace_Q=trace(Q), hs_sq_P_plus_Q=hs_sq(add(P, Q)),
            positive_index_Q_bound=1, hypotheses=hyp)
        self.assertEqual(cert.status, "INCONCLUSIVE")
        self.assertIn("P_positive_semidefinite", cert.blocker)

    def test_a_deliberately_violated_positive_index_bound_is_detectable(self):
        """b must bound Q's positive index; the engine reports what Q really has."""
        Q = C.congruence(C.diagonal([1, 1, -1]), C.unimodular(3, 8))
        hyp = hypotheses_from_matrices(C.diagonal([1, 1, 1]), Q)
        real_index = hyp["_positive_index_Q"]
        self.assertEqual(real_index, 2)
        self.assertEqual(exact_inertia(Q).signature, (2, 1, 0))
        # A caller claiming b = 1 is contradicted by the certified inertia.
        self.assertGreater(real_index, 1)

    def test_a_non_symmetric_q_fails_the_hermitian_hypothesis(self):
        Q = [[F(1), F(2)], [F(3), F(1)]]
        hyp = hypotheses_from_matrices(C.diagonal([1, 1]), Q, exact=True)
        self.assertFalse(hyp["Q_hermitian"]["verified"])


if __name__ == "__main__":
    unittest.main(verbosity=1)
