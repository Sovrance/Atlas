#!/usr/bin/env python3
"""ATLAS-RH-ENG-006 §3/§4 — inertia engine and exact congruence regression.

The engine is checked against an oracle that shares none of its code: the LDL
elimination is compared with an inertia read off the characteristic polynomial
via Descartes' rule of signs, which is exact for the real-rooted polynomial of a
symmetric matrix. Agreement between two routes that could fail independently is
the point; a self-consistent engine proves nothing.

The singular and near-singular cases are the ones that matter most. An interval
enclosure cannot witness an exact zero, so those must come back INCONCLUSIVE --
a test that let them return a signature would be certifying a guess.
"""
from __future__ import annotations

import json
import os
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import inertia.congruence as C  # noqa: E402
from inertia.certificate import (  # noqa: E402
    KIND_INERTIA,
    KIND_STRATIFICATION,
    build_inertia_certificate,
    build_stratification_certificate,
    satisfies_psd_requirement,
    validate_against_schema,
)
from inertia.ldl import exact_inertia, interval_inertia  # noqa: E402
from inertia.stratify import certify_inertia_family  # noqa: E402

SCHEMA_DIR = ROOT / "inertia" / "schemas"


def _flint():
    try:
        from flint import arb, ctx

        ctx.prec = 120
        return arb
    except ImportError:  # pragma: no cover
        return None


def _balls(A, rad=0):
    arb = _flint()
    return [[arb(str(x), str(rad)) if rad else arb(str(x)) for x in row] for row in A]


# --------------------------------------------------------------------------- #
# Exact engine vs an independent oracle                                        #
# --------------------------------------------------------------------------- #
class ExactInertiaAgreesWithCharpoly(unittest.TestCase):
    CASES = {
        "identity": [[1, 0], [0, 1]],
        "signature_1_1": [[1, 0], [0, -1]],
        # No usable diagonal pivot at all: needs the 2x2 pivot path.
        "antidiagonal": [[0, 1], [1, 0]],
        "zero_2x2": [[0, 0], [0, 0]],
        "rank_one_psd": [[1, 2], [2, 4]],
        "diag_1_m1_0": [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
        "zero_diagonal_3x3": [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
        "indefinite_3x3": [[1, 2, 3], [2, 1, 4], [3, 4, 1]],
        "hilbert_3x3": [[F(1), F(1, 2), F(1, 3)],
                        [F(1, 2), F(1, 3), F(1, 4)],
                        [F(1, 3), F(1, 4), F(1, 5)]],
        "negative_definite": [[-2, 1], [1, -2]],
        "wide_scale": [[F(10 ** 8), F(1)], [F(1), F(1, 10 ** 8)]],
    }

    def test_every_case_matches_the_independent_oracle(self):
        for name, A in self.CASES.items():
            with self.subTest(name):
                r = exact_inertia(A)
                self.assertEqual(r.status, "PASS", name)
                self.assertEqual(r.signature, C.charpoly_inertia(A), name)

    def test_signature_counts_sum_to_the_dimension(self):
        for name, A in self.CASES.items():
            with self.subTest(name):
                r = exact_inertia(A)
                self.assertEqual(sum(r.signature), len(A), name)

    def test_a_matrix_with_no_usable_diagonal_pivot_still_resolves(self):
        """[[0,1],[1,0]] has inertia (1,1,0) and no nonzero diagonal entry."""
        r = exact_inertia([[0, 1], [1, 0]])
        self.assertEqual(r.signature, (1, 1, 0))
        self.assertTrue(any(p.kind == "2x2" for p in r.pivots))


# --------------------------------------------------------------------------- #
# §4 congruence / Sylvester regression                                         #
# --------------------------------------------------------------------------- #
class CongruenceInvariance(unittest.TestCase):
    SIGNATURES = ([1, 1], [1, -1], [-1, -1], [1, 0], [1, 1, -1], [1, -1, -1],
                  [1, 1, 1], [1, -1, 0], [0, 0, 0], [1, 1, -1, -1])

    def test_inertia_is_invariant_under_exact_rational_congruence(self):
        for signs in self.SIGNATURES:
            D = C.diagonal(signs)
            want = exact_inertia(D).signature
            for seed in (1, 7, 42, 99, 1234):
                with self.subTest(signs=tuple(signs), seed=seed):
                    S = C.unimodular(len(signs), seed)
                    self.assertTrue(C.is_invertible(S), "shear product must be invertible")
                    A = C.congruence(D, S)
                    self.assertEqual(exact_inertia(A).signature, want)
                    self.assertEqual(C.charpoly_inertia(A), want)

    def test_congruence_actually_changes_the_matrix(self):
        """Guard against a vacuous pass from S = I."""
        D = C.diagonal([1, -1, 1])
        moved = 0
        for seed in (1, 7, 42, 99, 1234):
            if C.congruence(D, C.unimodular(3, seed)) != D:
                moved += 1
        self.assertGreaterEqual(moved, 4, "congruences must not be near-identity")

    def test_rank_one_psd_has_a_genuine_zero_direction(self):
        for v in ([1, 2], [3, -1], [F(1, 3), F(5, 7)]):
            A = [[v[i] * v[j] for j in range(2)] for i in range(2)]
            with self.subTest(v=str(v)):
                self.assertEqual(exact_inertia(A).signature, (1, 0, 1))

    def test_signature_1_1_blocks_stay_indefinite_under_congruence(self):
        D = C.diagonal([1, -1])
        for seed in range(1, 9):
            A = C.congruence(D, C.unimodular(2, seed))
            with self.subTest(seed=seed):
                self.assertEqual(exact_inertia(A).signature, (1, 1, 0))
                self.assertLess(C.determinant(A), 0)


# --------------------------------------------------------------------------- #
# Interval engine: agreement, and failing closed                               #
# --------------------------------------------------------------------------- #
@unittest.skipIf(_flint() is None, "python-flint not installed")
class IntervalInertia(unittest.TestCase):
    def test_agrees_with_exact_on_nonsingular_congruences(self):
        for signs in ([1, 1], [1, -1], [-1, -1], [1, 1, -1], [1, -1, -1]):
            for seed in range(1, 8):
                D = C.diagonal(signs)
                A = C.congruence(D, C.unimodular(len(signs), seed))
                with self.subTest(signs=tuple(signs), seed=seed):
                    want = exact_inertia(A).signature
                    got = interval_inertia(_balls([[float(x) for x in r] for r in A]))
                    self.assertEqual(got.status, "PASS")
                    self.assertEqual(got.signature, want)

    def test_singular_matrices_must_be_inconclusive_not_guessed(self):
        """§14.3: an interval cannot witness exact zero."""
        for name, A in {"zero_2x2": [[0, 0], [0, 0]],
                        "rank_one_psd": [[1, 2], [2, 4]],
                        "diag_1_0": [[1, 0], [0, 0]]}.items():
            with self.subTest(name):
                r = interval_inertia(_balls(A))
                self.assertEqual(r.status, "INCONCLUSIVE", name)
                self.assertIsNone(r.signature)
                self.assertTrue(r.blocker)

    def test_a_box_straddling_singularity_must_be_inconclusive(self):
        r = interval_inertia(_balls([[1, 1], [1, 1]], rad="1e-9"))
        self.assertEqual(r.status, "INCONCLUSIVE")

    def test_boxes_strictly_either_side_of_a_crossing_do_resolve(self):
        arb = _flint()
        for t, want in (("1", (2, 0, 0)), ("-1", (1, 1, 0))):
            M = [[arb("1"), arb("0")], [arb("0"), arb(t, "0.5")]]
            with self.subTest(t=t):
                r = interval_inertia(M)
                self.assertEqual(r.status, "PASS")
                self.assertEqual(r.signature, want)

    def test_interval_path_never_claims_an_exact_zero(self):
        for signs in ([1, 1], [1, -1], [1, 1, -1]):
            A = C.congruence(C.diagonal(signs), C.unimodular(len(signs), 3))
            r = interval_inertia(_balls([[float(x) for x in row] for row in A]))
            with self.subTest(signs=tuple(signs)):
                self.assertEqual(r.n_zero, 0)


# --------------------------------------------------------------------------- #
# Stratification over a parameter cell                                         #
# --------------------------------------------------------------------------- #
@unittest.skipIf(_flint() is None, "python-flint not installed")
class Stratification(unittest.TestCase):
    def _family(self, shift):
        arb = _flint()
        from interval_backend import interval_box

        def fam(lo, hi):
            t = interval_box(lo, hi)
            return [[arb(1), arb(0)], [arb(0), t - arb(str(shift))]]

        return fam

    def test_a_known_crossing_is_bracketed_not_absorbed(self):
        st = certify_inertia_family(self._family("0.5"), (0.0, 1.0),
                                    subdivision_policy={"initial_cells": 4,
                                                        "max_depth": 14})
        self.assertEqual([s.signature for s in st.strata], [(1, 1, 0), (2, 0, 0)])
        self.assertTrue(st.transitions, "the crossing must be reported, not swallowed")
        t = st.transitions[0]
        self.assertLessEqual(t.lo, 0.5)
        self.assertGreaterEqual(t.hi, 0.5)
        self.assertFalse(st.is_constant)

    def test_strata_and_transitions_tile_the_cell_exactly(self):
        st = certify_inertia_family(self._family("0.5"), (0.0, 1.0),
                                    subdivision_policy={"initial_cells": 4,
                                                        "max_depth": 14})
        pieces = sorted([(s.lo, s.hi) for s in st.strata]
                        + [(t.lo, t.hi) for t in st.transitions])
        self.assertEqual(pieces[0][0], 0.0)
        self.assertEqual(pieces[-1][1], 1.0)
        for (_, h1), (l2, _) in zip(pieces, pieces[1:]):
            self.assertEqual(h1, l2, "tiling must have no gap and no overlap")

    def test_a_constant_family_yields_one_stratum_and_no_transitions(self):
        arb = _flint()
        from interval_backend import interval_box

        def fam(lo, hi):
            t = interval_box(lo, hi)
            return [[arb(2) + t, arb(0)], [arb(0), arb(-3) - t]]

        st = certify_inertia_family(fam, (0.0, 1.0),
                                    subdivision_policy={"initial_cells": 4})
        self.assertTrue(st.is_constant)
        self.assertEqual(st.signature_if_constant(), (1, 1, 0))
        self.assertEqual(st.transitions, [])


# --------------------------------------------------------------------------- #
# §11 certificate semantics                                                    #
# --------------------------------------------------------------------------- #
class CertificateSemantics(unittest.TestCase):
    def _schema(self, name):
        return json.loads((SCHEMA_DIR / name).read_text(encoding="utf-8"))

    def _cert(self, A, evidence_class="E0"):
        return build_inertia_certificate(
            exact_inertia(A), dimension=len(A), program="test", work_order="TEST",
            evidence_class=evidence_class, normalization_certificate_id="norm_test")

    def test_certificates_validate_against_the_shipped_schema(self):
        schema = self._schema("inertia_certificate.schema.json")
        for A in ([[1, 0], [0, 1]], [[1, 0], [0, -1]], [[1, 2], [2, 4]]):
            with self.subTest(A=str(A)):
                self.assertEqual(validate_against_schema(self._cert(A), schema), [])

    def test_an_inertia_certificate_never_satisfies_a_psd_requirement(self):
        """§11 is categorical and binds on the content kind, not the signature.

        Reported on PR #10: a passing WEIL_INERTIA_CERTIFICATE with (2,0,0) used
        to satisfy a PSD consumer while its own body said psd_claim: false -- the
        predicate contradicting the certificate it was reading. A definite
        inertia artifact must be refused exactly like an indefinite one; the
        consumer should be handed something that claims positivity.
        """
        definite = self._cert([[2, 0], [0, 3]])
        self.assertEqual(definite["n_negative"], 0)
        self.assertIs(definite["psd_claim"], False)
        self.assertFalse(satisfies_psd_requirement(definite))

    def test_an_indefinite_block_never_satisfies_a_psd_requirement(self):
        self.assertFalse(satisfies_psd_requirement(self._cert([[1, 0], [0, -1]])))

    def test_a_positivity_certificate_may_satisfy_a_psd_requirement(self):
        """The other half of the rule: a positivity claim is allowed to answer."""
        positivity = {
            "content_kind": "WEIL_DEGREE3_POSITIVITY_CERTIFICATE",
            "status": "PASS", "evidence_class": "E1", "psd_claim": True,
            "n_positive": 2, "n_negative": 0, "n_zero": 0, "rh_proof_claim": False,
        }
        self.assertTrue(satisfies_psd_requirement(positivity))

    def test_a_positivity_claim_still_needs_a_signature_that_backs_it(self):
        for bad in ({"n_negative": 1, "n_zero": 0},
                    {"n_negative": 0, "n_zero": None},
                    {"n_negative": None, "n_zero": 0}):
            body = {
                "content_kind": "WEIL_DEGREE3_POSITIVITY_CERTIFICATE",
                "status": "PASS", "evidence_class": "E1", "psd_claim": True,
                "n_positive": 1, "rh_proof_claim": False, **bad,
            }
            with self.subTest(**bad):
                self.assertFalse(satisfies_psd_requirement(body))

    def test_psd_licensing_requires_an_explicit_claim(self):
        """Positivity is declared by the producer, never inferred here."""
        body = {
            "content_kind": "WEIL_DEGREE3_POSITIVITY_CERTIFICATE",
            "status": "PASS", "evidence_class": "E1",
            "n_positive": 2, "n_negative": 0, "n_zero": 0, "rh_proof_claim": False,
        }
        self.assertFalse(satisfies_psd_requirement(body), "no psd_claim field")
        self.assertFalse(satisfies_psd_requirement(dict(body, psd_claim=False)))
        self.assertTrue(satisfies_psd_requirement(dict(body, psd_claim=True)))

    def test_an_inconclusive_result_never_satisfies_a_psd_requirement(self):
        arb = _flint()
        if arb is None:
            self.skipTest("python-flint not installed")
        cert = build_inertia_certificate(
            interval_inertia(_balls([[1, 2], [2, 4]])), dimension=2, program="test",
            work_order="TEST", evidence_class="E1",
            normalization_certificate_id="norm_test")
        self.assertEqual(cert["status"], "INCONCLUSIVE")
        self.assertFalse(satisfies_psd_requirement(cert))

    def test_an_e3_diagnostic_never_satisfies_a_psd_requirement(self):
        """§14.4: a floating scan cannot promote an E1 claim."""
        body = {
            "content_kind": "WEIL_DEGREE3_POSITIVITY_CERTIFICATE",
            "status": "PASS", "evidence_class": "E3", "psd_claim": True,
            "n_positive": 2, "n_negative": 0, "n_zero": 0, "rh_proof_claim": False,
        }
        self.assertFalse(satisfies_psd_requirement(body))
        self.assertTrue(satisfies_psd_requirement(dict(body, evidence_class="E1")))

    def test_a_stratification_never_satisfies_a_psd_requirement(self):
        arb = _flint()
        if arb is None:
            self.skipTest("python-flint not installed")
        from interval_backend import interval_box

        def fam(lo, hi):
            t = interval_box(lo, hi)
            return [[arb(2) + t, arb(0)], [arb(0), arb(1) + t]]

        st = certify_inertia_family(fam, (0.0, 1.0),
                                    subdivision_policy={"initial_cells": 4})
        cert = build_stratification_certificate(
            st, dimension=2, program="test", work_order="TEST", evidence_class="E1",
            normalization_certificate_id="norm_test")
        # Every stratum is positive definite here...
        self.assertEqual(st.signature_if_constant(), (2, 0, 0))
        # ...and it still must not satisfy a PSD consumer.
        self.assertFalse(satisfies_psd_requirement(cert))
        self.assertEqual(cert["content_kind"], KIND_STRATIFICATION)

    def test_no_certificate_claims_rh_or_psd(self):
        for cert in (self._cert([[1, 0], [0, 1]]), self._cert([[1, 0], [0, -1]])):
            self.assertIs(cert["rh_proof_claim"], False)
            self.assertIs(cert["psd_claim"], False)
            self.assertEqual(cert["content_kind"], KIND_INERTIA)

    def test_schema_validation_actually_rejects_bad_bodies(self):
        schema = self._schema("inertia_certificate.schema.json")
        good = self._cert([[1, 0], [0, 1]])
        self.assertEqual(validate_against_schema(good, schema), [])
        for mutate, why in (
            (lambda c: c.pop("dimension"), "missing required"),
            (lambda c: c.update(rh_proof_claim=True), "const"),
            (lambda c: c.update(status="MAYBE"), "enum"),
            (lambda c: c.update(claim_scope="something_else"), "const"),
            (lambda c: c.update(dimension=-1), "minimum"),
        ):
            bad = json.loads(json.dumps(good))
            mutate(bad)
            with self.subTest(why):
                self.assertTrue(validate_against_schema(bad, schema),
                                f"validator missed a {why} violation")


if __name__ == "__main__":
    unittest.main(verbosity=1)
