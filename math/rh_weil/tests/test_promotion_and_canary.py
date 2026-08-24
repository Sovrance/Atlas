"""Promotion guards, quarantine persistence and the scalar canary (ENG-004 §8/§11).

Everything here is about *refusing* to promote. A certificate reaches a claim
only when the central predicate says so, and the rules that matter are the ones
that are easy to lose accidentally:

* a quarantined certificate cannot promote itself;
* repeated writes preserve the **original** ``prior_state``;
* a rigorous certificate without the active normalization id is refused;
* a stale normalization id or a stale source hash is refused;
* an authorised explicit release works, and releases nothing else;
* the scalar canary promotes only after a real regeneration;
* every other disputed E1 stays quarantined.
"""
from __future__ import annotations

import copy
import json
import math
import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import certificate_io  # noqa: E402
import normalization as N  # noqa: E402
import promotion  # noqa: E402

CERT_DIR = ROOT / "certificates"
SCALAR = "e1_scalar_log3_log4.json"

#: Disputed E1 artifacts that ENG-004 keeps quarantined (everything but SCALAR).
STILL_QUARANTINED = tuple(c for c in N.QUARANTINED_CERTIFICATES if c != SCALAR)


def _load(name: str) -> dict:
    return json.loads((CERT_DIR / name).read_text(encoding="utf-8"))


def _rigorous_stub(**over) -> dict:
    body = {
        "evidence_class": "E1",
        "rigorous": True,
        "status": "PASS",
        "hard_constraints_certified": True,
        "promotion_state": promotion.PROMOTED_STATE,
        promotion.NORMALIZATION_ID_FIELD: promotion.active_normalization_id(),
    }
    body.update(over)
    return body


# --------------------------------------------------------------------------- #
# Active normalization id                                                      #
# --------------------------------------------------------------------------- #
class ActiveNormalizationId(unittest.TestCase):
    def test_read_from_the_adjudication_artifact_not_a_filename(self):
        artifact = _load("normalization_adjudication.json")
        self.assertEqual(promotion.active_normalization_id(),
                         artifact["active_normalization_id"])

    def test_artifact_id_matches_the_frozen_definition(self):
        ok, why = promotion.normalization_id_consistent()
        self.assertTrue(ok, why)

    def test_strict_mode_raises_when_the_artifact_is_unreadable(self):
        original = promotion.CERT_DIR
        try:
            promotion.CERT_DIR = ROOT / "certificates" / "does-not-exist"
            with self.assertRaises(promotion.NormalizationUnavailable):
                promotion.active_normalization_id(strict=True)
            self.assertIsNone(promotion.active_normalization_id())
        finally:
            promotion.CERT_DIR = original


# --------------------------------------------------------------------------- #
# The predicate                                                                #
# --------------------------------------------------------------------------- #
class PromotionPredicate(unittest.TestCase):
    def test_a_clean_rigorous_certificate_may_promote(self):
        self.assertIsNone(promotion.promotion_refusal(_rigorous_stub()))

    def test_quarantined_certificate_cannot_self_promote(self):
        body = _rigorous_stub(promotion_state=N.QUARANTINE_STATE)
        reason = promotion.promotion_refusal(body)
        self.assertIsNotNone(reason)
        self.assertIn("QUARANTINED", reason)

    def test_missing_normalization_id_blocks_e1_promotion(self):
        body = _rigorous_stub()
        del body[promotion.NORMALIZATION_ID_FIELD]
        reason = promotion.promotion_refusal(body)
        self.assertIsNotNone(reason)
        self.assertIn(promotion.NORMALIZATION_ID_FIELD, reason)

    def test_mismatched_normalization_id_blocks_e1_promotion(self):
        body = _rigorous_stub(**{promotion.NORMALIZATION_ID_FIELD: "norm_sha256_deadbeef"})
        reason = promotion.promotion_refusal(body)
        self.assertIsNotNone(reason)
        self.assertIn("stale", reason)

    def test_e1_without_hard_constraints_is_refused(self):
        body = _rigorous_stub(hard_constraints_certified=False)
        self.assertIsNotNone(promotion.promotion_refusal(body))

    def test_stale_source_hash_blocks_promotion(self):
        body = _rigorous_stub(
            dependencies={"source_hashes": {"src/pole.py": "0" * 64}}
        )
        reason = promotion.promotion_refusal(body)
        self.assertIsNotNone(reason)
        self.assertIn("stale dependency hashes", reason)

    def test_current_source_hash_passes(self):
        body = _rigorous_stub(
            dependencies={"source_hashes": promotion.source_hashes(["src/pole.py"])}
        )
        self.assertIsNone(promotion.promotion_refusal(body))

    def test_missing_dependency_file_blocks_promotion(self):
        body = _rigorous_stub(
            dependencies={"source_hashes": {"src/not_a_module.py": "0" * 64}}
        )
        self.assertIn("missing", promotion.promotion_refusal(body))

    def test_non_rigorous_certificates_need_no_normalization_id(self):
        """E0/E2/E3 artifacts are not numerical claims bound to the pole."""
        self.assertIsNone(promotion.promotion_refusal(
            {"evidence_class": "E0", "status": "REGENERATED"}))


# --------------------------------------------------------------------------- #
# Quarantine persistence                                                       #
# --------------------------------------------------------------------------- #
class QuarantinePersistence(unittest.TestCase):
    def test_repeated_writes_preserve_the_original_prior_state(self):
        name = STILL_QUARANTINED[0]
        body = {"status": "ORIGINAL", "hard_constraints_certified": True}
        certificate_io._enforce_quarantine(name, body)
        first = copy.deepcopy(body["quarantine"]["prior_state"])
        for _ in range(3):
            body["hard_constraints_certified"] = False
            certificate_io._enforce_quarantine(name, body)
        self.assertEqual(body["quarantine"]["prior_state"], first)
        self.assertTrue(first["hard_constraints_certified"])

    def test_legacy_certifier_reruns_cannot_erase_quarantine(self):
        """The pre-ENG-003 certify_*.py bodies claim E1; the writer overrides."""
        for name in STILL_QUARANTINED:
            body = {"evidence_class": "E1", "status": "REGENERATED_BY_LEGACY_SCRIPT",
                    "hard_constraints_certified": True}
            certificate_io._enforce_quarantine(name, body)
            self.assertEqual(body["promotion_state"], N.QUARANTINE_STATE, name)
            self.assertFalse(body["hard_constraints_certified"], name)
            self.assertTrue(body["quarantine"]["prior_state"]["hard_constraints_certified"])

    def test_downstream_e1_remains_quarantined_on_disk(self):
        for name in STILL_QUARANTINED:
            cert = _load(name)
            self.assertEqual(cert.get("promotion_state"), N.QUARANTINE_STATE, name)
            self.assertFalse(cert.get("hard_constraints_certified"), name)
            self.assertIsNotNone(promotion.promotion_refusal(cert), name)

    def test_released_scalar_still_fails_closed_on_an_unauthorised_write(self):
        """Release is not permanent immunity.

        ``scripts/certify_scalar_e1.py`` still exists and still writes this file
        from the old sampling path. It carries no ``quarantine_released`` marker,
        so the write boundary must re-quarantine it rather than let a weaker
        artifact sit where a promoted one was. Verified end to end by running that
        script; reproduced here on its body shape.
        """
        body = {
            "certificate_version": "0.3",
            "evidence_class": "E1_SAMPLES_PLUS_E0_CURVATURE",
            "status": "ABSOLUTE_G00_REGENERATED_PENDING_FULL_CELL_COVER",
            "hard_constraints_certified": True,
        }
        certificate_io._enforce_quarantine(SCALAR, body)
        self.assertEqual(body["promotion_state"], N.QUARANTINE_STATE)
        self.assertFalse(body["hard_constraints_certified"])

    def test_a_released_body_survives_the_quarantine_pass(self):
        """...but the authorised release is not undone on every runner pass."""
        released = _load(SCALAR)
        self.assertIn("quarantine_released", released)
        self.assertIsNone(promotion.promotion_refusal(released))
        self.assertIn(SCALAR, N.RELEASED_CERTIFICATES)
        self.assertIn(SCALAR, N.QUARANTINED_CERTIFICATES,
                      "must stay registered so unauthorised writes still fail closed")

    def test_release_requires_the_explicit_flag(self):
        """Only ``allow_quarantine_change=True`` may lift a marker."""
        import inspect

        sig = inspect.signature(certificate_io.write_certificate)
        self.assertIn("allow_quarantine_change", sig.parameters)
        self.assertFalse(sig.parameters["allow_quarantine_change"].default)


# --------------------------------------------------------------------------- #
# The scalar canary                                                            #
# --------------------------------------------------------------------------- #
class ScalarCanary(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = _load(SCALAR) if (CERT_DIR / SCALAR).exists() else None

    def test_certificate_exists(self):
        self.assertIsNotNone(self.cert)

    def test_bound_is_positive_and_uniform(self):
        self.assertEqual(self.cert["status"], "PASS", self.cert.get("status"))
        self.assertGreater(float(self.cert["certified_lower_bound"]), 0.0)
        self.assertEqual(self.cert["domain"]["L_interval"], ["log(3)", "log(4)"])

    def test_promoted_under_the_active_normalization_id(self):
        self.assertEqual(self.cert["promotion_state"], promotion.PROMOTED_STATE)
        self.assertTrue(self.cert["hard_constraints_certified"])
        self.assertEqual(self.cert[promotion.NORMALIZATION_ID_FIELD],
                         promotion.active_normalization_id())
        self.assertIsNone(promotion.promotion_refusal(self.cert))

    def test_records_what_an_e1_claim_needs(self):
        for field in ("certified_lower_bound", "domain", "backend", "precision_bits",
                      "T", "subdivision_statistics", "dependencies", "claim_scope",
                      "convexity_certificate", "tail_lemma"):
            self.assertIn(field, self.cert, field)
        self.assertEqual(self.cert["claim_scope"], "finite_dimensional_weil_compression")
        self.assertIs(self.cert["rh_proof_claim"], False)
        self.assertEqual(self.cert["pole_candidate"], "A")

    def test_no_mpmath_produced_this_certificate(self):
        self.assertIs(self.cert["mpmath_used"], False)
        self.assertIn("flint", json.dumps(self.cert["backend"]).lower())

    def test_dependency_hashes_are_current(self):
        self.assertEqual(promotion.stale_dependencies(self.cert), [])

    def test_regression_review_contains_the_notebook_minimum(self):
        """Candidate A recovers a figure Candidate B put out of reach.

        Regression evidence only -- the acceptance gate is the positive bound.
        """
        review = self.cert["regression_review"]
        self.assertEqual(review["status"], "CONTAINED", review)
        self.assertTrue(review["containing_grid_enclosures"])

    def test_subdivision_statistics_are_recorded(self):
        stats = self.cert["subdivision_statistics"]
        self.assertGreaterEqual(stats["grid_points"], 3)
        self.assertGreaterEqual(stats["tangent_certificates"], 1)
        self.assertEqual(sorted(stats["prime_powers_in_cell"]), [2, 3])


class ScalarCanaryMath(unittest.TestCase):
    """The convexity engine, checked directly rather than through the artifact."""

    @classmethod
    def setUpClass(cls):
        try:
            from interval_backend import require_flint, set_precision_bits
        except ImportError:  # pragma: no cover
            raise unittest.SkipTest("interval backend unavailable")
        try:
            _, cls.arb, cls.acb, _ = require_flint()
        except Exception:
            raise unittest.SkipTest("python-flint unavailable")
        set_precision_bits(200)

    def test_assembled_curvature_equals_the_e0_formula(self):
        """``G00'' = 4cosh(L/2) - e^{L/2}/sinh L = 2(r^3-r-1)/(sqrt r (r^2-1))``."""
        import scalar_canary as SC

        for L in (math.log(3), 1.2, 1.2828, math.log(4)):
            La = self.arb(repr(L))
            assembled = SC.curvature_from_assembly(La)
            e0 = SC.curvature_e0_formula(La)
            self.assertTrue((assembled - e0).contains(0), L)

    def test_candidate_B_cannot_reproduce_the_e0_curvature(self):
        """The ``4cosh(L/2)`` term is Candidate A's doing and nothing else's.

        Candidate B multiplies the pole by ``(sqrt(3)/2)cosh(L/2)``, so its second
        derivative differs; assembling with it breaks the certified E0 curvature.
        """
        import rejected_pole as RP
        import scalar_canary as SC

        h = 1e-5
        for L in (1.2, 1.2828, math.log(4)):
            b2 = (RP.legacy_pole_entry("one", "one", L + h)
                  + RP.legacy_pole_entry("one", "one", L - h)
                  - 2 * RP.legacy_pole_entry("one", "one", L)) / (h * h)
            a2 = float(SC.pole.pole_scalar_g00_second_derivative(self.arb(repr(L))))
            self.assertGreater(abs(b2 - a2), 1e-3 * abs(a2), L)

    def test_tail_is_non_negative_and_small(self):
        import scalar_canary as SC

        T = 200_000
        self.assertGreater(float(SC.h_plus_at(T, self.arb, self.acb).lower()), 0.0)
        self.assertLess(float(SC.tail_bound(T, self.arb, self.acb).upper()), 1e-3)

    def test_lemma_A_constant_bounds_the_series(self):
        import scalar_canary as SC

        ok, rows = SC.lemma_A_numeric_check(self.arb)
        self.assertTrue(ok, rows)

    def test_prime_breakpoints_are_the_cell_endpoints(self):
        """Why ``Gp'' = 0`` inside the cell."""
        import scalar_canary as SC

        self.assertEqual([q for q, _ in SC.prime_powers_below(1.2)], [2, 3])
        self.assertEqual([q for q, _ in SC.prime_powers_below(math.log(4) - 1e-12)], [2, 3])


if __name__ == "__main__":
    unittest.main()
