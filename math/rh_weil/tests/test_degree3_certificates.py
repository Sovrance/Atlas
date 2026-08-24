#!/usr/bin/env python3
"""ATLAS-RH-ENG-006 §9/§10/§11 — the degree-3 artifacts and their semantics.

Checks the certificates that were actually emitted, not the ones that were
hoped for. §9 allows three outcomes and §17 counts two of them as success, so
these tests are written to accept a positivity result or an inertia
stratification and to insist only that whichever arrived is internally
consistent, correctly scoped, and honest about what it licenses.
"""
from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import promotion  # noqa: E402
from inertia.certificate import (  # noqa: E402
    KIND_INERTIA,
    KIND_STRATIFICATION,
    satisfies_psd_requirement,
    validate_against_schema,
)

CERT_DIR = ROOT / "certificates"
SCHEMA_DIR = ROOT / "inertia" / "schemas"
E1_ALTERNATIVES = ("e1_degree3_odd_positivity_log3_log4.json",
                   "e1_degree3_odd_inertia_log3_log4.json")


def _load(name):
    path = CERT_DIR / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _e1():
    for name in E1_ALTERNATIVES:
        cert = _load(name)
        if cert is not None:
            return name, cert
    return None, None


class DegreeThreeE1(unittest.TestCase):
    def setUp(self):
        self.name, self.cert = _e1()
        if self.cert is None:
            self.skipTest("degree-3 E1 certificate not generated")
        if self.cert.get("quick_mode"):
            self.skipTest("degree-3 certificate is a quick smoke test")

    def test_exactly_one_of_the_two_outcomes_was_emitted(self):
        present = [n for n in E1_ALTERNATIVES if (CERT_DIR / n).exists()]
        self.assertEqual(len(present), 1,
                         f"the outcome must be unambiguous, found {present}")

    def test_the_filename_matches_the_outcome_it_records(self):
        outcome = self.cert["outcome"]
        if outcome == "A_POSITIVE_DEFINITE":
            self.assertEqual(self.name, E1_ALTERNATIVES[0])
        elif outcome == "B_INERTIA_STRATIFICATION":
            self.assertEqual(self.name, E1_ALTERNATIVES[1])
        else:
            self.fail(f"unexpected outcome {outcome}")

    def test_it_is_promoted_and_binds_the_active_normalization(self):
        self.assertEqual(self.cert["status"], "PASS")
        self.assertEqual(self.cert["promotion_state"], promotion.PROMOTED_STATE)
        self.assertIsNone(promotion.promotion_refusal(self.cert))
        self.assertEqual(self.cert[promotion.NORMALIZATION_ID_FIELD],
                         promotion.active_normalization_id())

    def test_the_signature_sums_to_the_block_dimension(self):
        i = self.cert["inertia"]
        if i["n_positive"] is None:
            self.skipTest("inertia is not constant on the cell")
        self.assertEqual(i["n_positive"] + i["n_negative"] + i["n_zero"], 2)

    def test_a_positive_definite_outcome_carries_both_named_bounds(self):
        """§9 Outcome A names O1 > 0 and det > 0 explicitly."""
        if self.cert["outcome"] != "A_POSITIVE_DEFINITE":
            self.skipTest("not outcome A")
        bounds = self.cert["uniform_bounds"]
        self.assertEqual(set(bounds), {"O1", "det_odd3"})
        for key, b in bounds.items():
            with self.subTest(key):
                self.assertGreater(float(b["certified_lower_bound"]), 0.0)
                self.assertIn("no convexity or monotonicity assumed",
                              b["cover"]["topology_proved"])

    def test_the_certified_cell_is_the_whole_closed_cell(self):
        lo, hi = (float(x) for x in self.cert["certified_cell"])
        self.assertAlmostEqual(lo, math.log(3.0), places=15)
        self.assertAlmostEqual(hi, math.log(4.0), places=15)

    def test_the_strata_and_transitions_tile_the_certified_cell(self):
        strat = self.cert["inertia_stratification"]
        pieces = [[float(x) for x in s["interval"]] for s in strat["strata"]]
        pieces += [[float(x) for x in t["interval"]]
                   for t in strat["transition_regions"]]
        pieces.sort()
        lo, hi = (float(x) for x in self.cert["certified_cell"])
        self.assertEqual(pieces[0][0], lo)
        self.assertEqual(pieces[-1][1], hi)
        for (_, h1), (l2, _) in zip(pieces, pieces[1:]):
            self.assertEqual(h1, l2, "tiling must have no gap and no overlap")

    def test_the_scan_is_recorded_as_a_clue_not_a_warrant(self):
        self.assertIn("warrant", self.cert["strategy_chosen_from"])
        scan = _load("e3_degree3_odd_scan_log3_log4.json")
        self.assertIsNotNone(scan)
        self.assertEqual(scan["evidence_class"], "E3")
        self.assertFalse(scan["hard_constraints_certified"])

    def test_no_eigenvalue_solver_is_claimed_in_the_method(self):
        self.assertIn("no eigenvalue solver", self.cert["method"])
        self.assertIn("no termwise PSD domination", self.cert["method"])

    def test_the_nested_inertia_object_never_satisfies_a_psd_requirement(self):
        """§11, stated precisely: the rule binds the *inertia* certificate.

        The nested stratification is an inertia artifact and must refuse a PSD
        consumer whatever its signature says -- that is the whole point of the
        rule, and it has to hold here where the signature is in fact (2,0,0).
        """
        nested = self.cert["inertia_stratification"]
        self.assertIn(nested["content_kind"], (KIND_INERTIA, KIND_STRATIFICATION))
        self.assertFalse(satisfies_psd_requirement(nested))
        self.assertIs(nested["psd_claim"], False)

    def test_psd_licensing_of_the_outer_artifact_matches_what_it_proved(self):
        """An Outcome-A artifact is a positivity certificate and may say so."""
        i = self.cert["inertia"]
        allowed = satisfies_psd_requirement(self.cert)
        if self.cert["outcome"] == "A_POSITIVE_DEFINITE":
            self.assertEqual(i["n_negative"], 0)
            self.assertTrue(allowed,
                            "a certified positive definite block must satisfy a "
                            "PSD requirement")
            self.assertIs(self.cert["psd_claim"], True)
            self.assertNotIn(self.cert["content_kind"],
                             (KIND_INERTIA, KIND_STRATIFICATION))
        else:
            self.assertFalse(allowed)
            self.assertIs(self.cert["psd_claim"], False)

    def test_the_top_level_signature_matches_the_nested_one(self):
        i = self.cert["inertia"]
        for key in ("n_positive", "n_negative", "n_zero"):
            self.assertEqual(self.cert[key], i[key], key)

    def test_it_makes_no_rh_claim_and_stays_in_scope(self):
        self.assertIs(self.cert["rh_proof_claim"], False)
        self.assertEqual(self.cert["claim_scope"],
                         "finite_dimensional_weil_compression")

    def test_the_stratification_body_validates_against_the_schema(self):
        schema = json.loads((SCHEMA_DIR / "inertia_certificate.schema.json")
                            .read_text(encoding="utf-8"))
        self.assertEqual(
            validate_against_schema(self.cert["inertia_stratification"], schema), [])


class DegreeThreeMoments(unittest.TestCase):
    def setUp(self):
        self.cert = _load("e1_degree3_odd_moments_log3_log4.json")
        if self.cert is None:
            self.skipTest("moment certificate not generated")
        if self.cert.get("quick_mode"):
            self.skipTest("moment certificate is a quick smoke test")

    def test_it_validates_against_the_spectral_moment_schema(self):
        schema = json.loads((SCHEMA_DIR / "spectral_moment_certificate.schema.json")
                            .read_text(encoding="utf-8"))
        self.assertEqual(validate_against_schema(self.cert, schema), [])

    def test_every_point_reports_all_four_moments(self):
        for row in self.cert["points"]:
            ms = row["moment_analysis"]["moments"]
            with self.subTest(row["label"]):
                self.assertEqual(set(ms), {"m1", "m2", "m3", "m4"})
                for k, v in ms.items():
                    self.assertLessEqual(float(v["lo"]), float(v["hi"]))

    def test_even_moments_are_non_negative_and_no_sanity_violation_fired(self):
        for row in self.cert["points"]:
            ms = row["moment_analysis"]["moments"]
            with self.subTest(row["label"]):
                self.assertGreaterEqual(float(ms["m2"]["lo"]), 0.0)
                self.assertGreaterEqual(float(ms["m4"]["lo"]), 0.0)
                self.assertEqual(row["moment_analysis"]["sanity_violations"], [])

    def test_moments_recover_the_inertia_at_n_equals_2(self):
        """The 2x2 case is where m1 and m2 invert to the spectrum."""
        for row in self.cert["points"]:
            q = next(q for q in row["moment_analysis"]["b1_queries"]
                     if q["query"] == "inertia_determined_by_moments")
            with self.subTest(row["label"]):
                self.assertEqual(q["status"], "CONCLUSIVE")
                self.assertTrue(q["determined"])

    def test_insufficient_information_is_recorded_where_it_applies(self):
        """§6: it is a certified outcome, so it has to actually appear."""
        statuses = {q["status"] for row in self.cert["points"]
                    for q in row["moment_analysis"]["b1_queries"]}
        self.assertIn("INSUFFICIENT_INFORMATION", statuses,
                      "the adapter must not be claiming more than moments give")

    def test_rank_trace_names_its_theorem_and_hypotheses(self):
        for row in self.cert["points"]:
            rt = row["rank_trace"]
            with self.subTest(row["label"]):
                self.assertEqual(rt["theorem_id"], "rank_trace_hs_v1")
                self.assertIn("hypotheses", rt)
                for name, h in rt["hypotheses"].items():
                    self.assertIn("statement", h)
                    self.assertIn("verified", h)
                if rt["status"] == "PASS":
                    for name, h in rt["hypotheses"].items():
                        self.assertTrue(h["verified"], f"{name} must be verified")

    def test_the_rank_trace_bound_never_exceeds_the_true_rank(self):
        """A PD 2x2 has rank 2; a bound above that would be a broken theorem."""
        for row in self.cert["points"]:
            rt = row["rank_trace"]
            if rt["status"] != "PASS":
                continue
            with self.subTest(row["label"]):
                self.assertLessEqual(rt["result"]["certified_rank_lower_bound"], 2)

    def test_the_normalization_hypothesis_is_checked_not_declared(self):
        for row in self.cert["points"]:
            h = row["rank_trace"]["hypotheses"]["shared_normalization"]
            with self.subTest(row["label"]):
                self.assertIn("requirement", h["evidence"])
                self.assertIn("lambda_max_upper", h["evidence"])

    def test_it_makes_no_rh_claim(self):
        self.assertIs(self.cert["rh_proof_claim"], False)


if __name__ == "__main__":
    unittest.main(verbosity=1)
