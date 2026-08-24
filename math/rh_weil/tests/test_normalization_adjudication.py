"""WO-RH-17 / WO-RH-18 regression tests.

Locks the adjudicated normalization in place:

* the adopted pole is reproduced by an **independent** real-space route;
* the parity lemma and the parity block-diagonality hold;
* the rejected ``sqrt(3)/2`` block is exactly ``(sqrt(3)/2)cosh(L/2)`` times the
  adopted one, hence agrees at ``L = log 3`` and nowhere else;
* the repository's *odd* pivot already equalled the adopted pole;
* every certificate that depended on the disputed block is quarantined and the
  PIR guard refuses to promote it;
* no adopted-path quantity depends on the fitted constant.
"""
import json
import math
import os
import sys
import unittest

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import normalization as N  # noqa: E402
import rejected_pole as RP  # noqa: E402  (archival; this file audits it)

CERT_DIR = os.path.join(ROOT, "certificates")
L_POINTS = [math.log(3.0), 1.1059498113, 1.20, math.log(4.0)]
LOG3 = math.log(3.0)


class PoleDerivationTests(unittest.TestCase):
    def test_pole_matches_independent_realspace_route(self):
        """G0_ij = int_0^L K_ij(a) 2cosh(a/2) da — the same K_ij the prime block uses."""
        try:
            import mpmath as mp
        except ImportError:  # pragma: no cover
            self.skipTest("mpmath unavailable")
        mp.mp.dps = 40
        for L in L_POINTS:
            for i in N.BASIS_NAMES:
                for j in N.BASIS_NAMES:
                    closed = N.pole_entry(i, j, L)
                    quad = float(
                        mp.quad(lambda a: N.kernel_K(i, j, float(a), L) * 2 * mp.cosh(a / 2), [0, L])
                    )
                    self.assertLess(
                        abs(closed - quad) / max(1.0, abs(closed)), 1e-12, (i, j, L)
                    )

    def test_parity_lemma(self):
        """h(L-x) = ±h(x)  =>  E^- = ± e^{-L/2} E^+."""
        for L in L_POINTS:
            for name in N.BASIS_NAMES:
                sign = 1 if N.basis_parity(name) == "even" else -1
                lhs = N.E_pm(name, L, -1)
                rhs = sign * math.exp(-L / 2) * N.E_pm(name, L, 1)
                self.assertLess(abs(lhs - rhs), 1e-12 * max(1.0, abs(rhs)), (name, L))

    def test_pole_is_parity_block_diagonal(self):
        for L in L_POINTS:
            for even in ("one", "b"):
                self.assertLess(abs(N.pole_entry(even, "q1", L)), 1e-12, (even, L))

    def test_even_block_is_rank_one(self):
        """Delta = E0+Eb- - E0-Eb+ vanishes, so both candidates are rank 1.

        This is why a determinant/rank regression could not discriminate them.
        """
        for L in L_POINTS:
            delta = N.E_pm("one", L, 1) * N.E_pm("b", L, -1) - N.E_pm("one", L, -1) * N.E_pm("b", L, 1)
            self.assertLess(abs(delta), 1e-12, L)

    def test_odd_pivot_already_matched_the_repository(self):
        """E_q1^+ = 2 e^{L/4} A(L)  =>  G0[q1,q1] = -8 A^2 (the value already shipped)."""
        for L in L_POINTS:
            A = L * math.cosh(L / 4) - 4 * math.sinh(L / 4)
            self.assertLess(abs(N.E_pm("q1", L, 1) - 2 * math.exp(L / 4) * A), 1e-12 * max(1.0, abs(A)))
            self.assertLess(abs(N.pole_entry("q1", "q1", L) - (-8 * A * A)), 1e-12 * max(1.0, A * A))


class LegacyRejectionTests(unittest.TestCase):
    def test_quotient_is_sqrt3_over_2_cosh(self):
        for L in L_POINTS:
            ratio = RP.legacy_pole_entry("one", "one", L) / N.pole_entry("one", "one", L)
            self.assertLess(abs(ratio - RP.legacy_over_adopted_ratio(L)), 1e-12, L)

    def test_agreement_only_at_log3(self):
        self.assertLess(abs(RP.legacy_over_adopted_ratio(LOG3) - 1.0), 1e-14)
        for L in (1.1059498113, 1.20, math.log(4.0)):
            self.assertGreater(abs(RP.legacy_over_adopted_ratio(L) - 1.0), 1e-4, L)

    def test_discrepancy_is_L_dependent_so_not_a_basis_change(self):
        """A change of basis is a constant congruence; it cannot produce a factor
        that varies with L. The quotient does vary, so B is not a renormalised A."""
        r1 = RP.legacy_over_adopted_ratio(1.20)
        r2 = RP.legacy_over_adopted_ratio(math.log(4.0))
        self.assertGreater(abs(r1 - r2), 1e-3)

    def test_adopted_path_carries_no_fitted_constant(self):
        """The adopted pole is a pure E^± product: rescaling h rescales it
        bilinearly, with no residual constant left over."""
        L = 1.20
        direct = N.E_pm("one", L, 1) * N.E_pm("b", L, -1) + N.E_pm("one", L, -1) * N.E_pm("b", L, 1)
        self.assertLess(abs(N.pole_entry("one", "b", L) - direct), 1e-15 * max(1.0, abs(direct)))

    def test_legacy_model_reproduces_the_shipped_block(self):
        """Production ships Candidate A; the archival model still models B.

        Before ENG-004 this test pinned ``finite_weil.g0_even_block`` against the
        rejected block, because that is what production computed. ENG-004 moved
        production onto Candidate A, so the check splits in two: production must
        now equal A, and the archival module must still reproduce the historical
        B it exists to audit.
        """
        try:
            from interval_backend import require_flint

            _flint, arb, _acb, _ctx = require_flint()
            import finite_weil as FW
            import pole
        except Exception:  # pragma: no cover
            self.skipTest("python-flint unavailable")
        for L in (LOG3, 1.20, math.log(4.0)):
            g00, g0b, gbb = FW.g0_even_block(arb(L), arb)
            for shipped, adopted in (
                (g00, pole.pole_gram_entry("one", "one", L)),
                (g0b, pole.pole_gram_entry("one", "b", L)),
                (gbb, pole.pole_gram_entry("b", "b", L)),
            ):
                self.assertLess(abs(float(shipped.mid()) - adopted),
                                1e-9 * max(1.0, abs(adopted)), L)

    def test_archival_model_still_reproduces_the_rejected_block(self):
        """``(sqrt(3)/2)(E_i^+E_j^+ + E_i^-E_j^-)`` — the historical expression."""
        import pole

        root3_over_2 = math.sqrt(3.0) / 2.0
        for L in (LOG3, 1.20, math.log(4.0), 3.5):
            for i, j in (("one", "one"), ("one", "b"), ("b", "b")):
                want = root3_over_2 * (
                    pole.laplace_plus(i, L) * pole.laplace_plus(j, L)
                    + pole.laplace_minus(i, L) * pole.laplace_minus(j, L)
                )
                self.assertAlmostEqual(RP.legacy_pole_entry(i, j, L), want,
                                       delta=1e-12 * max(1.0, abs(want)), msg=(i, j, L))


class QuarantineTests(unittest.TestCase):
    AFFECTED = list(N.QUARANTINED_CERTIFICATES)

    def _load(self, name):
        with open(os.path.join(CERT_DIR, name), encoding="utf-8") as fh:
            return json.load(fh)

    def test_affected_certificates_are_quarantined(self):
        """Every disputed certificate stays quarantined unless a work order
        explicitly released it after a rigorous regeneration (ENG-004 §4)."""
        for name in self.AFFECTED:
            cert = self._load(name)
            self.assertFalse(cert.get("rh_proof_claim", False), name)
            if name in N.RELEASED_CERTIFICATES:
                self.assertIn("quarantine_released", cert, name)
                continue
            self.assertEqual(cert.get("promotion_state"), "QUARANTINED_NORMALIZATION_ADJUDICATION", name)
            self.assertFalse(cert.get("hard_constraints_certified"), name)

    def test_quarantine_preserves_history(self):
        """Nothing is deleted and the claimed evidence class is not rewritten."""
        for name in self.AFFECTED:
            if name in N.RELEASED_CERTIFICATES:
                continue  # released after a rigorous regeneration (ENG-004 §4)
            cert = self._load(name)
            q = cert.get("quarantine", {})
            self.assertIn("prior_state", q, name)
            self.assertIn("reason", q, name)
            self.assertIsNotNone(cert.get("evidence_class"), name)

    def test_pir_guard_refuses_quarantined_certificates(self):
        import pir_bridge

        refused = {r["certificate_file"] for r in pir_bridge.refused_promotions()}
        for name in self.AFFECTED:
            if name in N.RELEASED_CERTIFICATES:
                self.assertNotIn(name, refused, name)
                continue
            self.assertIn(name, refused, name)

    def test_registry_matches_the_certificates_on_disk(self):
        """The registry is the single source of truth; guard against it drifting."""
        self.assertEqual(len(N.QUARANTINED_CERTIFICATES), 5)
        for name in N.QUARANTINED_CERTIFICATES:
            self.assertTrue(os.path.exists(os.path.join(CERT_DIR, name)), name)

    def test_writer_reasserts_quarantine_on_regeneration(self):
        """A re-run of a certify script must not be able to lift the quarantine.

        The certify_*.py entrypoints predate the adjudication and rebuild these
        bodies from the REJECTED even pole block with hard_constraints_certified
        set from their own gates. Writing such a body must come back quarantined.
        """
        import certificate_io

        for name in N.QUARANTINED_CERTIFICATES:
            body = {
                "certificate_version": "0.1",
                "evidence_class": "E1_SAMPLES_PLUS_E0_CURVATURE",
                "status": "REGENERATED_BY_LEGACY_SCRIPT",
                "hard_constraints_certified": True,   # what the old script would claim
            }
            certificate_io._enforce_quarantine(name, body)
            self.assertEqual(body["promotion_state"], N.QUARANTINE_STATE, name)
            self.assertFalse(body["hard_constraints_certified"], name)
            self.assertFalse(body["rh_proof_claim"], name)
            # the contrary claim is preserved, not deleted
            self.assertTrue(body["quarantine"]["prior_state"]["hard_constraints_certified"], name)
            self.assertEqual(body["evidence_class"], "E1_SAMPLES_PLUS_E0_CURVATURE", name)

    def test_writer_leaves_unaffected_certificates_alone(self):
        import certificate_io

        body = {"status": "REGENERATED", "hard_constraints_certified": True}
        certificate_io._enforce_quarantine("e0_exact_identities.json", body)
        self.assertNotIn("promotion_state", body)
        self.assertTrue(body["hard_constraints_certified"])

    def test_quarantine_is_idempotent_and_keeps_original_prior_state(self):
        """Re-writing an already-quarantined body must not overwrite prior_state."""
        import certificate_io

        name = N.QUARANTINED_CERTIFICATES[0]
        body = {"status": "ORIGINAL", "hard_constraints_certified": True}
        certificate_io._enforce_quarantine(name, body)
        first_prior = dict(body["quarantine"]["prior_state"])
        body["hard_constraints_certified"] = False  # now reflects the quarantine
        certificate_io._enforce_quarantine(name, body)
        self.assertEqual(body["quarantine"]["prior_state"], first_prior)

    def test_certify_scripts_write_through_the_guarded_writer(self):
        """No certify entrypoint may bypass certificate_io.write_certificate."""
        scripts = os.path.join(os.path.dirname(CERT_DIR), "scripts")
        for fn in sorted(os.listdir(scripts)):
            if not fn.startswith("certify_"):
                continue
            with open(os.path.join(scripts, fn), encoding="utf-8") as fh:
                src = fh.read()
            self.assertIn("write_certificate", src, fn)
            self.assertNotIn("out.write_text(", src, fn)

    def test_rejected_pole_scale_is_marked_in_finite_weil(self):
        """The module still shipping candidate B must say so."""
        import finite_weil

        self.assertEqual(finite_weil.POLE_EVEN_SCALE_STATUS, "REJECTED_WO_RH_17")
        self.assertIn("REJECTED", finite_weil.__doc__)

    #: Every order the WO-RH-17 adjudication quarantined.
    QUARANTINED_BY_WO_RH_17 = ["WO-RH-05"] + [f"WO-RH-{n:02d}" for n in range(9, 16)]

    def test_quarantined_orders_are_either_still_flagged_or_explicitly_recovered(self):
        """An order may leave quarantine only by saying so, never by going blank.

        WO-RH-09 was recovered by ENG-004 and the rest by ENG-005. The guard that
        matters is not that they stay quarantined — they are meant to be recovered
        — but that no order can quietly lose the marker: the status has to name
        the work order that did the recovering.
        """
        st = self._load("work_order_status.json")
        for wo in self.QUARANTINED_BY_WO_RH_17:
            if wo not in st["orders"]:
                continue
            state = st["orders"][wo]
            if state == "quarantined_pending_WO-RH-17":
                continue
            self.assertTrue(state.startswith("recovered_"),
                            f"{wo} left quarantine without saying so: {state!r}")
            self.assertRegex(state, r"^recovered_(ENG-00[45])_",
                             f"{wo} does not name the recovering work order")

    def test_pre_quarantine_claims_are_retained_verbatim(self):
        """WO-RH-17 forbids deleting the contrary evidence a recovery supersedes."""
        st = self._load("work_order_status.json")
        for wo, state in st["orders"].items():
            if str(state).startswith("recovered_"):
                self.assertIn(wo, st["pre_quarantine_orders"],
                              f"{wo} was recovered but its pre-quarantine claim is gone")


class AdjudicationCertificateTests(unittest.TestCase):
    def _cert(self, name):
        path = os.path.join(CERT_DIR, name)
        if not os.path.exists(path):
            self.skipTest(f"{name} not generated yet")
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)

    def test_adjudication_certificate_shape(self):
        c = self._cert("normalization_adjudication.json")
        for field in (
            "certificate_version", "program", "status", "active_normalization_id",
            "fourier_convention", "tilde_convention", "pole_formula", "prime_formula",
            "archimedean_formula", "crosschecks", "rh_proof_claim",
        ):
            self.assertIn(field, c, field)
        self.assertEqual(c["program"], "RH/Weil normalization adjudication")
        self.assertEqual(c["status"], "ADJUDICATED")
        self.assertIs(c["rh_proof_claim"], False)

    def test_dispositions_are_recorded(self):
        c = self._cert("normalization_adjudication.json")
        d = c["candidate_dispositions"]
        self.assertEqual(d["candidate_A_explicit_formula"]["disposition"], "ADOPTED")
        self.assertEqual(d["candidate_B_repo_sqrt3_over_2"]["disposition"], "REJECTED")

    def test_normalization_id_is_stable_and_active(self):
        c = self._cert("normalization_adjudication.json")
        self.assertEqual(c["active_normalization_id"], N.normalization_id())
        self.assertEqual(N.normalization_id(), N.normalization_id())

    def test_crosscheck_certificate_has_no_disagreements(self):
        c = self._cert("normalization_crosscheck.json")
        self.assertEqual(c["disagreements"], [], c["disagreements"][:2])
        self.assertEqual(c["status"], "AGREE")
        self.assertGreater(c["coverage"]["cells_compared"], 0)


if __name__ == "__main__":
    unittest.main()
