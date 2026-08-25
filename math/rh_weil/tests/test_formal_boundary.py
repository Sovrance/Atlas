#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §10/§12 — the formal boundary and what it does not license.

The interesting assertions here are negative ones. A formal theorem is a strong
thing to have and an easy thing to over-read, so most of what follows checks
that having one changes nothing it should not change: no numeric warrant
appears, no PSD requirement is satisfied, no certificate's evidence class moves,
and nothing in the manifest is about the Riemann hypothesis.
"""
from __future__ import annotations

import json
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import formal_evidence as F  # noqa: E402
import promotion  # noqa: E402
from inertia.certificate import satisfies_psd_requirement  # noqa: E402

CERT_DIR = ROOT / "certificates"
FORMAL = ROOT / "formal"
MANIFEST = FORMAL / "manifests" / "theorem_manifest.json"
CERT_NAME = "formal_theorem_certificate.json"


def _load(name):
    p = CERT_DIR / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


class Manifest(unittest.TestCase):
    def setUp(self):
        self.m = F.load_manifest()
        if self.m is None:
            self.skipTest("no theorem manifest")

    def test_it_pins_an_exact_toolchain_and_mathlib_commit(self):
        self.assertRegex(self.m["lean_toolchain"], r"^leanprover/lean4:")
        self.assertRegex(self.m["mathlib_commit"], r"^[0-9a-f]{40}$")

    def test_the_lakefile_pins_a_commit_not_a_branch(self):
        text = (FORMAL / "lakefile.toml").read_text(encoding="utf-8")
        self.assertIsNotNone(
            re.search(r'name\s*=\s*"mathlib".*?rev\s*=\s*"[0-9a-f]{40}"', text, re.S),
            "mathlib must be pinned to a 40-hex commit; a floating dependency "
            "would let the meaning of a proved theorem change with no commit here",
        )

    def test_manifest_id_is_the_hash_of_its_own_body(self):
        sys.path.insert(0, str(ROOT / "scripts"))
        import check_formal_manifest as C

        self.assertEqual(self.m["manifest_id"], C.canonical_id(self.m))

    def test_every_theorem_uses_only_the_enumerated_axioms(self):
        allowed = set(self.m["allowed_axioms"])
        self.assertEqual(allowed, {"Classical.choice", "Quot.sound", "propext"})
        for t in self.m["theorems"]:
            self.assertLessEqual(set(t["axioms"]), allowed, t["id"])

    def test_it_makes_no_rh_proof_claim(self):
        self.assertIs(self.m["rh_proof_claim"], False)
        blob = json.dumps(self.m).lower()
        for word in ("riemann", "rh is true", "proves rh"):
            self.assertNotIn(word, blob)

    def test_unproved_statements_carry_no_warrant(self):
        # The general rank-trace inequality is recorded, not proved. Recording
        # it is what stops it being quietly assumed; a warrant on it would undo
        # exactly that.
        self.assertTrue(self.m["unproved_statements"])
        for u in self.m["unproved_statements"]:
            self.assertIsNone(u["warrant"], u["id"])
            self.assertEqual(u["status"], "EXTERNAL_THEOREM_PENDING_FORMAL_PROOF")

    def test_the_general_rank_trace_statement_has_no_proof_in_the_library(self):
        # It is a `def ... : Prop`. If it ever becomes a `theorem`, this test
        # fails and the manifest has to be updated deliberately.
        text = (FORMAL / "AtlasRH" / "RankTrace.lean").read_text(encoding="utf-8")
        self.assertIn("def RankTraceGeneralStatement : Prop", text)
        self.assertNotIn("theorem RankTraceGeneralStatement", text)


class NoSorryNoLocalAxioms(unittest.TestCase):
    def setUp(self):
        sys.path.insert(0, str(ROOT / "scripts"))
        import check_formal_manifest as C

        self.C = C

    def test_no_sorry_in_the_promoted_proof_library(self):
        self.assertEqual(self.C.check_no_sorry(), [])

    def test_no_project_local_axiom_declarations(self):
        self.assertEqual(self.C.check_no_local_axioms(), [])

    def test_the_statement_layer_imports_no_proof_module(self):
        self.assertEqual(self.C.check_statements_imports(), [])

    def test_the_scanner_actually_fires(self):
        # A gate nobody has watched fail is a gate nobody should trust.
        code = self.C.strip_lean_comments(
            "/- the word sorry in prose -/\ntheorem t : True := by sorry\n"
        )
        self.assertNotIn("sorry", code.splitlines()[0])
        self.assertIn("sorry", code.splitlines()[1])


class FormalCertificate(unittest.TestCase):
    def setUp(self):
        self.cert = _load(CERT_NAME)
        if self.cert is None:
            self.skipTest("no formal certificate")

    def test_it_is_internally_consistent(self):
        self.assertEqual(F.formal_certificate_problems(self.cert), [])

    def test_it_carries_no_numeric_warrant(self):
        self.assertIsNone(self.cert["numeric_warrant"])
        self.assertEqual(self.cert["logical_implication_warrant"], "FORMAL")

    def test_it_never_satisfies_a_psd_requirement(self):
        # It proves that positive bounds imply definiteness. It does not say the
        # bounds hold; only the E1 certificate says that.
        self.assertFalse(satisfies_psd_requirement(self.cert))

    def test_it_is_not_a_rigorous_numeric_artifact(self):
        self.assertFalse(promotion.is_rigorous(self.cert))

    def test_it_may_be_promoted(self):
        self.assertIsNone(promotion.promotion_refusal(self.cert))

    def test_its_lean_sources_are_current(self):
        self.assertEqual(promotion.stale_dependencies(self.cert), [])

    def test_it_binds_the_live_manifest(self):
        self.assertEqual(self.cert["formal_manifest_id"], F.manifest_id())


class BackingDoesNotUpgrade(unittest.TestCase):
    """§12: a formal implication never converts interval evidence to FORMAL."""

    def test_backed_certificates_keep_their_numeric_warrant(self):
        for name, ids in F.FORMAL_BACKING.items():
            cert = _load(name)
            if cert is None:
                continue
            declared = cert.get("evidence_class")
            block = F.formal_block(name, declared)
            self.assertIsNotNone(block, name)
            self.assertEqual(block["numeric_warrant"], declared, name)
            self.assertEqual(block["logical_implication_warrant"], "FORMAL", name)
            # And the certificate on disk is untouched by any of this.
            self.assertEqual(cert.get("evidence_class"), declared, name)
            self.assertNotIn("logical_implication_warrant", cert, name)

    def test_every_backing_claim_names_a_proved_theorem(self):
        proved = set(F.proved_theorem_ids())
        if not proved:
            self.skipTest("no theorem manifest")
        for name, ids in F.FORMAL_BACKING.items():
            self.assertLessEqual(set(ids), proved, name)

    def test_a_backing_claim_on_an_unproved_theorem_is_refused(self):
        original = dict(F.FORMAL_BACKING)
        try:
            F.FORMAL_BACKING["fake_certificate.json"] = ("rank_trace_general",)
            with self.assertRaises(ValueError):
                F.formal_block("fake_certificate.json", "E1")
        finally:
            F.FORMAL_BACKING.clear()
            F.FORMAL_BACKING.update(original)

    def test_an_unbacked_certificate_gains_nothing(self):
        self.assertIsNone(F.formal_block("e1_scalar_log3_log4.json", "E1"))


class ExternalZeta23(unittest.TestCase):
    """§11: the external reference stays a reference until it deliberately isn't.

    The failure mode this guards is the one §11 names: an imported theorem
    quietly becoming provenance-free E0. So the assertions are about the
    record's honesty, not about its content.
    """

    PATH = ROOT / "external" / "zeta23" / "theorem_manifest.json"

    def setUp(self):
        if not self.PATH.exists():
            self.skipTest("no zeta23 reference manifest")
        self.m = json.loads(self.PATH.read_text(encoding="utf-8"))

    def test_it_is_reference_only_and_carries_no_warrant(self):
        self.assertEqual(self.m["dependency_status"], "REFERENCE_ONLY")
        self.assertIsNone(self.m["warrant"])
        self.assertIs(self.m["rh_proof_claim"], False)
        for t in self.m["theorems"]:
            self.assertFalse(t["imported_by_atlas"], t["id"])

    def test_the_upstream_commit_is_pinned_exactly(self):
        self.assertRegex(self.m["upstream"]["commit"], r"^[0-9a-f]{40}$")
        for t in self.m["theorems"]:
            self.assertEqual(t["upstream_commit"], self.m["upstream"]["commit"], t["id"])

    def test_the_license_chain_is_recorded(self):
        up = self.m["upstream"]
        self.assertEqual(up["license"], "Apache-2.0")
        self.assertTrue(up["copyright"])
        self.assertTrue(up["notice"])

    def test_source_hashes_are_labelled_as_source_hashes(self):
        # A source-text hash is a pin. Calling it a statement comparison would
        # be the exact overclaim §11 is written against.
        self.assertEqual(self.m["hash_kind"], "source_text")
        for t in self.m["theorems"]:
            self.assertRegex(t["statement_source_hash"], r"^sha256:[0-9a-f]{64}$")
            self.assertIsNone(t["elaborated_statement_hash"], t["id"])
            self.assertIsNone(t["axiom_report"], t["id"])

    def test_toolchain_incompatibility_is_stated_not_glossed(self):
        self.assertIs(self.m["toolchain_compatible"], False)
        self.assertNotEqual(
            self.m["upstream"]["mathlib_commit"], self.m["atlas"]["mathlib_commit"]
        )

    def test_the_atlas_pins_in_the_record_match_the_live_project(self):
        live = F.load_manifest()
        if live is None:
            self.skipTest("no theorem manifest")
        self.assertEqual(self.m["atlas"]["lean_toolchain"], live["lean_toolchain"])
        self.assertEqual(self.m["atlas"]["mathlib_commit"], live["mathlib_commit"])

    def test_no_upstream_source_is_vendored(self):
        self.assertEqual(self.m["source_copied"], [])
        for lean in (ROOT / "formal").rglob("*.lean"):
            if ".lake" in lean.parts:
                continue
            text = lean.read_text(encoding="utf-8")
            self.assertNotIn("import Zeta23", text, str(lean))

    def test_the_lakefile_requires_nothing_from_upstream(self):
        text = (ROOT / "formal" / "lakefile.toml").read_text(encoding="utf-8")
        self.assertNotIn("zeta-23-lean", text)
        self.assertNotIn("Zeta23", text)


class PirVisibility(unittest.TestCase):
    """§12: the distinction must be visible in PIR, not merely true."""

    def setUp(self):
        sys.path.insert(0, str(ROOT / "src"))
        import pir_bridge

        self.bridge = pir_bridge
        if not pir_bridge.available():
            self.skipTest("pir package unavailable")

    def test_the_formal_content_kind_is_registered(self):
        self.assertIn(F.KIND_FORMAL, self.bridge.CONTENT_KINDS)

    def test_facts_separate_the_two_warrants(self):
        facts = self.bridge.certs_to_facts()
        by_file = {f.content["certificate_file"]: f for f in facts}
        formal = by_file.get(CERT_NAME)
        if formal is None:
            self.skipTest("formal certificate not promoted")
        self.assertIsNone(formal.content["numeric_warrant"])
        self.assertEqual(formal.content["logical_implication_warrant"], "FORMAL")
        self.assertFalse(formal.content["satisfies_psd_requirement"])
        for name in F.FORMAL_BACKING:
            f = by_file.get(name)
            if f is None:
                continue
            self.assertEqual(
                f.content["numeric_warrant"], f.content["evidence_class_declared"], name
            )
            self.assertEqual(f.content["logical_implication_warrant"], "FORMAL", name)

    def test_no_fact_claims_an_rh_proof(self):
        for f in self.bridge.certs_to_facts():
            self.assertIs(f.content["rh_proof_claim"], False)


if __name__ == "__main__":
    unittest.main(verbosity=2)
