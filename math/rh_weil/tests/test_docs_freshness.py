#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §14 — documentation freshness as a testable invariant.

`scripts/check_docs.py` is a gate, and a gate nobody has watched fail is a gate
nobody should trust. These tests exercise its scanners on synthetic input, so a
future edit that loosens one of them fails here rather than silently passing
everything.

The value matcher gets the most attention. It deliberately accepts a rounded
quote of a certified number, which means it could also accept a *wrong* one if
its tolerance were sloppy — and a README quietly disagreeing with a certificate
is exactly the defect this work order exists to close.
"""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import check_docs as C  # noqa: E402


class Status(unittest.TestCase):
    def setUp(self):
        self.status = C.load_status()

    def test_docs_status_agrees_with_the_machine_status(self):
        self.assertEqual(C.check_status_agreement(self.status), [])

    def test_it_names_the_current_and_previous_work_orders(self):
        self.assertEqual(self.status["current_work_order"], "ATLAS-RH-ENG-007")
        self.assertEqual(self.status["latest_completed_work_order"], "ATLAS-RH-ENG-006")

    def test_every_canonical_doc_exists_and_is_not_marked_historical(self):
        for doc in self.status["canonical_docs"]:
            path = C.resolve(doc)
            self.assertTrue(path.exists(), doc)
            self.assertFalse(C.is_historical(path, self.status), doc)

    def test_every_historical_doc_exists_and_is_marked(self):
        for doc in self.status["historical_docs"]:
            path = C.resolve(doc)
            self.assertTrue(path.exists(), doc)
            self.assertTrue(C.is_historical(path, self.status), doc)

    def test_it_makes_no_rh_proof_claim(self):
        self.assertIs(self.status["rh_proof_claim"], False)


class ScannersFire(unittest.TestCase):
    """Each scanner, on synthetic input, in both directions."""

    def setUp(self):
        self.status = C.load_status()
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, text: str) -> Path:
        p = self.dir / "doc.md"
        p.write_text(text, encoding="utf-8")
        return p

    def test_superseded_scanner_catches_a_stale_claim(self):
        p = self.write("# doc\n\nWO-RH-05 is still open, so do not proceed.\n")
        problems = C.check_superseded("doc.md", p, self.status)
        self.assertTrue(any("wo_rh_05_open" in s for s in problems), problems)

    def test_superseded_scanner_catches_a_degree_three_block(self):
        p = self.write("# doc\n\nDo not start degree 3 yet.\n")
        problems = C.check_superseded("doc.md", p, self.status)
        self.assertTrue(any("degree3_blocked" in s for s in problems), problems)

    def test_a_fenced_code_block_is_not_prose(self):
        p = self.write("# doc\n\n```\nWO-RH-05 is still open\n```\n")
        self.assertEqual(C.check_superseded("doc.md", p, self.status), [])

    def test_a_marked_quotation_region_is_exempt(self):
        p = self.write(
            "# doc\n\n"
            f"{C.QUOTE_OPEN}\nWO-RH-05 is still open\n{C.QUOTE_CLOSE}\n"
        )
        self.assertEqual(C.check_superseded("doc.md", p, self.status), [])

    def test_an_unbalanced_quotation_marker_does_not_silently_exempt(self):
        # Losing the closing marker must not turn the rest of the file into a
        # blanket exemption.
        p = self.write(f"# doc\n\n{C.QUOTE_OPEN}\nWO-RH-05 is still open\n")
        self.assertTrue(C.check_superseded("doc.md", p, self.status))

    def test_quotation_regions_are_counted(self):
        text = f"a\n{C.QUOTE_OPEN}\nx\n{C.QUOTE_CLOSE}\nb\n{C.QUOTE_OPEN}\ny\n{C.QUOTE_CLOSE}\n"
        _, used = C.strip_quoted_regions(text)
        self.assertEqual(used, 2)

    def test_boundary_scanner_catches_a_doc_with_no_claim_boundary(self):
        p = self.write("# doc\n\nNothing about scope here.\n")
        self.assertTrue(C.check_boundary("doc.md", p, self.status))

    def test_boundary_scanner_accepts_a_doc_that_states_it(self):
        p = self.write("# doc\n\nNo RH proof claim is made.\n")
        self.assertEqual(C.check_boundary("doc.md", p, self.status), [])

    def test_link_scanner_catches_a_broken_local_link(self):
        p = self.write("# doc\n\nSee [x](./does-not-exist.md).\n")
        self.assertTrue(C.check_links("doc.md", p))

    def test_link_scanner_ignores_external_links_and_anchors(self):
        p = self.write("# doc\n\n[a](https://example.invalid) [b](#section)\n")
        self.assertEqual(C.check_links("doc.md", p), [])


class ValueMatcher(unittest.TestCase):
    """The certified-value check: rounding is fine, a different number is not."""

    def test_the_exact_value_matches(self):
        self.assertTrue(C._prefix_match("3.4251152511218656e-06",
                                        "bound 3.4251152511218656e-06 on the cell"))

    def test_a_rounded_quote_matches(self):
        self.assertTrue(C._prefix_match("1.073120529992708e-06", "det >= 1.0731e-06"))

    def test_a_different_number_does_not_match(self):
        # One part in 1e3 off: close enough to look right in a table, far enough
        # to be a different certified bound.
        self.assertFalse(C._prefix_match("1.073120529992708e-06", "det >= 1.0742e-06"))

    def test_a_wrong_exponent_does_not_match(self):
        self.assertFalse(C._prefix_match("1.073120529992708e-06", "det >= 1.0731e-05"))

    def test_absent_value_does_not_match(self):
        self.assertFalse(C._prefix_match("3.4251152511218656e-06", "no numbers here"))

    def test_every_recorded_certified_value_is_live(self):
        status = C.load_status()
        for spec in status["certified_values"]:
            value, err = C.cert_value(spec)
            self.assertIsNone(err, err)
            self.assertIsNotNone(value, spec)


class AtlasVersionDetection(unittest.TestCase):
    def test_the_current_atlas_is_detected_not_hardcoded(self):
        status = C.load_status()
        current = C.current_atlas_doc(status)
        self.assertIsNotNone(current)
        others = sorted(C.REPO.glob(status["root_atlas_doc_glob"]))
        self.assertGreater(len(others), 1, "expected several atlas versions on disk")
        # Whatever it picked must be the newest by version, not by name order.
        import re

        def version(p):
            m = re.search(r"v(\d+)\.(\d+)", p.name)
            return tuple(int(g) for g in m.groups())

        self.assertEqual(version(current), max(version(p) for p in others))

    def test_the_root_readme_links_the_detected_current_atlas(self):
        self.assertEqual(C.check_root_readme(C.load_status()), [])


class WorkOrderStatus(unittest.TestCase):
    def setUp(self):
        self.wo = json.loads(
            (ROOT / "certificates" / "work_order_status.json").read_text(encoding="utf-8")
        )

    def test_it_records_eng007_as_current(self):
        self.assertEqual(self.wo["current_work_order"], "ATLAS-RH-ENG-007")

    def test_the_orders_this_work_order_closed_are_recorded(self):
        for order in ("WO-RH-37", "WO-RH-38", "WO-RH-39", "WO-RH-41",
                      "WO-RH-42", "WO-RH-43", "WO-RH-44", "WO-RH-45", "WO-RH-46"):
            self.assertIn(order, self.wo["orders"], order)
            self.assertTrue(str(self.wo["orders"][order]).startswith("done"), order)

    def test_wo_rh_40_is_recorded_as_partial(self):
        # The general rank-trace case is not proved, and the status must say so.
        self.assertTrue(str(self.wo["orders"]["WO-RH-40"]).startswith("partial"))

    def test_the_pre_quarantine_values_are_still_there(self):
        # WO-RH-17 forbids deleting contrary evidence.
        self.assertIn("WO-RH-05", self.wo["pre_quarantine_orders"])

    def test_it_makes_no_rh_proof_claim(self):
        self.assertIs(self.wo["rh_proof_claim"], False)


if __name__ == "__main__":
    unittest.main(verbosity=2)
