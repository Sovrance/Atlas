"""Production must not import the REJECTED Candidate-B pole (ENG-004 §1).

CI fails here if the rejected block creeps back into a production module. The
archival module itself, the adjudication scripts and the tests are the only
legitimate places it may appear -- they *audit* it, they do not compute with it.
"""
from __future__ import annotations

import ast
import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SRC))

#: The archival module. Nothing under src/ but this file may import it.
REJECTED_MODULE = "rejected_pole"

#: Scripts whose whole job is to audit the rejected candidate.
AUDIT_SCRIPTS = {"derive_normalization.py", "run_normalization_crosscheck.py"}

#: Textual fingerprints of the rejected calibration.
REJECTED_FINGERPRINTS = ("sqrt(3)/2", "sqrt(3) / 2", "arb(3).sqrt() / 2", "sqrt(3)/4")


def _imported_names(path: Path) -> set[str]:
    """Top-level module names imported by ``path``, including inside functions."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module and node.level == 0:
                names.add(node.module.split(".")[0])
    return names


def _executable_source(path: Path) -> str:
    """``path``'s source with every string constant and comment removed."""
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    segments = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            seg = ast.get_source_segment(text, node)
            if seg:
                segments.append(seg)
    for seg in sorted(segments, key=len, reverse=True):
        text = text.replace(seg, '""')
    return "\n".join(
        ln.split("#", 1)[0] for ln in text.splitlines() if not ln.strip().startswith("#")
    )


def _production_modules() -> list[Path]:
    return sorted(
        p for p in SRC.glob("*.py")
        if p.name != f"{REJECTED_MODULE}.py" and not p.name.startswith("__")
    )


class ProductionImportScan(unittest.TestCase):
    def test_no_production_module_imports_the_rejected_pole(self):
        offenders = [
            p.name for p in _production_modules()
            if REJECTED_MODULE in _imported_names(p)
        ]
        self.assertEqual(
            offenders, [],
            f"production modules import the REJECTED Candidate-B pole: {offenders}. "
            "Route the pole through src/pole.py (Candidate A) instead.",
        )

    def test_production_modules_do_not_spell_the_rejected_scale(self):
        """A copied-in ``sqrt(3)/2`` would evade the import scan.

        Only *executable* code is scanned. Production modules are expected to
        describe the rejected calibration in prose -- the quarantine reason and
        the adjudication notes both quote it -- so every string constant and
        comment is blanked before the fingerprints are matched.
        """
        offenders = []
        for path in _production_modules():
            if any(fp in _executable_source(path) for fp in REJECTED_FINGERPRINTS):
                offenders.append(path.name)
        self.assertEqual(
            offenders, [],
            f"the rejected sqrt(3)/2 pole scale is computed in production: {offenders}",
        )

    def test_the_archival_module_is_reachable_for_audits(self):
        """The rejected candidate is archived, not deleted -- audits still need it."""
        import rejected_pole

        self.assertEqual(rejected_pole.LEGACY_STATUS, "REJECTED_FITTED_CALIBRATION")
        self.assertTrue(callable(rejected_pole.legacy_pole_entry))

    def test_only_audit_scripts_import_the_rejected_pole(self):
        importers = {
            p.name for p in sorted(SCRIPTS.glob("*.py"))
            if REJECTED_MODULE in _imported_names(p)
        }
        self.assertTrue(
            importers <= AUDIT_SCRIPTS,
            f"non-audit scripts import the rejected pole: {sorted(importers - AUDIT_SCRIPTS)}",
        )

    def test_certify_scripts_do_not_import_it(self):
        for p in sorted(SCRIPTS.glob("certify_*.py")):
            self.assertNotIn(REJECTED_MODULE, _imported_names(p), p.name)


if __name__ == "__main__":
    unittest.main()
