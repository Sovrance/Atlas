#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 -- the ``rh-formal`` gate.

    python3 scripts/ci_formal.py [--with-lean]

Runs, in order:

  1. the theorem manifest gate (``check_formal_manifest.py``) -- source hashes,
     no ``sorry``, no project-local ``axiom``, manifest self-hash, pinned
     toolchain and Mathlib commit, and, when ``lake`` is available, the Lean
     build, the statement comparator and the axiom audit;
  2. the formal certificate's own consistency (``certify_formal_boundary.py
     --check``);
  3. the boundary tests (``tests/test_formal_boundary.py``), which assert the
     things a reader would otherwise have to take on trust: that a formal
     certificate never satisfies a PSD requirement, never carries a numeric
     warrant, and never upgrades the numeric warrant of anything it backs.

The Lean layer is opportunistic by default -- set ``ATLAS_LEAN_BIN`` to a Lean
toolchain's ``bin`` directory, or pass ``--with-lean`` to make its absence a
failure rather than a skip.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _env() -> dict:
    e = dict(os.environ)
    parts = [str(ROOT / "src"), str(ROOT)]
    if e.get("PYTHONPATH"):
        parts.append(e["PYTHONPATH"])
    e["PYTHONPATH"] = ":".join(parts)
    return e


def run(label: str, argv: list) -> int:
    print(f"\n--- {label} ---", flush=True)
    return subprocess.run(argv, cwd=str(ROOT), env=_env()).returncode


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--with-lean", action="store_true")
    args = ap.parse_args()
    extra = ["--with-lean"] if args.with_lean else []

    print("=== rh-formal ===")
    steps = [
        ("theorem manifest gate",
         [sys.executable, str(ROOT / "scripts" / "check_formal_manifest.py")] + extra),
        ("formal certificate consistency",
         [sys.executable, str(ROOT / "scripts" / "certify_formal_boundary.py"), "--check"] + extra),
        ("formal boundary tests",
         [sys.executable, str(ROOT / "tests" / "test_formal_boundary.py")]),
    ]
    for label, argv in steps:
        rc = run(label, argv)
        if rc:
            print(f"\nFAIL (rh-formal): {label}", file=sys.stderr)
            return rc
    print("\nPASS (rh-formal)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
