#!/usr/bin/env python3
"""Regenerate the Connes-CvS external cross-validation certificate.

Optional research runner. Requires:
  pip install 'connes-cvs==0.3.1' python-flint mpmath
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "external"))
sys.path.insert(0, str(ROOT / "src"))

import crosschecks  # noqa: E402


def main() -> int:
    results = crosschecks.run_acceptance_suite()
    path = crosschecks.write_certificate(results=results)
    gate = [r for r in results if r.test_id in {"XC-01", "XC-02", "XC-03"}]
    print(f"wrote {path}")
    for r in results:
        print(f"  {r.test_id} dps={r.precision_dps}: {r.status} ({r.evidence_class})")
    if not all(r.status == "pass" for r in gate):
        print("ACCEPTANCE GATE FAILED", file=sys.stderr)
        return 1
    print("ACCEPTANCE GATE PASSED (XC-01..03)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
