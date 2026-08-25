#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §12 (WO-RH-44) -- emit the FORMAL_THEOREM_CERTIFICATE.

Runs the manifest gate first, then writes
``certificates/formal_theorem_certificate.json`` from the manifest.

The certificate carries no numeric bound. It records which finite theorems are
proved, under which pinned toolchain, with which axioms, and which numeric
certificates they back. The numeric warrants of those certificates are
unchanged by it -- that is the whole point of keeping the two warrants in
separate fields.

Usage:

    python3 scripts/certify_formal_boundary.py            # gate offline, then write
    python3 scripts/certify_formal_boundary.py --with-lean  # require the Lean layer
    python3 scripts/certify_formal_boundary.py --check      # verify, do not write
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for extra in (ROOT, ROOT / "src"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))

import formal_evidence as F  # noqa: E402
from certificate_io import write_certificate  # noqa: E402

CERT_NAME = "formal_theorem_certificate.json"


def run_gate(with_lean: bool) -> int:
    cmd = [sys.executable, str(ROOT / "scripts" / "check_formal_manifest.py")]
    if with_lean:
        cmd.append("--with-lean")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    return proc.returncode


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--with-lean", action="store_true",
                    help="fail unless the Lean layer of the manifest gate can run")
    ap.add_argument("--check", action="store_true",
                    help="verify the committed certificate instead of rewriting it")
    args = ap.parse_args()

    rc = run_gate(args.with_lean)
    if rc:
        print("certify_formal_boundary: manifest gate failed", file=sys.stderr)
        return rc

    if args.check:
        path = ROOT / "certificates" / CERT_NAME
        if not path.exists():
            print(f"certify_formal_boundary: missing {CERT_NAME}", file=sys.stderr)
            return 1
        cert = json.loads(path.read_text(encoding="utf-8"))
        problems = F.formal_certificate_problems(cert)
        if problems:
            print("certify_formal_boundary: FAIL", file=sys.stderr)
            for p in problems:
                print(f"  - {p}", file=sys.stderr)
            return 1
        print(f"  {CERT_NAME}: {len(cert.get('formal_theorem_ids') or ())} theorems, "
              f"manifest {cert.get('formal_manifest_id')}")
        return 0

    body = F.build_formal_certificate()
    problems = F.formal_certificate_problems(body)
    if problems:
        print("certify_formal_boundary: refusing to write an inconsistent certificate",
              file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    path = write_certificate(CERT_NAME, body)
    print(f"wrote {path}")
    print(f"  manifest: {body['formal_manifest_id']}")
    print(f"  theorems: {len(body['formal_theorem_ids'])}  "
          f"axioms: {', '.join(body['axioms'])}")
    print(f"  recorded unproved: {len(body['unproved_statements'])}")
    print("  numeric warrant: none (this artifact grades an implication, not a number)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
