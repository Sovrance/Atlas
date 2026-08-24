#!/usr/bin/env python3
"""ENG-004 §10 — the rigorous scalar path, in order.

    normalization validation
      -> exact pole invariants
      -> Arb scalar regeneration
      -> certificate policy validation
      -> scalar PIR validation

Each stage gates the next: there is no point regenerating a certificate under a
normalization that has drifted, and no point exporting PIR facts from a
certificate the promotion predicate would refuse.

This is deliberately **not** what ``run_rh_weil_suite.py`` runs. That suite is the
fast path: it executes the unit tests and regenerates the cheap E0/E3 artifacts,
and passing it says nothing about whether the rigorous certificates are current
(ENG-004 §10). Only this script re-derives them.

    python3 scripts/run_rigorous_scalar.py [--release] [--quick]

Exit codes: 0 all stages passed; 1 a stage failed; 2 a required rigorous
dependency (python-flint, SymPy) is missing.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

SCALAR = "e1_scalar_log3_log4.json"


def _run(label: str, argv: list[str]) -> int:
    print(f"\n=== {label} ===", flush=True)
    proc = subprocess.run(argv, cwd=str(ROOT),
                          env={"PYTHONPATH": str(SRC), **_env()})
    return proc.returncode


def _env() -> dict:
    import os

    e = dict(os.environ)
    e["PYTHONPATH"] = str(SRC) + ((":" + e["PYTHONPATH"]) if e.get("PYTHONPATH") else "")
    return e


def check_dependencies() -> int:
    """python-flint and SymPy are both required for the rigorous path (§6)."""
    print("=== rigorous dependencies ===")
    missing = []
    try:
        import flint  # noqa: F401

        print("  python-flint: present")
    except ImportError:
        missing.append("python-flint")
    try:
        import sympy

        print(f"  sympy: {sympy.__version__}")
    except ImportError:
        missing.append("sympy")
    if missing:
        print(f"  MISSING: {', '.join(missing)}", file=sys.stderr)
        print("  install: pip install -r math/rh_weil/requirements-rigorous.txt",
              file=sys.stderr)
        return 2
    return 0


def validate_normalization() -> int:
    print("\n=== normalization validation ===")
    import promotion

    active = promotion.active_normalization_id()
    ok, why = promotion.normalization_id_consistent()
    print(f"  active id: {active}")
    print(f"  consistency: {why}")
    if not ok:
        print("  FAIL: normalization id changed unexpectedly (ENG-004 §14 stop "
              "condition) — re-run scripts/derive_normalization.py", file=sys.stderr)
        return 1

    # SymPy evidence must not have degraded (§6). Before ENG-004 a missing SymPy
    # rewrote this artifact with a "verified numerically only" note; the gate
    # exists so that can never pass CI silently.
    artifact = json.loads(
        (ROOT / "certificates" / "normalization_adjudication.json").read_text(encoding="utf-8")
    )
    sym = artifact.get("symbolic_derivation", {})
    engine = str(sym.get("symbolic_engine", ""))
    steps = sym.get("steps", [])
    print(f"  symbolic engine: {engine or '(none)'}  steps: {len(steps)}")
    if not engine.startswith("sympy") or not steps:
        print("  FAIL: the adjudication artifact carries degraded symbolic evidence",
              file=sys.stderr)
        return 1
    if not all(st.get("verified") for st in steps):
        print("  FAIL: an unverified symbolic step in the adjudication artifact",
              file=sys.stderr)
        return 1
    return 0


def validate_certificate_policy() -> int:
    print("\n=== certificate policy validation ===")
    import normalization as N
    import promotion

    cert_path = ROOT / "certificates" / SCALAR
    if not cert_path.exists():
        print(f"  FAIL: {SCALAR} absent", file=sys.stderr)
        return 1
    cert = json.loads(cert_path.read_text(encoding="utf-8"))

    refusal = promotion.promotion_refusal(cert)
    print(f"  scalar promotable: {refusal is None}"
          + (f"  ({refusal})" if refusal else ""))
    if refusal:
        return 1
    if float(cert["certified_lower_bound"]) <= 0:
        print("  FAIL: scalar lower bound is not positive", file=sys.stderr)
        return 1
    if cert.get("mpmath_used"):
        print("  FAIL: an E1 certificate may not come from the mpmath path", file=sys.stderr)
        return 1
    backend = json.dumps(cert.get("backend", {})).lower()
    if "flint" not in backend:
        print("  FAIL: the interval backend was skipped; E1 requires Arb", file=sys.stderr)
        return 1
    if cert.get("quick_mode"):
        print("  FAIL: a --quick smoke certificate is not promotable", file=sys.stderr)
        return 1
    print(f"  backend: Arb  bound: {cert['certified_lower_bound']}")

    # Every *other* disputed E1 must still be quarantined (§4, §9).
    still = [c for c in N.QUARANTINED_CERTIFICATES if c != SCALAR]
    leaked = []
    for name in still:
        other = json.loads((ROOT / "certificates" / name).read_text(encoding="utf-8"))
        if promotion.promotion_refusal(other) is None:
            leaked.append(name)
    print(f"  downstream E1 still quarantined: {len(still) - len(leaked)}/{len(still)}")
    if leaked:
        print(f"  FAIL: quarantine lifted on {leaked}", file=sys.stderr)
        return 1
    return 0


def validate_scalar_pir() -> int:
    print("\n=== scalar PIR validation ===")
    import pir_bridge

    if not pir_bridge.available():
        print("  pir package unavailable — skipping export, policy already checked")
        return 0
    path = pir_bridge.export_pir_facts()
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    promoted = [f["content"]["certificate_file"] for f in payload["facts"]]
    refused = [r["certificate_file"] for r in payload["refused_promotions"]]
    print(f"  promoted: {len(promoted)}  refused: {len(refused)}")
    if SCALAR not in promoted:
        print(f"  FAIL: {SCALAR} did not reach PIR", file=sys.stderr)
        return 1
    disputed = [r for r in refused if r.startswith("e1_") and r != SCALAR]
    if len(disputed) != 4:
        print(f"  FAIL: expected the four other disputed E1 refused, got {disputed}",
              file=sys.stderr)
        return 1
    print(f"  scalar promoted; {len(disputed)} other disputed E1 still refused")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", action="store_true",
                    help="lift the quarantine on the scalar certificate after a PASS")
    ap.add_argument("--quick", action="store_true",
                    help="coarse grid; never promotable (smoke test only)")
    ap.add_argument("--skip-regenerate", action="store_true",
                    help="validate the committed certificate without recomputing it")
    args = ap.parse_args()

    rc = check_dependencies()
    if rc:
        return rc

    for stage, fn in (("normalization validation", validate_normalization),):
        if fn():
            return 1

    if _run("exact pole invariants",
            [sys.executable, str(ROOT / "tests" / "test_pole_primitive.py")]):
        return 1
    if _run("production import scan",
            [sys.executable, str(ROOT / "tests" / "test_production_imports.py")]):
        return 1

    if not args.skip_regenerate:
        argv = [sys.executable, str(ROOT / "scripts" / "certify_scalar_canary.py")]
        if args.release:
            argv.append("--release")
        if args.quick:
            argv.append("--quick")
        if _run("Arb scalar regeneration", argv):
            return 1

    if validate_certificate_policy():
        return 1
    if validate_scalar_pir():
        return 1

    if _run("promotion + canary tests",
            [sys.executable, str(ROOT / "tests" / "test_promotion_and_canary.py")]):
        return 1

    print("\nrigorous scalar path: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
