#!/usr/bin/env python3
"""ENG-006 §13 — the two inertia CI gates.

This repository has no workflow files; its CI is script-shaped
(``ci/run_all_certified.py``, ``scripts/run_rigorous_chain.py``), so the two
gates §13 names are exposed the same way rather than as YAML that nothing would
run.

    python3 scripts/ci_inertia.py --gate fast       # no python-flint needed
    python3 scripts/ci_inertia.py --gate rigorous   # python-flint required

**fast** covers what is exact and cheap: congruence/Sylvester regression, the E0
degree-3 kernel identities, schema validation, and PIR content-kind validation.

**rigorous** requires python-flint and covers the interval engines, the degree-3
E1 result, interval moment extraction, the rank-trace runtime checks, and
dependency freshness.

The last check in the rigorous gate is the one worth stating plainly: §13 says
"no floating eigenvalue computation may satisfy a rigorous gate", so the gate
greps the rigorous modules for a floating eigenvalue solver rather than trusting
that none crept in. The odd degree-3 block is 2x2 and its spectrum has a closed
form, so there is no reason for a solver to appear -- which is exactly why an
appearance would be worth catching.

Exit codes: 0 pass, 1 a check failed, 2 a required dependency is missing.
"""
from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

CERT_DIR = ROOT / "certificates"
SCHEMA_DIR = ROOT / "inertia" / "schemas"

#: Modules that may support a rigorous claim. None of them may reach a floating
#: eigenvalue solver.
RIGOROUS_MODULES = (
    "src/degree3.py", "src/pole.py", "src/core.py", "src/weil_entries.py",
    "src/archimedean_realspace.py", "src/interval_cover.py",
    "src/interval_backend.py",
    "inertia/ldl.py", "inertia/stratify.py", "inertia/congruence.py",
    "inertia/certificate.py",
    "moments/spectral_moments.py", "moments/feasible_spectrum.py",
    "moments/adapter.py", "ranktrace/theorem.py",
)

#: Names that would mean a floating eigenvalue solver is in the path.
FORBIDDEN_IN_RIGOROUS = ("numpy", "scipy", "eigvalsh", "eigvals", "eigh", "eig")


def _run(label: str, argv: list) -> int:
    print(f"\n--- {label} ---", flush=True)
    return subprocess.run(argv, cwd=str(ROOT)).returncode


def check_schemas() -> int:
    print("\n--- certificate schemas ---")
    from inertia.certificate import (
        build_inertia_certificate,
        validate_against_schema,
    )
    from inertia.ldl import exact_inertia

    schema = json.loads((SCHEMA_DIR / "inertia_certificate.schema.json")
                        .read_text(encoding="utf-8"))
    cert = build_inertia_certificate(
        exact_inertia([[1, 0], [0, -1]]), dimension=2, program="ci",
        work_order="ENG-006", evidence_class="E0",
        normalization_certificate_id="norm_ci")
    errs = validate_against_schema(cert, schema)
    if errs:
        print(f"  FAIL: {errs}", file=sys.stderr)
        return 1
    # The validator must also reject: a permissive validator passes everything.
    bad = dict(cert, rh_proof_claim=True)
    if not validate_against_schema(bad, schema):
        print("  FAIL: validator accepted rh_proof_claim=true", file=sys.stderr)
        return 1
    for name in ("inertia_certificate.schema.json",
                 "rank_trace_certificate.schema.json",
                 "spectral_moment_certificate.schema.json"):
        json.loads((SCHEMA_DIR / name).read_text(encoding="utf-8"))
    print("  ok: 3 schemas parse; a valid body validates and an invalid one is rejected")
    return 0


def check_pir_kinds() -> int:
    print("\n--- PIR content kinds ---")
    import pir_bridge
    from inertia.certificate import satisfies_psd_requirement

    expected = {"WEIL_INERTIA_CERTIFICATE", "WEIL_INERTIA_STRATIFICATION",
                "WEIL_RANK_TRACE_CERTIFICATE", "WEIL_SPECTRAL_MOMENT_CERTIFICATE"}
    if set(pir_bridge.CONTENT_KINDS) != expected:
        print(f"  FAIL: content kinds are {pir_bridge.CONTENT_KINDS}", file=sys.stderr)
        return 1
    # §11: an inertia certificate must never satisfy a PSD-requiring consumer.
    indefinite = {"content_kind": "WEIL_INERTIA_CERTIFICATE", "status": "PASS",
                  "evidence_class": "E1", "n_negative": 1, "n_zero": 0,
                  "rh_proof_claim": False}
    strat_psd = {"content_kind": "WEIL_INERTIA_STRATIFICATION", "status": "PASS",
                 "evidence_class": "E1", "n_negative": 0, "n_zero": 0,
                 "rh_proof_claim": False}
    if satisfies_psd_requirement(indefinite):
        print("  FAIL: an indefinite inertia certificate satisfied a PSD consumer",
              file=sys.stderr)
        return 1
    if satisfies_psd_requirement(strat_psd):
        print("  FAIL: a stratification satisfied a PSD consumer", file=sys.stderr)
        return 1
    print(f"  ok: {len(expected)} content kinds; PSD gate refuses inertia and "
          "stratification certificates")
    return 0


def check_no_floating_eigensolver() -> int:
    print("\n--- no floating eigenvalue solver in a rigorous path ---")
    bad = []
    for rel in RIGOROUS_MODULES:
        path = ROOT / rel
        if not path.exists():
            print(f"  FAIL: {rel} is missing", file=sys.stderr)
            return 1
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    if a.name.split(".")[0] in FORBIDDEN_IN_RIGOROUS:
                        bad.append(f"{rel}: import {a.name}")
            elif isinstance(node, ast.ImportFrom) and node.module:
                if node.module.split(".")[0] in FORBIDDEN_IN_RIGOROUS:
                    bad.append(f"{rel}: from {node.module} import ...")
            elif isinstance(node, ast.Attribute) and node.attr in FORBIDDEN_IN_RIGOROUS:
                bad.append(f"{rel}: attribute .{node.attr}")
    if bad:
        for b in bad:
            print(f"  FAIL {b}", file=sys.stderr)
        return 1
    print(f"  ok: {len(RIGOROUS_MODULES)} rigorous modules, no eigenvalue solver reachable")
    return 0


def check_dependencies_current() -> int:
    print("\n--- degree-3 certificate freshness ---")
    import promotion

    names = ["e1_degree3_odd_moments_log3_log4.json"]
    for cand in ("e1_degree3_odd_positivity_log3_log4.json",
                 "e1_degree3_odd_inertia_log3_log4.json"):
        if (CERT_DIR / cand).exists():
            names.append(cand)
    if len(names) == 1:
        print("  FAIL: no degree-3 E1 certificate exists", file=sys.stderr)
        return 1
    rc = 0
    for name in names:
        cert = json.loads((CERT_DIR / name).read_text(encoding="utf-8"))
        stale = promotion.stale_dependencies(cert)
        refusal = promotion.promotion_refusal(cert)
        if stale or refusal:
            print(f"  FAIL {name}: {refusal or 'stale ' + str(stale)}", file=sys.stderr)
            rc = 1
        else:
            print(f"  ok {name}")
    return rc


def gate_fast() -> int:
    print("=== rh-inertia-fast ===")
    rc = 0
    for label, test in (
        ("exact congruence / Sylvester", "test_inertia_engine.py"),
        ("exact E0 degree-3 kernels", "test_degree3_exact.py"),
        ("rank-trace theorem runtime", "test_ranktrace.py"),
    ):
        rc = rc or _run(label, [sys.executable, str(ROOT / "tests" / test)])
    rc = rc or check_schemas()
    rc = rc or check_pir_kinds()
    return rc


def gate_rigorous() -> int:
    print("=== rh-inertia-rigorous ===")
    try:
        import flint  # noqa: F401
    except ImportError:
        print("ERROR: python-flint is required for the rigorous gate", file=sys.stderr)
        return 2
    rc = 0
    for label, test in (
        ("interval LDL / inertia", "test_inertia_engine.py"),
        ("spectral moments + B1 adapter", "test_moments_adapter.py"),
        ("rank-trace theorem runtime", "test_ranktrace.py"),
        ("degree-3 exact identities", "test_degree3_exact.py"),
    ):
        rc = rc or _run(label, [sys.executable, str(ROOT / "tests" / test)])
    rc = rc or _run("degree-3 E1 certificate checks",
                    [sys.executable, str(ROOT / "tests" / "test_degree3_certificates.py")])
    rc = rc or check_no_floating_eigensolver()
    rc = rc or check_dependencies_current()
    return rc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", choices=["fast", "rigorous", "both"], default="both")
    args = ap.parse_args()
    rc = 0
    if args.gate in ("fast", "both"):
        rc = rc or gate_fast()
    if args.gate in ("rigorous", "both"):
        rc = rc or gate_rigorous()
    print("\nPASS" if rc == 0 else "\nFAIL", f"(rh-inertia-{args.gate})")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
