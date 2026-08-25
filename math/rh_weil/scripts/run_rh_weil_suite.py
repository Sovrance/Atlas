#!/usr/bin/env python3
"""Dedicated RH/Weil runner (WO-RH-07) — the **fast** path.

Does NOT expand ci/run_all_certified.py. Runs the unit tests and regenerates the
cheap E0/E3 artifacts. Never promotes imported notebook claims.

ENG-004 §10: passing this suite does **not** mean the rigorous certificates are
current. It does not re-derive the scalar E1 canary (minutes of Arb quadrature)
and it does not check the promotion policy end to end. Use
``scripts/run_rigorous_scalar.py`` for that, and read its exit code -- not this
one -- before believing an E1 claim is fresh.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
TESTS = ROOT / "tests"
sys.path.insert(0, str(SRC))


def run_unittest_file(path: Path) -> int:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC) + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run(
        [sys.executable, str(path)],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
    )
    sys.stdout.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    return proc.returncode


def regenerate_certificates() -> None:
    import certificate_io
    import fourier
    import mpmath as mp
    import scalar

    e0 = certificate_io.build_e0_exact_certificate()
    certificate_io.write_certificate("e0_exact_identities.json", e0)

    report = scalar.verify_scalar_cell().to_dict()
    certificate_io.write_certificate(
        "e0_scalar_cell_log3_log4.json",
        certificate_io.build_e0_scalar_cell_certificate(report),
    )

    Ls = [mp.log(3), mp.mpf("1.1059498113"), mp.mpf("1.20"), mp.log(4)]
    scan = {
        "T": 84,
        "rows": fourier.scan_E2_probe(Ls, T=84, dps=25),
        "forms_implemented": ["H0", "Hb", "H0_L_jet", "Hb_L_jet"],
    }
    # ENG-005 §6 owns ``e3_fourier_T84_scan.json``: that file is the fresh
    # Candidate-A topology scan, built from exact jets under the adopted pole and
    # read by the T=84 certifier to choose its starting bracket. This probe is the
    # legacy mpmath one -- lower evidence, different shape, and Candidate-A only by
    # accident of the L values it happens to sample. Writing it to that name would
    # silently replace the scan the E1 stages depend on with weaker content, on
    # every fast-path run. It gets its own name.
    certificate_io.write_certificate(
        "e3_fourier_T84_probe_mpmath.json",
        certificate_io.build_e3_fourier_scan_certificate(scan),
    )
    certificate_io.write_certificate(
        "work_order_status.json",
        certificate_io.build_work_order_status(),
    )

    # WO-RH-17: the adjudication certificate is regenerated from the derivation,
    # and the quarantine is re-asserted after work_order_status.json is rewritten.
    import subprocess as _sp

    for script in ("derive_normalization.py", "quarantine_normalization.py"):
        rc = _sp.run([sys.executable, str(ROOT / "scripts" / script)],
                     cwd=str(ROOT), capture_output=True, text=True)
        if rc.returncode != 0:  # pragma: no cover
            raise RuntimeError(f"{script} failed: {rc.stderr[-400:]}")

    # Imported notebook state must remain pending — never silently rewrite to E1.
    imported = ROOT / "certificates" / "imported_notebook_state.json"
    data = json.loads(imported.read_text(encoding="utf-8"))
    assert data.get("status") == "IMPORTED_PENDING_REGENERATION"
    assert data.get("hard_constraints_certified") is False


def main() -> int:
    test_files = [
        TESTS / "test_exact_identities.py",
        TESTS / "test_scalar_verifier.py",
        TESTS / "test_fourier_forms.py",
        TESTS / "test_connes_cvs_adapter_contract.py",
        TESTS / "test_connes_cvs_crosschecks.py",
        TESTS / "test_eng002_parity.py",
        TESTS / "test_normalization_adjudication.py",
        TESTS / "test_pole_primitive.py",
        TESTS / "test_production_imports.py",
        TESTS / "test_promotion_and_canary.py",
        TESTS / "test_eng005_recovery.py",
    ]
    failed = 0
    for tf in test_files:
        if not tf.exists():
            continue
        print(f"=== {tf.name} ===")
        rc = run_unittest_file(tf)
        if rc != 0:
            failed += 1
    print("=== regenerate certificates ===")
    try:
        regenerate_certificates()
        print("certificates regenerated under math/rh_weil/certificates/")
    except Exception as exc:
        print(f"certificate regeneration failed: {exc}", file=sys.stderr)
        failed += 1
    # ENG-007 §3.3: the documentation gate runs in the fast gate. It is cheap,
    # it needs no rigorous dependency, and a stale instruction is the one defect
    # class that propagates into *future* mathematics rather than staying put.
    print("=== rh-docs ===")
    import subprocess as _sp

    rc = _sp.run([sys.executable, str(ROOT / "scripts" / "check_docs.py")],
                 cwd=str(ROOT))
    if rc.returncode != 0:
        failed += 1
    if run_unittest_file(TESTS / "test_docs_freshness.py") != 0:
        failed += 1

    print("NOTE (ENG-004 §10 / ENG-005 §12): this fast suite does not re-derive any "
          "rigorous certificate — not the scalar canary, not degree-1/degree-2, not "
          "the T=84 chain. Run scripts/run_rigorous_chain.py for that and read its "
          "exit code before believing an E1 claim is current.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
