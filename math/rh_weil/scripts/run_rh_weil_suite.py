#!/usr/bin/env python3
"""Dedicated RH/Weil runner (WO-RH-07).

Does NOT expand ci/run_all_certified.py. Regenerates E0 certificates and an
E3 Fourier probe scan. Never promotes imported notebook claims.
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
    certificate_io.write_certificate(
        "e3_fourier_T84_scan.json",
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
        TESTS / "test_normalization_adjudication.py",
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
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
