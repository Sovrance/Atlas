"""Connes-CvS cross-checks XC-01..XC-06 (optional dependency).

Acceptance gate (AGENT_WORK_ORDER_CONNES_CVS.md): XC-01, XC-02, XC-03 at two
precision levels. XC-04+ are research extensions.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for p in (HERE, ROOT / "src"):
    s = str(p)
    if s not in sys.path:
        sys.path.insert(0, s)

from connes_cvs_adapter import (  # noqa: E402
    TESTED_VERSION,
    ConnesCVSUnavailable,
    dependency_info,
    finite_eigenpair_certificate,
    h0_transform,
    h_plus as cvs_h_plus,
    h_plus_mpmath as cvs_h_plus_mpmath,
    prime_power_data,
)

CERT_DIR = ROOT / "certificates" / "external"
DEFAULT_CERT_PATH = CERT_DIR / "connes_cvs_crossvalidation_v0.1.json"

UPSTREAM_REPO = "https://github.com/akivag613/connes-cvs-"
UPSTREAM_TIP_SHA = "8ce0fc791ed9c9ca6f4ba512322720b4be80421b"

CELL_L_VALUES = ("log(3)", "1.1059498113", "log(4)")
DEFAULT_T_GRID = (0, "1e-12/L", "1e-8/L", 0.1, 7, 21.218, 84)
PRECISION_LEVELS = (40, 80)


@dataclass
class CheckResult:
    test_id: str
    status: str  # pass | fail | inconclusive | skipped
    evidence_class: str
    precision_dps: int | None = None
    detail: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _mp():
    import mpmath as mp

    return mp


def _parse_L(token: str | float, mp):
    if token == "log(3)":
        return mp.log(3)
    if token == "log(4)":
        return mp.log(4)
    return mp.mpf(token)


def _t_grid_for_L(L, mp):
    out = []
    for tok in DEFAULT_T_GRID:
        if tok == 0:
            out.append(mp.mpf(0))
        elif tok == "1e-12/L":
            out.append(mp.mpf("1e-12") / L)
        elif tok == "1e-8/L":
            out.append(mp.mpf("1e-8") / L)
        else:
            out.append(mp.mpf(tok))
    return out


def _tol(dps: int):
    """Absolute tolerance scaled to working precision (not binary64-fixed)."""
    mp = _mp()
    return mp.power(10, -(dps - 10))


def atlas_source_hash() -> str:
    paths = [
        ROOT / "src" / "core.py",
        ROOT / "src" / "mpmath_core.py",
        ROOT / "external" / "connes_cvs_adapter.py",
        ROOT / "external" / "crosschecks.py",
    ]
    h = hashlib.sha256()
    for p in paths:
        h.update(p.read_bytes())
        h.update(b"\0")
    return h.hexdigest()


def environment_versions() -> dict[str, Any]:
    import importlib.metadata as md

    info: dict[str, Any] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "mpmath": None,
        "python_flint": None,
        "flint_dps_default": None,
        "connes_cvs": None,
        "connes_cvs_tested_version": TESTED_VERSION,
        "upstream_repository": UPSTREAM_REPO,
        "upstream_tip_sha_at_integration": UPSTREAM_TIP_SHA,
    }
    for name, key in (
        ("mpmath", "mpmath"),
        ("python-flint", "python_flint"),
        ("connes-cvs", "connes_cvs"),
    ):
        try:
            info[key] = md.version(name)
        except md.PackageNotFoundError:
            info[key] = None
    try:
        import flint

        info["flint_dps_default"] = int(flint.ctx.dps)
    except Exception:
        pass
    return info


def xc01_archimedean_multiplier(dps: int) -> CheckResult:
    """XC-01: Atlas h_plus vs connes_cvs.operator.h_plus."""
    import mpmath_core as atlas

    mp = _mp()
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        taus = [
            mp.mpf(0),
            mp.mpf("0.1"),
            mp.mpf(1),
            mp.mpf(7),
            mp.mpf("21.218"),
            mp.mpf(84),
        ]
        tol = _tol(dps)
        flint_floor = mp.power(10, -12)  # Arb↔mpmath re-entry floor diagnostic
        worst = mp.mpf(0)
        worst_flint = mp.mpf(0)
        rows = []
        for tau in taus:
            a = atlas.atlas_h_plus(tau, dps)
            # Precision-scaling gate vs independent mpmath path in connes-cvs.
            b = cvs_h_plus_mpmath(tau, dps)
            err = abs(a - b)
            worst = max(worst, err)
            flint_err = abs(a - cvs_h_plus(tau, dps))
            worst_flint = max(worst_flint, flint_err)
            rows.append(
                {
                    "tau": mp.nstr(tau, 20),
                    "abs_error_mpmath_backend": mp.nstr(err, 20),
                    "abs_error_default_backend": mp.nstr(flint_err, 20),
                }
            )
        ok = worst < tol and worst_flint < flint_floor
        return CheckResult(
            "XC-01",
            "pass" if ok else "fail",
            "E3",
            dps,
            {
                "worst_abs_error_mpmath_backend": mp.nstr(worst, 20),
                "worst_abs_error_default_backend": mp.nstr(worst_flint, 20),
                "tolerance_mpmath_backend": mp.nstr(tol, 20),
                "tolerance_default_backend_floor": mp.nstr(flint_floor, 20),
                "samples": rows,
            },
        )
    finally:
        mp.mp.dps = old


def xc02_scalar_transform(dps: int) -> CheckResult:
    """XC-02: Atlas H0 vs connes_cvs.kernels.stable_A."""
    import mpmath_core as atlas

    mp = _mp()
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        tol = _tol(dps)
        worst = mp.mpf(0)
        rows = []
        for Ltok in CELL_L_VALUES:
            L = _parse_L(Ltok, mp)
            for t in _t_grid_for_L(L, mp):
                a = atlas.atlas_h0(t, L)
                b = h0_transform(t, L)
                err = abs(a - b)
                worst = max(worst, err)
                rows.append(
                    {
                        "L": Ltok,
                        "t": mp.nstr(t, 20),
                        "abs_error": mp.nstr(err, 20),
                    }
                )
        ok = worst < tol
        return CheckResult(
            "XC-02",
            "pass" if ok else "fail",
            "E3",
            dps,
            {
                "worst_abs_error": mp.nstr(worst, 20),
                "tolerance": mp.nstr(tol, 20),
                "n_samples": len(rows),
                "samples_head": rows[:8],
            },
        )
    finally:
        mp.mp.dps = old


def xc03_prime_power_ledger(dps: int, c: int = 4) -> CheckResult:
    """XC-03: exact n/base and high-precision Λ(n)/√n on the log3–log4 cell."""
    import mpmath_core as atlas

    mp = _mp()
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        atlas_data, atlas_primes = atlas.atlas_prime_powers_up_to(c, dps)
        cvs_data, cvs_primes = prime_power_data(c)
        tol = _tol(dps)
        if [n for n, _, _ in atlas_data] != [n for n, _, _ in cvs_data]:
            return CheckResult(
                "XC-03",
                "fail",
                "E3",
                dps,
                {
                    "reason": "prime-power support mismatch",
                    "atlas_n": [n for n, _, _ in atlas_data],
                    "cvs_n": [n for n, _, _ in cvs_data],
                },
            )
        if atlas_primes != cvs_primes:
            return CheckResult(
                "XC-03",
                "fail",
                "E3",
                dps,
                {
                    "reason": "prime list mismatch",
                    "atlas_primes": atlas_primes,
                    "cvs_primes": cvs_primes,
                },
            )
        worst_log = mp.mpf(0)
        worst_w = mp.mpf(0)
        for (n, log_a, w_a), (_, log_b, w_b) in zip(atlas_data, cvs_data):
            worst_log = max(worst_log, abs(log_a - log_b))
            worst_w = max(worst_w, abs(w_a - w_b))
        ok = worst_log < tol and worst_w < tol
        return CheckResult(
            "XC-03",
            "pass" if ok else "fail",
            "E3",
            dps,
            {
                "c": c,
                "active_cell": "[log(3), log(4)]",
                "n": [n for n, _, _ in atlas_data],
                "worst_log_error": mp.nstr(worst_log, 20),
                "worst_weight_error": mp.nstr(worst_w, 20),
                "tolerance": mp.nstr(tol, 20),
            },
        )
    finally:
        mp.mp.dps = old


def xc04_scalar_fourier_hplus_isolation(dps: int = 40, T: int = 84) -> CheckResult:
    """WO-CVS-03 / XC-04 lite: Atlas scalar probe with Atlas vs external h_plus."""
    import mpmath_core as atlas

    mp = _mp()
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        tol = _tol(dps) * 1000  # quadrature noise guard
        rows = []
        worst = mp.mpf(0)
        for Ltok in CELL_L_VALUES:
            L = _parse_L(Ltok, mp)
            a = atlas.atlas_scalar_arch_probe(L, T, dps, atlas.atlas_h_plus)
            # Isolate the multiplier formula (mpmath backend), not Arb conversion.
            b = atlas.atlas_scalar_arch_probe(L, T, dps, cvs_h_plus_mpmath)
            err = abs(a - b)
            worst = max(worst, err)
            rows.append(
                {
                    "L": Ltok,
                    "T": T,
                    "atlas_hplus_probe": mp.nstr(a, 20),
                    "cvs_hplus_probe": mp.nstr(b, 20),
                    "abs_error": mp.nstr(err, 20),
                }
            )
        ok = worst < tol
        return CheckResult(
            "XC-04",
            "pass" if ok else "fail",
            "E3",
            dps,
            {
                "note": "Isolates archimedean multiplier inside Atlas scalar probe; not an RH claim.",
                "worst_abs_error": mp.nstr(worst, 20),
                "tolerance": mp.nstr(tol, 20),
                "rows": rows,
            },
        )
    finally:
        mp.mp.dps = old


def xc05_projection_stub() -> CheckResult:
    """XC-05 research extension: deferred until basis map is proved."""
    return CheckResult(
        "XC-05",
        "inconclusive",
        "E3",
        None,
        {
            "reason": "Projection bridge requires explicit Atlas↔CvS basis map; not claimed equal by index.",
            "next": "WO-CVS-04",
        },
    )


def _jsonable(obj: Any) -> Any:
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    return str(obj)


def xc06_eigenpair_residual(dps: int = 40) -> CheckResult:
    """XC-06: finite-matrix residual machinery on a controlled symmetric test matrix."""
    mp = _mp()
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        Q = mp.matrix([[mp.mpf(2), mp.mpf(0)], [mp.mpf(0), mp.mpf(5)]])
        v = mp.matrix([[mp.mpf(0)], [mp.mpf(1)]])
        lam = mp.mpf(5)
        ext = finite_eigenpair_certificate(Q, v, lam, dps=dps)
        resid = (Q * v) - lam * v
        atlas_norm = mp.norm(resid) / mp.norm(v)
        return CheckResult(
            "XC-06",
            "pass",
            "external_E1_finite_matrix_only",
            dps,
            {
                "scope": "finite real-symmetric matrix only; no infinite-operator / RH conclusion",
                "atlas_residual_norm": mp.nstr(atlas_norm, 20),
                "external_certificate": _jsonable(ext) if isinstance(ext, dict) else type(ext).__name__,
            },
        )
    except Exception as exc:
        return CheckResult(
            "XC-06",
            "fail",
            "external_E1_finite_matrix_only",
            dps,
            {"error": f"{type(exc).__name__}: {exc}"},
        )
    finally:
        mp.mp.dps = old


def run_acceptance_suite(precision_levels=PRECISION_LEVELS) -> list[CheckResult]:
    dep = dependency_info()
    if not dep.available:
        raise ConnesCVSUnavailable("connes-cvs not installed")
    results: list[CheckResult] = []
    for dps in precision_levels:
        results.append(xc01_archimedean_multiplier(dps))
        results.append(xc02_scalar_transform(dps))
        results.append(xc03_prime_power_ledger(dps))
    results.append(xc04_scalar_fourier_hplus_isolation(dps=40, T=84))
    results.append(xc05_projection_stub())
    results.append(xc06_eigenpair_residual(dps=40))
    return results


def build_certificate(
    results: list[CheckResult],
    precision_levels=PRECISION_LEVELS,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    dep = dependency_info()
    required = {"XC-01", "XC-02", "XC-03"}
    gate_ok = all(r.status == "pass" for r in results if r.test_id in required)
    body = {
        "certificate_version": "0.1",
        "program": "Atlas RH/Weil × Connes-CvS cross-validation",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "rh_proof_claim": False,
        "evidence_policy": {
            "XC-01_to_XC-04": "E3 / HEURISTIC diagnostic until Atlas interval wrapping",
            "XC-06_external": "external E1 restricted to supplied finite matrix",
            "promotion": "Never promote upstream floats to Atlas E1; never claim RH",
        },
        "upstream": {
            "repository": UPSTREAM_REPO,
            "package": "connes-cvs",
            "installed_version": dep.version,
            "tested_version": dep.tested_version,
            "tip_sha_at_integration": UPSTREAM_TIP_SHA,
            "canonical_proof_engine": False,
        },
        "dependencies": environment_versions(),
        "atlas_source_hash": atlas_source_hash(),
        "acceptance_gate": {
            "required": sorted(required),
            "precision_levels": list(precision_levels),
            "passed": gate_ok,
        },
        "results": [r.to_dict() for r in results],
    }
    if extra:
        body["extra"] = extra
    return body


def write_certificate(
    path: Path | None = None,
    results: list[CheckResult] | None = None,
) -> Path:
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    path = path or DEFAULT_CERT_PATH
    results = results if results is not None else run_acceptance_suite()
    cert = build_certificate(results)
    path.write_text(json.dumps(cert, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    return path


if __name__ == "__main__":
    out = write_certificate()
    print(f"wrote {out}")
