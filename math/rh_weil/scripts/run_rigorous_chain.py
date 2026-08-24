#!/usr/bin/env python3
"""ENG-005 §12 — the full rigorous chain, in the canonical order.

    normalization
      -> E0
      -> scalar
      -> degree1
      -> compact degree2
      -> Candidate-A T=84 E3 scan
      -> T=84 point E1
      -> T=84 interior minimum (locate the minimiser)
      -> T=84 uniform E1
      -> inertia / rank-trace / moment self-tests   (ENG-006 §12)
      -> degree-3 exact identities
      -> degree-3 E3 scan
      -> degree-3 E1 positivity or inertia stratification
      -> degree-3 moments
      -> PIR
      -> clean-tree / hash validation

Any stage failing blocks everything downstream: there is no point certifying a
degree-2 determinant under a normalization that has drifted, and no point
exporting PIR facts from certificates the promotion predicate would refuse.

This supersedes ``run_rigorous_scalar.py`` as the full-chain entry point; that
script remains useful for exercising the scalar stage alone.
``run_rh_weil_suite.py`` is still the fast path and still does not re-derive any
of this.

    python3 scripts/run_rigorous_chain.py [--release] [--quick] [--skip-regenerate]

Exit codes: 0 all stages passed; 1 a stage failed; 2 a rigorous dependency is
missing.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
CERT_DIR = ROOT / "certificates"
sys.path.insert(0, str(SRC))

#: (certificate, human label) for every rigorous artifact the chain produces.
RIGOROUS_CERTS = [
    ("e1_scalar_log3_log4.json", "scalar"),
    ("e1_degree1_log3_log4.json", "degree1"),
    ("e1_degree2_compact_log3_log4.json", "compact degree2"),
    ("e1_fourier_T84_points.json", "T84 points"),
    ("e1_fourier_T84_interior_minimum.json", "T84 interior minimum"),
    ("e1_fourier_T84_uniform_degree2.json", "T84 uniform degree2"),
]

#: ENG-006 artifacts. The degree-3 E1 result lands under whichever filename the
#: certification reached, so the pair is given as alternatives: exactly one must
#: exist, and which one it is *is* the outcome (§9).
DEGREE3_E1_ALTERNATIVES = (
    "e1_degree3_odd_positivity_log3_log4.json",
    "e1_degree3_odd_inertia_log3_log4.json",
)
ENG006_CERTS = [
    ("e1_degree3_odd_moments_log3_log4.json", "degree3 moments"),
]


def degree3_e1_certificate():
    """The degree-3 E1 artifact that exists, with its filename."""
    for name in DEGREE3_E1_ALTERNATIVES:
        path = CERT_DIR / name
        if path.exists():
            return name, json.loads(path.read_text(encoding="utf-8"))
    return None, None


def _env() -> dict:
    e = dict(os.environ)
    e["PYTHONPATH"] = str(SRC) + ((":" + e["PYTHONPATH"]) if e.get("PYTHONPATH") else "")
    return e


def _run(label: str, argv: list) -> int:
    print(f"\n=== {label} ===", flush=True)
    return subprocess.run(argv, cwd=str(ROOT), env=_env()).returncode


def stage_dependencies() -> int:
    print("=== rigorous dependencies ===")
    missing = []
    for mod in ("flint", "sympy"):
        try:
            __import__(mod)
            print(f"  {mod}: present")
        except ImportError:
            missing.append(mod)
    if missing:
        print(f"  MISSING: {', '.join(missing)}", file=sys.stderr)
        print("  install: pip install -r math/rh_weil/requirements-rigorous.txt",
              file=sys.stderr)
        return 2
    return 0


def stage_normalization() -> int:
    print("\n=== normalization ===")
    import promotion

    ok, why = promotion.normalization_id_consistent()
    print(f"  active id: {promotion.active_normalization_id()}")
    print(f"  consistency: {why}")
    if not ok:
        print("  FAIL: normalization id moved (ENG-005 §15 stop condition)", file=sys.stderr)
        return 1

    artifact = json.loads((CERT_DIR / "normalization_adjudication.json").read_text(encoding="utf-8"))
    sym = artifact.get("symbolic_derivation", {})
    if not str(sym.get("symbolic_engine", "")).startswith("sympy") or not sym.get("steps"):
        print("  FAIL: degraded symbolic evidence in the adjudication artifact",
              file=sys.stderr)
        return 1
    print(f"  symbolic engine: {sym['symbolic_engine']} ({len(sym['steps'])} steps)")
    return 0


def stage_curvature_derivation() -> int:
    """§1: the Ginf'' derivation must reproduce, and must not overclaim."""
    print("\n=== curvature derivation (ENG-005 §1) ===")
    import curvature_derivation as CD

    report = CD.derivation_report()
    print(f"  status: {report['status']}  engine: {report['symbolic_engine']}")
    for st in report["steps"]:
        print(f"    [{'ok' if st['verified'] else '!!'}] {st['step']}")
    if report["status"] != "VERIFIED":
        print("  FAIL: the curvature derivation did not reproduce", file=sys.stderr)
        return 1
    if report["evidence_class"] != "E0" or not report["analytic_hypotheses"]:
        print("  FAIL: derivation must stay E0 with its interchange hypothesis named",
              file=sys.stderr)
        return 1
    print(f"  analytic hypotheses recorded: {len(report['analytic_hypotheses'])} "
          "(interchange NOT machine-verified — by design)")
    return 0


def stage_tail_lemma() -> int:
    """§3: the invalid t h_+' <= 1 assumption must still be rejected."""
    print("\n=== tail lemma (ENG-005 §3) ===")
    from interval_backend import require_flint, set_precision_bits
    import scalar_canary as SC

    _, arb, _acb, _ = require_flint()
    set_precision_bits(200)
    row = SC.invalid_assumption_is_rejected(arb)
    print(f"  t h_+'(2) in [{row['t_h_plus_prime_enclosure'][0][:12]}, "
          f"{row['t_h_plus_prime_enclosure'][1][:12]}] -> {row['verdict']}")
    if row["verdict"] != "REJECTED":
        print("  FAIL: the invalid assumption was not rejected", file=sys.stderr)
        return 1
    return 0


def _headline(cert: dict) -> str:
    """The number a certificate actually certifies, for the policy line.

    Not every rigorous certificate carries a whole-cell ``certified_lower_bound``:
    the point certificate is point-scoped by design, and the interior-minimum
    certificate bounds ``E2`` on the interval its argument governs. Falling back
    to "point-scoped" for the latter would misreport it.
    """
    if "certified_lower_bound" in cert:
        return str(cert["certified_lower_bound"])
    interior = cert.get("interior_minimum", {})
    bound = interior.get("basin_bound", {}).get("certified_lower_bound")
    if bound:
        gov = interior.get("governed_interval", ["?", "?"])
        return f"{bound} on [{gov[0]}, {gov[1]}]"
    if cert.get("point_scoped"):
        return "point-scoped"
    return "no bound field"


def stage_policy() -> int:
    print("\n=== certificate policy validation ===")
    import promotion

    bad = 0
    d3_name, _ = degree3_e1_certificate()
    to_check = list(RIGOROUS_CERTS) + list(ENG006_CERTS)
    if d3_name:
        to_check.append((d3_name, "degree3 E1"))
    for name, label in to_check:
        path = CERT_DIR / name
        if not path.exists():
            print(f"  FAIL: {name} absent", file=sys.stderr)
            bad += 1
            continue
        cert = json.loads(path.read_text(encoding="utf-8"))
        refusal = promotion.promotion_refusal(cert)
        marks = []
        if cert.get("mpmath_used"):
            marks.append("mpmath path")
        if cert.get("quick_mode"):
            marks.append("quick mode")
        if refusal or marks:
            print(f"  FAIL {label}: {refusal or ', '.join(marks)}", file=sys.stderr)
            bad += 1
        else:
            print(f"  ok {label}: {_headline(cert)}")
    return 1 if bad else 0


def stage_engines() -> int:
    """ENG-006 §12: the three new engines must pass their own self-tests."""
    print("\n=== inertia / rank-trace / moment engines (ENG-006) ===")
    for label, test in (
        ("inertia + congruence", "test_inertia_engine.py"),
        ("rank-trace theorem", "test_ranktrace.py"),
        ("spectral moments + B1", "test_moments_adapter.py"),
    ):
        if _run(label, [sys.executable, str(ROOT / "tests" / test)]):
            return 1
    return 0


def stage_degree3() -> int:
    """ENG-006 §9: report which outcome the degree-3 certification reached."""
    print("\n=== degree-3 odd block (ENG-006) ===")
    import promotion
    from inertia.certificate import satisfies_psd_requirement

    name, cert = degree3_e1_certificate()
    if cert is None:
        print(f"  FAIL: neither of {DEGREE3_E1_ALTERNATIVES} exists", file=sys.stderr)
        return 1
    inertia = cert.get("inertia", {})
    sig = tuple(inertia.get(k) for k in ("n_positive", "n_negative", "n_zero"))
    print(f"  artifact: {name}")
    print(f"  outcome : {cert.get('outcome')}  inertia {sig}")
    for k, v in cert.get("uniform_bounds", {}).items():
        print(f"    {k} >= {v['certified_lower_bound']}")
    refusal = promotion.promotion_refusal(cert)
    if refusal or cert.get("quick_mode"):
        print(f"  FAIL: {refusal or 'quick mode'}", file=sys.stderr)
        return 1
    # §11: an inertia certificate may satisfy a PSD consumer only when it
    # actually proved PSD. Report the answer rather than leaving it implied.
    print(f"  satisfies a PSD requirement: {satisfies_psd_requirement(cert)}")
    if cert.get("outcome") == "C_INCONCLUSIVE":
        print("  FAIL: degree-3 certification is inconclusive", file=sys.stderr)
        return 1
    return 0


def stage_pir() -> int:
    print("\n=== PIR ===")
    import pir_bridge

    if not pir_bridge.available():
        print("  pir package unavailable — policy already validated, export skipped")
        return 0
    payload = json.loads(Path(pir_bridge.export_pir_facts()).read_text(encoding="utf-8"))
    promoted = [f["content"]["certificate_file"] for f in payload["facts"]]
    refused = [r["certificate_file"] for r in payload["refused_promotions"]]
    print(f"  promoted: {len(promoted)}  refused: {len(refused)}")
    d3_name, _ = degree3_e1_certificate()
    expected = [n for n, _ in RIGOROUS_CERTS] + [n for n, _ in ENG006_CERTS]
    if d3_name:
        expected.append(d3_name)
    missing = [n for n in expected if n not in promoted]
    if missing:
        print(f"  FAIL: recovered certificates did not reach PIR: {missing}", file=sys.stderr)
        return 1
    for n in promoted:
        print(f"    + {n}")
    if refused:
        print(f"  still refused (expected: historical/quarantined only): {refused}")
    return 0


def stage_hashes() -> int:
    print("\n=== clean-tree / hash validation ===")
    import promotion

    bad = []
    d3_name, _ = degree3_e1_certificate()
    checked = [n for n, _ in RIGOROUS_CERTS] + [n for n, _ in ENG006_CERTS]
    if d3_name:
        checked.append(d3_name)
    for name in checked:
        cert = json.loads((CERT_DIR / name).read_text(encoding="utf-8"))
        stale = promotion.stale_dependencies(cert)
        if stale:
            bad.append((name, stale))
    if bad:
        for name, stale in bad:
            print(f"  FAIL {name}: stale {stale}", file=sys.stderr)
        return 1
    print(f"  all {len(checked)} rigorous certificates have current source hashes")

    proc = subprocess.run(["git", "status", "--porcelain", "math/rh_weil/certificates"],
                          cwd=str(ROOT.parents[1]), capture_output=True, text=True)
    dirty = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    print(f"  certificate tree: {'clean' if not dirty else str(len(dirty)) + ' modified'}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", action="store_true")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--skip-regenerate", action="store_true",
                    help="validate the committed certificates without recomputing")
    args = ap.parse_args()

    rc = stage_dependencies()
    if rc:
        return rc
    for stage in (stage_normalization, stage_curvature_derivation, stage_tail_lemma):
        if stage():
            return 1

    for label, test in (
        ("exact pole invariants", "test_pole_primitive.py"),
        ("production import scan", "test_production_imports.py"),
        ("E0 exact identities", "test_exact_identities.py"),
    ):
        if _run(label, [sys.executable, str(ROOT / "tests" / test)]):
            return 1

    if not args.skip_regenerate:
        common = (["--release"] if args.release else []) + (["--quick"] if args.quick else [])
        if _run("scalar E1", [sys.executable, str(ROOT / "scripts" / "certify_scalar_canary.py")] + common):
            return 1
        if _run("degree1 + compact degree2 E1",
                [sys.executable, str(ROOT / "scripts" / "certify_cutoff_free_e1.py")] + common):
            return 1
        if _run("T=84 scan, points, uniform",
                [sys.executable, str(ROOT / "scripts" / "certify_t84_e1.py")] + common):
            return 1
        if _run("degree-3 scan, E1, moments",
                [sys.executable, str(ROOT / "scripts" / "certify_degree3.py")] + common):
            return 1

    if _run("degree-3 exact identities (ENG-006 §7)",
            [sys.executable, str(ROOT / "tests" / "test_degree3_exact.py")]):
        return 1

    for stage in (stage_engines, stage_degree3, stage_policy, stage_pir, stage_hashes):
        if stage():
            return 1

    if _run("promotion + canary tests",
            [sys.executable, str(ROOT / "tests" / "test_promotion_and_canary.py")]):
        return 1
    if _run("ENG-005 recovery tests",
            [sys.executable, str(ROOT / "tests" / "test_eng005_recovery.py")]):
        return 1
    if _run("ENG-006 degree-3 certificate semantics",
            [sys.executable, str(ROOT / "tests" / "test_degree3_certificates.py")]):
        return 1
    if _run("ENG-006 information comparison report (§10)",
            [sys.executable, str(ROOT / "scripts" / "report_information_comparison.py")]):
        return 1

    print("\nrigorous chain: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
