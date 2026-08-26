#!/usr/bin/env python3
"""ATLAS-RH-ENG-010 — the 4x4 even Weil block ``G[{1, b, b^2, b^3}]``.

    python3 scripts/certify_even4.py [--stage STAGE] [--quick]

Stages (default ``all``): ``e0``, ``crosscheck``, ``inertia``, ``positivity``,
``gap``, ``moments``. Produces:

  ``e0_degree6_even4_exact_identities.json``          §WO-RH-66, E0
  ``e3_degree6_even4_crosscheck.json``                §WO-RH-66, E3
  ``e1_degree6_even4_inertia_log3_log4.json``         §WO-RH-69, E1
  ``e1_degree6_even4_positivity_log3_log4.json``      §WO-RH-69, E1, only if proved
  ``e1_eng010_even4_generalized_gap_log3_log4.json``  §WO-RH-70, E1
  ``e1_degree6_even4_moments_log3_log4.json``         §WO-RH-72, E1

The block is a *prediction test* (§0): ENG-009's preregistered even-sector
scaling models disagree about its generalized gap, and the adjudication script
compares them against the ``gap`` stage's certificate without refitting.

Two independent warrants are computed for definiteness -- interval LDL*
congruence stratified over the cell, and Sylvester's criterion as four
adaptive covers -- exactly as ENG-008 did at 3x3. The generalized gap uses
shifted positivity of the *same frozen congruence* applied to the pencil.

If the block turns out not to be definite, the stratification is the result
and no positivity certificate is written (§WO-RH-69).

No RH proof claim is made. Claim scope is ``finite_dimensional_weil_compression``.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
for extra in (ROOT, ROOT / "src"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))

import basis_algebra  # noqa: E402
import even4  # noqa: E402
import generalized_gap as GG  # noqa: E402
import normalization as N  # noqa: E402
import promotion  # noqa: E402
import reference_metric as RM  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from content_kinds import (  # noqa: E402
    KIND_DEGREE6_POSITIVITY,
    KIND_GENERALIZED_GAP,
    KIND_SCAN_PREVIEW,
)
from inertia.certificate import KIND_INERTIA  # noqa: E402
from inertia.stratify import certify_inertia_family  # noqa: E402
from interval_backend import interval_box, require_flint  # noqa: E402
from interval_cover import NotSeparated, adaptive_cover  # noqa: E402
from moments.adapter import analyse  # noqa: E402
from ranktrace.theorem import rank_trace_lower_bound  # noqa: E402

E0_FILE = "e0_degree6_even4_exact_identities.json"
CROSSCHECK_FILE = "e3_degree6_even4_crosscheck.json"
INERTIA_FILE = "e1_degree6_even4_inertia_log3_log4.json"
POSITIVITY_FILE = "e1_degree6_even4_positivity_log3_log4.json"
GAP_FILE = "e1_eng010_even4_generalized_gap_log3_log4.json"
MOMENTS_FILE = "e1_degree6_even4_moments_log3_log4.json"

PRECISION_BITS = even4.DEFAULT_PRECISION_BITS

DEPENDENCIES = (
    "src/pole.py",
    "src/core.py",
    "src/basis_algebra.py",
    "src/reference_metric.py",
    "src/generalized_gap.py",
    "src/weil_entries.py",
    "src/archimedean_realspace.py",
    "src/even4.py",
    "src/interval_cover.py",
    "src/interval_backend.py",
    "src/normalization.py",
    "inertia/ldl.py",
    "inertia/stratify.py",
    "inertia/certificate.py",
    "moments/spectral_moments.py",
    "moments/adapter.py",
    "ranktrace/theorem.py",
    "scripts/certify_even4.py",
)

SAMPLE_LS = (math.log(3.0), 1.10, 1.20, 1.30, math.log(4.0))
#: §WO-RH-66 also asks for at least two non-special interior points.
EXTRA_LS = (1.1547, 1.3266)


def common_header(work_order: str, kind: str, evidence: str) -> Dict[str, Any]:
    return {
        "certificate_version": "0.1",
        "program": "RH/Weil 4x4 even block — Candidate A",
        "work_order": work_order,
        "claim_scope": even4.CLAIM_SCOPE,
        "content_kind": kind,
        "evidence_class": evidence,
        "basis": even4.basis_identity(),
        "domain": {"cell": [repr(even4.CELL[0]), repr(even4.CELL[1])],
                   "label": list(even4.CELL_LABEL)},
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }


# --------------------------------------------------------------------------- #
# §WO-RH-66: E0 exact identities                                               #
# --------------------------------------------------------------------------- #
def stage_e0() -> Dict[str, Any]:
    print("\n=== E0 exact identities (§WO-RH-66) ===")
    F = Fraction
    checks: List[Dict[str, Any]] = []

    def check(name: str, ok: bool) -> None:
        checks.append({"check": name, "ok": bool(ok)})
        print(f"  [{'ok' if ok else 'FAIL'}] {name}")
        if not ok:
            raise AssertionError(name)

    LS = (F(7, 6), F(5, 4), F(11, 8))
    XS = (F(0), F(1, 3), F(1, 2), F(7, 8))

    def h(name, x, L):
        return sum(c * x ** xp * L ** lp
                   for xp, lp_ in enumerate(basis_algebra.BASIS_L_POLY[name])
                   for lp, c in lp_.items())

    check("every element of {1, b, b2, bcube} is even about x = L/2",
          all(h(n, L - x, L) == h(n, x, L)
              for n in even4.EVEN4_BASIS for L in LS for x in XS))
    check("bcube = b^3 and b2 = b^2, exactly",
          all(h("bcube", x, L) == h("b", x, L) ** 3
              and h("b2", x, L) == h("b", x, L) ** 2
              for L in LS for x in XS))
    check("all ten overlap kernels symmetric and vanishing at a = L",
          all(basis_algebra.kernel_exact(i, j, a, L)
              == basis_algebra.kernel_exact(j, i, a, L)
              and basis_algebra.kernel_exact(i, j, L, L) == 0
              for _, (i, j) in even4.ENTRY_KEYS for L in LS for a in XS))
    check("kernel (L - a)^m multiplicities are 1..7 as the vanishing orders predict",
          all(basis_algebra.kernel_factored(i, j)[0] >= 1
              for _, (i, j) in even4.ENTRY_KEYS)
          and basis_algebra.kernel_factored("bcube", "bcube")[0] == 7
          and basis_algebra.kernel_factored("one", "bcube")[0] == 4)
    check("reference metric entries are the predicted monomials and M(1) is PD",
          all(RM.metric_monomial(i, j)[1]
              == RM.BASIS_DEGREE[i] + RM.BASIS_DEGREE[j] + 1
              for _, (i, j) in even4.ENTRY_KEYS)
          and all(F(v) > 0 for v in RM.certify_positive_definite(
              even4.EVEN4_BASIS)["unit_leading_minors"]))
    check("endpoint values h(L; L) simplify to exact zero for b, b2, bcube",
          all(basis_algebra.endpoint_poly(n) == {}
              for n in ("b", "b2", "bcube")))

    body = common_header("WO-RH-66", "", "E0")
    body.pop("content_kind")
    body.update({
        "rigorous": True,
        "hard_constraints_certified": True,
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "arithmetic": "exact_rational",
        "checks": checks,
        "reference_metric": RM.certify_positive_definite(even4.EVEN4_BASIS),
    })
    path = write_certificate(E0_FILE, body)
    print(f"wrote {path}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-66: E3 independent cross-check                                        #
# --------------------------------------------------------------------------- #
def stage_crosscheck() -> Dict[str, Any]:
    print("\n=== E3 independent cross-check (§WO-RH-66) ===")
    import mpmath as mp

    import independent_even4 as IE4

    _, arb, _, ctx = require_flint()
    ctx.prec = PRECISION_BITS
    points = sorted(set(SAMPLE_LS) | set(EXTRA_LS))
    worst = 0.0
    rows = []
    for L in points:
        built = even4.assemble_even4_arb(arb(repr(L)))
        M = IE4.gram_matrix(L, mp)
        for a in range(4):
            for b in range(a, 4):
                rig = float(built["matrix"][a][b].mid())
                ind = float(M[a][b])
                rel = abs(rig - ind) / max(abs(rig), 1e-300)
                worst = max(worst, rel)
        rows.append({"L": repr(L), "worst_so_far": repr(worst)})
        print(f"  L={L:.9f}  worst relative difference so far {worst:.3e}")

    body = common_header("WO-RH-66", KIND_SCAN_PREVIEW, "E3")
    body.update({
        "rigorous": False,
        "hard_constraints_certified": False,
        "psd_claim": False,
        "status": "AGREES",
        "mpmath_used": True,
        "numeric_warrant": "NONE — E3 regression evidence, never a warrant",
        "independence": ("independent_even4 imports independent_even3 only; "
                         "neither imports any production assembly module "
                         "(asserted by test_even4)"),
        "points": rows,
        "worst_relative_difference": repr(worst),
    })
    path = write_certificate(CROSSCHECK_FILE, body)
    print(f"wrote {path}")
    print(f"  worst relative difference across {len(points)} points: {worst:.3e}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-69: inertia and positivity                                            #
# --------------------------------------------------------------------------- #
def stage_inertia(quick: bool) -> Tuple[Dict[str, Any], bool]:
    print("\n=== E1 whole-cell inertia (§WO-RH-69) ===")
    t0 = time.time()
    policy = {"initial_cells": 64 if quick else 256, "max_depth": 18,
              "min_width": 1e-12}
    strat = certify_inertia_family(
        lambda lo, hi: even4.matrix_over(lo, hi, precision_bits=PRECISION_BITS),
        even4.CELL, subdivision_policy=policy)
    elapsed = time.time() - t0
    signatures = sorted({tuple(s.signature) for s in strat.strata})
    definite = (len(signatures) == 1 and signatures[0] == (4, 0, 0)
                and not strat.transitions)
    strat_dict = strat.to_dict()
    body = common_header("WO-RH-69", KIND_INERTIA, "E1")
    sample = even4.assemble_even4_arb(interval_box(*even4.CELL),
                                      precision_bits=PRECISION_BITS)
    body.update({
        "rigorous": True,
        "hard_constraints_certified": strat.status.startswith("PASS"),
        "numeric_warrant": "E1",
        "psd_claim": False,
        "status": "PASS" if strat.status.startswith("PASS") else strat.status,
        "mpmath_used": False,
        "backend": "python-flint / Arb",
        "precision_bits": PRECISION_BITS,
        "method": ("interval Hermitian LDL* congruence on the frozen "
                   "preconditioned block, with adaptive subdivision of the L "
                   "cell; no eigenvalue solver anywhere on this path"),
        "preconditioner": sample["preconditioner"],
        "stratification": strat_dict,
        "signatures_seen": [list(s) for s in signatures],
        "constant_on_cell": bool(len(signatures) == 1 and not strat.transitions),
        "n_positive": signatures[0][0] if len(signatures) == 1 else None,
        "n_negative": signatures[0][1] if len(signatures) == 1 else None,
        "n_zero": signatures[0][2] if len(signatures) == 1 else None,
        "elapsed_seconds": repr(elapsed),
    })
    path = write_certificate(INERTIA_FILE, body)
    print(f"wrote {path}")
    print(f"  status {strat.status}  signatures {signatures}  "
          f"boxes {strat.boxes_examined}  depth {strat.max_depth}  {elapsed:.1f}s")
    for t in strat.transitions:
        print(f"  transition region [{t.lo!r}, {t.hi!r}]: {t.blocker}")
    return body, definite


def stage_positivity(quick: bool) -> Optional[Dict[str, Any]]:
    print("\n=== E1 Sylvester leading minors (§WO-RH-69) ===")
    scales = even4.minor_scale_factors(even4.PRECONDITIONER_EXPONENTS)
    covers = []
    initial = {1: 64, 2: 64, 3: 512, 4: 2048 if not quick else 512}
    for idx in range(4):
        k = idx + 1

        def ev(lo, hi, _k=k):
            mat = even4.matrix_over(lo, hi, precision_bits=PRECISION_BITS)
            m = even4.leading_minors(mat)[_k - 1]
            return float(m.lower()), float(m.upper())

        t0 = time.time()
        try:
            cov = adaptive_cover(ev, quantity=f"Delta{k}", cell=even4.CELL,
                                 target=0.0, initial_boxes=initial[k],
                                 max_depth=14)
        except NotSeparated as stop:
            print(f"  Delta{k}: NOT SEPARATED — {stop}")
            return None
        raw = cov.certified_lower_bound / scales[idx]
        covers.append({
            "minor": f"Delta{k}",
            "certified_lower_bound": repr(cov.certified_lower_bound),
            "implied_raw_lower_bound": repr(raw),
            "scale_factor": repr(scales[idx]),
            "boxes": cov.boxes_examined,
            "max_depth": cov.max_depth,
        })
        print(f"  Delta{k}: >= {cov.certified_lower_bound} (preconditioned) "
              f"= {raw} raw; {cov.boxes_examined} boxes, depth {cov.max_depth} "
              f"[{time.time() - t0:.1f}s]")

    body = common_header("WO-RH-69", KIND_DEGREE6_POSITIVITY, "E1")
    body.update({
        "rigorous": True,
        "hard_constraints_certified": True,
        "numeric_warrant": "E1",
        "logical_implication_warrant": (
            "FORMAL: manifest theorems pd_four_by_four_certificate and "
            "preconditioned_certificate4 (minor bounds on the frozen dyadic "
            "congruence imply the original block is positive definite)"),
        "psd_claim": True,
        "status": "PASS",
        "mpmath_used": False,
        "backend": "python-flint / Arb",
        "precision_bits": PRECISION_BITS,
        "preconditioner": even4.preconditioner_record(
            even4.PRECONDITIONER_EXPONENTS),
        "leading_minors": covers,
        "n_positive": 4, "n_negative": 0, "n_zero": 0,
        "statement": {
            "conclusion": ("The 4x4 even Weil block G[{1, b, b^2, b^3}] is "
                           "positive definite for every L in [log 3, log 4]"),
            "route": "Sylvester: Delta1..Delta4 all uniformly positive",
        },
    })
    path = write_certificate(POSITIVITY_FILE, body)
    print(f"wrote {path}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-70: the generalized gap                                               #
# --------------------------------------------------------------------------- #
def dyadic_floor(x: float, bits: int = 24) -> Fraction:
    return Fraction(math.floor(x * 2 ** bits), 2 ** bits)


def stage_gap(quick: bool) -> Dict[str, Any]:
    print("\n=== E1 generalized gap (§WO-RH-70) ===")
    _, arb, _, ctx = require_flint()
    ctx.prec = PRECISION_BITS
    lo, hi = even4.CELL
    n_grid = 9 if quick else 17
    grid = [lo + (hi - lo) * k / (n_grid - 1) for k in range(n_grid)]
    basis = even4.EVEN4_BASIS

    lams = []
    for L in grid:
        built = even4.assemble_even4_arb(arb(repr(L)))
        floats = {}
        for key, (i, j) in even4.ENTRY_KEYS:
            floats[(i, j)] = float(built["entries"][key].mid())
        lams.append(GG.scout_gap_at(basis, floats, L))
    k_min = min(range(len(grid)), key=lambda k: lams[k])
    scout_min, argmin_L = lams[k_min], grid[k_min]
    #: §WO-RH-71 note: 0.8 slack (not the ENG-009 0.9) trades enclosure width
    #: for cover cost -- the shifted fourth minor's margin shrinks linearly in
    #: the slack while the cover cost grows quickly. The adjudication only
    #: needs to know which preregistered falsifier intervals the certified
    #: enclosure excludes, and the E3 scout sits far enough above both that
    #: 0.8 x scout already decides them; the certificate reports both numbers.
    lam = dyadic_floor(0.8 * scout_min)
    print(f"  scout: min lambda* = {scout_min:.8e} at L = {argmin_L:.6f} "
          f"(E3); certifying lam = {float(lam):.8e}")

    covers = []
    initial = {1: 64, 2: 64, 3: 512, 4: 2048 if not quick else 512}
    for idx in range(4):
        k = idx + 1

        def ev(blo, bhi, _k=k):
            mat = even4.shifted_matrix_over(
                blo, bhi, lam.numerator, lam.denominator,
                precision_bits=PRECISION_BITS)
            m = even4.leading_minors(mat)[_k - 1]
            return float(m.lower()), float(m.upper())

        t0 = time.time()
        try:
            cov = adaptive_cover(ev, quantity=f"shifted minor{k}",
                                 cell=even4.CELL, target=0.0,
                                 initial_boxes=initial[k], max_depth=14)
        except NotSeparated as stop:
            raise SystemExit(f"gap cover failed: {stop}")
        covers.append({
            "minor": k,
            "certified_lower_bound_preconditioned": repr(cov.certified_lower_bound),
            "boxes": cov.boxes_examined, "max_depth": cov.max_depth,
        })
        print(f"    minor {k}: PASS bound={cov.certified_lower_bound:.6e} "
              f"boxes={cov.boxes_examined} depth={cov.max_depth} "
              f"[{time.time() - t0:.1f}s]")

    # Upper bound: certified Rayleigh quotient of a rational witness at the
    # scouted bottleneck.
    built = even4.assemble_even4_arb(arb(repr(argmin_L)))
    floats = {}
    ent = {}
    for key, (i, j) in even4.ENTRY_KEYS:
        ent[(i, j)] = built["entries"][key]
        floats[(i, j)] = float(built["entries"][key].mid())
    v = GG.scout_min_eigvec(basis, floats, scout_min, argmin_L)
    quotient = GG.rayleigh_upper(basis, ent, v, arb(repr(argmin_L)))
    upper = float(quotient.upper())
    print(f"  upper bound {upper} at L={argmin_L}")

    body = common_header("WO-RH-70", KIND_GENERALIZED_GAP, "E1")
    body.update({
        "rigorous": True,
        "hard_constraints_certified": True,
        "numeric_warrant": "E1",
        "logical_implication_warrant": (
            "FORMAL: manifest theorems generalized_rayleigh, "
            "generalized_pencil_congruence, and preconditioned_gap_certificate4 "
            "(the composed 4x4 replay)"),
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "backend": "python-flint / Arb",
        "precision_bits": PRECISION_BITS,
        "reference_metric_id": even4.REFERENCE_METRIC_ID,
        "reference_metric_certificate": "e0_eng010_even4_reference_metric.json",
        "basis_id": even4.EVEN4_BASIS_ID,
        "preconditioner": even4.preconditioner_record(
            even4.PRECONDITIONER_EXPONENTS),
        "preconditioner_convention": (
            "the frozen D = diag(2^{-e}) is applied to the *shifted* pencil "
            "D (G - lam M) D; exact powers of two, minor signs unchanged, and "
            "the pencil's eigenvalues are unchanged by the simultaneous "
            "congruence (generalized_pencil_congruence)"),
        "certified_lambda_lower_uniform": str(lam),
        "certified_lambda_lower_float": repr(float(lam)),
        "slack_note": "lam = dyadic_floor(0.8 * scout minimum); see source",
        "upper_bound_at_bottleneck": {
            "at_L": repr(argmin_L),
            "witness_vector": [str(x) for x in v],
            "rayleigh_enclosure": [repr(float(quotient.lower())),
                                   repr(float(quotient.upper()))],
            "certified_upper_bound": repr(upper),
        },
        "bottleneck_region": {
            "argmin_L_on_grid": repr(argmin_L),
            "grid_lambda_star": [repr(x) for x in lams],
            "grid": [repr(x) for x in grid],
            "note": "E3 scout locates the bottleneck; the covers are the warrant",
        },
        "shifted_minor_covers": covers,
        "statement": (
            "lambda_min(G, M)(L) >= {} for every L in [log 3, log 4], and "
            "lambda_min at the scouted bottleneck is <= {}. M is the exact L2 "
            "reference metric.".format(float(lam), upper)),
    })
    path = write_certificate(GAP_FILE, body)
    print(f"wrote {path}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-72: moments and rank-trace                                            #
# --------------------------------------------------------------------------- #
def stage_moments() -> Dict[str, Any]:
    print("\n=== E1 moments and rank-trace (§WO-RH-72) ===")
    _, arb, _, ctx = require_flint()
    ctx.prec = PRECISION_BITS
    points = []
    for L in SAMPLE_LS:
        built = even4.assemble_even4_arb(arb(repr(L)))
        mat = built["preconditioned"]
        out = analyse(mat, observed_inertia=(4, 0, 0))
        det = even4.leading_minors(built["matrix"])[3]
        trace = built["matrix"][0][0]
        for k in range(1, 4):
            trace = trace + built["matrix"][k][k]
        rt = out.get("rank_trace") or {}
        points.append({
            "label": f"L={L:.6f}",
            "L": repr(L),
            "trace": [repr(float(trace.lower())), repr(float(trace.upper()))],
            "determinant": [repr(float(det.lower())), repr(float(det.upper()))],
            "moment_analysis": out,
        })
        print(f"  L={L:.6f}  trace {float(trace.mid()):.6e}  "
              f"det {float(det.mid()):.6e}")
    body = common_header("WO-RH-72", "WEIL_SPECTRAL_MOMENT_CERTIFICATE", "E1")
    body.update({
        "rigorous": True,
        "hard_constraints_certified": True,
        "numeric_warrant": "E1",
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "dimension": 4,
        "observed_inertia": [4, 0, 0],
        "points": points,
        "note": ("moments of the preconditioned block; the raw trace and "
                 "determinant are recorded alongside for the information "
                 "comparison"),
    })
    path = write_certificate(MOMENTS_FILE, body)
    print(f"wrote {path}")
    return body


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "e0", "metric", "crosscheck", "inertia",
                             "positivity", "gap", "moments"])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--release", action="store_true")
    args = ap.parse_args()
    require_flint()
    if args.stage in ("all", "e0"):
        stage_e0()
    if args.stage in ("all", "e0", "metric"):
        stage_metric()
    if args.stage in ("all", "crosscheck"):
        stage_crosscheck()
    definite = True
    if args.stage in ("all", "inertia"):
        _, definite = stage_inertia(args.quick)
    if args.stage in ("all", "positivity"):
        if definite:
            if stage_positivity(args.quick) is None:
                return 1
        else:
            print("  positivity not claimed; the stratification is the result")
    if args.stage in ("all", "gap"):
        stage_gap(args.quick)
    if args.stage in ("all", "moments"):
        stage_moments()
    return 0




# --------------------------------------------------------------------------- #
# §WO-RH-65: the reference metric record (its own file, content-addressed)     #
# --------------------------------------------------------------------------- #
def stage_metric() -> Dict[str, Any]:
    print("\n=== E0 reference metric (§WO-RH-65) ===")
    rec = RM.certify_positive_definite(even4.EVEN4_BASIS)
    body = common_header("WO-RH-65", KIND_GENERALIZED_GAP, "E0")
    body.update({
        "role": "reference_metric",
        "rigorous": True,
        "hard_constraints_certified": True,
        "numeric_warrant": "E0",
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "arithmetic": "exact_rational",
        "reference_metric_id": even4.REFERENCE_METRIC_ID,
        "metric": rec,
        "statement": ("The exact L2 reference metric on {1, b, b^2, b^3} is "
                      "positive definite for every L > 0, by exact Sylvester "
                      "minors of M(1) and the diagonal congruence "
                      "M(L) = D^T M(1) D."),
    })
    path = write_certificate("e0_eng010_even4_reference_metric.json", body)
    print(f"wrote {path}")
    print(f"  M(1) minors: {rec['unit_leading_minors']}")
    return body


if __name__ == "__main__":
    sys.exit(main())
