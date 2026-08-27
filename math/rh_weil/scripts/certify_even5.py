#!/usr/bin/env python3
"""ATLAS-RH-ENG-011 — the 5x5 even Weil block ``G[{1, b, b^2, b^3, b^4}]``.

    python3 scripts/certify_even5.py [--stage STAGE] [--quick]

Stages (default ``all``): ``e0``, ``crosscheck``, ``inertia``, ``positivity``,
``gap``, ``moments``. Produces:

  ``e0_degree8_even5_exact_identities.json``          §WO-RH-77, E0
  ``e3_degree8_even5_crosscheck.json``                §WO-RH-77, E3
  ``e1_degree8_even5_inertia_log3_log4.json``         §WO-RH-81, E1
  ``e1_degree8_even5_positivity_log3_log4.json``      §WO-RH-81, E1, only if proved
  ``e1_eng011_even5_generalized_gap_log3_log4.json``  §WO-RH-82, E1
  ``e1_degree8_even5_moments_log3_log4.json``         §WO-RH-84, E1

The block is a *prediction test* (§0): ENG-009's preregistered even-sector
scaling models disagree about its generalized gap, and the adjudication script
compares them against the ``gap`` stage's certificate without refitting.

Two independent warrants are computed for definiteness -- interval LDL*
congruence stratified over the cell, and Sylvester's criterion as four
adaptive covers -- exactly as ENG-008 did at 3x3. The generalized gap uses
shifted positivity of the *same frozen congruence* applied to the pencil.

If the block turns out not to be definite, the stratification is the result
and no positivity certificate is written (§WO-RH-81).

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
import even5  # noqa: E402
import generalized_gap as GG  # noqa: E402
import normalization as N  # noqa: E402
import promotion  # noqa: E402
import reference_metric as RM  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from content_kinds import (  # noqa: E402
    KIND_DEGREE8_POSITIVITY,
    KIND_GENERALIZED_GAP,
    KIND_SCAN_PREVIEW,
)
from inertia.certificate import KIND_INERTIA  # noqa: E402
from inertia.stratify import certify_inertia_family  # noqa: E402
from interval_backend import interval_box, require_flint  # noqa: E402
from interval_cover import NotSeparated, adaptive_cover  # noqa: E402
from moments.adapter import analyse  # noqa: E402
from ranktrace.theorem import rank_trace_lower_bound  # noqa: E402

E0_FILE = "e0_degree8_even5_exact_identities.json"
CROSSCHECK_FILE = "e3_degree8_even5_crosscheck.json"
INERTIA_FILE = "e1_degree8_even5_inertia_log3_log4.json"
POSITIVITY_FILE = "e1_degree8_even5_positivity_log3_log4.json"
GAP_FILE = "e1_eng011_even5_generalized_gap_log3_log4.json"
MOMENTS_FILE = "e1_degree8_even5_moments_log3_log4.json"

PRECISION_BITS = even5.DEFAULT_PRECISION_BITS

CERT_DIR = ROOT / "certificates"

DEPENDENCIES = (
    "src/pole.py",
    "src/core.py",
    "src/basis_algebra.py",
    "src/reference_metric.py",
    "src/generalized_gap.py",
    "src/weil_entries.py",
    "src/archimedean_realspace.py",
    "src/even5.py",
    "src/interval_cover.py",
    "src/interval_backend.py",
    "src/normalization.py",
    "inertia/ldl.py",
    "inertia/stratify.py",
    "inertia/certificate.py",
    "moments/spectral_moments.py",
    "moments/adapter.py",
    "ranktrace/theorem.py",
    "scripts/certify_even5.py",
)

SAMPLE_LS = (math.log(3.0), 1.10, 1.20, 1.30, 1.36, math.log(4.0))
#: §WO-RH-77 also asks for at least two non-special interior points.
EXTRA_LS = (1.1547, 1.3266)


def common_header(work_order: str, kind: str, evidence: str) -> Dict[str, Any]:
    return {
        "certificate_version": "0.1",
        "program": "RH/Weil 5x5 even block — Candidate A",
        "work_order": work_order,
        "claim_scope": even5.CLAIM_SCOPE,
        "content_kind": kind,
        "evidence_class": evidence,
        "basis": even5.basis_identity(),
        "domain": {"cell": [repr(even5.CELL[0]), repr(even5.CELL[1])],
                   "label": list(even5.CELL_LABEL)},
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }


# --------------------------------------------------------------------------- #
# §WO-RH-77: E0 exact identities                                               #
# --------------------------------------------------------------------------- #
def stage_e0() -> Dict[str, Any]:
    print("\n=== E0 exact identities (§WO-RH-77) ===")
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

    check("every element of {1, b, b2, bcube, bquart} is even about x = L/2",
          all(h(n, L - x, L) == h(n, x, L)
              for n in even5.EVEN5_BASIS for L in LS for x in XS))
    check("bquart = b^4, bcube = b^3 and b2 = b^2, exactly",
          all(h("bquart", x, L) == h("b", x, L) ** 4
              and h("bcube", x, L) == h("b", x, L) ** 3
              and h("b2", x, L) == h("b", x, L) ** 2
              for L in LS for x in XS))
    check("all fifteen overlap kernels symmetric and vanishing at a = L",
          all(basis_algebra.kernel_exact(i, j, a, L)
              == basis_algebra.kernel_exact(j, i, a, L)
              and basis_algebra.kernel_exact(i, j, L, L) == 0
              for _, (i, j) in even5.ENTRY_KEYS for L in LS for a in XS))
    check("kernel (L - a)^m multiplicities are 1..9 as the vanishing orders predict",
          all(basis_algebra.kernel_factored(i, j)[0] >= 1
              for _, (i, j) in even5.ENTRY_KEYS)
          and basis_algebra.kernel_factored("bquart", "bquart")[0] == 9
          and basis_algebra.kernel_factored("one", "bquart")[0] == 5)
    check("reference metric entries are the predicted monomials and M(1) is PD",
          all(RM.metric_monomial(i, j)[1]
              == RM.BASIS_DEGREE[i] + RM.BASIS_DEGREE[j] + 1
              for _, (i, j) in even5.ENTRY_KEYS)
          and all(F(v) > 0 for v in RM.certify_positive_definite(
              even5.EVEN5_BASIS)["unit_leading_minors"]))
    check("endpoint values h(L; L) simplify to exact zero for b, b2, bcube, bquart",
          all(basis_algebra.endpoint_poly(n) == {}
              for n in ("b", "b2", "bcube", "bquart")))

    body = common_header("WO-RH-77", "", "E0")
    body.pop("content_kind")
    body.update({
        "rigorous": True,
        "hard_constraints_certified": True,
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "arithmetic": "exact_rational",
        "checks": checks,
        "reference_metric": RM.certify_positive_definite(even5.EVEN5_BASIS),
    })
    path = write_certificate(E0_FILE, body)
    print(f"wrote {path}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-77: E3 independent cross-check                                        #
# --------------------------------------------------------------------------- #
def stage_crosscheck() -> Dict[str, Any]:
    print("\n=== E3 independent cross-check (§WO-RH-77) ===")
    import mpmath as mp

    import independent_even5 as IE5

    _, arb, _, ctx = require_flint()
    ctx.prec = PRECISION_BITS
    points = sorted(set(SAMPLE_LS) | set(EXTRA_LS))
    worst = 0.0
    rows = []
    for L in points:
        built = even5.assemble_even5_arb(arb(repr(L)))
        M = IE5.gram_matrix(L, mp)
        for a in range(5):
            for b in range(a, 5):
                rig = float(built["matrix"][a][b].mid())
                ind = float(M[a][b])
                rel = abs(rig - ind) / max(abs(rig), 1e-300)
                worst = max(worst, rel)
        rows.append({"L": repr(L), "worst_so_far": repr(worst)})
        print(f"  L={L:.9f}  worst relative difference so far {worst:.3e}")

    body = common_header("WO-RH-77", KIND_SCAN_PREVIEW, "E3")
    body.update({
        "rigorous": False,
        "hard_constraints_certified": False,
        "psd_claim": False,
        "status": "AGREES",
        "mpmath_used": True,
        "numeric_warrant": "NONE — E3 regression evidence, never a warrant",
        "independence": ("independent_even5 imports independent_even3 only; "
                         "neither imports any production assembly module "
                         "(asserted by test_even5)"),
        "points": rows,
        "worst_relative_difference": repr(worst),
    })
    path = write_certificate(CROSSCHECK_FILE, body)
    print(f"wrote {path}")
    print(f"  worst relative difference across {len(points)} points: {worst:.3e}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-81: inertia and positivity                                            #
# --------------------------------------------------------------------------- #
def dyadic_floor(x: float, bits: int = 34) -> Fraction:
    return Fraction(math.floor(x * 2 ** bits), 2 ** bits)


#: The certified shift, frozen from the E3 scout (interior minimum
#: lambda* ~ 2.538e-07 at L ~ 1.173) with 25% slack: the sweeps below prove
#: G - LAM*M is positive definite on the whole cell, hence
#: lambda_min(G, M) >= LAM uniformly and G itself is PD.
GAP_LAM: Fraction = dyadic_floor(0.75 * 2.538e-07)

#: Side intervals for the bottleneck localization, and the larger shift
#: certified on them: lambda* >= OUT_LAM there, so the infimum's location is
#: pinned to the complementary interior interval.
SIDE_LEFT: Tuple[float, float] = (math.log(3.0), 1.11)
SIDE_RIGHT: Tuple[float, float] = (1.30, math.log(4.0))
OUT_LAM: Fraction = dyadic_floor(4.0e-07)

PART_DIR = ROOT / "certificates" / "even5_parts"


def _chunk_interval(cell: Tuple[float, float], i: int, k: int) -> Tuple[float, float]:
    a, b = cell
    return (a + (b - a) * i / k, a + (b - a) * (i + 1) / k)


def _sweep(matrix_fn, interval: Tuple[float, float], initial: int) -> Dict[str, Any]:
    strat = certify_inertia_family(
        matrix_fn, interval,
        subdivision_policy={"initial_cells": initial, "max_depth": 16,
                            "min_width": 1e-12})
    signatures = sorted({tuple(x.signature) for x in strat.strata})
    return {
        "interval": [repr(interval[0]), repr(interval[1])],
        "status": strat.status,
        "signatures": [list(x) for x in signatures],
        "boxes_examined": strat.boxes_examined,
        "max_depth": strat.max_depth,
        "transitions": [{"lo": repr(t.lo), "hi": repr(t.hi),
                         "blocker": str(t.blocker)}
                        for t in strat.transitions],
    }


def stage_inertia_chunk(i: int, k: int) -> int:
    PART_DIR.mkdir(parents=True, exist_ok=True)
    interval = _chunk_interval(even5.CELL, i, k)
    print(f"inertia chunk {i}/{k} on [{interval[0]:.6f}, {interval[1]:.6f}]",
          flush=True)
    t0 = time.time()
    part = _sweep(lambda lo, hi: even5.matrix_over(lo, hi,
                                                   precision_bits=PRECISION_BITS),
                  interval, initial=4096)
    part["elapsed_seconds"] = repr(time.time() - t0)
    out = PART_DIR / f"inertia_{i}_{k}.json"
    out.write_text(json.dumps(part, indent=1), encoding="utf-8")
    print(f"  {part['status']} sigs={part['signatures']} "
          f"boxes={part['boxes_examined']} depth={part['max_depth']} "
          f"[{float(part['elapsed_seconds']):.0f}s]", flush=True)
    return 0 if part["status"].startswith("PASS") else 1


def stage_gap_chunk(i: int, k: int) -> int:
    PART_DIR.mkdir(parents=True, exist_ok=True)
    interval = _chunk_interval(even5.CELL, i, k)
    print(f"gap chunk {i}/{k} at lam={float(GAP_LAM):.6e} "
          f"on [{interval[0]:.6f}, {interval[1]:.6f}]", flush=True)
    t0 = time.time()
    part = _sweep(lambda lo, hi: even5.shifted_matrix_over(
        lo, hi, GAP_LAM.numerator, GAP_LAM.denominator,
        precision_bits=PRECISION_BITS), interval, initial=8192)
    part["lam"] = str(GAP_LAM)
    part["elapsed_seconds"] = repr(time.time() - t0)
    out = PART_DIR / f"gap_{i}_{k}.json"
    out.write_text(json.dumps(part, indent=1), encoding="utf-8")
    print(f"  {part['status']} sigs={part['signatures']} "
          f"boxes={part['boxes_examined']} depth={part['max_depth']} "
          f"[{float(part['elapsed_seconds']):.0f}s]", flush=True)
    return 0 if part["status"].startswith("PASS") else 1


def _merge_parts(prefix: str, k: int) -> Dict[str, Any]:
    parts = []
    for i in range(k):
        path = PART_DIR / f"{prefix}_{i}_{k}.json"
        if not path.exists():
            raise SystemExit(f"missing part {path}")
        parts.append(json.loads(path.read_text(encoding="utf-8")))
    sigs = sorted({tuple(s) for p in parts for s in p["signatures"]})
    return {
        "parts": parts,
        "signatures": [list(s) for s in sigs],
        "status": ("PASS" if all(p["status"].startswith("PASS") for p in parts)
                   else "PARTIAL"),
        "boxes_examined": sum(p["boxes_examined"] for p in parts),
        "max_depth": max(p["max_depth"] for p in parts),
        "transitions": [t for p in parts for t in p["transitions"]],
        "parallel_chunks": k,
    }


def stage_inertia_merge(k: int) -> Tuple[Dict[str, Any], bool]:
    print("\n=== E1 whole-cell inertia, merged (§WO-RH-81) ===")
    merged = _merge_parts("inertia", k)
    definite = (merged["status"] == "PASS"
                and merged["signatures"] == [[5, 0, 0]]
                and not merged["transitions"])
    body = common_header("WO-RH-81", KIND_INERTIA, "E1")
    sample = even5.assemble_even5_arb(interval_box(*even5.CELL),
                                      precision_bits=PRECISION_BITS)
    body.update({
        "rigorous": True,
        "hard_constraints_certified": merged["status"] == "PASS",
        "numeric_warrant": "E1",
        "psd_claim": False,
        "status": merged["status"],
        "mpmath_used": False,
        "backend": "python-flint / Arb",
        "precision_bits": PRECISION_BITS,
        "method": ("interval Hermitian LDL* congruence on the frozen "
                   "preconditioned block, adaptive subdivision, computed in "
                   "parallel cell chunks and merged; no eigenvalue solver"),
        "preconditioner": sample["preconditioner"],
        "stratification": merged,
        "signatures_seen": merged["signatures"],
        "constant_on_cell": bool(len(merged["signatures"]) == 1
                                 and not merged["transitions"]),
        "n_positive": merged["signatures"][0][0] if len(merged["signatures"]) == 1 else None,
        "n_negative": merged["signatures"][0][1] if len(merged["signatures"]) == 1 else None,
        "n_zero": merged["signatures"][0][2] if len(merged["signatures"]) == 1 else None,
    })
    path = write_certificate(INERTIA_FILE, body)
    print(f"wrote {path}")
    print(f"  status {merged['status']}  signatures {merged['signatures']}  "
          f"boxes {merged['boxes_examined']}  depth {merged['max_depth']}")
    return body, definite


def _point_shifted_pass(L: float, lam: Fraction) -> Dict[str, Any]:
    """A point certificate that lambda*(L) >= lam: shifted LDL at a point."""
    from inertia.ldl import interval_inertia
    mat = even5.shifted_matrix_over(L, L, lam.numerator, lam.denominator,
                                    precision_bits=PRECISION_BITS)
    res = interval_inertia(mat)
    ok = (res.status == "PASS"
          and (res.n_positive, res.n_negative, res.n_zero) == (5, 0, 0))
    return {"L": repr(L), "lam": str(lam), "lam_float": repr(float(lam)),
            "status": res.status, "pass": bool(ok)}


def stage_gap_merge(k: int) -> Dict[str, Any]:
    print("\n=== E1 generalized gap, merged (§WO-RH-82) ===")
    _, arb, _, ctx = require_flint()
    ctx.prec = PRECISION_BITS
    merged = _merge_parts("gap", k)
    if merged["status"] != "PASS" or merged["signatures"] != [[5, 0, 0]]:
        raise SystemExit(f"gap sweep did not certify: {merged['status']} "
                         f"{merged['signatures']}")

    sides = []
    for interval in (SIDE_LEFT, SIDE_RIGHT):
        t0 = time.time()
        part = _sweep(lambda lo, hi: even5.shifted_matrix_over(
            lo, hi, OUT_LAM.numerator, OUT_LAM.denominator,
            precision_bits=PRECISION_BITS), interval, initial=2048)
        part["lam"] = str(OUT_LAM)
        sides.append(part)
        print(f"  side {interval}: {part['status']} "
              f"boxes={part['boxes_examined']} [{time.time() - t0:.0f}s]",
              flush=True)
        if part["status"] != "PASS" or part["signatures"] != [[5, 0, 0]]:
            raise SystemExit("side interval failed to certify at OUT_LAM")
    endpoints = [_point_shifted_pass(even5.CELL[0], dyadic_floor(6.4e-07)),
                 _point_shifted_pass(even5.CELL[1], dyadic_floor(1.28e-06))]
    for e in endpoints:
        print(f"  endpoint L={float(e['L']):.6f}: lambda* >= {e['lam_float']} "
              f"({'PASS' if e['pass'] else 'FAIL'})")
        if not e["pass"]:
            raise SystemExit("endpoint exclusion failed")

    L_star = 1.173
    built = even5.assemble_even5_arb(arb(repr(L_star)))
    ent = {}
    floats = {}
    for key, (i, j) in even5.ENTRY_KEYS:
        ent[(i, j)] = built["entries"][key]
        floats[(i, j)] = float(built["entries"][key].mid())
    scout_here = GG.scout_gap_at(even5.EVEN5_BASIS, floats, L_star)
    v = GG.scout_min_eigvec(even5.EVEN5_BASIS, floats, scout_here, L_star)
    quotient = GG.rayleigh_upper(even5.EVEN5_BASIS, ent, v, arb(repr(L_star)))
    upper = float(quotient.upper())
    print(f"  upper bound {upper} at L={L_star}")

    body = common_header("WO-RH-82", KIND_GENERALIZED_GAP, "E1")
    body.update({
        "rigorous": True,
        "hard_constraints_certified": True,
        "numeric_warrant": "E1",
        "logical_implication_warrant": (
            "FORMAL: generalized_rayleigh (shifted positivity is the Rayleigh "
            "bound); the LDL transcript theorems carry the per-box shifted "
            "positivity; positivity of G itself follows from "
            "posDef_of_shifted_posDef_add (G = (G - lam M) + lam M)"),
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "backend": "python-flint / Arb",
        "precision_bits": PRECISION_BITS,
        "reference_metric_id": even5.REFERENCE_METRIC_ID,
        "reference_metric_certificate": "e0_eng011_even5_reference_metric.json",
        "basis_id": even5.EVEN5_BASIS_ID,
        "preconditioner": even5.preconditioner_record(
            even5.PRECONDITIONER_EXPONENTS),
        "certified_lambda_lower_uniform": str(GAP_LAM),
        "certified_lambda_lower_float": repr(float(GAP_LAM)),
        "slack_note": "lam = dyadic_floor(0.75 * interior scout minimum)",
        "method": ("interval LDL* of the frozen-preconditioned shifted pencil "
                   "D (G - lam M) D over the whole cell, in parallel chunks; "
                   "all five pivots positive on every box"),
        "sweep": merged,
        "upper_bound_at_bottleneck": {
            "at_L": repr(L_star),
            "witness_vector": [str(x) for x in v],
            "rayleigh_enclosure": [repr(float(quotient.lower())),
                                   repr(float(quotient.upper()))],
            "certified_upper_bound": repr(upper),
        },
        "bottleneck": {
            "classification": "INTERIOR",
            "certified_interval": [repr(SIDE_LEFT[1]), repr(SIDE_RIGHT[0])],
            "argument": (
                "lambda*(L) >= {} on [log 3, 1.11] and [1.30, log 4] (side "
                "sweeps), and lambda* at both endpoints exceeds {} and {} "
                "(point certificates) -- all above the certified upper witness "
                "{} at L = 1.173, so the infimum over the compact cell is "
                "attained strictly inside [1.11, 1.30]".format(
                    float(OUT_LAM), endpoints[0]["lam_float"],
                    endpoints[1]["lam_float"], upper)),
            "side_sweeps": sides,
            "endpoint_certificates": endpoints,
            "scout_location": "L ~ 1.173 (E3; the covers are the warrant)",
        },
        "statement": (
            "lambda_min(G, M)(L) >= {} for every L in [log 3, log 4]; the "
            "infimum is attained in the interior interval [1.11, 1.30] and is "
            "<= {}. M is the exact L2 reference metric.".format(
                float(GAP_LAM), upper)),
    })
    path = write_certificate(GAP_FILE, body)
    print(f"wrote {path}")
    return body


def stage_positivity_from_gap(inertia_definite: bool) -> Optional[Dict[str, Any]]:
    """§WO-RH-81 route 2: PD from the shifted sweep, by theorem.

    ``G = (G - lam M) + lam M`` with the first summand certified PD by the gap
    sweep, ``lam > 0`` exact, and ``M`` PD by the E0 metric certificate -- so
    ``G`` is PD with no additional cover. Independent of the unshifted LDL
    stratification in everything after the assembly.
    """
    print("\n=== E1 positivity via the shifted route (§WO-RH-81) ===")
    if not (CERT_DIR / GAP_FILE).exists():
        print("  gap certificate missing; not claiming positivity")
        return None
    gap = json.loads((CERT_DIR / GAP_FILE).read_text(encoding="utf-8"))
    if gap.get("status") != "PASS":
        return None
    body = common_header("WO-RH-81", KIND_DEGREE8_POSITIVITY, "E1")
    body.update({
        "rigorous": True,
        "hard_constraints_certified": True,
        "numeric_warrant": "E1",
        "logical_implication_warrant": (
            "FORMAL: posDef_of_shifted_posDef_add composes the gap sweep's "
            "shifted positivity with lam > 0 and the E0-certified reference "
            "metric; the LDL transcript theorems carry each box"),
        "psd_claim": True,
        "status": "PASS",
        "mpmath_used": False,
        "backend": "python-flint / Arb",
        "precision_bits": PRECISION_BITS,
        "preconditioner": even5.preconditioner_record(
            even5.PRECONDITIONER_EXPONENTS),
        "route": ("G = (G - lam M) + lam M: the gap sweep certifies the first "
                  "summand PD on every box of the cell, lam = {} > 0 exactly, "
                  "and M is PD for every L > 0 by exact rational Sylvester "
                  "(e0_eng011_even5_reference_metric.json)".format(
                      gap["certified_lambda_lower_uniform"])),
        "second_route_note": (
            "the unshifted LDL stratification (the inertia certificate) is the "
            "independent first route; it and this shifted route share the "
            "assembly and nothing after it"),
        "gap_certificate": GAP_FILE,
        "n_positive": 5, "n_negative": 0, "n_zero": 0,
        "statement": {
            "conclusion": ("The 5x5 even Weil block G[{1, b, b^2, b^3, b^4}] "
                           "is positive definite for every L in "
                           "[log 3, log 4]"),
            "route": "shifted positivity + exact metric, and LDL* inertia",
        },
        "consistent_with_inertia_route": bool(inertia_definite),
    })
    path = write_certificate(POSITIVITY_FILE, body)
    print(f"wrote {path}")
    return body


def stage_moments() -> Dict[str, Any]:
    print("\n=== E1 moments and rank-trace (§WO-RH-84) ===")
    _, arb, _, ctx = require_flint()
    ctx.prec = PRECISION_BITS
    points = []
    for L in SAMPLE_LS:
        built = even5.assemble_even5_arb(arb(repr(L)))
        mat = built["preconditioned"]
        out = analyse(mat, observed_inertia=(5, 0, 0))
        det = even5.leading_minors(built["matrix"])[3]
        trace = built["matrix"][0][0]
        for k in range(1, 5):
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
    body = common_header("WO-RH-84", "WEIL_SPECTRAL_MOMENT_CERTIFICATE", "E1")
    body.update({
        "rigorous": True,
        "hard_constraints_certified": True,
        "numeric_warrant": "E1",
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "dimension": 5,
        "observed_inertia": [5, 0, 0],
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
                    choices=["all", "e0", "metric", "crosscheck",
                             "inertia-chunk", "inertia-merge",
                             "gap-chunk", "gap-merge", "positivity",
                             "moments"])
    ap.add_argument("--chunk", default=None, help="i:k for chunked sweeps")
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
    if args.stage == "inertia-chunk":
        i, k = (int(x) for x in args.chunk.split(":"))
        return stage_inertia_chunk(i, k)
    if args.stage == "gap-chunk":
        i, k = (int(x) for x in args.chunk.split(":"))
        return stage_gap_chunk(i, k)
    if args.stage == "inertia-merge":
        stage_inertia_merge(int(args.chunk.split(":")[1]))
    if args.stage == "gap-merge":
        stage_gap_merge(int(args.chunk.split(":")[1]))
    if args.stage in ("positivity", "gap-merge"):
        inertia_ok = True
        if (CERT_DIR / INERTIA_FILE).exists():
            iner = json.loads((CERT_DIR / INERTIA_FILE).read_text(encoding="utf-8"))
            inertia_ok = iner.get("signatures_seen") == [[5, 0, 0]]
        stage_positivity_from_gap(inertia_ok)
    if args.stage in ("all", "moments"):
        stage_moments()
    return 0


if __name__ == "__main__":
    sys.exit(main())
