#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §15 (WO-RH-46) -- E3 preview of the 3x3 even parity block.

Prepares, for ENG-008, the first genuinely three-dimensional parity block:
the even sector extended from ``{1, b}`` to ``{1, b, b2}`` with
``b2(x) = x^2 (L-x)^2``. See ``src/pilot3.py`` for why the even sector and not
the odd one.

This produces **E3 evidence only**. It runs in mpmath, uses a floating
eigenvalue solver for the conditioning report, and promotes nothing. §15
forbids a new E1 degree result in ENG-007 and none is claimed here.

    python3 scripts/preview_pilot3.py [--points N] [--T 84] [--dps 30]

The artifact answers the three questions ENG-008 has to make a plan from:

  1. is the block definite across [log 3, log 4], and does its inertia change?
  2. how badly conditioned is it, and does a symmetric rescaling fix that?
  3. does the third leading minor carry information the 2x2 determinant does
     not -- i.e. is this block actually worth the inertia and moment machinery?
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for extra in (ROOT, ROOT / "src"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))

import pilot3 as P  # noqa: E402
import promotion  # noqa: E402
from certificate_io import write_certificate  # noqa: E402

CERT_NAME = "e3_pilot3_even_conditioning_log3_log4.json"
DEPENDENCIES = ("src/pilot3.py", "scripts/preview_pilot3.py")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", type=int, default=9)
    ap.add_argument("--T", type=float, default=84.0)
    ap.add_argument("--dps", type=int, default=30)
    args = ap.parse_args()

    mp = P.require_mpmath()
    mp.mp.dps = args.dps
    lo, hi = P.CELL

    grid = [lo + (hi - lo) * k / (args.points - 1) for k in range(args.points)]
    rows = []
    inertias = set()
    minor3_min = None
    cond_max = 0.0
    rescaled_max = 0.0
    for L in grid:
        M = P.gram_matrix_mp(P.EVEN_BASIS, L, args.T, mp)
        minors = [float(v) for v in P.leading_minors(M)]
        rep = P.condition_report(M, mp)
        scaled, scales = P.jacobi_rescale(M, mp)
        srep = P.condition_report(scaled, mp)
        inertia = (rep["n_positive"], rep["n_negative"], rep["n_zero"])
        inertias.add(inertia)
        minor3_min = minors[2] if minor3_min is None else min(minor3_min, minors[2])
        cond_max = max(cond_max, rep["condition_number"] or 0.0)
        rescaled_max = max(rescaled_max, srep["condition_number"] or 0.0)
        rows.append({
            "L": repr(L),
            "leading_minors": [repr(v) for v in minors],
            "eigenvalues": [repr(v) for v in rep["eigenvalues"]],
            "inertia": list(inertia),
            "condition_number": repr(rep["condition_number"]),
            "diagonal_spread": repr(rep["diagonal_spread"]),
            "jacobi_scales": [repr(v) for v in scales],
            "condition_number_after_jacobi": repr(srep["condition_number"]),
            "inertia_after_jacobi": [srep["n_positive"], srep["n_negative"], srep["n_zero"]],
        })
        print(f"  L={L:.9f}  minors={['%.6e' % v for v in minors]}  "
              f"inertia={inertia}  cond={rep['condition_number']:.3e}  "
              f"after Jacobi={srep['condition_number']:.3e}")

    # The planning question: does the determinant already say everything? On a
    # 2x2 it does. Here the third minor and the smallest eigenvalue live on
    # different scales, which is exactly the gap inertia and moments would fill.
    smallest = min(float(r["eigenvalues"][0]) for r in rows)
    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil 3x3 even parity block, ENG-008 preparation",
        "work_order": "ATLAS-RH-ENG-007",
        "content_kind": "WEIL_PILOT_CONDITIONING_PREVIEW",
        "claim_scope": P.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "evidence_class": "E3",
        "rigorous": False,
        "hard_constraints_certified": False,
        "psd_claim": False,
        "status": "PREVIEW",
        "mpmath_used": True,
        "eigenvalue_solver_used": True,
        "basis": list(P.EVEN_BASIS),
        "basis_description": "one = 1, b = x(L-x), b2 = x^2 (L-x)^2; all even about x = L/2",
        "domain": {"cell": [repr(lo), repr(hi)], "label": list(P.CELL_LABEL)},
        "T": args.T,
        "dps": args.dps,
        "grid_points": args.points,
        "rows": rows,
        "summary": {
            "inertias_seen": sorted(list(i) for i in inertias),
            "inertia_constant_on_grid": len(inertias) == 1,
            "min_third_leading_minor": repr(minor3_min),
            "min_eigenvalue": repr(smallest),
            "max_condition_number": repr(cond_max),
            "max_condition_number_after_jacobi": repr(rescaled_max),
        },
        "eng008_notes": [
            "The Jacobi rescaling D M D is a congruence, so it changes no part of "
            "the inertia -- that is AtlasRH.posIndexAtLeast_congruence_iff and "
            "AtlasRH.rank_congruence, both proved. Conditioning is therefore free "
            "to fix and the fix costs no information.",
            "b2 is quadratic in L, unlike every basis element ENG-005 uses. "
            "pole.py's _laplace_d2L drops the second integral because d^2_L h = 0 "
            "for a linear-in-L element; that simplification is invalid for b2, so "
            "any E2''-style curvature argument on this block needs that machinery "
            "extended first.",
            "The third leading minor and the smallest eigenvalue sit on very "
            "different scales here, which is the gap a 2x2 block does not have: "
            "on a 2x2 the determinant fixes the spectrum given the trace, and "
            "ENG-006 found inertia and moments adding nothing as a result.",
            "This is E3. It is a plan input, never a warrant.",
        ],
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(CERT_NAME, body)
    print(f"\nwrote {path}")
    print(f"  inertias seen on grid: {sorted(list(i) for i in inertias)}")
    print(f"  min third leading minor: {minor3_min:.6e}")
    print(f"  condition number: {cond_max:.3e} -> {rescaled_max:.3e} after Jacobi rescaling")
    print("  evidence class E3 — heuristic preview, promoted as nothing")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
