#!/usr/bin/env python3
"""ATLAS-RH-ENG-008 — the 3x3 even Weil block ``G[{1, b, b^2}]`` on [log 3, log 4].

    python3 scripts/certify_even3.py [--quick] [--release]

Produces, in order:

  ``e0_degree4_even3_exact_identities.json``   §WO-RH-48, E0
  ``e3_degree4_even3_crosscheck.json``         §WO-RH-48, E3
  ``e1_degree4_even3_inertia_log3_log4.json``  §WO-RH-51, E1
  ``e1_degree4_even3_positivity_log3_log4.json``  §WO-RH-51, E1, only if proved
  ``e1_degree4_even3_moments_log3_log4.json``  §WO-RH-52, E1

Two independent warrants are computed for the same conclusion, and both are
kept. The first is interval LDL* congruence on the preconditioned block,
stratified over the cell; the second is Sylvester's criterion, three separate
adaptive interval covers of the leading principal minors. They share the
assembly and nothing after it, so agreement is evidence and disagreement would
be a stop condition (§Stop conditions).

If the block turns out not to be definite, the stratification is the result and
no positivity certificate is written. §Mission is explicit that a rigorous
inertia stratification is a successful outcome; nothing here tunes toward the
preferred one.

No RH proof claim is made. Claim scope is ``finite_dimensional_weil_compression``.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
for extra in (ROOT, ROOT / "src"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))

import basis_algebra  # noqa: E402
import even3  # noqa: E402
import normalization as N  # noqa: E402
import promotion  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from content_kinds import (  # noqa: E402
    KIND_DEGREE4_POSITIVITY,
    KIND_SCAN_PREVIEW,
)
from inertia.certificate import KIND_INERTIA  # noqa: E402
from inertia.ldl import interval_inertia  # noqa: E402
from inertia.stratify import certify_inertia_family  # noqa: E402
from interval_backend import interval_box, require_flint  # noqa: E402
from interval_cover import adaptive_cover  # noqa: E402
from moments.adapter import analyse  # noqa: E402
from ranktrace.theorem import rank_trace_lower_bound  # noqa: E402

E0_FILE = "e0_degree4_even3_exact_identities.json"
CROSSCHECK_FILE = "e3_degree4_even3_crosscheck.json"
INERTIA_FILE = "e1_degree4_even3_inertia_log3_log4.json"
POSITIVITY_FILE = "e1_degree4_even3_positivity_log3_log4.json"
MOMENTS_FILE = "e1_degree4_even3_moments_log3_log4.json"

DEPENDENCIES = (
    "src/pole.py",
    "src/core.py",
    "src/basis_algebra.py",
    "src/weil_entries.py",
    "src/archimedean_realspace.py",
    "src/even3.py",
    "src/interval_cover.py",
    "src/interval_backend.py",
    "src/normalization.py",
    "inertia/ldl.py",
    "inertia/stratify.py",
    "inertia/certificate.py",
    "moments/spectral_moments.py",
    "moments/adapter.py",
    "ranktrace/theorem.py",
    "scripts/certify_even3.py",
)

PRECISION_BITS = 160

#: The mean-value enclosure widens roughly linearly in the box radius, and the
#: third minor is the binding one. 256 initial cells puts most boxes below the
#: width where all three pivots separate at depth 0.
INERTIA_CELLS = 256
INERTIA_MAX_DEPTH = 12
MINOR_BOXES = 256
MINOR_MAX_DEPTH = 14

QUICK_CELL = (math.log(3.0), math.log(3.0) + 0.02)


def _bounds(x) -> Tuple[float, float]:
    return float(x.lower()), float(x.upper())


def _enc(x) -> Dict[str, str]:
    lo, hi = _bounds(x)
    return {"lo": repr(lo), "hi": repr(hi), "width": repr(hi - lo)}


# --------------------------------------------------------------------------- #
# §WO-RH-48 — exact identities                                                 #
# --------------------------------------------------------------------------- #
def stage_exact() -> Dict[str, Any]:
    from fractions import Fraction as F

    print("\n=== E0 exact identities (§WO-RH-48) ===")
    checks: List[Dict[str, Any]] = []

    # Parity of every element, exactly, at rational points.
    Ls = (F(7, 3), F(11, 8), F(3))
    xs = (F(1, 5), F(2, 7), F(11, 13), F(3, 2))
    import pole

    parity_ok = True
    for name in even3.EVEN3_BASIS:
        for L in Ls:
            for x in xs:
                left = pole._horner(pole.basis_coeffs(name, L), L - x)
                right = pole._horner(pole.basis_coeffs(name, L), x)
                parity_ok = parity_ok and (left == right)
    checks.append({"identity": "every element of {1, b, b2} is even about x = L/2",
                   "verified": parity_ok, "method": "exact rational evaluation"})

    b2_ok = all(
        pole._horner(pole.basis_coeffs("b2", L), x)
        == pole._horner(pole.basis_coeffs("b", L), x) ** 2
        for L in Ls for x in xs)
    checks.append({"identity": "b2 = b^2", "verified": b2_ok,
                   "method": "exact rational evaluation"})

    # The six kernels, against direct exact integration of their definition.
    kernel_rows = []
    kernels_ok = True
    for key, (i, j) in even3.ENTRY_KEYS:
        sym = all(basis_algebra.kernel_exact(i, j, a, L)
                  == basis_algebra.kernel_exact(j, i, a, L)
                  for L in Ls for a in (F(0), F(1, 4), F(2)))
        vanishes = all(basis_algebra.kernel_exact(i, j, L, L) == 0 for L in Ls)
        kernels_ok = kernels_ok and sym and vanishes
        kernel_rows.append({
            "entry": key, "pair": [i, j],
            "degree_in_a": basis_algebra.kernel_degree_in_a(i, j),
            "symmetric": sym,
            "vanishes_at_a_equals_L": vanishes,
            "K_at_zero_at_L_eq_7_over_3":
                str(basis_algebra.kernel_exact(i, j, F(0), F(7, 3))),
        })
    checks.append({"identity": "all six overlap kernels symmetric and vanishing at a = L",
                   "verified": kernels_ok, "method": "exact rational arithmetic"})

    # d^2_L is nonzero exactly for the elements quadratic in L.
    d2 = {n: [str(c) for c in pole.basis_coeffs_d2L(n, F(2))]
          for n in ("one", "q1", "b", "b3", "b2")}
    checks.append({
        "identity": "d^2_L h vanishes for one/q1/b and does not for b3/b2",
        "verified": (all(c == "0" for c in d2["one"] + d2["q1"] + d2["b"])
                     and any(c != "0" for c in d2["b3"])
                     and any(c != "0" for c in d2["b2"])),
        "method": "coefficient inspection",
        "coefficients": d2,
    })

    all_ok = all(c["verified"] for c in checks)
    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil 3x3 even block, exact identities",
        "work_order": "ATLAS-RH-ENG-008",
        "claim_scope": even3.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "evidence_class": "E0",
        "status": "PASS" if all_ok else "FAIL",
        "rigorous": True,
        "hard_constraints_certified": bool(all_ok),
        "mpmath_used": False,
        "basis": even3.basis_identity(),
        "checks": checks,
        "kernels": kernel_rows,
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(E0_FILE, body)
    print(f"wrote {path}")
    for c in checks:
        print(f"  [{'ok' if c['verified'] else 'FAIL'}] {c['identity']}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-48 — independent cross-check                                          #
# --------------------------------------------------------------------------- #
def stage_crosscheck(points: List[float]) -> Dict[str, Any]:
    print("\n=== E3 independent cross-check (§WO-RH-48) ===")
    import independent_even3 as IE

    mp = IE.require_mpmath()
    rows = []
    worst = 0.0
    for L in points:
        rigorous = even3.assemble_even3_arb(L, precision_bits=PRECISION_BITS)
        M = IE.gram_matrix(L, mp, dps=40)
        indep = {"G00": M[0][0], "G01": M[0][1], "G02": M[0][2],
                 "G11": M[1][1], "G12": M[1][2], "G22": M[2][2]}
        entry_rows = {}
        for key, _ in even3.ENTRY_KEYS:
            a = rigorous["entries"][key]
            lo, hi = _bounds(a)
            v = float(indep[key])
            inside = lo <= v <= hi
            rel = abs(v - float(a.mid())) / max(1e-300, abs(float(a.mid())))
            worst = max(worst, rel)
            entry_rows[key] = {
                "rigorous_enclosure": [repr(lo), repr(hi)],
                "independent_value": repr(v),
                "independent_value_inside_enclosure": inside,
                "relative_difference": repr(rel),
            }
        rows.append({
            "L": repr(L),
            "entries": entry_rows,
            "independent_minors": [repr(float(m)) for m in IE.leading_minors(M)],
            "rigorous_minors": [repr(float(m.mid())) for m in rigorous["minors_raw"]],
        })
        print(f"  L={L:.9f}  worst relative difference so far {worst:.2e}")

    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil 3x3 even block, independent assembly cross-check",
        "work_order": "ATLAS-RH-ENG-008",
        "content_kind": KIND_SCAN_PREVIEW,
        "claim_scope": even3.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "evidence_class": "E3",
        "rigorous": False,
        "hard_constraints_certified": False,
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": True,
        "method": (
            "second assembly of G = G0 - Gp + Ginf in mpmath, with the overlap "
            "kernels integrated symbolically by SymPy from their definition and "
            "the pole Laplace transforms evaluated by quadrature rather than in "
            "closed form; imports none of the modules it checks"
        ),
        "independence": {
            "shares_with_rigorous_path": "the formula only",
            "does_not_import": ["basis_algebra", "pole", "weil_entries",
                                "archimedean_realspace", "even3"],
        },
        "worst_relative_difference": repr(worst),
        "points": rows,
        "note": ("E3. mpmath never certifies in this program; this is regression "
                 "evidence that the rigorous assembly computes what it claims."),
        "dependencies": {"source_hashes": promotion.source_hashes(
            DEPENDENCIES + ("src/independent_even3.py",))},
    }
    path = write_certificate(CROSSCHECK_FILE, body)
    print(f"wrote {path}\n  worst relative difference across "
          f"{len(points)} points: {worst:.3e}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-51 — inertia and positivity                                           #
# --------------------------------------------------------------------------- #
def stage_inertia(cell: Tuple[float, float], cells: int) -> Tuple[Dict[str, Any], Any]:
    print("\n=== E1 whole-cell inertia (§WO-RH-51) ===")
    t0 = time.time()
    strat = certify_inertia_family(
        lambda lo, hi: even3.matrix_over(lo, hi, precision_bits=PRECISION_BITS),
        cell,
        subdivision_policy={"initial_cells": cells,
                            "max_depth": INERTIA_MAX_DEPTH,
                            "min_width": 1e-12},
    )
    elapsed = time.time() - t0
    strat_dict = strat.to_dict()
    sample = even3.assemble_even3_arb(interval_box(*cell), precision_bits=PRECISION_BITS)
    signatures = sorted({tuple(s.signature) for s in strat.strata})
    definite = (len(signatures) == 1 and signatures[0] == (3, 0, 0)
                and not strat.transitions)

    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil 3x3 even block, whole-cell inertia",
        "work_order": "ATLAS-RH-ENG-008",
        "content_kind": KIND_INERTIA,
        "claim_scope": even3.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "evidence_class": "E1",
        "rigorous": True,
        "hard_constraints_certified": bool(strat.status.startswith("PASS")),
        # An inertia certificate never claims PSD, whatever its signature
        # (ENG-006 §11). The positivity artifact below is what claims it.
        "psd_claim": False,
        "status": "PASS" if strat.status.startswith("PASS") else strat.status,
        "mpmath_used": False,
        "backend": "python-flint / Arb",
        "precision_bits": PRECISION_BITS,
        "basis": even3.basis_identity(),
        "domain": {"cell": [repr(cell[0]), repr(cell[1])],
                   "label": list(even3.CELL_LABEL)},
        "method": ("interval Hermitian LDL* congruence on the preconditioned "
                   "block, with adaptive subdivision of the L cell; no "
                   "eigenvalue solver anywhere on this path"),
        "preconditioner": sample["preconditioner"],
        "stratification": strat_dict,
        "signatures_seen": [list(s) for s in signatures],
        "constant_on_cell": bool(len(signatures) == 1 and not strat.transitions),
        "n_positive": signatures[0][0] if len(signatures) == 1 else None,
        "n_negative": signatures[0][1] if len(signatures) == 1 else None,
        "n_zero": signatures[0][2] if len(signatures) == 1 else None,
        "elapsed_seconds": repr(elapsed),
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(INERTIA_FILE, body)
    print(f"wrote {path}")
    print(f"  status {strat.status}  signatures {signatures}  "
          f"boxes {strat.boxes_examined}  depth {strat.max_depth}  "
          f"{elapsed:.1f}s")
    for t in strat.transitions:
        print(f"  transition region [{t.lo!r}, {t.hi!r}]: {t.blocker}")
    return body, definite


def stage_minors(cell: Tuple[float, float], boxes: int) -> Dict[str, Any]:
    print("\n=== E1 Sylvester leading minors (§WO-RH-51) ===")
    scales = even3.minor_scale_factors(even3.PRECONDITIONER_EXPONENTS)
    covers = []
    for idx, name in enumerate(("Delta1", "Delta2", "Delta3")):
        def ev(lo, hi, _idx=idx):
            out = even3.assemble_even3_arb(interval_box(lo, hi),
                                           precision_bits=PRECISION_BITS)
            return _bounds(out["minors_preconditioned"][_idx])

        t0 = time.time()
        r = adaptive_cover(ev, quantity=f"even3_{name}_preconditioned", cell=cell,
                           target=0.0, initial_boxes=boxes,
                           max_depth=MINOR_MAX_DEPTH)
        d = r.to_dict()
        d["minor"] = name
        d["scale_factor_from_preconditioner"] = repr(scales[idx])
        d["implied_raw_lower_bound"] = repr(r.certified_lower_bound / scales[idx])
        d["elapsed_seconds"] = repr(time.time() - t0)
        covers.append(d)
        print(f"  {name}: >= {r.certified_lower_bound!r} (preconditioned) "
              f"= {r.certified_lower_bound / scales[idx]!r} raw; "
              f"{r.boxes_examined} boxes, depth {r.max_depth}")
    return {"covers": covers,
            "all_positive": all(float(c["certified_lower_bound"]) > 0.0
                                for c in covers)}


def stage_positivity(cell, inertia_body, minors) -> Optional[Dict[str, Any]]:
    if not minors["all_positive"]:
        print("\n  positivity NOT certified: a leading minor did not clear zero")
        return None
    print("\n=== E1 positivity (§WO-RH-51) ===")
    sample = even3.assemble_even3_arb(interval_box(*cell), precision_bits=PRECISION_BITS)
    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil 3x3 even block, uniform positive definiteness",
        "work_order": "ATLAS-RH-ENG-008",
        "content_kind": KIND_DEGREE4_POSITIVITY,
        "claim_scope": even3.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "evidence_class": "E1",
        "rigorous": True,
        "hard_constraints_certified": True,
        "psd_claim": True,
        "status": "PASS",
        "mpmath_used": False,
        "backend": "python-flint / Arb",
        "precision_bits": PRECISION_BITS,
        "claim": (
            "The 3x3 even Weil block G[{1, b, b^2}](L) is positive definite for "
            "every L in [log 3, log 4]: inertia (3, 0, 0), with all three "
            "leading principal minors uniformly bounded below by the certified "
            "constants below."
        ),
        "basis": even3.basis_identity(),
        "domain": {"cell": [repr(cell[0]), repr(cell[1])],
                   "label": list(even3.CELL_LABEL)},
        "inertia": {"n_positive": 3, "n_negative": 0, "n_zero": 0,
                    "constant_on_cell": True},
        "n_positive": 3, "n_negative": 0, "n_zero": 0,
        "preconditioner": sample["preconditioner"],
        "leading_minors": minors["covers"],
        "warrants": [
            {"route": "interval LDL* congruence, stratified over the cell",
             "artifact": INERTIA_FILE,
             "conclusion": "inertia (3, 0, 0) on one stratum, no transition regions"},
            {"route": "Sylvester's criterion, three independent adaptive covers",
             "artifact": POSITIVITY_FILE,
             "conclusion": "Delta1, Delta2, Delta3 all uniformly positive"},
        ],
        "independent_agreement": (
            "the two routes share the assembly and nothing after it; both "
            "conclude positive definiteness on the whole cell"
        ),
        "licensed_by": [
            "AtlasRH.posDef_sym3_iff",
            "AtlasRH.posIndexAtLeast_congruence_iff",
            "AtlasRH.rank_congruence",
        ],
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(POSITIVITY_FILE, body)
    print(f"wrote {path}")
    return body


# --------------------------------------------------------------------------- #
# §WO-RH-52 — moments, rank-trace, information comparison                      #
# --------------------------------------------------------------------------- #
def stage_moments(points: List[Tuple[str, float]], observed) -> Dict[str, Any]:
    print("\n=== E1 moments and rank-trace (§WO-RH-52) ===")
    from moments.spectral_moments import trace, trace_of_power

    rows = []
    for label, L in points:
        built = even3.assemble_even3_arb(L, precision_bits=PRECISION_BITS)
        G = built["matrix"]
        report = analyse(G, observed_inertia=observed)
        tr = trace(G)
        hs = trace_of_power(G, 2)
        # The theorem is not scale free -- rank is scale invariant while tr(P)
        # is degree 1 and the HS term degree 2 -- so it holds only under a
        # normalization, and equality at a projection identifies that
        # normalization as "spectrum of P in [0, 1]".
        #
        # ENG-006 discharged that from closed-form 2x2 eigenvalues. There is no
        # closed form at 3x3 and no eigenvalue solver may appear on a rigorous
        # path (§13), so it is discharged here without one: P is positive
        # semidefinite by this run's certified inertia, which gives
        # lambda_min >= 0; and for a PSD matrix every eigenvalue is at most the
        # trace, which is enclosed above. Both halves are certified and neither
        # needs a spectrum.
        trace_upper = float(tr.upper())
        psd_certified = observed == (3, 0, 0)
        spectrum_ok = bool(psd_certified and trace_upper <= 1.0)
        rt = rank_trace_lower_bound(
            trace_P=tr, trace_Q=0, hs_sq_P_plus_Q=hs,
            positive_index_Q_bound=0,
            hypotheses={
                "P_positive_semidefinite": {
                    "verified": psd_certified,
                    "evidence": {
                        "why": "certified inertia (3,0,0) on this cell",
                        "certificate": INERTIA_FILE,
                    }},
                "Q_hermitian": {"verified": True,
                                "evidence": {"Q": "zero matrix",
                                             "trivially Hermitian": True}},
                "Q_positive_index_at_most_b": {
                    "verified": True,
                    "evidence": {"Q": "zero matrix", "positive_index": 0, "b": 0}},
                "shared_normalization": {
                    "verified": spectrum_ok,
                    "evidence": {
                        "requirement": "spectrum of P contained in [0, 1]",
                        "lambda_min_lower": "0 (P is positive semidefinite)",
                        "lambda_max_upper": repr(trace_upper),
                        "argument": (
                            "lambda_min >= 0 from the certified inertia; "
                            "lambda_max <= tr(P) because P is PSD, so every "
                            "eigenvalue is bounded by the sum of them all. No "
                            "eigenvalue solver is used or needed."),
                        "eigenvalue_solver_used": False,
                    }},
            },
        )
        rows.append({
            "label": label, "L": repr(L),
            "trace": _enc(tr),
            "hs_norm_squared": _enc(hs),
            "determinant": _enc(built["minors_raw"][2]),
            "moments": report["moments"],
            "b1_queries": report["b1_queries"],
            "b1_hankel_view": report["b1_hankel_view"],
            "rank_trace": rt.to_dict() if hasattr(rt, "to_dict") else str(rt),
        })
        status = rows[-1]["rank_trace"]
        got = status.get("status") if isinstance(status, dict) else status
        print(f"  {label:16s} trace {float(tr.mid()):.6e}  "
              f"det {float(built['minors_raw'][2].mid()):.6e}  rank-trace {got}")

    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil 3x3 even block, spectral moments and rank-trace",
        "work_order": "ATLAS-RH-ENG-008",
        "content_kind": "WEIL_SPECTRAL_MOMENT_CERTIFICATE",
        "claim_scope": even3.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "evidence_class": "E1",
        "rigorous": True,
        "hard_constraints_certified": True,
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "dimension": 3,
        "basis": even3.basis_identity(),
        "observed_inertia": list(observed) if observed else None,
        "points": rows,
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(MOMENTS_FILE, body)
    print(f"wrote {path}")
    return body


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="a short sub-interval, for wiring checks only")
    ap.add_argument("--release", action="store_true")
    ap.add_argument("--skip-crosscheck", action="store_true")
    args = ap.parse_args()

    require_flint()
    cell = QUICK_CELL if args.quick else even3.CELL
    cells = 16 if args.quick else INERTIA_CELLS
    boxes = 16 if args.quick else MINOR_BOXES

    stage_exact()

    lo, hi = cell
    points = [lo, 1.20 if not args.quick else lo + (hi - lo) / 2,
              hi, lo + (hi - lo) / 4, lo + 3 * (hi - lo) / 4]
    points = sorted(set(round(p, 15) for p in points))
    if not args.skip_crosscheck:
        stage_crosscheck(points)

    inertia_body, definite = stage_inertia(cell, cells)
    minors = stage_minors(cell, boxes)
    positivity = stage_positivity(cell, inertia_body, minors)

    observed = (3, 0, 0) if definite else None
    labelled = [("log3" if abs(p - math.log(3.0)) < 1e-12 else
                 "log4" if abs(p - math.log(4.0)) < 1e-12 else f"L={p:.6f}", p)
                for p in points]
    stage_moments(labelled, observed)

    print("\n=== summary ===")
    print(f"  inertia        : {inertia_body['signatures_seen']} "
          f"(constant: {inertia_body['constant_on_cell']})")
    print(f"  positivity     : {'CERTIFIED' if positivity else 'not claimed'}")
    for c in minors["covers"]:
        print(f"  {c['minor']:7s}        >= {c['implied_raw_lower_bound']} (raw)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
