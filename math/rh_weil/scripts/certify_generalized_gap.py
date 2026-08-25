#!/usr/bin/env python3
"""ATLAS-RH-ENG-009 §WO-RH-58 — reference metric and generalized gap certificates.

    python3 scripts/certify_generalized_gap.py [--quick]

Produces, in order:

  ``e0_eng009_reference_metric.json``            E0 — exact PD proof of M
  ``e1_eng009_generalized_gap_log3_log4.json``   E1 — per-block gap enclosures

The measured object is ``lambda_min(G, M)(L)`` for the pencil of each certified
cutoff-free block against the exact L^2 reference metric. Lower bounds are
uniform over the cell by shifted positivity (adaptive interval covers of the
leading minors of ``G - lam M``, exactly preconditioned); upper bounds are
certified Rayleigh quotients of rational witness vectors at the scouted
bottleneck. The scout (float Sylvester bisection) is E3 and recorded as such;
no eigensolver appears anywhere, including the scout.

The T=84 finite-cutoff family is deliberately absent: §Anti-overclaim is
explicit that finite-T and cutoff-free families are not one scaling sequence.

No RH proof claim is made. Claim scope is ``finite_dimensional_weil_compression``.
"""
from __future__ import annotations

import argparse
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

import archimedean_realspace as AR  # noqa: E402
import even3  # noqa: E402
import generalized_gap as GG  # noqa: E402
import normalization as N  # noqa: E402
import promotion  # noqa: E402
import reference_metric as RM  # noqa: E402
import weil_entries as WE  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from content_kinds import KIND_GENERALIZED_GAP  # noqa: E402
from interval_backend import interval_box, require_flint  # noqa: E402
from interval_cover import NotSeparated, adaptive_cover  # noqa: E402

METRIC_FILE = "e0_eng009_reference_metric.json"
GAP_FILE = "e1_eng009_generalized_gap_log3_log4.json"

PRECISION_BITS = 160

#: The cutoff-free blocks the pencil is certified for, with the preconditioner
#: applied to the *shifted* matrix. The even3 exponents are the ENG-008 frozen
#: ones; the small blocks need none.
BLOCKS: Tuple[Dict[str, Any], ...] = (
    {"name": "scalar", "basis": ("one",), "exponents": None},
    {"name": "degree1_odd", "basis": ("q1",), "exponents": None},
    {"name": "degree2_even", "basis": ("one", "b"), "exponents": None},
    {"name": "degree3_odd", "basis": ("q1", "b3"), "exponents": None},
    # even3 stores exponents for d_i = 2^(-e_i); GG.precondition multiplies by
    # 2^(+e), so the sign flips here to reproduce the same scale-up congruence.
    {"name": "degree4_even3", "basis": ("one", "b", "b2"),
     "exponents": tuple(-e for e in even3.PRECONDITIONER_EXPONENTS)},
)

DEPENDENCIES = (
    "src/pole.py",
    "src/core.py",
    "src/basis_algebra.py",
    "src/reference_metric.py",
    "src/generalized_gap.py",
    "src/weil_entries.py",
    "src/archimedean_realspace.py",
    "src/even3.py",
    "src/interval_cover.py",
    "src/interval_backend.py",
    "src/normalization.py",
    "scripts/certify_generalized_gap.py",
)


def cell_bounds() -> Tuple[float, float]:
    return math.log(3.0), math.log(4.0)


def assemble_entries(basis: Sequence[str], box: Any) -> Dict[Tuple[str, str], Any]:
    """The block's Gram entries over an L-box, centred mean-value form."""
    _, arb, acb, _ = require_flint()
    mid = float(box.mid()) if hasattr(box, "mid") else float(box)
    primes = WE.prime_powers_below(mid)
    out: Dict[Tuple[str, str], Any] = {}
    for a, i in enumerate(basis):
        for b, j in enumerate(basis):
            if b < a:
                continue
            val, _rec = AR.gram_entry_centred(i, j, box, arb, acb,
                                              prime_powers=primes)
            out[(i, j)] = val
    return out


# --------------------------------------------------------------------------- #
# Stage 1: the reference metric, E0                                            #
# --------------------------------------------------------------------------- #
def stage_metric() -> Dict[str, Any]:
    print("\n=== E0 reference metric (§WO-RH-58) ===")
    records = {}
    for spec in BLOCKS:
        records[spec["name"]] = RM.certify_positive_definite(spec["basis"])
        minors = records[spec["name"]]["unit_leading_minors"]
        print(f"  [ok] {spec['name']}: M(1) minors {minors}")
    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil generalized gap — Candidate A",
        "work_order": "WO-RH-58",
        "claim_scope": "finite_dimensional_weil_compression",
        "content_kind": "WEIL_GENERALIZED_GAP_CERTIFICATE",
        "role": "reference_metric",
        "evidence_class": "E0",
        "numeric_warrant": "E0",
        "rigorous": True,
        "hard_constraints_certified": True,
        "psd_claim": False,
        "status": "PASS",
        "mpmath_used": False,
        "arithmetic": "exact_rational",
        "statement": (
            "The L^2 reference metric M_ij(L) = int_0^L h_i h_j dx is positive "
            "definite on every listed block for every L > 0, by exact Sylvester "
            "minors of M(1) and the diagonal congruence M(L) = D^T M(1) D, "
            "D = diag(L^(d_i + 1/2))."),
        "blocks": records,
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(METRIC_FILE, body)
    print(f"wrote {path}")
    return body


# --------------------------------------------------------------------------- #
# Stage 2: per-block gap enclosures, E1                                        #
# --------------------------------------------------------------------------- #
def dyadic_floor(x: float, bits: int = 24) -> Fraction:
    return Fraction(math.floor(x * 2 ** bits), 2 ** bits)


def scout_block(spec: Dict[str, Any], cell: Tuple[float, float],
                n_grid: int) -> Dict[str, Any]:
    """E3: locate min_L lambda_min(G, M)(L) on a grid. Never a warrant."""
    _, arb, _, _ = require_flint()
    basis = spec["basis"]
    lo, hi = cell
    points = [lo + (hi - lo) * k / (n_grid - 1) for k in range(n_grid)]
    lams = []
    for L in points:
        entries = assemble_entries(basis, arb(repr(L)))
        floats = {k: float(v.mid()) for k, v in entries.items()}
        lams.append(GG.scout_gap_at(basis, floats, L))
    k_min = min(range(len(points)), key=lambda k: lams[k])
    return {
        "grid": [repr(p) for p in points],
        "lambda_star": [repr(v) for v in lams],
        "min_lambda": lams[k_min],
        "argmin_L": points[k_min],
        "evidence_class": "E3",
        "warrant_note": "scouting only; the covers below are the warrant",
    }


def certify_lower(spec: Dict[str, Any], lam: Fraction,
                  cell: Tuple[float, float], initial_boxes: int) -> List[Dict[str, Any]]:
    """E1: adaptive covers proving every leading minor of G - lam M positive."""
    basis = spec["basis"]
    n = len(basis)
    covers = []
    for k in range(1, n + 1):
        def ev(blo, bhi, _k=k):
            box = interval_box(blo, bhi)
            minors = GG.shifted_minors_over(
                basis, lambda b: assemble_entries(basis, b), lam, box, box,
                exponents=spec["exponents"])
            m = minors[_k - 1]
            return float(m.lower()), float(m.upper())

        try:
            cov = adaptive_cover(ev, quantity=f"minor{k}(G - lam*M)", cell=cell,
                                 target=0.0, initial_boxes=initial_boxes,
                                 max_depth=12)
        except NotSeparated as stop:
            covers.append({"minor": k, "status": "NOT_SEPARATED",
                           "detail": str(stop)})
            print(f"    minor {k}: NOT_SEPARATED — {stop}")
            break
        covers.append({
            "minor": k,
            "status": "PASS",
            "certified_lower_bound_preconditioned": repr(cov.certified_lower_bound),
            "boxes": cov.boxes_examined, "max_depth": cov.max_depth,
        })
        print(f"    minor {k}: PASS bound={cov.certified_lower_bound:.6e} "
              f"boxes={cov.boxes_examined} depth={cov.max_depth}")
    return covers


def certify_upper(spec: Dict[str, Any], scout: Dict[str, Any]) -> Dict[str, Any]:
    """E1: a certified Rayleigh quotient at the scouted bottleneck point."""
    _, arb, _, _ = require_flint()
    basis = spec["basis"]
    L = scout["argmin_L"]
    box = arb(repr(L))
    entries = assemble_entries(basis, box)
    floats = {k: float(v.mid()) for k, v in entries.items()}
    if len(basis) == 1:
        v = [Fraction(1)]
    else:
        v = GG.scout_min_eigvec(basis, floats, scout["min_lambda"], L)
    quotient = GG.rayleigh_upper(basis, entries, v, box)
    return {
        "at_L": repr(L),
        "witness_vector": [str(x) for x in v],
        "rayleigh_enclosure": [repr(float(quotient.lower())),
                               repr(float(quotient.upper()))],
        "certified_upper_bound": repr(float(quotient.upper())),
    }


def stage_gap(cell: Tuple[float, float], *, quick: bool) -> Dict[str, Any]:
    print("\n=== E1 generalized gap enclosures (§WO-RH-58) ===")
    n_grid = 9 if quick else 17
    initial_boxes = 32 if quick else 64
    blocks = []
    all_pass = True
    for spec in BLOCKS:
        t0 = time.time()
        print(f"  block {spec['name']} {tuple(spec['basis'])}")
        scout = scout_block(spec, cell, n_grid)
        lam = dyadic_floor(0.9 * scout["min_lambda"])
        print(f"    scout min lambda* = {scout['min_lambda']:.8e} at "
              f"L = {scout['argmin_L']:.6f}; certifying lam = {float(lam):.8e}")
        covers = certify_lower(spec, lam, cell, initial_boxes)
        ok = all(c["status"] == "PASS" for c in covers) and len(covers) == len(spec["basis"])
        upper = certify_upper(spec, scout)
        all_pass = all_pass and ok
        blocks.append({
            "block": spec["name"],
            "basis": list(spec["basis"]),
            "dimension": len(spec["basis"]),
            "reference_metric": "l2_gram_on_support",
            "preconditioner_exponents": (list(spec["exponents"])
                                         if spec["exponents"] else None),
            "certified_lambda_lower_uniform": str(lam) if ok else None,
            "certified_lambda_lower_float": repr(float(lam)) if ok else None,
            "status": "PASS" if ok else "FAIL",
            "shifted_minor_covers": covers,
            "upper_bound_at_bottleneck": upper,
            "scout": scout,
            "elapsed_seconds": repr(time.time() - t0),
        })
        print(f"    upper bound {upper['certified_upper_bound']} at "
              f"L={upper['at_L']}  [{time.time() - t0:.1f}s]")
    body = {
        "certificate_version": "0.1",
        "program": "RH/Weil generalized gap — Candidate A",
        "work_order": "WO-RH-58",
        "claim_scope": "finite_dimensional_weil_compression",
        "content_kind": KIND_GENERALIZED_GAP,
        "role": "gap_enclosures",
        "evidence_class": "E1",
        "numeric_warrant": "E1",
        "logical_implication_warrant": (
            "FORMAL: manifest theorems generalized_rayleigh (shifted positivity "
            "is the Rayleigh bound), generalized_pencil_congruence (the bound "
            "is invariant under simultaneous congruence), and "
            "preconditioned_gap_certificate3 (the composed 3x3 replay)"),
        "rigorous": True,
        "hard_constraints_certified": all_pass,
        "psd_claim": False,
        "status": "PASS" if all_pass else "PARTIAL",
        "mpmath_used": False,
        "backend": "python-flint / Arb",
        "precision_bits": PRECISION_BITS,
        "statement": (
            "For each block, lambda_min(G, M)(L) >= the stated lower bound for "
            "every L in [log 3, log 4], by positivity of all leading minors of "
            "G - lam M under adaptive interval covers; and lambda_min at the "
            "scouted bottleneck is <= the stated Rayleigh upper bound. M is the "
            "exact L^2 reference metric of e0_eng009_reference_metric.json."),
        "preconditioner_convention": (
            "exponents e mean the congruence D (G - lam M) D with "
            "D = diag(2^e); exact powers of two, minor signs unchanged"),
        "invariance_note": (
            "generalized eigenvalues of (G, M) are unchanged under simultaneous "
            "congruence G -> S^T G S, M -> S^T M S; raw eigenvalues and raw "
            "determinants are not"),
        "domain": {"cell": [repr(cell[0]), repr(cell[1])],
                   "label": ["log 3", "log 4"]},
        "blocks": blocks,
        "reference_metric_certificate": METRIC_FILE,
        "normalization_certificate_id": N.normalization_id(),
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
    }
    path = write_certificate(GAP_FILE, body)
    print(f"wrote {path}")
    return body


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--release", action="store_true")
    args = ap.parse_args()
    require_flint()
    cell = cell_bounds()
    stage_metric()
    gap = stage_gap(cell, quick=args.quick)
    print("\n=== summary ===")
    for blk in gap["blocks"]:
        lo = blk["certified_lambda_lower_float"]
        up = blk["upper_bound_at_bottleneck"]["certified_upper_bound"]
        print(f"  {blk['block']:16s} lambda_min in [{lo}, {up}]  ({blk['status']})")
    return 0 if gap["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
