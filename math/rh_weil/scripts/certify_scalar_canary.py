#!/usr/bin/env python3
"""ENG-004 §5 — regenerate the scalar E1 canary under Candidate A.

Emits ``certificates/e1_scalar_log3_log4.json`` with a **uniform** rigorous
positive lower bound over ``[log 3, log 4]``, bound to the active normalization
id and to the hashes of the sources it was produced from.

This is the one certificate ENG-004 authorises to leave quarantine, and only
after a successful rigorous regeneration: the release is explicit
(``--release``), it happens after the bound is certified, and every other
disputed E1 stays quarantined.

    python3 scripts/certify_scalar_canary.py [--release] [--T 200000] [--quick]

``--quick`` shrinks the grid and cutoff for smoke-testing; it refuses to write a
promotable certificate.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import promotion  # noqa: E402
import scalar as _scalar  # noqa: E402
import scalar_canary as SC  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from interval_backend import FlintUnavailable  # noqa: E402

#: Sources whose bytes the certificate is bound to. A change to any of them
#: invalidates the certificate through ``promotion.stale_dependencies``.
DEPENDENCIES = (
    "src/pole.py",
    "src/scalar_canary.py",
    "src/scalar.py",
    "src/normalization.py",
    "src/interval_backend.py",
    "scripts/certify_scalar_canary.py",
)

#: Historical notebook figure. Regression evidence only — never an acceptance
#: constant (ENG-004 §5). The gate is ``certified_lower_bound > 0``.
NOTEBOOK_REGRESSION_TARGET = "0.0753795566117244"


def build(*, T: int, precision_bits: int, quick: bool) -> dict:
    a, b = math.log(3.0), math.log(4.0)
    grid = SC.default_grid(a, b, coarse=5, refine=2) if quick else None
    result = SC.certify_scalar_canary(T=T, precision_bits=precision_bits, grid=grid)

    bound = result.certified_lower_bound
    passed = bound > 0.0
    structural = _scalar.verify_scalar_cell().to_dict()

    # Regression review (§5): does the recovered geometry still contain the
    # historical notebook minimum? Under the rejected Candidate B the entry never
    # dropped below ~0.1276 on this cell, so the notebook figure was unreachable.
    target = float(NOTEBOOK_REGRESSION_TARGET)
    containing = [
        {"L": repr(gp.L), "enclosure": [repr(gp.lower), repr(gp.upper)]}
        for gp in result.grid
        if gp.lower <= target <= gp.upper
    ]

    body = {
        "certificate_version": "1.0",
        "program": "RH/Weil scalar cell — Candidate-A E1 canary",
        "work_order": "ATLAS-RH-ENG-004",
        "evidence_class": "E1",
        "rigorous": True,
        "status": "PASS" if passed else "FAIL_SCALAR_NOT_SEPARATED_FROM_ZERO",
        "hard_constraints_certified": bool(passed),
        "promotion_state": promotion.PROMOTED_STATE if passed else "REFUSED",
        promotion.NORMALIZATION_ID_FIELD: promotion.active_normalization_id(strict=True),
        "claim_scope": SC.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "domain": {
            "L_interval": list(SC.CELL_LABEL),
            "L_left": repr(a),
            "L_right": repr(b),
            "closed": True,
            "prime_power_breakpoints": "log 3 and log 4 are the cell endpoints; none interior",
        },
        "claim": (
            "G00(L) = G0 - Gp + Ginf >= certified_lower_bound > 0 for every L in "
            "[log 3, log 4], with G0 the adjudicated Candidate-A pole."
        ),
        "certified_lower_bound": repr(bound),
        "anchor_L": repr(result.anchor_L),
        "pole_candidate": "A",
        "pole_formula": SC.pole.POLE_FORMULA,
        "convexity_certificate": result.convexity,
        "tail_lemma": {
            "statement": "0 <= R_T(L) <= (4/pi)(h_+(T) + kappa(T))/T for every L > 0",
            "upper_bound": repr(result.tail_upper),
            "sign": "non-negative: h_+ increasing and h_+(T) > 0, so dropping the "
                    "tail can only lower the certified bound",
            **result.lemma_A,
        },
        "structural_report": structural,
        "grid": [gp.to_dict() for gp in result.grid],
        "tangent_certificates": result.tangent_bounds,
        "subdivision_statistics": result.stats,
        "backend": result.stats["backend"],
        "precision_bits": result.precision_bits,
        "T": result.T,
        "mpmath_used": False,
        "regression_review": {
            "notebook_target_G00": NOTEBOOK_REGRESSION_TARGET,
            "status": "CONTAINED" if containing else "NOT_CONTAINED",
            "containing_grid_enclosures": containing,
            "observed_min_lower_on_grid": repr(result.stats["observed_min_lower_on_grid"]),
            "observed_argmin_on_grid": repr(result.stats["observed_argmin_on_grid"]),
            "note": (
                "Regression evidence only, never an acceptance constant. Recorded "
                "because the rejected Candidate B kept this entry above ~0.1276 on "
                "the cell, which put the historical minimum out of reach; Candidate "
                "A recovers it."
            ),
        },
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
        "quick_mode": quick,
    }
    if quick:
        body["promotion_state"] = "REFUSED"
        body["hard_constraints_certified"] = False
        body["status"] = "QUICK_SMOKE_TEST_NOT_PROMOTABLE"
    return body


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=int, default=SC.DEFAULT_T)
    ap.add_argument("--precision-bits", type=int, default=SC.DEFAULT_PRECISION_BITS)
    ap.add_argument("--quick", action="store_true",
                    help="coarse grid for smoke tests; never promotable")
    ap.add_argument("--release", action="store_true",
                    help="lift the WO-RH-17 quarantine on the scalar certificate "
                         "after a successful rigorous regeneration (ENG-004 §4)")
    args = ap.parse_args()

    try:
        body = build(T=args.T, precision_bits=args.precision_bits, quick=args.quick)
    except FlintUnavailable as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        print("The scalar E1 canary requires python-flint/Arb; there is no mpmath "
              "fallback that may emit E1 (ENG-004 §5).", file=sys.stderr)
        return 2

    passed = body["status"] == "PASS"
    if args.release and not passed:
        print("ERROR: refusing to release the quarantine without a PASS", file=sys.stderr)
        return 1

    if args.release:
        # Explicit, authorised release of *this* certificate only (§4). The
        # writer's quarantine guard is bypassed for exactly this write.
        import normalization as N

        body["quarantine_released"] = {
            "by": "ATLAS-RH-ENG-004 §4 scalar canary",
            "previous_state": N.QUARANTINE_STATE,
            "authorised_because": "rigorous Candidate-A regeneration passed",
        }
        path = write_certificate("e1_scalar_log3_log4.json", body,
                                 allow_quarantine_change=True)
    else:
        path = write_certificate("e1_scalar_log3_log4.json", body)

    print(f"wrote {path}")
    print(f"status={body['status']}  bound={body['certified_lower_bound']}  "
          f"promotion_state={body['promotion_state']}")
    print(f"regression review: {body['regression_review']['status']} "
          f"(notebook {NOTEBOOK_REGRESSION_TARGET})")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
