#!/usr/bin/env python3
"""ENG-005 §4/§5 — recover degree-1 and compact degree-2 E1 under Candidate A.

Emits ``certificates/e1_degree1_log3_log4.json`` (odd pivot ``O1``) and
``certificates/e1_degree2_compact_log3_log4.json`` (compact even determinant
``E2``), each with a uniform rigorous positive lower bound on ``[log 3, log 4]``.

Both are cutoff-free: the archimedean term uses the exact real-space form, so
there is no frequency truncation to bound away.

    python3 scripts/certify_cutoff_free_e1.py [--release] [--quick]

``--release`` lifts the WO-RH-17 quarantine on these two certificates, and only
after a PASS. Nothing else is released.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import e1_cutoff_free as CF  # noqa: E402
import promotion  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from interval_backend import FlintUnavailable  # noqa: E402

DEPENDENCIES = (
    "src/pole.py",
    "src/weil_entries.py",
    "src/archimedean_realspace.py",
    "src/e1_cutoff_free.py",
    "src/curvature_derivation.py",
    "src/normalization.py",
    "scripts/certify_cutoff_free_e1.py",
)

SPECS = {
    "degree1": {
        "file": "e1_degree1_log3_log4.json",
        "quantity": "O1",
        "program": "RH/Weil degree-1 odd pivot — Candidate-A recovery",
        "work_order": "ATLAS-RH-ENG-005 §4",
        "basis": "q1 = x - L/2",
        "prime_kernel": "K_q1q1(a;L) = (L-a)(L^2-2La-2a^2)/6",
        "claim": "O1(L) = G[q1,q1] >= certified_lower_bound > 0 on [log 3, log 4]",
        # O1 clears every starting box without subdivision, so a target here only
        # costs time; the bound is already ~1.5e-2.
        "target": 0.0,
        "initial_boxes": 24,
    },
    "degree2": {
        "file": "e1_degree2_compact_log3_log4.json",
        "quantity": "E2",
        "program": "RH/Weil compact degree-2 even determinant — Candidate-A recovery",
        "work_order": "ATLAS-RH-ENG-005 §5",
        "basis": "{1, b}, b = x(L-x)",
        "prime_kernel": ("K00 = 2(L-a); K0b = (L-a)^2(L+2a)/3; "
                         "Kbb = (L-a)^3(L^2+3La+a^2)/15"),
        "claim": "E2(L) = G00 Gbb - G0b^2 >= certified_lower_bound > 0 on [log 3, log 4]",
        # E2 runs ~4.6e-6 at the left endpoint, so ask the cover to refine until
        # it clears 2e-6 rather than merely clearing 0 -- a bound 30x below the
        # true minimum would be a weak certificate for a quantity this small.
        "target": 2e-6,
        "initial_boxes": 32,
    },
}


def build(name: str, *, quick: bool) -> dict:
    spec = SPECS[name]
    result = CF.certify_positive(
        spec["quantity"],
        initial_boxes=8 if quick else spec["initial_boxes"],
        target=0.0 if quick else spec["target"],
        max_depth=6 if quick else CF.DEFAULT_MAX_DEPTH,
    )
    passed = result.certified_lower_bound > 0.0
    identities = CF.parity_identities()
    identities_ok = all(r["D2_matches_E2_plus_L2_G00_O1"] and r["det_matches_O1_times_E2"]
                        for r in identities)

    body = {
        "certificate_version": "1.0",
        "program": spec["program"],
        "work_order": spec["work_order"],
        "evidence_class": "E1",
        "rigorous": True,
        "status": "PASS" if (passed and identities_ok) else "FAIL",
        "hard_constraints_certified": bool(passed and identities_ok),
        "promotion_state": promotion.PROMOTED_STATE if (passed and identities_ok) else "REFUSED",
        promotion.NORMALIZATION_ID_FIELD: promotion.active_normalization_id(strict=True),
        "claim_scope": CF.CLAIM_SCOPE,
        "rh_proof_claim": False,
        "domain": {
            "L_interval": list(CF.CELL_LABEL),
            "L_left": repr(math.log(3.0)),
            "L_right": repr(math.log(4.0)),
            "closed": True,
        },
        "basis": spec["basis"],
        "prime_kernel": spec["prime_kernel"],
        "pole_candidate": "A",
        "claim": spec["claim"],
        "certified_lower_bound": repr(result.certified_lower_bound),
        "cutoff_free": True,
        "archimedean_route": {
            "method": "exact real-space kernel transform (no frequency cutoff)",
            "formula": ("Ginf_ij(L) = (K(0)/2)h_+(0) + int_0^L [K(0)-K(u)]w(u)du "
                        "+ K(0)S(L),  w(u)=e^{-u/2}/(1-e^{-2u}), "
                        "S(L)=sum_n e^{-(2n+1/2)L}/(2n+1/2)"),
            "analytic_hypothesis": (
                "termwise transform of the digamma series; see "
                "src/curvature_derivation.py INTERCHANGE_HYPOTHESIS"
            ),
            "independent_cross_check": (
                "frequency-space panel integration in src/weil_entries.py agrees "
                "to 2e-16 on the fast-decaying entries and differs on the "
                "slow-decaying ones by exactly the expected T-truncation tail"
            ),
        },
        "subdivision_statistics": result.to_dict(),
        "parity_identities": {
            "checked": identities,
            "all_hold": identities_ok,
            "note": ("D2 = E2 + L^2 G00 O1 and det(G_deg<=2) = O1 E2 both follow "
                     "from the pole and prime blocks being parity block diagonal"),
        },
        "dependencies": {"source_hashes": promotion.source_hashes(DEPENDENCIES)},
        "backend": result.to_dict().get("precision_bits"),
        "precision_bits": result.precision_bits,
        "mpmath_used": False,
        "quick_mode": quick,
        "historical_values_note": (
            "Historical degree-1/degree-2 numbers were produced under the REJECTED "
            "Candidate-B pole and are regression evidence only, never a warrant."
        ),
    }
    if quick:
        body["promotion_state"] = "REFUSED"
        body["hard_constraints_certified"] = False
        body["status"] = "QUICK_SMOKE_TEST_NOT_PROMOTABLE"
    return body


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", action="store_true")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", choices=sorted(SPECS), default=None)
    args = ap.parse_args()

    names = [args.only] if args.only else list(SPECS)
    rc = 0
    for name in names:
        try:
            body = build(name, quick=args.quick)
        except FlintUnavailable as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            print("Cutoff-free E1 requires python-flint/Arb (ENG-005).", file=sys.stderr)
            return 2
        except ValueError as exc:
            print(f"ERROR ({name}): {exc}", file=sys.stderr)
            return 1

        passed = body["status"] == "PASS"
        if args.release and not passed:
            print(f"ERROR: refusing to release {name} without a PASS", file=sys.stderr)
            return 1

        if args.release:
            import normalization as N

            body["quarantine_released"] = {
                "by": f"{SPECS[name]['work_order']} Candidate-A recovery",
                "previous_state": N.QUARANTINE_STATE,
                "authorised_because": "rigorous cutoff-free Candidate-A regeneration passed",
            }
            path = write_certificate(SPECS[name]["file"], body,
                                     allow_quarantine_change=True)
        else:
            path = write_certificate(SPECS[name]["file"], body)

        print(f"wrote {path}")
        print(f"  {name}: status={body['status']} bound={body['certified_lower_bound']} "
              f"promotion_state={body['promotion_state']}")
        stats = body["subdivision_statistics"]
        print(f"  boxes={stats['boxes_examined']} depth={stats['max_subdivision_depth']} "
              f"parity_identities={body['parity_identities']['all_hold']}")
        rc = rc or (0 if passed else 1)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
