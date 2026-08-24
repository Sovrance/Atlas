#!/usr/bin/env python3
"""WO-RH-18 — four-way normalization cross-check.

Compares the ``pole``, ``prime``, ``archimedean`` and ``total`` components across
four independent providers at the mandated points

    L in {log 3, 1.1059498113, 1.20, log 4},  basis {1, q1 = x-L/2, b = x(L-x)},

by **interval overlap**, and writes ``certificates/normalization_crosscheck.json``.

The run also records the legacy-vs-adopted pole audit: the rejected
``(sqrt(3)/2)`` even block against the adopted explicit-formula pole, showing the
``(sqrt(3)/2)cosh(L/2)`` discrepancy that vanishes only at ``L = log 3``.

No RH proof claim. Usage::

    python3 scripts/run_normalization_crosscheck.py [--T 84] [--arch-full] [--no-arch]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import cross_validation as CV  # noqa: E402
import normalization as N  # noqa: E402
import providers as P  # noqa: E402
from certificate_io import write_certificate  # noqa: E402

L_POINTS = [
    ("log3", math.log(3.0)),
    ("1.1059498113", 1.1059498113),
    ("1.20", 1.20),
    ("log4", math.log(4.0)),
]
BASIS = list(N.BASIS_NAMES)


def _pairs():
    for a_i, i in enumerate(BASIS):
        for j in BASIS[a_i:]:
            yield i, j


def run(T: float = 84.0, arch_full: bool = False, with_arch: bool = True) -> dict:
    provs = P.all_providers()
    rows = []
    disagreements = []

    for lname, L in L_POINTS:
        for i, j in _pairs():
            for component in ("pole", "prime", "archimedean", "total"):
                if component == "archimedean" and not with_arch:
                    continue
                # the expensive independent archimedean route: diagonal only unless --arch-full
                slow_ok = arch_full or (i == j)
                meas = {}
                for pr in provs:
                    if component in ("archimedean", "total") and not slow_ok and pr.name == "DirectFourierProvider":
                        meas[pr.name] = None
                        continue
                    try:
                        if component == "pole":
                            m = pr.pole_entry(i, j, L)
                        elif component == "prime":
                            m = pr.prime_entry(i, j, L)
                        elif component == "archimedean":
                            m = pr.arch_entry(i, j, L, T=T)
                        else:
                            m = pr.gram_entry(i, j, L, T=T) if with_arch else None
                    except Exception as exc:  # provider failure is reported, never hidden
                        m = None
                        disagreements.append(
                            {"L": lname, "entry": f"{i},{j}", "component": component,
                             "provider": pr.name, "error": f"{type(exc).__name__}: {exc}"}
                        )
                    meas[pr.name] = m
                pairs = CV.compare_all(meas)
                summary = CV.summarize(pairs)
                rows.append(
                    {
                        "L_label": lname,
                        "L": L,
                        "entry": [i, j],
                        "component": component,
                        "measurements": {k: (v.to_dict() if v else None) for k, v in meas.items()},
                        "pairs": [p.to_dict() for p in pairs],
                        "summary": summary,
                    }
                )
                if summary["status"] == CV.DISAGREE:
                    disagreements.append(
                        {"L": lname, "entry": f"{i},{j}", "component": component,
                         "detail": [p.to_dict() for p in pairs if p.status == CV.DISAGREE]}
                    )

    # legacy (rejected) vs adopted pole audit on the even sector
    legacy_audit = []
    for lname, L in L_POINTS:
        for i, j in (("one", "one"), ("one", "b"), ("b", "b")):
            adopted = N.pole_entry(i, j, L)
            legacy = N.legacy_pole_entry(i, j, L)
            ratio = legacy / adopted if adopted != 0 else None
            legacy_audit.append(
                {
                    "L_label": lname,
                    "L": L,
                    "entry": [i, j],
                    "adopted_candidate_A": adopted,
                    "legacy_candidate_B": legacy,
                    "ratio_B_over_A": ratio,
                    "predicted_ratio_sqrt3_over_2_cosh_L_over_2": N.legacy_over_adopted_ratio(L),
                    "relative_discrepancy": (ratio - 1.0) if ratio is not None else None,
                }
            )

    compared = [r for r in rows if r["summary"]["status"] != CV.UNAVAILABLE]
    agreeing = [r for r in compared if r["summary"]["status"] == CV.AGREE]
    return {
        "certificate_version": "0.1",
        "program": "RH/Weil normalization cross-check",
        "work_order": "WO-RH-18",
        "evidence_class": "E2_numeric_crosscheck",
        "status": "AGREE" if len(agreeing) == len(compared) and compared else "DISAGREE",
        "hard_constraints_certified": False,
        "rh_proof_claim": False,
        "normalization_id": N.normalization_id(),
        "normalization": N.normalization_content(),
        "providers": [
            {"name": p.name, "description": p.description,
             "external": bool(getattr(p, "external", False)),
             "available": bool(getattr(p, "available", lambda: True)())}
            for p in provs
        ],
        "settings": {"T": T, "arch_full": arch_full, "with_arch": with_arch,
                     "flint_available": P.FLINT},
        "coverage": {
            "L_points": [n for n, _ in L_POINTS],
            "basis": BASIS,
            "cells_total": len(rows),
            "cells_compared": len(compared),
            "cells_all_unavailable": len(rows) - len(compared),
            "note": "cells with a single available provider cannot be cross-checked and "
                    "are reported UNAVAILABLE rather than counted as agreement",
        },
        "results": rows,
        "legacy_pole_audit": legacy_audit,
        "disagreements": disagreements,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=84.0)
    ap.add_argument("--arch-full", action="store_true")
    ap.add_argument("--no-arch", action="store_true")
    args = ap.parse_args()

    body = run(T=args.T, arch_full=args.arch_full, with_arch=not args.no_arch)
    path = write_certificate("normalization_crosscheck.json", body)
    print(f"wrote {path}")
    print(f"status={body['status']}  compared={body['coverage']['cells_compared']}"
          f"/{body['coverage']['cells_total']}  disagreements={len(body['disagreements'])}")
    for row in body["legacy_pole_audit"]:
        if row["entry"] == ["one", "one"]:
            print(f"  legacy/adopted at L={row['L_label']:>14}: "
                  f"{row['ratio_B_over_A']:.12f}  (discrepancy {100*row['relative_discrepancy']:+.4f} %)")
    return 0 if body["status"] == "AGREE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
