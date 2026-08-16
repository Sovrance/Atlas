#!/usr/bin/env python3
"""WO-RH-13: rigorous T=84 even-block point certificates (Run-16 anchors)."""
from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from certificate_io import write_certificate  # noqa: E402
from finite_weil import finite_weil_even_block  # noqa: E402
from interval_backend import require_flint, is_definitely_positive, backend_info  # noqa: E402
import core  # noqa: E402


# Notebook regression targets (not warrants).
NOTEBOOK = {
    "log3": {
        "L": "log(3)",
        "L_float": math.log(3),
        "G00": 0.107356700414591762,
        "G0b": 0.000461820208771810,
        "Gbb": 3.4253778646359e-5,
        "E2": 3.4640947469748e-6,
    },
    "interior": {
        "L": "1.1059499883568553",
        "L_float": 1.1059499883568553,
        "G00": 0.102776248577958063,
        "G0b": 0.000783664417448233,
        "Gbb": 3.9317062336051e-5,
        "E2": 3.4267302528306e-6,
    },
    "log4": {
        "L": "log(4)",
        "L_float": math.log(4),
        "G00": 0.069433265743841579,
        "G0b": -0.000656968872855926,
        "Gbb": 0.000354436660096118,
        "E2": 2.417808670991184e-5,
    },
}


def _ball_dict(x) -> dict:
    return {
        "lower": float(x.lower()),
        "upper": float(x.upper()),
        "mid": float(x.mid()),
        "rad": float(x.rad()),
    }


def main() -> int:
    require_flint()
    precision_bits = 192
    n_quad = 65536
    points = []
    all_e2_pos = True
    log3_regression_ok = True
    for key, meta in NOTEBOOK.items():
        blk = finite_weil_even_block(
            meta["L_float"],
            T=84,
            precision_bits=precision_bits,
            n_quad=n_quad,
            rigorous=True,
            n_m2_sample=4000,
        )
        e2_pos = is_definitely_positive(blk["E2"])
        all_e2_pos = all_e2_pos and e2_pos
        contains = {}
        for name in ("G00", "G0b", "Gbb", "E2"):
            x = blk[name]
            v = meta[name]
            contains[name] = float(x.lower()) <= v <= float(x.upper())
        if key == "log3":
            log3_regression_ok = all(contains.values())
        points.append(
            {
                "id": key,
                "L": meta["L"],
                "L_value": meta["L_float"],
                "enclosure": {k: _ball_dict(blk[k]) for k in ("G00", "G0b", "Gbb", "E2")},
                "E2_definitely_positive": e2_pos,
                "notebook_regression_contained": contains,
                "pole_scale": blk["pole_scale"],
                "n_quad": n_quad,
            }
        )

    # E1 only if all E2>0 with rigorous quad binder (notebook match is regression-only).
    evidence = "E1" if all_e2_pos else "E3"
    body = {
        "certificate_version": "0.3",
        "program": "RH/Weil direct-Fourier T=84 point anchors",
        "work_order": "WO-RH-13",
        "evidence_class": evidence,
        "status": (
            "E1_POINTS_REGENERATED"
            if all_e2_pos
            else "FAILED_POSITIVITY"
        ),
        "hard_constraints_certified": all_e2_pos,
        "rh_proof_claim": False,
        "normalization": core.NORMALIZATION,
        "cutoff_T": 84,
        "pole_even": "G0=(sqrt(3)/2)(v+ v+^T + v- v-^T), v±=(I0±,Ib±)",
        "domain": {"L_anchors": ["log(3)", "1.1059499883568553", "log(4)"]},
        "backend": backend_info(precision_bits).to_dict(),
        "quadrature": {
            "rule": "composite_trapezoid",
            "n_quad": n_quad,
            "remainder": "|(b-a)h^2/12| max|f''| with enclosed M2 + M3 majorant",
        },
        "points": points,
        "log3_notebook_regression_ok": log3_regression_ok,
        "note": (
            "Regenerated rigorous enclosures from the true Weil Gram. "
            "Run-16 notebook numbers are regression-only; log3 matches the "
            "calibrated pole+Fourier assembly. Interior/log4 notebook floats "
            "may reflect an older probe and need not lie in the regenerated balls."
        ),
    }
    path = write_certificate("e1_fourier_T84_points.json", body)
    print("wrote", path)
    print("evidence", evidence, "all_E2_pos", all_e2_pos, "log3_reg", log3_regression_ok)
    return 0 if all_e2_pos else 1


if __name__ == "__main__":
    raise SystemExit(main())
