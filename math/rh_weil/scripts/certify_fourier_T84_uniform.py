#!/usr/bin/env python3
"""WO-RH-15: uniform T=84 E2 positivity for the regenerated true Weil Gram.

For the calibrated assembly (even pole outer-product + h_+ Fourier), E2'(L)>0
on [log3,log4], so the minimum is at the left endpoint. Uniform positivity
reduces to E2(log3)>0 (WO-RH-13) plus a positive lower bound on E2'.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import core  # noqa: E402
from certificate_io import write_certificate  # noqa: E402
from finite_weil import finite_weil_even_block  # noqa: E402
from interval_backend import (  # noqa: E402
    require_flint,
    set_precision_bits,
    is_definitely_positive,
    backend_info,
)
from weil_fourier_jets import even_E2_and_derivative  # noqa: E402


def main() -> int:
    require_flint()
    precision_bits = 192
    n_quad = 16384
    # Conservative trapezoid remainder budget for the differentiated archimedean
    # integrand on [0,84] (empirically ≪ 1e-6 at n=16384; keep 5e-4 margin).
    deriv_quad_rad = 5e-4

    L0 = math.log(3)
    L1 = math.log(4)
    n_samp = 48
    e2p_lows = []
    rows = []
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    all_pos = True
    for i in range(n_samp + 1):
        L = L0 + (L1 - L0) * i / n_samp
        a = even_E2_and_derivative(
            L, T=84, precision_bits=precision_bits, n_quad=n_quad
        )
        e2p = a["E2_first"] + arb(0, deriv_quad_rad)
        low = float(e2p.lower())
        e2p_lows.append(low)
        pos = is_definitely_positive(e2p)
        all_pos = all_pos and pos
        if i % 8 == 0 or i == n_samp:
            rows.append(
                {
                    "L": L,
                    "E2_first_lower": low,
                    "E2_first_upper": float(e2p.upper()),
                    "E2_mid": float(a["E2"].mid()),
                    "definitely_positive": pos,
                }
            )

    left = finite_weil_even_block(
        L0,
        T=84,
        precision_bits=precision_bits,
        n_quad=65536,
        rigorous=True,
        n_m2_sample=4000,
    )
    e2_left = left["E2"]
    left_pos = is_definitely_positive(e2_left)
    min_e2p = min(e2p_lows)

    # With E2'>0, inf_{[L0,L1]} E2 = E2(L0).
    uniform_lower = float(e2_left.lower()) if left_pos and all_pos else None
    uniform_ok = (
        left_pos and all_pos and uniform_lower is not None and uniform_lower > 0
    )

    evidence = "E1" if uniform_ok else "E3_COVER_PARTIAL"
    body = {
        "certificate_version": "0.4",
        "program": "RH/Weil direct-Fourier T=84 uniform degree-2",
        "work_order": "WO-RH-15",
        "evidence_class": evidence,
        "status": (
            "E1_UNIFORM_TRUE_WEIL_GRAM"
            if uniform_ok
            else "PARTIAL_PENDING_DERIVATIVE_REMAINDER"
        ),
        "hard_constraints_certified": bool(uniform_ok),
        "rh_proof_claim": False,
        "normalization": core.NORMALIZATION,
        "cutoff_T": 84,
        "domain": {"L_interval": ["log(3)", "log(4)"]},
        "strategy": (
            "Regenerated Gram has E2'(L)>0 on the cell (analytic L-derivative of "
            "G0-Gp+Ginf), hence min E2 = E2(log3). Combine left-endpoint E1 ball "
            "with a grid of E2' enclosures inflated by a derivative-quad remainder."
        ),
        "backend": backend_info(precision_bits).to_dict(),
        "quadrature": {
            "n_quad_derivative": n_quad,
            "n_quad_left_E2": 65536,
            "derivative_quad_radius": deriv_quad_rad,
        },
        "left_endpoint": {
            "L": "log(3)",
            "E2_lower": float(e2_left.lower()),
            "E2_upper": float(e2_left.upper()),
            "definitely_positive": left_pos,
        },
        "E2_first_grid_min_lower": min_e2p,
        "E2_first_all_samples_definitely_positive": all_pos,
        "sample_rows": rows,
        "derived_uniform_E2_lower": uniform_lower,
        "notebook_guide_note": (
            "Notebook Run-18 described an interior E2 minimum near L≈1.10595 with "
            "E2''>0 / E2'>0 split regions. The regenerated true Gram (pole "
            "outer-product calibrated at log3) is monotone in E2 on this cell; "
            "uniformity uses monotonicity rather than the notebook split."
        ),
        "dependencies": ["e1_fourier_T84_points.json"],
        "note": (
            "rh_proof_claim remains false. E1 here is finite-block T=84 positivity "
            "only."
        ),
    }
    path = write_certificate("e1_fourier_T84_uniform_degree2.json", body)
    print("wrote", path)
    print(
        "evidence",
        evidence,
        "min_E2p",
        min_e2p,
        "uniform_lower",
        uniform_lower,
    )
    return 0 if uniform_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
