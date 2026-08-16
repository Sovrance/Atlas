"""Compact real-space degree-2 even block (WO-RH-11)."""
from __future__ import annotations

from typing import Any, Dict

import core
from finite_weil import (
    gp_even_block,
    g0_even_block,
    pole_even_helpers,
    ginf_even_block_quad,
    archimedean_tail_even,
)
from interval_backend import require_flint, set_precision_bits, is_definitely_positive, backend_info
from weil_degree1 import assemble_O1


def assemble_E2_compact(
    L, precision_bits: int = 256, T_cap: int = 4000, n: int = 8000
) -> Dict[str, Any]:
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    L_a = arb(L)
    gp00, gp0b, gpbb = gp_even_block(L_a, arb)
    g000, g00b, g0bb = g0_even_block(L_a, arb)
    ebp, ebm = pole_even_helpers(L_a, arb)
    gi00, gi0b, gibb = ginf_even_block_quad(L_a, T_cap, arb, n=n)
    t00, t0b, tbb = archimedean_tail_even(L_a, float(T_cap), arb)
    # Enclose cutoff-free Ginf ∈ trunc ± tail (absolute majorant).
    gi00 = gi00 + arb(0, float(t00.upper()))
    gi0b = gi0b + arb(0, float(t0b.upper()))
    gibb = gibb + arb(0, float(tbb.upper()))
    G00 = g000 - gp00 + gi00
    G0b = g00b - gp0b + gi0b
    Gbb = g0bb - gpbb + gibb
    E2 = G00 * Gbb - G0b * G0b
    pole_det = g000 * g0bb - g00b * g00b
    return {
        "G00": G00,
        "G0b": G0b,
        "Gbb": Gbb,
        "E2": E2,
        "Eb_plus": ebp,
        "Eb_minus": ebm,
        "pole_det": pole_det,
        "normalization": core.NORMALIZATION,
        "evidence_class": "E3_PENDING_TIGHT_TAIL",
        "rh_proof_claim": False,
        "backend": backend_info(precision_bits).to_dict(),
        "notebook_regression_target": "4.07220283229438138e-6",
        "T_cap": T_cap,
    }


def certify_degree2_compact_e1(precision_bits: int = 256) -> Dict[str, Any]:
    import math

    samples = [math.log(3), 1.1059498113, 1.20, math.log(4)]
    rows = []
    all_pos = True
    min_lower = None
    pole_ok = True
    for L in samples:
        a = assemble_E2_compact(L, precision_bits=precision_bits)
        e2 = a["E2"]
        low = float(e2.lower())
        o1 = assemble_O1(L, precision_bits=precision_bits)["O1"]
        D2 = a["E2"] + (require_flint()[1](L) ** 2) * a["G00"] * o1
        full = o1 * a["E2"]
        # Pole det ≈ 0 (rank ≤ 1)
        pd = abs(float(a["pole_det"].mid()))
        pole_ok = pole_ok and pd < 1e-20
        rows.append(
            {
                "L": L,
                "E2_lower": low,
                "D2_lower": float(D2.lower()),
                "full_det_lower": float(full.lower()),
                "pole_det_abs": pd,
            }
        )
        all_pos = all_pos and is_definitely_positive(e2)
        min_lower = low if min_lower is None else min(min_lower, low)
    return {
        "certificate_version": "0.3",
        "program": "RH/Weil compact degree-2 even block",
        "work_order": "WO-RH-11",
        "evidence_class": "E3",
        "status": "POLE_WIRED_PENDING_TIGHT_TAIL_FOR_E1",
        "hard_constraints_certified": False,
        "rh_proof_claim": False,
        "normalization": core.NORMALIZATION,
        "domain": {"L_interval": ["log(3)", "log(4)"]},
        "samples": rows,
        "min_sample_E2_lower": min_lower,
        "all_samples_definitely_positive": all_pos,
        "pole_rank1_regression_ok": pole_ok,
        "backend": backend_info(precision_bits).to_dict(),
        "dependencies": ["e0_scalar_cell_log3_log4.json", "WO-RH-10 assembly"],
        "note": (
            "Even pole outer-product (√3/2)(v₊v₊ᵀ+v₋v₋ᵀ) wired. "
            "Cutoff-free archimedean uses truncated Fourier + absolute tail majorant; "
            "tighten tail / quad before E1 promotion."
        ),
    }
