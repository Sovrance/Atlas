"""Scalar f0 verifier for localized Weil positivity (WO-RH-02).

Algebraic positivity of W00'' on [log 3, log 4] is E0.
Jump bookkeeping and uniqueness of an interior critical point are structural.
No RH claim.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from math import exp, log, sqrt
from typing import Any, Dict, List

import core
from cells import Cell, current_research_cell, split_cells


@dataclass(frozen=True)
class ScalarCellReport:
    cell: tuple[float, float]
    w00_second_positive: bool
    algebraic_reason: str
    left_jump: float | None
    right_jump: float | None
    at_most_one_interior_minimizer: bool
    evidence_class: str
    rh_proof_claim: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def w00_second_positive_on_r_interval(r_lo: float, r_hi: float) -> tuple[bool, str]:
    """Prove W00''(r) > 0 for r in [r_lo, r_hi] when r_lo >= plastic-ish bound.

    For r > 1: denominator √r (r²-1) > 0.
    Numerator r³ - r - 1 is increasing for r > 0 and positive for r >= 3
    (since 27 - 3 - 1 = 23 > 0). Hence on [3, 4] positivity is algebraic.
    """
    if r_lo <= 1:
        return False, "r_lo must exceed 1 for the cell formula denominator"
    # Evaluate numerator at left endpoint; monotonicity for r>1/sqrt(3).
    num_lo = r_lo**3 - r_lo - 1
    if num_lo <= 0:
        return False, f"numerator non-positive at r_lo={r_lo}"
    if r_hi < r_lo:
        return False, "empty interval"
    return True, (
        f"For r∈[{r_lo},{r_hi}], r³-r-1 ≥ {num_lo} > 0 and √r(r²-1)>0 "
        "⇒ W00''(L)=2(r³-r-1)/(√r(r²-1)) > 0 algebraically."
    )


def verify_scalar_cell(cell: Cell | None = None) -> ScalarCellReport:
    cell = cell or current_research_cell()
    r_lo, r_hi = exp(cell.L_left), exp(cell.L_right)
    ok, reason = w00_second_positive_on_r_interval(r_lo, r_hi)

    left_jump = None
    right_jump = None
    if cell.left_break is not None:
        q, p, _ = cell.left_break
        left_jump = -core.von_mangoldt_jump(q, p)
    if cell.right_break is not None:
        q, p, _ = cell.right_break
        right_jump = -core.von_mangoldt_jump(q, p)

    # W00'' > 0 ⇒ W00' strictly increasing on the open cell ⇒ ≤1 root.
    at_most_one = ok

    return ScalarCellReport(
        cell=(cell.L_left, cell.L_right),
        w00_second_positive=ok,
        algebraic_reason=reason,
        left_jump=left_jump,
        right_jump=right_jump,
        at_most_one_interior_minimizer=at_most_one,
        evidence_class="E0" if ok else "FAIL",
        rh_proof_claim=False,
    )


def sample_curvature_grid(cell: Cell | None = None, n: int = 21) -> List[float]:
    """Floating diagnostic samples (E3); not used for E0/E1 promotion."""
    cell = cell or current_research_cell()
    vals = []
    for i in range(n):
        t = i / (n - 1)
        L = cell.L_left * (1 - t) + cell.L_right * t
        # Stay slightly inside endpoints to avoid breakpoint ambiguity.
        if i == 0:
            L = cell.L_left + 1e-12
        if i == n - 1:
            L = cell.L_right - 1e-12
        vals.append(core.scalar_curvature(L))
    return vals
