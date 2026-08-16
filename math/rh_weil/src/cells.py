"""Prime-power cell splitting for the scalar Weil verifier (WO-RH-02).

Stdlib-only. No RH claim: cell geometry and jump bookkeeping only.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import log
from typing import Iterable, List, Sequence, Tuple


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def prime_powers_in_log_interval(
    L_min: float, L_max: float
) -> List[Tuple[int, int, float]]:
    """Return (q, p, log q) for prime powers with log q in [L_min, L_max]."""
    if L_max < L_min:
        raise ValueError("L_max < L_min")
    # q <= exp(L_max); use a safe integer ceiling.
    from math import ceil, exp

    c = max(2, int(ceil(exp(L_max) + 1e-12)))
    out: List[Tuple[int, int, float]] = []
    for p in range(2, c + 1):
        if not is_prime(p):
            continue
        q = p
        while q <= c:
            lq = log(q)
            if L_min - 1e-15 <= lq <= L_max + 1e-15:
                out.append((q, p, lq))
            if q > c // p:
                break
            q *= p
    out.sort(key=lambda t: t[2])
    return out


@dataclass(frozen=True)
class Cell:
    L_left: float
    L_right: float
    left_break: Tuple[int, int, float] | None  # (q,p,logq) or None
    right_break: Tuple[int, int, float] | None

    @property
    def open_interval(self) -> Tuple[float, float]:
        return (self.L_left, self.L_right)


def split_cells(L_min: float, L_max: float) -> List[Cell]:
    """Split [L_min, L_max] at prime-power log breakpoints."""
    breaks = prime_powers_in_log_interval(L_min, L_max)
    # Include endpoints even if not breakpoints.
    points = [L_min]
    tagged: List[Tuple[float, Tuple[int, int, float] | None]] = [(L_min, None)]
    for q, p, lq in breaks:
        if lq <= L_min + 1e-15 or lq >= L_max - 1e-15:
            # Endpoint breaks still annotate the adjacent cell edges.
            continue
        tagged.append((lq, (q, p, lq)))
        points.append(lq)
    tagged.append((L_max, None))
    points.append(L_max)

    # Rebuild with endpoint annotations.
    end_breaks = {lq: (q, p, lq) for q, p, lq in breaks}
    coords = sorted(set(points))
    cells: List[Cell] = []
    for i in range(len(coords) - 1):
        left, right = coords[i], coords[i + 1]
        if right - left <= 1e-15:
            continue
        lb = end_breaks.get(left)
        rb = end_breaks.get(right)
        # Match floats carefully
        if lb is None:
            for q, p, lq in breaks:
                if abs(lq - left) <= 1e-12:
                    lb = (q, p, lq)
                    break
        if rb is None:
            for q, p, lq in breaks:
                if abs(lq - right) <= 1e-12:
                    rb = (q, p, lq)
                    break
        cells.append(Cell(left, right, lb, rb))
    return cells


def current_research_cell() -> Cell:
    """The imported notebook cell L ∈ [log 3, log 4]."""
    return Cell(log(3), log(4), (3, 3, log(3)), (4, 2, log(4)))
