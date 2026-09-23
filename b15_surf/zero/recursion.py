"""Route B — factorization recursion for Tr(phi^3) tree amplitudes.

Seeds (transcribed): A_3 = 1,  A_4 = 1/X_{1,3} + 1/X_{2,4}
(arXiv:2401.05483 §I; arXiv:2312.16282 §2).

Recursion: for the polygon with cyclically ordered vertex list V = (v_0,...,v_m),
the triangle containing the edge (v_0, v_m) has apex v_k, 0 < k < m; the two
sub-polygons (v_0..v_k) and (v_k..v_m) are joined across the chords (v_0,v_k)
and (v_k,v_m) (each contributing 1/X when it is a genuine chord, i.e. not an
edge). This is the propagator-factorization form of the amplitude, organised
as a memoised recursion on vertex tuples — an implementation independent of
the brute-force chord enumeration of Route A.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import Callable, Tuple

from .kinematics import Point, X


def amplitude_B(pt: Point, n: int) -> Fraction:
    Xf = lambda i, j: X(pt, n, i, j)
    return _amp(tuple(range(1, n + 1)), Xf)


def amplitude_B_general(Xf: Callable[[int, int], Fraction], n: int) -> Fraction:
    """Route B over an arbitrary accessor X(i,j) (used by delta_shift / splits)."""
    return _amp(tuple(range(1, n + 1)), Xf)


def _amp(V: Tuple[int, ...], Xf) -> Fraction:
    @lru_cache(maxsize=None)
    def rec(vs: Tuple[int, ...]) -> Fraction:
        m = len(vs)
        if m <= 3:
            return Fraction(1)          # A_3 = 1 (and degenerate 2-gons)
        total = Fraction(0)
        v0, vm = vs[0], vs[-1]
        for k in range(1, m - 1):
            vk = vs[k]
            left = vs[: k + 1]
            right = vs[k:]
            term = rec(left) * rec(right)
            if len(left) >= 3:
                term /= Xf(v0, vk)
            if len(right) >= 3:
                term /= Xf(vk, vm)
            total += term
        return total
    return rec(V)
