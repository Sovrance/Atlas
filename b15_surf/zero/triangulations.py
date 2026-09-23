"""Route A — Tr(phi^3) tree amplitude as a sum over triangulations.

A_n = sum_{T} prod_{chords (i,j) in T} 1 / X_{i,j}          (arXiv:2312.16282 §2;
                                                              e.g. eq. (2.6) of
                                                              arXiv:2405.09608 for n=5,
                                                              eq. (2.3) for n=6)
Triangulations are enumerated by brute force as the non-crossing (n-3)-subsets of
chords; their number must be the Catalan number C_{n-2} (5, 14, 42, 132 for
n = 5, 6, 7, 8). This route is deliberately independent of the recursion in
``recursion.py`` (Route B).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import List, Tuple

from .kinematics import Chord, Point, chords


def crossing(a: Chord, b: Chord) -> bool:
    (i, j), (k, l) = a, b
    return (i < k < j < l) or (k < i < l < j)


def triangulations(n: int) -> List[Tuple[Chord, ...]]:
    ch = chords(n)
    out = []
    for sub in combinations(ch, n - 3):
        ok = True
        for x in range(len(sub)):
            for y in range(x + 1, len(sub)):
                if crossing(sub[x], sub[y]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            out.append(tuple(sub))
    return out


def catalan(m: int) -> int:
    from math import comb
    return comb(2 * m, m) // (m + 1)


def amplitude_A(pt: Point, n: int) -> Fraction:
    """Route A: A_n = sum_T prod 1/X."""
    total = Fraction(0)
    for T in triangulations(n):
        term = Fraction(1)
        for chord in T:
            term /= pt[chord]
        total += term
    return total


def amplitude_A_shifted(pt_shifted, n: int) -> Fraction:
    """Same sum with an arbitrary (i,j) -> value accessor (used by delta_shift)."""
    total = Fraction(0)
    for T in triangulations(n):
        term = Fraction(1)
        for chord in T:
            term /= pt_shifted[chord]
        total += term
    return total
