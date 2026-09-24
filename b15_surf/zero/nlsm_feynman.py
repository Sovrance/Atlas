"""Independent NLSM flavour-ordered tree amplitudes from Feynman rules.

Lagrangian (transcribed, arXiv:2312.16282 eq. (7.4)):
    L_NLSM = (1/(8 lambda^2)) Tr(d_mu U^dag d^mu U),   U = (I + lambda Phi)(I - lambda Phi)^{-1}
(Cayley parametrisation). For U unitary, Phi is anti-Hermitian; writing
Phi = i pi with pi Hermitian and using dU = 2 lambda (1-lambda Phi)^{-1} dPhi (1-lambda Phi)^{-1},
    Tr(dU^dag dU) = 4 lambda^2 Tr[ dpi (1 + lambda^2 pi^2)^{-1} dpi (1 + lambda^2 pi^2)^{-1} ],
so, with lambda = 1 (the paper's "units with f_pi = 1", 2401.05483 eq. 2) and the
mostly-plus signature under which the letter's A^NLSM_4 = -(X13+X24) is
reproduced (declared stipulation S1 of prereg-002 §ZERO; the signature only
fixes the OVERALL sign of every even-point amplitude),
    L = -(1/2) sum_{a,b>=0} (-1)^{a+b} Tr[ dpi pi^{2a} dpi pi^{2b} ].

Flavour-ordered Feynman rules (Tr(T^a T^b) = delta^{ab}, iM = sum of diagrams,
d_mu -> -i p_mu incoming, propagator -i/P^2): the (2m+2)-point vertex is i*v with
    v(P_1..P_{2m+2}) = (1/2) (-1)^m sum_{a+b=m} sum_{r=0}^{2m+1} P_r . P_{r+2a+1}
(one term per cyclic placement of the two derivatives in Tr[dpi pi^{2a} dpi pi^{2b}]),
and the tree amplitude is  M = sum_{planar trees} prod v / prod P^2 .
Trees are summed by Berends-Giele flavour-ordered currents, exactly in Fraction;
all dot products come from the X's via 2 p_a.p_b = -c_{a,b} (eq. 2.3/2.4) and
P^2 of a consecutive block (i..j) is X_{i,j+1} (eq. 2.2).
Check: v_4 = -(s+t) = -(X13+X24) = A^NLSM_4  (matches 2401.05483 §I).
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import List, Tuple

from .kinematics import Point, X, dot

Block = Tuple[int, int]     # consecutive legs (i..j), 1-based, i <= j


def _gram(pt: Point, n: int) -> List[List[Fraction]]:
    return [[dot(pt, n, a, b) for b in range(1, n + 1)] for a in range(1, n + 1)]


def _block_dot(G, A: List[int], B: List[int]) -> Fraction:
    s = Fraction(0)
    for a in A:
        for b in B:
            s += G[a - 1][b - 1]
    return s


def vertex(G, legs: List[List[int]]) -> Fraction:
    """v for legs given as lists of external indices (momentum = sum of block;
    the off-shell leg is passed as the COMPLEMENT block, i.e. -P by momentum
    conservation is represented by the remaining external legs)."""
    m2 = len(legs)
    assert m2 % 2 == 0 and m2 >= 4
    mm = (m2 - 2) // 2
    total = Fraction(0)
    for a in range(0, mm + 1):
        step = 2 * a + 1
        for r in range(m2):
            total += _block_dot(G, legs[r], legs[(r + step) % m2])
    return Fraction(1, 2) * (-1) ** mm * total


def _odd_partitions(i: int, j: int):
    """Partitions of consecutive legs i..j into k >= 3 (odd) consecutive blocks."""
    def rec(start, parts):
        if start > j:
            if len(parts) >= 3 and len(parts) % 2 == 1:
                yield list(parts)
            return
        for end in range(start, j + 1):
            parts.append((start, end))
            yield from rec(end + 1, parts)
            parts.pop()
    yield from rec(i, [])


def nlsm_amplitude(pt: Point, n: int) -> Fraction:
    """Flavour-ordered NLSM tree amplitude A_n (n even) from Feynman rules."""
    if n % 2:
        raise ValueError("odd-point NLSM amplitudes vanish / do not exist")
    G = _gram(pt, n)
    all_legs = list(range(1, n + 1))

    @lru_cache(maxsize=None)
    def J(i: int, j: int) -> Fraction:
        """Off-shell current for legs i..j (amputated leg = -P)."""
        if i == j:
            return Fraction(1)
        P2 = X(pt, n, i, j + 1)
        total = Fraction(0)
        for parts in _odd_partitions(i, j):
            legs = [list(range(a, b + 1)) for (a, b) in parts]
            comp = [l for l in all_legs if l < i or l > j]
            v = vertex(G, legs + [comp])
            prod = v
            for (a, b) in parts:
                prod *= J(a, b)
            total += prod
        return total / P2

    total = Fraction(0)
    for parts in _odd_partitions(1, n - 1):
        legs = [list(range(a, b + 1)) for (a, b) in parts]
        v = vertex(G, legs + [[n]])
        prod = v
        for (a, b) in parts:
            prod *= J(a, b)
        total += prod
    return total
