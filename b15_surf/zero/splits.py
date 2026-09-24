"""Tree-level 2-splits of Tr(phi^3) amplitudes (transcribed).

Source: arXiv:2405.09608 (R37), §2.
  eq. (1.7)   A_S[X_{S = S1 (x) S2}] = A_{S1}(x) * A_{S2}(y)
  §2.3        two sub-surfaces overlapping on the triangle tau = (i, j, k):
              P = (i, i+1, ..., j, ..., k)  and  Q = (i, j, k, k+1, ..., i-1),
              index blocks B = (i..j-1), C = (j..k-1), A = (k..i-1)   [eq. 2.9]
  eq. (2.10)  X_{a,b} -> x_{k,b} + y_{a,i},   X_{a,c} -> x_{c,k} + y_{j,a}
              (a in A, b in B, c in C); chords living only in P (resp. Q) map
              trivially to their x (resp. y) variable.
  eq. (2.13)  on split kinematics  c_{a-1,c-1} = 0 and c_{b-1,c-1} = 0
              (the locus is a collection of vanishing mesh variables).
  eq. (2.7)   n=5 example, tau=(1,3,5) written with S1=(1,2,3,5), S2=(1,3,4,5):
              X_{1,3}->y_{1,3}, X_{1,4}->y_{1,3}+x_{1,4}, X_{2,4}->x_{1,4},
              X_{2,5}->y_{2,5}, X_{3,5}->x_{3,5};  eq. (2.8) A_5 -> (1/y13+1/y25)(1/x14+1/x35)
  eq. (2.1)   n=6 example (4-pt (x) 5-pt): X_{1,4} -> x_{1,3}+y_{1,4},
              X_{2,4}->y_{1,4}, X_{3,6}->y_{3,5}+x_{3,6}, X_{4,6}->x_{3,6}, rest trivial.
  §2.5        factorization near zeros (2312.16282 eq. 3.8) = sequential splits.

Transcription note (recorded, not repaired): the paper's letters x/y are used
for S2/S1 in §2.2 but the general formula (2.10) is written with x for the
(i..k) surface; we fix x := kinematics of P = (i..k) and y := kinematics of Q,
and VERIFY the fixed convention against the explicit examples (2.1) and (2.7)
in tests/test_b15_zero.py::T5 — that check is the transcription guard.

Near-zero factorization (2312.16282 eq. 3.3, transcribed): at n=5 with
c_{1,4} = 0,   A_5 -> (1/X_{1,3} + 1/X_{2,5}) (1/X_{1,4} + 1/X_{3,5}).
"""

from __future__ import annotations

import random
from fractions import Fraction
from typing import Dict, List, Tuple

from .kinematics import Chord, Point, X, c, chords, is_chord, norm, random_point


def surfaces(n: int, i: int, j: int, k: int) -> Tuple[List[int], List[int]]:
    """Vertex lists (cyclic order) of P = (i..k) and Q = (i, j, k, k+1, .., i-1)."""
    def arc(a, b):
        out = [a]
        while out[-1] != b:
            out.append(out[-1] % n + 1)
        return out
    P = arc(i, k)
    if j not in P[1:-1]:
        raise ValueError("j must lie strictly inside the arc i..k")
    Q = [i, j] + arc(k, (i - 2) % n + 1)
    return P, Q


def _sub_chords(V: List[int]) -> List[Chord]:
    """Chords of the sub-polygon with vertex list V, as pairs of ORIGINAL labels."""
    m = len(V)
    out = []
    for a in range(m):
        for b in range(a + 2, m):
            if a == 0 and b == m - 1:
                continue
            out.append((V[a], V[b]))
    return out


def split_point(n: int, i: int, j: int, k: int, x: Dict[Chord, Fraction],
                y: Dict[Chord, Fraction]) -> Point:
    """Build the big-surface point X from sub-surface kinematics (eq. 2.10).
    ``x`` is keyed by original-label pairs of chords of P, ``y`` of Q."""
    P, Q = surfaces(n, i, j, k)
    setP, setQ = set(P), set(Q)
    posP = {v: t for t, v in enumerate(P)}
    posQ = {v: t for t, v in enumerate(Q)}

    def xv(a, b):   # x_{a,b}: 0 unless (a,b) is a chord of P
        if a not in setP or b not in setP:
            raise KeyError((a, b))
        pa, pb = sorted((posP[a], posP[b]))
        if pb - pa <= 1 or (pa == 0 and pb == len(P) - 1):
            return Fraction(0)
        return x[(a, b)] if (a, b) in x else x[(b, a)]

    def yv(a, b):
        if a not in setQ or b not in setQ:
            raise KeyError((a, b))
        pa, pb = sorted((posQ[a], posQ[b]))
        if pb - pa <= 1 or (pa == 0 and pb == len(Q) - 1):
            return Fraction(0)
        return y[(a, b)] if (a, b) in y else y[(b, a)]

    B = P[: P.index(j)]          # i .. j-1
    C = P[P.index(j): -1]        # j .. k-1
    A = Q[2:]                    # k .. i-1   (contains k)
    setA, setBC = set(A), set(B) | set(C)
    pt: Point = {}
    for (u, v) in chords(n):
        if u in setA and v in setA:
            pt[(u, v)] = yv(u, v)                      # chord of Q only
        elif u in setBC and v in setBC:
            pt[(u, v)] = xv(u, v)                      # chord of P only
        else:
            a, other = (u, v) if u in setA else (v, u)
            if other in B:
                pt[(u, v)] = xv(k, other) + yv(a, i)   # eq. (2.10), first rule
            else:
                pt[(u, v)] = xv(other, k) + yv(j, a)   # eq. (2.10), second rule
    return pt


def sub_point(V: List[int], vals: Dict[Chord, Fraction]) -> Point:
    """Relabel a sub-surface's kinematics to a standard 1..m point."""
    m = len(V)
    pos = {v: t + 1 for t, v in enumerate(V)}
    out: Point = {}
    for (a, b), val in vals.items():
        out[norm(m, pos[a], pos[b])] = val
    return out


def random_split_kinematics(n: int, i: int, j: int, k: int, rng: random.Random):
    """Seeded (x, y) with every sub-chord nonzero; returns (X, x_point, y_point, |P|, |Q|)."""
    P, Q = surfaces(n, i, j, k)
    from .kinematics import random_rational
    while True:
        x = {ch: random_rational(rng) for ch in _sub_chords(P)}
        y = {ch: random_rational(rng) for ch in _sub_chords(Q)}
        Xp = split_point(n, i, j, k, x, y)
        if all(v != 0 for v in Xp.values()):      # x + y sums may vanish by chance
            return Xp, sub_point(P, x), sub_point(Q, y), len(P), len(Q)


def split_locus_pairs(n: int, i: int, j: int, k: int) -> List[Chord]:
    """Mesh variables that vanish IDENTICALLY on split kinematics.

    Computed exactly from the linear map (2.10): c_{a,b}(X(x,y)) is a linear
    form in (x, y); it is evaluated on every basis vector e_x, e_y and the pair
    is reported iff every coefficient is zero. This is a certificate, not a
    sample test. (Our literal reading of the index ranges in eq. (2.13) does
    NOT reproduce the explicit examples (2.5)/(2.8) — see edit-011 erratum
    E1; the registered locus is this exact derivation, which does reproduce
    (2.5): {c_{1,3}, c_{3,5}} for n=6, tau=(5,1,3).)"""
    P, Q = surfaces(n, i, j, k)
    basis = []
    for ch in _sub_chords(P):
        basis.append(({ch: Fraction(1)}, {}))
    for ch in _sub_chords(Q):
        basis.append(({}, {ch: Fraction(1)}))
    out = []
    for (a, b) in chords(n):
        vanishes = True
        for (x, y) in basis:
            xx = {ch: x.get(ch, Fraction(0)) for ch in _sub_chords(P)}
            yy = {ch: y.get(ch, Fraction(0)) for ch in _sub_chords(Q)}
            if c(split_point(n, i, j, k, xx, yy), n, a, b) != 0:
                vanishes = False
                break
        if vanishes:
            out.append((a, b))
    return out


def registered_splits(n: int) -> Dict[str, Tuple[int, int, int]]:
    """Registered split triangles per n (prereg-002 §ZERO)."""
    table = {
        5: {"S5-a": (3, 5, 1), "S5-b": (1, 3, 4)},
        6: {"S6-a": (5, 1, 3), "S6-b": (1, 3, 5), "S6-c": (2, 4, 6)},
        8: {"S8-a": (1, 3, 5), "S8-b": (1, 4, 7), "S8-c": (2, 5, 7)},
    }
    return table[n]
