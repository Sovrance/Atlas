"""Planar kinematics for n-point colour-ordered tree amplitudes (transcribed).

Source: arXiv:2312.16282 §2 (R37).
  eq. (2.2)  X_{i,j} = (p_i + ... + p_{j-1})^2            (planar variables)
  eq. (2.3)  c_{i,j} = -2 p_i . p_j                       (mesh / non-planar)
  eq. (2.4)  c_{i,j} = X_{i,j} + X_{i+1,j+1} - X_{i,j+1} - X_{i+1,j}
with X_{i,i+1} = p_i^2 = 0 (massless), X_{i,j} = X_{j,i} (Mobius symmetry of
the mesh, §2.1), and exactly n(n-3)/2 independent X's.

A kinematic *point* is a dict {(i,j): Fraction} over the chords of the n-gon
(1 <= i < j <= n, j - i >= 2, (i,j) != (1,n)). Everything is exact.
"""

from __future__ import annotations

import random
from fractions import Fraction
from typing import Dict, List, Tuple

Chord = Tuple[int, int]
Point = Dict[Chord, Fraction]


def chords(n: int) -> List[Chord]:
    """The n(n-3)/2 chords (i,j) of the n-gon in lexicographic order."""
    out = [(i, j) for i in range(1, n + 1) for j in range(i + 2, n + 1)
           if not (i == 1 and j == n)]
    assert len(out) == n * (n - 3) // 2
    return out


def norm(n: int, i: int, j: int) -> Chord:
    """Normalise a (possibly wrapped / reversed) label pair to a canonical
    (i,j) with 1 <= i < j <= n; boundary pairs return (i,i) style markers
    that :func:`X` maps to 0."""
    i = (i - 1) % n + 1
    j = (j - 1) % n + 1
    if i > j:
        i, j = j, i
    return (i, j)


def is_chord(n: int, i: int, j: int) -> bool:
    i, j = norm(n, i, j)
    if i == j:
        return False
    if j - i == 1 or (i == 1 and j == n):
        return False
    return True


def X(pt: Point, n: int, i: int, j: int) -> Fraction:
    """X_{i,j} with wrap-around and X_{i,i+1} = X_{i,i} = 0 (eq. 2.2 conventions)."""
    if not is_chord(n, i, j):
        return Fraction(0)
    return pt[norm(n, i, j)]


def c(pt: Point, n: int, i: int, j: int) -> Fraction:
    """Mesh variable c_{i,j} = X_{i,j} + X_{i+1,j+1} - X_{i,j+1} - X_{i+1,j}  (eq. 2.4).
    Valid for every pair (adjacent pairs give c_{i,i+1} = -X_{i,i+2} = -2 p_i.p_{i+1})."""
    return (X(pt, n, i, j) + X(pt, n, i + 1, j + 1)
            - X(pt, n, i, j + 1) - X(pt, n, i + 1, j))


def dot(pt: Point, n: int, a: int, b: int) -> Fraction:
    """p_a . p_b = -c_{a,b}/2 (eq. 2.3); p_a . p_a = 0 (massless)."""
    if (a - 1) % n == (b - 1) % n:
        return Fraction(0)
    return -c(pt, n, a, b) / 2


def c_rows(n: int, pairs: List[Chord]) -> List[List[Fraction]]:
    """Linear forms c_{a,b} as rows over the chord basis (eq. 2.4 is linear in X)."""
    idx = {ch: k for k, ch in enumerate(chords(n))}
    rows = []
    for (a, b) in pairs:
        row = [Fraction(0)] * len(idx)
        for (i, j, s) in ((a, b, 1), (a + 1, b + 1, 1), (a, b + 1, -1), (a + 1, b, -1)):
            if is_chord(n, i, j):
                row[idx[norm(n, i, j)]] += s
        rows.append(row)
    return rows


def random_rational(rng: random.Random, num_max: int = 60, den_max: int = 9) -> Fraction:
    """Nonzero rational with |numerator| <= num_max, 1 <= denominator <= den_max."""
    while True:
        p = rng.randint(-num_max, num_max)
        if p != 0:
            return Fraction(p, rng.randint(1, den_max))


def random_point(n: int, rng: random.Random) -> Point:
    """Seeded generic rational point (all X nonzero)."""
    return {ch: random_rational(rng) for ch in chords(n)}


def point_from_vector(n: int, v: List[Fraction]) -> Point:
    return {ch: Fraction(v[k]) for k, ch in enumerate(chords(n))}


def vector_from_point(n: int, pt: Point) -> List[Fraction]:
    return [pt[ch] for ch in chords(n)]


def relabel(pt: Point, n: int, perm: Dict[int, int]) -> Point:
    """Apply a vertex relabelling i -> perm[i] to a point: X'_{perm(i),perm(j)} = X_{i,j}."""
    out: Point = {}
    for (i, j), v in pt.items():
        out[norm(n, perm[i], perm[j])] = v
    return out
