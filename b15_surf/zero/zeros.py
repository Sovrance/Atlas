"""Registered hidden-zero loci of Tr(phi^3) tree amplitudes (transcribed).

Source: arXiv:2312.16282 (R37).
  §3.2 (p.16) "Zeros": pick a planar variable X_B and the causal diamond anchored
  on it; setting all c_{i,j} inside the diamond to zero makes A_n vanish.
  eq. (6.14) — explicit index ranges for the diamond anchored at X_{1,i} in the
  ray-like triangulation (1,3),...,(1,n-1) (field-theory case n_{a,b} = 0):
        c_{a,b} = 0   for  1 <= a <= i-2 ,  i <= b <= n-1 .
  §6.1 (p.26): "the total number of zeros for the Tr(phi^3) field theory
  amplitude is n(n-3)/2, equal to the number of X_{i,j}'s".
  §4.1: the i = 3 diamond is the "skinny rectangle" c_{1,j} = 0 for all j not
  adjacent to 1 (the zero that implies the NLSM Adler zero).

The remaining loci are the cyclic images (i -> i+k) of (6.14); the mesh's
Mobius symmetry c_{i,j} = c_{j,i} makes some images coincide, and the
de-duplicated count must equal n(n-3)/2 (asserted by the tests).

Points on a locus are produced EXACTLY: the constraints are linear in the X's
(eq. 2.4), so the locus is the null space of the constraint matrix, obtained
from the frozen ``pir.symbolic.linear.solve`` (status UNDERDETERMINED) and
sampled with seeded rational coefficients.
"""

from __future__ import annotations

import random
from fractions import Fraction
from typing import Dict, List, Tuple

from pir.symbolic.linear import solve, verify_solution

from .kinematics import (Chord, Point, c_rows, chords, norm, point_from_vector,
                         random_rational)


def diamond_6_14(n: int, i: int) -> List[Chord]:
    """c-pairs of the causal diamond anchored at X_{1,i}: eq. (6.14), n_{a,b}=0."""
    if not (3 <= i <= n - 1):
        raise ValueError("i must satisfy 3 <= i <= n-1")
    return [(a, b) for a in range(1, i - 1) for b in range(i, n)]


def _canon_pairs(n: int, pairs: List[Chord]) -> Tuple[Chord, ...]:
    return tuple(sorted({norm(n, a, b) for (a, b) in pairs}))


def registered_loci(n: int) -> Dict[str, List[Chord]]:
    """All n(n-3)/2 zero loci as {locus_id: [(a,b) with c_{a,b}=0]}.
    locus_id = "Z{n}-i{i}-k{k}" names the (6.14) anchor i and the cyclic shift k
    of the FIRST representative found for that locus."""
    seen: Dict[Tuple[Chord, ...], str] = {}
    out: Dict[str, List[Chord]] = {}
    for i in range(3, n):
        base = diamond_6_14(n, i)
        for k in range(n):
            shifted = [(a + k, b + k) for (a, b) in base]
            key = _canon_pairs(n, shifted)
            if key in seen:
                continue
            lid = f"Z{n}-i{i}-k{k}"
            seen[key] = lid
            out[lid] = list(key)
    return out


def locus_matrix(n: int, pairs: List[Chord]) -> List[List[Fraction]]:
    return c_rows(n, pairs)


def locus_null_space(n: int, pairs: List[Chord]) -> List[List[Fraction]]:
    A = locus_matrix(n, pairs)
    res = solve(A, [Fraction(0)] * len(A))
    if res["status"] == "UNIQUE":
        return []
    if res["status"] == "INCONSISTENT":      # cannot happen for b = 0
        raise RuntimeError("homogeneous system reported inconsistent")
    return res["null_space"]


def sample_on_locus(n: int, pairs: List[Chord], rng: random.Random,
                    max_tries: int = 200) -> Point:
    """Exact point on the locus with every X nonzero (seeded)."""
    null = locus_null_space(n, pairs)
    A = locus_matrix(n, pairs)
    for _ in range(max_tries):
        coeffs = [random_rational(rng) for _ in null]
        v = [Fraction(0)] * len(chords(n))
        for cf, basis in zip(coeffs, null):
            for k in range(len(v)):
                v[k] += cf * basis[k]
        if all(x != 0 for x in v):
            assert verify_solution(A, [Fraction(0)] * len(A), v)
            return point_from_vector(n, v)
    raise RuntimeError("could not sample a locus point with all X nonzero")


def random_constraints(n: int, m: int, rng: random.Random) -> List[List[Fraction]]:
    """Negative control: m random rational linear constraints on the X's (dense,
    so no single chord is forced to zero)."""
    k = len(chords(n))
    return [[random_rational(rng) for _ in range(k)] for _ in range(m)]


def sample_on_constraints(n: int, A: List[List[Fraction]], rng: random.Random,
                          max_tries: int = 200) -> Point:
    res = solve(A, [Fraction(0)] * len(A))
    null = res.get("null_space", [])
    for _ in range(max_tries):
        v = [Fraction(0)] * len(chords(n))
        for basis in null:
            cf = random_rational(rng)
            for k in range(len(v)):
                v[k] += cf * basis[k]
        if all(x != 0 for x in v):
            return point_from_vector(n, v)
    raise RuntimeError("could not sample a control point with all X nonzero")
