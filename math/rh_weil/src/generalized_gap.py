"""Generalized spectral gap of the pencil ``(G, M)`` (ENG-009 §WO-RH-58).

The object measured here is

    lambda_min(G, M)(L) = min { v^T G(L) v / v^T M(L) v : v != 0 },

the smallest generalized eigenvalue of ``G v = lambda M v`` with ``M`` the
exact L^2 reference metric of :mod:`reference_metric`. Unlike raw eigenvalues
or raw determinants, this number does not move under a change of basis applied
to both forms: ``det(S^T G S - lam S^T M S) = det(S)^2 det(G - lam M)``, so the
pencil's roots are invariant (that identity is checked in exact arithmetic in
the tests, and the Rayleigh-quotient form of the same fact is proved in Lean).

Everything rigorous reduces to *shifted positivity*:

* **Lower bounds.** If ``G - lam M`` is positive semidefinite then
  ``v^T G v >= lam * v^T M v`` for every ``v``, so ``lambda_min >= lam``. The
  runtime certifies shifted positivity the same way ENG-008 certified ``G``
  itself -- Sylvester leading minors of the (exactly preconditioned) shifted
  block under adaptive interval covers. No eigensolver anywhere.

* **Upper bounds.** Any single vector ``v`` gives
  ``lambda_min <= v^T G v / v^T M v``; a rational ``v`` evaluated on an interval
  carrier makes that a certified bound. The scouting phase proposes the vector,
  the certificate never trusts how it was found.

The float scouting phase locates the candidate ``lam`` by bisection on the
*sign pattern of leading minors* -- the same Sylvester logic, run in floating
point. It is E3, feeds only the choice of question, and is recorded as such.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import reference_metric as RM

Entries = Dict[Tuple[str, str], Any]


def entry(entries: Entries, i: str, j: str) -> Any:
    return entries[(i, j)] if (i, j) in entries else entries[(j, i)]


def shifted_matrix(basis: Sequence[str], entries: Entries, lam: Any,
                   L: Any) -> List[List[Any]]:
    """``G - lam * M`` on the caller's carrier."""
    out: List[List[Any]] = []
    for i in basis:
        row = []
        for j in basis:
            row.append(entry(entries, i, j) - lam * RM.metric_value(i, j, L))
        out.append(row)
    return out


def leading_minors(matrix: Sequence[Sequence[Any]]) -> List[Any]:
    """Leading principal minors by cofactor expansion, any dimension.

    Division-free, so interval carriers do not widen through pivoting. The
    blocks here are at most 3x3; the general form exists so the next block does
    not need a new function.
    """
    out = []
    for k in range(1, len(matrix) + 1):
        out.append(_det([row[:k] for row in matrix[:k]]))
    return out


def _det(m: Sequence[Sequence[Any]]) -> Any:
    n = len(m)
    if n == 1:
        return m[0][0]
    total = None
    for col in range(n):
        minor = [row[:col] + row[col + 1:] for row in m[1:]]
        term = m[0][col] * _det(minor)
        if col % 2:
            term = -term
        total = term if total is None else total + term
    return total


def precondition(matrix: Sequence[Sequence[Any]],
                 exponents: Sequence[int]) -> List[List[Any]]:
    """``D A D`` with ``D = diag(2^e)`` -- an exact congruence on any carrier."""
    return [[matrix[a][b] * (2.0 ** (exponents[a] + exponents[b]))
             for b in range(len(matrix))] for a in range(len(matrix))]


# --------------------------------------------------------------------------- #
# E3 scouting: locate the gap by float Sylvester bisection                     #
# --------------------------------------------------------------------------- #
def _floats(basis: Sequence[str], entries: Entries) -> Entries:
    return {k: float(v) for k, v in entries.items()}


def _sylvester_pd(basis: Sequence[str], entries: Entries, lam: float,
                  L: float) -> bool:
    return all(m > 0 for m in leading_minors(shifted_matrix(basis, entries, lam, L)))


def scout_gap_at(basis: Sequence[str], entries: Entries, L: float,
                 *, tol: float = 1e-12) -> float:
    """Float bisection for the largest ``lam`` with ``G - lam M`` PD at ``L``.

    Sign checks on leading minors only -- no eigensolver even in the scout, so
    the rigorous path and the scout disagree about arithmetic, never about
    method.
    """
    entries = _floats(basis, entries)
    lo = 0.0
    if not _sylvester_pd(basis, entries, lo, L):
        return 0.0
    hi = 1.0
    while _sylvester_pd(basis, entries, hi, L):
        hi *= 2.0
        if hi > 1e6:  # pragma: no cover - the pencil is bounded in practice
            raise AssertionError("scout runaway: G - lam M stayed PD past 1e6")
    while hi - lo > tol * max(1.0, hi):
        mid = 0.5 * (lo + hi)
        if _sylvester_pd(basis, entries, mid, L):
            lo = mid
        else:
            hi = mid
    return lo


def scout_min_eigvec(basis: Sequence[str], entries: Entries, lam: float,
                     L: float) -> List[Fraction]:
    """A rational near-kernel vector of ``G - lam M`` at the scouted gap.

    At ``lam`` just past the crossing the shifted matrix is nearly singular;
    the adjugate's largest column is a numerically fine kernel direction. The
    vector is *proposed* here and *verified* by the certified Rayleigh quotient
    -- a bad proposal weakens the upper bound, it cannot make one wrong.
    """
    entries = _floats(basis, entries)
    m = shifted_matrix(basis, entries, lam, L)
    n = len(basis)
    cols = []
    for b in range(n):
        col = []
        for a in range(n):
            minor = [[m[r][c] for c in range(n) if c != b]
                     for r in range(n) if r != a]
            cof = _det(minor) if minor else 1.0
            col.append(cof * (-1.0) ** (a + b))
        cols.append(col)
    best = max(cols, key=lambda c: sum(x * x for x in c))
    norm = max(abs(x) for x in best) or 1.0
    return [Fraction(x / norm).limit_denominator(10 ** 6) for x in best]


# --------------------------------------------------------------------------- #
# E1: certified bounds                                                        #
# --------------------------------------------------------------------------- #
def shifted_minors_over(basis: Sequence[str],
                        assemble: Callable[[Any], Entries],
                        lam: Fraction, box: Any, L_carrier: Any,
                        exponents: Optional[Sequence[int]] = None
                        ) -> List[Any]:
    """Enclosures of the leading minors of ``D (G - lam M) D`` over a box.

    ``lam`` is an exact dyadic rational, applied as ``* p / q`` so the shift
    itself adds no rounding; ``exponents`` is the frozen preconditioner (powers
    of two, exact congruence, minor signs unchanged).
    """
    entries = assemble(box)
    lamc = (L_carrier * 0 + 1) * lam.numerator / lam.denominator
    m = shifted_matrix(basis, entries, lamc, L_carrier)
    if exponents is not None:
        m = precondition(m, exponents)
    return leading_minors(m)


def rayleigh_upper(basis: Sequence[str], entries: Entries,
                   v: Sequence[Fraction], L: Any) -> Any:
    """Certified enclosure of ``v^T G v / v^T M v`` -- an upper bound carrier.

    Its upper endpoint bounds ``lambda_min(G, M)`` from above at that ``L``.
    """
    num = None
    den = None
    for a, i in enumerate(basis):
        for b, j in enumerate(basis):
            c = Fraction(v[a]) * Fraction(v[b])
            if not c:
                continue
            gterm = entry(entries, i, j) * c.numerator / c.denominator
            mterm = RM.metric_value(i, j, L) * c.numerator / c.denominator
            num = gterm if num is None else num + gterm
            den = mterm if den is None else den + mterm
    return num / den
