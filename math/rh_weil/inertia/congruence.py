"""Exact congruence algebra and independent inertia oracles (§4, WO-RH-29).

Two jobs.

**Congruence.** ``A -> S* A S`` for invertible ``S`` preserves inertia
(Sylvester). This module builds such transforms exactly over the rationals so
the regression layer can hit :mod:`inertia.ldl` with matrices whose signature is
known by construction but whose entries look nothing like the original.

**An independent oracle.** :func:`charpoly_inertia` computes the inertia from
the characteristic polynomial instead of by elimination, so it shares no code
path with the LDL engine. For a real symmetric matrix every eigenvalue is real,
and for a real-rooted polynomial Descartes' rule of signs is exact rather than
an upper bound: the number of sign changes in the coefficient sequence *equals*
the number of positive roots. That turns a signature into a finite exact count
with no root finding, no floating point, and no shared assumptions with the
elimination it is checking.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import Any, List, Sequence, Tuple

Mat = List[List[Fraction]]


def to_fraction_matrix(A: Sequence[Sequence[Any]]) -> Mat:
    return [[Fraction(x) for x in row] for row in A]


def matmul(A: Mat, B: Mat) -> Mat:
    n, m, p = len(A), len(B), len(B[0])
    return [[sum((A[i][k] * B[k][j] for k in range(m)), Fraction(0))
             for j in range(p)] for i in range(n)]


def transpose(A: Mat) -> Mat:
    return [list(col) for col in zip(*A)]


def congruence(A: Mat, S: Mat) -> Mat:
    """``S^T A S`` — exact, and inertia-preserving when ``S`` is invertible."""
    return matmul(matmul(transpose(S), A), S)


def determinant(A: Mat) -> Fraction:
    """Exact determinant by fraction-free-ish Gaussian elimination."""
    n = len(A)
    M = [list(row) for row in A]
    det = Fraction(1)
    for k in range(n):
        piv = next((i for i in range(k, n) if M[i][k] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != k:
            M[k], M[piv] = M[piv], M[k]
            det = -det
        det *= M[k][k]
        for i in range(k + 1, n):
            if M[i][k] == 0:
                continue
            f = M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
    return det


def is_invertible(A: Mat) -> bool:
    return determinant(A) != 0


def principal_minor_sums(A: Mat) -> List[Fraction]:
    """``[e_1, ..., e_n]`` where ``e_k`` sums the ``k x k`` principal minors.

    These are the elementary symmetric functions of the eigenvalues, so they are
    exactly the characteristic-polynomial coefficients up to sign.
    """
    n = len(A)
    out: List[Fraction] = []
    for k in range(1, n + 1):
        total = Fraction(0)
        for idx in combinations(range(n), k):
            sub = [[A[i][j] for j in idx] for i in idx]
            total += determinant(sub)
        out.append(total)
    return out


def charpoly_coeffs(A: Mat) -> List[Fraction]:
    """Coefficients of ``det(xI - A)``, constant term first, monic last."""
    n = len(A)
    e = principal_minor_sums(A)
    # det(xI - A) = x^n - e1 x^{n-1} + e2 x^{n-2} - ... + (-1)^n e_n
    coeffs = [Fraction(0)] * (n + 1)
    coeffs[n] = Fraction(1)
    for k in range(1, n + 1):
        coeffs[n - k] = (-1) ** k * e[k - 1]
    return coeffs


def _sign_changes(seq: Sequence[Fraction]) -> int:
    changes = 0
    last = 0
    for v in seq:
        s = (v > 0) - (v < 0)
        if s == 0:
            continue
        if last != 0 and s != last:
            changes += 1
        last = s
    return changes


def charpoly_inertia(A: Sequence[Sequence[Any]]) -> Tuple[int, int, int]:
    """``(n_positive, n_negative, n_zero)`` of an exact symmetric matrix.

    Valid because the matrix is symmetric, hence real-rooted, which is exactly
    the condition under which Descartes' rule of signs is an equality. The zero
    eigenvalue's multiplicity is the number of trailing zero coefficients.
    """
    M = to_fraction_matrix(A)
    n = len(M)
    for i in range(n):
        for j in range(n):
            if M[i][j] != M[j][i]:
                raise ValueError("charpoly_inertia requires a symmetric matrix")
    c = charpoly_coeffs(M)          # index i is the coefficient of x^i
    n_zero = 0
    while n_zero < len(c) and c[n_zero] == 0:
        n_zero += 1
    reduced = c[n_zero:]            # divide out x^n_zero; remaining roots nonzero
    ascending = list(reduced)
    n_pos = _sign_changes(list(reversed(ascending)))
    # Negative roots of p(x) are positive roots of p(-x): flip odd coefficients.
    flipped = [v if i % 2 == 0 else -v for i, v in enumerate(ascending)]
    n_neg = _sign_changes(list(reversed(flipped)))
    return n_pos, n_neg, n_zero


def inertia_2x2(trace: Fraction, det: Fraction) -> Tuple[int, int, int]:
    """Closed-form inertia of a symmetric 2x2 from its trace and determinant.

    The odd degree-3 Weil block is 2x2, so this is the whole story there and it
    needs no elimination at all: the eigenvalue product is ``det`` and their sum
    is ``trace``.
    """
    t, d = Fraction(trace), Fraction(det)
    if d < 0:
        return (1, 1, 0)
    if d > 0:
        return (2, 0, 0) if t > 0 else (0, 2, 0)
    if t > 0:
        return (1, 0, 1)
    if t < 0:
        return (0, 1, 1)
    return (0, 0, 2)


def diagonal(signs: Sequence[int]) -> Mat:
    return [[Fraction(signs[i]) if i == j else Fraction(0)
             for j in range(len(signs))] for i in range(len(signs))]


def unimodular(n: int, seed: int, steps: int = 6) -> Mat:
    """A deterministic invertible rational ``S`` built from shear operations.

    Shears have determinant 1, so the product is invertible by construction --
    no determinant check has to succeed for the test to be meaningful, and the
    inertia of ``S^T A S`` is forced to equal that of ``A``.
    """
    S = [[Fraction(1) if i == j else Fraction(0) for j in range(n)] for i in range(n)]
    state = seed
    for _ in range(steps):
        state = (state * 1103515245 + 12345) % (1 << 31)
        i = state % n
        state = (state * 1103515245 + 12345) % (1 << 31)
        j = state % n
        if i == j:
            j = (j + 1) % n
        state = (state * 1103515245 + 12345) % (1 << 31)
        num = (state % 9) - 4 or 1
        state = (state * 1103515245 + 12345) % (1 << 31)
        den = (state % 4) + 1
        f = Fraction(num, den)
        for c in range(n):
            S[i][c] += f * S[j][c]
    return S
