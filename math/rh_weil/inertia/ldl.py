"""Verified Hermitian LDL* / congruence reduction (ATLAS-RH-ENG-006 §3).

The engine answers one question about a finite Hermitian matrix ``A``::

    Inertia(A) = (n_positive, n_negative, n_zero)

by symmetric (congruence) elimination. Sylvester's law of inertia is what makes
this rigorous: for any invertible ``S``, ``S* A S`` has the same inertia as
``A``. Symmetric row/column elimination and symmetric permutation are both
congruences, so the signature of the pivots *is* the signature of ``A``.

Why this is sound on intervals
------------------------------
Run the elimination on a matrix of balls and record the pivot order. For any
point matrix inside the box, the *same* pivot order produces point pivots lying
inside the corresponding pivot balls -- interval arithmetic encloses the point
computation step by step. So if every pivot ball has a determined sign, every
point matrix in the box has that same signature, and the box has one constant
inertia. If any pivot ball straddles zero the point signs are not determined by
these data and the answer is ``INCONCLUSIVE`` -- never a guess.

Two pivot kinds are needed, not one
-----------------------------------
1x1 diagonal pivots alone are not enough: ``[[0, 1], [1, 0]]`` has inertia
``(1, 1, 0)`` and no usable diagonal entry at all. So when no diagonal pivot has
a determined sign, the engine looks for a symmetric 2x2 block whose determinant
is *definitely negative*. Such a block has inertia exactly ``(1, 1)`` whatever
its diagonal does -- its two eigenvalues have product < 0 -- and it can be
eliminated as a unit via the Schur complement. Without this, ordinary indefinite
matrices would report INCONCLUSIVE for no better reason than pivot order.

Exact zero is not numerical zero (§14.3)
----------------------------------------
The interval path never reports a nonzero ``n_zero``. A ball containing zero is
not a proof that the underlying entry *is* zero, and a zero-radius ball produced
by arithmetic is not one either. Exact zero multiplicity is reported only on the
exact rational path, where a vanishing pivot with a vanishing row is a structural
fact rather than a measurement.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Callable, List, Optional, Sequence, Tuple

#: Returned instead of an inertia when the arithmetic cannot determine a sign.
INCONCLUSIVE = "INCONCLUSIVE"


class SignOracle:
    """How to read the sign of a scalar, and whether exact zero is knowable."""

    name = "abstract"
    #: True when a proven-zero scalar is meaningful (exact arithmetic only).
    exact_zero_available = False

    def sign(self, x) -> Optional[int]:
        """``+1`` / ``-1`` / ``0``, or ``None`` when undetermined."""
        raise NotImplementedError

    def is_zero(self, x) -> bool:
        raise NotImplementedError

    def describe(self, x) -> Any:
        raise NotImplementedError


class ExactSignOracle(SignOracle):
    """Exact rational arithmetic: every sign is decidable, zero included."""

    name = "exact_rational"
    exact_zero_available = True

    def sign(self, x) -> Optional[int]:
        return (x > 0) - (x < 0)

    def is_zero(self, x) -> bool:
        return x == 0

    def describe(self, x) -> Any:
        return str(x)


class BallSignOracle(SignOracle):
    """Arb balls: a sign is readable only when the ball excludes zero."""

    name = "arb_interval"
    exact_zero_available = False

    def sign(self, x) -> Optional[int]:
        if x.lower() > 0:
            return 1
        if x.upper() < 0:
            return -1
        return None

    def is_zero(self, x) -> bool:
        # Deliberately never true: see the module docstring. An interval cannot
        # witness exact zero, so claiming it would manufacture evidence.
        return False

    def describe(self, x) -> Any:
        return [repr(float(x.lower())), repr(float(x.upper()))]


@dataclass
class PivotRecord:
    step: int
    kind: str                 # "1x1" or "2x2"
    indices: Tuple[int, ...]
    sign: Tuple[int, ...]     # signature contributed by this pivot
    value: Any                # pivot (or 2x2 determinant) as reported


@dataclass
class InertiaResult:
    status: str                       # "PASS" or "INCONCLUSIVE"
    n_positive: Optional[int] = None
    n_negative: Optional[int] = None
    n_zero: Optional[int] = None
    method: str = ""
    pivots: List[PivotRecord] = field(default_factory=list)
    blocker: Optional[str] = None
    oracle: str = ""

    @property
    def signature(self) -> Optional[Tuple[int, int, int]]:
        if self.status != "PASS":
            return None
        return (self.n_positive, self.n_negative, self.n_zero)

    def to_dict(self) -> dict:
        return {
            "status": self.status,
            "n_positive": self.n_positive,
            "n_negative": self.n_negative,
            "n_zero": self.n_zero,
            "method": self.method,
            "sign_oracle": self.oracle,
            "pivot_intervals": [
                {"step": p.step, "kind": p.kind, "indices": list(p.indices),
                 "signature": list(p.sign), "value": p.value}
                for p in self.pivots
            ],
            **({"blocker": self.blocker} if self.blocker else {}),
        }


def _copy(A: Sequence[Sequence[Any]]) -> List[List[Any]]:
    return [list(row) for row in A]


def _mirror(A: List[List[Any]], lo: int, hi: int) -> None:
    """Copy the upper triangle of the trailing block onto the lower.

    The elimination only writes ``c >= i``, so the lower triangle is stale
    afterwards. Mathematically the two halves agree; mirroring keeps that true
    of the stored balls as well, and makes the upper triangle authoritative
    rather than leaving the result dependent on update order.
    """
    for i in range(lo, hi):
        for c in range(i + 1, hi):
            A[c][i] = A[i][c]


def _swap_symmetric(A: List[List[Any]], i: int, j: int) -> None:
    """Symmetric permutation — a congruence, so it cannot change the inertia."""
    if i == j:
        return
    A[i], A[j] = A[j], A[i]
    for row in A:
        row[i], row[j] = row[j], row[i]


def ldl_inertia(A: Sequence[Sequence[Any]], oracle: SignOracle) -> InertiaResult:
    """Inertia of Hermitian ``A`` by congruence reduction under ``oracle``."""
    n = len(A)
    for row in A:
        if len(row) != n:
            raise ValueError("matrix must be square")
    W = _copy(A)
    _mirror(W, 0, n)  # upper triangle is authoritative
    pos = neg = zer = 0
    pivots: List[PivotRecord] = []
    k = 0
    step = 0

    while k < n:
        step += 1
        # --- 1x1 pivot: any remaining diagonal entry with a determined sign.
        chosen = None
        for j in range(k, n):
            s = oracle.sign(W[j][j])
            if s in (1, -1):
                chosen = (j, s)
                break

        if chosen is not None:
            j, s = chosen
            _swap_symmetric(W, k, j)
            d = W[k][k]
            pivots.append(PivotRecord(step, "1x1", (k,), (s,), oracle.describe(d)))
            pos += 1 if s > 0 else 0
            neg += 1 if s < 0 else 0
            # Multipliers come first, from the *pre-update* column. Computing
            # them inside the update loop would let an already-eliminated row
            # feed back into a later one and subtract the pivot row twice.
            mult = [W[i][k] / d for i in range(k + 1, n)]
            for i in range(k + 1, n):
                f = mult[i - k - 1]
                if oracle.is_zero(f):
                    continue
                for c in range(i, n):
                    W[i][c] = W[i][c] - f * W[k][c]
            _mirror(W, k + 1, n)
            k += 1
            continue

        # --- exact path only: a structurally zero row is a null direction.
        if oracle.exact_zero_available:
            zj = None
            for j in range(k, n):
                if oracle.is_zero(W[j][j]) and all(
                    oracle.is_zero(W[j][c]) for c in range(k, n)
                ):
                    zj = j
                    break
            if zj is not None:
                _swap_symmetric(W, k, zj)
                pivots.append(PivotRecord(step, "1x1", (k,), (0,),
                                          oracle.describe(W[k][k])))
                zer += 1
                k += 1
                continue

        # --- 2x2 pivot: a symmetric block with definitely negative determinant
        #     contributes exactly (1, 1) whatever its diagonal looks like.
        found = None
        for a in range(k, n):
            for b in range(a + 1, n):
                det = W[a][a] * W[b][b] - W[a][b] * W[b][a]
                if oracle.sign(det) == -1:
                    found = (a, b, det)
                    break
            if found:
                break

        if found is None:
            return InertiaResult(
                status=INCONCLUSIVE,
                method="hermitian_ldl_congruence",
                pivots=pivots,
                oracle=oracle.name,
                blocker=(
                    f"no 1x1 pivot has a determined sign and no 2x2 sub-block has a "
                    f"definitely negative determinant, on the trailing {n - k}x{n - k} "
                    f"block at step {step}"
                ),
            )

        a, b, _ = found
        # b > a >= k, so the first swap never moves b: it exchanges k with a,
        # and b equals neither.
        _swap_symmetric(W, k, a)
        _swap_symmetric(W, k + 1, b)
        p11, p12, p22 = W[k][k], W[k][k + 1], W[k + 1][k + 1]
        det = p11 * p22 - p12 * W[k + 1][k]
        pivots.append(PivotRecord(step, "2x2", (k, k + 1), (1, -1),
                                  oracle.describe(det)))
        pos += 1
        neg += 1
        # Schur complement against the 2x2 block: S = A22 - A21 P^-1 A12,
        # with P^-1 = [[p22, -p12], [-p12, p11]] / det.
        # Same rule as the 1x1 case: all multipliers before any update.
        gs = []
        for i in range(k + 2, n):
            u, v = W[i][k], W[i][k + 1]
            gs.append(((u * p22 - v * p12) / det, (v * p11 - u * p12) / det))
        for i in range(k + 2, n):
            gu, gv = gs[i - k - 2]
            for c in range(i, n):
                W[i][c] = W[i][c] - (gu * W[k][c] + gv * W[k + 1][c])
        _mirror(W, k + 2, n)
        k += 2

    return InertiaResult(
        status="PASS",
        n_positive=pos,
        n_negative=neg,
        n_zero=zer,
        method="hermitian_ldl_congruence",
        pivots=pivots,
        oracle=oracle.name,
    )


def exact_inertia(A: Sequence[Sequence[Fraction]]) -> InertiaResult:
    """Inertia of an exact rational symmetric matrix. Zeros are real zeros."""
    M = [[Fraction(x) for x in row] for row in A]
    return ldl_inertia(M, ExactSignOracle())


def interval_inertia(A: Sequence[Sequence[Any]]) -> InertiaResult:
    """Inertia of a matrix of Arb balls, valid for every point in the box."""
    return ldl_inertia(A, BallSignOracle())
