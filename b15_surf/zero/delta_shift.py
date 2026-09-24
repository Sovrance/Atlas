"""The delta-shift Tr(phi^3) -> NLSM (transcribed).

Source: arXiv:2401.05483 §I (R37), eqs. (1)-(2):
    X~_{i,j} = X_{i,j} + delta_{i,j},   delta_{e,e} = -delta_{o,o} = delta,  delta_{o,e} = 0,
i.e.  X_{e,e} -> X_{e,e} + delta,   X_{o,o} -> X_{o,o} - delta,   X_{o,e} unchanged
("e"/"o" = even/odd index; the shift preserves every c_{i,j}, eq. (1)), and
    lim_{delta -> inf}  delta^{2n-2}  A^{delta}_{2n}(X~)  =  A^{NLSM}_{2n}(X)        (2)
"in units with f_pi = 1"; the letter states it is non-trivial that the
first non-vanishing order is O(delta^0) (all O(delta^{-m}) contributions cancel).
Check (transcribed, 4 points):  delta^2 (X13+X24)/((X13-delta)(X24+delta)) -> -(X13+X24) = A^NLSM_4.
Cross-check (arXiv:2312.16282 eq. 7.20, 6 points):
    A^{Tr(phi^3)}_6(X -> X +- delta) -> -(1/delta^4) [ -(X13+X24)(X15+X46)/X14 + (cyclic, i->i+2)
                                                  + X13+X35+X15+X24+X46+X26 ] + O(1/delta^5).

Implementation is EXACT: with eps = 1/delta, each shifted propagator is
1/(X + s*delta) = eps / (s + X*eps) for s = +-1 (or 1/X for s = 0); the amplitude
becomes a truncated power series in eps over Fraction, and the limit (2) is the
coefficient of eps^{2n-2}. The vanishing of every lower coefficient is asserted
(that is the letter's "non-trivial cancellation", verified rather than assumed).
"""

from __future__ import annotations

from fractions import Fraction
from typing import Dict, List, Tuple

from .kinematics import Chord, Point, X, chords
from .triangulations import triangulations

Series = List[Fraction]     # coefficients of eps^0 .. eps^K


def parity_shift(i: int, j: int) -> int:
    """delta_{i,j}/delta: +1 (even,even), -1 (odd,odd), 0 (mixed)."""
    if i % 2 == 0 and j % 2 == 0:
        return 1
    if i % 2 == 1 and j % 2 == 1:
        return -1
    return 0


def shifted_point(pt: Point, n: int, delta: Fraction) -> Point:
    """Finite-delta shifted kinematics X~ (used for the registered shifted dataset)."""
    return {ch: pt[ch] + parity_shift(*ch) * delta for ch in chords(n)}


def _ser_mul(a: Series, b: Series, K: int) -> Series:
    out = [Fraction(0)] * (K + 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j in range(0, K + 1 - i):
            out[i + j] += x * b[j]
    return out


def _inv_linear(s: int, x: Fraction, K: int) -> Series:
    """Series of eps/(s + x*eps) (s = +-1) truncated at eps^K."""
    # 1/(s + x eps) = (1/s) sum_k (-x eps / s)^k
    out = [Fraction(0)] * (K + 1)
    coef = Fraction(1, s)
    for k in range(0, K):
        out[k + 1] = coef          # times eps
        coef *= Fraction(-x, s)
    return out


def shifted_amplitude_series(pt: Point, n: int) -> Series:
    """A^{delta}_n as a truncated series in eps = 1/delta, up to eps^{n-2}."""
    K = n - 2
    total = [Fraction(0)] * (K + 1)
    for T in triangulations(n):
        term = [Fraction(0)] * (K + 1)
        term[0] = Fraction(1)
        for (i, j) in T:
            s = parity_shift(i, j)
            if s == 0:
                term = [t / pt[(i, j)] for t in term]
            else:
                term = _ser_mul(term, _inv_linear(s, pt[(i, j)], K), K)
        total = [a + b for a, b in zip(total, term)]
    return total


def nlsm_from_shift(pt: Point, n: int) -> Tuple[Fraction, List[Fraction]]:
    """lim delta^{n-2} A^delta_n  (eq. 2 with 2n -> n). Returns (limit, lower
    coefficients that must all vanish)."""
    if n % 2:
        raise ValueError("NLSM amplitudes exist for even n only")
    ser = shifted_amplitude_series(pt, n)
    return ser[n - 2], ser[: n - 2]


def nlsm6_eq_7_20(pt: Point, n: int = 6) -> Fraction:
    """Closed form of A^NLSM_6 transcribed from arXiv:2312.16282 eq. (7.20)
    (combined with eq. (2) of 2401.05483: A^NLSM = lim delta^4 A^delta):
       -[ -(X13+X24)(X15+X46)/X14 + (cyclic i->i+2) + X13+X35+X15+X24+X46+X26 ]."""
    def Xv(i, j):
        return X(pt, n, i, j)
    pole = Fraction(0)
    for s in (0, 2, 4):
        a, b, c_, d, e = [((v - 1 + s) % 6) + 1 for v in (1, 2, 3, 4, 5)]
        f = ((6 - 1 + s) % 6) + 1
        # -(X_{1,3}+X_{2,4})(X_{1,5}+X_{4,6})/X_{1,4} shifted by s
        pole += -(Xv(a, c_) + Xv(b, d)) * (Xv(a, e) + Xv(d, f)) / Xv(a, d)
    contact = (Xv(1, 3) + Xv(3, 5) + Xv(1, 5) + Xv(2, 4) + Xv(4, 6) + Xv(2, 6))
    return -(pole + contact)
