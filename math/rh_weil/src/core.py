"""Exact/algebraic core for Atlas RH/Weil notebook integration.

No RH claim is made by this module. It contains finite-block identities only.
Stdlib-only by design so exact regression tests can run everywhere.
"""
from __future__ import annotations
from math import comb, sqrt
from fractions import Fraction


def overlap_c(i: int, j: int, a, L):
    """C_ij(a;L) for numeric/symbolic-like scalar types supporting arithmetic."""
    d = L - a
    total = 0
    for r in range(j + 1):
        total += comb(j, r) * (a ** (j-r)) * (d ** (i+r+1)) / (i+r+1)
    return total


def kernel_00(a, L):
    return 2*(L-a)


def kernel_01(a, L):
    return L*(L-a)


def kernel_11(a, L):
    d=L-a
    return d*d*(2*L+a)/3


def kernel_q1q1(a, L):
    d=L-a
    return d*(L*L-2*L*a-2*a*a)/6


def kernel_0b(a, L):
    d=L-a
    return d*d*(L+2*a)/3


def kernel_bb(a, L):
    d=L-a
    return d**3*(L*L+3*L*a+a*a)/15


def kernel_bubble_det(a, L):
    d=L-a
    return d**4*(L*L-2*L*a-14*a*a)/45


def scalar_curvature(L):
    from math import exp
    r=exp(L)
    return 2*(r**3-r-1)/(sqrt(r)*(r*r-1))


def q1_sign_threshold():
    return (sqrt(3)-1)/2


def bubble_det_threshold():
    return (sqrt(15)-1)/14


def degree2_raw_det(g00, odd_pivot, even_det, L):
    return even_det + L*L*g00*odd_pivot


def degree2_full_det(odd_pivot, even_det):
    return odd_pivot*even_det
