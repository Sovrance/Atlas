"""Forward-limit moment / Hankel construction for the registered EFT-hedron bound.

Source: Arkani-Hamed, Huang & Huang, "The EFT-Hedron", arXiv:2012.15849 (R26).
  eq. (2.16)/(3.1)  M_IR(s,t) = {massless poles} + sum_{k,q} a_{k,q} s^{k-q} t^q
                    (k = total Mandelstam degree, q = degree in t; s = (p1+p2)^2,
                    t = (p2+p3)^2, u = -s-t; a_{k,q} <-> dimension 2k+4 operators).
  eq. (7.17)        forward limit (q = 0):  a_{k,0} = sum_a p'_a (x_a)^k,   p'_a > 0,
                    x_a = 1/m_a^2, valid for k >= 2.
  eq. (7.18)-(7.19) the vector (1, a_{3,0}/a_{2,0}, a_{4,0}/a_{2,0}, ...) lies in the
                    convex hull of the half moment curve (1, x, x^2, ...), x > 0, so
                    its Hankel matrix K[a~_0] is totally positive.
  eq. (7.34)-(7.35) with a mass gap M_Gap (x_a <= 1 in units M_Gap = 1):
                    a_{2,0} >= M_Gap^2 a_{3,0} >= ... >= M_Gap^{2(k-2)} a_{k,0} >= 0.

Registered target (prereg-002 §POS, awaiting §9-OI-1 confirmation) — the
lowest-order TWO-SIDED forward-limit statement, in units a_{2,0} = 1, M_Gap = 1:
    mu_0 = 1,  mu_1 = a_{3,0}/a_{2,0},  mu_2 = a_{4,0}/a_{2,0}
    lower wall (Hankel 2x2 minor, eq. 7.19):     mu_0 mu_2 - mu_1^2 >= 0   i.e.  mu_2 >= mu_1^2
    upper wall (gap chain, eq. 7.35):            mu_1 - mu_2 >= 0          i.e.  mu_2 <= mu_1
so at the registered slice mu_1 = 1/2:   1/4 <= mu_2 <= 1/2 .

Truncation: we keep moments through a_{4,0} (Hankel order t = 1, the first
order at which a two-sided bound exists). The truncation order is a HEURISTIC
choice and is carried as a located warning; every arithmetic step is exact.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Dict, List, Tuple

from b1_moment_solver.exact import hankel

Mat = List[List[Fraction]]

# Registered slice and points (prereg-002 §POS).
MU1_SLICE = Fraction(1, 2)
PUBLISHED_LOWER_WALL = MU1_SLICE ** 2           # 1/4  (eq. 7.19 at mu_1 = 1/2)
PUBLISHED_UPPER_WALL = MU1_SLICE                # 1/2  (eq. 7.35 at mu_1 = 1/2)
INTERIOR_POINT = (MU1_SLICE, Fraction(3, 8))
EXTERIOR_POINT_HANKEL = (MU1_SLICE, Fraction(1, 5))     # below the Hankel wall
EXTERIOR_POINT_GAP = (MU1_SLICE, Fraction(3, 5))        # above the gap wall
PERTURBATION = Fraction(-1, 4)                          # T4 negative control on mu_2
TRUNCATION_ORDER = 1


def moments(mu1: Fraction, mu2: Fraction) -> List[Fraction]:
    """(mu_0, mu_1, mu_2) = (1, a30/a20, a40/a20)."""
    return [Fraction(1), Fraction(mu1), Fraction(mu2)]


def hankel_matrix(mu: List[Fraction]) -> Mat:
    """K[a~_0] truncated at order t = 1: [[mu0, mu1],[mu1, mu2]]  (eq. 7.18)."""
    return hankel(mu, TRUNCATION_ORDER)


def gap_matrices(mu: List[Fraction]) -> List[Mat]:
    """The gap chain (7.35) as 1x1 'shifted Hankel' blocks: [mu_j - mu_{j+1}] >= 0."""
    return [[[mu[j] - mu[j + 1]]] for j in range(len(mu) - 1)]


def warnings_for_truncation() -> List[Dict[str, str]]:
    return [{"location": "b15_surf/pos/dispersive.py::TRUNCATION_ORDER",
             "text": "moment sequence truncated at a_{4,0} (Hankel order t=1); "
                     "higher-order EFT-hedron walls are not tested (HEURISTIC choice)"}]
