"""B15-POS2 — PERMITTED / REJECTED certification at Hankel order t = 2 (prereg-004).

Registered system (even truncated Hausdorff problem, degree 4, on [0, 1]):
    H = [mu_{i+j}]_{i,j=0..2}              >= 0   (Hankel wall)
    B = [mu_{i+j+1} - mu_{i+j+2}]_{i,j=0,1} >= 0  (gap wall; localizing weight x(1-x))
Warrant: Krein-Nudel'man Thm III.2.3 as stated in Curto-Fialkow 1991, Remark 4.4.

Primal (SCHUR_PIVOT_EXACT): ``psd_certificate`` on H and on B; PERMITTED iff both
are PSD/PD-certified; witness = both pivot chains.

Dual (DUAL_EXCLUSION_FUNCTIONAL rev. 2): on REJECTED, a rational functional y
over (mu_0..mu_4) with y.mu < 0 and the Gram-form nonnegativity witness
    sum_j y_j x^j = z2^T Q0 z2 + x(1-x) z1^T Q1 z1,   z2 = (1, x, x^2), z1 = (1, x),
checked by coefficient matching through the frozen ``verify_solution`` and with
Q0, Q1 >= 0 checked by ``psd_certificate``. The functional is pivot-derived: if
the first negative Schur pivot of M (M = H or B) sits at index k, the rational
vector v with v_k = 1, v_j = 0 (j > k) and M[:k,:k] v[:k] = -M[:k,k] satisfies
v^T M v = that pivot < 0, and Q = v v^T is the rank-1 Gram witness.
No Farkas companion at t = 2 (the Hankel-wall boundary measure has irrational
atoms, so the equality system is not over Q) — registered in prereg-004 T2.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Dict, List, Optional

from b1_moment_solver.exact import psd_certificate
from pir.symbolic.linear import solve, verify_solution

from ..exact import fmt

Vec = List[Fraction]
Mat = List[List[Fraction]]

OP = "DUAL_EXCLUSION_FUNCTIONAL"
OP_REV = 2


def moments(mu1, mu2, mu3, mu4) -> Vec:
    return [Fraction(1), Fraction(mu1), Fraction(mu2), Fraction(mu3), Fraction(mu4)]


def hankel3(mu: Vec) -> Mat:
    return [[mu[i + j] for j in range(3)] for i in range(3)]


def localizing2(mu: Vec) -> Mat:
    return [[mu[i + j + 1] - mu[i + j + 2] for j in range(2)] for i in range(2)]


def _ok(status: str) -> bool:
    return status in ("PSD_CERTIFIED", "PD_CERTIFIED")


# --------------------------------------------------------------------------- #
# primal                                                                       #
# --------------------------------------------------------------------------- #
def primal(mu: Vec) -> Dict:
    chains = {}
    ok = True
    for name, M in (("hankel", hankel3(mu)), ("gap", localizing2(mu))):
        st, piv, rank = psd_certificate(M)
        chains[name] = {"status": st, "pivots": [fmt(p) for p in piv], "rank": rank}
        ok = ok and _ok(st)
    return {"feasible": ok, "chains": chains, "mu": [fmt(m) for m in mu]}


# --------------------------------------------------------------------------- #
# dual exclusion functional, rev. 2 (Gram witness)                             #
# --------------------------------------------------------------------------- #
def negative_direction(M: Mat) -> Optional[Vec]:
    """Rational v with v^T M v = first negative Schur pivot of M, or None if M is
    PSD-certified. Requires the leading k x k block to be PD (true whenever the
    first failing pivot is negative after k positive pivots)."""
    st, piv, _ = psd_certificate(M)
    if _ok(st):
        return None
    k = len(piv) - 1
    if piv[k] >= 0:
        raise ValueError("indefiniteness via zero pivot is outside the registered design")
    n = len(M)
    v = [Fraction(0)] * n
    v[k] = Fraction(1)
    if k > 0:
        res = solve([row[:k] for row in M[:k]], [-M[i][k] for i in range(k)])
        assert res["status"] == "UNIQUE", res
        v[:k] = res["solution"]
    return v


def _quad(v: Vec, M: Mat) -> Fraction:
    return sum(v[i] * M[i][j] * v[j] for i in range(len(v)) for j in range(len(v)))


def _outer(v: Vec) -> Mat:
    return [[a * b for b in v] for a in v]


def gram_system(Q0: Mat, Q1: Mat):
    """Linear system A q = y expressing the coefficients of
    z2^T Q0 z2 + x(1-x) z1^T Q1 z1 (degree 4) in the unknowns
    q = (Q0 upper triangle, Q1 upper triangle); returns (A, q)."""
    cols, q = [], []
    for i in range(3):
        for j in range(i, 3):
            c = [Fraction(0)] * 5
            c[i + j] += 1 if i == j else 2
            cols.append(c)
            q.append(Q0[i][j])
    for i in range(2):
        for j in range(i, 2):
            m = 1 if i == j else 2
            c = [Fraction(0)] * 5
            c[i + j + 1] += m        # x * x^{i+j}
            c[i + j + 2] -= m        # -x^2 * x^{i+j}
            cols.append(c)
            q.append(Q1[i][j])
    A = [[cols[k][r] for k in range(len(cols))] for r in range(5)]
    return A, q


def verify_gram_witness(y: Vec, Q0: Mat, Q1: Mat) -> Dict:
    A, q = gram_system(Q0, Q1)
    identity = verify_solution(A, y, q)
    s0, p0, _ = psd_certificate(Q0)
    s1, p1, _ = psd_certificate(Q1)
    return {"identity_verified": identity,
            "Q0_status": s0, "Q0_pivots": [fmt(p) for p in p0],
            "Q1_status": s1, "Q1_pivots": [fmt(p) for p in p1],
            "certified": identity and _ok(s0) and _ok(s1)}


def dual_functional(mu: Vec) -> Optional[Dict]:
    zero3 = [[Fraction(0)] * 3 for _ in range(3)]
    zero2 = [[Fraction(0)] * 2 for _ in range(2)]
    v = negative_direction(hankel3(mu))
    if v is not None:
        Q0, Q1, wall = _outer(v), zero2, "hankel (3x3; arXiv:2012.15849 eqs. 7.19/7.46)"
        pivot = _quad(v, hankel3(mu))
    else:
        w = negative_direction(localizing2(mu))
        if w is None:
            return None
        Q0, Q1, wall = zero3, _outer(w), "gap (2x2 x(1-x)-localizing; eqs. 7.45/7.46)"
        pivot = _quad(w, localizing2(mu))
    A, q = gram_system(Q0, Q1)
    y = [sum(A[r][k] * q[k] for k in range(len(q))) for r in range(5)]
    value = sum(a * b for a, b in zip(y, mu))
    assert value == pivot < 0, (value, pivot)
    chk = verify_gram_witness(y, Q0, Q1)
    return {"op": OP, "op_rev": OP_REV, "y": [fmt(t) for t in y], "evaluation": fmt(value),
            "wall": wall, "direction": [fmt(t) for t in (v if v is not None else w)],
            "gram": {"Q0": [[fmt(t) for t in r] for r in Q0], "Q1": [[fmt(t) for t in r] for r in Q1]},
            **chk}


# --------------------------------------------------------------------------- #
# verdict + bracket                                                            #
# --------------------------------------------------------------------------- #
def certify_point(mu1, mu2, mu3, mu4) -> Dict:
    mu = moments(mu1, mu2, mu3, mu4)
    pr = primal(mu)
    out = {"point": {"mu1": fmt(mu[1]), "mu2": fmt(mu[2]), "mu3": fmt(mu[3]), "mu4": fmt(mu[4])},
           "primal": pr, "witness": {"pivots": {k: c["pivots"] for k, c in pr["chains"].items()}}}
    if pr["feasible"]:
        out["verdict"] = "PERMITTED"
        out["impossibility_certificate"] = None
    else:
        out["verdict"] = "REJECTED"
        out["impossibility_certificate"] = {"dual_functional": dual_functional(mu),
                                            "farkas_atom": None}
    return out


def bisect_mu4(slice3, inner: Fraction, outer: Fraction, width: Fraction) -> Dict:
    """Exact bisection on mu_4 at fixed (mu1, mu2, mu3) between a feasible
    ``inner`` and an infeasible ``outer`` value."""
    feas = lambda m4: primal(moments(*slice3, m4))["feasible"]
    assert feas(inner) and not feas(outer)
    steps = 0
    while abs(outer - inner) > width:
        mid = (inner + outer) / 2
        if feas(mid):
            inner = mid
        else:
            outer = mid
        steps += 1
    return {"slice": [fmt(x) for x in slice3],
            "certified_inner_interval": [fmt(min(inner, outer)), fmt(max(inner, outer))],
            "inner_feasible": fmt(inner), "outer_infeasible": fmt(outer),
            "width": fmt(abs(outer - inner)), "steps": steps}


# --------------------------------------------------------------------------- #
# independent routes                                                           #
# --------------------------------------------------------------------------- #
def gauss_legendre2_mu4() -> Fraction:
    """mu_4 of the 2-point Gauss-Legendre rule on [0, 1] (weights 1/2; nodes are the
    roots of x^2 - x + 1/6), computed exactly by reducing x^4 modulo that quadratic:
    x^4 = a x + b  =>  sum_i w_i x_i^4 = a (x_1 + x_2)/2 + b = a/2 + b."""
    # coefficients of x^k mod (x^2 - x + 1/6) as (a, b) meaning a x + b
    a, b = Fraction(1), Fraction(0)          # x^1
    for _ in range(3):                        # multiply by x three times -> x^4
        # x (a x + b) = a x^2 + b x = a (x - 1/6) + b x
        a, b = a + b, -a / 6
    return a / 2 + b


def simpson_mu4() -> Fraction:
    nodes = (Fraction(0), Fraction(1, 2), Fraction(1))
    weights = (Fraction(1, 6), Fraction(2, 3), Fraction(1, 6))
    return sum(w * x ** 4 for w, x in zip(weights, nodes))


def _det(M: Mat) -> Fraction:
    if len(M) == 1:
        return M[0][0]
    return sum((-1) ** j * M[0][j] * _det([r[:j] + r[j + 1:] for r in M[1:]]) for j in range(len(M)))


def tower_minors(mu: Vec) -> List[Dict]:
    """Every contiguous Hankel minor defined through mu_4 of the eq. (7.46) tower:
    the sequence, its first differences, its second differences."""
    d1 = [mu[i] - mu[i + 1] for i in range(4)]
    d2 = [d1[i] - d1[i + 1] for i in range(3)]
    out = []
    for name, seq in (("mu", mu), ("delta1", d1), ("delta2", d2)):
        n = len(seq)
        for size in range(1, (n + 1) // 2 + 1):
            for s in range(0, n - 2 * size + 2):
                M = [[seq[s + i + j] for j in range(size)] for i in range(size)]
                out.append({"sequence": name, "start": s, "size": size, "minor": fmt(_det(M))})
    return out
