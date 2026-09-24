"""PERMITTED / REJECTED certification of a forward-limit point with primal AND
dual certificates, plus an exact bisection bracket along a registered ray.

Primal (SCHUR_PIVOT_EXACT, frozen op): ``b1_moment_solver.exact.psd_certificate``
on the Hankel block and on each gap block; PERMITTED iff every chain is
PSD/PD-certified; ``witness = {pivots}``.

Dual (exclusion functional): on REJECTED we return an exact rational linear
functional  ell(mu) = y . mu  which is >= 0 on the WHOLE truncated moment cone
and < 0 at the point — the bootstrap "dual certificate" (R26: He–Kruczenski
dual convex problem; EFT-hedron walls). Nonnegativity of ell on the cone is
itself certified by an exact polynomial identity (Markov–Lukács form on [0,1]):
    sum_j y_j x^j  =  (v0 + v1 x)^2 * s0  +  x (1 - x) * s1 ,   s0, s1 >= 0,
so that ell(mu) = integral of a nonnegative polynomial against the (positive)
spectral measure. The identity is checked by coefficient matching through the
frozen ``pir.symbolic.linear`` bridge (``verify_solution``); the evaluation
y . mu_ext < 0 is exact.

NOTE (reconciliation, recorded not repaired): the work order asked for
``verify_farkas`` on this dual. Farkas' lemma in ``pir/symbolic/linear.py``
certifies inconsistency of a LINEAR-EQUALITY system; the moment-cone exclusion
is a conic (SDP-type) dual, so the honest exact primitive is the
polynomial-identity check above. ``verify_farkas`` IS exercised, as a genuine
certificate, on the equality system "mu_ext is the moment vector of the single
atom forced by the flat (rank-1) boundary" (see ``farkas_atom_certificate``):
for a Hankel-wall exterior point that system is inconsistent and the Farkas
vector is an exact witness. Both are stored in ``impossibility_certificate``.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Dict, List, Optional, Tuple

from b1_moment_solver.exact import psd_certificate
from pir.symbolic.linear import solve, verify_farkas, verify_solution

from ..exact import fmt
from .dispersive import gap_matrices, hankel_matrix, moments

Vec = List[Fraction]


# --------------------------------------------------------------------------- #
# primal                                                                       #
# --------------------------------------------------------------------------- #
def primal(mu1: Fraction, mu2: Fraction) -> Dict:
    mu = moments(mu1, mu2)
    H = hankel_matrix(mu)
    st, piv, rank = psd_certificate(H)
    chains = {"hankel": {"status": st, "pivots": [fmt(p) for p in piv], "rank": rank}}
    ok = st in ("PSD_CERTIFIED", "PD_CERTIFIED")
    for j, G in enumerate(gap_matrices(mu)):
        gst, gpiv, grank = psd_certificate(G)
        chains[f"gap_{j}"] = {"status": gst, "pivots": [fmt(p) for p in gpiv], "rank": grank}
        ok = ok and gst in ("PSD_CERTIFIED", "PD_CERTIFIED")
    return {"feasible": ok, "chains": chains, "mu": [fmt(m) for m in mu]}


# --------------------------------------------------------------------------- #
# dual exclusion functional                                                    #
# --------------------------------------------------------------------------- #
def _poly_mul(a: Vec, b: Vec) -> Vec:
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def dual_functional(mu1: Fraction, mu2: Fraction) -> Optional[Dict]:
    """Exclusion functional y (over mu_0, mu_1, mu_2) with an exact
    nonnegativity witness, or None if the point is feasible."""
    mu = moments(mu1, mu2)
    H = hankel_matrix(mu)
    st, piv, _ = psd_certificate(H)
    if st == "NOT_PSD_CERTIFIED":
        # negative Schur pivot d = mu2 - mu1^2/mu0 < 0 ; direction v = (-mu1/mu0, 1):
        # v^T H v = d.  ell(mu) = v0^2 mu0 + 2 v0 v1 mu1 + v1^2 mu2 = moment of (v0+v1 x)^2.
        v0, v1 = -mu[1] / mu[0], Fraction(1)
        y = [v0 * v0, 2 * v0 * v1, v1 * v1]
        witness = {"form": "square", "v": [fmt(v0), fmt(v1)], "s0": "1/1", "s1": "0/1"}
        wall = "hankel (eq. 7.19)"
    else:
        # gap wall: ell(mu) = mu_j - mu_{j+1} = moment of x^j (1 - x)
        bad = None
        for j, G in enumerate(gap_matrices(mu)):
            gst, _, _ = psd_certificate(G)
            if gst == "NOT_PSD_CERTIFIED":
                bad = j
                break
        if bad is None:
            return None
        y = [Fraction(0)] * 3
        y[bad] += 1
        y[bad + 1] -= 1
        witness = {"form": "interval", "v": ["0/1", "0/1"], "s0": "0/1", "s1": "1/1",
                   "x_power": bad}
        wall = "gap chain (eq. 7.35)"
    value = sum(yi * mi for yi, mi in zip(y, mu))
    assert value < 0
    return {"y": [fmt(t) for t in y], "evaluation": fmt(value), "wall": wall,
            "nonnegativity_witness": witness,
            "identity_verified": verify_nonnegativity_identity(y, witness)}


def verify_nonnegativity_identity(y: Vec, witness: Dict) -> bool:
    """Check  sum_j y_j x^j == s0 (v0 + v1 x)^2 + s1 x^p (1 - x)  by coefficient
    matching, expressed as the linear system A c = b in c = (s0, s1) and checked
    with the frozen ``verify_solution``."""
    y = [Fraction(t) for t in y]
    v0, v1 = (Fraction(t) for t in witness["v"])
    s0, s1 = Fraction(witness["s0"]), Fraction(witness["s1"])
    p = int(witness.get("x_power", 0))
    sq = _poly_mul([v0, v1], [v0, v1])                       # (v0 + v1 x)^2
    iv = [Fraction(0)] * p + [Fraction(1), Fraction(-1)]     # x^p (1 - x)
    deg = max(len(sq), len(iv), len(y))
    sq += [Fraction(0)] * (deg - len(sq))
    iv += [Fraction(0)] * (deg - len(iv))
    yy = y + [Fraction(0)] * (deg - len(y))
    A = [[sq[k], iv[k]] for k in range(deg)]
    return s0 >= 0 and s1 >= 0 and verify_solution(A, yy, [s0, s1])


def farkas_atom_certificate(mu1: Fraction, mu2: Fraction) -> Optional[Dict]:
    """Genuine Farkas usage: on the Hankel wall the moment problem is FLAT
    (rank 1) and the representing measure is a single atom at x* = mu1/mu0
    with weight mu0. The linear system  [w, w x*, w x*^2] = mu  in the unknown
    weight w is inconsistent for a Hankel-exterior point; ``solve`` returns the
    Farkas vector y with y^T A = 0, y^T b != 0, re-checked by ``verify_farkas``."""
    mu = moments(mu1, mu2)
    xs = mu[1] / mu[0]
    A = [[Fraction(1)], [xs], [xs * xs]]
    res = solve(A, mu)
    if res["status"] != "INCONSISTENT":
        return None
    y = res["farkas"]
    ok = verify_farkas(A, mu, y)
    return {"system": "single-atom moment equations at x*=mu1/mu0",
            "A": [[fmt(a) for a in row] for row in A], "b": [fmt(m) for m in mu],
            "farkas_y": [fmt(t) for t in y], "verify_farkas": ok}


# --------------------------------------------------------------------------- #
# verdict + bracket                                                            #
# --------------------------------------------------------------------------- #
def certify_point(mu1: Fraction, mu2: Fraction) -> Dict:
    pr = primal(mu1, mu2)
    out = {"point": {"mu1": fmt(mu1), "mu2": fmt(mu2)}, "primal": pr}
    if pr["feasible"]:
        out["verdict"] = "PERMITTED"
        out["witness"] = {"pivots": {k: v["pivots"] for k, v in pr["chains"].items()}}
        out["impossibility_certificate"] = None
    else:
        out["verdict"] = "REJECTED"
        out["witness"] = {"pivots": {k: v["pivots"] for k, v in pr["chains"].items()}}
        dual = dual_functional(mu1, mu2)
        # The Farkas (equality-form) companion is meaningful only on the Hankel
        # wall, where the boundary measure is flat (rank 1); on the gap wall the
        # exclusion rests on the dual functional alone.
        farkas = farkas_atom_certificate(mu1, mu2) if dual["wall"].startswith("hankel") else None
        out["impossibility_certificate"] = {"dual_functional": dual, "farkas_atom": farkas}
    return out


def bisect_ray(mu1: Fraction, inner: Fraction, outer: Fraction,
               width: Fraction) -> Dict:
    """Exact bisection on mu_2 between a feasible ``inner`` and an infeasible
    ``outer`` value at fixed mu_1; returns the certified bracket
    [inner feasible, outer infeasible] with |outer - inner| <= width."""
    assert primal(mu1, inner)["feasible"] and not primal(mu1, outer)["feasible"]
    steps = 0
    while abs(outer - inner) > width:
        mid = (inner + outer) / 2
        if primal(mu1, mid)["feasible"]:
            inner = mid
        else:
            outer = mid
        steps += 1
    return {"mu1": fmt(mu1), "certified_inner_interval": [fmt(min(inner, outer)), fmt(max(inner, outer))],
            "inner_feasible": fmt(inner), "outer_infeasible": fmt(outer),
            "width": fmt(abs(outer - inner)), "steps": steps}
