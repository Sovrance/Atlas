"""The odd degree-3 Weil block (ATLAS-RH-ENG-006 §7, WO-RH-32).

Midpoint parity splits the basis. About ``x = L/2``::

    one(x) = 1                   even
    b(x)   = x(L-x)              even
    q1(x)  = x - L/2             odd
    b3(x)  = x(L-x)(x-L/2)       odd

so the Gram matrix is block diagonal and the odd degree-3 block is the 2x2

    [ G[q1,q1]  G[q1,b3] ]
    [ G[q1,b3]  G[b3,b3] ]

assembled the same way as everything else in the program,
``G = G0 - Gp + Ginf``: the Candidate-A pole through :mod:`pole`, the prime
block from the exact kernels, and the archimedean term from the ENG-005
real-space form, which is what makes an interval ``L`` usable at all.

Being 2x2 is a genuine simplification, not just a small case. The inertia of a
symmetric 2x2 is fixed by its trace and determinant with no elimination
(:func:`inertia.congruence.inertia_2x2`), its eigenvalues have a closed form, and
its first two spectral moments invert to the spectrum exactly. So every channel
ENG-006 builds -- inertia, moments, rank-trace -- has a sharp answer here rather
than a bound.

What is deliberately *not* attempted
------------------------------------
§7 records that the individual active prime-shift block is indefinite on this
cell, and forbids termwise PSD domination. That is preserved as a regression
test (:func:`prime_shift_block`) rather than worked around: the assembled entry
can be positive while a single shift's contribution is not, and any argument
that bounds the assembly by dominating it term by term is unsound here. The
indefiniteness is a fact about the decomposition, not an obstacle to route
around.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import archimedean_realspace as AR
import pole
import weil_entries as WE

CELL = (math.log(3.0), math.log(4.0))
CELL_LABEL = ("log(3)", "log(4)")
CLAIM_SCOPE = "finite_dimensional_weil_compression"

#: The odd block, in the order its entries appear in the 2x2.
ODD_KEYS: Tuple[Tuple[str, Tuple[str, str]], ...] = (
    ("Oqq", ("q1", "q1")),
    ("Oqb", ("q1", "b3")),
    ("Obb", ("b3", "b3")),
)

#: Parity of each basis element about ``x = L/2`` (§7).
PARITY = {"one": "even", "b": "even", "q1": "odd", "b3": "odd"}

DEFAULT_PRECISION_BITS = 160
DEFAULT_INTEGRAL_OPTIONS = {"rel_tol": 1e-20}


def basis_value(name: str, x, L):
    """``h_name(x; L)`` from the pole primitive's own coefficients.

    Deliberately reuses :func:`pole.basis_coeffs` rather than restating the
    polynomials: a second copy of a basis definition is exactly the kind of
    duplication that lets two parts of the program drift apart.
    """
    coeffs = pole.basis_coeffs(name, L)
    total = 0 * x
    power = 0 * x + 1
    for c in coeffs:
        total = total + c * power
        power = power * x
    return total


def odd_block(L, arb, acb, *, prime_powers=None, options=None) -> Dict[str, Any]:
    """The 2x2 odd block plus its trace and determinant.

    ``L`` may be an Arb ball; every entry is then an enclosure valid for every
    ``L`` in it, so the result feeds the inertia engine directly.
    """
    L_a = arb(L) if not hasattr(L, "mid") else L
    if prime_powers is None:
        mid = float(L_a.mid()) if hasattr(L_a, "mid") else float(L_a)
        prime_powers = WE.prime_powers_below(mid)
    if options is None:
        options = DEFAULT_INTEGRAL_OPTIONS

    out: Dict[str, Any] = {}
    records: Dict[str, Any] = {}
    for key, (i, j) in ODD_KEYS:
        val, rec = AR.gram_entry_realspace(i, j, L_a, arb, acb,
                                           prime_powers=prime_powers,
                                           options=options)
        out[key] = val
        records[key] = rec

    out["trace"] = out["Oqq"] + out["Obb"]
    out["det"] = out["Oqq"] * out["Obb"] - out["Oqb"] * out["Oqb"]
    out["_quadrature"] = records
    return out


def odd_matrix(block: Dict[str, Any]) -> List[List[Any]]:
    """The block as a 2x2 matrix, ready for the inertia and moment engines."""
    return [[block["Oqq"], block["Oqb"]], [block["Oqb"], block["Obb"]]]


def odd_matrix_at(L, arb, acb, *, prime_powers=None, options=None) -> List[List[Any]]:
    return odd_matrix(odd_block(L, arb, acb, prime_powers=prime_powers,
                                options=options))


# --------------------------------------------------------------------------- #
# §7 regression: the prime-shift block is indefinite, and stays that way       #
# --------------------------------------------------------------------------- #
def prime_shift_block(q: int, L, arb) -> List[List[Any]]:
    """The single-shift contribution ``K_ij(log q; L)`` as a 2x2 symmetric block.

    This is one term of ``Gp``, not an entry of ``G``. §7 records that it is
    indefinite on this cell and forbids termwise PSD domination; the test suite
    pins that so nobody re-derives a bound by dominating the assembly term by
    term.
    """
    a = arb(q).log()
    kqq = WE.kernel("q1", "q1", a, L)
    kqb = WE.kernel("q1", "b3", a, L)
    kbb = WE.kernel("b3", "b3", a, L)
    return [[kqq, kqb], [kqb, kbb]]


def prime_shift_determinant(q: int, L, arb):
    M = prime_shift_block(q, L, arb)
    return M[0][0] * M[1][1] - M[0][1] * M[1][0]


def active_shifts(L_value: float) -> List[Tuple[int, int]]:
    """The prime powers active on the open cell at ``L_value``."""
    return WE.prime_powers_below(L_value)


# --------------------------------------------------------------------------- #
# Exact identities (E0)                                                        #
# --------------------------------------------------------------------------- #
def parity_identities() -> Dict[str, Any]:
    """The block-diagonal statement the odd block relies on (§7)."""
    return {
        "basis_parity": dict(PARITY),
        "statement": ("h(L - x) = +h(x) for the even elements and -h(x) for the odd "
                      "ones, so an even/odd Gram entry is an integral of an odd "
                      "function about L/2 and vanishes"),
        "odd_block": [k for k, _ in ODD_KEYS],
        "even_block": ["G00", "G0b", "Gbb"],
        "consequence": ("the degree-3 Gram is block diagonal; the odd block is 2x2 "
                        "and can be treated on its own"),
    }


# --------------------------------------------------------------------------- #
# §8 E3 scan                                                                   #
# --------------------------------------------------------------------------- #
def _f(x) -> float:
    return float(x.mid()) if hasattr(x, "mid") else float(x)


def eigenvalues_2x2(trace, det, arb):
    """``(lambda_min, lambda_max)`` in closed form -- no eigenvalue solver.

    ``lambda = (T +- sqrt(T^2 - 4D)) / 2``. For a real symmetric 2x2 the
    discriminant is non-negative; a computed enclosure can still dip below zero,
    and the lower end is re-seated at 0 outward-safely rather than producing a
    complex result from a real-rooted problem.
    """
    disc = trace * trace - 4 * det
    if disc.lower() < 0:
        from interval_backend import interval_box

        disc = interval_box(0.0, max(0.0, float(disc.upper())))
    root = disc.sqrt()
    return (trace - root) / 2, (trace + root) / 2


def scan_row(Lv: float, arb, acb, *, prime_powers=None, options=None) -> Dict[str, Any]:
    """One E3 scan row: entries, spectrum, inertia diagnostic, moments (§8)."""
    from inertia.congruence import inertia_2x2

    L = arb(repr(Lv))
    blk = odd_block(L, arb, acb, prime_powers=prime_powers, options=options)
    tr, det = blk["trace"], blk["det"]
    lo, hi = eigenvalues_2x2(tr, det, arb)
    lo_f, hi_f = _f(lo), _f(hi)
    # m1..m4 of a 2x2 from its eigenvalues -- exact in terms of trace and det,
    # so no matrix powers are needed here.
    m1 = lo_f + hi_f
    m2 = lo_f ** 2 + hi_f ** 2
    m3 = lo_f ** 3 + hi_f ** 3
    m4 = lo_f ** 4 + hi_f ** 4
    cond = abs(hi_f / lo_f) if lo_f != 0 else float("inf")
    return {
        "L": repr(Lv),
        "Oqq": repr(_f(blk["Oqq"])),
        "Oqb": repr(_f(blk["Oqb"])),
        "Obb": repr(_f(blk["Obb"])),
        "trace": repr(_f(tr)),
        "det": repr(_f(det)),
        "lambda_min": repr(lo_f),
        "lambda_max": repr(hi_f),
        "floating_inertia_diagnostic": list(
            inertia_2x2(round(_f(tr), 15), round(_f(det), 15))),
        "condition_number": repr(cond),
        "m1": repr(m1), "m2": repr(m2), "m3": repr(m3), "m4": repr(m4),
    }


def topology_scan(arb, acb, *, n_points: int = 41, cell=None,
                  options=None) -> Dict[str, Any]:
    """Fresh Candidate-A scan of the odd degree-3 block (§8). E3 only."""
    a, b = cell if cell else CELL
    prime_powers = WE.prime_powers_below((a + b) / 2)
    rows = [scan_row(a + (b - a) * k / (n_points - 1), arb, acb,
                     prime_powers=prime_powers, options=options)
            for k in range(n_points)]

    def crossings(key: str):
        out = []
        for r0, r1 in zip(rows, rows[1:]):
            v0, v1 = float(r0[key]), float(r1[key])
            if v0 == 0.0 or v1 == 0.0 or (v0 < 0) != (v1 < 0):
                out.append({"between": [r0["L"], r1["L"]],
                            "values": [r0[key], r1[key]]})
        return out

    dets = [float(r["det"]) for r in rows]
    lams = [float(r["lambda_min"]) for r in rows]
    i_min = min(range(len(rows)), key=lambda k: dets[k])
    sigs = {tuple(r["floating_inertia_diagnostic"]) for r in rows}
    return {
        "cell": [repr(a), repr(b)],
        "pole_candidate": "A",
        "n_points": n_points,
        "prime_shifts_active": [list(t) for t in prime_powers],
        "rows": rows,
        "det_zero_crossings": crossings("det"),
        "lambda_min_zero_crossings": crossings("lambda_min"),
        "Oqb_sign_changes": crossings("Oqb"),
        "det_min_on_grid": {"L": rows[i_min]["L"], "det": rows[i_min]["det"],
                            "at_endpoint": i_min in (0, len(rows) - 1)},
        "lambda_min_min_on_grid": repr(min(lams)),
        "distinct_floating_inertias": [list(s) for s in sorted(sigs)],
        "apparent_constant_inertia": len(sigs) == 1,
    }
