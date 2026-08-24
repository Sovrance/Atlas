"""Candidate-A direct-Fourier T=84 objects and exact jets (ENG-005 §6/§7/§9).

The T=84 finite Weil matrix is a *different object* from the cutoff-free entries
of §4/§5: its archimedean term is truncated at ``T = 84`` by definition, not
approximated. So the exact real-space transform (which sums the whole half-line)
does not apply here, and the frequency-space route is the definition. That is
fine at ``T = 84`` -- the panel integrator needs only 8 panels and the oscillation
``cos(Lt)`` has ``t <= 84``, so even an ``L``-box of width ``1e-3`` moves the
phase by under 0.1 and interval evaluation stays usable.

Exact support-length jets (§9)
------------------------------
No finite differences anywhere. Differentiating ``H_i(t;L) = int_0^L h_i(x;L)
e^{itx} dx`` in the support length ``L``:

    d_L^n H0  = (it)^{n-1} e^{itL}                       (n >= 1)
    d_L   Hb  = int_0^L x e^{itx} dx
    d_L^2 Hb  = L e^{itL}
    d_L^n Hb  = e^{itL} [ L(it)^{n-2} + (n-2)(it)^{n-3} ] (n >= 3)
    d_L   Hq1 = (L/2) e^{itL} - H0/2
    d_L^2 Hq1 = (L/2)(it) e^{itL}

``Hq1``'s jets follow from ``h_q1 = x - L/2``: the moving endpoint contributes
``h_q1(L;L) e^{itL} = (L/2)e^{itL}`` and the moving integrand ``-1/2``.

Gram-entry jets come from binomial convolution on the *analytic* product
``P_ij(z) = [H_i(-z)H_j(z) + H_i(z)H_j(-z)]/2``:

    d_L^n P_ij = sum_k C(n,k) [ d^k H_i(-z) d^{n-k} H_j(z)
                              + d^k H_i(z)  d^{n-k} H_j(-z) ] / 2

Candidate-B derivative code is not reachable from here: the pole jets come from
``pole.pole_gram_entry_dL`` / ``_d2L`` and the prime jets from the exact kernel
coefficient expansions.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from math import comb
from typing import Any, Dict, List, Optional, Tuple

import archimedean_realspace as AR
import pole
import weil_entries as WE
from rigorous_integration import PANELS_T84, rigorous_panel_integral

T84 = 84.0
BASIS_KEYS = (("G00", ("one", "one")), ("G0b", ("one", "b")),
              ("Gbb", ("b", "b")), ("O1", ("q1", "q1")))

#: The even block alone -- all E2 needs. Skipping the odd pivot cuts a quarter of
#: the integrals off every box in the uniform cover.
EVEN_KEYS = BASIS_KEYS[:3]

#: Default integrator tolerance. Arb otherwise targets full working precision
#: (~1e-60 at 200 bits), which costs ~4x for accuracy nothing here consumes: the
#: quantities being bounded are ~1e-5 and the claims need ~1e-10. At 140 bits with
#: rel_tol 1e-25 the enclosures still come back with radius ~1e-26.
DEFAULT_INTEGRAL_OPTIONS = {"rel_tol": 1e-25}
DEFAULT_PRECISION_BITS = 140


# --------------------------------------------------------------------------- #
# Exact jets of H_i in the support length                                      #
# --------------------------------------------------------------------------- #
def H_jet(name: str, order: int, z, L, acb):
    """``d_L^order H_name(z; L)`` — exact, analytic in ``z``."""
    iz = acb(0, 1) * z
    Lc = acb(L)
    if order == 0:
        return pole.poly_exp_integral([acb(c) for c in pole.basis_coeffs(name, L)], iz, Lc)
    e = (iz * Lc).exp()
    if name == "one":
        # d_L^n H0 = (it)^{n-1} e^{itL}
        return (iz ** (order - 1)) * e
    if name == "q1":
        if order == 1:
            h0 = pole.poly_exp_integral([acb(1)], iz, Lc)
            return Lc / 2 * e - h0 / 2
        # d_L^2 Hq1 = (L/2)(it)e^{itL}; higher orders by differentiating that
        return Lc / 2 * (iz ** (order - 1)) * e + (order - 2) / 2 * (iz ** (order - 2)) * e \
            if order >= 3 else Lc / 2 * iz * e
    if name == "b":
        if order == 1:
            return pole.poly_exp_integral([acb(0), acb(1)], iz, Lc)
        if order == 2:
            return Lc * e
        return e * (Lc * (iz ** (order - 2)) + (order - 2) * (iz ** (order - 3)))
    raise KeyError(f"unknown basis element {name!r}")


def spectral_product_jet(i: str, j: str, order: int, z, L, acb):
    """``d_L^order`` of ``[H_i(-z)H_j(z) + H_i(z)H_j(-z)]/2`` by binomial convolution."""
    total = acb(0)
    for k in range(order + 1):
        c = comb(order, k)
        total += c * (H_jet(i, k, -z, L, acb) * H_jet(j, order - k, z, L, acb)
                      + H_jet(i, k, z, L, acb) * H_jet(j, order - k, -z, L, acb))
    return total / 2


def arch_jet(i: str, j: str, order: int, L, arb, acb, *, T: float = T84,
             panels=None, options=None):
    """``d_L^order Ginf_ij(L; T)`` — the truncated archimedean term's jet."""
    log_pi = arb.pi().log()
    L_a = arb(L) if not hasattr(L, "mid") else L

    def integrand(z, _analytic):
        return WE.h_plus_analytic(z, arb, acb, log_pi) * spectral_product_jet(
            i, j, order, z, L_a, acb)

    total, record = rigorous_panel_integral(
        integrand, T, acb, panels=panels if panels is not None else list(PANELS_T84),
        options=options)
    return total.real / arb.pi(), record


# --------------------------------------------------------------------------- #
# Prime-block jets (exact polynomial)                                          #
# --------------------------------------------------------------------------- #
def kernel_coeffs_d2L_in_u(i: str, j: str, L, acb) -> List[Any]:
    """Coefficients of ``d^2/dL^2 K_ij(u; L)`` as a polynomial in ``u``."""
    Lc = acb(L)
    if (i, j) == ("one", "one"):
        return [acb(0), acb(0)]
    if {i, j} == {"one", "b"}:
        return [2 * Lc, acb(0), acb(0), acb(0)]
    if (i, j) == ("b", "b"):
        return [4 * Lc**3 / 3, acb(0), -2 * Lc, acb(2) / 3, acb(0), acb(0)]
    if (i, j) == ("q1", "q1"):
        # Kq1q1 = L^3/6 - (L^2/2)u + u^3/3, so d^2/dL^2 = [L, -1, 0, 0].
        # This read L/2 first time round -- the derivative of L^2/2 is L, not L/2 --
        # which threw d^2 O1 off by exactly sum_q w_q (L - log q) = 0.7028 at L=1.25.
        # The exact-jet-vs-finite-difference check caught it; nothing else would have.
        return [Lc, acb(-1), acb(0), acb(0)]
    raise KeyError(f"no d2L kernel expansion for {(i, j)!r}")


def prime_jet(i: str, j: str, order: int, L, arb, acb, prime_powers=None):
    """``d_L^order Gp_ij(L)`` — exact; the prime set is fixed on the open cell."""
    if prime_powers is None:
        mid = float(L.mid()) if hasattr(L, "mid") else float(L)
        prime_powers = WE.prime_powers_below(mid)
    if order == 0:
        return WE.prime_entry(i, j, L, arb, prime_powers)
    coeffs = (AR.kernel_coeffs_dL_in_u(i, j, L, acb) if order == 1
              else kernel_coeffs_d2L_in_u(i, j, L, acb))
    total = arb(0)
    for q, p in prime_powers:
        a = arb(q).log()
        val = acb(0)
        power = acb(1)
        for c in coeffs:
            val += c * power
            power = power * acb(a)
        total += arb(p).log() / arb(q).sqrt() * val.real
    return total


def pole_jet(i: str, j: str, order: int, L):
    """``d_L^order G0_ij(L)`` — Candidate A only."""
    if order == 0:
        return pole.pole_gram_entry(i, j, L)
    if order == 1:
        return pole.pole_gram_entry_dL(i, j, L)
    if order == 2:
        return pole.pole_gram_entry_d2L(i, j, L)
    raise ValueError(f"pole jets implemented to order 2, got {order}")


def entry_jet(i: str, j: str, order: int, L, arb, acb, *, T: float = T84,
              prime_powers=None, panels=None, options=None):
    """``d_L^order G_ij(L; T) = d^n G0 - d^n Gp + d^n Ginf``."""
    arch, record = arch_jet(i, j, order, L, arb, acb, T=T, panels=panels,
                            options=options)
    return (pole_jet(i, j, order, L) - prime_jet(i, j, order, L, arb, acb, prime_powers)
            + arch), record


# --------------------------------------------------------------------------- #
# The T=84 block and its topology                                              #
# --------------------------------------------------------------------------- #
def block_t84(L, arb, acb, *, order: int = 0, T: float = T84, prime_powers=None,
              panels=None, options=None, keys=None) -> Dict[str, Any]:
    """The four entries at derivative ``order``, plus ``E2`` and its jets.

    ``E2 = G00 Gbb - G0b^2``, so by the product rule
    ``E2' = G00' Gbb + G00 Gbb' - 2 G0b G0b'`` and
    ``E2'' = G00'' Gbb + 2 G00' Gbb' + G00 Gbb'' - 2(G0b'^2 + G0b G0b'')``.
    """
    L_a = arb(L) if not hasattr(L, "mid") else L
    if prime_powers is None:
        mid = float(L_a.mid()) if hasattr(L_a, "mid") else float(L_a)
        prime_powers = WE.prime_powers_below(mid)

    if options is None:
        options = DEFAULT_INTEGRAL_OPTIONS
    jets: Dict[str, List[Any]] = {}
    records: Dict[str, Any] = {}
    for key, (i, j) in (keys if keys is not None else BASIS_KEYS):
        vals = []
        for n in range(order + 1):
            v, rec = entry_jet(i, j, n, L_a, arb, acb, T=T,
                               prime_powers=prime_powers, panels=panels,
                               options=options)
            vals.append(v)
            records[f"{key}_d{n}"] = rec
        jets[key] = vals

    out: Dict[str, Any] = {f"{k}": v[0] for k, v in jets.items()}
    for k, v in jets.items():
        for n in range(1, order + 1):
            out[f"{k}_d{n}"] = v[n]

    if not all(k in jets for k in ("G00", "G0b", "Gbb")):
        out["_quadrature"] = records
        return out
    g00, g0b, gbb = jets["G00"], jets["G0b"], jets["Gbb"]
    out["E2"] = g00[0] * gbb[0] - g0b[0] ** 2
    if order >= 1:
        out["E2_d1"] = g00[1] * gbb[0] + g00[0] * gbb[1] - 2 * g0b[0] * g0b[1]
    if order >= 2:
        out["E2_d2"] = (g00[2] * gbb[0] + 2 * g00[1] * gbb[1] + g00[0] * gbb[2]
                        - 2 * (g0b[1] ** 2 + g0b[0] * g0b[2]))
    if "O1" in out:
        out["det_deg2"] = out["O1"] * out["E2"]
    out["_quadrature"] = records
    return out


def topology_scan(arb, acb, *, n_points: int = 25, T: float = T84,
                  cell: Optional[Tuple[float, float]] = None) -> Dict[str, Any]:
    """Fresh Candidate-A scan of the T=84 objects and their first two jets (§6).

    Emphatically **not** a reuse of the Candidate-B monotonicity topology: every
    number here comes from the adopted pole. Reports sign changes of ``E2'`` (the
    stationary points) and of ``E2''`` (curvature changes) as *apparent* features
    located on a grid -- the uniform certificate is what proves anything.
    """
    a, b = cell if cell else (math.log(3.0), math.log(4.0))
    primes = WE.prime_powers_below((a + b) / 2)
    rows: List[Dict[str, Any]] = []
    for k in range(n_points):
        L = a + (b - a) * k / (n_points - 1)
        blk = block_t84(arb(repr(L)), arb, acb, order=2, T=T, prime_powers=primes)
        rows.append({
            "L": repr(L),
            "G00": repr(float(blk["G00"])), "G0b": repr(float(blk["G0b"])),
            "Gbb": repr(float(blk["Gbb"])), "O1": repr(float(blk["O1"])),
            "E2": repr(float(blk["E2"])),
            "E2_d1": repr(float(blk["E2_d1"])), "E2_d2": repr(float(blk["E2_d2"])),
            "G00_d1": repr(float(blk["G00_d1"])), "G00_d2": repr(float(blk["G00_d2"])),
        })

    def sign_changes(key: str):
        out = []
        for p, q in zip(rows, rows[1:]):
            if float(p[key]) * float(q[key]) < 0:
                out.append({"between": [p["L"], q["L"]],
                            "values": [p[key], q[key]]})
        return out

    e2 = [float(r["E2"]) for r in rows]
    argmin = min(range(len(rows)), key=lambda k: e2[k])
    return {
        "T": T,
        "pole_candidate": "A",
        "n_points": n_points,
        "rows": rows,
        "E2_stationary_points": sign_changes("E2_d1"),
        "E2_curvature_changes": sign_changes("E2_d2"),
        "E2_min_on_grid": {"L": rows[argmin]["L"], "E2": rows[argmin]["E2"],
                           "at_endpoint": argmin in (0, len(rows) - 1)},
        "E2_positive_on_grid": all(v > 0 for v in e2),
        "note": ("Apparent features located on a grid; E3 topology evidence only. "
                 "The Candidate-B monotonicity topology is NOT reused."),
    }
