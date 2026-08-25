"""Real-space form of the archimedean Weil term (ATLAS-RH-ENG-005 §4/§5/§8).

The frequency-space definition

    Ginf_ij(L) = (1/pi) int_0^inf h_+(t) Re(conj(H_i) H_j) dt

is an oscillatory integral over a half-line. It is fine at a *point* ``L`` --
that is what ENG-004's canary does -- but hopeless on an ``L``-interval: the
factor ``cos(Lt)`` has argument width ``rad(L)*t``, so a box of any useful width
destroys the enclosure long before ``t`` reaches the cutoff. That is why ENG-004
had to reach for convexity instead of subdividing.

Transforming to real space removes the oscillation entirely. Using
``Re(conj(H_i)H_j) = int_0^L K_ij(u) cos(tu) du`` (the same ``K_ij`` as the prime
block) together with the transform proved in ``curvature_derivation``,

    int_0^inf cos(Lt) a/(a^2+t^2/4) dt = pi e^{-2aL},

and the digamma series ``h_+(t) = -gamma - log pi + sum_n [1/(n+1) -
a_n/(a_n^2+t^2/4)]``, ``a_n = n+1/4``, one gets

    Ginf_ij(L) = (K(0)/2) h_+(0)
               + int_0^L [K(0) - K(u)] w(u) du
               + K(0) * S(L),

    w(u) = e^{-u/2}/(1 - e^{-2u}),   S(L) = sum_{n>=0} e^{-(2n+1/2)L}/(2n+1/2).

Every piece is benign on an ``L``-interval: ``K`` is an exact polynomial, the
integral runs over the *compact* ``[0, L]`` with no oscillation, and ``S``
converges geometrically (ratio ``e^{-2L} < 0.09`` on this cell).

Why the naive rearrangement fails
---------------------------------
Interchanging directly gives ``-int_0^L K(u) w(u) du``, which **diverges**:
``w(u) ~ 1/(2u)`` while ``K(0) != 0``. The constant part of ``h_+`` contributes a
delta at ``u = 0`` that the naive swap drops. Keeping the ``F(0)`` terms is what
produces the convergent ``K(0) - K(u)`` numerator above, and that difference
vanishes linearly at the origin exactly cancelling the ``1/u``.

The same interchange hypothesis recorded in ``curvature_derivation`` applies
here; this module is E0 algebra conditional on it, and the frequency-space route
in ``weil_entries`` remains as an independent numerical cross-check.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial
from typing import Any, Dict, List, Optional, Tuple

import core
import basis_algebra
import weil_entries as WE

#: Bernoulli numbers B_n^+ (B_1 = +1/2), for u/(1-e^{-u}) = sum_n B_n^+ u^n/n!.
_BERNOULLI_PLUS: Tuple[Fraction, ...] = (
    Fraction(1), Fraction(1, 2), Fraction(1, 6), Fraction(0), Fraction(-1, 30),
    Fraction(0), Fraction(1, 42), Fraction(0), Fraction(-1, 30), Fraction(0),
    Fraction(5, 66), Fraction(0), Fraction(-691, 2730), Fraction(0), Fraction(7, 6),
    Fraction(0), Fraction(-3617, 510), Fraction(0), Fraction(43867, 798),
)
_BERNOULLI_RADIUS = 1.0  # |u| below which the series branch is used


def _u_over_one_minus_exp_neg(u, acb):
    """``u / (1 - e^{-u})``, entire at ``u = 0`` where it equals 1.

    Series branch for ``|u| <= 1`` with an explicit remainder. Using
    ``|B_n| <= 4 n!/(2pi)^n`` for ``n >= 2``, the tail past ``N`` terms is at most
    ``4 (|u|/2pi)^{N+1} / (1 - |u|/2pi)``.
    """
    import pole

    mag = pole._mag(u)
    if mag > _BERNOULLI_RADIUS:
        return u / (1 - (-u).exp())
    total = acb(0)
    power = acb(1)
    for n, b in enumerate(_BERNOULLI_PLUS):
        if b != 0:
            total += acb(float(b.numerator)) / acb(float(b.denominator)) * power
        power = power * u / (n + 1)
    n_terms = len(_BERNOULLI_PLUS)
    ratio = mag / 6.283185307179586
    tail = 4.0 * (ratio ** n_terms) / max(1e-300, 1.0 - ratio)
    return total + pole._ball(acb(0), tail)


def kernel_difference_over_u(i: str, j: str, u, L, acb):
    """``[K_ij(0; L) - K_ij(u; L)] / u`` — an exact polynomial, no singularity.

    ``K_ij(0) - K_ij(u)`` vanishes at ``u = 0``, so the quotient is a polynomial
    and needs no series branch. Computed by exact synthetic division on the
    coefficients rather than by dividing two numbers, which would be a ``0/0``
    ball at the origin.
    """
    coeffs = kernel_coeffs_in_u(i, j, L, acb)
    # p(u) = K(u); we need (p(0) - p(u))/u = -(c1 + c2 u + c3 u^2 + ...)
    out = acb(0)
    power = acb(1)
    for c in coeffs[1:]:
        out += c * power
        power = power * u
    return -out


def kernel_coeffs_in_u(i: str, j: str, L, acb) -> List[Any]:
    """Coefficients of ``K_ij(u; L)`` as a polynomial in ``u`` (``L`` fixed).

    ENG-008 §WO-RH-48: derived from the basis coefficients by
    :mod:`basis_algebra`, not tabulated. The hand-written table this replaces
    had one entry per pair and raised ``KeyError`` for anything else, so every
    new basis element meant deriving and pasting six polynomial identities into
    three modules. ``tests/test_kernel_algebra.py`` pins the derived values
    against that table's entries, which ENG-005 and ENG-006 had verified against
    SymPy, so the generalization is checked rather than assumed.
    """
    return basis_algebra.kernel_coeffs_in_a(i, j, acb(L))


def kernel_at_zero(i: str, j: str, L, arb):
    """``K_ij(0; L)``."""
    return WE.kernel(i, j, arb(0) if not hasattr(L, "mid") else 0 * L, L)


def geometric_tail(L, arb, terms: int = 64):
    """``S(L) = sum_{n>=0} e^{-(2n+1/2)L}/(2n+1/2)``, with a rigorous remainder.

    Terms decay by ``e^{-2L} < 0.09`` on this cell, so the truncation ball is
    ``<= e^{-(2N+1/2)L}/((2N+1/2)(1-e^{-2L}))``.
    """
    total = arb(0)
    for n in range(terms):
        c = 2 * arb(n) + arb(1) / 2
        total += (-c * L).exp() / c
    cN = 2 * arb(terms) + arb(1) / 2
    ratio = (-2 * L).exp()
    rem = (-cN * L).exp() / (cN * (1 - ratio))
    return total + arb(0, float(rem.upper()))


def h_plus_at_zero(arb, acb):
    """``h_+(0) = psi(1/4) - log pi``."""
    return acb(arb("0.25")).digamma().real - arb.pi().log()


def arch_entry_realspace(i: str, j: str, L, arb, acb, *, options=None,
                         series_terms: int = 64):
    """``Ginf_ij(L)`` exactly, with no frequency cutoff and no oscillation.

    ``L`` may be an Arb **ball**; the result is a valid enclosure for every ``L``
    in it. Returns ``(value, record)``.
    """
    L_a = arb(L) if not hasattr(L, "mid") else L
    K0 = kernel_at_zero(i, j, L_a, arb)

    # Substitute u = L*s so the limits are the fixed [0, 1] and L appears only
    # inside the integrand. Integrating to sup(L) instead would be wrong on an
    # L-ball: for L' < sup(L) the true integral stops at L', so the enclosure
    # would include mass the real integral never sees.
    L_c = acb(L_a)

    def integrand(s, _analytic):
        u = L_c * s
        pref = kernel_difference_over_u(i, j, u, L_a, acb)
        return pref * (-u / 2).exp() * _u_over_one_minus_exp_neg(2 * u, acb) / 2

    integral = L_c * acb.integral(integrand, 0, 1, **(options or {}))
    if not integral.is_finite():
        raise ValueError(f"real-space archimedean integral did not converge for {(i, j)}")

    S = geometric_tail(L_a, arb, terms=series_terms)
    value = K0 / 2 * h_plus_at_zero(arb, acb) + integral.real + K0 * S
    record = {
        "method": "realspace_kernel_transform",
        "substitution": "u = L*s, s in [0,1] (fixed limits; L-interval safe)",
        "compact_domain": [0.0, 1.0],
        "geometric_terms": series_terms,
        "integral_radius": float(integral.real.rad()),
        "oscillatory": False,
    }
    return value, record


def gram_entry_realspace(i: str, j: str, L, arb, acb, *, prime_powers=None,
                         options=None):
    """``G_ij(L) = G0_ij - Gp_ij + Ginf_ij`` with the exact archimedean term."""
    import pole

    L_a = arb(L) if not hasattr(L, "mid") else L
    arch, record = arch_entry_realspace(i, j, L_a, arb, acb, options=options)
    g0 = pole.pole_gram_entry(i, j, L_a)
    gp = WE.prime_entry(i, j, L_a, arb, prime_powers)
    return g0 - gp + arch, record


def block_realspace(L, arb, acb, *, prime_powers=None, options=None) -> Dict[str, Any]:
    """``G00``, ``G0b``, ``Gbb``, ``O1`` and the derived ``E2`` / ``det``."""
    L_a = arb(L) if not hasattr(L, "mid") else L
    if prime_powers is None:
        mid = float(L_a.mid()) if hasattr(L_a, "mid") else float(L_a)
        prime_powers = WE.prime_powers_below(mid)
    entries: Dict[str, Any] = {}
    records: Dict[str, Any] = {}
    for key, (i, j) in (("G00", ("one", "one")), ("G0b", ("one", "b")),
                        ("Gbb", ("b", "b")), ("O1", ("q1", "q1"))):
        val, rec = gram_entry_realspace(i, j, L_a, arb, acb,
                                        prime_powers=prime_powers, options=options)
        entries[key] = val
        records[key] = rec
    entries["E2"] = entries["G00"] * entries["Gbb"] - entries["G0b"] ** 2
    entries["det_deg2"] = entries["O1"] * entries["E2"]
    entries["_quadrature"] = records
    return entries


# --------------------------------------------------------------------------- #
# Exact L-derivative and the centred enclosure                                 #
# --------------------------------------------------------------------------- #
# Differentiating the real-space form in L:
#
#   Ginf'(L) = (K0'/2) h_+(0)
#            + [K0(L) - K(L;L)] w(L)                      <- boundary term
#            + int_0^L [K0' - d_L K(u;L)] w(u) du
#            + K0' S(L) + K0 S'(L).
#
# Every kernel carries a factor (L-a), so K(L;L) = 0 and the boundary term is
# K0(L) w(L); and S'(L) = -sum_n e^{-c_n L} = -w(L), so K0 S'(L) = -K0 w(L).
# The two cancel **exactly**, leaving the same shape as the value itself with
# K replaced by d_L K:
#
#   Ginf'(L) = (K0'/2) h_+(0) + int_0^L [K0' - d_L K(u;L)] w(u) du + K0' S(L).
#
# The bracket still vanishes at u = 0 (since d_L K(0;L) = K0'), so the 1/u in w
# is cancelled just as before. No finite differences anywhere (ENG-005 §9).
def kernel_coeffs_dL_in_u(i: str, j: str, L, acb) -> List[Any]:
    """Coefficients of ``d/dL K_ij(u; L)`` as a polynomial in ``u``.

    Differentiated exactly, from the same bivariate table the value comes from
    (ENG-008 §WO-RH-48). No finite differences anywhere (ENG-005 §9).
    """
    return basis_algebra.kernel_dL_coeffs_in_a(i, j, acb(L))


def _difference_over_u(coeffs: List[Any], u, acb):
    """``[p(0) - p(u)]/u`` for the polynomial with those coefficients."""
    out = acb(0)
    power = acb(1)
    for c in coeffs[1:]:
        out += c * power
        power = power * u
    return -out


def arch_entry_dL_realspace(i: str, j: str, L, arb, acb, *, options=None,
                            series_terms: int = 64):
    """``d/dL Ginf_ij(L)`` — exact, same shape as the value with ``K -> d_L K``."""
    L_a = arb(L) if not hasattr(L, "mid") else L
    L_c = acb(L_a)
    dcoeffs = kernel_coeffs_dL_in_u(i, j, L_a, acb)
    K0p = dcoeffs[0]

    def integrand(s, _analytic):
        u = L_c * s
        pref = _difference_over_u(dcoeffs, u, acb)
        return pref * (-u / 2).exp() * _u_over_one_minus_exp_neg(2 * u, acb) / 2

    integral = L_c * acb.integral(integrand, 0, 1, **(options or {}))
    if not integral.is_finite():
        raise ValueError(f"real-space derivative integral did not converge for {(i, j)}")
    S = geometric_tail(L_a, arb, terms=series_terms)
    return (K0p / 2 * acb(h_plus_at_zero(arb, acb)) + integral
            + K0p * acb(S)).real


def arch_entry_centred(i: str, j: str, box, arb, acb, *, options=None,
                       series_terms: int = 64):
    """``Ginf_ij`` on an ``L``-**box**, via the mean-value (centred) form.

    A direct interval evaluation is badly lossy here: the three terms
    ``(K0/2)h_+(0)``, the integral and ``K0 S`` are each of order 1 and largely
    cancel (the ``q1`` entry is ~0.04 from terms of ~0.9), so their widths add
    while the true variation does not. Amplification is ~5000x, which no
    practical subdivision recovers.

    The centred form evaluates the value at the midpoint -- a *point*, so no
    dependency loss at all -- and adds ``[-r, r]`` times an enclosure of the
    derivative over the box. The derivative enclosure may itself be loose; it is
    only multiplied by the radius.
    """
    mid = arb(box.mid()) if hasattr(box, "mid") else arb(box)
    rad = float(box.rad()) if hasattr(box, "rad") else 0.0
    centre, record = arch_entry_realspace(i, j, mid, arb, acb, options=options,
                                          series_terms=series_terms)
    if rad == 0.0:
        record["form"] = "point"
        return centre, record
    deriv = arch_entry_dL_realspace(i, j, box, arb, acb, options=options,
                                    series_terms=series_terms)
    slope_mag = max(abs(float(deriv.lower())), abs(float(deriv.upper())))
    record["form"] = "centred_mean_value"
    record["radius"] = rad
    record["derivative_enclosure"] = [repr(float(deriv.lower())), repr(float(deriv.upper()))]
    return centre + arb(0, rad * slope_mag), record


def prime_entry_dL(i: str, j: str, L, arb, acb, prime_powers=None):
    """``d/dL Gp_ij(L) = sum_q (log p/sqrt q) d_L K_ij(log q; L)``.

    Exact: the prime set is fixed on the open cell (its breakpoints are the
    endpoints), so only the kernels depend on ``L``.
    """
    if prime_powers is None:
        mid = float(L.mid()) if hasattr(L, "mid") else float(L)
        prime_powers = WE.prime_powers_below(mid)
    dcoeffs = kernel_coeffs_dL_in_u(i, j, L, acb)
    total = arb(0)
    for q, p in prime_powers:
        a = arb(q).log()
        val = acb(0)
        power = acb(1)
        for c in dcoeffs:
            val += c * power
            power = power * acb(a)
        total += arb(p).log() / arb(q).sqrt() * val.real
    return total


def gram_entry_dL(i: str, j: str, L, arb, acb, *, prime_powers=None, options=None):
    """``d/dL G_ij = d_L G0 - d_L Gp + d_L Ginf`` — every piece exact (§9)."""
    import pole

    return (pole.pole_gram_entry_dL(i, j, L)
            - prime_entry_dL(i, j, L, arb, acb, prime_powers)
            + arch_entry_dL_realspace(i, j, L, arb, acb, options=options))


def gram_entry_centred(i: str, j: str, box, arb, acb, *, prime_powers=None,
                       options=None):
    """``G_ij`` on an ``L``-box via the mean-value form applied to the **whole** entry.

    Centring only the archimedean term is not enough. The three blocks
    ``G0``, ``-Gp`` and ``Ginf`` are individually of order 1 and cancel heavily --
    ``Gbb`` is ~3.6e-5 built from pieces of order 0.1 -- so evaluating them
    separately on a box makes their widths add while the true variation stays
    tiny. ``d/dL Gbb`` is ~1e-3, four orders below the per-block slopes.

    So the centre value is a single *point* evaluation of the assembled entry (no
    dependency loss), and the radius term uses an enclosure of the assembled
    derivative, where the same cancellation is taken analytically rather than by
    interval arithmetic.
    """
    import pole

    mid = arb(box.mid()) if hasattr(box, "mid") else arb(box)
    rad = float(box.rad()) if hasattr(box, "rad") else 0.0
    if prime_powers is None:
        prime_powers = WE.prime_powers_below(float(mid))

    arch_c, record = arch_entry_realspace(i, j, mid, arb, acb, options=options)
    centre = (pole.pole_gram_entry(i, j, mid)
              - WE.prime_entry(i, j, mid, arb, prime_powers) + arch_c)
    if rad == 0.0:
        record["form"] = "point"
        return centre, record

    deriv = gram_entry_dL(i, j, box, arb, acb, prime_powers=prime_powers,
                          options=options)
    slope = max(abs(float(deriv.lower())), abs(float(deriv.upper())))
    record["form"] = "centred_mean_value_assembled"
    record["radius"] = rad
    record["assembled_derivative_enclosure"] = [
        repr(float(deriv.lower())), repr(float(deriv.upper()))
    ]
    return centre + arb(0, rad * slope), record


def block_centred(box, arb, acb, *, prime_powers=None, options=None) -> Dict[str, Any]:
    """``G00``, ``G0b``, ``Gbb``, ``O1``, ``E2``, ``det`` on an ``L``-box."""
    if prime_powers is None:
        mid = float(box.mid()) if hasattr(box, "mid") else float(box)
        prime_powers = WE.prime_powers_below(mid)
    entries: Dict[str, Any] = {}
    records: Dict[str, Any] = {}
    for key, (i, j) in (("G00", ("one", "one")), ("G0b", ("one", "b")),
                        ("Gbb", ("b", "b")), ("O1", ("q1", "q1"))):
        val, rec = gram_entry_centred(i, j, box, arb, acb,
                                      prime_powers=prime_powers, options=options)
        entries[key] = val
        records[key] = rec
    entries["E2"] = entries["G00"] * entries["Gbb"] - entries["G0b"] ** 2
    entries["det_deg2"] = entries["O1"] * entries["E2"]
    entries["_quadrature"] = records
    return entries
