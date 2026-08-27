"""Canonical adjudicated pole primitive (ATLAS-RH-ENG-004 §1).

This module is the **single** implementation of the pole contribution to the
finite Weil Gram. WO-RH-17 adjudicated it as **Candidate A**, derived from the
explicit formula:

    G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i^+ E_j^- + E_i^- E_j^+,
    E_i^± = int_0^L h_i(x) e^{±x/2} dx.

Candidate B -- the legacy ``(sqrt(3)/2)(v+v+^T + v-v-^T)`` block -- is **not**
implemented here and must not be implemented here. It equals Candidate A times
``(sqrt(3)/2)cosh(L/2)``, a factor that is 1 only at ``L = log 3``; it therefore
cannot arise from a change of basis (a constant congruence) and is a calibration
fitted at one test point. It survives only in ``rejected_pole.py``, an archival
module that production code is forbidden to import (see
``tests/test_production_imports.py``).

Carriers
--------
Every routine is carrier-generic: it works with ``float``, ``mpmath``, ``sympy``
and ``flint.arb`` arguments. On an interval carrier the result is an **outward**
enclosure -- no midpoint is ever taken to narrow a value. The one place where a
midpoint is consulted (``_mag``) only chooses between two branches that are both
valid enclosures, so the choice cannot invalidate the result.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from math import factorial

import basis_algebra
from typing import Any, Dict, List, Sequence, Tuple

POLE_FORMULA = (
    "G0_ij = E_i^+ E_j^- + E_i^- E_j^+, E_i^± = int_0^L h_i(x) e^{±x/2} dx"
)
POLE_CANDIDATE = "A"
POLE_STATUS = "ADOPTED_WO_RH_17"

# Basis elements on [0, L], as monomial coefficient tuples (c0, c1, ...).
BASIS_NAMES: Tuple[str, ...] = ("one", "q1", "b", "b3", "b2", "bcube", "bquart")

#: Parity about the cell midpoint x = L/2. Drives ``E^- = ± e^{-L/2} E^+``.
BASIS_PARITY: Dict[str, str] = {"one": "even", "q1": "odd", "b": "even",
                                "b3": "odd", "b2": "even", "bcube": "even",
                                "bquart": "even"}

# Below this |a*L| the endpoint closed form cancels and the series branch (with a
# rigorously bounded remainder) is used instead.
#
# This was 1/1024, which is far too conservative once ``a`` is complex. The
# endpoint form carries terms ~1/a^{k+1}: for the bubble basis (n = 2) that is
# 1/(iz)^3, which at |z| ~ 1e-3 is ~1e9 while the answer is O(L^3). With an exact
# L the cancellation is harmless at working precision, but on an L-**ball** each
# term's width is amplified by that same 1e9 factor and the enclosure explodes --
# measured widths of 4e3 for a Gbb entry whose true value is 3.5e-4.
#
# The series has no cancellation at all and its remainder bound
# |z|^M e^{|z|}/(M!(n+M+1)) is valid for |z| <= 1 (M = 24 gives ~6e-27 there), so
# the cutoff belongs at 1, not 1/1024. Above it the endpoint form's worst term is
# 1/|a|^{n+1} <= 1, which is benign.
_SERIES_CUTOFF = 1.0
_SERIES_TERMS = 24


def basis_coeffs(name: str, L: Any) -> Tuple[Any, ...]:
    """Monomial coefficients of a basis element on ``[0, L]``."""
    if name == "one":
        return (1,)
    if name == "q1":
        return (-L / 2, 1)
    if name == "b":  # x(L - x) = L*x - x^2
        return (0 * L, L, -1)
    if name == "b3":
        # ENG-006 §7: b3(x) = x(L-x)(x-L/2) = -x^3 + (3L/2)x^2 - (L^2/2)x.
        # Odd about x = L/2, like q1, so it joins the odd block.
        return (0 * L, -L * L / 2, 3 * L / 2, 0 * L - 1)
    if name == "b2":
        # ENG-008 §WO-RH-47: b2(x) = b(x)^2 = x^2(L-x)^2 = L^2 x^2 - 2L x^3 + x^4.
        # Even about x = L/2, and the third element of the even sector: in
        # u = x - L/2 the even sector is span{1, u^2, u^4} and b2 supplies u^4.
        # It is the first basis element that is *quadratic* in L, which is why
        # the second-derivative machinery below had to be generalized.
        return (0 * L, 0 * L, L * L, -2 * L, 0 * L + 1)
    if name == "bcube":
        # ENG-009 §WO-RH-62: bcube(x) = b(x)^3 = x^3(L-x)^3
        #                             = L^3 x^3 - 3L^2 x^4 + 3L x^5 - x^6.
        # Even about x = L/2; the fourth even-sector element (u^6 in
        # u = x - L/2), cubic in L. E0-prepared for ENG-010; no production
        # certificate uses it yet.
        L2 = L * L
        return (0 * L, 0 * L, 0 * L, L2 * L, -3 * L2, 3 * L, 0 * L - 1)
    if name == "bquart":
        # ENG-011 §WO-RH-76: bquart(x) = b(x)^4 = x^4 (L - x)^4. Even; the
        # fifth even-sector element (u^8), quartic in L. E0-prepared; no
        # production certificate uses it until the even5 work below.
        L2 = L * L
        return (0 * L, 0 * L, 0 * L, 0 * L, L2 * L2, -4 * L2 * L, 6 * L2,
                -4 * L, 0 * L + 1)
    raise KeyError(f"unknown basis element {name!r}")


def basis_parity(name: str) -> str:
    """Parity of ``h_name`` about ``x = L/2``."""
    try:
        return BASIS_PARITY[name]
    except KeyError:
        raise KeyError(f"unknown basis element {name!r}") from None


# --------------------------------------------------------------------------- #
# Carrier helpers                                                              #
# --------------------------------------------------------------------------- #
def _mag(z: Any) -> float:
    """An **upper bound** on |z| as a float.

    This must over-estimate, never under-estimate: it both selects the branch in
    ``monomial_exp_integral`` and scales that branch's remainder ball. On a ball
    carrier ``float(abs(z))`` returns the *midpoint* of ``|z|``, which for a wide
    ball is far below its supremum -- e.g. 0.5 for a ball spanning ``[0, 1.5]``.
    Using that to size a ``|z|^24`` remainder would understate it by orders of
    magnitude. Harmless for the pole itself, where ``a = ±1/2`` and ``z = ±L/2``
    is never near 0, but the archimedean transforms evaluate this with ``a = it``
    on balls that straddle the origin.
    """
    m = abs(z)
    if hasattr(m, "upper"):  # flint arb (acb.__abs__ also returns arb)
        return float(m.upper())
    if hasattr(m, "mid"):  # pragma: no cover - ball carrier without .upper()
        return abs(float(m.mid())) + abs(float(m.rad()))
    return abs(float(m))


def _exp(z: Any) -> Any:
    if hasattr(z, "exp"):
        return z.exp()
    if isinstance(z, complex):
        import cmath

        return cmath.exp(z)
    from math import exp as _e

    return _e(z)


def _ball(carrier: Any, radius: float) -> Any:
    """A symmetric ball of the given radius on ``carrier``'s type.

    The two-argument constructors differ and getting them confused is silent:
    ``arb(mid, rad)`` builds a ball, but ``acb(re, im)`` builds a *rectangular*
    complex number, so ``acb(0, r)`` is ``r*i`` -- a purely imaginary offset that
    does not enclose a real remainder at all. A complex remainder therefore needs
    ``acb(arb(0, r), arb(0, r))``, a box containing the disc of radius ``r``.
    """
    if hasattr(carrier, "real") and hasattr(carrier, "imag") and hasattr(carrier, "mid"):
        # flint acb: build the ball on each component.
        from flint import arb as _arb

        return type(carrier)(_arb(0, radius), _arb(0, radius))
    if hasattr(carrier, "mid"):  # flint arb
        return type(carrier)(0, radius)
    return 0 * carrier  # exact carriers: the series is truncated exactly below


def _require_interval(L: Any, backend: str | None) -> None:
    if backend is None:
        return
    if backend != "flint":
        raise ValueError(f"backend={backend!r} unsupported; use 'flint' or None")
    if not hasattr(L, "mid"):
        from interval_backend import FlintUnavailable

        raise FlintUnavailable(
            "backend='flint' requires an Arb ball carrier for L; "
            f"got {type(L).__name__}"
        )


# --------------------------------------------------------------------------- #
# int_0^L x^n e^{a x} dx                                                       #
# --------------------------------------------------------------------------- #
def monomial_exp_integral(n: int, a: Any, L: Any) -> Any:
    """Exact ``int_0^L x^n e^{a x} dx``, outward-rounded on interval carriers.

    Two branches, both valid enclosures:

    * ``|a L| > _SERIES_CUTOFF``: the endpoint closed form
      ``e^{ax} sum_k (-1)^k n!/(n-k)! x^{n-k}/a^{k+1}`` evaluated at 0 and L.
      For the pole we always have ``a = ±1/2``, so ``1/a^{k+1} = (∓2)^{k+1}`` is
      benign -- there is no cancellation to fear in this regime.
    * ``|a L| <= _SERIES_CUTOFF``: the everywhere-convergent series
      ``L^{n+1} sum_m (aL)^m / (m! (n+m+1))``, truncated at ``_SERIES_TERMS``
      with an explicit remainder ball. The endpoint form carries terms
      ``~1/a^{k+1}`` and would cancel catastrophically here.
    """
    if a == 0:
        return (L ** (n + 1)) / (n + 1)

    z = a * L
    if _mag(z) <= _SERIES_CUTOFF:
        term = 0 * z + 1
        total = term / (n + 1)
        for m in range(1, _SERIES_TERMS):
            term = term * z / m
            total = total + term / (n + m + 1)
        # |remainder| <= sum_{m>=M} |z|^m/(m!(n+m+1)) <= |z|^M/(M!(n+M+1)) * e^|z|
        M = _SERIES_TERMS
        mag = _mag(z)
        rem = (mag ** M) / (factorial(M) * (n + M + 1)) * 2.72
        total = total + _ball(z, rem)
        return (L ** (n + 1)) * total

    total_L = 0 * z
    for k in range(n + 1):
        coeff = ((-1) ** k) * factorial(n) / factorial(n - k)
        total_L = total_L + coeff * (L ** (n - k)) / (a ** (k + 1))
    at_L = _exp(z) * total_L
    at_0 = ((-1) ** n) * factorial(n) / (a ** (n + 1))
    return at_L - at_0


def poly_exp_integral(coeffs: Sequence[Any], a: Any, L: Any) -> Any:
    """``int_0^L p(x) e^{a x} dx`` for the polynomial with those coefficients."""
    total = 0 * (L * a) if a != 0 else 0 * L
    for n, c in enumerate(coeffs):
        if c == 0:
            continue
        total = total + c * monomial_exp_integral(n, a, L)
    return total


# --------------------------------------------------------------------------- #
# The ENG-004 §1 API                                                           #
# --------------------------------------------------------------------------- #
def _half(L: Any) -> Any:
    """The constant 1/2 on ``L``'s carrier (keeps Arb balls exact)."""
    return (0 * L + 1) / 2 if hasattr(L, "mid") else 0.5


def laplace_plus(basis: str, L: Any, backend: str | None = None) -> Any:
    """``E_basis^+ = int_0^L h_basis(x) e^{+x/2} dx``."""
    _require_interval(L, backend)
    return poly_exp_integral(basis_coeffs(basis, L), _half(L), L)


def laplace_minus(basis: str, L: Any, backend: str | None = None) -> Any:
    """``E_basis^- = int_0^L h_basis(x) e^{-x/2} dx``."""
    _require_interval(L, backend)
    return poly_exp_integral(basis_coeffs(basis, L), -_half(L), L)


def pole_gram_entry(
    basis_i: str, basis_j: str, L: Any, backend: str | None = None
) -> Any:
    """``G0_ij = E_i^+ E_j^- + E_i^- E_j^+`` -- Candidate A, exactly."""
    _require_interval(L, backend)
    half = _half(L)
    ci, cj = basis_coeffs(basis_i, L), basis_coeffs(basis_j, L)
    eip = poly_exp_integral(ci, half, L)
    eim = poly_exp_integral(ci, -half, L)
    ejp = poly_exp_integral(cj, half, L)
    ejm = poly_exp_integral(cj, -half, L)
    return eip * ejm + eim * ejp


def pole_gram_matrix(
    basis: Sequence[str], L: Any, backend: str | None = None
) -> List[List[Any]]:
    """The symmetric Candidate-A pole Gram over ``basis``."""
    _require_interval(L, backend)
    names = list(basis)
    half = _half(L)
    plus = {n: poly_exp_integral(basis_coeffs(n, L), half, L) for n in names}
    minus = {n: poly_exp_integral(basis_coeffs(n, L), -half, L) for n in names}
    return [
        [plus[i] * minus[j] + minus[i] * plus[j] for j in names] for i in names
    ]


# --------------------------------------------------------------------------- #
# Closed forms used by the scalar canary (identical values, cheaper carriers)  #
# --------------------------------------------------------------------------- #
def pole_scalar_g00(L: Any) -> Any:
    """``G0[one,one] = 16 (cosh(L/2) - 1)``.

    ``E_one^± = ±2(e^{±L/2} - 1)`` gives
    ``G0 = 2 E^+ E^- = 8(e^{L/2}-1)(1-e^{-L/2}) = 16(cosh(L/2) - 1)``.
    Algebraically identical to ``pole_gram_entry("one","one",L)``; kept because
    the scalar canary evaluates it in a tight loop and this form avoids the
    polynomial machinery. ``test_pole_primitive`` pins the two together.
    """
    half = L / 2
    ch = half.cosh() if hasattr(half, "cosh") else __import__("math").cosh(half)
    return 16 * (ch - 1)


def pole_scalar_g00_second_derivative(L: Any) -> Any:
    """``d^2/dL^2 G0[one,one] = 4 cosh(L/2) = 2 (r + 1)/sqrt(r)``, ``r = e^L``."""
    half = L / 2
    ch = half.cosh() if hasattr(half, "cosh") else __import__("math").cosh(half)
    return 4 * ch


# --------------------------------------------------------------------------- #
# L-derivative of the pole block                                               #
# --------------------------------------------------------------------------- #
# The basis itself depends on L (``q1 = x - L/2``, ``b = x(L-x)``), so
# differentiating ``E_i^± = int_0^L h_i(x;L) e^{±x/2} dx`` picks up both the
# moving limit and the moving integrand:
#
#     d/dL E_i^± = h_i(L;L) e^{±L/2} + int_0^L (d/dL h_i)(x;L) e^{±x/2} dx.
#
# Keeping this next to the primitive is the point of ENG-004 §1: the jet module
# used to carry its own copy, and that copy was of the *rejected* candidate.
def _horner(coeffs: Sequence[Any], x: Any) -> Any:
    """``p(x)`` from ascending monomial coefficients."""
    out = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        out = out * x + c
    return out


def _derivative_coeffs(coeffs: Sequence[Any]) -> Tuple[Any, ...]:
    """Coefficients of ``dp/dx``. Empty input differentiates to zero."""
    if len(coeffs) <= 1:
        return (coeffs[0] * 0,) if coeffs else ()
    return tuple(n * c for n, c in enumerate(coeffs) if n >= 1)


def basis_at_right_endpoint(name: str, L: Any) -> Any:
    """``h_name(L; L)`` -- the integrand at the moving upper limit.

    Evaluated from an exact ``L``-polynomial that :mod:`basis_algebra` has
    already simplified, not by substituting ``x = L`` on the carrier.

    That distinction is not cosmetic. ``b``, ``b3`` and ``b2`` all vanish
    identically at ``x = L``, and letting the cancellation happen on an ``L``-ball
    turns an exact zero into a ball of radius comparable to the box: on a box of
    radius 1e-2, ``b(L; L) = L*L - L*L`` came back with radius 2.2e-2. The width
    propagates through every derivative bound built on it, and it moved the
    degree-1 and degree-2 certified bounds when it was first tried. Doing the
    cancellation in exact rational arithmetic instead costs nothing and returns
    exact zero.
    """
    return basis_algebra.evaluate_l_poly(basis_algebra.endpoint_poly(name), L)


def basis_endpoint_dL(name: str, L: Any) -> Any:
    """``(d/dL)[h_name(L; L)]`` -- the *total* derivative at the moving limit.

    By the chain rule this is ``(d_x h + d_L h)(L; L)``. Exact, and identically
    zero for every element of the current basis -- ``one`` is constant and the
    other four vanish at ``x = L`` for every ``L``. The general form is computed
    anyway, since that is an accident of this basis rather than a fact about the
    construction.
    """
    return basis_algebra.evaluate_l_poly(
        basis_algebra.endpoint_total_dL_poly(name), L)


def basis_coeffs_dL(name: str, L: Any) -> Tuple[Any, ...]:
    """Monomial coefficients of ``(d/dL) h_name(x; L)``, exactly."""
    return basis_algebra.basis_coeffs_dL_on(name, L)


def basis_coeffs_d2L(name: str, L: Any) -> Tuple[Any, ...]:
    """Monomial coefficients of ``(d^2/dL^2) h_name(x; L)``, exactly.

    ENG-008 §WO-RH-49. ``one``, ``q1`` and ``b`` are linear in ``L`` so this
    vanishes; ``b3`` is quadratic in ``L`` through its ``-(L^2/2) x`` term, and
    ``b2`` is quadratic through ``L^2 x^2``. Nothing here assumes a third
    derivative would vanish, even though for this basis it does.
    """
    return basis_algebra.basis_coeffs_d2L_on(name, L)


def _laplace_dL(name: str, L: Any, sign: int) -> Any:
    half = _half(L) * sign
    boundary = basis_at_right_endpoint(name, L) * _exp(half * L)
    return boundary + poly_exp_integral(basis_coeffs_dL(name, L), half, L)


def laplace_plus_dL(basis: str, L: Any, backend: str | None = None) -> Any:
    """``d/dL E_basis^+``."""
    _require_interval(L, backend)
    return _laplace_dL(basis, L, 1)


def laplace_minus_dL(basis: str, L: Any, backend: str | None = None) -> Any:
    """``d/dL E_basis^-``."""
    _require_interval(L, backend)
    return _laplace_dL(basis, L, -1)


def pole_gram_entry_dL(
    basis_i: str, basis_j: str, L: Any, backend: str | None = None
) -> Any:
    """``d/dL (E_i^+ E_j^- + E_i^- E_j^+)`` -- Candidate A, differentiated."""
    _require_interval(L, backend)
    eip, eim = laplace_plus(basis_i, L), laplace_minus(basis_i, L)
    ejp, ejm = laplace_plus(basis_j, L), laplace_minus(basis_j, L)
    dip, dim = _laplace_dL(basis_i, L, 1), _laplace_dL(basis_i, L, -1)
    djp, djm = _laplace_dL(basis_j, L, 1), _laplace_dL(basis_j, L, -1)
    return dip * ejm + eip * djm + dim * ejp + eim * djp


# --------------------------------------------------------------------------- #
# Second L-derivative of the pole block                                        #
# --------------------------------------------------------------------------- #
# Differentiating F(L) = int_0^L h(x; L) e^{s x/2} dx twice, with H(L) = h(L; L):
#
#   F'  = H(L) e^{sL/2} + int_0^L h_L(x; L) e^{sx/2} dx
#   F'' = H'(L) e^{sL/2} + (s/2) H(L) e^{sL/2} + h_L(L; L) e^{sL/2}
#         + int_0^L h_LL(x; L) e^{sx/2} dx
#
# ENG-008 §WO-RH-49 replaced a per-element table of closed forms with this, the
# general expression, evaluated from the coefficient tables. The table was not
# wrong -- its four entries are reproduced exactly, and
# ``tests/test_pole_primitive.py`` pins them -- but it had to be extended by
# hand for every new basis element, and it silently raised ``KeyError`` for one
# it did not know. ``b2`` was that element.
#
# The one part worth stating: the integral term is not optional. ``one``, ``q1``
# and ``b`` are linear in ``L`` so their ``h_LL`` vanishes and it drops, which is
# why the old table could write them as pure boundary terms; ``b3`` and ``b2``
# are quadratic in ``L`` and it does not.
def _laplace_d2L(name: str, L: Any, sign: int) -> Any:
    """``d^2/dL^2 E_name^sign``, from the coefficient tables."""
    half = _half(L) * sign
    e = _exp(half * L)
    boundary = (
        basis_endpoint_dL(name, L)
        + half * basis_at_right_endpoint(name, L)
        + _horner(basis_coeffs_dL(name, L), L)
    )
    return boundary * e + poly_exp_integral(basis_coeffs_d2L(name, L), half, L)


def laplace_plus_d2L(basis: str, L: Any, backend: str | None = None) -> Any:
    """``d^2/dL^2 E_basis^+``."""
    _require_interval(L, backend)
    return _laplace_d2L(basis, L, 1)


def laplace_minus_d2L(basis: str, L: Any, backend: str | None = None) -> Any:
    """``d^2/dL^2 E_basis^-``."""
    _require_interval(L, backend)
    return _laplace_d2L(basis, L, -1)


def pole_gram_entry_d2L(
    basis_i: str, basis_j: str, L: Any, backend: str | None = None
) -> Any:
    """``d^2/dL^2 (E_i^+ E_j^- + E_i^- E_j^+)`` -- Candidate A, twice differentiated."""
    _require_interval(L, backend)
    ip, im = laplace_plus(basis_i, L), laplace_minus(basis_i, L)
    jp, jm = laplace_plus(basis_j, L), laplace_minus(basis_j, L)
    dip, dim = _laplace_dL(basis_i, L, 1), _laplace_dL(basis_i, L, -1)
    djp, djm = _laplace_dL(basis_j, L, 1), _laplace_dL(basis_j, L, -1)
    d2ip, d2im = _laplace_d2L(basis_i, L, 1), _laplace_d2L(basis_i, L, -1)
    d2jp, d2jm = _laplace_d2L(basis_j, L, 1), _laplace_d2L(basis_j, L, -1)
    return (d2ip * jm + 2 * dip * djm + ip * d2jm
            + d2im * jp + 2 * dim * djp + im * d2jp)
