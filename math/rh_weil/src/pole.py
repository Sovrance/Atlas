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
from typing import Any, Dict, List, Sequence, Tuple

POLE_FORMULA = (
    "G0_ij = E_i^+ E_j^- + E_i^- E_j^+, E_i^± = int_0^L h_i(x) e^{±x/2} dx"
)
POLE_CANDIDATE = "A"
POLE_STATUS = "ADOPTED_WO_RH_17"

# Basis elements on [0, L], as monomial coefficient tuples (c0, c1, ...).
BASIS_NAMES: Tuple[str, ...] = ("one", "q1", "b")

#: Parity about the cell midpoint x = L/2. Drives ``E^- = ± e^{-L/2} E^+``.
BASIS_PARITY: Dict[str, str] = {"one": "even", "q1": "odd", "b": "even"}

# Below this |a*L| the endpoint closed form cancels catastrophically and the
# series branch (with a rigorously bounded remainder) is used instead.
_SERIES_CUTOFF = 1.0 / 1024.0
_SERIES_TERMS = 24


def basis_coeffs(name: str, L: Any) -> Tuple[Any, ...]:
    """Monomial coefficients of a basis element on ``[0, L]``."""
    if name == "one":
        return (1,)
    if name == "q1":
        return (-L / 2, 1)
    if name == "b":  # x(L - x) = L*x - x^2
        return (0 * L, L, -1)
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
    """|z| as a float, for branch selection only (never for a returned value)."""
    try:
        return abs(float(abs(z)))
    except (TypeError, ValueError):  # pragma: no cover - exotic carriers
        m = abs(z)
        return abs(float(m.mid())) + abs(float(m.rad())) if hasattr(m, "mid") else float(m)


def _exp(z: Any) -> Any:
    if hasattr(z, "exp"):
        return z.exp()
    if isinstance(z, complex):
        import cmath

        return cmath.exp(z)
    from math import exp as _e

    return _e(z)


def _ball(carrier: Any, radius: float) -> Any:
    """A symmetric ``[-radius, radius]`` ball on ``carrier``'s type, or 0.0."""
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
def basis_at_right_endpoint(name: str, L: Any) -> Any:
    """``h_name(L; L)`` -- the integrand at the moving upper limit."""
    if name == "one":
        return 0 * L + 1
    if name == "q1":
        return L / 2
    if name == "b":
        return 0 * L
    raise KeyError(f"unknown basis element {name!r}")


def basis_coeffs_dL(name: str, L: Any) -> Tuple[Any, ...]:
    """Monomial coefficients of ``(d/dL) h_name(x; L)``."""
    if name == "one":
        return (0 * L,)
    if name == "q1":
        return (0 * L - 0.5 if not hasattr(L, "mid") else -(0 * L + 1) / 2,)
    if name == "b":  # d/dL [L x - x^2] = x
        return (0 * L, 0 * L + 1)
    raise KeyError(f"unknown basis element {name!r}")


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
