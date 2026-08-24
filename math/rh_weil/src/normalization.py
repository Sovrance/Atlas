"""Frozen normalization for the finite Weil quadratic form (WO-RH-17).

This module is the single mathematical source of truth for the finite Weil
Gram assembly. Every convention is stated explicitly; nothing here is fitted.

Adjudication summary (see ``docs/NORMALIZATION_ADJUDICATION_v0.1.md``):

* **Candidate A (ADOPTED)** — the pole block obtained from the explicit formula,

      G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i^+ E_j^- + E_i^- E_j^+,
      E_i^± = int_0^L h_i(x) e^{±x/2} dx.

* **Candidate B (REJECTED)** — the legacy even block ``(sqrt(3)/2)(v+v+^T +
  v-v-^T)``. On the even sector it equals Candidate A multiplied by the
  *L-dependent* factor ``(sqrt(3)/2)cosh(L/2)``, which is 1 **only** at
  ``L = log 3``. It is therefore a multiplicative calibration fitted at a single
  test point, not a normalization convention, and is retained here only so the
  audit is reproducible.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from math import cosh, exp, factorial, log, sqrt
from typing import Any, Callable, Dict, Sequence, Tuple

NORMALIZATION_VERSION = "v0.1"

# --------------------------------------------------------------------------- #
# 1. Conventions (all seven axes the work order requires to be pinned)         #
# --------------------------------------------------------------------------- #
FOURIER_CONVENTION = "Fhat(xi) = int_R F(x) e^{-i*xi*x} dx"
TILDE_CONVENTION = "tilde_h(x) = conj(h(-x)); real basis => tilde_h(x) = h(-x)"
CONVOLUTION_CONVENTION = "F_ij = h_i * tilde_h_j, supported on [-L, L]"
POLE_FORMULA = (
    "G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i^+ E_j^- + E_i^- E_j^+, "
    "E_i^± = int_0^L h_i(x) e^{±x/2} dx"
)
PRIME_FORMULA = (
    "Gp_ij = sum_{q=p^k, log q < L} (log p / sqrt(q)) * K_ij(log q; L), "
    "K_ij(a;L) = int_0^{L-a} [h_i(x)h_j(x+a) + h_j(x)h_i(x+a)] dx"
)
ARCH_FORMULA = (
    "Ginf_ij(T) = (1/pi) int_0^T h_+(t) Re(conj(H_i(t;L)) H_j(t;L)) dt, "
    "one-sided frequency integral with the even/real basis"
)
ASSEMBLY = "G = G0 - Gp + Ginf"
FREQUENCY_INTEGRAL = "one-sided (0,T], real part taken; even basis"
BASIS_NORMALIZATION = "unnormalised polynomials on [0,L]: 1, q1 = x - L/2, b = x(L-x)"

CLAIM_BOUNDARY = "finite_block_only_no_rh_proof"
RH_PROOF_CLAIM = False

BASIS_NAMES: Tuple[str, ...] = ("one", "q1", "b")


def basis_coeffs(name: str, L: Any) -> Tuple[Any, ...]:
    """Monomial coefficients (c0, c1, ...) of a basis element on [0, L]."""
    if name == "one":
        return (1,)
    if name == "q1":
        return (-L / 2, 1)
    if name == "b":  # x(L-x) = L*x - x^2
        return (0 * L, L, -1)
    raise KeyError(f"unknown basis element {name!r}")


def basis_parity(name: str) -> str:
    """Parity about the cell midpoint x = L/2 (drives the E^- relation)."""
    return {"one": "even", "b": "even", "q1": "odd"}[name]


# --------------------------------------------------------------------------- #
# 2. E^± integrals — exact closed form for any polynomial basis element        #
# --------------------------------------------------------------------------- #
def monomial_exp_integral(n: int, a: Any, L: Any) -> Any:
    """Exact ``int_0^L x^n e^{a x} dx`` for any (possibly complex) exponent ``a``.

    Closed form ``int x^n e^{ax} = e^{ax} sum_k (-1)^k n!/(n-k)! x^{n-k}/a^{k+1}``
    evaluated at the endpoints. Works for float / mpmath / flint-arb carriers, and
    for complex ``a`` (used by the Fourier-side ``H_i(t;L)``). The removable
    ``a -> 0`` case integrates the monomial directly.
    """
    if a == 0:
        return (L ** (n + 1)) / (n + 1)
    # Stability: the endpoint closed form carries terms ~ 1/a^{k+1}, so for small
    # |a*L| it cancels catastrophically (the result is only O(L^{n+1})). Use the
    # everywhere-convergent series there instead:
    #     int_0^L x^n e^{ax} dx = L^{n+1} * sum_m (aL)^m / (m! (n+m+1)).
    z = a * L
    if _abs_like(z) <= 1.0:
        term = 0 * z + 1.0
        total = term / (n + 1)
        for m in range(1, 64):
            term = term * z / m
            delta = term / (n + m + 1)
            total = total + delta
            if _abs_like(delta) <= 1e-24 * max(_abs_like(total), 1e-300):
                break
        return (L ** (n + 1)) * total
    total_L = 0 * z
    for k in range(n + 1):
        coeff = ((-1) ** k) * factorial(n) / factorial(n - k)
        total_L = total_L + coeff * (L ** (n - k)) / (a ** (k + 1))
    at_L = _exp_like(z) * total_L
    at_0 = ((-1) ** n) * factorial(n) / (a ** (n + 1))
    return at_L - at_0


def _abs_like(z: Any) -> float:
    """|z| as a float for float / complex / mpmath / arb carriers."""
    try:
        return float(abs(z))
    except TypeError:  # pragma: no cover - arb balls
        return float(abs(z).mid())


def _exp_like(z: Any) -> Any:
    """exp(z) for float / complex / mpmath / flint-arb style carriers."""
    if hasattr(z, "exp"):
        return z.exp()
    if isinstance(z, complex):
        import cmath

        return cmath.exp(z)
    return exp(z)


def poly_exp_integral(coeffs: Sequence[Any], a: Any, L: Any) -> Any:
    """``int_0^L p(x) e^{a x} dx`` for the polynomial with the given coefficients."""
    total = 0 * (L * a) if a != 0 else 0 * L
    for n, c in enumerate(coeffs):
        if c == 0:
            continue
        total = total + c * monomial_exp_integral(n, a, L)
    return total


def E_pm(name: str, L: Any, sign: int) -> Any:
    """``E^±_i = int_0^L h_i(x) e^{±x/2} dx`` in closed form."""
    if sign not in (1, -1):
        raise ValueError("sign must be +1 or -1")
    half = (0 * L + 1) / 2 if hasattr(L, "mid") else 0.5
    return poly_exp_integral(basis_coeffs(name, L), sign * half, L)


def H_transform(name: str, t: Any, L: Any) -> Any:
    """``H_i(t;L) = int_0^L h_i(x) e^{i t x} dx`` in closed form (complex)."""
    return poly_exp_integral(basis_coeffs(name, L), 1j * t, L)


# --------------------------------------------------------------------------- #
# 3. Candidate A — the ADOPTED pole block (derived from the explicit formula)  #
# --------------------------------------------------------------------------- #
def pole_entry(i: str, j: str, L: Any) -> Any:
    """``G0_ij = E_i^+ E_j^- + E_i^- E_j^+`` (Candidate A, adopted)."""
    Eip, Eim = E_pm(i, L, 1), E_pm(i, L, -1)
    Ejp, Ejm = E_pm(j, L, 1), E_pm(j, L, -1)
    return Eip * Ejm + Eim * Ejp


# --------------------------------------------------------------------------- #
# 4. Candidate B — the REJECTED legacy even block, kept for audit only         #
# --------------------------------------------------------------------------- #
LEGACY_POLE_SCALE = "sqrt(3)/2"
LEGACY_STATUS = "REJECTED_FITTED_CALIBRATION"


def legacy_pole_entry(i: str, j: str, L: Any) -> Any:
    """``(sqrt(3)/2)(v+v+^T + v-v-^T)`` entry — REJECTED; audit use only."""
    s = sqrt(3) / 2 if isinstance(L, float) else _sqrt3_over_2(L)
    Eip, Eim = E_pm(i, L, 1), E_pm(i, L, -1)
    Ejp, Ejm = E_pm(j, L, 1), E_pm(j, L, -1)
    return s * (Eip * Ejp + Eim * Ejm)


def _sqrt3_over_2(carrier: Any) -> Any:
    three = 0 * carrier + 3
    return (three.sqrt() / 2) if hasattr(three, "sqrt") else (sqrt(3) / 2)


def legacy_over_adopted_ratio(L: Any) -> Any:
    """Even-sector ratio ``B/A = (sqrt(3)/2) cosh(L/2)``.

    Derivation: on the even sector ``h(L-x) = h(x)`` gives ``E^- = e^{-L/2}E^+``,
    hence ``A = 2 e^{-L/2} v+ v+^T`` and ``B = (sqrt(3)/2)(1 + e^{-L}) v+ v+^T``.
    The quotient is ``(sqrt(3)/4)(e^{L/2} + e^{-L/2}) = (sqrt(3)/2)cosh(L/2)``,
    which equals 1 exactly at ``L = log 3`` and nowhere else.
    """
    if isinstance(L, float):
        return (sqrt(3) / 2) * cosh(L / 2)
    half = L / 2
    ch = half.cosh() if hasattr(half, "cosh") else cosh(float(half))
    return _sqrt3_over_2(L) * ch


CALIBRATION_FIXED_POINT = "L = log 3"


# --------------------------------------------------------------------------- #
# 5. Prime block — K_ij kernel (shared real-space correlation)                 #
# --------------------------------------------------------------------------- #
def _poly_mul(p: Sequence[Any], q: Sequence[Any]) -> list:
    out = [0 * (p[0] + q[0])] * (len(p) + len(q) - 1)
    out = list(out)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] = out[i + j] + a * b
    return out


def _poly_shift(p: Sequence[Any], a: Any) -> list:
    """Coefficients of p(x + a)."""
    from math import comb

    out = [0 * (p[0] + a)] * len(p)
    out = list(out)
    for n, c in enumerate(p):
        for r in range(n + 1):
            out[r] = out[r] + c * comb(n, r) * (a ** (n - r))
    return out


def _poly_integral(p: Sequence[Any], upper: Any) -> Any:
    total = 0 * upper
    for n, c in enumerate(p):
        total = total + c * (upper ** (n + 1)) / (n + 1)
    return total


def kernel_K(i: str, j: str, a: Any, L: Any) -> Any:
    """``K_ij(a;L) = int_0^{L-a}[h_i(x)h_j(x+a) + h_j(x)h_i(x+a)]dx`` (exact)."""
    hi, hj = basis_coeffs(i, L), basis_coeffs(j, L)
    upper = L - a
    term1 = _poly_integral(_poly_mul(hi, _poly_shift(hj, a)), upper)
    term2 = _poly_integral(_poly_mul(hj, _poly_shift(hi, a)), upper)
    return term1 + term2


def prime_powers_below(L_value: float):
    """Exact ``(q, p, log q)`` for prime powers with ``log q < L``."""
    import math

    c = max(2, int(math.floor(math.exp(float(L_value)) + 1e-12)))
    sieve = [True] * (c + 1)
    if c >= 0:
        sieve[0] = False
    if c >= 1:
        sieve[1] = False
    for i in range(2, int(c**0.5) + 1):
        if sieve[i]:
            for m in range(i * i, c + 1, i):
                sieve[m] = False
    out = []
    for p in (i for i in range(2, c + 1) if sieve[i]):
        q = p
        while q <= c:
            lq = math.log(q)
            if lq < float(L_value):
                out.append((q, p, lq))
            if q > c // p:
                break
            q *= p
    return sorted(out)


def prime_entry(i: str, j: str, L: Any) -> Any:
    """``Gp_ij = sum_q (log p / sqrt q) K_ij(log q; L)``."""
    total = 0 * L
    for q, p, _ in prime_powers_below(float(L)):
        w = log(p) / sqrt(q)
        total = total + w * kernel_K(i, j, log(q), L)
    return total


# --------------------------------------------------------------------------- #
# 6. Identity used by the independent real-space pole route                    #
# --------------------------------------------------------------------------- #
POLE_REALSPACE_IDENTITY = (
    "G0_ij = int_0^L K_ij(a;L) * 2*cosh(a/2) da  "
    "(same K_ij as the prime block; equals E_i^+E_j^- + E_i^-E_j^+)"
)


def normalization_content() -> Dict[str, Any]:
    """The frozen, content-addressable normalization definition."""
    return {
        "normalization_version": NORMALIZATION_VERSION,
        "fourier_convention": FOURIER_CONVENTION,
        "tilde_convention": TILDE_CONVENTION,
        "convolution_convention": CONVOLUTION_CONVENTION,
        "pole_formula": POLE_FORMULA,
        "pole_realspace_identity": POLE_REALSPACE_IDENTITY,
        "prime_formula": PRIME_FORMULA,
        "archimedean_formula": ARCH_FORMULA,
        "assembly": ASSEMBLY,
        "frequency_integral": FREQUENCY_INTEGRAL,
        "basis_normalization": BASIS_NORMALIZATION,
        "claim_boundary": CLAIM_BOUNDARY,
        "rh_proof_claim": RH_PROOF_CLAIM,
    }


def normalization_id() -> str:
    """Content-addressed id of the frozen normalization."""
    import hashlib
    import json

    blob = json.dumps(normalization_content(), sort_keys=True, separators=(",", ":"))
    return "norm_sha256_" + hashlib.sha256(blob.encode("utf-8")).hexdigest()[:32]


# --------------------------------------------------------------------------- #
# 6. Quarantine registry (WO-RH-17 §3.4) — single source of truth              #
# --------------------------------------------------------------------------- #
# Certificates whose numbers were produced under the REJECTED even pole block.
# Kept here rather than in ``scripts/quarantine_normalization.py`` so that the
# certificate writer can enforce the quarantine at the point of write: a certify
# script that regenerates one of these files must not be able to silently drop
# the marker (the same hazard already fixed for ``work_order_status.json``).
QUARANTINE_STATE = "QUARANTINED_NORMALIZATION_ADJUDICATION"

QUARANTINED_CERTIFICATES: Tuple[str, ...] = (
    "e1_scalar_log3_log4.json",
    "e1_degree1_log3_log4.json",
    "e1_degree2_compact_log3_log4.json",
    "e1_fourier_T84_points.json",
    "e1_fourier_T84_uniform_degree2.json",
)

QUARANTINE_REASON = (
    "Depends on the even pole block under adjudication (WO-RH-17). The repository "
    "block (sqrt(3)/2)(v+v+^T+v-v-^T) equals the explicit-formula pole times "
    "(sqrt(3)/2)cosh(L/2), which is 1 only at L = log 3; values at any other L "
    "carry that factor. Regenerate under the frozen normalization (WO-RH-19/20) "
    "before any promotion."
)


def is_quarantined_certificate(name: str) -> bool:
    """True if ``name`` must carry the WO-RH-17 quarantine marker when written."""
    return name in QUARANTINED_CERTIFICATES


def quarantine_block(body: Dict[str, Any]) -> Dict[str, Any]:
    """Build the ``quarantine`` sub-object, preserving ``body``'s claim as prior state."""
    return {
        "reason": QUARANTINE_REASON,
        "work_order": "WO-RH-17",
        "adjudication_certificate": "normalization_adjudication.json",
        "active_normalization_id": normalization_id(),
        "prior_state": {
            "hard_constraints_certified": body.get("hard_constraints_certified"),
            "status": body.get("status"),
        },
        "pre_quarantine_content_hash": body.get("content_hash"),
        "evidence_class_preserved": body.get("evidence_class"),
        "note": "not relabelled E3 on purpose; the historical claim is preserved as evidence",
        "numerically_unaffected_at": "L = log 3 (the calibration fixed point) — still "
                                     "not promotable until regenerated",
    }
