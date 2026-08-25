"""ATLAS-RH-ENG-007 §15 (WO-RH-46) -- the first genuinely 3-dimensional parity block.

ENG-006's "degree-3" odd block is 2x2, and that is a real limitation on what
ENG-006 could learn. In two dimensions the trace and determinant already fix the
spectrum, so inertia, the spectral moments and the rank-trace bound all collapse
onto information the determinant carried anyway -- as the ENG-006 information
report says in as many words. To find out whether those channels add anything,
they have to be pointed at a block where the determinant is not the whole story.

This module prepares that block. It does not certify it. §15 is explicit that
ENG-007 runs only exact (E0) identities and a heuristic (E3) conditioning and
topology preview, and promotes no new E1 degree result.

The block
---------
Write ``u = x - L/2``. Then ``b(x) = x(L-x) = L^2/4 - u^2``, so about the
midpoint::

    one = 1                    even    span{1}
    b   = x(L-x)               even    span{1, u^2}
    b2  = x^2 (L-x)^2 = b^2    even    span{1, u^2, u^4}
    q1  = x - L/2              odd     span{u}
    b3  = x(L-x)(x-L/2)        odd     span{u, u^3}

The odd sector through degree 3 is exactly ``span{u, u^3}`` -- two-dimensional,
which is why ENG-006 stopped there and why ``q1^3`` would have added nothing.
Extending the **even** sector with ``b2`` gives the first 3x3 parity block, and
it is the sector worth extending: it contains the constant function, so its Gram
is the Weil form restricted to a three-dimensional space containing the scalar
case, and its 2x2 leading block is exactly the ``E2 = G00 Gbb - G0b^2`` that
ENG-004 and ENG-005 certified.

What this module provides
-------------------------
* exact prime-overlap kernels for the three new pairs, in closed form and as
  coefficient lists in ``a`` (both routes, so each checks the other);
* the pole block for the extended basis, from the same Candidate-A primitive;
* an E3 preview assembly in mpmath -- deliberately *not* interval arithmetic, so
  that nothing here can be mistaken for a certificate;
* the conditioning diagnostics ENG-008 needs in order to choose a basis scaling.

A caution recorded here rather than discovered later
----------------------------------------------------
Every basis element ENG-005 works with is **linear in L**, and ``pole.py``'s
second-derivative machinery is built on that: ``_laplace_d2L`` drops the
second integral because ``d^2_L h_i = 0``. ``b2`` is *quadratic* in ``L``
(``L^2 x^2 - 2L x^3 + x^4``), so that simplification fails for it. Any ENG-008
work that wants an ``E2''``-style curvature argument on this block has to extend
``_laplace_d2L`` first. This module does not use or need it.

No RH proof claim is made here. Nothing in this module is promoted, and the
preview it produces is E3 evidence, which is never a warrant.
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Any, Callable, Dict, List, Sequence, Tuple

CLAIM_SCOPE = "finite_dimensional_weil_compression"
EVIDENCE_NOTE = "E0 identities and an E3 preview only; no E1 result is claimed"

CELL = (math.log(3.0), math.log(4.0))
CELL_LABEL = ("log(3)", "log(4)")

#: The extended even sector, in Gram order.
EVEN_BASIS: Tuple[str, ...] = ("one", "b", "b2")

#: The odd sector, unchanged from ENG-006.
ODD_BASIS: Tuple[str, ...] = ("q1", "b3")

#: Parity about ``x = L/2`` of every element this module knows.
PARITY: Dict[str, str] = {
    "one": "even", "b": "even", "b2": "even", "q1": "odd", "b3": "odd",
}


# --------------------------------------------------------------------------- #
# Basis elements as monomial coefficient lists                                 #
# --------------------------------------------------------------------------- #
def basis_coeffs(name: str, L: Any) -> Tuple[Any, ...]:
    """Monomial coefficients of ``h_name(x; L)``, ascending in ``x``.

    Deliberately a local copy rather than an extension of ``pole.basis_coeffs``.
    Editing that function would change the source hash of ``pole.py``, which
    every promoted E1 certificate binds -- a preparatory pilot must not be able
    to invalidate a certified result.
    """
    one = L * 0 + 1
    zero = L * 0
    if name == "one":
        return (one,)
    if name == "q1":
        return (-L / 2, one)
    if name == "b":  # L x - x^2
        return (zero, L, -one)
    if name == "b3":  # -x^3 + (3L/2) x^2 - (L^2/2) x
        return (zero, -L * L / 2, 3 * L / 2, -one)
    if name == "b2":  # (L x - x^2)^2 = L^2 x^2 - 2L x^3 + x^4
        return (zero, zero, L * L, -2 * L, one)
    raise KeyError(f"unknown basis element {name!r}")


def evaluate(name: str, x: Any, L: Any) -> Any:
    """``h_name(x; L)`` by Horner on the coefficient list."""
    coeffs = basis_coeffs(name, L)
    out = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        out = out * x + c
    return out


def is_even_about_midpoint(name: str, L: Fraction, samples: Sequence[Fraction]) -> bool:
    """Exact parity check: ``h(L - x) == h(x)`` at rational sample points."""
    return all(evaluate(name, L - x, L) == evaluate(name, x, L) for x in samples)


def is_odd_about_midpoint(name: str, L: Fraction, samples: Sequence[Fraction]) -> bool:
    return all(evaluate(name, L - x, L) == -evaluate(name, x, L) for x in samples)


# --------------------------------------------------------------------------- #
# Exact prime-overlap kernels                                                  #
# --------------------------------------------------------------------------- #
# K_fg(a; L) = int_0^{L-a} [ f(x) g(x+a) + g(x) f(x+a) ] dx.
#
# The three new closed forms were derived symbolically and are checked in
# tests/test_pilot3_exact.py against exact rational integration of the
# integrand -- two independent routes, neither trusting the other's algebra.
def kernel_one_b2(a: Any, L: Any) -> Any:
    d = L - a
    return d**3 * (L * L + 3 * L * a + 6 * a * a) / 15


def kernel_b_b2(a: Any, L: Any) -> Any:
    d = L - a
    return d**4 * (3 * L**3 + 12 * L**2 * a + 16 * L * a**2 + 4 * a**3) / 210


def kernel_b2_b2(a: Any, L: Any) -> Any:
    d = L - a
    return d**5 * (L**4 + 5 * L**3 * a + 9 * L**2 * a**2 + 5 * L * a**3 + a**4) / 315


#: Only the pairs this module adds. The existing pairs stay in ``core``.
NEW_KERNELS: Dict[Tuple[str, str], Callable[[Any, Any], Any]] = {
    ("one", "b2"): kernel_one_b2,
    ("b", "b2"): kernel_b_b2,
    ("b2", "b2"): kernel_b2_b2,
}

#: The same three kernels as coefficient lists in ``a`` (ascending), with ``L``
#: symbolic. The real-space archimedean route consumes this form, and having
#: both lets each check the other exactly.
def kernel_coeffs_in_a(i: str, j: str, L: Any) -> List[Any]:
    key = (i, j) if (i, j) in NEW_KERNELS else (j, i)
    zero = L * 0
    if key == ("one", "b2"):
        return [L**5 / 15, zero, zero, -2 * L**2 / 3, L, -(L * 0 + 2) / 5]
    if key == ("b", "b2"):
        return [L**7 / 70, zero, -L**5 / 15, zero, L**3 / 6, -2 * L**2 / 15, zero,
                (L * 0 + 2) / 105]
    if key == ("b2", "b2"):
        return [L**9 / 315, zero, -2 * L**7 / 105, zero, L**5 / 15, -L**4 / 15,
                zero, 2 * L**2 / 105, zero, -(L * 0 + 1) / 315]
    raise KeyError(f"no coefficient expansion for {(i, j)!r}")


def kernel(i: str, j: str, a: Any, L: Any) -> Any:
    """``K_ij(a; L)`` for any pair over ``EVEN_BASIS + ODD_BASIS``."""
    key = (i, j) if (i, j) in NEW_KERNELS else (j, i)
    if key in NEW_KERNELS:
        return NEW_KERNELS[key](a, L)
    import core

    legacy = {
        ("one", "one"): core.kernel_00,
        ("one", "b"): core.kernel_0b,
        ("b", "b"): core.kernel_bb,
        ("q1", "q1"): core.kernel_q1q1,
        ("q1", "b3"): core.kernel_q1b3,
        ("b3", "b3"): core.kernel_b3b3,
    }
    key = (i, j) if (i, j) in legacy else (j, i)
    if key not in legacy:
        raise KeyError(f"no kernel for the pair {(i, j)!r}")
    return legacy[key](a, L)


def kernel_by_quadrature(i: str, j: str, a: Fraction, L: Fraction) -> Fraction:
    """``K_ij(a; L)`` by exact rational integration of the polynomial integrand.

    The independent route. Polynomial coefficients are combined exactly and the
    integral is evaluated termwise, so this involves no floating point and no
    closed form -- if it agrees with :func:`kernel`, the closed form is right.
    """
    ci = [Fraction(c) for c in basis_coeffs(i, L)]
    cj = [Fraction(c) for c in basis_coeffs(j, L)]
    # (x + a)^n expanded once, reused for both orderings.
    def shifted(coeffs: List[Fraction]) -> List[Fraction]:
        out = [Fraction(0)] * len(coeffs)
        for n, c in enumerate(coeffs):
            for r in range(n + 1):
                out[r] += c * Fraction(math.comb(n, r)) * a ** (n - r)
        return out

    prod: Dict[int, Fraction] = {}
    for left, right in ((ci, shifted(cj)), (cj, shifted(ci))):
        for m, cm in enumerate(left):
            for n, cn in enumerate(right):
                prod[m + n] = prod.get(m + n, Fraction(0)) + cm * cn
    d = L - a
    return sum(c * d ** (k + 1) / (k + 1) for k, c in prod.items())


def kernel_from_coeffs(i: str, j: str, a: Fraction, L: Fraction) -> Fraction:
    """``K_ij(a; L)`` evaluated from :func:`kernel_coeffs_in_a`."""
    coeffs = kernel_coeffs_in_a(i, j, Fraction(L))
    out = Fraction(0)
    power = Fraction(1)
    for c in coeffs:
        out += Fraction(c) * power
        power *= a
    return out


# --------------------------------------------------------------------------- #
# Blocks                                                                       #
# --------------------------------------------------------------------------- #
def prime_powers_below(L_value: float) -> List[Tuple[int, int]]:
    """``(q, p)`` for prime powers ``q = p^k`` with ``log q < L``."""
    cap = int(math.floor(math.exp(L_value)))
    out = []
    for p in range(2, cap + 1):
        if any(p % div == 0 for div in range(2, int(p**0.5) + 1)):
            continue
        q = p
        while q <= cap and math.log(q) < L_value:
            out.append((q, p))
            q *= p
    return sorted(out)


def pole_entry(i: str, j: str, L: Any, exp_fn: Callable[[Any], Any]) -> Any:
    """``G0_ij = E_i^+ E_j^- + E_i^- E_j^+`` for the Candidate-A pole.

    ``E_i^± = int_0^L h_i(x) e^{±x/2} dx``, evaluated in closed form from the
    monomial coefficients. Same primitive as ``pole.pole_gram_entry``; the
    coefficients come from this module's own table.
    """
    def laplace(name: str, sign: int) -> Any:
        coeffs = basis_coeffs(name, L)
        s = (L * 0 + sign) / 2
        # int_0^L x^n e^{s x} dx by repeated integration by parts:
        #   I_n = (L^n e^{sL} - n I_{n-1}) / s,  I_0 = (e^{sL} - 1)/s
        e = exp_fn(s * L)
        integrals = [(e - 1) / s]
        for n in range(1, len(coeffs)):
            integrals.append((L**n * e - n * integrals[n - 1]) / s)
        return sum(c * I for c, I in zip(coeffs, integrals))

    return (laplace(i, 1) * laplace(j, -1)) + (laplace(i, -1) * laplace(j, 1))


def prime_entry(i: str, j: str, L: Any, *, log_fn, sqrt_fn,
                prime_powers: Sequence[Tuple[int, int]] | None = None) -> Any:
    """``Gp_ij(L) = sum_q (log p / sqrt q) K_ij(log q; L)``."""
    if prime_powers is None:
        prime_powers = prime_powers_below(float(L))
    total = L * 0
    for q, p in prime_powers:
        total += (log_fn(p) / sqrt_fn(q)) * kernel(i, j, log_fn(q), L)
    return total


def pole_matrix(basis: Sequence[str], L: Any, exp_fn) -> List[List[Any]]:
    names = list(basis)
    return [[pole_entry(i, j, L, exp_fn) for j in names] for i in names]


def cross_parity_pole_entries(L: Any, exp_fn) -> Dict[Tuple[str, str], Any]:
    """Every even/odd pole entry, which parity says must vanish."""
    return {
        (i, j): pole_entry(i, j, L, exp_fn)
        for i in EVEN_BASIS for j in ODD_BASIS
    }


# --------------------------------------------------------------------------- #
# E3 preview (mpmath)                                                          #
# --------------------------------------------------------------------------- #
# Everything below is *heuristic*. It uses mpmath, not Arb, and mpmath never
# certifies anything in this program (ENG-004 §5). The preview exists to answer
# ENG-008's planning questions -- is the block definite on the cell, how badly
# conditioned is it, does a basis rescaling help -- not to establish a bound.
def require_mpmath():
    import mpmath  # noqa: F401

    return mpmath


def h_plus_mp(t, mp):
    """``Re psi(1/4 + i t / 2) - log pi``."""
    return mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * t / 2)) - mp.log(mp.pi)


def fourier_transform_mp(name: str, t, L, mp):
    """``H_i(t; L) = int_0^L h_i(x) e^{i t x} dx``, in closed form.

    Closed form rather than quadrature: the integrand oscillates at rate ``t``
    up to ``t = T``, and quadrature of it inside another quadrature is both slow
    and inaccurate exactly where the preview needs to be trusted least.

    Two branches, because the obvious one is unstable. The recursion

        I_n = (L^n e^{sL} - n I_{n-1}) / s,   I_0 = (e^{sL} - 1)/s

    divides by ``s`` at every step, and for small ``|s|`` the numerator is a
    difference of nearly equal numbers. mpmath's default tanh-sinh quadrature
    clusters nodes exponentially close to the panel endpoints, so ``t`` really
    does reach ``1e-25`` and smaller; at that point the recursion amplifies
    rounding by ``|s|^{-n}`` and returns values around ``1e30``. Below
    ``|sL| = 1/2`` the Taylor series is used instead, where every term is
    positive-definite in magnitude and nothing cancels.
    """
    coeffs = basis_coeffs(name, mp.mpf(L))
    Lm = mp.mpf(L)
    s = 1j * t
    if abs(s) * Lm < mp.mpf("0.5"):
        # I_n = sum_k s^k L^{n+k+1} / (k! (n+k+1))
        integrals = []
        eps = mp.mpf(10) ** (-(mp.mp.dps + 5))
        for n in range(len(coeffs)):
            total = mp.mpc(0)
            term = Lm ** (n + 1)
            k = 0
            while True:
                contrib = term / (n + k + 1)
                total += contrib
                if k > 0 and abs(contrib) < eps * (abs(total) + eps):
                    break
                k += 1
                term = term * s * Lm / k
                if k > 4 * mp.mp.dps + 40:
                    break
            integrals.append(total)
        return sum(c * I for c, I in zip(coeffs, integrals))
    e = mp.e ** (s * Lm)
    integrals = [(e - 1) / s]
    for n in range(1, len(coeffs)):
        integrals.append((Lm ** n * e - n * integrals[n - 1]) / s)
    return sum(c * I for c, I in zip(coeffs, integrals))


def arch_entry_mp(i: str, j: str, L, T: float, mp, *, maxdegree: int = 8):
    """``(1/pi) int_0^T h_+(t) Re(conj(H_i) H_j) dt`` -- E3, mpmath quadrature."""
    def integrand(t):
        hi = fourier_transform_mp(i, t, L, mp)
        hj = fourier_transform_mp(j, t, L, mp)
        return h_plus_mp(t, mp) * mp.re(mp.conj(hi) * hj)

    # Subdivide: h_+ varies slowly but the transforms oscillate, so a single
    # panel over [0, T] would be quadrature theatre.
    panels = max(16, int(T / 4))
    edges = [mp.mpf(T) * k / panels for k in range(panels + 1)]
    total = mp.quad(integrand, edges, maxdegree=maxdegree)
    return total / mp.pi


def gram_matrix_mp(basis: Sequence[str], L: float, T: float, mp, *,
                   prime_powers=None, maxdegree: int = 8) -> List[List[Any]]:
    """``G = G0 - Gp + Ginf`` over ``basis`` at a single real ``L``. E3 only."""
    names = list(basis)
    Lm = mp.mpf(L)
    if prime_powers is None:
        prime_powers = prime_powers_below(float(L))
    out = [[None] * len(names) for _ in names]
    for m, i in enumerate(names):
        for n, j in enumerate(names):
            if n < m:
                out[m][n] = out[n][m]
                continue
            g0 = pole_entry(i, j, Lm, lambda z: mp.e**z)
            gp = prime_entry(i, j, Lm, log_fn=lambda v: mp.log(mp.mpf(v)),
                             sqrt_fn=lambda v: mp.sqrt(mp.mpf(v)),
                             prime_powers=prime_powers)
            gi = arch_entry_mp(i, j, Lm, T, mp, maxdegree=maxdegree)
            out[m][n] = g0 - gp + gi
    return out


def leading_minors(M: Sequence[Sequence[Any]]) -> List[Any]:
    """The leading principal minors, by exact cofactor expansion on the entries.

    Sylvester's criterion reads these; ``AtlasRH.posDef_sym3_iff`` is the proved
    statement that all three positive is equivalent to positive definiteness.
    """
    n = len(M)
    out = []
    for k in range(1, n + 1):
        sub = [[M[a][b] for b in range(k)] for a in range(k)]
        out.append(_det(sub))
    return out


def _det(M: Sequence[Sequence[Any]]):
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    total = None
    for c in range(n):
        minor = [[M[r][b] for b in range(n) if b != c] for r in range(1, n)]
        term = M[0][c] * _det(minor)
        term = term if c % 2 == 0 else -term
        total = term if total is None else total + term
    return total


def condition_report(M: Sequence[Sequence[Any]], mp) -> Dict[str, Any]:
    """Eigenvalues, condition number and diagonal scale spread. E3.

    ``mp.eigsy`` is a floating eigenvalue solver and is therefore banned from
    every rigorous path in this program. It is used here, and only here,
    because a *preview* of conditioning is exactly the sort of question a
    non-rigorous tool is for.
    """
    n = len(M)
    A = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            A[a, b] = M[a][b]
    eig = mp.eigsy(A, eigvals_only=True)
    vals = sorted(float(eig[k]) for k in range(n))
    diag = [float(M[k][k]) for k in range(n)]
    amax = max(abs(v) for v in vals)
    amin = min(abs(v) for v in vals)
    return {
        "eigenvalues": vals,
        "n_positive": sum(1 for v in vals if v > 0),
        "n_negative": sum(1 for v in vals if v < 0),
        "n_zero": sum(1 for v in vals if v == 0),
        "condition_number": (amax / amin) if amin > 0 else None,
        "diagonal": diag,
        "diagonal_spread": (max(diag) / min(diag)) if min(diag) > 0 else None,
    }


def jacobi_rescale(M: Sequence[Sequence[Any]], mp) -> Tuple[List[List[Any]], List[float]]:
    """Symmetric diagonal rescaling ``D M D`` with ``D = diag(1/sqrt(M_kk))``.

    A congruence, so by ``AtlasRH.posIndexAtLeast_congruence_iff`` and
    ``AtlasRH.rank_congruence`` it changes no part of the inertia. That is the
    point: it is a free change of conditioning, and the proved theorem is what
    says it is free. Requires a positive diagonal, which is the case worth
    scaling anyway.
    """
    n = len(M)
    scales = [1.0 / math.sqrt(float(M[k][k])) for k in range(n)]
    out = [[M[a][b] * scales[a] * scales[b] for b in range(n)] for a in range(n)]
    return out, scales
