"""Rigorous scalar E1 canary under Candidate A (ATLAS-RH-ENG-004 §5).

Certifies a **uniform** positive lower bound for the cutoff-free scalar Weil
entry on the cell ``[log 3, log 4]``, using python-flint/Arb throughout. There is
no mpmath path here: an E1 claim must come from interval arithmetic (§5).

The assembled entry
-------------------
``G00(L) = G0(L) - Gp(L) + Ginf(L)`` with, for ``h = 1`` on ``[0, L]``:

* ``G0(L)  = E^+ E^- + E^- E^+ = 16(cosh(L/2) - 1)``  (Candidate A, exact)
* ``Gp(L)  = sum_{q=p^k, log q < L} (2 log p / sqrt q)(L - log q)``  (exact)
* ``Ginf(L) = (1/pi) int_0^inf h_+(t) (2 - 2 cos(L t))/t^2 dt``

Why the cell is convex -- and why that is Candidate A's doing
-------------------------------------------------------------
``Gp`` is piecewise linear with breakpoints exactly at ``log q``; the cell
``(log 3, log 4)`` contains none in its interior, so ``Gp'' = 0`` there. The
pole term contributes ``G0'' = 4 cosh(L/2)`` exactly. For the archimedean term,
the series ``Re psi(1/4 + it/2) = -gamma + sum_n [1/(n+1) - a_n/(a_n^2+t^2/4)]``
with ``a_n = n + 1/4`` and the standard transform
``int_0^inf cos(Lt) a/(a^2+t^2/4) dt = pi e^{-2aL}`` give, for ``L > 0``,

    Ginf''(L) = (2/pi) int_0^inf h_+(t) cos(Lt) dt
              = -2 sum_{n>=0} e^{-2(n+1/4)L}
              = -2 e^{-L/2}/(1 - e^{-2L}) = -e^{L/2}/sinh(L).

Hence, with ``r = e^L``,

    G00''(L) = 4 cosh(L/2) - e^{L/2}/sinh(L) = 2(r^3 - r - 1)/(sqrt(r)(r^2-1)),

which is **exactly** the repository's E0-algebraic curvature ``W00''`` (see
``core.scalar_curvature`` / ``scalar.w00_second_positive_on_r_interval``), proved
positive on ``[log 3, log 4]`` by an algebraic argument. So ``G00`` is convex on
the cell and the uniform bound needs only point evaluations plus convexity --
no interval-``L`` quadrature of an oscillatory integrand.

This is also an independent check on WO-RH-17: the ``4 cosh(L/2)`` term is
produced by Candidate A's pole and by nothing else. The rejected Candidate B
pole ``(sqrt(3)/2)(E^{+2} + E^{-2})`` has a different second derivative, so it
cannot reproduce the certified E0 curvature. ``test_scalar_canary`` pins both
directions.

The bound
---------
``Ginf`` is evaluated on ``[0, T]`` by Arb's rigorous adaptive integrator
(``acb.integral``) after continuing the integrand analytically, and the discarded
tail is **non-negative**, so dropping it only lowers the bound:

    0 <= R_T(L) <= (4/pi)(h_+(T) + kappa(T))/T                (Lemma T below)

Convexity then turns point enclosures into a global bound: at any interior grid
point ``p`` the tangent lies below ``G00`` everywhere, and convexity brackets the
slope between the neighbouring secants, which are computable from the point
enclosures alone.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pole
from interval_backend import backend_info, interval_box, require_flint, set_precision_bits
from rigorous_integration import panel_schedule, rigorous_panel_integral

CELL_LABEL = ("log(3)", "log(4)")
CLAIM_SCOPE = "finite_dimensional_weil_compression"
RH_PROOF_CLAIM = False

#: Frequency cutoff for the rigorous archimedean quadrature.
DEFAULT_T = 200_000
DEFAULT_PRECISION_BITS = 200

#: Series length for the entire factor phi(u) = (2 - 2 cos u)/u^2 near u = 0.
_PHI_TERMS = 14


# --------------------------------------------------------------------------- #
# Analytic continuation of the integrand                                       #
# --------------------------------------------------------------------------- #
def _phi(u, acb):
    """``(2 - 2 cos u)/u^2`` — entire, ``phi(0) = 1``.

    The quotient form cancels catastrophically (and is a 0/0 ball) near the
    origin, so an alternating series with an explicit remainder ball is used
    there. ``acb.integral`` evaluates on complex balls, hence the branch is
    taken on ``|u|``'s upper bound, which is sound for either branch.
    """
    au = u.abs_upper()
    if au < 1:
        total = acb(0)
        power = acb(1)
        u2 = u * u
        for k in range(1, _PHI_TERMS):
            total += ((-1) ** (k + 1)) * (acb(2) / acb(math.factorial(2 * k))) * power
            power = power * u2
        # Alternating with decreasing terms: |remainder| <= first omitted term.
        rem = float(2 / math.factorial(2 * _PHI_TERMS)) * (au ** (2 * _PHI_TERMS - 2))
        from flint import arb as _arb

        return total + acb(_arb(0, rem), _arb(0, rem))
    return (2 - 2 * u.cos()) / (u * u)


def _h_plus_analytic(z, arb, acb, log_pi):
    """``h_+`` continued analytically: ``(psi(1/4+iz/2)+psi(1/4-iz/2))/2 - log pi``.

    For real ``z`` this is ``Re psi(1/4+iz/2) - log pi`` because ``psi`` is
    real-analytic, and unlike ``Re`` it is analytic in ``z`` -- which
    ``acb.integral`` requires. Its poles sit on the imaginary axis, away from the
    real integration path.
    """
    quarter = acb(arb("0.25"))
    i = acb(0, 1)
    return (( quarter + i * z / 2).digamma() + (quarter - i * z / 2).digamma()) / 2 - log_pi


# --------------------------------------------------------------------------- #
# Components of the assembled entry                                            #
# --------------------------------------------------------------------------- #
def pole_term(L):
    """``G0(L) = 16(cosh(L/2) - 1)`` — Candidate A, via the canonical primitive."""
    return pole.pole_scalar_g00(L)


def prime_powers_below(L_value: float) -> List[Tuple[int, int]]:
    """``(q, p)`` for prime powers ``q = p^k`` with ``log q < L``."""
    cap = int(math.floor(math.exp(L_value)))
    out = []
    for p in range(2, cap + 1):
        if any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        q = p
        while q <= cap and math.log(q) < L_value:
            out.append((q, p))
            q *= p
    return sorted(out)


def prime_term(L, arb, prime_powers: Optional[Sequence[Tuple[int, int]]] = None):
    """``Gp(L) = sum_q (2 log p / sqrt q)(L - log q)`` with ``K_00(a;L) = 2(L-a)``."""
    if prime_powers is None:
        prime_powers = prime_powers_below(float(L))
    total = arb(0)
    for q, p in prime_powers:
        total += 2 * arb(p).log() / arb(q).sqrt() * (L - arb(q).log())
    return total


def arch_term_truncated(L, T, arb, acb, *, options: Optional[Dict[str, Any]] = None,
                        panels=None):
    """``(1/pi) int_0^T h_+(t)(2-2cos(Lt))/t^2 dt`` as a rigorous Arb enclosure.

    Delegates to the canonical panel integrator (ENG-005 §2). This module used to
    carry its own decade-edge schedule, which had a coverage bug: for any ``T``
    that did not exceed its last hard-coded edge it integrated *past* ``T`` --
    e.g. ``T = 20000`` was integrated over ``[0, 100000]``. The shipped canary
    used ``T = 200000``, where the schedule happened to be exact, so the ENG-004
    certificate was unaffected; but nothing had checked that. The canonical
    integrator validates that the schedule covers ``[0, T]`` exactly and refuses
    otherwise.

    Returns ``(value, quadrature_record)``.
    """
    log_pi = arb.pi().log()
    L_a = arb(L) if not hasattr(L, "mid") else L
    L_c = acb(L_a)

    def integrand(z, _analytic):
        return _h_plus_analytic(z, arb, acb, log_pi) * (L_c**2) * _phi(L_c * z, acb)

    total, record = rigorous_panel_integral(integrand, T, acb, panels=panels,
                                            options=options)
    return total.real / arb.pi(), record


# --------------------------------------------------------------------------- #
# Lemma T — the discarded tail is non-negative and small                       #
# --------------------------------------------------------------------------- #
def h_plus_at(t, arb, acb):
    """``h_+(t) = Re psi(1/4 + it/2) - log pi`` (Arb)."""
    return acb(arb("0.25"), arb(t) / 2).digamma().real - arb.pi().log()


#: The bound ``t h_+'(t) <= 1`` looks natural and is **false**. Kept as a named
#: constant so the regression test can assert it is rejected rather than quietly
#: re-adopted (ENG-005 §3).
INVALID_TAIL_ASSUMPTION = "t * h_+'(t) <= 1 for t >= 2"
INVALID_NEAR_T = 2.0


def lemma_A_constant(T: float) -> float:
    """``kappa(T)`` with ``t h_+'(t) <= kappa(T)`` for all ``t >= T`` (Lemma A').

    From the series, ``h_+'(t) = (t/2) sum_n a_n/(a_n^2+c^2)^2`` with
    ``a_n = n+1/4`` and ``c = t/2`` -- every term positive. The summand
    ``f(a) = a/(a^2+c^2)^2`` is unimodal with peak ``f(c/sqrt 3) = 9/(16 sqrt3 c^3)``,
    so sampling at spacing 1 gives

        sum_{n>=0} f(a_n) <= int_0^inf f + max f = 1/(2c^2) + 9/(16 sqrt3 c^3),

    hence ``t h_+'(t) <= (t^2/2)(2/t^2 + 2.59808/t^3) = 1 + 1.29904/t``, and on
    ``[T, inf)`` that is at most ``1 + 1.3/T``.

    The domain matters. The bound tends to 1 **from above**, and the tempting
    ``t h_+' <= 1`` is false near the low end: numerically ``t h_+'(2) = 1.0601``.
    ENG-005 §3 requires the constant be carried with a derivation covering the
    exact tail domain, which is why ``kappa`` depends on ``T`` rather than being
    replaced by 1. :func:`invalid_assumption_is_rejected` regression-guards it.
    """
    if T <= 0:
        raise ValueError(f"tail domain must start at T > 0, got {T!r}")
    return 1.0 + 1.3 / float(T)


def invalid_assumption_is_rejected(arb, at: float = INVALID_NEAR_T,
                                   terms: int = 20000) -> Dict[str, Any]:
    """Demonstrate that ``t h_+'(t) <= 1`` fails near ``t = 2`` (ENG-005 §3).

    Computes ``t h_+'(t)`` from the positive series with a rigorous truncation
    ball, and reports that it exceeds 1 while remaining under ``kappa(t)``.
    """
    t_a = arb(repr(at))
    c = t_a / 2
    total = arb(0)
    for n in range(terms):
        a = arb(n) + arb("0.25")
        total += a / ((a * a + c * c) ** 2)
    aN = arb(terms) + arb("0.25")
    total += arb(0, float((1 / (2 * (aN * aN + c * c))).upper()))
    value = (t_a * t_a / 2) * total
    lo, hi = float(value.lower()), float(value.upper())
    return {
        "assumption": INVALID_TAIL_ASSUMPTION,
        "t": at,
        "t_h_plus_prime_enclosure": [repr(lo), repr(hi)],
        "exceeds_one": lo > 1.0,
        "kappa_at_t": lemma_A_constant(at),
        "within_kappa": hi <= lemma_A_constant(at),
        "verdict": "REJECTED" if lo > 1.0 else "NOT_REJECTED",
    }


def lemma_A_numeric_check(arb, ts=(2.0, 10.0, 1e3, 1e5, 2e5), terms: int = 4000):
    """Corroborate Lemma A' by summing the positive series directly.

    ``h_+'(t) = (t/2) sum_n a_n/(a_n^2+c^2)^2`` with every term positive, so the
    truncation is one-sided and its remainder is bounded by the tail integral.
    """
    rows = []
    for t in ts:
        t_a = arb(repr(t))
        c = t_a / 2
        total = arb(0)
        for n in range(terms):
            a = arb(n) + arb("0.25")
            total += a / ((a * a + c * c) ** 2)
        aN = arb(terms) + arb("0.25")
        total += arb(0, float((1 / (2 * (aN * aN + c * c))).upper()))
        t_hp = (t_a * t_a / 2) * total
        rows.append({"t": t,
                     "t_h_plus_prime_upper": repr(float(t_hp.upper())),
                     "kappa_bound": repr(lemma_A_constant(t)),
                     "holds": float(t_hp.upper()) <= lemma_A_constant(t)})
    return all(r["holds"] for r in rows), rows


def tail_bound(T, arb, acb):
    """Rigorous ``0 <= R_T(L) <= (4/pi)(h_+(T) + kappa(T))/T`` for every ``L > 0``.

    *Sign.* ``h_+`` is increasing: in the series ``Re psi(1/4+it/2) = -gamma +
    sum_n [1/(n+1) - a_n/(a_n^2+t^2/4)]`` every summand increases with ``t``. So
    ``h_+ >= h_+(T) > 0`` on ``[T, inf)`` once ``h_+(T) > 0`` is checked with Arb,
    and since ``2 - 2cos >= 0`` the discarded tail is non-negative -- dropping it
    can only *lower* the certified bound, which is the direction we need.

    *Size.* Lemma A' gives ``h_+(t) <= h_+(T) + kappa log(t/T)``, so with
    ``0 <= 2-2cos <= 4`` and ``int_T^inf log(t/T)/t^2 dt = 1/T``,

        R_T <= (4/pi) int_T^inf [h_+(T) + kappa log(t/T)]/t^2 dt
             = (4/pi)(h_+(T) + kappa)/T.
    """
    hT = h_plus_at(T, arb, acb)
    if not (float(hT.lower()) > 0):
        raise ValueError(f"h_+({T}) is not certified positive; tail sign unknown")
    return (4 / arb.pi()) * (hT + lemma_A_constant(T)) / arb(T)


# --------------------------------------------------------------------------- #
# Curvature — the convexity engine                                             #
# --------------------------------------------------------------------------- #
def curvature_from_assembly(L):
    """``G00''(L) = 4 cosh(L/2) - e^{L/2}/sinh(L)``, from the three components."""
    half = L / 2
    return 4 * half.cosh() - half.exp() / L.sinh()


def curvature_e0_formula(L):
    """The repository's E0-algebraic ``W00''(L) = 2(r^3-r-1)/(sqrt r (r^2-1))``."""
    r = L.exp()
    return 2 * (r**3 - r - 1) / (r.sqrt() * (r * r - 1))


# --------------------------------------------------------------------------- #
# The certification                                                            #
# --------------------------------------------------------------------------- #
@dataclass
class GridPoint:
    L: float
    lower: float
    upper: float

    def to_dict(self) -> Dict[str, Any]:
        return {"L": repr(self.L), "G00_lower": repr(self.lower), "G00_upper": repr(self.upper)}


@dataclass
class CanaryResult:
    certified_lower_bound: float
    anchor_L: float
    grid: List[GridPoint]
    tangent_bounds: List[Dict[str, Any]]
    tail_upper: float
    T: int
    precision_bits: int
    convexity: Dict[str, Any]
    lemma_A: Dict[str, Any]
    stats: Dict[str, Any] = field(default_factory=dict)


def default_grid(a: float, b: float, *, coarse: int = 9, refine_halfwidth: float = 0.02,
                 refine: int = 4, focus: float = 1.2828) -> List[float]:
    """Uniform cover of the cell plus a refinement around the expected minimizer.

    The refinement is what tightens the slope bracket (its width scales like
    ``curvature * spacing + tail / spacing``); the coarse points make the
    subdivision a genuine cover rather than a single lucky sample. The focus is a
    *starting guess only* -- correctness does not depend on it, because every
    tangent bound is independently valid and the certificate takes their max.
    """
    pts = [a + (b - a) * k / (coarse - 1) for k in range(coarse)]
    for k in range(-refine, refine + 1):
        p = focus + refine_halfwidth * k / refine
        if a < p < b:
            pts.append(p)
    return sorted(set(pts))


def certify_scalar_canary(
    *,
    T: int = DEFAULT_T,
    precision_bits: int = DEFAULT_PRECISION_BITS,
    grid: Optional[Sequence[float]] = None,
    integral_options: Optional[Dict[str, Any]] = None,
) -> CanaryResult:
    """Certify ``inf_{[log3,log4]} G00 >= certified_lower_bound > 0`` with Arb."""
    import scalar as _scalar

    _, arb, acb, _ = require_flint()
    set_precision_bits(precision_bits)

    a, b = math.log(3.0), math.log(4.0)
    pts = list(grid) if grid is not None else default_grid(a, b)

    # --- convexity -------------------------------------------------------- #
    # Positivity of the curvature is the repository's E0 result and is
    # *algebraic*, not numeric: on r in [3,4], r^3-r-1 >= 23 > 0 and
    # sqrt(r)(r^2-1) > 0. Interval evaluation of that expression over the whole
    # cell is far looser than the algebraic argument, so the algebra is the
    # warrant and the interval cover below only pins the two implementations
    # against each other.
    e0_positive, e0_reason = _scalar.w00_second_positive_on_r_interval(3.0, 4.0)
    if not e0_positive:
        raise ValueError(f"E0 curvature not certified positive: {e0_reason}")

    panels = 32
    curvature_rows = []
    worst_lower = None
    for k in range(panels):
        lo = a + (b - a) * k / panels
        hi = a + (b - a) * (k + 1) / panels
        ball = interval_box(lo, hi)
        assembled = curvature_from_assembly(ball)
        e0 = curvature_e0_formula(ball)
        if not (assembled - e0).contains(0):
            raise ValueError(
                f"assembled curvature disagrees with the E0 formula on [{lo},{hi}] — "
                "ENG-004 §14 stop condition (scalar geometry contradicts the "
                "adjudicated derivation)"
            )
        lower = float(assembled.lower())
        worst_lower = lower if worst_lower is None else min(worst_lower, lower)
        curvature_rows.append({"L_lo": lo, "L_hi": hi, "assembled_lower": lower})
    if worst_lower is None or worst_lower <= 0:
        raise ValueError("assembled curvature not positive on every cover panel")

    # --- tail lemma -------------------------------------------------------- #
    lemA_ok, lemA_rows = lemma_A_numeric_check(arb)
    if not lemA_ok:
        raise ValueError("Lemma A' numeric corroboration failed")
    tail = tail_bound(T, arb, acb)
    tail_up = float(tail.upper())

    # --- point enclosures --------------------------------------------------- #
    # The prime set is constant across the open cell: its breakpoints are exactly
    # the endpoints log 3 and log 4, which is what makes Gp'' vanish inside.
    primes = prime_powers_below((a + b) / 2)
    gridpoints: List[GridPoint] = []
    panel_radius_max = 0.0
    for L in pts:
        L_a = arb(L)
        arch, qrecord = arch_term_truncated(L_a, T, arb, acb, options=integral_options)
        panel_radius_max = max(panel_radius_max, qrecord["max_panel_radius"])
        g = pole_term(L_a) - prime_term(L_a, arb, primes) + arch
        # The discarded tail is non-negative, so ``g`` already lower-bounds G00.
        gridpoints.append(
            GridPoint(L=L, lower=float(g.lower()), upper=float(g.upper()) + tail_up)
        )

    # --- convexity turns point enclosures into a global bound ---------------- #
    # At an interior grid point p the tangent lies below G00 everywhere, and
    # convexity brackets G00'(p) between the neighbouring secants. Each k yields
    # an independently valid global bound; the certificate takes the best.
    tangent_bounds: List[Dict[str, Any]] = []
    best: Optional[Tuple[float, float]] = None
    for k in range(1, len(gridpoints) - 1):
        prev, cur, nxt = gridpoints[k - 1], gridpoints[k], gridpoints[k + 1]
        s_lo = (cur.lower - prev.upper) / (cur.L - prev.L)
        s_hi = (nxt.upper - cur.lower) / (nxt.L - cur.L)
        offsets = (a - cur.L, b - cur.L)
        worst = min(sl * d for sl in (s_lo, s_hi) for d in offsets)
        bound = cur.lower + worst
        tangent_bounds.append(
            {
                "L": repr(cur.L),
                "slope_lower": repr(s_lo),
                "slope_upper": repr(s_hi),
                "tangent_lower_bound": repr(bound),
            }
        )
        if best is None or bound > best[0]:
            best = (bound, cur.L)

    if best is None:
        raise ValueError("grid too small: need at least three points")

    observed_min = min(gp.lower for gp in gridpoints)
    argmin = min(gridpoints, key=lambda gp: gp.lower).L

    return CanaryResult(
        certified_lower_bound=best[0],
        anchor_L=best[1],
        grid=gridpoints,
        tangent_bounds=tangent_bounds,
        tail_upper=tail_up,
        T=T,
        precision_bits=precision_bits,
        convexity={
            "identity": "G00'' = 4cosh(L/2) - e^{L/2}/sinh(L) = 2(r^3-r-1)/(sqrt(r)(r^2-1))",
            "e0_algebraic_reason": e0_reason,
            "pole_second_derivative": "4 cosh(L/2)  (Candidate A)",
            "prime_second_derivative": "0 (piecewise linear; breakpoints are the cell endpoints)",
            "arch_second_derivative": "-e^{L/2}/sinh(L)  (cosine transform of h_+)",
            "cover_panels": panels,
            "assembled_curvature_min_lower": worst_lower,
            "assembled_matches_e0_on_every_panel": True,
            "evidence_class": "E0",
        },
        lemma_A={
            "statement": "t h_+'(t) <= 1 + 1.29904/t, hence <= kappa(T) = 1 + 1.3/T on [T, inf)",
            "kappa": lemma_A_constant(T),
            "numeric_corroboration": lemA_rows,
        },
        stats={
            "grid_points": len(gridpoints),
            "tangent_certificates": len(tangent_bounds),
            "quadrature": qrecord,
            "quadrature_panels_per_point": qrecord["n_panels"],
            "quadrature_panel_radius_max": panel_radius_max,
            "prime_powers_in_cell": [q for q, _ in primes],
            "observed_min_lower_on_grid": observed_min,
            "observed_argmin_on_grid": argmin,
            "integral_options": dict(integral_options or {}),
            "backend": backend_info(precision_bits).to_dict(),
        },
    )
