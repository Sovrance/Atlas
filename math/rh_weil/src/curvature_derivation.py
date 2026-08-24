"""Machine-reproducible derivation of ``Ginf''`` (ATLAS-RH-ENG-005 §1).

ENG-004's scalar canary rests on

    G00''(L) = 4 cosh(L/2) - e^{L/2}/sinh(L) = 2(r^3 - r - 1)/(sqrt(r)(r^2-1)),

whose archimedean half is

    Ginf''(L) = (2/pi) int_0^inf h_+(t) cos(Lt) dt = -e^{L/2}/sinh(L).

ENG-004 derived that by hand and checked it numerically. §1 asks for the
derivation itself to be reproducible, so each algebraic step below is re-derived
by SymPy at run time and its verdict recorded in the certificate.

What is and is not proved here
------------------------------
Five steps are **exact symbolic identities**, re-derived by SymPy on every run:

1. ``Re[1/(n + 1/4 + it/2)] = a_n/(a_n^2 + t^2/4)`` with ``a_n = n + 1/4``,
   giving the digamma series for ``Re psi(1/4 + it/2)``.
2. ``int_0^inf cos(Lt) a/(a^2 + t^2/4) dt = pi e^{-2aL}`` for ``a, L > 0``.
3. the geometric sum, as an exact *finite* partial-sum identity plus the
   elementary limit ``N -> inf`` valid because ``e^{-2L} < 1`` for ``L > 0``.
4. ``-2 e^{-L/2}/(1 - e^{-2L}) = -e^{L/2}/sinh(L)``.
5. ``4cosh(L/2) - e^{L/2}/sinh(L) = 2(r^3-r-1)/(sqrt(r)(r^2-1))``, ``r = e^L``.

One step is **not** machine-verified and is not claimed to be: interchanging the
sum with the cosine transform in

    (2/pi) int_0^inf h_+(t) cos(Lt) dt  =  (2/pi) sum_n int_0^inf [...] cos(Lt) dt.

The integral is not absolutely convergent — it is an oscillatory/distributional
Fourier transform of a logarithmically growing function — so the interchange
needs a summability argument (Abel or Cesaro) that is stated as a hypothesis,
not discharged by SymPy. :data:`INTERCHANGE_HYPOTHESIS` records it, and
:func:`derivation_report` reports ``evidence_class`` E0 for the algebra with the
interchange listed under ``analytic_hypotheses``. Calling the whole thing "exact"
would overclassify it, which §1 forbids.

Independent regression: :func:`series_regression` evaluates the *convergent*
series form ``-2 sum_n e^{-2(n+1/4)L}`` directly and compares it to the closed
form, at high precision, without SymPy.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List

CLAIM = "Ginf''(L) = -e^{L/2}/sinh(L) for L > 0"

INTERCHANGE_HYPOTHESIS = (
    "Termwise Fourier cosine transform of the digamma series. The integral "
    "int_0^inf h_+(t) cos(Lt) dt is only conditionally/distributionally "
    "convergent (h_+ grows like log t), so exchanging sum and integral requires "
    "an Abel or Cesaro summability argument. Stated as a hypothesis; NOT "
    "discharged symbolically here."
)


class SymPyRequired(RuntimeError):
    """The symbolic derivation cannot run without SymPy (ENG-004 §6)."""


def _sympy():
    try:
        import sympy as sp
    except ImportError as exc:  # pragma: no cover - exercised by the CI gate
        raise SymPyRequired(
            "sympy is required for the curvature derivation (ENG-005 §1); "
            "pip install -r math/rh_weil/requirements-rigorous.txt"
        ) from exc
    return sp


# --------------------------------------------------------------------------- #
# The five exact steps                                                         #
# --------------------------------------------------------------------------- #
def symbolic_steps() -> List[Dict[str, Any]]:
    """Re-derive each algebraic step with SymPy and report its verdict."""
    sp = _sympy()
    t, L, a, n = sp.symbols("t L a n", positive=True)
    steps: List[Dict[str, Any]] = []

    # 1. real part of a single digamma term
    a_n = n + sp.Rational(1, 4)
    lhs = sp.re(1 / (a_n + sp.I * t / 2))
    rhs = a_n / (a_n**2 + t**2 / 4)
    steps.append({
        "step": "digamma term real part",
        "statement": "Re[1/(n+1/4+it/2)] = a_n/(a_n^2+t^2/4),  a_n = n+1/4",
        "verified": bool(sp.simplify(lhs - rhs) == 0),
    })

    # 2. the cosine transform. SymPy returns it in sinh/cosh form
    #    (pi(cosh(2La) - sinh(2La))), so both sides are rewritten to exp before
    #    comparison -- simplify() will not bridge the two families on its own.
    integral = sp.integrate(sp.cos(L * t) * a / (a**2 + t**2 / 4), (t, 0, sp.oo))
    closed = sp.pi * sp.exp(-2 * a * L)
    steps.append({
        "step": "cosine transform",
        "statement": "int_0^inf cos(Lt) a/(a^2+t^2/4) dt = pi e^{-2aL}",
        "sympy_result": str(sp.simplify(integral))[:200],
        "verified": bool(sp.simplify((integral - closed).rewrite(sp.exp)) == 0),
    })

    # 3. the geometric sum, as an exact finite identity plus an elementary limit.
    #    SymPy will not evaluate the infinite sum directly (it does not assume
    #    |e^{-2L}| < 1), and asserting the closed form without the ratio
    #    condition would be the kind of unearned step §1 warns about. The finite
    #    partial sum is exact algebra; the limit needs only e^{-2L} < 1, which
    #    holds for every L > 0.
    k, N = sp.symbols("k N", integer=True, nonnegative=True)
    partial = sp.summation(sp.exp(-2 * (k + sp.Rational(1, 4)) * L), (k, 0, N - 1))
    partial_closed = sp.exp(-L / 2) * (1 - sp.exp(-2 * N * L)) / (1 - sp.exp(-2 * L))
    # SymPy returns a Piecewise guarding the degenerate ratio e^{-2L} = 1. That
    # branch is L = 0, excluded here since L > 0, so the generic branch is the
    # whole story on our domain -- take it explicitly rather than leaving a
    # Piecewise that never compares equal to 0.
    generic = partial
    if isinstance(partial, sp.Piecewise):
        generic = partial.args[-1][0]
    elif partial.has(sp.Piecewise):
        generic = sp.piecewise_fold(partial)
        if isinstance(generic, sp.Piecewise):
            generic = generic.args[-1][0]
        else:
            generic = partial.replace(
                lambda e: isinstance(e, sp.Piecewise), lambda e: e.args[-1][0]
            )
    finite_ok = bool(sp.simplify((generic - partial_closed).rewrite(sp.exp)) == 0)
    limit_expr = sp.limit(partial_closed, N, sp.oo)
    limit_ok = bool(
        sp.simplify((limit_expr - sp.exp(-L / 2) / (1 - sp.exp(-2 * L))).rewrite(sp.exp)) == 0
    )
    steps.append({
        "step": "geometric sum over n",
        "statement": ("finite: sum_{k<N} e^{-2(k+1/4)L} = e^{-L/2}(1-e^{-2NL})/(1-e^{-2L}); "
                      "limit N->inf with e^{-2L} < 1 gives e^{-L/2}/(1-e^{-2L})"),
        "finite_sum_verified": finite_ok,
        "limit_verified": limit_ok,
        "degenerate_branch_excluded": "e^{-2L} = 1 requires L = 0, outside L > 0",
        "verified": finite_ok and limit_ok,
    })

    # 4. closed form of Ginf''
    geom = sp.exp(-L / 2) / (1 - sp.exp(-2 * L))
    ginf2 = -2 * geom
    target = -sp.exp(L / 2) / sp.sinh(L)
    steps.append({
        "step": "Ginf'' closed form",
        "statement": "-2 e^{-L/2}/(1-e^{-2L}) = -e^{L/2}/sinh(L)",
        "verified": bool(sp.simplify((ginf2 - target).rewrite(sp.exp)) == 0),
    })

    # 5. the assembled curvature equals the E0 algebraic formula
    r = sp.exp(L)
    assembled = 4 * sp.cosh(L / 2) + target
    e0 = 2 * (r**3 - r - 1) / (sp.sqrt(r) * (r**2 - 1))
    steps.append({
        "step": "assembled curvature",
        "statement": "4cosh(L/2) - e^{L/2}/sinh(L) = 2(r^3-r-1)/(sqrt(r)(r^2-1)), r = e^L",
        "verified": bool(sp.simplify((assembled - e0).rewrite(sp.exp)) == 0),
    })
    return steps


# --------------------------------------------------------------------------- #
# Independent regression, no SymPy                                             #
# --------------------------------------------------------------------------- #
def series_regression(L_values=(math.log(3.0), 1.2, 1.2828, math.log(4.0)),
                      terms: int = 400) -> List[Dict[str, Any]]:
    """Compare ``-2 sum_n e^{-2(n+1/4)L}`` to ``-e^{L/2}/sinh L`` numerically.

    This series converges geometrically, so it is an honest independent check of
    the closed form — it shares no code with the SymPy path above.
    """
    rows = []
    for L in L_values:
        total = 0.0
        for k in range(terms):
            total += math.exp(-2 * (k + 0.25) * L)
        series = -2 * total
        closed = -math.exp(L / 2) / math.sinh(L)
        rows.append({
            "L": repr(L),
            "series": repr(series),
            "closed_form": repr(closed),
            "abs_diff": repr(abs(series - closed)),
            "agrees": abs(series - closed) < 1e-12 * max(1.0, abs(closed)),
        })
    return rows


def derivation_report() -> Dict[str, Any]:
    """The artifact the scalar certificate embeds and the runner gates on."""
    import sympy as sp  # noqa: F401  (SymPyRequired already raised if missing)

    steps = symbolic_steps()
    regression = series_regression()
    all_symbolic = all(s["verified"] for s in steps)
    all_numeric = all(r["agrees"] for r in regression)
    return {
        "claim": CLAIM,
        "symbolic_engine": f"sympy {sp.__version__}",
        "steps": steps,
        "all_symbolic_steps_verified": all_symbolic,
        "independent_series_regression": regression,
        "all_regression_rows_agree": all_numeric,
        # Deliberately E0 for the algebra, with the analytic gap named. The
        # program's evidence ladder is only useful if a step that is argued
        # rather than machine-checked is visible as such.
        "evidence_class": "E0",
        "analytic_hypotheses": [INTERCHANGE_HYPOTHESIS],
        "overclaim_guard": (
            "The five algebraic steps are symbolic identities. The termwise "
            "interchange is a stated hypothesis, not a machine-checked step. This "
            "report is E0 for the algebra only."
        ),
        "status": "VERIFIED" if (all_symbolic and all_numeric) else "FAILED",
        "rh_proof_claim": False,
    }
