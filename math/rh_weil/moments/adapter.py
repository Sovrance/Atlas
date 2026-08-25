"""Atlas B1 adapter for Weil spectral moments (§6, WO-RH-31).

Answers the five queries §6 names, each with an explicit status, and refuses to
upgrade a necessary condition into a sufficient one. The queries:

* minimum feasible positive eigenvalue count
* maximum feasible negative eigenvalue count
* bounds on the smallest eigenvalue
* whether the moments force PSD
* whether an observed inertia is uniquely determined by the moments

The last two are where an adapter like this usually goes wrong. A truncated
localizing matrix being PSD is *necessary* for the spectrum to sit in
``[0, inf)`` and not sufficient, so "the localizing matrix is PSD" must be
reported as ``INSUFFICIENT_INFORMATION``, never as "the moments force PSD". The
conclusive direction is the other one: a localizing matrix that provably fails
to be PSD proves a negative eigenvalue exists.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

_PROGRAM = Path(__file__).resolve().parents[1]
if str(_PROGRAM) not in sys.path:
    sys.path.insert(0, str(_PROGRAM))

from .feasible_spectrum import (  # noqa: E402
    CONCLUSIVE,
    INSUFFICIENT,
    b1_available,
    eigenvalue_bounds_from_two_moments,
    hankel_matrix,
    localizing_matrix,
    psd_status,
    spectrum_from_two_moments_n2,
)
from .spectral_moments import KIND_SPECTRAL_MOMENT, moment_report, spectral_moments  # noqa: E402


def _moment_sequence(n: int, moments: Dict[str, Any]) -> List[Any]:
    """``[m_0, m_1, ...]`` with ``m_0 = n``, the count of eigenvalues."""
    seq: List[Any] = [n]
    k = 1
    while f"m{k}" in moments:
        seq.append(moments[f"m{k}"])
        k += 1
    return seq


def query_forces_psd(seq: Sequence[Any]) -> Dict[str, Any]:
    """Do these moments prove the spectrum is non-negative?

    Almost always no, and saying so is the point. The localizing matrix test can
    *refute* non-negativity but cannot establish it from finitely many moments
    without a flat extension.
    """
    if len(seq) < 4:
        return {"query": "moments_force_psd", "status": INSUFFICIENT,
                "reason": "need m_0..m_3 to form the degree-1 localizing matrix"}
    L = localizing_matrix(seq, 1)
    st = psd_status(L)
    if st["definitely_not_psd"]:
        return {
            "query": "moments_force_psd",
            "status": CONCLUSIVE,
            "answer": False,
            "conclusion": ("the localizing matrix for x is certifiably not PSD, so no "
                           "representing measure lives on [0, inf): the matrix has a "
                           "negative eigenvalue"),
            "localizing_matrix_status": st,
        }
    return {
        "query": "moments_force_psd",
        "status": INSUFFICIENT,
        "answer": None,
        "reason": ("the degree-1 localizing matrix is not refuted, but PSD-ness of a "
                   "truncated localizing matrix is necessary and not sufficient for "
                   "support in [0, inf); establishing it would need a flat extension "
                   "(Curto-Fialkow) that these moments do not supply"),
        "localizing_matrix_status": st,
    }


def query_negative_count(seq: Sequence[Any]) -> Dict[str, Any]:
    """Lower bound on the number of negative eigenvalues."""
    out: Dict[str, Any] = {"query": "minimum_negative_eigenvalue_count"}
    if len(seq) < 4:
        out.update(status=INSUFFICIENT, reason="need m_0..m_3")
        return out
    st = psd_status(localizing_matrix(seq, 1))
    if st["definitely_not_psd"]:
        out.update(status=CONCLUSIVE, minimum=1,
                   conclusion="at least one eigenvalue is negative",
                   localizing_matrix_status=st)
    else:
        out.update(status=INSUFFICIENT, minimum=0,
                   reason="no localizing obstruction found; 0 negatives stays feasible",
                   localizing_matrix_status=st)
    return out


def query_positive_count(seq: Sequence[Any]) -> Dict[str, Any]:
    """Lower bound on the number of positive eigenvalues."""
    out: Dict[str, Any] = {"query": "minimum_positive_eigenvalue_count"}
    if len(seq) < 4:
        out.update(status=INSUFFICIENT, reason="need m_0..m_3")
        return out
    # Support in (-inf, 0] would need the localizing matrix for -x to be PSD.
    st = psd_status(localizing_matrix(seq, 1, flip=True))
    if st["definitely_not_psd"]:
        out.update(status=CONCLUSIVE, minimum=1,
                   conclusion="at least one eigenvalue is positive",
                   localizing_matrix_status=st)
    else:
        out.update(status=INSUFFICIENT, minimum=0,
                   reason="no localizing obstruction found; 0 positives stays feasible",
                   localizing_matrix_status=st)
    return out


def query_smallest_eigenvalue(n: int, seq: Sequence[Any]) -> Dict[str, Any]:
    """Two-sided bounds on ``lambda_min`` from ``m1`` and ``m2``."""
    if len(seq) < 3:
        return {"query": "smallest_eigenvalue_bounds", "status": INSUFFICIENT,
                "reason": "need m_1 and m_2"}
    b = eigenvalue_bounds_from_two_moments(n, seq[1], seq[2])
    return {
        "query": "smallest_eigenvalue_bounds",
        "status": CONCLUSIVE,
        "lambda_min": b["lambda_min"],
        "lambda_max": b["lambda_max"],
        "tight": b["tight"],
        "method": b["method"],
        "note": b["note"],
    }


def query_inertia_determined(n: int, seq: Sequence[Any],
                             observed: Optional[Sequence[int]] = None) -> Dict[str, Any]:
    """Do the moments determine the inertia, or merely permit the observed one?

    At ``n = 2`` they determine it: two moments, two eigenvalues, and the map is
    invertible. Above that they generally do not, and the answer is
    ``INSUFFICIENT_INFORMATION`` rather than a shrug dressed as a result.
    """
    out: Dict[str, Any] = {"query": "inertia_determined_by_moments", "n": int(n)}
    if n == 2 and len(seq) >= 3:
        spec = spectrum_from_two_moments_n2(seq[1], seq[2])
        if spec["status"] != CONCLUSIVE:
            out.update(status=INSUFFICIENT, reason=spec.get("blocker"))
            return out
        lo_lo, lo_hi = (float(x) for x in spec["lambda_1"])
        hi_lo, hi_hi = (float(x) for x in spec["lambda_2"])
        sig: Optional[List[int]] = None
        if lo_hi < 0 < hi_lo:
            sig = [1, 1, 0]
        elif lo_lo > 0:
            sig = [2, 0, 0]
        elif hi_hi < 0:
            sig = [0, 2, 0]
        out.update(
            status=CONCLUSIVE if sig else INSUFFICIENT,
            determined=True,
            spectrum=spec,
            implied_inertia=sig,
            reason=(None if sig else
                    "the eigenvalue enclosures straddle zero; the spectrum is "
                    "determined but its signature is not resolved at this precision"),
            conclusion=("m1 and m2 invert to the spectrum of a 2x2 matrix, so the "
                        "inertia is a consequence of the moments, not an extra "
                        "observation"),
        )
        if observed is not None and sig is not None:
            out["matches_observed"] = list(observed) == sig
        return out
    out.update(
        status=INSUFFICIENT,
        determined=False,
        reason=(f"for n = {n} the map from a spectrum to (m1..m4) is not injective; "
                "the moments constrain the inertia but do not fix it"),
    )
    return out


def b1_hankel_view(seq: Sequence[Any]) -> Dict[str, Any]:
    """Hankel rank and flatness through B1, when the moments are exact.

    Reported as structural context: ``rank M_t`` lower-bounds the number of
    distinct eigenvalues, and flatness is what would license a unique
    representing measure at all.
    """
    exact = all(isinstance(x, (int, Fraction)) for x in seq)
    if not exact:
        return {"available": False,
                "reason": "B1 is exact-rational; these moments are enclosures"}
    if not b1_available():
        return {"available": False, "reason": "b1_moment_solver not importable"}
    from b1_moment_solver.exact import hankel_rank, is_flat

    view: Dict[str, Any] = {"available": True,
                            "engine": "b1_moment_solver.exact"}
    ms = [Fraction(x) for x in seq]
    for t in (1, 2):
        if len(ms) >= 2 * t + 1:
            view[f"rank_M{t}"] = hankel_rank(ms, t)
            view[f"flat_at_{t}"] = bool(is_flat(ms, t))
    if "rank_M2" in view:
        view["distinct_eigenvalues_at_least"] = view["rank_M2"]
    return view


def analyse(G: Sequence[Sequence[Any]], *,
            observed_inertia: Optional[Sequence[int]] = None,
            upto: int = 4) -> Dict[str, Any]:
    """Full moment analysis of a Hermitian matrix, ready to certify."""
    n = len(G)
    ms = spectral_moments(G, upto=upto)
    seq = _moment_sequence(n, ms)
    report = moment_report(G, upto=upto)
    queries = [
        query_forces_psd(seq),
        query_negative_count(seq),
        query_positive_count(seq),
        query_smallest_eigenvalue(n, seq),
        query_inertia_determined(n, seq, observed_inertia),
    ]
    return {
        "content_kind": KIND_SPECTRAL_MOMENT,
        "dimension": n,
        "moments": report["moments"],
        "moment_method": report["method"],
        "sanity_violations": report["sanity_violations"],
        "b1_hankel_view": b1_hankel_view(seq),
        "b1_queries": queries,
        "claim_scope": "finite_dimensional_weil_compression",
        "rh_proof_claim": False,
        "note": ("INSUFFICIENT_INFORMATION is a certified outcome (§6): it records "
                 "that these moments do not decide the question, which is different "
                 "from the question being undecided."),
    }
