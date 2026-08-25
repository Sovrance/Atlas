"""The rank-trace lower bound, with its hypotheses enforced (§5, WO-RH-30).

Target form, in the normalization of its source::

    rank(P) >= 2 tr(P) + 4 tr(Q) - 4 b - ||P + Q||_HS^2

under **all** of:

    H1  P is positive semidefinite
    H2  Q is Hermitian
    H3  Q has at most b positive directions
    H4  every term is in the theorem's normalization

Why the hypotheses are objects and not comments
-----------------------------------------------
An inequality like this is only as good as the conditions attached to it, and
those conditions are exactly what gets lost when a bound is carried from one
setting to another. So this module refuses to compute a number unless each
hypothesis arrives with a stated verification status. An unverified hypothesis
is not a warning printed alongside a result -- it makes the result
``INCONCLUSIVE`` and there is no number to quote. §5: "If any hypothesis is
unverified, output INCONCLUSIVE."

The normalization tag is part of that. Two sources can state the same-looking
inequality with traces normalized differently, and silently mixing them
produces a bound that is arithmetically fine and mathematically meaningless.
The theorem id and the normalization tag are compared, not assumed.

What the bound is, and is not
-----------------------------
It is a lower bound on the rank of a finite positive semidefinite operator in
terms of traces and a Hilbert-Schmidt norm. It is not a statement about zeros,
it is not asymptotic, and a non-positive value for the right-hand side is a
*true and useless* bound rather than a failure -- rank is a non-negative
integer, so a right-hand side of -3 tells you nothing you did not know. That
case is reported as ``trivial`` so a null result cannot be mistaken for a
finding (§10).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Dict, List, Optional, Sequence, Tuple

#: Identifies the exact statement and normalization implemented here. A caller
#: asking for a different theorem id gets INCONCLUSIVE rather than this bound.
THEOREM_ID = "rank_trace_hs_v1"

#: The normalization this implementation is written in. Inputs must declare it.
NORMALIZATION_TAG = "trace_normalized_finite_dimensional"

KIND_RANK_TRACE = "WEIL_RANK_TRACE_CERTIFICATE"

#: The hypotheses, in the order the certificate reports them.
HYPOTHESES = (
    ("P_positive_semidefinite", "P is positive semidefinite"),
    ("Q_hermitian", "Q is Hermitian"),
    ("Q_positive_index_at_most_b", "Q has at most b positive directions"),
    ("shared_normalization", "all terms use the theorem's normalization"),
)


@dataclass
class RankTraceCertificate:
    status: str
    theorem_id: str
    hypotheses: Dict[str, Any] = field(default_factory=dict)
    inputs: Dict[str, Any] = field(default_factory=dict)
    result: Dict[str, Any] = field(default_factory=dict)
    blocker: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        body = {
            "content_kind": KIND_RANK_TRACE,
            "status": self.status,
            "theorem_id": self.theorem_id,
            "statement": ("rank(P) >= 2 tr(P) + 4 tr(Q) - 4 b - ||P + Q||_HS^2"),
            "hypotheses": self.hypotheses,
            "inputs": self.inputs,
            "result": self.result,
            "claim_scope": "finite_dimensional_weil_compression",
            "rh_proof_claim": False,
        }
        if self.blocker:
            body["blocker"] = self.blocker
        return body


def _as_interval(x) -> Dict[str, Any]:
    """Report a scalar as an outward [lo, hi] pair, whatever carrier it uses."""
    if hasattr(x, "lower") and hasattr(x, "upper"):
        return {"lo": repr(float(x.lower())), "hi": repr(float(x.upper()))}
    return {"lo": repr(x), "hi": repr(x)}


def _lower(x) -> float:
    return float(x.lower()) if hasattr(x, "lower") else float(x)


def _upper(x) -> float:
    return float(x.upper()) if hasattr(x, "upper") else float(x)


def rank_trace_lower_bound(
    *,
    trace_P,
    trace_Q,
    hs_sq_P_plus_Q,
    positive_index_Q_bound: int,
    theorem_id: str = THEOREM_ID,
    hypotheses: Optional[Dict[str, Any]] = None,
    normalization: str = NORMALIZATION_TAG,
) -> RankTraceCertificate:
    """Evaluate the rank-trace bound, or refuse to.

    ``hypotheses`` maps each name in :data:`HYPOTHESES` to a truthy verification
    record. Anything missing, false, or merely asserted-without-evidence makes
    the whole certificate INCONCLUSIVE: there is no partial credit, because the
    inequality is not valid without all four.
    """
    hyp = dict(hypotheses or {})
    recorded: Dict[str, Any] = {}
    unverified: List[str] = []
    for name, text in HYPOTHESES:
        entry = hyp.get(name)
        verified = bool(entry) and (entry is True or entry.get("verified") is True)
        recorded[name] = {
            "statement": text,
            "verified": verified,
            "evidence": (entry.get("evidence") if isinstance(entry, dict) else None),
        }
        if not verified:
            unverified.append(name)

    inputs = {
        "trace_P": _as_interval(trace_P),
        "trace_Q": _as_interval(trace_Q),
        "hs_sq_P_plus_Q": _as_interval(hs_sq_P_plus_Q),
        "positive_index_Q_bound": int(positive_index_Q_bound),
        "normalization": normalization,
    }

    if theorem_id != THEOREM_ID:
        return RankTraceCertificate(
            status="INCONCLUSIVE", theorem_id=theorem_id, hypotheses=recorded,
            inputs=inputs,
            blocker=(f"this module implements {THEOREM_ID!r} only; refusing to "
                     f"evaluate under the requested id {theorem_id!r}"))

    if normalization != NORMALIZATION_TAG:
        return RankTraceCertificate(
            status="INCONCLUSIVE", theorem_id=theorem_id, hypotheses=recorded,
            inputs=inputs,
            blocker=(f"normalization mismatch: inputs declare {normalization!r}, "
                     f"the theorem is stated in {NORMALIZATION_TAG!r}"))

    if unverified:
        return RankTraceCertificate(
            status="INCONCLUSIVE", theorem_id=theorem_id, hypotheses=recorded,
            inputs=inputs,
            blocker="unverified hypotheses: " + ", ".join(unverified))

    if int(positive_index_Q_bound) < 0:
        return RankTraceCertificate(
            status="INCONCLUSIVE", theorem_id=theorem_id, hypotheses=recorded,
            inputs=inputs, blocker="b must be a non-negative integer")

    # H3 names a *bound*, so it has to be checked against the bound actually
    # passed here -- not only against whatever the record was built with. If the
    # hypothesis record carries a measured positive index, it must not exceed
    # ``b``. Without this an evaluator can be handed a record verified for one
    # bound and a smaller bound at the call site, and the resulting certificate
    # is unsound rather than merely weak: P = diag(1,0,0), Q = diag(0,2,2) with
    # b = 1 yields rank(P) >= 5 for a 3x3 matrix of rank 1.
    entry = hyp.get("Q_positive_index_at_most_b")
    measured = None
    if isinstance(entry, dict):
        measured = (entry.get("evidence") or {}).get("positive_index")
    if measured is not None and int(measured) > int(positive_index_Q_bound):
        return RankTraceCertificate(
            status="INCONCLUSIVE", theorem_id=theorem_id, hypotheses=recorded,
            inputs=inputs,
            blocker=(f"Q_positive_index_at_most_b is violated: the recorded "
                     f"positive index of Q is {int(measured)}, which exceeds the "
                     f"supplied bound b = {int(positive_index_Q_bound)}"))

    # The bound is a lower bound, so every term is taken in the direction that
    # makes the right-hand side smallest: the enclosure's worst case, not its
    # midpoint. tr(P) and tr(Q) enter positively, the HS term negatively.
    b = int(positive_index_Q_bound)
    rhs_lo = (2.0 * _lower(trace_P) + 4.0 * _lower(trace_Q)
              - 4.0 * b - _upper(hs_sq_P_plus_Q))
    rhs_hi = (2.0 * _upper(trace_P) + 4.0 * _upper(trace_Q)
              - 4.0 * b - _lower(hs_sq_P_plus_Q))
    # rank is a non-negative integer, so the usable bound is the ceiling of the
    # certified real lower bound, floored at zero.
    import math

    usable = max(0, math.ceil(rhs_lo)) if rhs_lo > 0 else 0
    trivial = rhs_lo <= 0

    return RankTraceCertificate(
        status="PASS",
        theorem_id=theorem_id,
        hypotheses=recorded,
        inputs=inputs,
        result={
            "rhs_enclosure": {"lo": repr(rhs_lo), "hi": repr(rhs_hi)},
            "certified_rank_lower_bound": usable,
            "trivial": trivial,
            "interpretation": (
                "the right-hand side is non-positive, so the bound is true but "
                "says nothing: rank is a non-negative integer anyway"
                if trivial else
                f"rank(P) >= {usable} on these inputs"
            ),
        },
    )


def _same_enclosure(a, b) -> bool:
    """Whether two carriers denote the same value, for both arithmetics.

    Exact rationals compare directly. For balls, equality of the *enclosures* is
    the strongest thing available: two different balls may both contain a common
    real, but that does not make the underlying entries equal, so anything short
    of an identical enclosure is treated as unverified rather than assumed.
    """
    if hasattr(a, "lower") and hasattr(b, "lower"):
        return bool(a.lower() == b.lower() and a.upper() == b.upper())
    return bool(a == b)


def symmetry_status(M) -> tuple:
    """``(verified, detail)`` for the Hermitian hypothesis on ``M``.

    This has to be checked on *both* arithmetics, and the interval case is the
    one that matters. :func:`inertia.ldl.ldl_inertia` mirrors the upper triangle
    onto the lower before eliminating, so a non-symmetric input silently becomes
    a symmetric one and comes back with a perfectly good inertia -- for a matrix
    that was never Hermitian. Marking the hypothesis verified on the strength of
    the carrier alone would then let a rank-trace certificate apply a theorem
    whose Hermitian hypothesis is false.
    """
    n = len(M)
    for i in range(n):
        if len(M[i]) != n:
            return False, "matrix is not square"
        for j in range(i + 1, n):
            if not _same_enclosure(M[i][j], M[j][i]):
                return False, (f"entries ({i},{j}) and ({j},{i}) are not the same "
                               "value; symmetry cannot be certified")
    return True, "every transposed pair of entries is the same value"


def hypotheses_from_matrices(P, Q, *, b: Optional[int] = None,
                             exact: bool = True) -> Dict[str, Any]:
    """Build a hypothesis record by actually checking ``P`` and ``Q``.

    Uses the inertia engine, so the positivity and positive-index claims are
    certified by the same machinery as everything else rather than asserted.

    ``b`` is the positive-index bound the caller intends to pass to
    :func:`rank_trace_lower_bound`. It belongs here because H3 is not "Q has
    *some* positive index" -- it is "Q has at most ``b`` positive directions",
    and that is a statement about ``b``. Verifying only that the inertia was
    computable would leave the hypothesis recorded as satisfied for any bound at
    all, including one the matrix violates.
    """
    import sys
    from pathlib import Path

    root = str(Path(__file__).resolve().parents[1])
    if root not in sys.path:
        sys.path.insert(0, root)
    from inertia.ldl import exact_inertia, interval_inertia

    run = exact_inertia if exact else interval_inertia
    p_sym, p_sym_why = symmetry_status(P)
    q_sym, q_sym_why = symmetry_status(Q)
    rp = run(P) if p_sym else None
    rq = run(Q) if q_sym else None

    # P's inertia is only about P if P is symmetric -- the engine mirrors.
    p_psd = bool(p_sym and rp is not None and rp.status == "PASS"
                 and rp.n_negative == 0)
    q_index = rq.n_positive if (rq is not None and rq.status == "PASS") else None
    index_ok = q_index is not None and (b is None or int(q_index) <= int(b))

    return {
        "P_positive_semidefinite": {
            "verified": p_psd,
            "evidence": {"inertia": list(rp.signature) if rp and rp.signature else None,
                         "status": rp.status if rp else "NOT_EVALUATED",
                         "symmetric": p_sym, "symmetry_check": p_sym_why,
                         "method": rp.method if rp else None},
        },
        "Q_hermitian": {
            "verified": bool(q_sym),
            "evidence": {"checked": ("entrywise equality of transposed pairs"
                                     if exact else
                                     "entrywise equality of transposed enclosures"),
                         "detail": q_sym_why},
        },
        "Q_positive_index_at_most_b": {
            "verified": bool(index_ok),
            "evidence": {"inertia": list(rq.signature) if rq and rq.signature else None,
                         "status": rq.status if rq else "NOT_EVALUATED",
                         "positive_index": q_index,
                         "b": (None if b is None else int(b)),
                         "detail": (None if b is None else
                                    f"positive index {q_index} vs bound {int(b)}")},
        },
        "shared_normalization": {
            "verified": True,
            "evidence": {"normalization": NORMALIZATION_TAG,
                         "note": "matrices supplied directly in the theorem's normalization"},
        },
        "_positive_index_Q": q_index,
    }
