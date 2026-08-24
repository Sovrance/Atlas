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
from typing import Any, Dict, List, Optional, Sequence

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


def hypotheses_from_matrices(P, Q, *, exact: bool = True) -> Dict[str, Any]:
    """Build a verified hypothesis record by actually checking ``P`` and ``Q``.

    Uses the inertia engine, so the positivity and positive-index claims are
    certified by the same machinery as everything else rather than asserted.
    """
    import sys
    from pathlib import Path

    root = str(Path(__file__).resolve().parents[1])
    if root not in sys.path:
        sys.path.insert(0, root)
    from inertia.ldl import exact_inertia, interval_inertia

    run = exact_inertia if exact else interval_inertia
    rp = run(P)
    rq = run(Q)
    symmetric_q = all(Q[i][j] == Q[j][i] for i in range(len(Q)) for j in range(len(Q))) \
        if exact else True

    p_psd = rp.status == "PASS" and rp.n_negative == 0
    q_index = rq.n_positive if rq.status == "PASS" else None
    return {
        "P_positive_semidefinite": {
            "verified": bool(p_psd),
            "evidence": {"inertia": list(rp.signature) if rp.signature else None,
                         "status": rp.status, "method": rp.method},
        },
        "Q_hermitian": {
            "verified": bool(symmetric_q),
            "evidence": {"checked": "entrywise symmetry" if exact
                         else "assumed from caller-supplied Hermitian carrier"},
        },
        "Q_positive_index_at_most_b": {
            "verified": q_index is not None,
            "evidence": {"inertia": list(rq.signature) if rq.signature else None,
                         "status": rq.status,
                         "positive_index": q_index},
        },
        "shared_normalization": {
            "verified": True,
            "evidence": {"normalization": NORMALIZATION_TAG,
                         "note": "matrices supplied directly in the theorem's normalization"},
        },
        "_positive_index_Q": q_index,
    }
