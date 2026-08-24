"""Rank-trace / Hilbert-Schmidt certificates (ATLAS-RH-ENG-006 §5).

The theorem this implements is a *finite-dimensional* inequality with named,
checkable hypotheses. It is implemented in the exact normalization of its
source and nothing is generalized by analogy: §5 is explicit that constants may
not be adjusted to fit, and §14.5 requires the conclusion to name the theorem
and its hypotheses. Every hypothesis is verified or the result is INCONCLUSIVE.

No RH proof claim is made by this package.
"""
from .theorem import (  # noqa: F401
    KIND_RANK_TRACE,
    THEOREM_ID,
    RankTraceCertificate,
    hypotheses_from_matrices,
    rank_trace_lower_bound,
)
