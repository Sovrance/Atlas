"""The 3x3 odd parity pilot for ENG-008 (ATLAS-RH-ENG-007 §15, WO-RH-46).

**Preparation only. This module promotes nothing.** §15: "run only E0 identities and an E3
conditioning/topology preview; do not promote a new E1 degree result in ENG-007 unless it is
essentially free." It is not free, so nothing here is E1.

Why a third odd element
-----------------------
ENG-006's odd degree-3 block is the pair `{q1, b3}`, which is `2x2`. At `2x2` the
determinant and trace already encode nearly all of the spectral behaviour, so the inertia and
spectral-moment channels ENG-006 built cannot say much that `det > 0` does not. They only
start earning their keep at `3x3`, where a signature `(2,0,1)` and a moment sequence carry
information a single determinant cannot.

So the pilot extends the **odd** sector -- the one that already has two elements -- with the
next odd polynomial under the midpoint reflection `sigma(x) = L - x`:

    q1 = x - L/2                    odd, degree 1
    b3 = x(L-x)(x - L/2)            odd, degree 3
    q3 = (x - L/2)^3                odd, degree 3   <- added here

All three are odd, so the block stays inside one parity sector and the even/odd
block-diagonalisation is untouched. This matters: mixing parities would destroy the
factorization `det G = O1 * E2` that the existing certificates rely on.

Conditioning
------------
`q1` and `q3` are both powers of the same midpoint coordinate, so the Gram block is closer to
a Hilbert-like matrix than the `{q1, b3}` pair was, and its condition number grows fast. The
`conditioning_preview` below reports that explicitly (E3, diagnostic only) so ENG-008 chooses
a basis scaling deliberately rather than discovering an ill-conditioned block after committing
to it. A rescaling of `q3` by a power of `L` is the obvious lever; the preview reports the
effect so the choice can be made on evidence.

Kernels
-------
Symmetrized prime-overlap convention, matching `core`:

    K_ij(a, L) = integral_0^(L-a) [ f_i(x) f_j(x+a) + f_i(x+a) f_j(x) ] dx

Derived exactly with SymPy and pinned by `tests/test_degree3_pilot_exact.py`, which
re-derives them from the basis rather than trusting these closed forms.
"""
from __future__ import annotations

from typing import Any, Dict, List

#: The odd basis of the pilot block, in index order.
PILOT_BASIS = ("q1", "b3", "q3")

#: E0 only. Nothing in this module is a rigorous interval result.
EVIDENCE_CLASS = "E0"

CONTENT_KIND_PILOT = "WEIL_DEGREE3_PILOT_PREVIEW"


def kernel_q1q3(a, L):
    """``K_q1q3(a; L)`` -- exact, verified against SymPy in the test suite."""
    d = L - a
    return d * (L**4 - 4 * L**3 * a + 6 * L**2 * a**2 - 4 * L * a**3 - 4 * a**4) / 40


def kernel_b3q3(a, L):
    """``K_b3q3(a; L)`` -- exact, verified against SymPy in the test suite."""
    d = L - a
    return (d * d
            * (L**5 + 2 * L**4 * a - 4 * L**3 * a**2 - 10 * L**2 * a**3
               - 16 * L * a**4 - 8 * a**5) / 560)


def kernel_q3q3(a, L):
    """``K_q3q3(a; L)`` -- exact, verified against SymPy in the test suite."""
    d = L - a
    return (d
            * (5 * L**6 - 30 * L**5 * a + 54 * L**4 * a**2 - 16 * L**3 * a**3
               - 16 * L**2 * a**4 - 16 * L * a**5 - 16 * a**6) / 1120)


#: The three new prime-overlap kernels the pilot needs, beyond ENG-006's `{q1, b3}` pair.
PILOT_KERNELS = {
    ("q1", "q3"): kernel_q1q3,
    ("b3", "q3"): kernel_b3q3,
    ("q3", "q3"): kernel_q3q3,
}


def scaled_basis_note(scale_exponent: int = 2) -> Dict[str, Any]:
    """The basis-scaling strategy, stated rather than left implicit.

    `q3` carries three powers of the midpoint coordinate where `q1` carries one, so on a cell
    of width `L` their natural magnitudes differ by `L^2`. Rescaling `q3 -> q3 / L^2` puts the
    two on comparable footing. Scaling a basis element is a **congruence** by a diagonal
    matrix, so by `AtlasRH.congruence_posDef_iff` it changes neither positive definiteness nor
    inertia -- only the conditioning. That is exactly why it is a safe lever.
    """
    return {
        "strategy": "diagonal_rescaling",
        "element": "q3",
        "factor": f"L^-{scale_exponent}",
        "justification": (
            "diagonal congruence; preserves inertia and definiteness "
            "(AtlasRH.congruence_posDef_iff), affects only conditioning"
        ),
    }


def conditioning_preview(L: float, entries: Dict[str, float]) -> Dict[str, Any]:
    """E3 conditioning/topology preview for the pilot block.

    Floating point, diagnostic only, and labelled as such: §15 asks for a preview, and the
    evidence policy is explicit that a floating scan never certifies. The output exists so
    ENG-008 can pick a basis scaling on evidence instead of discovering the problem later.

    `entries` maps the six independent Gram entries of the symmetric 3x3 odd block, keyed
    ``"q1q1"``, ``"q1b3"``, ``"q1q3"``, ``"b3b3"``, ``"b3q3"``, ``"q3q3"``.
    """
    import math

    g = [
        [entries["q1q1"], entries["q1b3"], entries["q1q3"]],
        [entries["q1b3"], entries["b3b3"], entries["b3q3"]],
        [entries["q1q3"], entries["b3q3"], entries["q3q3"]],
    ]

    # Leading principal minors: the Sylvester criterion the runtime already uses, extended to
    # 3x3. Reported, not certified -- these are floats.
    m1 = g[0][0]
    m2 = g[0][0] * g[1][1] - g[0][1] * g[0][1]
    m3 = (g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[1][2])
          - g[0][1] * (g[0][1] * g[2][2] - g[1][2] * g[0][2])
          + g[0][2] * (g[0][1] * g[1][2] - g[1][1] * g[0][2]))

    frob = math.sqrt(sum(v * v for row in g for v in row))
    diag_ratio = (max(abs(g[i][i]) for i in range(3))
                  / max(1e-300, min(abs(g[i][i]) for i in range(3))))

    return {
        "content_kind": CONTENT_KIND_PILOT,
        "evidence_class": "E3",
        "rigorous": False,
        "certifies": False,
        "promotion_state": "NOT_PROMOTED_PREVIEW_ONLY",
        "rh_proof_claim": False,
        "claim_scope": "finite_dimensional_weil_compression",
        "basis": list(PILOT_BASIS),
        "parity_sector": "odd",
        "L": L,
        "leading_principal_minors": {"m1": m1, "m2": m2, "m3": m3},
        "sylvester_all_positive_float": bool(m1 > 0 and m2 > 0 and m3 > 0),
        "frobenius_norm": frob,
        "diagonal_magnitude_ratio": diag_ratio,
        "scaling_strategy": scaled_basis_note(),
        "note": (
            "E3 preview. Floating point never certifies (evidence policy); the positive "
            "minors below are a topology hint for ENG-008, not a positivity result. A "
            "rigorous claim requires the interval engine and an E1 certificate."
        ),
    }


def pilot_summary() -> Dict[str, Any]:
    """What ENG-008 inherits from this preparation."""
    return {
        "work_order": "WO-RH-46",
        "status": "PREPARED_NOT_CERTIFIED",
        "basis": list(PILOT_BASIS),
        "parity_sector": "odd",
        "new_kernels": [f"K_{i}{j}" for (i, j) in PILOT_KERNELS],
        "evidence_class": EVIDENCE_CLASS,
        "e1_promoted": False,
        "why_3x3": (
            "ENG-006's odd block is 2x2, where determinant and trace already encode most "
            "spectral behaviour; inertia and moments only add information at 3x3."
        ),
        "open_for_eng008": [
            "rigorous interval assembly of the pilot block",
            "E1 positivity or inertia certificate on the cell",
            "basis scaling decision informed by the conditioning preview",
        ],
    }
