"""The 5x5 even Weil block ``G[{1, b, b^2, b^3, b^4}]`` (ATLAS-RH-ENG-011 §WO-RH-76/80).

The first block certified as a *prediction test*: ENG-009's two exploratory
even-sector scaling models disagree by roughly 8x about this block's
generalized gap, and §WO-RH-71 adjudicates them against the certified result
before any refit. The block itself extends the even sector by one direction::

    e0 = 1                          span{1}
    e1 = b     = x(L - x)           adds u^2   (u = x - L/2)
    e2 = b2    = x^2 (L - x)^2      adds u^4
    e3 = bcube = x^3 (L - x)^3      adds u^6

``bcube`` was E0-prepared in ENG-009 §WO-RH-62: its kernels, endpoint
polynomials, derivatives and reference-metric column all *derive* from
``basis_algebra.BASIS_L_POLY`` and were swept by the generic exact tests
(including direct symbolic integration) when it was added. Nothing about the
assembly here is new machinery; what is new is that there are ten independent
entries instead of six, and that the reference metric and the shifted pencil
``G - lam M`` are first-class outputs, because the generalized gap -- not the
raw determinant -- is this block's headline observable (ENG-009 §WO-RH-59).

Preconditioning follows ENG-008 exactly: a diagonal of exact powers of two,
frozen for the whole cell, applied without rounding, licensed by the ENG-007
congruence theorems. The raw diagonal spans ~1e5, so the raw fourth minor
sits ~1e-19 below entries of order 1e-1; the frozen dyadic congruence brings
the preconditioned minors to O(1)-O(1e-3) without adding any width.

No RH proof claim is made by this module. Everything here concerns one finite
block on one interval under one normalization.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import archimedean_realspace as AR
import basis_algebra
import pole
import reference_metric as RM
import weil_entries as WE

CLAIM_SCOPE = "finite_dimensional_weil_compression"

_names = ("one", "b", "b2", "bcube", "bquart")

#: The frozen basis, in Gram order (§WO-RH-65).
EVEN5_BASIS: Tuple[str, ...] = ("one", "b", "b2", "bcube", "bquart")

#: A stable identifier for that basis, so a certificate names what it certified.
EVEN5_BASIS_ID = "weil_even5_one_b_b2_bcube_bquart_v1"

#: And for the reference metric it is measured against (§WO-RH-65).
REFERENCE_METRIC_ID = "l2_gram_on_support"

CELL: Tuple[float, float] = (math.log(3.0), math.log(4.0))
CELL_LABEL: Tuple[str, str] = ("log(3)", "log(4)")

#: The ten independent entries, in the order a symmetric 4x4 stores them.
ENTRY_KEYS: Tuple[Tuple[str, Tuple[str, str]], ...] = tuple(
    (f"G{a}{b}", (_names[a], _names[b]))
    for a in range(5) for b in range(a, 5)
)

#: Default working precision. The raw fifth minor is far below 1e-19 while the entries
#: are O(1e-1), so the determinant loses ~18 digits to cancellation before any
#: interval widening is counted; 192 bits keeps the point enclosures far below
#: that floor.
DEFAULT_PRECISION_BITS = 224

#: The preconditioner exponents, **frozen for the cell** (§WO-RH-80), in the
#: ENG-008 convention ``D = diag(2^{-e_k})``: ``round(log2(sqrt(G_kk)))`` at
#: the cell midpoint. The first four agree with the ENG-010 even4 choice; the
#: bquart axis lands at -13 (its diagonal is ~2.4e-8).
PRECONDITIONER_EXPONENTS: Tuple[int, int, int, int, int] = (-2, -6, -10, -10, -13)


def basis_identity() -> Dict[str, Any]:
    """What the basis *is*, for a certificate to record rather than imply."""
    return {
        "basis_id": EVEN5_BASIS_ID,
        "elements": list(EVEN5_BASIS),
        "definitions": {
            "one": "1",
            "b": "x(L - x)",
            "b2": "x^2 (L - x)^2 = b(x)^2",
            "bcube": "x^3 (L - x)^3 = b(x)^3",
            "bquart": "x^4 (L - x)^4 = b(x)^4",
        },
        "parity_about_midpoint": {n: pole.basis_parity(n) for n in EVEN5_BASIS},
        "spans_in_u": "span{1, u^2, u^4, u^6, u^8} with u = x - L/2",
        "kernel_degrees_in_a": {
            f"{i}_{j}": basis_algebra.kernel_degree_in_a(i, j)
            for _, (i, j) in ENTRY_KEYS
        },
        "reference_metric_id": REFERENCE_METRIC_ID,
    }


# --------------------------------------------------------------------------- #
# Preconditioner (ENG-008 machinery, one axis wider)                           #
# --------------------------------------------------------------------------- #
def _exponent_for(value: Any) -> int:
    mid = abs(float(value.mid()) if hasattr(value, "mid") else float(value))
    if mid == 0.0 or not math.isfinite(mid):
        return 0
    return int(round(math.log2(math.sqrt(mid))))


def preconditioner_exponents(diagonal: Sequence[Any]) -> List[int]:
    return [_exponent_for(v) for v in diagonal]


def apply_preconditioner(matrix: Sequence[Sequence[Any]],
                         exponents: Sequence[int]) -> List[List[Any]]:
    """``D^T A D`` with ``D = diag(2^{-e_k})``, exactly (no width added)."""
    n = len(matrix)
    scales = [2.0 ** (-e) for e in exponents]
    return [[matrix[i][j] * scales[i] * scales[j] for j in range(n)]
            for i in range(n)]


def preconditioner_record(exponents: Sequence[int]) -> Dict[str, Any]:
    return {
        "form": "diagonal_dyadic",
        "frozen_for_cell": True,
        "definition": "D = diag(2^{-e_k}), e_k = round(log2(sqrt(G_kk)))",
        "exponents": [int(e) for e in exponents],
        "diagonal": [repr(2.0 ** (-int(e))) for e in exponents],
        "inverse_diagonal": [repr(2.0 ** int(e)) for e in exponents],
        "exactly_representable": True,
        "invertible": True,
        "invertibility_argument": (
            "diagonal with entries 2^{-e_k}, e_k a finite integer, so every "
            "entry is a nonzero dyadic rational and det D = 2^{-sum e_k} != 0 "
            "exactly; no enclosure is involved"
        ),
        "width_added": "none: scaling an Arb ball by a power of two is exact",
        "licensed_by": [
            "AtlasRH.posIndexAtLeast_congruence_iff",
            "AtlasRH.negIndexAtLeast_congruence_iff",
            "AtlasRH.rank_congruence",
            "AtlasRH.generalized_pencil_congruence",
        ],
        "claim_effect": (
            "none -- congruence by an invertible matrix preserves the positive "
            "index, the negative index and the rank; applied to G and M "
            "simultaneously it also preserves every generalized eigenvalue of "
            "the pencil, which is the ENG-009 invariance theorem"
        ),
    }


# --------------------------------------------------------------------------- #
# Minors                                                                       #
# --------------------------------------------------------------------------- #
def leading_minors(matrix: Sequence[Sequence[Any]]) -> List[Any]:
    """``Delta1..Delta4`` by division-free cofactor expansion.

    The 3x3 sub-expansion is exactly :func:`even3.leading_minors`'s, so the
    fourth minor is one cofactor layer on top of arithmetic the ENG-008
    certificates already exercised.
    """
    out = []
    for k in range(1, len(matrix) + 1):
        out.append(_det([row[:k] for row in matrix[:k]]))
    return out


def _det(m: Sequence[Sequence[Any]]) -> Any:
    n = len(m)
    if n == 1:
        return m[0][0]
    total = None
    for col in range(n):
        minor = [row[:col] + row[col + 1:] for row in m[1:]]
        term = m[0][col] * _det(minor)
        if col % 2:
            term = -term
        total = term if total is None else total + term
    return total


def minor_scale_factors(exponents: Sequence[int]) -> List[float]:
    """``Delta~_k = Delta_k * prod_{i<=k} 2^{-2 e_i}`` -- exact powers of two."""
    out, running = [], 1.0
    for e in exponents:
        running *= (2.0 ** (-e)) ** 2
        out.append(running)
    return out


# --------------------------------------------------------------------------- #
# Assembly (§WO-RH-68)                                                         #
# --------------------------------------------------------------------------- #
def _as_matrix(entries: Dict[str, Any]) -> List[List[Any]]:
    def get(a: int, b: int) -> Any:
        key = f"G{min(a, b)}{max(a, b)}"
        return entries[key]
    return [[get(a, b) for b in range(5)] for a in range(5)]


def reference_matrix(L: Any) -> List[List[Any]]:
    """The exact L2 reference metric on the caller's carrier.

    Each entry is a single rational monomial in ``L`` (checked at import of
    :mod:`reference_metric`), so the interval evaluation is one power and one
    scalar multiple -- no dependency problem exists to manage.
    """
    return RM.metric_matrix_over(EVEN5_BASIS, L)


def shifted_matrix(raw: Sequence[Sequence[Any]], m_ref: Sequence[Sequence[Any]],
                   lam_num: int, lam_den: int) -> List[List[Any]]:
    """``G - lam M`` with ``lam = lam_num / lam_den`` exact on the carrier."""
    n = len(raw)
    return [[raw[i][j] - m_ref[i][j] * lam_num / lam_den for j in range(n)]
            for i in range(n)]


def assemble_even5_arb(
    L_interval: Any,
    *,
    precision_bits: int = DEFAULT_PRECISION_BITS,
    precondition: bool = True,
    exponents: Optional[Sequence[int]] = None,
    prime_powers: Optional[Sequence[Tuple[int, int]]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """The rigorous 4x4 even block over an ``L``-interval (§WO-RH-68).

    Same centred mean-value route as every certified entry since ENG-005; the
    additions over :func:`even3.assemble_even3_arb` are the exact reference
    matrix and a shifted-pencil builder, because the generalized gap is this
    block's primary observable.
    """
    from interval_backend import require_flint

    _, arb, acb, ctx = require_flint()
    previous = ctx.prec
    ctx.prec = int(precision_bits)
    try:
        box = L_interval if hasattr(L_interval, "mid") else arb(L_interval)
        mid = float(box.mid()) if hasattr(box, "mid") else float(box)
        if prime_powers is None:
            prime_powers = WE.prime_powers_below(mid)

        entries: Dict[str, Any] = {}
        records: Dict[str, Any] = {}
        for key, (i, j) in ENTRY_KEYS:
            val, rec = AR.gram_entry_centred(
                i, j, box, arb, acb, prime_powers=prime_powers, options=options)
            entries[key] = val
            records[key] = rec

        raw = _as_matrix(entries)
        m_ref = reference_matrix(box)
        if precondition:
            chosen = list(PRECONDITIONER_EXPONENTS if exponents is None else exponents)
        else:
            chosen = [0, 0, 0, 0, 0]
        conditioned = apply_preconditioner(raw, chosen) if precondition else raw

        return {
            "basis": basis_identity(),
            "L": repr(mid),
            "L_radius": repr(float(box.rad()) if hasattr(box, "rad") else 0.0),
            "precision_bits": int(precision_bits),
            "entries": entries,
            "matrix": raw,
            "reference_matrix": m_ref,
            "preconditioned": conditioned,
            "preconditioner": preconditioner_record(chosen),
            "minor_scale_factors": minor_scale_factors(chosen),
            "preconditioned_applied": bool(precondition),
            "minors_raw": leading_minors(raw),
            "minors_preconditioned": leading_minors(conditioned),
            "prime_powers": [(int(q), int(p)) for q, p in prime_powers],
            "quadrature": records,
        }
    finally:
        ctx.prec = previous


def matrix_over(lo: float, hi: float, *, precision_bits: int = DEFAULT_PRECISION_BITS,
                precondition: bool = True,
                exponents: Optional[Sequence[int]] = None,
                options: Optional[Dict[str, Any]] = None) -> List[List[Any]]:
    """The preconditioned block enclosing the family over ``[lo, hi]``."""
    from interval_backend import interval_box, require_flint

    _, arb, _, _ = require_flint()
    box = interval_box(lo, hi)
    built = assemble_even5_arb(box, precision_bits=precision_bits,
                               precondition=precondition, exponents=exponents,
                               options=options)
    return built["preconditioned"] if precondition else built["matrix"]


def shifted_matrix_over(lo: float, hi: float, lam_num: int, lam_den: int,
                        *, precision_bits: int = DEFAULT_PRECISION_BITS,
                        exponents: Optional[Sequence[int]] = None,
                        options: Optional[Dict[str, Any]] = None) -> List[List[Any]]:
    """``D (G - lam M) D`` over a box -- what the gap covers evaluate.

    ``G`` and ``M`` are shifted *before* preconditioning, and the same frozen
    ``D`` is applied to the pencil, so the certified minors describe one fixed
    congruence of ``G - lam M`` over the whole cell.
    """
    from interval_backend import interval_box, require_flint

    _, arb, _, _ = require_flint()
    box = interval_box(lo, hi)
    built = assemble_even5_arb(box, precision_bits=precision_bits,
                               precondition=False, options=options)
    shifted = shifted_matrix(built["matrix"], built["reference_matrix"],
                             lam_num, lam_den)
    chosen = list(PRECONDITIONER_EXPONENTS if exponents is None else exponents)
    return apply_preconditioner(shifted, chosen)
