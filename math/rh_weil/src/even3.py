"""The 3x3 even Weil block ``G[{1, b, b^2}]`` (ATLAS-RH-ENG-008 §WO-RH-47/50).

The first block in this program where the determinant is not the whole story.
On a 2x2 the trace and determinant fix the spectrum, so inertia, the spectral
moments and the rank-trace bound all collapse onto information the determinant
already carried -- which is what the ENG-006 information report found and said.
Here they have room to differ.

The basis
---------
Writing ``u = x - L/2``, the even sector about the cell midpoint is
``span{1, u^2, u^4}`` and the three elements supply it in order::

    e0 = 1                        span{1}
    e1 = b  = x(L - x)            adds u^2   (b = L^2/4 - u^2)
    e2 = b2 = x^2 (L - x)^2       adds u^4

The odd sector through degree 3 is exactly ``span{u, u^3}``, and
``q1^3 = (L^2/4) q1 - b3`` lies in it, so extending *there* at degree 3 is
impossible. ENG-007 §15 established both facts; ENG-008 §WO-RH-47 freezes this
basis.

Assembly
--------
Every entry is the canonical ``G = G0 - Gp + Ginf`` under the adjudicated
Candidate-A pole, with the prime block from the exact overlap kernels and the
archimedean term from the ENG-005 real-space form -- the same route, the same
modules and the same panel rules as every other certified entry in this
program. Nothing about the assembly is new; what is new is that there are six
independent entries instead of three.

Preconditioning
---------------
The raw block is badly scaled: ENG-007's preview measured a condition number up
to 1.3e5 across the cell, almost all of it the diagonal spread. §WO-RH-47 allows
a diagonal preconditioner provided the congruence is certified and the
mathematical claim is unchanged.

The preconditioner here is diagonal with **exact powers of two**::

    D = diag(2^{-e_0}, 2^{-e_1}, 2^{-e_2}),   e_k = round(log2(sqrt(G_kk)))

Two properties follow, and both matter. It is exactly invertible -- a diagonal
matrix of nonzero dyadic rationals, with ``D^{-1} = diag(2^{e_k})`` exactly, no
enclosure involved. And ``G~ = D^T G D`` adds **no width at all**: scaling an
Arb ball by a power of two is exact, so the preconditioned entries are the
original entries with their exponents shifted. A general Jacobi scaling by
``1/sqrt(G_kk)`` would have to round, and would inflate every entry it touched.

By ``AtlasRH.posIndexAtLeast_congruence_iff`` and ``AtlasRH.rank_congruence``
-- proved in ENG-007 against a pinned Mathlib -- congruence by an invertible
matrix preserves the positive index, the negative index and the rank. So the
inertia read off ``G~`` is the inertia of ``G``. That is the whole justification
for preconditioning at all, and it is a theorem rather than an assumption.

No RH proof claim is made by this module. Everything here concerns one finite
block on one interval under one normalization.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import archimedean_realspace as AR
import basis_algebra
import pole
import weil_entries as WE

CLAIM_SCOPE = "finite_dimensional_weil_compression"

#: The frozen basis, in Gram order (§WO-RH-47).
EVEN3_BASIS: Tuple[str, ...] = ("one", "b", "b2")

#: A stable identifier for that basis, so a certificate names what it certified.
EVEN3_BASIS_ID = "weil_even3_one_b_b2_v1"

CELL: Tuple[float, float] = (math.log(3.0), math.log(4.0))
CELL_LABEL: Tuple[str, str] = ("log(3)", "log(4)")

#: The six independent entries, in the order a symmetric 3x3 stores them.
ENTRY_KEYS: Tuple[Tuple[str, Tuple[str, str]], ...] = (
    ("G00", ("one", "one")),
    ("G01", ("one", "b")),
    ("G02", ("one", "b2")),
    ("G11", ("b", "b")),
    ("G12", ("b", "b2")),
    ("G22", ("b2", "b2")),
)

#: Default working precision. The block's third minor is ~1e-11 while its
#: entries are O(1e-1), so the determinant loses roughly eleven digits to
#: cancellation before any interval widening is counted.
DEFAULT_PRECISION_BITS = 160

#: The preconditioner exponents, **frozen for the cell** (§WO-RH-47).
#:
#: Choosing them per box from that box's own diagonal would give a slightly
#: better-scaled matrix, and would also mean the certified numbers were about a
#: different matrix on each box -- the sign conclusion would survive, since
#: every admissible ``D`` is a positive diagonal scaling, but a uniform
#: numerical bound would not be a statement about any single object. Freezing
#: them makes ``G~ = D^T G D`` one matrix family over the cell, congruent to
#: ``G`` by one fixed invertible ``D``.
#:
#: The values are ``round(log2(sqrt(G_kk)))`` at the cell midpoint. Across the
#: cell the per-box choice ranges over ``(-1,-7,-9)``, ``(-2,-7,-9)``,
#: ``(-2,-6,-9)`` and ``(-2,-6,-10)``, so the frozen choice is never worse than
#: one binary order of magnitude per axis anywhere on the cell.
PRECONDITIONER_EXPONENTS: Tuple[int, int, int] = (-2, -6, -9)


def basis_identity() -> Dict[str, Any]:
    """What the basis *is*, for a certificate to record rather than imply."""
    return {
        "basis_id": EVEN3_BASIS_ID,
        "elements": list(EVEN3_BASIS),
        "definitions": {
            "one": "1",
            "b": "x(L - x)",
            "b2": "x^2 (L - x)^2 = b(x)^2",
        },
        "parity_about_midpoint": {n: pole.basis_parity(n) for n in EVEN3_BASIS},
        "spans_in_u": "span{1, u^2, u^4} with u = x - L/2",
        "kernel_degrees_in_a": {
            f"{i}_{j}": basis_algebra.kernel_degree_in_a(i, j)
            for _, (i, j) in ENTRY_KEYS
        },
    }


# --------------------------------------------------------------------------- #
# Preconditioner                                                               #
# --------------------------------------------------------------------------- #
def _exponent_for(value: Any) -> int:
    """``round(log2(sqrt(v)))`` from a ball's midpoint, as a plain int.

    Read off the midpoint on purpose: the preconditioner is a *choice*, not a
    claim. Any nonzero diagonal would leave the inertia unchanged, so nothing is
    certified about this number -- it only has to be a good choice, and it has
    to be exactly representable, which a power of two is.
    """
    mid = abs(float(value.mid()) if hasattr(value, "mid") else float(value))
    if mid == 0.0 or not math.isfinite(mid):
        return 0
    return int(round(math.log2(math.sqrt(mid))))


def preconditioner_exponents(diagonal: Sequence[Any]) -> List[int]:
    return [_exponent_for(v) for v in diagonal]


def apply_preconditioner(matrix: Sequence[Sequence[Any]],
                         exponents: Sequence[int]) -> List[List[Any]]:
    """``D^T A D`` with ``D = diag(2^{-e_k})``, exactly.

    Multiplying an Arb ball by a power of two shifts its exponent and changes
    nothing else -- neither the midpoint's significand nor the radius' -- so
    this is a congruence performed without rounding. The enclosure of
    ``D^T A D`` is exactly the image of the enclosure of ``A``.
    """
    n = len(matrix)
    scales = [2.0 ** (-e) for e in exponents]
    return [[matrix[i][j] * scales[i] * scales[j] for j in range(n)]
            for i in range(n)]


def preconditioner_record(exponents: Sequence[int]) -> Dict[str, Any]:
    """Everything a reader needs to reconstruct ``D`` and check it is invertible."""
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
        ],
        "claim_effect": (
            "none -- congruence by an invertible matrix preserves the positive "
            "index, the negative index and the rank, so the inertia read off "
            "D^T G D is the inertia of G"
        ),
    }


# --------------------------------------------------------------------------- #
# Assembly                                                                     #
# --------------------------------------------------------------------------- #
def _as_matrix(entries: Dict[str, Any]) -> List[List[Any]]:
    return [
        [entries["G00"], entries["G01"], entries["G02"]],
        [entries["G01"], entries["G11"], entries["G12"]],
        [entries["G02"], entries["G12"], entries["G22"]],
    ]


def leading_minors(matrix: Sequence[Sequence[Any]]) -> List[Any]:
    """``Delta1, Delta2, Delta3`` -- the leading principal minors.

    Written out rather than looped so the 2x2 minor is literally the ``E2`` the
    degree-2 certificates bound, and the 3x3 is a cofactor expansion with no
    pivoting: on an interval carrier a division would widen, and there is
    nothing here that needs one.
    """
    a, b, c = matrix[0]
    _, d, e = matrix[1]
    _, _, f = matrix[2]
    d1 = a
    d2 = a * d - b * b
    d3 = (a * (d * f - e * e) - b * (b * f - e * c) + c * (b * e - d * c))
    return [d1, d2, d3]


def minor_scale_factors(exponents: Sequence[int]) -> List[float]:
    """How much each leading minor is scaled by ``D^T . D``.

    ``Delta~_k = Delta_k * prod_{i<k} d_i^2`` with ``d_i = 2^{-e_i}``, so every
    factor is an exact positive power of two: the sign of a minor is unchanged,
    and a certified bound on ``Delta~_k`` converts to one on ``Delta_k`` by an
    exact division rather than a second cover.
    """
    out, running = [], 1.0
    for e in exponents:
        running *= (2.0 ** (-e)) ** 2
        out.append(running)
    return out


def assemble_even3_arb(
    L_interval: Any,
    *,
    precision_bits: int = DEFAULT_PRECISION_BITS,
    precondition: bool = True,
    exponents: Optional[Sequence[int]] = None,
    prime_powers: Optional[Sequence[Tuple[int, int]]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """The rigorous 3x3 even block over an ``L``-interval (§WO-RH-50).

    ``L_interval`` may be a point or an Arb ball; a ball uses the centred
    mean-value form, which is what makes an interval ``L`` tractable at all.
    Returns the entries, the raw and preconditioned matrices, the leading
    minors of whichever matrix the caller will read, and a record of how it was
    built.
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
        if precondition:
            chosen = list(PRECONDITIONER_EXPONENTS if exponents is None else exponents)
        else:
            chosen = [0, 0, 0]
        conditioned = apply_preconditioner(raw, chosen) if precondition else raw

        return {
            "basis": basis_identity(),
            "L": repr(mid),
            "L_radius": repr(float(box.rad()) if hasattr(box, "rad") else 0.0),
            "precision_bits": int(precision_bits),
            "entries": entries,
            "matrix": raw,
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
    """The preconditioned block enclosing the family over ``[lo, hi]``.

    The shape :func:`inertia.stratify.certify_inertia_family` expects.
    """
    from interval_backend import interval_box, require_flint

    _, arb, _, _ = require_flint()
    box = interval_box(lo, hi)
    built = assemble_even3_arb(box, precision_bits=precision_bits,
                               precondition=precondition, exponents=exponents,
                               options=options)
    return built["preconditioned"] if precondition else built["matrix"]
