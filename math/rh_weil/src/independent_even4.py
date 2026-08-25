"""An independent assembly of the 4x4 even block (ATLAS-RH-ENG-010 §WO-RH-66).

The same construction as :mod:`independent_even3` -- SymPy derives the overlap
kernels straight from their defining integral, mpmath does the quadrature --
extended by the ``bcube`` column. It deliberately reuses that module's
machinery rather than the production one: :mod:`independent_even3` imports none
of the production assembly modules (a test asserts it), so neither does this.

Evidence class is **E3**. mpmath never certifies anything in this program, and
nothing here is promoted; the job is regression against the rigorous assembly.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from typing import Any, List, Tuple

import independent_even3 as IE3

CLAIM_SCOPE = "finite_dimensional_weil_compression"
EVIDENCE_CLASS = "E3"

BASIS: Tuple[str, ...] = ("one", "b", "b2", "bcube")

CELL: Tuple[float, float] = (math.log(3.0), math.log(4.0))


def gram_matrix(L_val, mp, *, dps: int = 40) -> List[List[Any]]:
    previous = mp.mp.dps
    mp.mp.dps = dps
    try:
        pp = IE3.prime_powers_below(float(L_val))
        n = len(BASIS)
        out = [[None] * n for _ in range(n)]
        for a in range(n):
            for b in range(a, n):
                v = IE3.gram_entry(BASIS[a], BASIS[b], L_val, mp, pp)
                out[a][b] = v
                out[b][a] = v
        return out
    finally:
        mp.mp.dps = previous


def leading_minors(M) -> List[Any]:
    """``Delta1..Delta4`` by division-free cofactor expansion."""
    out = []
    for k in range(1, len(M) + 1):
        out.append(_det([row[:k] for row in M[:k]]))
    return out


def _det(m) -> Any:
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
