"""An independent assembly of the 5x5 even block (ATLAS-RH-ENG-011 §WO-RH-77).

Same construction and same independence argument as :mod:`independent_even4`:
the machinery is :mod:`independent_even3`'s (SymPy kernels from the defining
integral, mpmath quadrature), which imports no production assembly module.

Evidence class is **E3**; nothing here is promoted. No RH proof claim.
"""
from __future__ import annotations

import math
from typing import Any, List, Tuple

import independent_even3 as IE3

CLAIM_SCOPE = "finite_dimensional_weil_compression"
EVIDENCE_CLASS = "E3"

BASIS: Tuple[str, ...] = ("one", "b", "b2", "bcube", "bquart")

CELL: Tuple[float, float] = (math.log(3.0), math.log(4.0))


def gram_matrix(L_val, mp, *, dps: int = 45) -> List[List[Any]]:
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
