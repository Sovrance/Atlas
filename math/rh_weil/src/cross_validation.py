"""Interval cross-validation across Weil providers (WO-RH-18).

Agreement is decided by **interval overlap**, not by decimal closeness: two
measurements agree only when their enclosures ``[v±r]`` intersect. The report
carries both the absolute overlap (intersection length, or a negative number
equal to the gap when disjoint) and a scale-free relative figure, so a
disagreement can never hide behind a small absolute number.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Any, Dict, List, Optional, Sequence

from providers import Measurement

DISAGREE = "DISAGREE"
AGREE = "AGREE"
UNAVAILABLE = "UNAVAILABLE"


@dataclass(frozen=True)
class PairComparison:
    left: str
    right: str
    status: str
    overlap_abs: Optional[float] = None
    overlap_rel: Optional[float] = None
    separation: Optional[float] = None  # |Δ| / (r_left + r_right); ≤1 ⇔ overlap
    delta: Optional[float] = None
    ratio: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "left": self.left,
            "right": self.right,
            "status": self.status,
            "overlap_abs": self.overlap_abs,
            "overlap_rel": self.overlap_rel,
            "separation": self.separation,
            "delta": self.delta,
            "ratio": self.ratio,
        }


def compare(name_a: str, a: Optional[Measurement], name_b: str, b: Optional[Measurement]) -> PairComparison:
    if a is None or b is None:
        return PairComparison(name_a, name_b, UNAVAILABLE)
    lo = max(a.lo(), b.lo())
    hi = min(a.hi(), b.hi())
    overlap_abs = hi - lo  # negative when disjoint (= -gap)
    widths = max(a.hi() - a.lo(), b.hi() - b.lo(), 0.0)
    scale = max(abs(a.value), abs(b.value), 1.0)
    overlap_rel = (overlap_abs / widths) if widths > 0 else (0.0 if overlap_abs < 0 else 1.0)
    denom = a.rad + b.rad
    delta = a.value - b.value
    separation = (abs(delta) / denom) if denom > 0 else (0.0 if delta == 0 else float("inf"))
    ratio = (a.value / b.value) if b.value != 0 else None
    status = AGREE if overlap_abs >= 0 else DISAGREE
    # Guard against two razor-thin intervals that miss only by float noise.
    if status == DISAGREE and abs(delta) <= 1e-12 * scale:
        status = AGREE
    return PairComparison(name_a, name_b, status, overlap_abs, overlap_rel, separation, delta, ratio)


def compare_all(measurements: Dict[str, Optional[Measurement]]) -> List[PairComparison]:
    """Pairwise comparison over every provider that supplied a value."""
    out: List[PairComparison] = []
    for a, b in combinations(sorted(measurements), 2):
        out.append(compare(a, measurements[a], b, measurements[b]))
    return out


def summarize(pairs: Sequence[PairComparison]) -> Dict[str, Any]:
    considered = [p for p in pairs if p.status != UNAVAILABLE]
    disagreements = [p for p in considered if p.status == DISAGREE]
    return {
        "pairs_total": len(pairs),
        "pairs_compared": len(considered),
        "pairs_unavailable": len(pairs) - len(considered),
        "disagreements": len(disagreements),
        "status": AGREE if considered and not disagreements else (
            UNAVAILABLE if not considered else DISAGREE
        ),
        "worst_separation": max((p.separation for p in considered), default=None),
    }
