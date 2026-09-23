"""B15 -> PIR fact emission (b13_cdl pir-bridge pattern, work-order §6).

B15-POS and B15-ZERO facts live in the ``domain`` namespace at L2/E0 on the
DOMAIN layer; the GID candidate forest is ``pir.candidates.lattice_fact``
re-homed to ``analyst`` (see b15_surf/gid/identify.py). Assumption taint
(``asm:...``) rides on the fact and propagates by the existing invalidation
traversal — no bespoke logic.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from pir.models import AnalyzerRef, Fact, Warning_
from pir.types import PassTag


def domain_fact(*, benchmark: str, content: Dict, verdict: str, evidence_level: str,
                soundness: str, witness: Optional[Dict], impossibility_certificate: Optional[Dict],
                assumptions: List[str], warnings: List[Dict[str, str]],
                measurement_interface: str, version: str = "0.1.0") -> Fact:
    analyzer = AnalyzerRef(id=f"b15_surf.{benchmark.split('-')[1].lower()}", version=version,
                           tag=PassTag.SOUND if soundness == "SOUND" else PassTag.HEURISTIC)
    fid = Fact.compute_id(content, analyzer, assumptions=tuple(assumptions))
    return Fact(
        fact_id=fid, pir_level="L2", evidence_level=evidence_level, layer="DOMAIN",
        namespace="domain", status="SUPPORTED", analyzer=analyzer, content=content,
        created_at="1970-01-01T00:00:00Z", assumptions=tuple(assumptions),
        source_spans=({"artifact_id": f"certificates/{benchmark.lower()}_certificate.json",
                       "span": "results"},),
        measurement_interface=(measurement_interface,),
        warnings=tuple(Warning_(w["location"], w["text"]) for w in warnings),
        verdict=verdict, witness=witness, impossibility_certificate=impossibility_certificate,
    )
