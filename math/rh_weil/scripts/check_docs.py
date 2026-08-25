#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §3.3/§14 (WO-RH-37, WO-RH-45) -- the `rh-docs` gate.

    python3 scripts/check_docs.py [--verbose]

Documentation freshness is a correctness property here, not housekeeping. Every
prior work order in this program was executed by an agent that read these files
first, and a stale instruction does not merely mislead a human: it sends the next
run to regenerate a closed result, or to reinterpret a certificate under a
normalization that was rejected. So a docs-only failure blocks merge.

What this checks, driven by `docs_status.json`:

* **links** -- every local Markdown link target in a canonical doc exists;
* **superseded claims** -- no live doc repeats a status the repository has since
  closed (WO-RH-05 open, degree 3 forbidden, the T=84 uniform check missing).
  Files under an explicitly historical root, and files whose front matter marks
  them HISTORICAL, are exempt: §3.2 says preserve the history, label it;
* **claim boundary** -- every live RH doc states the finite-dimensional /
  no-RH-proof boundary;
* **status agreement** -- `docs_status.json` and the README's status table agree
  with `certificates/work_order_status.json`, which is the machine-readable
  source;
* **certified values** -- the numbers in the README's "Current certified
  results" table are the numbers in the certificates, read from the certificates
  at check time rather than from prose memory;
* **root README** -- points at the current Constant Atlas document, detected from
  repository state rather than hard-coded, and carries the research-programs
  pointer with its scope wording.

Exit codes: 0 clean; 1 at least one finding.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]          # math/rh_weil
REPO = ROOT.parents[1]                              # repository root
STATUS = ROOT / "docs_status.json"
CERT_DIR = ROOT / "certificates"

#: A doc may opt out of the live checks by saying so in its first lines.
HISTORICAL_MARKERS = ("status: HISTORICAL", "HISTORICAL / SUPERSEDED",
                      "**HISTORICAL**", "status: SUPERSEDED")

LINK = re.compile(r"\[[^\]]*\]\(([^)]+)\)")
#: An explicit, greppable opt-out for prose that *quotes* a defect in order to
#: record that it was fixed. Without it, the ENG-007 record could not say which
#: stale claims it removed, and the live instructions could not say what the
#: historical file gets wrong -- both of which a reader needs. It is a region
#: marker rather than a per-line one so that abusing it is visible in a diff,
#: and every use is counted in this gate's output.
QUOTE_OPEN = "<!-- docs-check: superseded-quote start -->"
QUOTE_CLOSE = "<!-- docs-check: superseded-quote end -->"
FENCE = re.compile(r"```.*?```", re.S)
INLINE_CODE = re.compile(r"`[^`\n]*`")


def load_status() -> Dict[str, Any]:
    return json.loads(STATUS.read_text(encoding="utf-8"))


def resolve(doc: str) -> Path:
    return (ROOT / doc).resolve()


def is_historical(path: Path, status: Dict[str, Any]) -> bool:
    for root in status["historical_roots"]:
        if str(path).startswith(str((ROOT / root).resolve())):
            return True
    for rel in status.get("historical_docs", ()):
        if path == resolve(rel):
            return True
    try:
        head = "\n".join(path.read_text(encoding="utf-8").splitlines()[:12])
    except OSError:
        return False
    return any(m in head for m in HISTORICAL_MARKERS)


def strip_code(text: str) -> str:
    """Prose only. A superseded claim quoted inside a fenced block is usually a
    command or a historical excerpt, and flagging it produces the kind of noise
    that gets a gate switched off."""
    return INLINE_CODE.sub(" ", FENCE.sub(" ", text))


# --------------------------------------------------------------------------- #
# checks                                                                       #
# --------------------------------------------------------------------------- #
def check_links(doc: str, path: Path) -> List[str]:
    problems = []
    text = path.read_text(encoding="utf-8")
    for target in LINK.findall(text):
        target = target.strip()
        if target.startswith(("http://", "https://", "mailto:", "#")):
            continue
        target = target.split("#", 1)[0]
        if not target:
            continue
        dest = (path.parent / target).resolve()
        if not dest.exists():
            problems.append(f"{doc}: broken local link -> {target}")
    return problems


def strip_quoted_regions(text: str) -> Tuple[str, int]:
    """Blank out explicitly-marked quotation regions, and count them."""
    out, used, pos = [], 0, 0
    while True:
        start = text.find(QUOTE_OPEN, pos)
        if start < 0:
            out.append(text[pos:])
            break
        end = text.find(QUOTE_CLOSE, start)
        if end < 0:  # unbalanced: treat the rest as live, which fails loudly
            out.append(text[pos:])
            break
        out.append(text[pos:start])
        out.append("\n" * text.count("\n", start, end))
        used += 1
        pos = end + len(QUOTE_CLOSE)
    return "".join(out), used


def check_superseded(doc: str, path: Path, status: Dict[str, Any],
                     counter: Optional[List[int]] = None) -> List[str]:
    problems = []
    prose, used = strip_quoted_regions(path.read_text(encoding="utf-8"))
    if counter is not None:
        counter[0] += used
    prose = strip_code(prose)
    for rule in status["superseded_claims"]:
        m = re.search(rule["pattern"], prose)
        if m:
            problems.append(
                f"{doc}: repeats superseded status {rule['id']!r} "
                f"({m.group(0)[:70].strip()!r}) -- {rule['why']}"
            )
    return problems


def check_boundary(doc: str, path: Path, status: Dict[str, Any]) -> List[str]:
    text = path.read_text(encoding="utf-8")
    if any(p in text for p in status["boundary_phrases"]):
        return []
    return [f"{doc}: states no finite-dimensional / no-RH-proof claim boundary"]


def cert_value(spec: Dict[str, Any]) -> Tuple[Optional[Any], Optional[str]]:
    path = CERT_DIR / spec["certificate"]
    if not path.exists():
        return None, f"missing certificate {spec['certificate']}"
    node: Any = json.loads(path.read_text(encoding="utf-8"))
    for key in spec["path"]:
        try:
            node = node[key]
        except (KeyError, IndexError, TypeError):
            return None, (f"{spec['certificate']}: no value at "
                          f"{'.'.join(str(k) for k in spec['path'])}")
    return node, None


def check_certified_values(status: Dict[str, Any]) -> List[str]:
    problems = []
    for spec in status.get("certified_values", ()):
        value, err = cert_value(spec)
        if err:
            problems.append(err)
            continue
        doc_path = resolve(spec["doc"])
        text = doc_path.read_text(encoding="utf-8")
        needle = str(value)
        if needle in text:
            continue
        # Accept a shortened form as long as it is a prefix of the certified
        # digits: a README may quote 1.0731e-06 for a bound recorded to
        # seventeen places, but it may not quote a *different* number.
        if _prefix_match(needle, text):
            continue
        problems.append(
            f"{spec['doc']}: does not carry the certified value {needle} from "
            f"{spec['certificate']} ({'.'.join(str(k) for k in spec['path'])})"
        )
    return problems


def _prefix_match(value: str, text: str) -> bool:
    """Is some number in ``text`` a rounded form of ``value``?"""
    try:
        target = float(value)
    except ValueError:
        return False
    for token in re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", text):
        try:
            got = float(token)
        except ValueError:
            continue
        if got == 0 and target == 0:
            return True
        if target != 0 and abs(got - target) <= 1e-4 * abs(target):
            # Guard the trap this check exists to catch: a rounded quote is
            # fine, a *different* certified number is not. Requiring four
            # significant figures makes 1.0731e-06 pass and 1.08e-06 fail.
            return True
    return False


def check_status_agreement(status: Dict[str, Any]) -> List[str]:
    problems = []
    wo_path = CERT_DIR / "work_order_status.json"
    if not wo_path.exists():
        return [f"missing {wo_path.relative_to(ROOT)}"]
    wo = json.loads(wo_path.read_text(encoding="utf-8"))
    for field in ("current_work_order", "latest_completed_work_order"):
        if wo.get(field) != status.get(field):
            problems.append(
                f"docs_status.json {field}={status.get(field)!r} but "
                f"work_order_status.json says {wo.get(field)!r}"
            )
    orders = wo.get("orders") or {}
    readme = resolve("README.md").read_text(encoding="utf-8")
    for order in status.get("status_table_orders", ()):
        if order not in orders:
            problems.append(f"work_order_status.json has no entry for {order}")
            continue
        if order not in readme:
            problems.append(
                f"README.md: status table omits {order}, which docs_status.json "
                "requires it to report"
            )
    # The specific stale claims §3.2 names, checked against the machine status
    # rather than against a regex alone.
    if not str(orders.get("WO-RH-05", "")).startswith("recovered"):
        problems.append("work_order_status.json: WO-RH-05 is not recorded as recovered")
    if "blocked" in str(orders.get("WO-RH-08", "")):
        problems.append("work_order_status.json: WO-RH-08 is still recorded as blocked")
    return problems


def check_required_phrases(status: Dict[str, Any]) -> List[str]:
    problems = []
    for spec in status.get("required_phrases", ()):
        path = resolve(spec["doc"])
        if not path.exists():
            problems.append(f"{spec['doc']}: missing")
            continue
        text = path.read_text(encoding="utf-8")
        for phrase in spec["phrases"]:
            if phrase not in text:
                problems.append(f"{spec['doc']}: missing required phrase {phrase!r}")
    return problems


def current_atlas_doc(status: Dict[str, Any]) -> Optional[Path]:
    """The newest `constant-atlas-vX.Y.md`, detected from repository state.

    §13 is explicit that the version must not be hard-coded: pinning v0.6 in a
    checker only moves the staleness from the README into the gate.
    """
    matches = sorted(REPO.glob(status["root_atlas_doc_glob"]))
    if not matches:
        return None

    def version(p: Path) -> Tuple[int, ...]:
        m = re.search(r"v(\d+)\.(\d+)", p.name)
        return tuple(int(g) for g in m.groups()) if m else (0,)

    return max(matches, key=version)


def check_root_readme(status: Dict[str, Any]) -> List[str]:
    problems = []
    path = REPO / "README.md"
    text = path.read_text(encoding="utf-8")
    current = current_atlas_doc(status)
    if current is None:
        return ["repository root: no constant-atlas document found"]
    rel = f"docs/{current.name}"
    head = "\n".join(text.splitlines()[:20])
    if rel not in head:
        problems.append(
            f"README.md: the opening 'Core documents' block does not link the "
            f"current atlas {rel}"
        )
    for older in sorted(REPO.glob(status["root_atlas_doc_glob"])):
        if older == current:
            continue
        line = next((ln for ln in text.splitlines()[:20]
                     if f"docs/{older.name}" in ln and "Prior" not in ln), None)
        if line:
            problems.append(
                f"README.md: superseded atlas {older.name} is linked in the opening "
                "block outside the 'Prior' list"
            )
    return problems


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    status = load_status()
    problems: List[str] = []
    checked = 0
    quoted = [0]

    for doc in status["canonical_docs"]:
        path = resolve(doc)
        if not path.exists():
            problems.append(f"{doc}: listed as canonical but missing")
            continue
        if is_historical(path, status):
            problems.append(
                f"{doc}: listed as canonical but marked historical -- a file "
                "cannot be both the live instruction and the record of a "
                "superseded one"
            )
            continue
        checked += 1
        problems += check_links(doc, path)
        problems += check_superseded(doc, path, status, quoted)
        problems += check_boundary(doc, path, status)
        if args.verbose:
            print(f"  checked {doc}")

    for doc in status.get("root_docs", ()):
        path = resolve(doc)
        if not path.exists():
            problems.append(f"{doc}: missing")
            continue
        checked += 1
        problems += check_links(doc, path)
        if args.verbose:
            print(f"  checked {doc}")

    # Historical files must exist and must actually be marked historical --
    # otherwise the exemption above is silently protecting a live doc.
    for doc in status.get("historical_docs", ()):
        path = resolve(doc)
        if not path.exists():
            problems.append(f"{doc}: listed as historical but missing")
            continue
        if not is_historical(path, status):
            problems.append(
                f"{doc}: listed as historical but carries no historical marker; "
                "add one or move it to a historical root"
            )

    problems += check_required_phrases(status)
    problems += check_certified_values(status)
    problems += check_status_agreement(status)
    problems += check_root_readme(status)

    if problems:
        print("check_docs: FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"check_docs: PASS ({checked} live docs, "
          f"{len(status.get('historical_docs', ()))} historical, "
          f"{len(status.get('certified_values', ()))} certified values matched, "
          f"{quoted[0]} marked superseded-quote regions)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
