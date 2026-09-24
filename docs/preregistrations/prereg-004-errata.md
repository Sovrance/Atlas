# prereg-004 — errata (recorded, not silently repaired)

prereg-004 is frozen at `1ccea83`, and its content hash (`prereg_sha256 8d68d570fa13020a…`) is
bound into certificate `b15-pos2-6facf90a799a`. Errata are therefore filed here. The frozen file
and the certificate are left exactly as they were generated. No registered rule, value, point,
threshold or verdict is affected by anything below.

## E1 — Kreĭn–Nudel′man theorem number (citation transcription)
- **Frozen text:** "Kreĭn & Nudel′man … **Theorem III.2.3** … as stated in … Curto–Fialkow 1991,
  **Remark 4.4**". It adds that a text extraction of the scan read "11.2.3" and interprets this
  as III.2.3 by analogy with Remark 4.2's "III.2.4".
- **What the paper prints** (checked 2026-09-24 by rendering the scanned pages of Houston J. Math.
  17(4) as images and reading them, not via text extraction):
  - p. 624, **Remark 4.4** (even case, the criterion prereg-004 relies on): "Krein and Nudel'man
    proved in **[KrN, Theorem II.2.3]** that (i) is equivalent to (v) A(k) ≥ 0 and
    (a + b)B(k − 1) ≥ abA(k − 1) + C, where C := (γ_{i+j})_{i,j=1}^k".
  - p. 622, **Remark 4.2** (odd case): "[KrN, Theorem **III**.2.4]".
- **Correct citation:** Curto & Fialkow, Houston J. Math. 17 (1991) 603–635, Remark 4.4 (p. 624),
  citing **[KrN, Theorem II.2.3]**. The text extraction was right, and the "III" reading was an
  over-correction.
- **Open:** Remarks 4.2 and 4.4 cite different chapters of KrN (III vs II). Only KrN itself can
  settle whether that is correct or a misprint in Curto–Fialkow; it has not been checked against
  KrN. Until it is, cite the theorem as "KrN Theorem II.2.3 as cited in Curto–Fialkow 1991,
  Remark 4.4". The mathematical content (positivity of the two matrices is necessary and
  sufficient for the even truncated Hausdorff problem) is unaffected, because it is stated
  verbatim in Remark 4.4.
- **Where the frozen/certified wording appears** (left as generated): prereg-004 (§"Registered
  target" and P4-OI-1); certificate `b15-pos2-6facf90a799a` fields `ground_truth_route` and
  `results.T7_tower_consistency.role`; the string literals in `tests/test_b15_pos_t2.py` that
  generate those two fields (changing them would change the certified artifact).
- **Corrected in place** (unfrozen, not certified): the `certify_t2.py` and
  `tests/test_b15_pos_t2.py` module docstrings, and `docs/sprint-B15-POS2-report.md`.
