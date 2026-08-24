# Normalization Adjudication v0.1 (WO-RH-17)

**Work order:** ATLAS-RH-ENG-003 / WO-RH-17 (P0)
**Status:** ADJUDICATED
**Certificate:** `certificates/normalization_adjudication.json`
**Cross-check:** `certificates/normalization_crosscheck.json` (WO-RH-18)
**Claim boundary:** finite-dimensional Weil-form research only. **No RH proof claim.**

---

## 0. Verdict

| Candidate | Formula | Disposition |
|---|---|---|
| **A** — explicit-formula convolution form | `G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i⁺E_j⁻ + E_i⁻E_j⁺` | **ADOPTED** |
| **B** — repository even block | `G0_even = (√3/2)(v₊v₊ᵀ + v₋v₋ᵀ)` | **REJECTED** |

Candidate B is not a differently-normalized Candidate A. On the even sector

$$B = \frac{\sqrt3}{2}\cosh\!\left(\frac{L}{2}\right)\cdot A,$$

and that factor equals **1 at exactly one point, `L = log 3`**. It is a
multiplicative calibration fitted at a single test point, which the P0 runbook
explicitly forbids.

---

## 1. Conventions frozen (the seven axes)

| Axis | Frozen choice |
|---|---|
| Fourier transform | `Fhat(ξ) = ∫_ℝ F(x) e^{-iξx} dx` |
| reflection / tilde | `tilde_h(x) = conj(h(-x))`; real basis ⇒ `tilde_h(x) = h(-x)` |
| convolution | `F_ij = h_i * tilde_h_j`, supported on `[-L, L]` |
| pole term | `Fhat_ij(i/2) + Fhat_ij(-i/2)` |
| prime term | `Gp_ij = Σ_{q=p^k, log q < L} (log p/√q) · K_ij(log q; L)`, subtracted |
| archimedean | `Ginf_ij(T) = (1/π)∫_0^T h₊(t) Re(conj(H_i)H_j) dt`, one-sided |
| basis | unnormalised `1`, `q1 = x − L/2`, `b = x(L−x)` on `[0,L]` |

Assembly: `G = G0 − Gp + Ginf`. The content-addressed id of this frozen set is
recorded as `active_normalization_id` and is bound into every emitted PIR fact.

---

## 2. Derivation of the pole term (Candidate A)

With `Fhat(ξ) = ∫F(x)e^{-iξx}dx`, the evaluation points `ξ = ±i/2` give real
exponential weights, `e^{-i(\pm i/2)x} = e^{\pm x/2}`:

$$\widehat F_{ij}(i/2) = \int F_{ij}(x)\,e^{x/2}\,dx .$$

Insert `F_ij = h_i * tilde_h_j`, i.e. `F_ij(x) = ∫ h_i(y) h_j(y-x)\,dy`, and
substitute `u = y − x` (so `x = y − u`). The exponential separates,
`e^{x/2} = e^{y/2}e^{-u/2}`, and the double integral factorises:

$$\widehat F_{ij}(i/2) = \Big(\int h_i(y)e^{y/2}dy\Big)\Big(\int h_j(u)e^{-u/2}du\Big) = E_i^{+}E_j^{-}.$$

Symmetrically `Fhat_ij(-i/2) = E_i⁻E_j⁺`. Hence

$$\boxed{G^0_{ij} = E_i^{+}E_j^{-} + E_i^{-}E_j^{+}},\qquad
E_i^{\pm} = \int_0^L h_i(x)e^{\pm x/2}dx .$$

Every sign and factor above is derived; nothing is fitted.

### 2.1 Independent confirmation (real-space route)

Because `F_ij(a) + F_ij(−a) = K_ij(a)` — the *same* correlation kernel the prime
block uses — the pole term also equals

$$G^0_{ij} = \int_0^L K_{ij}(a)\,2\cosh(a/2)\,da .$$

This route shares no code path with the closed form. The two agree to
`< 1e-12` relative on all nine `(i,j)` entries at all four cross-check points.

### 2.2 Parity lemma

For `h(L−x) = ±h(x)` the substitution `x → L−x` gives

$$E^{-} = \pm\,e^{-L/2}E^{+}.$$

The basis splits: `1` and `b` are **even** about `L/2`, `q1` is **odd**. Two
consequences follow immediately.

* The pole matrix is **parity block diagonal** — even–odd cross entries vanish
  identically (verified numerically to `1e-12`).
* On the even sector `v₋ = e^{-L/2}v₊`, so `Δ = E_1⁺E_b⁻ − E_1⁻E_b⁺ ≡ 0` and the
  even pole block is **rank 1 under both candidates**.

That last point explains how the conflict survived earlier review: the repository
carries a *rank-1 / det-0 regression* for the pole block, and **both** candidates
pass it. A rank or determinant test cannot discriminate them; only the scale can.

---

## 3. Audit of Candidate B

The repository ships (`src/finite_weil.py::g0_even_block`, `POLE_EVEN_SCALE = "sqrt(3)/2"`):

```
G0_even = (√3/2)(v₊v₊ᵀ + v₋v₋ᵀ),   v± = (I₀±, I_b±)
```

**The building blocks are correct.** `I₀±` and `I_b±` were verified to equal
`E_1^±` and `E_b^±` exactly (symbolically). The defect is in the *assembly*:
Candidate A pairs opposite signs (`E⁺E⁻ + E⁻E⁺`), Candidate B pairs like signs
(`E⁺E⁺ + E⁻E⁻`) and then rescales.

Using `v₋ = e^{-L/2}v₊`:

* `A = 2e^{-L/2}\,v₊v₊ᵀ`
* `B = (√3/2)(1 + e^{-L})\,v₊v₊ᵀ`

so

$$\frac{B}{A} = \frac{\sqrt3}{4}\big(e^{L/2}+e^{-L/2}\big) = \frac{\sqrt3}{2}\cosh(L/2).$$

Solving `(√3/2)cosh(L/2) = 1` gives the **unique** solution `L = log 3`
(`cosh(log3 / 2) = 2/√3`).

| `L` | `B/A` | discrepancy |
|---|---|---|
| `log 3` | 1.000000000000 | 0.0000 % |
| 1.1059498113 | 1.001841114685 | +0.1841 % |
| 1.20 | 1.026642994301 | +2.6643 % |
| `log 4` | 1.082531754731 | **+8.2532 %** |

### 3.1 Why this is not a basis normalization

A change of basis acts on the Gram matrix by a **constant congruence**
`G ↦ SᵀGS`. Such a transformation cannot introduce a factor that *varies with
L*. The observed quotient `(√3/2)cosh(L/2)` does vary with `L`; therefore no
invertible basis transformation relates B to A, and B cannot be rescued as "A in
another normalization". No search for such an `S` can succeed, so none was
fitted. Candidate B is rejected.

### 3.2 Provenance of the constant

No derivation of `√3/2` from the explicit formula exists in the repository or its
history; the source comment records it as a *"Run-16 calibration"*. Combined with
§3 the reading is unambiguous: the constant was chosen so that the even block
reproduced a preferred value at `L = log 3`. Both the fitting and the preference
are forbidden by the runbook.

### 3.3 The odd sector was already correct

For the odd element, `E_{q1}^{+} = 2e^{L/4}A(L)` with
`A(L) = L cosh(L/4) − 4 sinh(L/4)`, hence

$$G^0[q_1,q_1] = 2E^{+}E^{-} = -2e^{-L/2}(E^{+})^2 = -8A^2,$$

which is exactly the value the repository already shipped
(`finite_weil.pole_odd_A` / `g0_odd_pivot`). **The defect is confined to the even
block.** Candidate A reproduces the odd sector with no change at all — itself a
strong consistency check on the adopted convention.

---

## 4. Consequences for existing certificates

The following are marked `promotion_state = QUARANTINED_NORMALIZATION_ADJUDICATION`,
with `hard_constraints_certified = false`, by
`scripts/quarantine_normalization.py`:

```
e1_scalar_log3_log4.json
e1_degree1_log3_log4.json
e1_degree2_compact_log3_log4.json
e1_fourier_T84_points.json
e1_fourier_T84_uniform_degree2.json
```

Nothing is deleted and **no certificate is relabelled E3** — the historical claim
is itself evidence. Prior state is preserved under `quarantine.prior_state`, and
work orders WO-RH-05 and WO-RH-09…15 read `quarantined_pending_WO-RH-17`.

Note the numerical consequence: values computed **at `L = log 3` are unaffected**
(the calibration fixed point), while values at `L = log 4` carry an ~8.25 % error
in the pole block. This is consistent with the conflict recorded in the
engineering spec — a curve whose interior minimum near `L ≈ 1.10595` disappeared
once the even block was rescaled by an `L`-dependent factor. Such results must
still be regenerated, not reinterpreted, under WO-RH-19/20.

The PIR bridge now refuses to promote any quarantined or stale-normalization
certificate (`pir_bridge.promotion_refusal`), and every emitted fact carries
`asm:normalization_id:<active id>`.

### 4.1 The quarantine has to survive a re-run

Marking the five files once is not enough. The `scripts/certify_*.py`
entrypoints predate this adjudication: they rebuild those same bodies from the
**rejected** even pole block in `src/finite_weil.py` and set
`hard_constraints_certified` from their own gates. Re-running one therefore
overwrote the marker and restored a promotable-looking certificate — the same
failure mode already noted for the hardcoded `work_order_status.json` builder,
left open for the certificates themselves. Observed directly:

```
$ python3 scripts/certify_scalar_e1.py
$ jq '.promotion_state, .hard_constraints_certified' certificates/e1_scalar_log3_log4.json
null
true
```

The quarantine is therefore enforced at the single point of write rather than by
one stamping pass:

* `src/normalization.py` owns the registry (`QUARANTINED_CERTIFICATES`,
  `QUARANTINE_STATE`, `QUARANTINE_REASON`, `quarantine_block`) — one source of
  truth, previously duplicated in the script and again in the tests.
* `certificate_io.write_certificate` re-asserts the marker for any registered
  file, recording whatever claim the writer supplied as `quarantine.prior_state`
  so the contrary evidence is preserved rather than discarded.
* The three certify scripts that wrote raw JSON (`certify_scalar_e1.py`,
  `certify_degree1_e1.py`, `certify_degree2_compact_e1.py`) now write through
  that function, so no entrypoint bypasses the guard.
* Only `scripts/quarantine_normalization.py --release` may lift a marker, via an
  explicit `allow_quarantine_change=True`. Releasing is legitimate only after the
  WO-RH-19/20 regeneration.
* `finite_weil.POLE_EVEN_SCALE_STATUS == "REJECTED_WO_RH_17"` marks the module
  that still assembles Candidate B, so it cannot be read as current. Replacing
  that assembly remains WO-RH-19/20 and is **not** done here.

Regression tests cover all of it (`QuarantineTests`); they fail against the
pre-fix tree.

---

## 5. Reproduction

```bash
cd math/rh_weil
python3 scripts/derive_normalization.py           # derivation + adjudication certificate
python3 scripts/run_normalization_crosscheck.py   # four-way cross-check (add --no-arch to skip the slow route)
python3 scripts/quarantine_normalization.py       # idempotent quarantine
PYTHONPATH=src python3 tests/test_normalization_adjudication.py
```

## 6. Acceptance gate (§3.6)

| Requirement | Status |
|---|---|
| every sign and factor derived | ✅ §2 |
| real-space and Fourier formulations analytically equivalent | ✅ §2.1 (`∫K·2cosh` route) |
| selected formula reproduces low-degree symbolic identities | ✅ §3.3 (odd pivot `−8A²`) |
| any basis transformation explicit and invertible | ✅ §3.1 (none exists; argued, not fitted) |
| no fitted scale remains | ✅ adopted path is a pure `E^±` product |
| two independent implementations agree at test points | ✅ WO-RH-18, `< 1e-12` |
| `normalization_adjudication.json` has `"status": "ADJUDICATED"` | ✅ |

---

## 7. Limitations (stated, not hidden)

* The adjudication is **symbolic and high-precision numeric**. Emission of a new
  E1 certificate still requires the interval backend (`python-flint`) plus a
  certified tail bound; this document promotes nothing.
* The archimedean cross-check is `high_precision_numeric`, not interval
  certified. Only rows labelled `interval_certified` may ever support an E1
  claim.
* `ConnesCvSProjectedProvider` remains an **external diagnostic**: its
  projection/truncation error is not certified, so it returns no certifying value.
* A first draft of the closed-form `H_i(t;L)` used the endpoint formula
  everywhere and cancelled catastrophically as `t → 0` (terms `~1/t^{n+1}`). The
  four-way harness caught it — the two archimedean routes disagreed by up to
  `1e+65` — and it was fixed with a series branch for `|aL| ≤ 1`. Recorded here
  because a cross-check that has never caught anything has not been tested.
