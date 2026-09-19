# DRESP-06A spatial ALSDA comparison

Status: the bounded Fe execution gate passes, but the material campaign remains `CONVERGENCE_PENDING` / `FAIL-B`. The run is useful as an auditable implementation and oracle result; it is not a claim that the finite-q covariance and GF closures are converged.

## Scope and frozen state

The comparison uses one accepted 4 x 4 x 4 reciprocal k-space SCF state for bcc Fe, with Gamma, +q, -q, and a second finite q. All routes reuse the same eigenpairs, occupations, Fermi level, radial state, product basis, and response eta. `source/exchange.f90` and the `lr_lmto_turek_*` family remain frozen. No Goldstone/BES/GCR correction, kernel rescaling, or empirical eigenvalue shift is applied.

The ALSDA route is solved in the complete weighted-orthonormal LMTO product basis and only then observed in the DRESP site basis:

```text
chi_prod^ALSDA = (I - chi_prod^KS Kxc)^(-1) chi_prod^KS
chi_site       = F^H chi_prod F
```

For Fe/spd with `response_lmax=4`, the retained product dimension is 232. `Kxc` follows the existing KXC-01 contract, `Bxc=(vxc_up-vxc_down)/2`, divided by the Pauli magnetization, with the documented mixed `P <- SR` action.

## Implemented checks

- `lr_dresp06a_bridge_mod` provides the exact DRESP `F^H chi F` bridge, loss projection, and an independent reconstruct-pointwise-multiply-project Kxc action oracle.
- DRESP-05 Juelich static selection uses the minimum eta and the second-smallest distinct eta as holdout. The static chi inputs remain complex.
- Mills, Juelich, and ALSDA frequency refinements use route-specific coarse peak windows.
- The one-site downfolded diagnostic records `K_eff = chi0_site^-1 - chi_ALSDA_site^-1` only when both scalar inversions are directly conditioned; no regularization is used.
- Product-GF and Lehmann are compared at a representative finite-q/nonzero-frequency spot when `gf_closure_audit=.true.`.
- q/−q covariance uses the established compact circular angular/radial transport. The direct independently-compressed site conjugation mismatch is retained separately as a diagnostic.

The focused oracle is `UnitDresp06aProjection`; the independent finite-H known-U and eta/complex-holdout checks are in `UnitLrProjectedJuelichInteraction`. Existing compact Dyson and GF oracles remain separate tests.

## Bounded Fe ledger

Artifact: `/tmp/dresp06a_fe_4k.dat` (deliberately not committed).

| Diagnostic | Measured value | Interpretation |
|---|---:|---|
| Product dimension | 232 | Complete one-site Fe/spd product space |
| Kxc point identity, absolute / relative | 1.39e-17 / 5.72e-18 | Pass |
| Independent Kxc action oracle | 1.33e-15 | Pass |
| Compact-to-site bare closure | 3.07e-12 | Pass |
| Loss projection commutation | 1.42e-13 | Pass |
| Magnetization projection residual / relative | 8.02e-4 / 1.17e-3 | Reported finite product-space projection error |
| Raw ALSDA Gamma Ward residual / relative | 1.61e-1 / 2.35e-1 | Not a Goldstone-corrected pass |
| Mills U | -2.55976e-2 Ry | Scalar projected comparison |
| Juelich real U | -3.14865e-2 Ry | Eta-limited local-sum-rule comparison |
| Juelich selected eta | 0.01 Ry | Minimum ladder eta |
| Juelich holdout relative residual | 1.80e-1 | Eta stability is false |
| GF spot relative difference | 7.98e-1 | 101-point bounded spot; not quadrature-converged |
| Compact interacting covariance max residual | 8.68e-1 | `FAIL-B`; requires product/channel investigation |
| Projected covariance residual | 1.02e1 | Same failure survives the canonical site projection |
| Direct site conjugation mismatch | 1.02e1 | Diagnostic only; independent SVD/circular-coordinate comparison is not the canonical transport |

The route refinements were emitted independently for `spd_Mills`, `spd_Juelich`, and `ALSDA`. The scalar ALSDA comparison is serialized alongside the bare, Mills, and Juelich site rows for every requested eta/q/frequency.

## Verdict and next bounded work

`PASS` for representation algebra, independent Kxc action, exact site projection, loss commutation, finite-H Juelich recovery, eta selection, and end-to-end execution.

`FAIL-B / CONVERGENCE_PENDING` for material acceptance: the raw Gamma Ward residual is not small, Juelich eta stability is false, the bounded product-GF spot is underresolved, and interacting q/−q covariance fails at the current accepted-state/product-channel seam. No correction was applied to make these diagnostics pass.

The next campaign should isolate the covariance defect at the compact bare-response and Kxc transport levels, then repeat the GF spot with a controlled integration ladder. Only after those are closed should a material comparison be promoted beyond this bounded DRESP-06A gate.

## DRESP-06A-FINAL closure

The historical first bounded-run ledger above is retained unchanged.  The final
closure reran the accepted Fe state on the existing 4 x 4 x 4 mesh with the
commensurate pair `q=(+0.25,0,0)` and `-q=(-0.25,0,0)`.  The old `q=0.03`
result is therefore classified as a **NONCOMMENSURATE DISCRETE-QUADRATURE
DIAGNOSTIC**, not as a covariance implementation failure.  No 8 x 8 x 8 mesh
repeat was run; the optional mesh ladder remains outside this bounded closure.

### Commensurate covariance

The hard covariance gate passed at every one of the nine dynamic frequencies.
The maximum residuals were:

| Stage | Maximum residual or difference |
|---|---:|
| Bare compact response | 5.44e-15 |
| Kxc transported representation | 2.60e-11 relative Frobenius |
| Dyson denominator matrix | 1.99e-13 |
| Denominator minimum-singular-value difference | 1.83e-15 |
| Denominator condition-number difference | 2.40e-14 |
| Interacting compact response | 2.37e-14 |
| Compact loss | 3.11e-15 |
| Projected site response | 4.32e-13 |
| Projected site loss | 6.04e-14 |

The output records `q/channel covariance status = PASS`, so the final
classification is **DRESP-06A COVARIANCE = CLOSED**.  The new focused
`UnitDresp06aCommensurateCovariance` fixture exercises bare response, Kxc
transport, Dyson/interacting response, site projection, and loss covariance on
a genuinely finite, grid-commensurate q pair.

### Controlled product-GF closure

The representative spot uses finite `q=+0.25`, nonzero frequency index 2, and
`integration_eta=0.001 Ry`.  The energy interval is
`[-1.3219963339, 2.4984082194] Ry`.

| Quadrature | N_E | h (Ry) | h/eta_int | norm(chi_Lehmann) | norm(chi_GF) | dF | relative dF | dInf |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Base | 9553 | 3.999586e-4 | 0.399959 | 4.156046 | 4.160141 | 3.508615e-2 | 8.433885e-3 | 6.266478e-3 |
| Fine | 19105 | 1.999793e-4 | 0.199979 | 4.156046 | 4.160210 | 3.494005e-2 | 8.398628e-3 | 6.206559e-3 |

The base-to-fine GF difference is `dF=2.738490e-4`, relative `6.582578e-5`,
and `dInf=8.050153e-5`.  The controlled quadrature is stable and remains at
the previously accepted order of 1e-2 Lehmann relative discrepancy; it is not
an O(1) quadrature failure.

### Raw Ward eta ladder

For the complete 232-dimensional product space, the static denominator is
`D=I-chi0(0,0;eta) Kxc`.  The raw, uncorrected results are:

| eta (Ry) | Ward residual | relative | min sv | max sv | condition | min_abs_eig |
|---:|---:|---:|---:|---:|---:|---:|
| 0.040 | 2.274534e-1 | 3.309710e-1 | 1.860285e-1 | 1.335681 | 7.179983 | 2.519447e-1 |
| 0.020 | 1.775909e-1 | 2.584154e-1 | 9.636790e-2 | 1.342879 | 13.934916 | 1.317607e-1 |
| 0.010 | 1.614350e-1 | 2.349067e-1 | 4.997206e-2 | 1.344961 | 26.914267 | 6.865039e-2 |
| 0.005 | 1.570135e-1 | 2.284729e-1 | 2.738843e-2 | 1.345519 | 49.127296 | 3.768718e-2 |

The retained magnetization projection diagnostic is `R_m=1.167380e-3`
(absolute residual `8.022597e-4`), far below the Ward relative residual at
every eta.  The one-site unregularized downfolding remains descriptive only:
`K_eff^ALSDA(Gamma,omega=0)` is `-0.0334690`, `-0.0326048`, `-0.0321024`,
and `-0.0319324 Ry` along the eta ladder, compared with
`U_Mills=-0.0255976 Ry` and `U_Juelich=-0.0313139 Ry`.  Its scalar inversions
remain conditioned, but it is not used to repair the Ward identity.

### Denominator eigenmodes and SVD

The tracked candidate mode is mode 1 at every eta.  Its eigenvalue evolves as
`0.053740-0.246147i`, `0.003003-0.131726i`,
`-0.011960-0.067600i`, and `-0.016079-0.034085i Ry` for eta
`0.04, 0.02, 0.01, 0.005`, respectively.  Its magnitude decreases toward
zero, while its right-vector overlap with the compact magnetization stays
`0.9609--0.9612`; its biorthogonal weight is `0.9331--0.9372`, and its
magnetization fraction is `0.8707--0.8783`.

At eta `0.005`, the next reported eigenvalue magnitudes are approximately
`0.5964, 0.5964, 0.7086, 0.7086, 0.7276`; their magnetization overlaps are
negligible except for the last mode (`0.398`).  The smallest-right-singular-
vector overlap with `m` is `0.95425`, `0.95969`, `0.96111`, and `0.96146`
along the same eta ladder.  The eigenmode reconstruction and `Dm` modal
reconstruction residuals are at 1e-15--1e-14, validating the non-Hermitian
biorthogonal decomposition.

The decisive `Dm` decomposition is:

| eta (Ry) | candidate fraction | next four | remainder |
|---:|---:|---:|---:|
| 0.040 | 5.025e-1 | 6.32e-31 | 4.975e-1 |
| 0.020 | 1.898e-1 | 2.11e-29 | 8.102e-1 |
| 0.010 | 6.685e-2 | 1.02e-29 | 9.331e-1 |
| 0.005 | 1.334e-2 | 1.05e-30 | 9.867e-1 |

Thus the magnetization-like eigenmode becomes increasingly rigid, but it does
not carry the Ward residual as eta decreases.  The residual is distributed
over the remaining spectrum rather than concentrated in the candidate mode.
The raw Ward-defect classification is therefore
**DISTRIBUTED_KERNEL_INCONSISTENCY**, not a one-mode Halle-like defect.

### Final disposition

The final campaign is **PASS-B — NUMERICAL CLOSURE WITH MATERIAL LIMITATION**:
all representation, commensurate covariance, projected-response, loss, GF,
and Ward-algebra gates pass, while the raw physical/kernel consistency defect
remains large and only weakly eta-convergent.  The output metadata explicitly
remains:

```text
goldstone_correction = false
BES = false
GCR = false
kernel_rescaling = false
eigenvalue_pinning = false
```

The Halle/Buczek-Ernst-Sandratskii logic would require
`B_xc=chi0^{-1}(0)m=K_xc m`, hence `D m=0`, and a correction would act in the
Goldstone eigenspace.  This run shows a magnetization-aligned near-zero mode,
but the measured `Dm` defect is distributed.  DRESP-07 should therefore begin
with an ALSDA ground-state/response consistency audit, especially the mixed
`P <- SR` kernel contract, before considering any one-mode BES/Halle
correction.  No correction is implemented in DRESP-06A-FINAL.
