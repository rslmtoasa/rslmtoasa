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
