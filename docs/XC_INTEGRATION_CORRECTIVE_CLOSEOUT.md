# XCF-05 — XC integration corrective closeout

## Scope

This patch closes two post-reconciliation integration defects without changing
the reconciled XC formulas, density conventions, LMTO potential parameters, or
material-specific behavior:

1. the legacy `VXC0SP` to `PBEGGA` origin handoff now passes raw radial
   derivatives; and
2. the historical atomic-SCF control residual is separate from the integrated
   residual retained for diagnostics.

## `VXC0SP` → `PBEGGA` origin contract

The radial derivative producer in `VXC0SP` returns

```text
RHOP  = dn/dr
RHOPP = d2n/dr2
```

`PBEGGA` consumes those raw derivatives.  For `R > 0` it constructs

```text
laplacian(n) = d2n/dr2 + 2*(dn/dr)/R
```

and at the regular origin it constructs `laplacian(n) = 3*d2n/dr2`.
Consequently, `VXC0SP` must not pass a pre-formed spherical Laplacian.

The former origin path multiplied `RHOPP` by 3 for two spin channels.  Since
`PBEGGA` applies the origin limit itself, this applied the spherical factor
twice.  The corrected handoff passes raw channel derivatives in the
historical down/up argument order:

```text
NSP=2: RHODD(down) = RHOPP(down), RHODD(up) = RHOPP(up)
NSP=1: RHODD(up) = RHODD(down) = 0.5*RHOPP(total)
```

The first derivative is explicitly zero at the origin for a regular spherical
density.  The factor of 3 is therefore applied exactly once, inside
`PBEGGA`.

## Production-path regression

`UnitVxc0spGgaOrigin` calls the actual `VXC0SP` entry point with stored radial
density data generated from

```text
n(r) = n0 + a*r**2 + b*r**4
```

for both `NSP=2` and `NSP=1`.  It reconstructs the same radial derivatives,
compares the production origin potentials against direct `PBEGGA` evaluation
with raw `n''`, and checks the explicit spin conversion.  The polarized case
also compares against the historical factor-three input and requires a
different result.

## Atomic-SCF residual contract

`ATOMSC` now maintains two named quantities:

```text
DRHO_CONTROL    = sum(WGT * abs(RHO_new - RHO_old))
DRHO_INTEGRATED = sum(WGT * DRDI * abs(RHO_new - RHO_old))
```

`DRHO_CONTROL` is the historical unweighted mesh-coordinate L1 residual.  It
drives every existing nonlinear control decision:

| decision | historical condition | action |
|---|---|---|
| radial solver tolerance | `ITER >= 2` and `DRHO_CONTROL > 2` | use `TL = 1e-3`; otherwise retain `TOLRSQ = 1e-8` |
| `BETA1` switch | `MOD(ITER,3) == 2` and `DRHO_CONTROL < 1` | use `BETA1 = 0.5`; otherwise retain `BETA = 0.3` |
| inner termination | `DRHO_CONTROL < TOL`, with `TOL = 1e-6`, or `ITER == NITER-1` | mark the final pass |

The diagnostic residual file writes both quantities under the explicit column
names `drho_control` and `drho_integrated`.  `DRHO_INTEGRATED` is the
Jacobian/radial-quadrature-weighted physical electron-number residual.  It is
an observable for calibration and later analysis; it is not used for solver
tolerance selection, `BETA1`, termination, or the ordinary diagnostic control
branch.  No production threshold was recalibrated.

`UnitAtomscResidualControls` uses a deterministic radial update with distinct
unweighted and integrated residuals, checks values around the historical
thresholds 2, 1, and `1e-6`, and verifies that changing the integrated value
does not change any control result.

## Verification

The following focused checks were run in the configured GNU Fortran/libXC
build:

| test | result |
|---|---|
| `UnitVxc0spGgaOrigin` | PASS |
| `UnitAtomscResidualControls` | PASS |
| `UnitRadialGga` | PASS |
| `UnitLegacyPbeGga` | PASS |
| `UnitLibxcGgaRadial` | PASS |
| `UnitPbeGgaComparison` | PASS |
| `UnitXcSelectorSemantics` | PASS |
| `UnitXcLdaLegacyKernel` | PASS |
| `UnitLibxcXcBaseline` | PASS |
| `UnitXcLdaReconciliation` | PASS |
| `Example_bulk_bccFe_sd_smoke` | PASS |

The complete configured `ctest -L unit` run also passed all 11 XC-focused
tests and the two new corrective tests.  Its only failure was the unrelated
existing `AccP1ReciprocalBatchedSource` source-contract assertion about
`tests/benchmarks/accp1_real_material.py`; no XCF-05 test failed.

No extensive fcc-Fe seed sweep was run.  No threshold retuning or
material-specific physics tuning was performed.

## Completion checklist

- [x] raw `n''` passed from `VXC0SP` to `PBEGGA` at the origin
- [x] factor 3 applied exactly once
- [x] `NSP=1` origin convention corrected
- [x] `NSP=2` origin convention corrected
- [x] production-contract GGA origin test added
- [x] legacy GGA tests pass
- [x] libXC GGA tests pass
- [x] LDA tests unchanged
- [x] every `DRHO` control use mapped
- [x] historical control residual restored
- [x] integrated residual retained separately
- [x] `BETA1` behavior follows the historical metric
- [x] radial solver tolerance follows the historical metric
- [x] convergence termination follows the historical metric
- [x] no threshold retuning performed
- [x] no material-specific physics tuning performed
- [x] documentation completed
