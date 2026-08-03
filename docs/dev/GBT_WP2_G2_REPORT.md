# GBT WP2 / G2 report

Date: 2026-08-03.

## Post-report architecture clarification

The numerical WP2 implementation is now classified as **G2E: PASS**. The
production operator slice remains **G2O: OPEN/FAIL**, so final G2 is still
**FAIL**.

The agreed representation is a collinear rotating-frame potential with the
GBT endpoint link applied only to primitive directed `S`/`Sdot` blocks. In
that representation a one-sublattice common cone at q=0 correctly produces
the identical collinear Hamiltonian; laboratory-frame moments are reconstructed
for output. The original failure below was produced by the superseded double
representation (`ham0m_nc` rotation followed by reciprocal `h0/bz`
reconstruction).

An uncommitted experiment commenting the q/cone branches out of `ham0m_nc`
restores q=0 equality and preserves the user's MX/MY/MZ periodic-NC energies,
but the current auto-probe angles then reach neither `ham0m_nc` nor the
structure layer. The diagnostic finite-q auto run consequently gives exactly
zero for q=`0`, `0.025`, and `0.05`. That is not G2O acceptance: the S/Sdot
endpoint link and a non-no-op finite-q or relative-sublattice assertion are
still required. The next permitted work is the documented WP3+WP4 S/Sdot
gate repair; later work packages remain blocked.

## Prerequisite

G1 passed before WP2 started. `UnitGbtOracles` reports a maximum algebraic
relative error of `1.780e-15`, below the `1e-12` G1 threshold. See
`docs/dev/GBT_WP1_G1_REPORT.md`.

## Canonical contract

`reciprocal_occupations.f90` now owns the projection-free occupation path:

\[
N = \frac{\sum_{kn} w_k f(\epsilon_{nk},E_F,T)}{\sum_k w_k},
\qquad
E_{\rm band} =
\frac{\sum_{kn}w_k f(\epsilon_{nk},E_F,T)\epsilon_{nk}}
     {\sum_k w_k}.
\]

The Fermi-level root and both sums use the same stable Fermi-Dirac function and
the same effective `kT` floor. Raw multiplicities and normalized weights are
both supported because the global raw weight sum is accumulated explicitly and
both observables are divided by it.

For a distributed k mesh, each rank accumulates its owned eigenvalues and the
three raw scalars `(weight, electron count, energy)` are combined in one
`MPI_ALLREDUCE`. A replicated mesh is not reduced. Tetrahedron and Blöchl DOS
modes first establish the full uniform mesh and use this full-mesh
Fermi-Dirac path as the canonical energy oracle.

`frozen_magnon_probe_energy` now performs
`build -> diagonalize -> calculate_canonical_band_energy`; it no longer builds
a DOS or reads projected moments. The k-space SCF band energy uses the same
canonical result. Total-DOS integration and projected `m0*m1` are computed and
logged only as cross-checks.

Every Hamiltonian rebuild calls `invalidate_spectral_cache`. Hamiltonian,
overlap, eigenpair, DOS, projection, tetrahedron, and canonical-energy state is
discarded, while the reusable k mesh, weights, and configured DOS grid are
retained.

## Tests and budgets

The production evaluator is exercised by `UnitKspaceOccupations` in both
serial and two-rank MPI builds. Its reference calculation is independent of the
production Fermi root/accumulator.

| Check | Observed error | Budget | Result |
| --- | ---: | ---: | --- |
| Electron conservation | `1.155e-12` electron | `2e-10` | PASS |
| Fermi level vs independent root | `1.020e-13 Ry` | `2e-10 Ry` | PASS |
| Canonical EBAND vs independent sum | `5.434e-13 Ry` | `2e-10 Ry` | PASS |
| Raw-weight rescaling, electron count | `0` | `2e-10` | PASS |
| Raw-weight rescaling, EBAND | `8.327e-17 Ry` | `2e-10 Ry` | PASS |
| DOS window/grid change, EBAND | `0` | `2e-10 Ry` | PASS |
| Dense total-DOS integral vs EBAND | `4.087e-6 Ry` | `2e-5 Ry` | PASS |
| Dense SU(2) q=0 evaluator rotation, electron count | `0` | `2e-13` | PASS |
| Dense SU(2) q=0 evaluator rotation, EBAND | `0` | `2e-13 Ry` | PASS |
| Two-rank MPI decomposition | same assertions | same budgets | PASS |

The production bcc Fe narrow and wide DOS fixtures use the same full `4^3`
k mesh. Both frozen-magnon probes report canonical
`EBAND=-1.98567348 Ry`; their difference is zero at the `1e-8 Ry` output
resolution. All canonical probe electron-count errors satisfy
`|dN| <= 4.744e-12` electron.

The deliberately truncated DOS window shows why it is diagnostic only:

| Production DOS diagnostic | Narrow window | Wide/converging window |
| --- | ---: | ---: |
| Total-DOS EBAND minus canonical | `+2.3161e-1 Ry` | `-1.4508e-3 Ry` |
| Projected `m0*m1` minus canonical | `+2.3161e-1 Ry` | `-1.4496e-3 Ry` |
| DOS electron-count error | `-3.9227e-1` | `+1.3294e-2` |

The dense unit oracle, rather than the still-coarse production DOS grid, is the
converged total-DOS agreement test and passes its integration-error budget.

Validation commands:

```bash
cmake --build build_13 -j2
ctest --test-dir build_13 --output-on-failure -L unit
cmake --build build_wp2_mpi --target UnitKspaceOccupations -j2
ctest --test-dir build_wp2_mpi --output-on-failure -R UnitKspaceOccupations_mpi
/Users/andersb/envs/p311/bin/python tests/run_test.py --binary build_13/bin/rslmto.x --cases-json tests/gbt_wp2/cases.json --case-name wp2_q0_theta20_narrow_dos --scratch-root /private/tmp/rslmto_wp2_final --force-serial
/Users/andersb/envs/p311/bin/python tests/run_test.py --binary build_13/bin/rslmto.x --cases-json tests/gbt_wp2/cases.json --case-name wp2_q0_theta20_wide_dos --scratch-root /private/tmp/rslmto_wp2_final --force-serial
/Users/andersb/envs/p311/bin/python tests/run_test.py --binary build_13/bin/rslmto.x --cases-json tests/gbt_wp2/cases.json --case-name wp2_q0_fixed_potential_rotation --scratch-root /private/tmp/rslmto_wp2_final --force-serial
```

The complete unit label passes `14/14`; the two-rank MPI test passes. All three
production fixtures run successfully, but process success is not the G2 physics
criterion.

## q=0 production result and gate decision

The independent dense SU(2) evaluator test is exactly invariant, but the full
fixed-potential bcc Fe q=0 rotation is not:

- collinear reference canonical EBAND: `-1.99664163 Ry`;
- globally rotated canonical EBAND: `-1.98567348 Ry`;
- absolute EBAND change: `1.096815e-2 Ry`;
- reported acoustic `omega`: `9.15373585e-2 Ry`, versus the `1e-10 Ry`
  acceptance budget.

The canonical evaluator conserves eight electrons in both spectra, so this is
not an occupation, normalization, MPI, DOS, projection, or stale-cache error.
The discrepancy enters upstream: `hamiltonian_build::ham0m_nc` constructs
already rotated spin blocks, after which the legacy
`reciprocal_fourier::fourier_transform_gbt_array` re-extracts collinear `h0/bz`
components and applies another cone construction. Repairing that representation
boundary is the gated WP3/WP4 operator work; it is not hidden by weakening the
WP2 acceptance test.

**G2E: PASS. G2O: FAIL/OPEN. Final G2: FAIL.** WP3 and WP4 are permitted only
as the documented operator gate repair. WP5 and later work must not start until
the production S/Sdot path closes G2O and final G2 is rerun.

## Completion checklist

- [x] Canonical electron count is implemented.
- [x] Canonical eigenvalue EBAND is implemented.
- [x] MPI and k-weight normalization are tested.
- [x] Frozen-magnon MFT uses canonical EBAND.
- [x] DOS settings do not move canonical EBAND beyond budget.
- [x] Electron/energy cross-checks and error budgets are reported.
- [x] G2 PASS/FAIL is stated: **G2E PASS; G2O OPEN/FAIL; final G2 FAIL**.
