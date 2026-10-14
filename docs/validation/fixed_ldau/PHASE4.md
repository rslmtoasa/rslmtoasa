# Phase 4 — fixed-U convergence criterion

The fixed-Liechtenstein SCF stop now requires the existing charge/constraint
criteria plus a full complex occupation and potential update test. For each
site, the absolute P residual is the Frobenius norm of the combined up/down
`P_out-P_in` blocks; V is evaluated nonmutating at both snapshots and its
combined Frobenius residual is reported in Ry. Relative values divide by the
larger of one and the endpoint Frobenius norms. Nonfinite, non-Hermitian, or
out-of-range occupation snapshots reset the stable-update count. Three
consecutive accepted updates are required. This is a convergence criterion,
not an integration-noise floor or an energy-functional validation.

## Verification

- Source: `source/self.f90` SHA256
  `379d9c2e787a9e9038220efabf68f9e35bd23023b57d76ca4c754a389ec46499`.
- Current-source Debug executable: SHA256
  `080d803f339461237fc4da944138567b7478f2e031632f5badf2cac91056131e`.
  `ctest --test-dir /tmp/rslmto-ldau-remediation-debug -L unit
  --output-on-failure` passed all 113 unit-labeled tests (243.90 s); full
  transcript is `/tmp/rslmto-remediate-debug-unit-379d.log`.
- The adversarial runtime run uses fixed U=1.8 eV, one SCF update, and
  deliberately loose `conv_thr=1e6`: it reports `rp_abs=1.1028e-1`,
  `rv_abs=1.4590e-2 Ry`, `streak=0`, and remains not converged. A matched
  no-U control with U/J cleared from both input and restart converges under its
  original charge criterion and emits no fixed-U residual. Exact inputs,
  seeds, and logs are in `ENERGY_ACCOUNTING.md`.
- MPI build: `/tmp/rslmto-ldau-remediation-current-mpi`, current-source Debug,
  MPI Fortran 3.1, OpenMP OFF. Binary SHA256
  `60f7951a2844798af112b575081d5c4b4f4ee9e685a57aa226edc90c2c6e84aa`.
  `UnitLdaUFixedUResidualMPI` and
  `UnitLdaUFixedUResidualMPI_mpi_2` pass. The initial full MPI `-L unit` run
  recorded 113/119 pass and six failures in
  `/tmp/rslmto-remediate-mpi-unit-final.log`: five MPI launcher socket-bind
  failures (`UnitKspaceOccupations_mpi`, the fixed-U residual two-rank test,
  and `UnitParallelContext_mpi_{1,2,4}`) plus a real communicator teardown
  defect in `UnitLdaUProjectorCommonFrozenH`. After importing the MPI context
  assignment, targeted MPI CTest of the common-H, multisite, K-space, residual,
  and parallel-context tests passed 7/7 with local socket access; log:
  `/tmp/rslmto-remediate-mpi-targeted-final.log`. The earlier failure log is
  retained as provenance. The C residual test injects imaginary P and V
  changes on one rank, checks globally reduced maxima on every rank, and rejects
  a nonphysical occupation present on only one rank.

The finite-operator projector tests cover one and two inequivalent sites for
identity metric only. They do not establish a general nonorthogonal metric
contract or material convergence.

This gate verifies the stopping implementation and its scalar-collinear MPI
reduction. It does not close the physical energy accounting requirement in
R04, nor does it establish a converged U/J material campaign.
