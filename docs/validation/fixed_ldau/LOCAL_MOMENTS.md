# Fixed-U local-moment accounting

This record follows `cleanup/RS_LMTO_ASA_LDAU_Local_Moments_Protocol.md`,
which is user-owned and intentionally remains outside this worktree. Work is
limited to its Stage M0 until the M0 gate is measured. No M1 moment correction
or energy change is permitted during this stage.

## Scope and conventions

The supported test domain is scalar-relativistic, orthogonal
(`hoh=.false.`), collinear `nsp=1`, `lmax=2` (spd), fixed Liechtenstein U/J,
with RS and reciprocal KS backends. SOC, other fixed-U spin layouts, ACBN0,
SC-U, intersite V, nonorthogonal overlap, spdf, and GPU paths are out of scope.

Use the production operator convention `P_ab = sum_n f_n psi_(n,a)
conjg(psi_(n,b))`, spherical d order m=-2..2 at offset `l*l`, and
`spin_off=9`. The energy zero is the eigenvalue frame before subtracting
`vmad`; the corrected centre, if and only if M1 passes, is then shifted by
`vmad` as production does. Energies and matrices are in Ry; input U/J may be
specified in eV. `V_in` means the type-keyed d-block and spin-block actually
inserted into the Hamiltonian, captured from `hamiltonian%hubbard_u_pot`.
The tests must distinguish `P_out` from the density used to build that
potential. Temperature, Fermi level, and occupation weights are part of every
KS manifest.

The M0 instrumentation adds the energy-weighted occupation matrix
`M1_ab = sum_n f_n epsilon_n psi_(n,a) conjg(psi_(n,b))`. It is diagnostic
only: it must not feed production `gravity_center`, `ql`, `calculate_pl`,
`newrho`, `ETOT`, or the convergence decision. Existing physical values must
remain unchanged by enabling the trace.

For each atom/type, shell and spin, the independent checker evaluates

```
t       = Tr(V_in P_out)
M1c     = M1 - c P_out
c_local = c - t/q0
q2_local= q2 - 2 Tr(V_in M1c) + Tr(V_in V_in P_out) - t*t/q0
```

The signed carrier statistic is fixed in advance as
`(SUMEV_spectral - SUMEV_local) / sum Tr(V_in P_out)`. The prediction under
test is approximately +1 within the independently measured U=0 radial
inversion/relaxation floor; the opposite numerator convention must be shown
with its sign. A value far from +1 is a measured contradiction to explain,
not a reason to alter tolerances or promote M1.

## M0 execution plan and gates

| Work item | Frozen acceptance criterion | Negative/control case | Status |
|---|---|---|---|
| M0.1 RS/KS M1 extraction | Complex matrix agrees with independent finite-H KS eigenpair sum at `5e-12 Ry` and refined analytic finite-window Lorentz RS integral at `3e-4 Ry`; same-call `Tr(M1)=q0*c` per atom/shell/spin at `1e-10 Ry` because both use identical quadrature | Conjugated/transposed and wrong-site M1 controls must mismatch by `>3e-3 Ry` (10x the accepted RS M1 oracle budget); coarse-to-fine mesh error must decrease | PASS, finite-H probes only |
| M0.2 no-U calibration | Converged nsp=1 bcc Fe in RS and true KS; record spectral `sum(q0*c)`, converged `SUMEV`, and one frozen-input-potential `newrho` replay; report charge threshold and measured floor | Verify both restart U/J arrays are zero; deterministic restart and 16³ refinement; no-U must not activate R05 P/V stopping | FAIL: KS scalar intervention exceeds the calibrated allowance; RS has no admissible floor because its run did not converge |
| M0.3 U=2 eV, J=0 carrier | Reuse only converged accepted RS/KS states that meet the existing R05 P/V criteria; independent correction arithmetic; production `calculate_pl` and one frozen-sphere `newrho` pass for legacy and local moments; report signed ratio, `q2_local>=0`, radial d charge, sphere spin moment and ASA energy components | `V=vI` gives exactly invariant q2 and `delta SUMEV=v*q0`; substituting `V_out` or `P^T` must violate the predeclared carrier criterion | NOT_RUN; stopped by the M0.2 failure |

M0 passes only when all three rows have admissible production evidence. The
previous one-step `-0.153 Ry` comparison and synthetic/calibrated radial tests
are ineligible. A missing converged state is NOT_RUN; a measured carrier
failure is reported and explained before any M1 implementation. Since M0.2
failed its scalar control, no U-active M0.3 campaign was started.

The RS independent-oracle allowance `3e-4 Ry` is preregistered from the
existing projector integral budget `1e-4` times the largest absolute energy
in its finite test window (`3 Ry`). This allowance is not reused for the
same-call trace identity, whose two scalar moments share the same grid and
therefore has the tighter `1e-10 Ry` assertion above.

The R05 Hubbard P/V residual is inactive when U is absent or U=J=0, so it is
not a no-U convergence criterion. M0.2 uses the actual charge threshold and
reports successive accepted density stability/repeatability as a diagnostic;
it must not turn on the fixed-U stopping path or modify the no-U calculation.
M0.3 uses the existing fixed-U R05 P/V convergence gate, without changing its
thresholds.

## Evidence ledger

### M0.1 finite-operator measurements

The common-H one-site and two-inequivalent-site probes passed 2/2 on the
frozen Debug source snapshot. Exact source/test/executable/log hashes,
commands, predeclared limits, and machine-readable values are in
[`evidence/m0_1_finite_operator.json`](evidence/m0_1_finite_operator.json).
The RS 6000-channel analytic Lorentz M1 errors were `2.894374e-7 Ry` and
`2.845306e-7 Ry`; KS exact-eigenpair M1 errors were `7.757919e-18 Ry` and
`2.797158e-17 Ry`. The maximum same-call `Tr(M1)-q0*c` residual was
`1.054712e-15 Ry`. The 750/1500/3000/6000 channel error sequences decreased
strictly. M1 conjugation/transpose mismatches were `0.069864–0.072836 Ry`,
and the two-site production site-swap mismatches were `0.047631–0.047690 Ry`,
above the preregistered `0.003 Ry` minimum.

The same frozen source/tests also passed 4/4 in an MPI-enabled Debug build:
serial one-site and multisite tests plus their two-rank variants. The MPI
result validates the M1 real/imaginary reduction for these probes; it does
not establish behavior in a full RS/KS production campaign. Exact build and
test executable hashes are recorded in the JSON evidence.

This is class-B finite-operator evidence only. It does not measure M0.2 or
M0.3 and does not pass the overall M0 gate.

The final diagnostic source snapshot is `self.f90` SHA256
`d436f4a21f3a0c2d8e1a06c68db1f3334455b254374dbd4b6262488a96bdfd7f`.
The serial Debug executable SHA256 is
`73a6b36fb34656fecc9515740c3c721078017704ab601716b4637bfa55a42ecc`;
the fresh MPI Debug executable SHA256 is
`59ed8719fc1ea350ec12889ce64645e27107d77fab670f0a6ce59bea105ad863`.
On that exact source, serial `ctest --test-dir /private/tmp/rslmto-localmom-debug
-LE 'validation|performance' --output-on-failure` passed 119/119 in 331.24 s
(log SHA256 `1d8768ba1f61cc11d308098e95771856bb82dc89613819e82551071a2184289d`).
The focused serial LDA+U tests passed 5/5. MPI-enabled focused residual and
one-/two-site M1 projector tests passed 3/3 in 21.32 s with
`ctest --test-dir /private/tmp/rslmto-localmom-mpi-current -R
'UnitLdaUFixedUResidualMPI_mpi_2|UnitLdaUProjector(CommonFrozenH|MultisiteFrozenH)_mpi_2'
--output-on-failure` (log SHA256
`a3aefa1d99c1dd2d006da23dc6b26aba9971826825b780adf14fe547e291e6d7`). MPI trace tests also
passed both a two-rank root-only CSV write and an invalid output-path case;
both produced the same `Fe_out`, and failed trace-file creation disabled
tracing without rank divergence. These checks validate the diagnostic path
and its collectives, not the M0 radial-energy hypothesis. The RS positive
case used input/seed prehashes
`9a0ffae2bbdde7ddf1c2b0b8396f5745a3dd944c5aaa71550d01436d9ce36a67` and
`fcbf9541ca8c60940f9662246fc7acb841a44b7832b141241e32e8fa9caa7df1`, CSV
SHA256 `07a1407ee5008558892b5fbb30bca25fe62765fd29cb52dbbaed5079428aa3bb`,
and log SHA256 `01d38893c7ed45dfc0fa2da8ad3dfa33fae8ba2868e8f8e86c1fd0505c2edf32`;
the invalid-path log SHA256 is
`7e5fbb5c48ef3392cdfad12919348d24bb3f64885940969127e8116446229953`.
Both have the same `Fe_out` SHA256
`7de5fcc993a8d67ca4c7473bbe7e64cfa93255076541ce541f6a487fb7b02909`.
For true KS, rolling-final selection (`RSLMTO_FIXED_U_TRACE_ITERATION=0`)
captured selected iteration 2 from an input with `use_kspace=.true.`; the CSV
contains 1,236,242 rows, 70 M1 rows, 1,728 k weights, 31,104 eigenvalues,
and 559,872 eigenvector/H rows. Its trace/plain `Fe_out` hashes match at
`be13516c…701b3`. Full per-run provenance is retained under
`/private/tmp/rslmto-localmom-mpi-open-smoke/{success,open_fail,ks_rolling}`.
Open MPI required approved local socket access; sandbox-only launches failed
before reaching the application.

### M0.2 no-U calibration outcome

The KS M1/eigenpair comparison, same-call `Tr(M1)=q0*c`, and exact frozen
NEWRHO replay checks pass. The preregistered deterministic repeat and 16³
refinement give a baseline `D = sum(q0*c)-SUMEV = -5.028166848219939e-6`
Ry, repeat `D=-5.023148958116508e-6` Ry, refinement
`D=-5.053452480030529e-6` Ry. Thus the no-U floor is `5.028166848219939e-6`
Ry, the separately measured uncertainty is `2.528563181058985e-8` Ry, and
the acceptance allowance is `5.053452480030529e-6` Ry.

For the frozen scalar intervention, `c_local=c_raw-0.01` Ry with `q0`,
`q2`, radial arrays, and seed held fixed. The expected raw-minus-local
`SUMEV` change is `+0.06498480421874168` Ry; production gives
`+0.06497426910958881` Ry. The absolute residual is `1.0535109152873257e-5`
Ry, or `2.08473498` times the frozen allowance. This is a measured KS
scalar-control **FAIL**, not an M1 quadrature error; the finite-H M1 tolerance
does not apply to this radial carrier test. Hash-bound details, including
the source, executable, input/seed, raw CSV, manifest, checker, and floor
hashes, are recorded in
[`evidence/m02_ks_scalar_v0_fail.json`](evidence/m02_ks_scalar_v0_fail.json).

The corresponding RS no-U run reached its iteration limit without satisfying
the `1e-9` charge criterion, so no accepted RS floor or carrier verdict is
available. The stopped continuation is explicitly cancelled and excluded.
The exact-source diagnostic plumbing smoke reproduced same-call LIVE_STATE
components, but its RAW/LOCAL comparison is not a converged M0 result. The
M0.2 failure blocks M0.3; no U-active carrier, M1 correction, or physical
energy correction was attempted.

### Frozen 3d7571a baseline

The source worktree `/private/tmp/rslmto-ldau-remediate-worktree` was clean at
`3d7571ad4ac82d1c8c48e82c44299290507a787b` before baseline configuration.
Fresh builds used GNU Fortran 16.2, Accelerate BLAS/LAPACK, unit tests enabled,
and MPI/OpenMP disabled. The Debug executable
`/private/tmp/rslmto-localmom-base-debug/bin/rslmto.x` has SHA256
`0e7a10a13289809605105154c9f33ab553020ee2b250c0c35aefdf04e476dd7a`.
`ctest -L '^quick$' --output-on-failure` passed 24/24 in 88.19 s (log SHA256
`7432305334c23ba4f1330d1300a8e267b8886f5157b6680caa50f27245375607`);
`ctest -LE 'validation|performance' --output-on-failure` passed 119/119 in
189.64 s (log SHA256
`f4a6a4d2bd3fa13f8683ce2e5559acf79b9a3b2b9cbff23843918ea485fa5bbc`).
The fresh Release executable
`/private/tmp/rslmto-localmom-base-release/bin/rslmto.x` has SHA256
`ca58fef96fe947a165d7ddc94b7ad9abd5ec7d6fdb7921c515d604a9885cf862`;
`ctest -L unit --output-on-failure` passed 114/114 in 97.07 s (log SHA256
`dcda8922ce87729f188d27a2d18031b18841387fed0bfafafba35ead6833ac42`).
All build/test outputs remain in `/private/tmp`; the original user worktree and
its `cleanup/` documents were not edited.

No M0 measurement is recorded yet. For every later result, preserve the exact
source and executable SHA256, command, immutable pre-run deck and restart
seed hashes, raw-output location/hash, selected iteration, measured value,
predeclared tolerance, and a negative control. Large raw data stays in
`/private/tmp`; compact machine-readable summaries belong under
`docs/validation/fixed_ldau/evidence/`.
