# Fixed-U local-moment accounting

This record follows `cleanup/RS_LMTO_ASA_LDAU_Local_Moments_Protocol.md`,
which is user-owned and intentionally remains outside this worktree. The
user-approved Protocol amendment A2 below supersedes A1.1 criterion (i) and
the floor-only M0.3 allowance; A1 remains the source of the other registered
M0 shape tests. Current work is still at Stage M0. No M1 moment correction or
production energy change has been made. A2 permits Stage M1 for KS only if
every required M0 gate passes, after which work stops at Gate M1.

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

### A1.1 measurement pre-registration

Before any new shift result is inspected, define the frozen-state difference
as `DeltaSUMEV(v) = SUMEV_RAW - SUMEV(c_raw-v)` for positive `v`, with the
same definition applied to signed `v`; therefore the expected first-order
response is `v*q0`. Compute `r(v)` from the measured difference and captured
`q0`, never from a producer-reported residual. All eight signed cases must
start from identical copied radial arrays, EC/EV seeds, `q0`, `q2`, and other
`calculate_pl` inputs; only the center may change, by exactly `-v`. Keep each
site/spin/type mapping explicit and reject missing or duplicate cases.

For each site/spin, calculate `e(v)=s(v)/q0-1` and fit the four positive
magnitudes to the unweighted even model `a0+a2*v^2+a4*v^4`. The extrapolation
criterion is the user-specified `abs(a0) < 1e-6`; report prediction changes
from a linear-in-`v^2` fit and leave-one-magnitude-out quartic fits as
sensitivity diagnostics, not as additional pass criteria. Measured
`abs(e)` must strictly decrease as `|v|` decreases through 0.02, 0.01, 0.005,
0.0025 Ry. For signed `r`, report the all-point unweighted cubic fit
`r3(v)=a0+a1*v+a2*v^2+a3*v^3` and the unweighted signed quadratic fit
`r2(v)=b0+b1*v+b2*v^2`. The primary log-log OLS fit is `log|r|` versus
`log|v|` over all eight signed samples; separate positive/negative fits and
adjacent pairwise orders are diagnostics. An exact zero makes the log fit
undefined and must be reported for A1.3 review, without a pseudocount. Require
`|r|` to strictly decrease as `|v|` halves on both sign branches, and report
all fit residuals and alternative-fit sensitivity. Assess a v-independent
plateau using the fitted constant together with the full measured `r(v)`
sequence and its shape; the criterion is no fitted plateau above `10*D0`,
where `D0` is the already measured KS no-U floor
`5.028166848219939e-6 Ry`. There is no additional cap on the smallest-shift
`|r|`. If M0.3 is reached, the
preregistered allowance is the fitted residual magnitude evaluated at each
actual per-spin shift `t_s/q0_s`, summed once over spins, plus the existing
total `D0` floor. Evaluate and label model extrapolation if a shift is outside
the measured `[-0.02,0.02] Ry` range; do not silently treat it as
interpolation. Report fit sensitivity separately, not as an added allowance
or gate. If measurements expose a wrong error model, stop under A1.3 and
propose a corrected criterion for user approval. These analysis choices may
not be changed after viewing the shift measurements.

| Work item | Frozen acceptance criterion | Negative/control case | Status |
|---|---|---|---|
| M0.1 RS/KS M1 extraction | Complex matrix agrees with independent finite-H KS eigenpair sum at `5e-12 Ry` and refined analytic finite-window Lorentz RS integral at `3e-4 Ry`; same-call `Tr(M1)=q0*c` per atom/shell/spin at `1e-10 Ry` because both use identical quadrature | Conjugated/transposed and wrong-site M1 controls must mismatch by `>3e-3 Ry` (10x the accepted RS M1 oracle budget); coarse-to-fine mesh error must decrease | PASS, finite-H probes only |
| M0.2 no-U calibration | Converged nsp=1 bcc Fe in RS and true KS; record spectral `sum(q0*c)`, converged `SUMEV`, and one frozen-input-potential `newrho` replay; report charge threshold and measured floor | Verify both restart U/J arrays are zero; deterministic restart and 16³ refinement; no-U must not activate R05 P/V stopping | Initial scalar-v0 verdict superseded by A1 as an ill-posed criterion. A1.1 scaling test NOT_RUN; RS has no admissible floor because its run did not converge |
| M0.3 U=2 eV, J=0 carrier | Reuse only converged accepted RS/KS states that meet the existing R05 P/V criteria; independent correction arithmetic; production `calculate_pl` and one frozen-sphere `newrho` pass for legacy and local moments; report signed ratio, `q2_local>=0`, radial d charge, sphere spin moment and ASA energy components | `V=vI` gives exactly invariant q2 and `delta SUMEV=v*q0`; substituting `V_out` or `P^T` must violate the predeclared carrier criterion | NOT_RUN; gated on A1.1 passing and its fitted error model/allowance being locked before U-state outcomes |

M0 passes only when all three rows have admissible production evidence. The
previous one-step `-0.153 Ry` comparison and synthetic/calibrated radial tests
are ineligible. A missing converged state is NOT_RUN. The original scalar-v0
criterion is superseded by A1; no U-active M0.3 campaign has been started.

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

For the historical frozen scalar intervention, `c_local=c_raw-0.01` Ry with `q0`,
`q2`, radial arrays, and seed held fixed. The expected raw-minus-local
`SUMEV` change is `+0.06498480421874168` Ry; production gives
`+0.06497426910958881` Ry. The absolute residual is `1.0535109152873257e-5`
Ry, or `2.08473498` times the then-frozen allowance. Under user-approved A1,
this former verdict is reclassified as **criterion ill-posed, superseded by
A1**: the zero-shift floor does not bound smooth, shift-dependent
`calculate_pl`/potential-parameter inversion error. The finite-H M1 tolerance
does not apply to this radial carrier test. Keep the original hash-bound details, including
the source, executable, input/seed, raw CSV, manifest, checker, and floor
hashes, are recorded in
[`evidence/m02_ks_scalar_v0_fail.json`](evidence/m02_ks_scalar_v0_fail.json).

The original evidence file is preserved unchanged as historical evidence. The
new A1.1 signed-shift scaling run used the same Release build/source, PRE deck,
and seed as the accepted M0.2 no-U KS state, and reached its same iteration 23.
Before inspecting shifts, the independent state check passed with maximum
differences of zero for P, M1, q0/c/q2, calculate_pl inputs/returns, and frozen
radial/replay inputs. Candidate CSV SHA256 is
`55ef74516393256f262de6f02ecf0f061536ff9bff62641d1e6e828149586ac4`; executable
SHA256 `7a98a68a11f6a4a062e1f8256ffb80a64b4dc8967c0165dab1ffb8910425a3fa`;
source self.f90 SHA256 `5f9cf847738803a49fc8d35c51be7b5dd3e7083c851a5c07d1785e07bc1dab89`.
The hash-bound equivalence record is
[`evidence/a11_accepted_state_equivalence.json`](evidence/a11_accepted_state_equivalence.json)
(SHA256 `1a7238ed4f98649feb8e654c04209b987d80cd715d03eea103b64b8a0a7ef005`).

The A1.1 fit reports total q0 `6.498480421874168`. Its symmetric-slope
relative errors at |v|=`0.0025, 0.005, 0.01, 0.02` Ry are respectively
`-1.76935e-5, -4.64255e-5, -1.62889e-4, -6.27649e-4`; the fitted v→0
intercept is `-7.86514e-6`, outside the required `1e-6`. The spin-up and
spin-down intercepts also fail, at `-7.58154e-6` and `-8.43595e-6`.
The total all-eight-point log-log order for |r| versus |v| is `2.72568`; the
signed residual branches are monotone and the fitted cubic constant is
`-5.11131e-8 Ry`, below `10 D0 = 5.02817e-5 Ry`. The total residual at
`v=+0.01 Ry` is `-1.0535109153e-5 Ry`. The historical verdict was
**FAIL_A1_1** because the extrapolated symmetric-slope error missed the
declared limit. Under user-approved A2, that limit is reclassified as a
criterion without an error model and is superseded; the original measurement
and JSON files remain unchanged. The measured spin intercepts are
`a0_up=-7.58154e-6` and `a0_down=-8.43595e-6`, giving carrier efficiencies
`1+a0` of approximately `0.999992` for both spins. The residual shape and
power-law checks passed, with no plateau; A2 replaces the unweighted relative
slope cutoff with a decision-relevant comparison against the U=0 floor. No A2
outcome has been inspected yet.

The complete signed r(v), slopes, spin-resolved fits, diagnostics, and run
provenance are in
[`evidence/a11_scaling_result_reviewed.json`](evidence/a11_scaling_result_reviewed.json)
(SHA256 `48b0cb79291b12fd5e3c027b92379f0f70b45e8418882348616a1ea1abed27ae`)
and the concise table is in
[`evidence/a11_scaling_summary.md`](evidence/a11_scaling_summary.md)
(SHA256 `91c4c9c89e4e206ea07cb0e9a9654ad5e6d51d128f650c94f9854ac10c1c2792`).
The earlier `a11_scaling_result.json` remains unchanged as the initial result
record; the reviewed artifact includes the explicit binary/source hashes and
aggregate-spin diagnostics.

Under A1, **M0.3 was not run**. No U-state carrier result, M1, M2, physical
energy correction, or SOC+U conclusion follows from the A1 evidence. The
corresponding RS no-U run reached its 120-step limit without satisfying the
`1e-9` charge criterion, so no accepted RS floor or carrier verdict is
available. Its full
residual history is retained in
[`evidence/rs_u0_residual_history_v0.json`](evidence/rs_u0_residual_history_v0.json)
(source run log SHA256 `088f83094298f66fd900cc6452ba64102e317e9ac19908868dd39eaf0373d9e1`,
artifact SHA256 `24ddd3c743d7c6780167e366876fbc47d764cffca10c48e05df75f115f54f4a5`).
The final residual was `1.57e-8`; the last-ten range was `2.9e-9` to
`1.81e-8`. No RS floor or stagnation threshold is inferred from this finite
history. The A1.2 restart continuation with more iterations and Broyden mixing
is a separate, not-yet-run work item; it is not part of the KS A1.1 gate.
The exact-source diagnostic plumbing smoke reproduced same-call LIVE_STATE
components, but its RAW/LOCAL comparison is not a converged M0 result. M0 was
**FAIL** under the superseded KS A1.1 criterion, with M0.3 **NOT_RUN**; the
A1 result does not validate a physical total-energy correction.

### Final-source verification

The A1 implementation source was frozen at self.f90 SHA256
`5f9cf847738803a49fc8d35c51be7b5dd3e7083c851a5c07d1785e07bc1dab89`.
The Debug binary SHA256 was
`87a939337a8773d0b9c8363c11fe123a94d4401ebae47f3648db4e6bb30925cb`; on that
source, `ctest --test-dir /private/tmp/rslmto-localmom-debug -L unit
--output-on-failure` passed 114/114 (log SHA256
`9e6a8be399fc7eda7c608a187e77d0c2bfabbf2dc6d77f99195fe35b763c5f41`), and the
required stage-end `ctest --test-dir /private/tmp/rslmto-localmom-debug -LE
'validation|performance' --output-on-failure` passed 119/119 (log SHA256
`6b7dda0ee81d4919b5f8acbfdd019750df25f52bdbb993cfb620e914117e6e00`). Both
logs are retained under `/private/tmp`; no production source changed after
these runs.

## Protocol amendment A1 (verbatim)

The following user-approved amendment is reproduced verbatim. It supersedes
the original scalar-control criterion, but does not edit or replace the
user-owned protocol file or the historical scalar-v0 evidence artifact.

```text
Protocol amendment A1 (user-approved, supersedes the M0 scalar-control criterion
in cleanup/RS_LMTO_ASA_LDAU_Local_Moments_Protocol.md; record this amendment
verbatim in LOCAL_MOMENTS.md and do not edit the user-owned protocol file).

Rationale: the zero-shift floor D measures the radial inversion offset at v=0.
A finite centre shift adds a smooth, shift-dependent calculate_pl/potential-
parameter inversion error that this floor cannot bound. The measured KS result
(ΔSUMEV ratio 0.99984 at v=0.01 Ry) is therefore not a refutation, and the
previous verdict is reclassified from FAIL to "criterion ill-posed, superseded
by A1". Keep the original evidence file unchanged and add the new result
alongside it.

A1.1 Scalar-shift scaling test (KS, frozen converged no-U state, same seed and
frozen-potential single newrho pass as before):
  - v in {±0.0025, ±0.005, ±0.01, ±0.02} Ry, applied to c with q0, q2 held fixed.
  - Report r(v) = ΔSUMEV(v) − v·q0 for each v.
  - Report the symmetric slope s(v) = [ΔSUMEV(+v) − ΔSUMEV(−v)]/(2v).
  - PASS if both hold:
      (i) |s(v)/q0 − 1| decreases with v and its v→0 extrapolation is below
          1e-6 relative;
      (ii) r(v) is smooth, and a log-log fit gives order ≥ 2 with no
          v-independent plateau above 10× the U=0 floor.
  - If r(v) is flat or non-monotonic, stop and investigate the radial solver
    tolerance (TOLRSQ/RSEQSR) before continuing. Report what you find; do not
    tune the tolerance to pass.
  - The fitted r(v) model, evaluated at the actual per-spin shifts t_s/q0_s of
    the U state plus the U=0 floor, is the preregistered allowance for M0.3.

A1.2 Backend decoupling:
  - KS proceeds through M0.3, M1 and M2 independently.
  - RS becomes a separate work item:
      * record the full residual history of the no-U run;
      * continue from restart with more iterations and Broyden mixing;
      * report the stagnation level.
  - If 1e-9 is unattainable, use the measured stagnation level as RS's
    convergence threshold for its own floor and Hellmann–Feynman noise budget.
  - RS gates are reported separately and must not be marked PASS by
    inference from KS.

A1.3 Governance: if any later criterion is found to rest on a wrong error
model, stop, write a short justification with the measurement that exposes
it, and propose a corrected criterion for user approval. Never self-amend.

Continue with A1.1 and then M0.3 (KS) on fix/ldau-local-moments. Stop at the
KS M0 gate and report.
```

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

At initial preregistration no M0 measurement was recorded. For each result,
preserve the exact source and executable SHA256, command, immutable pre-run
deck and restart seed hashes, raw-output location/hash, selected iteration,
measured value, predeclared tolerance, and a negative control. Large raw data
stays in `/private/tmp`; compact machine-readable summaries belong under
`docs/validation/fixed_ldau/evidence/`.

## Protocol amendment A2 (verbatim)

~~~text
Protocol amendment A2 (user-approved; supersedes A1.1 criterion (i) and the
M0.3 allowance rule; record verbatim in LOCAL_MOMENTS.md; do not edit the
user-owned protocol file or any historical evidence).

Rationale: the A1.1 limit |a0| < 1e-6 had no error model and protected no
decision. M0 exists to decide whether radial SUMEV carries a centre shift
accurately enough to discriminate the energy candidates (a), (b) and (c),
which differ by Tr(V_in P) ≈ 0.25 Ry at U = 2 eV. Measured A1.1 values:
  - linear carrier efficiency 1 + a0 = 0.999992, the same sign for both spins;
  - smooth residuals of order 2.73, with no plateau.
Reclassify FAIL_A1_1 as "criterion (i) without an error model, superseded
by A2". Keep the existing evidence files unchanged.

A2.1 Decision-relevant replacement for A1.1(i):
  Using the accepted converged KS U = 2 eV, J = 0 state, record t_s = Tr(V_in P_out)
  for each spin. PASS if sum_s |a0_s|·|t_s| ≤ D0 (5.028e-6 Ry). Recompute this
  from the existing A1.1 fit coefficients; do not refit. Criterion (ii) stands
  as already passed.

A2.2 Extended scan so that M0.3 is an interpolation, not an extrapolation:
  - Same frozen no-U KS state, same pre-registered analysis code, and the same
    fitting models as A1.
  - Add v ∈ {±0.04, ±0.08, ±0.12} Ry to the existing eight shifts.
  - Before running M0.3, confirm that every actual per-spin shift t_s/q0_s of
    the U state lies inside the scanned range.
  - The M0.3 allowance is sum_s |r_fit,s(t_s/q0_s)| + D0, where r_fit
    includes its linear term.
  - If the extended residuals break the smooth-power-law shape (a plateau, a
    sign change, or non-monotonic behaviour), stop under A1.3.

A2.3 M0.3 (KS): run as specified in the protocol with the A2.2 allowance. Report:
  - the signed carrier ratio;
  - q2_local ≥ 0;
  - the d radial charge;
  - the sphere spin moment;
  - all ASA energy components;
  - the scalar-shift, V_out and P^T controls.

A2.4 Non-blocking diagnostic, reported but not gated: identify the origin of a0.
  Map D(E) → EVAL directly: for E near E_nu, form D from the POTPAR
  parameters exactly as calculate_pl does, solve the radial equation, and
  report f(E) = EVAL − E and f'(E_nu). Check whether f'(E_nu) is consistent
  with the φ̇ numerical derivative accuracy or with the known GFAC
  normalization mismatch (~2e-4). Do not change production code here.

A2.5 RS (separate work item, not blocking KS):
  - Run the A1.2 restart continuation with Broyden mixing.
  - If the residual stagnates at about 1e-8, adopt the measured stagnation level
    (e.g. 2e-8) as the RS convergence threshold.
  - Measure the RS D0 at that threshold and record it. RS M0 is reported
    separately.

A2.6 Criterion governance, which applies to all future gates: every
  preregistered tolerance must state (1) the decision it protects, (2) the
  quantity it is compared against and why that magnitude matters, and (3) its
  error model. A proposed criterion lacking any of these goes to the user for
  review before any measurement is inspected.

If A2.1 passes and M0.3 passes, M0 (KS) is PASS: proceed to Stage M1 for KS
and stop at Gate M1. If M0.3 fails, stop and report.
~~~

## A2 execution plan and locked decisions

This section records how the approved amendment will be applied. It does not
change the A2 text above, the user-owned protocol, or any historical evidence.
No A2 outcomes were inspected when this plan was written.

### A2.1 accepted-state intercept test

Use only a fully converged KS `U=2 eV, J=0`, `nsp=1`, `lmax=2` state that
satisfies the existing charge and fixed-U P/V convergence gates. The selected
capture must be the accepted final iteration and must include the actual
Hamiltonian `V_in` and its spectral `P_out`; the matrix residuals and input
provenance will be checked before calculating the criterion. A prior
one-update trace seeded from a converged `Fe_out` is not an accepted-state
capture. Use the exact existing A1 spin intercepts, without refitting:
`a0_up=-7.58154e-6` and `a0_down=-8.43595e-6`. Independently compute
`t_s=Tr(V_in,s P_out,s)` from the captured matrices, then
`B_A2.1=|a0_up|*|t_up|+|a0_down|*|t_down|` in Ry. The gate passes exactly
when `B_A2.1 <= D0`, where `D0=5.028166848219939e-6 Ry` (reported also as
`5.028e-6 Ry` per the amendment). This protects whether the measured A1
intercept bias at the actual U-state shifts can exceed the no-U radial
inversion offset; the comparator is that independently measured floor, whose
magnitude is the baseline systematic in `SUMEV`. The error model is the fixed
A1 intercept multiplied by each measured `|t_s|`; no refit or extra
uncertainty band is added. If the accepted state or actual `V_in/P_out` pair
cannot be established, A2.1 is not assessed and later KS energy gates do not
start.

### A2.2 extended no-U scan

Only after A2.1 passes, repeat the A1 frozen-state no-U capture from its exact
accepted pre-run deck and seed and add the six signed shifts
`±0.04, ±0.08, ±0.12 Ry` to the original eight. Preserve the original A1
analyzer and its outputs unchanged; the new analysis is separately
hash-bound. Use the same models on all fourteen signed points: the signed
cubic `r_fit(v)`, all-point OLS fit of `log|r|` against `log|v|`, and the
existing even fit for symmetric-slope diagnostics on the seven positive
magnitudes. The smooth-power-law check reuses the A1 criteria: on each sign
branch, no sign reversal and strict monotonic decrease of `|r|` as `|v|` is
halved; the all-point log-log order is at least 2; the fitted constant term
does not exceed `10*D0`. Exact zero residuals make the log fit undefined and
fail closed for review. These are the existing A1 limits, not new cutoffs;
`D0` is not treated as a per-point noise band. A failure of any shape check
stops before M0.3 under A1.3.

Before M0.3, establish from the accepted converged U-state capture that each
`t_s/q0_s` lies within `[-0.12, +0.12] Ry`. The scan decision is whether the
same fitted response model is supported across the full interval containing
those target shifts. The fit’s signed residual sequence, branch behavior,
order, and intercept are the evidence; no tolerance is enlarged after seeing
the results.

### M0.3 KS allowance and controls

For each spin, evaluate the locked fourteen-point signed cubic residual model
at `x_s=t_s/q0_s`. The total energy allowance is exactly
`A=|r_fit,up(x_up)|+|r_fit,down(x_down)|+D0` Ry, with no additional fit
sensitivity term. The principal carrier comparison is
`|DeltaSUMEV - (t_up+t_down)| <= A`, where
`DeltaSUMEV=SUMEV_spectral-SUMEV_local`; report the signed ratio
`DeltaSUMEV/(t_up+t_down)` as well. This tests whether the local radial
replay removes the measured Hubbard first-moment contribution within the
shift-dependent no-U response plus the independent floor.

Run production `calculate_pl`/single frozen-sphere `newrho` RAW and LOCAL
cases from the same accepted U state, plus the scalar-shift, `V_out`, and
transposed-P controls. The scalar control must preserve `q2` algebraically and
its SUMEV response must meet the same registered allowance. The original
protocol requires `V_out` and transposed-P to violate the carrier criterion;
they count as discriminating only if their absolute deviation from the
correct `V_in/P_out` target exceeds `A`. A nondiscriminating control is not a
PASS and leaves the corresponding overall M0 gate incomplete. Report, without
adding fitted thresholds, `q2_local` (which must be nonnegative), d radial
charge, sphere spin moment, every ASA energy component, and all three
controls. The evidence must record the actual frozen density/potential state,
R05 residuals, input/restart hashes, source and executable hashes, and raw
capture hashes.

### Scope after M0

If A2.1, A2.2, and every required M0.3 gate/control pass, proceed through
Stage M1 for KS and stop at Gate M1; do not start M2. The A2 allowance above
also replaces the old M0.3 floor-only allowance for any Stage M1 KS check.
Otherwise stop at the failed or unassessed gate. The A2.4 radial-inversion
diagnostic is private and nonblocking; it may not change production code.
The A2.5 RS Broyden continuation and RS floor are reported separately and
cannot establish or refute a KS result.

The RS continuation starts from the existing 120-step no-U RS endpoint
`Fe_out.nml` (SHA256
`bb09ff1734cd50a24da7f1f99b64e3d4485d907f7f96abf6dcbfb2b90d134ad2`). Its
original physical input has `conv_thr=1e-9`, Broyden mixing, 2,500 channels,
300 K, `nsp=1/lmax=2`, and no U/J; the pre-run input hash is
`bdb8f3723e4e028bdb754015a8aa5117d0f647348d026c25f83c3fc107f0e53c`. The
continuation preserves those physical settings and the requested `1e-9`
charge threshold, raises only `nstep` to 240, and uses a separate scratch
directory and immutable Release executable. The full residual history is
retained. If the original threshold is reached, use that accepted state. If
not, report the full continuation history and measured tail; only a plateau
near `1e-8` justifies the amendment's RS-specific measured stagnation
threshold (about `2e-8` as its example). That threshold applies only to the RS
floor and cannot be carried over to KS. Measure and label the RS floor
separately.

### A2.1 and A2.5 measured execution

A2.1 passes on the accepted converged KS `U=2 eV, J=0` capture. The selected
iteration is 32; charge residual is `7.69e-8`, and the final three fixed-U
residuals satisfy the registered three-update streak. Independent full-matrix
reconstruction gives `t_up=-0.27138788113273177 Ry` with
`q0_up=4.527615147715704`, and `t_down=+0.016746867166353676 Ry` with
`q0_down=2.08410054869144`. Thus the target shifts are `-0.059940580698351495`
and `+0.00803553704588229 Ry`. The recomputed A1-intercept bias is
`B_A2.1=2.19881406993552e-6 Ry`, below `D0=5.028166848219939e-6 Ry`
(ratio `0.437299`). The hash-bound source and accepted-state provenance are
recorded in `/private/tmp/rslmto-localmom-a2-ks-u2-full/manifest.final.json`.

The first A2.2 SCF recapture used the exact no-U accepted PRE deck and seed,
but it is **not an accepted frozen-state scan**. It stopped at the
pre-outcome state-equivalence gate: the reference accepted iteration was 23
and the recapture selected iteration 25. Before inspecting any shifts, the
independent reviewer found differences in directly captured frozen operands:
maximum `|ΔPL input|=1.7177e-9`, `|Δ radial moment|=4.0072e-9 Ry`,
`|ΔV|=3.6679e-10 Ry`, `|Δrho_in|=3.6100e-9`,
`|ΔEVAL/energy seed|=2.6946e-9 Ry`, and `|ΔQL|=1.4969e-9`; mesh and index
sets match. The capture-only run with A2 scanning disabled produced exactly
the same `Fe_out` hash as the rejected fourteen-shift run, so the appended
private replays did not alter the converged SCF output. No shift values or
fits from that rejected capture have been used. The accepted state was later
rehydrated into private clones with all required replay operands hash-bound;
the original state-equivalence condition was retained, and no additional SCF
run was used to force a match. The two incomplete/rejected capture manifests
are `/private/tmp/rslmto-localmom-a2-no-u-14scan/manifest.json` and
`/private/tmp/rslmto-localmom-a2-parity-captureonly/manifest.json`.

Under A2.5, the 240-step RS continuation did not reach its requested
`1e-9` threshold. Its preceding 60-step window (121–180) ranged from
`2.3e-9` to `2.53e-8` (mean `9.4517e-9`); the final window (181–240) ranged
from `2.4e-9` to `2.67e-8` (mean `9.06e-9`), with a small positive trend.
The locked one-significant-digit upward-rounding rule therefore selected the
RS-only threshold `3e-8`; this describes a low-residual plateau, not
convergence at `1e-9`. A separate restart from that endpoint converged in two
steps at the registered threshold. The complete 240-step residual history and
tail calculations are preserved in
[`evidence/a2_rs_floor_3e8.json`](evidence/a2_rs_floor_3e8.json) (SHA256
`187772f351330001555e7669a0ffb856d528e7bc2410a4399a71aef6991795df`). At
selected iteration 2,
`ΣTr(M1)=-1.9975997502219833 Ry` and `SUMEV=-1.9975946937217053 Ry`, giving
signed `D=-5.0565002780e-6 Ry` and `D0=5.0565002780e-6 Ry`. This is a single
RS floor measurement; its repeat/refinement uncertainty is unmeasured and it
does not pass or imply any KS gate. Provenance is in
`/private/tmp/rslmto-localmom-a2-rs-floor3e8/manifest.json`.

### A2.2 accepted-state replay and extended fit

The original accepted no-U KS iteration-23 state was rehydrated from its
hash-bound CSV into private clones and replayed through the production
`calculate_pl` and unchanged `NEWRHO` routines. The zero-shift replay matched
the archived PL, radial density, SUMEC, SUMEV, and EV vector exactly. All
eight historical A1 shifts also matched their archived PL and per-spin SUMEV
exactly. The independent reviewer confirmed these state/replay gates before
the six appended shifts were run. This is a fixed-state radial replay, not a
new SCF run. The baseline/A1 report is
[`evidence/a2_iter23_raw_a1_replay.json`](evidence/a2_iter23_raw_a1_replay.json)
(SHA256 `549edf2183ce13fbec7556f2bcde3d23b4f1c78eb908c04e616641a4d694d179`),
and the six new case outputs are
[`evidence/a2_iter23_extra_six.json`](evidence/a2_iter23_extra_six.json)
(SHA256 `93236aeaa5c85afcbf8c0075104f78462888c45cdbcdb2a70b51fa066c86fc79`).
Build, executable, fixture, source and command provenance is in
[`evidence/a2_iter23_replay_provenance.json`](evidence/a2_iter23_replay_provenance.json)
(SHA256 `362ec22c98a488d292c7d5239fbc2269470e37e1af78fceace40982350633d45`).

The independently reviewed fourteen-point fit passes the existing residual
shape criteria for both physical spins and the aggregate diagnostic: each
signed branch has a constant nonzero sign, `|r|` decreases strictly toward
zero, the all-point log-log order is above 2, and the fitted constant is below
`10*D0`. The A2.1 target shifts remain inside the scanned interval. The
hash-bound result is
[`evidence/a2_iter23_extended_fit.json`](evidence/a2_iter23_extended_fit.json)
(SHA256 `0ad6a7703aa26eb132602e7a43773d7619f70651c6b90817e4f945ff259cb993`).
It gives the locked true-target M0.3 allowance `A=0.0017022383465738291 Ry`
and the scalar `v=+0.01 Ry` allowance `6.022982854018227e-5 Ry`. The
expected scalar SUMEV change from the accepted U-state spin populations is
`0.06611715696407144 Ry`.

### M0.3 accepted KS frozen-state replay

The accepted iteration-32 U=2 eV, J=0 capture passed the registered charge
and R05 gates and its LIVE_STATE baseline reproduced captured PL, density,
SUMEC, SUMEV, and EV within their preregistered limits. The five replay cases
(RAW, LOCAL, scalar +0.01 Ry, VOUT, and transposed-P) use the same accepted
capture and the production `calculate_pl`, `NEWRHO`, `POISS0`, and `VXC0SP`
fixture seams. The captured Hamiltonian `v_hub/P_OUT` rows remain `V_in`.
For the deliberately wrong `V(P_out)` pairing, the adapter independently
reconstructs the d-block from the fixed Liechtenstein J=0 expression
`U_Ry*(I/2-P_sigma)`, using the manifest-bound U=2 eV and the source
`ry2ev=13.605703976`. Its per-spin `t`, corrected center, `q0`, and `q2_local`
agree with the separate production VOUT_CONTROL snapshot within `2e-12`.
The adapter and reports therefore distinguish the true VOUT matrix from
`V_in` without changing the frozen source, executable, accepted capture, or
locked A2.2 model.

The accepted replay reports `DeltaSUMEV_RAW-LOCAL=-0.25305509271796467 Ry`
against `t_up+t_down=-0.2546410140151576 Ry`, signed ratio `0.9937719329962356`,
and absolute error `0.001585921297192927 Ry`, within the preregistered
`A=0.001702238347231377 Ry`. The scalar control error is
`1.0447191816392865e-5 Ry`, within its `6.022982854018227e-5 Ry` allowance,
and `q2_local` is nonnegative for both spins (minimum `0.009544969868103789
Ry^2`). The measured NEWRHO l=2 output-density charge changes, sphere spin
moment, and all ten ASA components for every case are included in the case
report; the energy components use frozen input rho, with output-rho
components retained as a separate diagnostic.

The VOUT and transposed-P controls are both **nondiscriminating** under the
locked rule: their absolute carrier deviations are respectively
`0.0015859303878082986 Ry` and `0.001585921297192927 Ry`, each below A. They
are not counted as passing controls. Thus the required M0.3 controls do not
all pass, overall M0 (KS) is **not PASS**, and Stage M1 is not authorized.
No production energy correction, M1, or M2 follows from these measurements.
The preregistration, five-case output, and raw replay CSV are preserved in
[`evidence/a2_m03_preregistration.json`](evidence/a2_m03_preregistration.json),
[`evidence/a2_m03_cases.json`](evidence/a2_m03_cases.json), and
[`evidence/a2_m03_cases.csv`](evidence/a2_m03_cases.csv). The reports bind the
accepted capture, manifest, frozen executable, driver, A1 intercept, and A2.2
fit hashes.

Under A1.3/A2.6, the measured response allowance exceeds both wrong-control
deviations: A is `0.001702238347 Ry`, while each deviation is about
`0.00158593 Ry`. The independently measured operator separation is small at
this accepted state: `t(Vout)-t(Vin)` is `+1.06977e-8 Ry` for spin up and
`-1.41968e-9 Ry` for spin down, and `P_out` is symmetric at the recorded
precision, so transposing it leaves the first moment unchanged. The smooth
radial response error therefore masks the energy effect of these
converged-state controls; they are nondiscriminating by measurement, not
passing controls.

Proposed future-only criterion for user review: require an energy control to
be discriminating only when its preregistered operator separation `|δ_k|`
exceeds `A_true + A_k + E_num`, where `E_num` is the measured numerical replay
guard. If it does not, require the independent operator/convention oracle and
provenance check, and label the energy control nondiscriminating. This protects
whether the wrong pairing can be resolved above combined response and replay
error; its error model sums the locked true/control response allowances and
the measured replay guard. This proposal has not been applied to the current
outcome. The user must approve any future criterion change before new
measurements are inspected.

### Final A2 source verification

The final production snapshot used by the A2 capture is `source/self.f90`
SHA256 `d91a468a4493da21710fc7f60e745b5f1bef8f20d5d667ae939c31d5ba3061ef`.
On that exact source, a fresh Debug build passed 119/119 tests with
`ctest --test-dir /private/tmp/rslmto-localmom-a2-debug -LE
'validation|performance' --output-on-failure` in 335.66 seconds. The complete
CTest log is `/private/tmp/rslmto-localmom-a2-debug/Testing/Temporary/LastTest.log`
(SHA256 `41cc132aba07b7fe415b8cc95d30268bb3af1db61ddbd68888c7a190904e4752`).
The targeted MPI residual and common-/multisite frozen-H projector tests
passed 3/3 on the same source; their log SHA256 is
`9a01d251c77677384764aec77f7e05ea528e6802643ef33329d926e60ab9f07a`.
These tests validate the diagnostic plumbing and production-source regression
surface; they do not override the nondiscriminating M0.3 control verdict.
