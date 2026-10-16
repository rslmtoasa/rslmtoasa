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

Proposed future-only criterion for user review: before future control outcomes
are inspected, compute each wrong-control separation `δ_k=Σ_s(t_{k,s}−t_s)`
from the frozen matrices. Let `A_true` be the locked allowance at the correct
shifts, `A_k` the same locked fit's allowance at the wrong-control shifts, and
`E_num` the predeclared replay roundoff guard. Require an energy-discrimination
verdict only when `|δ_k| > A_true + A_k + E_num`; an eligible control must then
violate `A_true` in production replay. If neither control is eligible, require
independent operator-convention and provenance oracles, and keep the gate
incomplete if those oracles fail. This protects whether the wrong pairing can
be resolved above combined model and numerical error; its error model sums the
locked true/control fit-plus-floor allowances and measured replay guard. Any
energy control must use an accepted physical state and production
`calculate_pl`/`newrho`. This proposal has not been applied to the current
outcome. The user must approve any future criterion change before new
measurements are inspected.

## Protocol amendment A3 (verbatim)

```text
Protocol amendment A3 (user-approved; supersedes the M0.3 negative-control
criterion; record verbatim in LOCAL_MOMENTS.md; do not edit the user-owned
protocol file or historical evidence).

Rationale: on the converged collinear cubic bcc-Fe KS state, both mandatory
controls are null by construction:
  - P^T = P by time-reversal and cubic symmetry, so Tr(V P^T) = Tr(V P);
  - V_out = V_in to the R05 tolerance at convergence.
A control that cannot fail on the tested state protects no decision. Codex's
proposed rule is adopted, with the substitutions below.

A3.1 Control applicability rule (all future gates):
  - Before measuring, compute each wrong control's predicted separation from the
    captured state.
  - A control is admissible only if that predicted separation exceeds the sum of
    the relevant carrier allowances plus the numerical guard.
  - An inadmissible control is recorded as "symmetry-/construction-null" together
    with its measured null quantity, and is replaced by an admissible control that
    tests the same convention. Where no state-level control exists, cite an
    independent oracle that tests that convention.

A3.2 Record the null quantities on the M0.3 state:
  - max|P - P^T| and max|Im P| per spin;
  - max|V_out - V_in| per spin.

A3.3 Substitute admissible controls on the same frozen M0.3 state. Run each
  through the same production calculate_pl/newrho pass and predeclare its
  predicted carrier statistic before running:
  (a) Sign control: c_wrong = c + t/q0. Predicted ratio ≈ -1.
  (b) Spin-provenance control: t_wrong,s = Tr(V_in,s' P_out,s), using the
      other spin's potential. Predicted separation |Tr((V_up - V_dn) P_s)|
      per spin; it must exceed the allowance.
  (c) Iteration-provenance control: V from a non-converged earlier SCF
      iteration of the same run (e.g. iteration 1) in place of V_in, with the
      predicted separation from the captured arrays. If this does not exceed
      the allowance, record it as null under A3.1 and rely on (b) together
      with the M1.4 finite-H V_out/V_in unit test.

  The transpose convention is covered by the M0.1 class-B noncubic oracles
  (conjugation/transpose mismatches 0.070-0.073 Ry, already PASS) and will
  be re-tested in the M1.4 unit tests. Cite these; do not re-run them here.

A3.4 KS M0 is PASS if:
  - the A2.1 and A2.2 results stand;
  - the M0.3 carrier, scalar and q2 checks stand;
  - the A3.2 nulls are recorded;
  - the admissible controls (a) and (b), and (c) if admissible, fail as predicted.

  Then proceed to Stage M1 for KS and stop at Gate M1.

A3.5 RS (separate, non-blocking):
  - Measure the D0_RS uncertainty with the same deterministic repeat plus
    refinement procedure used for KS, at the 3e-8 threshold.
  - Then run RS M0.3 on the converged RS U = 2 eV state under A2 and A3.
  - Report RS M0 separately.
```

### A3 control preregistration (before radial responses)

The pre-A3 evidence inventory is bound in
[`evidence/a3_pre_a3_evidence_snapshot.json`](evidence/a3_pre_a3_evidence_snapshot.json).
The protocol hash and every listed historical evidence hash were recorded
before A3 work; those files and the user-owned protocol remain unchanged.
The verbatim A3 amendment above was recorded before this preregistration.

For the accepted KS U=2 eV, J=0 capture (iteration 32; CSV SHA256
`41e8bec19e50ba9931ccab2077686647300e1142774eea45567d1f65602472bc`),
the locked A2 model remains unchanged (SHA256
`0ad6a7703aa26eb132602e7a43773d7619f70651c6b90817e4f945ff259cb993`).
Before any new radial response, the controls were evaluated from the frozen
matrices. The pre-response machine-readable record is
[`evidence/a3_m03_control_preregistration.json`](evidence/a3_m03_control_preregistration.json)
(SHA256 `1f26829dfebffd2586ac891100f485bb3eda3ec7f6655fd100d56c6267dc177a`).
It locks `E_num=2e-8 Ry`: the `RAW-control` difference propagates two
per-case `ENERGY_TOL=1e-8 Ry` guards; the existing matrix/control arithmetic
cross-check is at most `2e-12 Ry` and negligible at this scale. It also
conservatively assigns the full `D0=5.028166848219939e-6 Ry` to each side of
each spinwise applicability comparison. The true allowance is
`A_true=0.001702238347231377 Ry`.

The recorded per-spin nulls are `max|P-P^T|=(1.8015e-16,6.2830e-17)`,
`max|Im P|=(9.0854e-17,3.1442e-17)`, and
`max|V_out-V_in|=(3.9630e-9,1.4636e-9) Ry` (up, down). The trace differences
for the last quantity are also recorded in the JSON. For the sign control,
only the center changes (`c_wrong=c+t/q0`) and the correct `q2_local` is held
fixed; its predicted ratio is `-0.9941686941`, with separation
`0.5092820280 Ry` against an applicability bound of `0.0031921762 Ry`. For
the spin-provenance control, `V_in` is taken from the other spin and
`q2_local` is recomputed. Predicted separations are `0.3303599107` and
`0.1373697148 Ry` against conservative spinwise bounds `0.0017413743` and
`0.0005938364 Ry`; both are admissible. The transpose convention cites the
already-passing M0.1 Class-B noncubic oracle (`0.070–0.073 Ry` mismatch) and
is scheduled again in M1.4.

Iteration-1 potential provenance is absent from the accepted iteration-32
capture. Its applicability is MISSING/UNASSESSED, so no KS M0 PASS can be
claimed until an exact d91 replay binds iteration 1 to the accepted endpoint
and complete charge/R05 trajectory, or the user approves a scope resolution.
The sign and spin predictions are preregistered and released for production
response measurement; iteration control remains pending the exact d91 replay.
RS remains a separate, non-blocking A3.5 task.

The RS repeat/refinement protocol was preregistered separately before either
outcome in [`evidence/a3_rs_floor_preregistration.json`](evidence/a3_rs_floor_preregistration.json)
(SHA256 `ee9d52104fc3fbd8071a68d4f84dbefeb73e053c9cc76f9a2cc740ac5276eef6`).
It binds the accepted no-U RS floor's final restart (SHA256
`b6f3ce8058139380bd6fe44d74cf2889221da92ea8f30307944d3bfe515f998b`), the
same base executable/source, repeats the 2,500-channel run unchanged, and
refines only `channels_ldos` to 6,000. As in the KS procedure, both runs start
from the accepted base final restart, not its pre-floor continuation seed. The
6,000 endpoint is the existing M0.1 RS quadrature refinement endpoint. The
RS uncertainty is the maximum absolute signed-D deviation of repeat/refined
runs from the baseline and is reported separately. It is not silently added
to the A2 allowance. RS requires its own no-U scalar fit and allowance
`sum_s |r_fit,s(x_s)|+D0_RS`; no KS fit/allowance is imported. Runs use only
the RS-specific `3e-8` threshold.

The RS scalar-model preregistration for A2.2 is separately hash-bound in
[`evidence/a3_rs_a22_scalar_preregistration.json`](evidence/a3_rs_a22_scalar_preregistration.json)
(SHA256 `18e704b1fe9efe4ffedd9da3c54dd90d9e4f15d6a47eeacfce893b9cd135fd6d`).
It binds the accepted RS no-U iteration-2 CSV and the existing base plus eight
A1 signed responses, and preregisters only the six additional shifts
`±0.04, ±0.08, ±0.12 Ry` before those responses. It locks the same A2.2
14-point signed-cubic model and shape checks, but uses the RS floor and
RS-specific allowance `Σ|r_fit,s(t_s/q0_s)|+D0_RS`. The A2.5 repeat/refinement
uncertainty remains a separate diagnostic.

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

### A3 KS Gate M1

The independently reviewed A3.4 KS M0 gate passed, authorizing Stage M1 for KS
only. Work stopped at Gate M1; no M2 work was started. The M1 production source
is the frozen six-file snapshot recorded in the build/run manifests. The
matched Release build uses GNU Fortran 16.2 with Accelerate, MPI and OpenMP
disabled; its executable SHA256 is
`71ad3c0861996ca714da5d4d2e03af2611995eb3cf627f856ed8c649e8136184`.

The M1 local-moment KS U=2 eV, J=0 run converged at iteration 34 with final
SCF difference `2.04e-8`; the accepted fixed-U residual streak was three, with
final `RP=2.8483418472811053e-8` and `RV=4.1869823014390071e-9 Ry`. The
capture-only deterministic repeat used the identical input PRE deck, Fe seed,
executable, and physical settings. Both endpoint `Fe_out.nml` files have
SHA256 `23137d3165f92980da0a5dbe88940463850b2c8708e3e9acc4bce03dadbc57d7`;
all 102 normalized Hubbard-energy and convergence log lines match exactly
(SHA256 `55ee48877b44e47ddeac9af0d9a49dbbede5a2223b69f0675bdb233399daff66`).
The capture manifest is `/private/tmp/a3-ks-m1-local-u2-capture/manifest.json`
(SHA256 `474f142e976616fffc0ae2d7610e34a30143489e41e5f9bd94026530214e3502`);
its records CSV SHA256 is
`f10790d5b007d41cfebfe86df1f1a305140769e0c43053f8a7eaf932e1de1c7c`.

Before the paired frozen radial pass, the local and spectral P_OUT operands
were preregistered in
[`evidence/a3_m13_pout_local_spectral_dual_replay_preregistration.json`](evidence/a3_m13_pout_local_spectral_dual_replay_preregistration.json)
(SHA256 `4b8e182cfabe175f030103581082a41e20c8f330a5035a6ead0ca3b1f4436325`).
An earlier draft was superseded before its proposed radial responses because
it mixed P_OUT moments with a LIVE_STATE radial q2 operand; its immutable
supersession record is
[`evidence/a3_m13_dual_replay_preregistration_supersession.json`](evidence/a3_m13_dual_replay_preregistration_supersession.json).
Both cases use the same P_OUT d-shell occupation and the same captured radial
mesh, potential, density, core/valence seeds, other-shell moments, and PL
inputs. Only the d-shell gravity center and q2 differ: local P_OUT moments for
the local case, and independently reconstructed spectral moments for the
spectral case. The direct eigenpair reconstruction was independently checked
against P_OUT/M1_OUT and the locked inverse identity in the separately locked
[`evidence/a3_m13_spectral_reconstruction_preregistration.json`](evidence/a3_m13_spectral_reconstruction_preregistration.json)
(SHA256 `ee1184d6776ad239b6cb1b3dbe902ba93e371a5620511c9471eabda75a33b8dd`).
It does not use `m0_carrier_*` diagnostics or the misnormalized displayed
second-moment field. The earlier mixed LIVE_STATE replay is reported only as a
timing diagnostic and is not combined with the unmixed P_OUT moments.

The frozen dual response measured
`ΔΣSUMEV(spectral−local)=-0.2565860181649944 Ry`, with signed per-spin ratios
`0.993661677965` and `0.999792868074` against `t_s=Tr(V_in,s P_out,s)`. The
locked model predicted `-0.2565645199104135 Ry`; the `-2.14982546e-5 Ry`
difference is within the unchanged A2 allowance `0.00185641944033825 Ry` plus
the two-case `2e-8 Ry` numerical guard. `ΔΣSUMEC=0` exactly. Both local and
spectral q2 values are nonnegative. The P_OUT local candidate reproduces its
captured CALCULATE_PL output exactly; the separately captured mixed LIVE_STATE
anchor reproduces PL, SUMEC, SUMEV and all ten frozen-input-density components
exactly, with maximum rho difference `3.55e-15`.

For the P_OUT local case, the sphere spin moment is `2.51078366041631`; the
actual l=2 radial density changes integrate to `4.591832540887848` and
`1.929888414728295` for up and down. For the spectral case, these are
`2.510781639399223`, `4.591830496386397`, and `1.929888391243929`. The full
ten-component frozen-input-density and output-density diagnostic tables,
per-spin energies, charge integrals, response hashes, and mixed-anchor timing
differences are in
[`evidence/a3_m13_ks_m1_dual_replay_results.json`](evidence/a3_m13_ks_m1_dual_replay_results.json)
(SHA256 recorded with the evidence artifact). Component labels distinguish
the frozen input density used by POISS0/VXC0SP from output-density diagnostics.

On the same converged local endpoint, production reports `E_int=2.757238890068984 Ry`,
`E_DC=2.646759512588076 Ry`, `E_Hub=0.1104793774809081 Ry`,
`t_in_out=-0.2583779167349475 Ry`, and physical total energy
`-2541.972663505471 Ry`. Candidate b (`-2541.862184127990 Ry`) is emitted for
the local mode. Candidates a and c remain `UNAVAILABLE` in that production
run; the source marks all candidates `PROVISIONAL_ONLY`. The frozen radial
comparison is diagnostic and does not relabel the self-consistent physical
total energy. Applying the candidate formulas to the paired frozen P_OUT
responses gives `ETOT_local=-2541.972663473148 Ry` and
`ETOT_spectral=-2542.229249491313 Ry`, then
`a=ETOT_spectral-t+E_Hub=-2541.860392197098 Ry`,
`b=ETOT_local+E_Hub=-2541.862184095667 Ry`, and
`c=ETOT_spectral+E_Hub=-2542.118770113832 Ry`. The site and total values
coincide for this one-site Fe case. These are `PROVISIONAL_ONLY` frozen
response candidates, not self-consistent production-mode totals. Candidate b
differs from the separately logged live-mode b by `3.23e-8 Ry`, the recorded
P_OUT-versus-LIVE mixed-snapshot timing difference. The formulas, endpoints,
and hashes are in
[`evidence/a3_m13_pout_energy_candidates.json`](evidence/a3_m13_pout_energy_candidates.json)
(SHA256 `c0b561dcf99ebcf2481b93dd640dd2b3fb95a7b0b04b9b1c142cea49e26c0323`).
No M1 scalar `+0.01 Ry` response was preregistered or run; the scalar response
gate passed earlier under A3 M0.3. No M2 physical-energy determination is
claimed.

The matched baseline comparison against source commit `3d7571a` passed for
no-U, explicit U=J=0, and spectral U>0. All eight registered numerical output
artifacts matched bitwise in every pair; the M1 default absent keyword also
matched explicit `spectral`. Source review confirmed the remaining differences
were manifests/logs and elapsed `cpu_s` only. The comparison report is
`/private/tmp/a3-m14-baseline-comparisons/comparison.final.json`, SHA256
`b7d60bf6c1c9c4f47e1b1053c99d430d33878885051f799f0b75dc84be295d7a`.
The first serial CTest run was superseded after audit found its registry lacked
the newly added site-swap target. The refreshed Release registry contains 126
tests, including both local-moment targets; the final filtered suite passed
121/121 with validation/performance labels excluded. Its immutable log and
source/build hashes are recorded in
[`evidence/a3_m14_final_serial_ctest.json`](evidence/a3_m14_final_serial_ctest.json)
(SHA256 `1a81bdd0781e30fb9f96ca769a59e7cf76d1de75db6624e0d5b1d1fd810e70ef`);
the copied CTest log SHA256 is
`3eded0b86179284e77ff23bd455764d16a035ae8992a0e38d86d4fed536c7238`. Targeted
MPI checks passed 3/3 on the same frozen production source. These results
complete KS Gate M1 only.

### A3.5 RS floor uncertainty

The RS deterministic repeat and 6,000-channel refinement were measured from
the accepted base final restart, with the repeat retaining 2,500 channels and
the refinement changing only `channels_ldos`. Both runs converged under the
RS-only `3e-8 Ry` threshold. The signed floor values are
`D_base=-5.056500278044496e-6 Ry`, `D_repeat=-5.052056331500410e-6 Ry`, and
`D_refine6000=-4.990233713142445e-6 Ry`. Thus
`u_RS=max(|D_repeat-D_base|,|D_refine-D_base|)=6.626656490205107e-8 Ry`.
`D0_RS+u_RS=5.122766842946547e-6 Ry` is reported as a diagnostic; `u_RS` is
not added to the A2 primary allowance, which remains the RS-specific
`Σ|r_fit,s(t_s/q0_s)|+D0_RS` with `D0_RS` counted once. Full PRE/run manifests,
input and endpoint hashes, convergence markers, and signed-D inputs are in
[`evidence/a3_rs_floor_uncertainty_results.json`](evidence/a3_rs_floor_uncertainty_results.json)
(SHA256 `b705d83cb206028115475afdb19a68939c19d503c7a6a6475126ea4881ef5056`).
The uncertainty result is separate from the RS U=2 M0.3 assessment; it does
not change the KS Gate M1 PASS or authorize M2.

### A3.5 RS M0.3 result

The RS-specific A2.2 scalar fit used the accepted RS no-U state, its original
eight A1 shifts, and the six preregistered shifts. The source-aware fit report
and raw measurements are preserved in
[`evidence/a3_rs_a22_scalar_fit_results.json`](evidence/a3_rs_a22_scalar_fit_results.json)
(SHA256 `0f0b543f81f53bc1c7fdfc46d3b404adef1ca455b8366a90e87154afd2b0987c`),
[`evidence/a3_rs_a22_six_shift_responses.csv`](evidence/a3_rs_a22_six_shift_responses.csv)
(SHA256 `293db40ac192edb5f4428e7a87fff85ad99061c6c51d602c1291a67aac46cdb0`),
and [`evidence/a3_rs_a22_raw_anchor.csv`](evidence/a3_rs_a22_raw_anchor.csv)
(SHA256 `60f5a8ff67804adcdbc4fcc7c8c4a4ef3d24bce6cf826762bb28116e9f09560f`).
The six-point extensions passed the locked signed-cubic shape checks for both
spins. This RS fit supplies the RS-specific residual allowances; no KS fit is
used.

The converged RS U=2 eV, J=0 capture used the matched d91 executable and the
RS-only `3e-8 Ry` threshold. Iteration 6 was accepted at charge difference
`1.66e-8 Ry` with a three-step physical R05 streak. The hash-bound capture
manifest and CSV are `54d0c8d8b510ac73254c3cf10f6fd7f9bdbcb16b125dec1e18ddbfa621443054`
and `9982b7fc115e72d9f5236253717c8a6edcb8e9d2f95ade6a9766b7fd23189d4e`.
The capture-only same-run iteration-1 replay has matching endpoint and
trajectory identity (`e65908b9756558ab343ea3722639b0d2738a6b9d74e6c2a0981c5bfa6bea74ec`,
`b5f4fa1691e391d3f588137c809f39fa8e1e2fdcd19b8c7e93879b0b24fc776e`). Its
potential separation is below the locked applicability threshold, so the
iteration control is recorded as construction-null and no radial response was
run for it.

The first frozen LIVE anchor reproduced PL, SUMEC, SUMEV, EV, and all ten
frozen-input-density energy components exactly; maximum density difference
was `3.55e-15`. The measured carrier has signed per-spin ratios
`0.994089411811` and `0.999901088591`; its total absolute deviation from
`ΣTr(V_in P_out)` is `0.00161403350146 Ry`, within the RS allowance
`0.00173036300335 Ry`. The minimum corrected `q2_local` is
`0.00952358124507 Ry²`. The sign control has a `0.510786194545 Ry`
separation from the true carrier; the opposite-spin control has a
`0.197099739922 Ry` separation. Both responses agree with their
residual-inclusive predictions within the RS-specific wrong-control
allowances. The `+0.01 Ry` scalar response differs from `v q0` by
`1.04448396420e-5 Ry`, within its `6.02422395038e-5 Ry` allowance. The
measured `P-Pᵀ`, imaginary-P, and `V_out-V_in` null quantities are recorded
per spin in the results report; transpose, Vout, and iteration controls are
construction-null on this converged cubic state. A3’s transpose convention
also remains covered by the already-passed M0.1 class-B noncubic oracles.

The exact preregistration, prediction lock, anchor report and case report are
[`evidence/a3_rs_m03_preregistration.json`](evidence/a3_rs_m03_preregistration.json),
[`evidence/a3_rs_m03_prediction_lock.json`](evidence/a3_rs_m03_prediction_lock.json),
[`evidence/a3_rs_m03_live_anchor_results.json`](evidence/a3_rs_m03_live_anchor_results.json)
(SHA256 `f8f12836d9c33073d84b57512fb22d4a472c00fd3da5a62328721d7ccfa31e10`),
and [`evidence/a3_rs_m03_response_results.json`](evidence/a3_rs_m03_response_results.json)
(SHA256 `b7c7f2100f9c036949935d722c26ea0e4458fe481ef9261167b83b49ee883cb6`).
The raw anchor and response CSV hashes are `1f227d9e41fbd1cefd0c0a4836bce83951a291050ccf852277b7ae63e1aa4e46` and
`af6aa568b4a548430b359e1c633d24792e1a1d5365052576bdc64fc9f91ce7d7`.
Independent review passed RS M0 under A2/A3. This closes the separate RS M0
reporting item; no RS M1 or M2 work was started. All candidate energy values
remain diagnostic/provisional and do not establish a physical-energy
correction.

An earlier set of RS A3 runner/lock copies remains in the original user
workspace at `/Users/andersb/Jobb/rslmto_fable_v4b`:
`tests/validation/fixed_ldau/rs_a3_m03_runner.py`,
`docs/validation/fixed_ldau/evidence/a3_rs_m03_prediction_lock.json`, and
`docs/validation/fixed_ldau/evidence/a3_rs_m03_preregistration.json`. These
copies differ from the governing worktree versions and were left untouched.
Automatic approval review rejected their removal as destructive because they
are frozen runner/preregistration evidence; the active worktree versions and
hashes above govern this record.
