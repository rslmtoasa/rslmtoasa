# Fixed-U local-moment accounting

This record keeps the conventions, registered criteria and measured key
numbers of the fixed-U local-moment (M0 to M2) accounting. Run provenance
(hashes, run directories, per-run narratives, verbatim amendment texts) is
not kept here; it is at tag `archive/fable_v4b_u-2026-11-05`, and
`VALIDATION_SUMMARY.md` gives the stage summary.

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
| M0.1 RS/KS M1 extraction | Complex matrix agrees with independent finite-H KS eigenpair sum at `5e-12 Ry` and refined analytic finite-window Lorentz RS integral at `3e-4 Ry`; same-call `Tr(M1)=q0*c` per atom/shell/spin at `1e-10 Ry` because both use identical quadrature | Conjugated/transposed and wrong-site M1 controls must mismatch by `>3e-3 Ry` (10x the accepted RS M1 oracle budget); coarse-to-fine mesh error must decrease | Measured: see M0.1 below (finite-H probes only) |
| M0.2 no-U calibration | Converged nsp=1 bcc Fe in RS and true KS; record spectral `sum(q0*c)`, converged `SUMEV`, and one frozen-input-potential `newrho` replay; report charge threshold and measured floor | Verify both restart U/J arrays are zero; deterministic restart and 16³ refinement; no-U must not activate R05 P/V stopping | Scalar-v0 criterion replaced by A1; measured values below |
| M0.3 U=2 eV, J=0 carrier | Reuse only converged accepted RS/KS states that meet the existing R05 P/V criteria; independent correction arithmetic; production `calculate_pl` and one frozen-sphere `newrho` pass for legacy and local moments; report signed ratio, `q2_local>=0`, radial d charge, sphere spin moment and ASA energy components | `V=vI` gives exactly invariant q2 and `delta SUMEV=v*q0`; substituting `V_out` or `P^T` must violate the predeclared carrier criterion | See A2 and M0.3 results below |

M0 needs admissible production evidence for all three rows. The previous
one-step `-0.153 Ry` comparison and synthetic/calibrated radial tests are
ineligible. A missing converged state is NOT_RUN. The original scalar-v0
criterion was replaced by A1, and A1 by A2 and A3 (below).

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

## Protocol criteria (amendments A1 to A5, condensed)

The user-approved amendments A1 to A5 changed the registered criteria in the
order below. Only the criteria and their stated decisions are kept here; the
verbatim amendment texts and approval exchanges are at tag
`archive/fable_v4b_u-2026-11-05`.

**A1 (scalar-shift scaling, KS, frozen no-U state).** Shifts
`v = +-0.0025, +-0.005, +-0.01, +-0.02 Ry` applied to `c` with `q0`, `q2`
held fixed. Report `r(v) = DeltaSUMEV(v) - v*q0` and the symmetric slope
`s(v) = [DeltaSUMEV(+v) - DeltaSUMEV(-v)]/(2v)`. Shape criterion: `r(v)`
smooth, log-log order of `|r|` versus `|v|` at least 2, no
`v`-independent plateau above `10*D0`. Where `r(v)` is flat or
non-monotonic, investigate the radial solver tolerance (TOLRSQ/RSEQSR)
without tuning it. The fitted `r(v)` model evaluated at the actual per-spin
shifts `t_s/q0_s` plus the U=0 floor is the M0.3 allowance. A1.2: RS is a
separate work item (full no-U residual history, continuation with Broyden
mixing, measured stagnation level used as the RS threshold for its own floor);
RS gates are reported separately from KS. A1.3: a criterion found to rest on
a wrong error model is not self-amended; a corrected criterion goes to the
user before further measurement.

**A2 (supersedes A1.1 criterion (i) and the floor-only M0.3 allowance).**
- A2.1: with `t_s = Tr(V_in,s P_out,s)` from the accepted KS U=2 eV, J=0
  state and the A1 spin intercepts `a0_up = -7.58154e-6`,
  `a0_down = -8.43595e-6` (not refit), require
  `B = |a0_up||t_up| + |a0_down||t_down| <= D0`, with
  `D0 = 5.028166848219939e-6 Ry` (the no-U KS floor).
- A2.2: add shifts `+-0.04, +-0.08, +-0.12 Ry` to the A1 eight; fourteen-point
  signed cubic `r_fit(v)`, all-point OLS of `log|r|` against `log|v|`. Shape
  checks: constant sign on each branch, `|r|` strictly decreasing as `|v|`
  halves, order at least 2, fitted constant at most `10*D0`. Each `t_s/q0_s`
  must lie within `[-0.12, +0.12] Ry`.
- M0.3 allowance: `A = |r_fit,up(x_up)| + |r_fit,down(x_down)| + D0`, with
  `x_s = t_s/q0_s`. Carrier check:
  `|DeltaSUMEV - (t_up+t_down)| <= A`, with
  `DeltaSUMEV = SUMEV_spectral - SUMEV_local`; the signed ratio
  `DeltaSUMEV/(t_up+t_down)` is reported. `q2_local >= 0` is required. Also
  reported: d radial charge, sphere spin moment, all ASA energy components.
  Controls: scalar shift (`q2` preserved, same allowance), `V_out`,
  transposed P.
- A2.4 (reported, not a gate): map `D(E) -> EVAL` near `E_nu`, report
  `f(E) = EVAL - E` and `f'(E_nu)` against the `phi-dot` derivative accuracy
  and the GFAC normalization mismatch (about 2e-4).
- A2.6: every preregistered tolerance states the decision it protects, the
  quantity it is compared against, and its error model.

**A3 (supersedes the A2 negative-control criterion).** On the converged cubic
collinear state `P^T = P` and `V_out = V_in` to the R05 tolerance, so those
controls cannot fail. A3.1: a control is admissible only if its predicted
separation (computed from the captured state before measuring) exceeds the
sum of the relevant carrier allowances plus the numerical guard; otherwise it
is recorded as symmetry- or construction-null with its measured null quantity
(`max|P-P^T|`, `max|Im P|`, `max|V_out-V_in|` per spin) and replaced by an
admissible control or an independent oracle. Substitute controls on the same
frozen state: (a) sign control `c_wrong = c + t/q0` (predicted ratio about
-1); (b) spin-provenance control using the other spin's `V_in` (predicted
separation `|Tr((V_up - V_dn) P_s)|` per spin); (c) iteration-provenance
control with `V` from iteration 1, recorded as null if its separation does
not exceed the allowance. The transpose convention is covered by the M0.1
noncubic oracles. A3.5: RS floor uncertainty from a deterministic repeat and
a 6000-channel refinement at the RS threshold `3e-8`; RS M0.3 under A2 and A3
with an RS-specific fit and `D0_RS`.

**A4 (Stage M2).** A4.0: RS local-mode identity, SCF and paired frozen
`P_OUT` replay before any RS M2 measurement; KS Mermin entropy diagnostic
`S = -sum_k w_k sum_n [f ln f + (1-f) ln(1-f)]` with `F = E - T*S` reported
beside each KS energy candidate.
- A4.1 external-potential control: test-only insertion of `lambda*W` in the d
  block through the `hubbard_u_pot` path, energy term `lambda*Tr(WP)`;
  `W1 = identity`, `W2` Hermitian anisotropic with complex off-diagonals,
  `max|eig| = 1`; U=0, `lambda = 0, +-0.0025, +-0.005, +-0.01 Ry` in both
  modes. Expected `dF/dlambda|0 = Tr(W P(0))` with
  `F_local = ETOT_local + lambda*Tr(WP) - T*S`,
  `F_spectral = ETOT_spectral - T*S`. Central differences, h^2 scaling,
  repeat runs at 10x tighter tolerances. Later h^2-dominance rule for the
  resolved leading term: `|L| - uL > 16(|H| + uH)`.
- A4.2 main test: U0 = 2 eV, J = 0, runs at `U0 +- delta`, `U0 +- 2 delta`,
  `delta = 0.2, 0.1, 0.05 eV`; target is `dE_Hub/dU` at fixed converged
  `P(U0)` from `evaluate_hubbard_energy`. J test at U0 = 2 eV,
  `J = 0.5, 0.9 eV`. Candidates:
  `a = F_spectral - t + E_Hub`, `b = F_local + E_Hub`,
  `c = F_spectral + E_Hub`. Prediction: c misses by about `dt/dU`
  (about 0.13 Ry/eV), b agrees within the error floor, a is reported only.
  Decision rule: b agrees and c misses, adopt local + b; b misses, stop and
  report; c agrees, stop (contradicts the derivation).
- A4.3 like-for-like observables, spectral versus local, at
  `(U,J) = (0,0), (2,0), (4,0), (2,0.9) eV`: sphere spin moment, d `q0` per
  spin, eigenvalue bounds of P, d-centre exchange splitting, E_F, E_int, E_DC,
  E_Hub, t, all candidates, F and T*S, residuals.
- A4.4 (after the KS M2 comparison): local as the default for fixed
  Liechtenstein `nsp=1` (KS; RS only with its own M1 and M2), candidate b in
  `physical_total_energy`, restart and MPI reproducibility.

**A5 (governance simplification and U/J test).** A5.0: one short
preregistration per experiment (prediction and magnitude, error budget with
order-of-magnitude justification, what each outcome implies); hashes kept for
source, executable, inputs and outputs only. A5.1: the W1 mismatch is a
property of the base ASA functional (see Measured results). A5.2 diagnostics
(hoh at U=0 with a W-only bypass, temperature halved and doubled, denser
k-mesh) are optional and have no production change. A5.3: error budget is the
base-functional non-stationarity propagated to `dE/dU`. With spin-resolved
occupations `q_s`, a shift `lambda_a` applied to spin `a` gives
`R_s,a = dq_s/dlambda_a` and `B_a = dF/dlambda_a - q_a`; solve `R^T g = B`
for the residual gradients `g_s` (Ry/electron) under a charge-only
assumption, and use the allowance `sum_s |g_s dq_s/dU|` plus propagated
uncertainty (Ry/eV). Decision rule: c must miss by far more than the budget,
b must agree within it, a is reported only; if b agrees and c misses, adopt
local + b with the label "valid to the ASA base-functional accuracy measured
in A5.1". A5.4: W2 state-level complex-convention test (transposed-P and
conjugated-P frozen dual replays), slope reported, not gated.

## Measured results

Raw data for everything below: tag `archive/fable_v4b_u-2026-11-05`.
Summary of the stages: `VALIDATION_SUMMARY.md`.

### M0.1 finite-operator probes

- Measured: common-H one-site and two-inequivalent-site probes, 2/2 pass in
  serial and 4/4 in an MPI Debug build.
- Measured: RS 6000-channel analytic Lorentz M1 errors `2.894374e-7 Ry` and
  `2.845306e-7 Ry`; KS eigenpair M1 errors `7.757919e-18 Ry` and
  `2.797158e-17 Ry`; maximum same-call `Tr(M1) - q0*c` residual
  `1.054712e-15 Ry`; the 750/1500/3000/6000 channel error sequences decrease
  strictly.
- Measured: M1 conjugation/transpose mismatches `0.069864` to `0.072836 Ry`;
  two-site site-swap mismatches `0.047631` to `0.047690 Ry` (minimum required
  `0.003 Ry`).
- Measured: trace-on and trace-off runs give the same `Fe_out` (serial and
  two-rank MPI; failed trace-file creation disables tracing without rank
  divergence). This covers the diagnostic path only, not the radial-energy
  hypothesis.
- Measured: serial `ctest -LE 'validation|performance'` 119/119 and 114/114
  unit tests on the source used for the A1 and A2 captures.

### M0.2 no-U KS calibration

- Measured: baseline `D = sum(q0*c) - SUMEV = -5.028166848219939e-6 Ry`,
  repeat `-5.023148958116508e-6 Ry`, 16^3 refinement
  `-5.053452480030529e-6 Ry`. Floor `D0 = 5.028166848219939e-6 Ry`,
  uncertainty `2.528563181058985e-8 Ry`, allowance
  `5.053452480030529e-6 Ry`.
- Measured: historical scalar shift `c_local = c_raw - 0.01 Ry`: expected
  `DeltaSUMEV = +0.06498480421874168 Ry`, production
  `+0.06497426910958881 Ry`, residual `1.0535109152873257e-5 Ry` (2.08 times
  the floor-only allowance of that time). The floor does not bound the
  shift-dependent `calculate_pl` inversion error, which is why A1 replaced
  that criterion.
- Measured (A1.1 scaling, KS, total q0 `6.498480421874168`): symmetric-slope
  relative error at `|v| = 0.0025, 0.005, 0.01, 0.02 Ry` is `-1.76935e-5`,
  `-4.64255e-5`, `-1.62889e-4`, `-6.27649e-4`; v to 0 intercept
  `-7.86514e-6` (spin up `-7.58154e-6`, spin down `-8.43595e-6`), carrier
  efficiency `1 + a0` about `0.999992`; all-point log-log order `2.72568`;
  fitted cubic constant `-5.11131e-8 Ry` (limit `10*D0 = 5.02817e-5 Ry`);
  residual at `v = +0.01 Ry` `-1.0535109153e-5 Ry`.
- Measured: RS no-U run reached its 120-step limit without meeting the `1e-9`
  charge criterion (final residual `1.57e-8`, last-ten range `2.9e-9` to
  `1.81e-8`).

### A2 results (KS, then RS)

- Measured (A2.1): accepted KS U=2 eV, J=0 state at iteration 32, charge
  residual `7.69e-8`; `t_up = -0.27138788113273177 Ry` (`q0_up =
  4.527615147715704`), `t_down = +0.016746867166353676 Ry` (`q0_down =
  2.08410054869144`); target shifts `-0.059940580698351495` and
  `+0.00803553704588229 Ry`; `B = 2.19881406993552e-6 Ry`, ratio to `D0`
  `0.437299`.
- Measured (A2.2): fourteen-point fit satisfies the shape checks for both
  spins and the aggregate; allowance `A = 0.0017022383465738291 Ry`; scalar
  `v = +0.01 Ry` allowance `6.022982854018227e-5 Ry`; expected scalar
  `DeltaSUMEV = 0.06611715696407144 Ry`. A first SCF recapture selected a
  different iteration (25 against 23) with operand differences up to
  `4.0e-9`; it was not used, and the accepted iteration-23 state was replayed
  from its stored operands (zero-shift and the eight A1 shifts reproduced the
  archived values bit for bit).
- Measured (RS floor, threshold `3e-8`): the 240-step continuation did not
  reach `1e-9`; windows 121-180 and 181-240 ranged `2.3e-9` to `2.53e-8` and
  `2.4e-9` to `2.67e-8`. `D = D0_RS = 5.0565002780e-6 Ry`
  (`sum Tr(M1) = -1.9975997502219833 Ry`, `SUMEV = -1.9975946937217053 Ry`).
- Measured (RS floor uncertainty): `D_base = -5.056500278044496e-6`,
  `D_repeat = -5.052056331500410e-6`, `D_refine6000 =
  -4.990233713142445e-6 Ry`; `u_RS = 6.626656490205107e-8 Ry`;
  `D0_RS + u_RS = 5.122766842946547e-6 Ry` (diagnostic only).

### M0.3 KS frozen-state replay (iteration 32)

- Measured: `DeltaSUMEV(RAW-LOCAL) = -0.25305509271796467 Ry` against
  `t_up + t_down = -0.2546410140151576 Ry`; signed ratio `0.9937719329962356`;
  error `0.001585921297192927 Ry` against `A = 0.001702238347231377 Ry`.
- Measured: scalar control error `1.0447191816392865e-5 Ry` against
  `6.022982854018227e-5 Ry`; minimum `q2_local` `0.009544969868103789 Ry^2`.
- Measured: V_out and transposed-P controls deviate by `0.0015859303878082986`
  and `0.001585921297192927 Ry`, both below `A`, so they do not discriminate
  at this state. Operator separation `t(Vout) - t(Vin)` is `+1.06977e-8 Ry`
  (up) and `-1.41968e-9 Ry` (down); `P_out` is symmetric at the recorded
  precision. This is why A3 replaced these controls.
- Measured (A3 nulls): `max|P-P^T| = (1.8015e-16, 6.2830e-17)`,
  `max|Im P| = (9.0854e-17, 3.1442e-17)`, `max|V_out - V_in| =
  (3.9630e-9, 1.4636e-9) Ry` (up, down). Numerical guard `E_num = 2e-8 Ry`.
- Measured (A3 controls, predicted before the responses): sign control ratio
  `-0.9941686941`, separation `0.5092820280 Ry` against bound
  `0.0031921762 Ry`; spin-provenance separations `0.3303599107` and
  `0.1373697148 Ry` against bounds `0.0017413743` and `0.0005938364 Ry`.
  Iteration-1 potentials are absent from the iteration-32 capture, so
  control (c) was not assessed for KS.

### Gate M1 (KS, local moments)

- Measured: KS U=2 eV, J=0 local-mode run converged at iteration 34, final
  SCF difference `2.04e-8`, final `RP = 2.8483418472811053e-8`,
  `RV = 4.1869823014390071e-9 Ry`; a deterministic repeat gave identical
  `Fe_out.nml` and identical Hubbard-energy and convergence log lines.
- Measured (paired frozen P_OUT replay): `DeltaSUMEV(spectral - local) =
  -0.2565860181649944 Ry`; per-spin ratios `0.993661677965` and
  `0.999792868074`; locked-model prediction `-0.2565645199104135 Ry`,
  difference `-2.14982546e-5 Ry` against allowance
  `0.00185641944033825 Ry` plus `2e-8 Ry`; `DeltaSUMEC = 0`; both `q2`
  nonnegative.
- Measured: sphere spin moment `2.51078366041631` (local) and
  `2.510781639399223` (spectral); l=2 radial density charge `4.591832540887848`
  and `1.929888414728295` (up, down, local).
- Measured (local endpoint): `E_int = 2.757238890068984 Ry`,
  `E_DC = 2.646759512588076 Ry`, `E_Hub = 0.1104793774809081 Ry`,
  `t_in_out = -0.2583779167349475 Ry`, physical total energy
  `-2541.972663505471 Ry`; frozen-response candidates
  `a = -2541.860392197098`, `b = -2541.862184095667`,
  `c = -2542.118770113832 Ry`.
- Measured: no-U, explicit U=J=0 and spectral U>0 outputs are bitwise
  identical to the baseline source in all eight registered artifacts; filtered
  serial suite 121/121 with 126 registered tests.

### RS M0.3 and RS Gate M1

- Measured (RS M0.3): accepted iteration 6, charge difference `1.66e-8 Ry`;
  per-spin ratios `0.994089411811` and `0.999901088591`; carrier deviation
  `0.00161403350146 Ry` against allowance `0.00173036300335 Ry`; minimum
  `q2_local` `0.00952358124507 Ry^2`; sign control separation
  `0.510786194545 Ry`; opposite-spin separation `0.197099739922 Ry`; scalar
  `+0.01 Ry` error `1.04448396420e-5 Ry` against `6.02422395038e-5 Ry`.
  Transpose, V_out and iteration controls are construction-null on this
  cubic state.
- Measured (RS Gate M1, local U=2 eV, J=0): converged at iteration 103,
  final charge difference `4.8e-9 Ry`. Paired frozen replay:
  `DeltaSUMEV = -0.2590196642435809 Ry` against `sum(t) =
  -0.26088126108452275 Ry`; error `0.001861596840941837 Ry` against
  `A_RS = 0.0019179425667351765 Ry`; minimum `q2_local`
  `0.01075725460260394 Ry^2`, minimum `q2_spectral`
  `0.009745244577861018 Ry^2`. The extra point-prediction test
  `sum(t) + sum(r_fit)` within `2.0002e-8 Ry` was not met (deviation
  `1.0234967739219947e-5 Ry`); under the approved supplemental criterion
  (carrier rule and `q2 >= 0` only) this extra test is diagnostic.

### A4.0 housekeeping

- Measured: finite-H Fermi-Dirac entropy oracle, `S = 2.1062540069077236` per
  cell, `T*S = 0.0040021025481680765 Ry/cell`; production reports `T*S` and
  `F = E - T*S` beside each KS candidate. RS entropy is reported as
  unavailable, not zero.
- The stored onsite d block uses spherical channels `m = -2..+2` at indices
  5:9 after the Cartesian-to-spherical transform; the source comment was
  corrected accordingly (basis analysis: tag `archive/fable_v4b_u-2026-11-05`).

### A4.1 W1 identity (U=0, KS)

- Measured: 24 runs (8 lambda-zero standard/tight captures plus the nonzero
  W1 set) all converged. Spectral and local modes both gave an unresolved
  h^2-dominance result under the A4.1 rule, with the `R_fine - g0` bias
  outside the propagated uncertainty. Under A5.1 this is a measurement of the
  base ASA functional: the U=0 Hellmann-Feynman identity for a d-level shift
  is off by `+0.29 %` (spectral) and `-0.14 %` (local); finite-difference
  error at most `5e-5` against a bias of at least `9e-3`. Aggregate W1
  residuals `+0.0185609` (spectral) and `-0.00921488` (local).

### A4.4 current stage (KS, 2026-10-18)

- Measured (U=2 eV, J=0): local b derivative `+0.0550976939 Ry/eV` against
  the direct-P production target `+0.0552396891 Ry/eV`; residual
  `-0.0001419952 Ry/eV`, inside the fine-grid local allowance
  `[0.00102115, 0.00155303] Ry/eV`.
- Measured (J): local b residuals `-0.0003734601` (J=0.5 eV) and
  `+0.0001251389 Ry/eV` (J=0.9 eV), inside the allowances
  `[0.00118574, 0.00174694]` and `[0.00117606, 0.00258153] Ry/eV`.
  Spectral c is outside its allowance at U and both J values.
- Candidate b = `F_local + E_Hub(P_OUT)` is the selected candidate, scoped to
  "valid to the ASA base-functional accuracy measured in A5.1". No a-versus-b
  preference is claimed.
- Eligible route: native reciprocal KS (`use_kspace=.true.`,
  `calctype='B'`), clean bulk (`nbulk=0`, `nrec=ntype`), explicit fixed
  general Liechtenstein U/J, scalar-collinear `nsp=1`, `lmax=2`, no
  constraints. SOC/noncollinear, HOH, self-consistent U, intersite V,
  impurity U, impurity/surface hosts and nonzero test fields are excluded.
  RS stays on its guarded spectral route unless its own M1 and M2 are met.
- Measured: serial restart continuity and matched serial/MPI runs reproduce;
  a constraints-enabled input selects spectral mode before the expected
  zero-target setup fatal. This checks the eligibility guard only, not a
  successful constrained SCF route. Raw data: tag
  `archive/fable_v4b_u-2026-11-05`.
