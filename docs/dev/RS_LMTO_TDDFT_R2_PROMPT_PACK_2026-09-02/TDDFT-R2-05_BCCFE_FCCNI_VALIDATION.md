# TDDFT-R2-05 — Production bcc Fe + fcc Ni physics validation

Work on current `fable_v4`.

This is **TDDFT-R2-05: perform the first genuine production-level material validation of transverse TD-DFT using the existing `examples/susceptibility` bcc Fe and fcc Ni fixtures**.

This task is not another smoke test.

Its purpose is to decide whether the current TD-DFT implementation is scientifically trustworthy for the baseline case

\[
\text{collinear FM},\qquad SOC=0,\qquad \text{transverse response}.
\]

Do not modify the physics implementation merely because a result disagrees with expectations. Treat disagreement as evidence to diagnose.

## 1. Use the existing material fixtures

Start from the existing:

- `examples/susceptibility/...bccFe...`
- `examples/susceptibility/...fccNi...`

or their actual current names.

Preserve these as recognizable examples rather than creating unrelated synthetic decks.

It is acceptable to add dedicated high-accuracy validation variants beside them.

## 2. Establish a convergence hierarchy

For each material converge at least

\[
N_k,\qquad
\eta,\qquad
N_\omega/\Delta\omega,\qquad
\text{electronic energy/band window}.
\]

For native RGF also converge

\[
R_{\rm source},\qquad
R_{\max},\qquad
N_{\rm contour}
\]

where applicable.

Do not vary everything simultaneously.

Establish a documented reference calculation and then one-axis convergence sweeps.

## 3. Three-backend \(\chi^0\) equivalence

For converged representative points explicitly compare

\[
\chi^0_{\rm eig},
\qquad
\chi^0_{kGF},
\qquad
\chi^0_{RGF}.
\]

Compare the complex response matrices themselves before solving Dyson.

Required points:

- \(q=0,\omega=0\);
- a small nonzero \(q\);
- a representative magnon-region \(\omega\);
- a point in or near the Stoner continuum.

Report norms such as

\[
\frac{\|\chi^0_A-\chi^0_B\|}
{\max(\|\chi^0_A\|,\|\chi^0_B\|)}.
\]

## 4. Raw Goldstone/Ward validation

With SOC disabled, evaluate

\[
r =
[1-\chi^0(0,0)K_{\rm xc}]m_{\rm GS}.
\]

Record:

- raw Ward residual;
- relevant eigenvalue of \(\chi^0K_{\rm xc}\);
- overlap of its mode with the ground-state magnetization;
- raw acoustic gap.

Do this **without Goldstone correction first**.

The raw result is the scientific diagnostic.

If a correction mode is also demonstrated, report separately:

- raw result;
- corrected result;
- magnitude of the correction.

Never substitute the corrected value for the raw release evidence.

## 5. Small-\(q\) dispersion

Use a dedicated dense small-\(q\) sequence near \(\Gamma\), not only a conventional high-symmetry path.

Fit

\[
\omega(q)
=
\Delta + Dq^2 + Cq^4
\]

and also test a restricted quadratic fit

\[
\omega(q)\approx \Delta+Dq^2.
\]

Report:

- \(\Delta\);
- \(D\);
- uncertainty/sensitivity to fitting window;
- residuals;
- whether adding \(q^4\) is statistically/numerically necessary.

For SOC=0, the physical target is

\[
\Delta\rightarrow0
\]

under convergence.

## 6. Special emphasis on fcc Ni

The previous TD-DFT implementation showed an insufficiently convincing quadratic dispersion around \(\Gamma\).

Treat this as a targeted regression problem.

Establish whether the old behavior came from:

- \(q\)-path/backfolding conventions;
- insufficient \(k\)-sampling;
- broadening;
- frequency-grid resolution;
- vertex/kernel inconsistency;
- mode-tracking failure;
- some other cause.

Do not declare Ni fixed merely because the spectral plot looks smoother.

Demonstrate quantitatively that

\[
\omega(q)\propto q^2
\]

over a defensible low-\(q\) interval.

## 7. Independent adiabatic comparison

Where the existing RS-LMTO workflows allow it, compare the TD-DFT low-energy dispersion/stiffness against an independent magnetic reference:

\[
\text{LKAG}
\]

and, when the validated GBT path is practical,

\[
\text{GBT/frozen magnon}.
\]

The comparison need not be exact at higher energy because TD-DFT contains Stoner hybridization and Landau damping.

But near \(\Gamma\) the stiffness should be mutually consistent within understood approximations and numerical convergence.

## 8. Stoner continuum and damping

For several \(q\) values record:

- \(-\mathrm{Im}\,\chi^0\);
- interacting loss spectrum;
- magnon peak position;
- FWHM where a meaningful quasiparticle peak exists.

Verify that increasing overlap with the Stoner continuum produces physically sensible damping rather than numerical peak broadening.

Explicitly separate physical FWHM from artificial \(\eta\)-broadening by convergence analysis.

## 9. Required material-validation artifacts

Store a compact reproducible evidence bundle in the repository or in a deliberately referenced validation-data location containing:

- exact input decks;
- git commit;
- compiler/build information;
- MPI/OpenMP settings;
- all convergence tables;
- raw three-backend comparison data;
- raw Goldstone diagnostics;
- small-\(q\) fit data;
- spectra used for damping extraction;
- independent LKAG/GBT comparison where available.

Generate a concise Markdown report:

`docs/validation/TDDFT_BCCFE_FCCNI_VALIDATION.md`

or an equivalent project-consistent location.

The report must contain explicit

**PASS / PARTIAL / FAIL**

for Fe and Ni separately.

## Acceptance checklist

The completed evidence report is
[`docs/validation/TDDFT_BCCFE_FCCNI_VALIDATION.md`](../../validation/TDDFT_BCCFE_FCCNI_VALIDATION.md).
A checked item means that the evidence or investigation was performed; an
unchecked item remains an unmet acceptance gate.

### bcc Fe

- [x] high-accuracy fixture established.
- [x] \(k\)-mesh convergence measured; non-monotonic raw Ward behavior prevents a convergence claim.
- [x] \(\eta\) convergence measured; widths remain numerical-broadening diagnostics.
- [ ] frequency-grid convergence demonstrated.
- [ ] three-backend static \(\chi^0\) agreement; current Γ,0 matrix norms fail the `1e-8` gate.
- [ ] three-backend dynamic \(\chi^0\) agreement; current selected-point norms fail the gate.
- [x] raw Ward residual reported.
- [ ] raw \(\Gamma\)-point acoustic gap established from a converged interacting branch.
- [ ] quadratic small-\(q\) dispersion demonstrated; measured zero-intercept fit is rejected.
- [x] diagnostic spin-wave stiffness extracted.
- [x] independent adiabatic comparison performed where available; the Jij record is incomplete.
- [x] damping/Stoner behavior examined without claiming intrinsic damping.

### fcc Ni

- [x] high-accuracy diagnostic fixture established; the retained reference is not production-valid itinerant Ni.
- [x] \(k\)-mesh convergence measured; stable diagnostics do not validate the reference state.
- [x] \(\eta\) convergence measured.
- [ ] frequency-grid convergence demonstrated.
- [ ] three-backend static \(\chi^0\) agreement; native Γ,0 disagrees maximally with the reciprocal routes.
- [ ] three-backend dynamic \(\chi^0\) agreement.
- [x] raw Ward residual reported.
- [ ] raw \(\Gamma\)-point acoustic gap established from a converged interacting branch.
- [x] old non-\(q^2\) issue explicitly investigated.
- [x] quadratic small-\(q\) failure clearly documented.
- [ ] stiffness extracted.
- [ ] independent adiabatic comparison performed where available.
- [x] damping/Stoner behavior examined without claiming intrinsic damping.

### Release decision

- [x] Fe assigned **PARTIAL**.
- [x] Ni assigned **FAIL**.
- [x] no synthetic fixture used as substitute for the material evidence.
- [x] no Goldstone correction used to hide a poor raw Ward result.

## Commit

`Validate TDDFT for bcc Fe and fcc Ni`
