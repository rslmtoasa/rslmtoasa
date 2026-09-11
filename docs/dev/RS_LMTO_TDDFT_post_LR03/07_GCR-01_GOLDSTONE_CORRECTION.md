# GCR-01 — Implement the published BES Goldstone eigenvalue correction

## Prerequisites

- LR-06 PASS.
- KXC-01 PASS.
- LR-04 metric/operator algebra available.
- LR-03 rigid-rotation vector convention available.

## Goal

Implement the Buczek-Ernst-Sandratskii numerical Goldstone correction as a
strictly optional, explicitly selected route.

Do not enable it by default.

Do not use it to conceal an unvalidated bare susceptibility or a grossly wrong
kernel.

## Object to correct

Construct the exact discrete form of

\[
D=I-\chi_{KS}(0)K_{xc}
\]

using the LR-04 operator algebra.

Do not multiply raw pointwise matrices with ordinary `matmul` unless LR-04 says
that representation is already metric free.

## Eigenproblem

Determine from the actual static representation whether `D` is:

- Hermitian;
- metric-Hermitian;
- complex non-Hermitian.

Use the corresponding mathematically correct eigensolver.

Do not call a Hermitian solver merely because the physical zero mode is expected
to be real.

If non-Hermitian, reconstruct the correction using the appropriate right/left or
inverse-eigenvector machinery required by the published method.

## Goldstone-mode identification

A small eigenvalue alone is insufficient.

Require:

1. proximity to zero relative to the rest of the spectrum;
2. strong overlap with the LR-03 rigid-rotation vector using the LR-04 metric;
3. consistency across modest eta/k-mesh changes;
4. no physical SOC/external/constraining field that should gap the mode.

If assignment is ambiguous, REFUSE correction.

## Correction

Follow the published relation exactly:

- change only the identified Goldstone eigenvalue;
- preserve all other eigenvalues/eigenvectors as prescribed;
- reconstruct the corrected denominator/kernel according to BES;
- retain both raw and corrected objects for diagnostics.

Do not apply a global scalar rescaling.

## Refusal criteria

Do not correct when:

- multiple unrelated small modes are present;
- rigid-rotation overlap is weak;
- raw Goldstone error is comparable with non-Goldstone spectral defects;
- LR-06 convergence is poor;
- KXC units/provenance mismatch;
- finite external field/SOC makes zero mode unphysical.

## Tests

- exact Goldstone manufactured system;
- controlled perturbation of its single zero mode;
- wrong-small-eigenvalue fixture with weak rigid overlap;
- two-small-eigenvalue ambiguous fixture;
- verify all untouched eigenvalues remain unchanged;
- verify metric-aware rigid mode is restored;
- verify correction remains optional.

## Acceptance checklist

- [x] BES `D=I-chiKS(0)Kxc` is formed with LR-04 canonical operator composition;
- [x] correction is disabled by default and requires explicit selection;
- [x] Hermitian, metric-Hermitian, and complex non-Hermitian representations select the corresponding eigensolver;
- [x] non-Hermitian reconstruction retains left/right vectors and uses an explicit right-vector inverse;
- [x] Goldstone assignment requires isolated proximity, LR-03 metric overlap, eta/k-mesh consistency, and no gapping field;
- [x] weak-overlap and two-small-mode fixtures refuse correction;
- [x] only the identified eigenvalue is changed and untouched eigenvalues remain unchanged;
- [x] BES corrected denominator and corrected Kxc are reconstructed and verified;
- [x] raw and corrected denominator/kernel objects are retained;
- [x] KXC units/provenance mismatch refuses correction;
- [x] exact, perturbed, metric, non-Hermitian, ambiguity, optionality, and reconstruction tests pass;
- [x] `docs/GOLDSTONE_EIGENVALUE_CORRECTION.md` created;
- [x] no global rescaling, bare-susceptibility modification, default enablement, Dyson solve, or material-validation claim added.

## Implementation outcome

`source/lr_goldstone_correction.f90` implements the optional BES route in the
positive-measure LR-04 response space.  It forms the canonical static
denominator, classifies its actual metric, applies the correct LAPACK
eigensolver, enforces explicit campaign/provenance/field/consistency gates,
changes only the selected eigenvalue, reconstructs the corrected kernel through
the published relation, and retains raw/corrected diagnostics.

`UnitLrGoldstoneCorrection` passes the exact, controlled-perturbation,
weak-overlap, ambiguous, metric-Hermitian, complex non-Hermitian, optionality,
provenance, untouched-spectrum, rigid-mode, and kernel-reconstruction checks.

## Deliverable

`docs/GOLDSTONE_EIGENVALUE_CORRECTION.md`

## Commit

`td-dft: add optional Goldstone eigenvalue correction`
