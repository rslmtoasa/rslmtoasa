# TDDY-01 — Implement enhanced susceptibility and loss matrix

## Prerequisites

- LR-06 PASS.
- LR-04 PASS.
- At least one interaction route PASS:
  - KXC-01 direct ALSDA; or
  - GSR-01 sum-rule interaction.
- GCR-01 is optional.

## Goal

Implement the interacting transverse susceptibility and loss matrix without
introducing new physics choices.

## Dyson equation

Use the exact LR-03 convention and LR-04 operator representation.

Schematically,

\[
\chi=(I-\chi_{KS}K)^{-1}\chi_{KS},
\]

but derive the discrete production equation from LR-04 before coding.

No duplicated radial weights.

No implicit conversion between raw and weighted matrices.

## Interaction routes

Expose explicit choices:

- `direct_alsda`;
- `goldstone_sumrule`;
- `direct_alsda_goldstone_corrected` if GCR-01 is selected.

Do not create an automatic hybrid.

Store/report which route generated every result.

## Numerical solve

Prefer solving the linear system over explicitly forming an inverse.

Report conditioning of the Dyson denominator.

Distinguish:

- physically meaningful near-singular collective pole;
- numerically singular/ill-conditioned solve.

Do not regularize a physical pole away.

## Frequency handling

Use LR-03 retarded convention.

Each result must retain:

- q;
- omega;
- eta;
- interaction route;
- electronic-state provenance;
- response-space metadata.

## Loss matrix

Implement the anti-Hermitian/loss convention exactly as fixed in LR-03.

Do not negate a circular sector merely to make plots positive.

Do not force `chi+` and `chi-` onto the same frequency half-axis if the convention
does not support that.

## Tests

### Analytic scalar Dyson

Compare with closed form.

### Small matrix Dyson

Use a noncommuting `chi0` and `K` with an independently evaluated solve.

### Metric-sensitive fixture

Demonstrate that the production result matches explicit continuum-discretized
quadrature and differs from an intentionally wrong weight-free multiplication.

### Goldstone fixture

For an interaction route with an exact static Goldstone identity, demonstrate the
expected q=0 singular behavior under controlled finite eta.

### Covariance

Verify LR-03 q/frequency/spin-reversal covariance for the enhanced response.

### Interaction separation

For the same fixture, run different interaction routes and verify metadata and
objects remain separate.

## No mode finder

Do not:

- fit Lorentzians;
- track peaks;
- assign magnon branches;
- compute stiffness.

Those belong to validation/post-processing.

## Completion checklist

- [x] Dyson equation exact in the LR-04 canonical representation
- [x] Explicit direct ALSDA, Goldstone sum-rule, and corrected-ALSDA routes
- [x] No pair-Xi enhancement route
- [x] LR-03 anti-Hermitian/loss convention exact, with the LR-04 metric adjoint
- [x] Denominator conditioning and singular/near-pole status reported
- [x] q/frequency/circular and spin-reversal covariance tests pass
- [x] No mode finder or stiffness calculation added
- [x] `docs/TDDFT_DYSON_AND_LOSS.md` completed

## Deliverable

`docs/TDDFT_DYSON_AND_LOSS.md`

## Commit

`td-dft: add enhanced susceptibility and loss matrix`
