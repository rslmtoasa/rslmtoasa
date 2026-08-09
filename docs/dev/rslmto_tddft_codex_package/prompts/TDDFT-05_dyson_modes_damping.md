# Codex Prompt — TDDFT-05: Dyson Enhancement, Self-Enhancement, Modes, and Damping

## Mission

Turn `chi_KS` into the interacting TDDFT susceptibility and implement mode diagnostics suitable for itinerant magnons.

Read first:

- `README.md`
- `00_MASTER_BLUEPRINT.md`
- `01_PHYSICS_CONVENTIONS.md`

## Dyson equation

For every `(q,omega)`, solve

\[
[I-\chi_{\rm KS}K]\chi=\chi_{\rm KS}.
\]

Do not explicitly form a matrix inverse in production.

Also retain

\[
\Xi=\chi_{\rm KS}K.
\]

## Tasks

- [x] Implement robust complex LAPACK linear solves.
- [x] Expose/store/stream:
  - [x] `chi_KS`;
  - [x] Xi;
  - [x] enhanced chi.
- [x] Define a documented spectral/loss matrix from the anti-Hermitian part of chi.
- [x] Compute site-resolved and trace spectral weights.
- [x] Diagonalize the loss matrix when requested.
- [x] Diagonalize/analyze Xi.
- [x] Identify candidate collective modes through Xi eigenvalues near/crossing unity.
- [x] Track mode branches along q using eigenvector overlap, not only nearest peak frequency.
- [x] Implement controlled peak/linewidth fitting for isolated modes.
- [x] Report fit residual/quality.
- [x] Keep FWHM versus HWHM convention explicit.
- [x] Keep numerical eta distinct from intrinsic linewidth.
- [x] Support multi-eta checks or extrapolation workflows.
- [x] Do not force strongly non-Lorentzian/Stoner-dominated structures into a magnon fit.
- [x] Output enough information to classify:
  - [x] coherent collective magnon;
  - [x] strongly Landau-damped/incoherent enhanced mode;
  - [x] predominantly Stoner feature.

## Toy collective-mode test

Construct a tiny deterministic `chi_KS` + kernel fixture with a known enhanced pole.

Require:

- [x] enhanced pole absent in bare KS spectrum;
- [x] pole present in chi;
- [x] corresponding Xi unity mode;
- [x] correct response eigenvector.

## Completion protocol

Tick boxes and commit once:

`TDDFT-05: add enhanced susceptibility and magnon mode analysis`
