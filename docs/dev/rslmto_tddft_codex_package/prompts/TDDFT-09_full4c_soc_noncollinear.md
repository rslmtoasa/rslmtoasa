# Codex Prompt — TDDFT-09: Full Four-Component, SOC, and Non-Collinear Response

## Mission

Generalize the already validated transverse/longitudinal engine to the full charge-spin response tensor.

Do not begin this task before TDDFT-03 through TDDFT-08 have stable tests.

## Part A — full response tensor

- [x] Promote active channels to site x `{0,x,y,z}`.
- [x] Build all `chi_KS^{mu,nu}` from generalized spinor vertices.
- [x] Add the charge-charge Hartree Coulomb kernel.
- [x] Add the full ALSDA XC derivative matrix from the common kernel provider.
- [x] Verify the collinear no-SOC limit reduces to the previously validated decoupled blocks.

## Part B — local frames for non-collinearity

- [x] Define each site's local frame from the ground-state magnetization.
- [x] Build the local ALSDA kernel there.
- [x] Rotate the kernel consistently into the global response frame.
- [x] Add global spin-rotation covariance tests.
- [x] In the no-SOC limit, verify the expected rigid-rotation zero modes.

## Part C — SOC

- [x] Allow charge/transverse/longitudinal mixing.
- [x] Preserve physical anisotropy/external-field gaps.
- [x] Replace zero-gap forcing by generalized sum-rule/zero-mode diagnostics.
- [x] Add small SOC fixtures before realistic material tests.

## Part D — compatibility

- [x] Original collinear transverse outputs must be unchanged when computed through the generalized route within tolerance.
- [x] Original longitudinal outputs must also be recovered in the decoupled limit.

## Non-goals

Do not implement the GF `chi_KS` backend in this task unless the provider interface requires a tiny preparatory change.

## Completion protocol

Tick boxes and commit once:

`TDDFT-09: generalize magnetic response to four-component LR-TDDFT`
