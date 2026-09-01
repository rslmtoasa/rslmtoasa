# TDDFT-08 — production energy integration: mixed contour first, Halle analytic continuation second

## Agent mission

Replace expensive brute-force real-axis GF integration with a validated complex-energy strategy, prioritizing a Lounis-style mixed contour decomposition and retaining Halle-style analytic continuation as an optional later mode.

## Context

Both Halle and Jülich exploit complex energies to move GF evaluation away from one-electron poles. Lounis separates an analytic contour-friendly part from a small near-real interval; Halle evaluates susceptibility at complex frequency and continues to the real axis. The latter can be very efficient but analytic continuation is ill-conditioned.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Keep the direct real-axis integration as the reference mode.
2. Derive the Lounis/Mills decomposition in the code’s exact convention and full LMTO response basis.
3. Implement the analytic part on a deformable complex contour and the nonanalytic part on its finite near-Fermi interval.
4. Validate contour deformation by varying nodes/shape and comparing with direct integration on analytic and Fe/Ni fixtures.
5. Cache/batch GF evaluations where this improves measured cost without excessive memory.
6. Only after the mixed-contour path is stable, prototype the Halle finite-complex-frequency susceptibility + analytic continuation mode behind an explicit experimental option.
7. For analytic continuation, implement at least two stability checks (e.g. sampling/order variation and comparison with direct/mixed results), reject noncausal continuations, and report uncertainty.
8. Never interpret numerical eta/gamma as physical magnon damping.

## Acceptance checklist

- [x] Direct reference path retained.
- [x] Mixed contour reproduces direct chi0.
- [x] Convergence versus contour nodes is documented.
- [x] GF evaluations demonstrably move away from problematic real-axis poles.
- [x] No analytic-continuation mode is exposed; Halle remains explicitly deferred until it can be stability-checked.
- [x] Causality and the existing response normalization/spectral conventions remain satisfied on the TDDFT-08 fixture and regression set.
- [x] Measured runtime improvement is reported.

## Required evidence / deliverables

Write `docs/TDDFT_GF_INTEGRATION.md` with derivation, quadrature choices, convergence tables and performance data for both G(k) and G(R) where applicable.

TDDFT-08 evidence is recorded in `docs/TDDFT_GF_INTEGRATION.md`.  The current
native G(R) source exposes sampled real-axis blocks only, so the mixed contour
is implemented and measured for G(k); selecting it for G(R) is explicitly
rejected until a genuine complex-energy native source exists.  No Fe/Ni
fixture is available in this unit-test milestone.

## One-line commit message

`tddft: add contour integration for GF response`
