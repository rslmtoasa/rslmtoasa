# TDDFT-03 — stabilize the explicit eigenpair chi0 baseline

## Agent mission

Make the current explicit transition/eigenpair backend a trustworthy reference implementation before adding new GF backends.

## Context

All later backend tests need one transparent baseline. The eigenpair route should remain even if a GF route is faster because it makes individual transitions and q mapping inspectable.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Use the audit and convention documents to repair any remaining k/k+q provenance, occupation, spin/channel or projection inconsistencies.
2. Verify the q-path coordinates and reciprocal-basis mapping, especially around Gamma.
3. Handle degeneracies and small energy denominators robustly with the correct finite-temperature/static limits.
4. Separate numerical broadening eta from physical linewidth reporting.
5. Add a true static chi0(q,0) evaluation path required by Ward tests.
6. Add deterministic output/provenance for k mesh, q, omega, eta, Fermi level and projection.
7. Reproduce small Fe/Ni spectra and document convergence trends; do not attempt to “fit” literature yet.
8. Profile the backend and record the dominant loops, but do not optimize materially in this prompt.

## Acceptance checklist

- [x] Eigenpair chi0 passes analytic fixtures.
- [x] q mapping is unit/functional tested.
- [x] Static limit is explicit and stable.
- [x] No endpoint provenance is lost.
- [ ] Fe/Ni runs are reproducible under fixed inputs.
- [x] Broadening is clearly labelled numerical.
- [x] Profile evidence is recorded.

## Required evidence / deliverables

Write a short baseline report with representative chi0 matrix elements and spectra that later GF prompts must reproduce.

## One-line commit message

`tddft: stabilize eigenpair chi0 baseline`
