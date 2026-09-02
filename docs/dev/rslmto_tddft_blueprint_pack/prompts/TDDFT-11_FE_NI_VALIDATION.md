# TDDFT-11 — bcc Fe and fcc Ni production validation gate

## Agent mission

Perform the decisive collinear SOC-free production validation on bcc Fe and fcc Ni, including the previously problematic Ni long-wavelength regime.

## Context

These systems test itinerant metallic magnons, Stoner damping, Goldstone behavior and small-q stiffness. The results must also close the low-energy consistency triangle with independent LKAG and frozen-magnon/GBT data from the same ground state.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Prepare reproducible Fe and Ni input decks with documented ground-state convergence.
2. For selected q,w points, demonstrate converged agreement of all three chi0 backends.
3. At q=0/near q=0, document the raw Ward/Goldstone residual and any optional numerical cleanup magnitude.
4. Use a dense, explicitly documented small-q path and fit `omega=D q^2` only inside a demonstrated quadratic window.
5. Compare D and low-energy dispersion with LKAG-derived and frozen-magnon/GBT references computed consistently in the same code.
6. For Ni, explicitly test reciprocal-coordinate/path correctness and rule out backfolding/path discontinuities.
7. Map the KS Stoner continuum and enhanced magnon linewidth trend.
8. Perform convergence sweeps in k mesh/Rmax, energy contour, eta/gamma and frequency mesh.
9. Produce a release-gate report; do not claim validation if one backend disagrees or q^2 is absent.

## Acceptance checklist

- [ ] Fe passes backend-equivalence gate.
- [ ] Ni passes backend-equivalence gate.
- [ ] SOC=0 acoustic mode is gapless within converged numerical uncertainty.
- [ ] Small-q dispersion is demonstrably quadratic.
- [ ] Stiffness agrees with independent low-energy references within established uncertainty.
- [ ] Ni path/backfolding concerns are explicitly resolved.
- [ ] Damping/Stoner trends are stable under numerical convergence.
- [x] Release-gate report states **FAIL** without cosmetic tuning (see `docs/TDDFT_FE_NI_VALIDATION.md`).

## Required evidence / deliverables

Commit `docs/TDDFT_FE_NI_VALIDATION.md`, input decks and machine-readable convergence/evidence data.

## One-line commit message

`tddft: validate transverse response in Fe and Ni`
