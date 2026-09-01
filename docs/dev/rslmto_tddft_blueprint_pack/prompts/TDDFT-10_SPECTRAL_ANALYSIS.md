# TDDFT-10 — loss-matrix modes, linewidths, and Stoner diagnostics

## Agent mission

Harden the common post-chi0 physics layer so magnon energies, mode character and Landau damping are extracted from the matrix response in a way consistent with the Halle literature.

## Context

A visually bright spectral ridge is not enough. The Halle approach analyzes the anti-Hermitian/loss matrix, collective eigenmodes, Stoner continuum and linewidths. This is particularly important for Ni and for future multi-sublattice/chiral systems.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Define the loss matrix exactly in the project convention and verify its Hermiticity/sign.
2. Diagonalize it at each q,w and retain eigenvalues/eigenvectors.
3. Track modes along q using overlap of eigenvectors plus energy continuity; avoid branch swapping at crossings by using mode character.
4. Fit peak energy/FWHM only when a resonance passes objective quasiparticle-quality criteria; otherwise label it overdamped/continuum-like.
5. Retain the KS loss spectrum as the Stoner-continuum reference.
6. Add optional Landau-map style output for debugging q-resolved Stoner weight.
7. Ensure all post-processing is backend-independent.
8. Add synthetic two-mode crossing and broadened-resonance tests.

## Acceptance checklist

- [x] Loss-matrix convention documented/tested.
- [x] Mode eigenvectors retained.
- [x] Branch tracking survives a synthetic crossing.
- [x] FWHM fitting rejects non-Lorentzian/overdamped cases appropriately.
- [x] KS Stoner and enhanced spectra are both output.
- [x] Results are identical across chi0 backends when fed the same chi0.

## Required evidence / deliverables

Document spectral definitions and mode-quality criteria in `docs/TDDFT_SPECTRAL_ANALYSIS.md`.

Completion evidence: the loss convention, retained mode data, branch assignment,
quality gates, KS/enhanced output, and backend-independent analysis are documented
in `docs/TDDFT_SPECTRAL_ANALYSIS.md`. They are covered by
`UnitTddftDysonModes` and the TDDFT response-convention/backend-equivalence tests,
including synthetic two-mode crossing and broadened-resonance cases.

## One-line commit message

`tddft: harden loss-matrix magnon analysis`
