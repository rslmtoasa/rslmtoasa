# GCR-01 — Implement the published Goldstone eigenvalue correction as a separate numerical option

## Prerequisites
Direct ALSDA Kxc and raw q=0 response must work.

## Goal
Implement the Buczek–Ernst–Sandratskii 2011 numerical correction exactly as a distinct option.

Do not enable it automatically.

## Suggested file
`source/tddft_goldstone_correction.f90`

## Algorithm
From BES-04:
\[
D=I-\chi_{\rm KS}(0)K_{\rm xc}.
\]

1. diagonalize D;
2. identify the small eigenvalue associated with the rigid-rotation vector;
3. require strong physical overlap with ground-state magnetization;
4. set only that eigenvalue to zero;
5. leave all other eigenvalues unchanged;
6. reconstruct D_corr;
7. derive corrected Kxc from the published relation;
8. use corrected Kxc for all q,omega only when this option is selected.

## Refusal criteria
Do not correct if:
- the expected Goldstone eigenvector cannot be identified;
- multiple unrelated small eigenvalues make assignment ambiguous;
- the raw error is not plausibly numerical;
- SOC/external field physically gaps the mode.

First implementation: SOC-off collinear only.

## Tests
- controlled numerical perturbation of an exact Goldstone system;
- verify only intended eigenvalue changes;
- verify other spectrum unchanged;
- q=0 pole restored;
- deliberately wrong kernel rejected if eigenvector character is wrong.

## Checklist
- [ ] follows BES-04
- [ ] optional, not default
- [ ] rigid-rotation character gate
- [ ] only intended eigenvalue changed
- [ ] wrong-physics fixture rejected
- [ ] corrected Kxc reused consistently for all q,omega only when selected

## Commit
`td-dft: add published Goldstone eigenvalue correction`
