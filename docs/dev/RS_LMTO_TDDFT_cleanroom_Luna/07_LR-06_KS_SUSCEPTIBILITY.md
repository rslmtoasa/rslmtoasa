# LR-06 — Implement collinear transverse Kohn-Sham susceptibility

## Prerequisites
LR-01 through LR-05 PASS.

## Goal
Implement the bare retarded KS susceptibility in the approved radial/angular response basis.

No Kxc, Dyson enhancement, Goldstone correction or mode fitter.

## Suggested files
- `source/lr_conventions.f90`
- `source/lr_ks_susceptibility.f90`

Only create files actually needed.

## Physics
Use only the retarded susceptibility fixed in LR-03.

Do not copy old `tddft_chi0` signs, circular conventions or endpoint ordering.

## Backend
Implement only a backend whose exact population of the approved spatial response was proved in LR-02.

No site-projected fallback.

## Requirements
- explicit occupation-state provenance;
- exact q/k+q phase;
- exact spin operator;
- radial/angular matrix elements;
- full response-space matrix;
- finite eta convention from LR-03.

## Validation before material magnons
1. q->-q covariance;
2. omega/circular covariance from LR-03;
3. global spin reversal;
4. raw static sum-rule left side using independent ground-state Bxc;
5. radial convergence;
6. k-mesh convergence;
7. backend cross-validation if a second exact backend exists.

## No Goldstone enforcement
Report raw q=0 Goldstone/Ward discrepancy. Do not fix it.

## Checklist
- [ ] approved spatial basis used
- [ ] occupation provenance explicit
- [ ] q phase exact
- [ ] spin factors exact
- [ ] q covariance passes
- [ ] spin reversal passes
- [ ] radial convergence
- [ ] k convergence
- [ ] raw static Ward evidence
- [ ] no Kxc/Dyson

## Commit
`linear-response: implement radial transverse KS susceptibility`
