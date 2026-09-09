# TDDY-01 — Implement enhanced susceptibility and loss matrix

## Prerequisites
LR-06 and at least KXC-01 PASS.

## Goal
Implement the enhanced susceptibility only after chiKS and an interaction route are independently validated.

## Suggested file
`source/tddft_dyson.f90`

## Physics
Use BES-01:
\[
\chi=(I-\chi_{\rm KS}K)^{-1}\chi_{\rm KS}.
\]

No pair-Xi route.

Supported interaction routes remain explicit and separate:
- direct radial ALSDA;
- radial sum-rule interaction;
- optional eigenvalue-corrected ALSDA.

## Loss
Implement the anti-Hermitian/loss convention exactly as fixed in LR-03.

Do not force opposite circular sectors onto the same frequency half-axis or alter signs to make all spectra positive.

## Numerical rules
- use stable linear solves where appropriate;
- report conditioning of the Dyson denominator;
- distinguish near-singular physics from numerical failure;
- no mode-fitting heuristics in this task.

## Tests
- scalar analytic Dyson;
- small matrix analytic Dyson;
- q=0 Goldstone divergence;
- finite-q nonsingular case;
- covariance under q/omega/spin reversal;
- interaction routes compared but never conflated.

## Checklist
- [ ] Dyson equation exact
- [ ] interaction routes separate
- [ ] no pair-Xi
- [ ] loss convention exact
- [ ] conditioning reported
- [ ] covariance tests pass
- [ ] no mode finder

## Commit
`td-dft: add enhanced susceptibility and loss matrix`
