# GSR-01 — Implement the published ASA Goldstone sum-rule interaction as an independent route

## Prerequisites
KXC-01 and all LR tasks PASS.

## Goal
Implement the Lounis et al. ASA sum-rule construction as a separate published interaction route.

Do not merge it with direct ALSDA or the eigenvalue correction.

## Suggested file
`source/tddft_goldstone_sumrule.f90`

## Physics
Use LCMM-01 through LCMM-04 after the exact LR-03 normalization mapping.

Operate in the radial/angular representation corresponding to the published sum rule.

## Critical normalization
The published 4pi factor belongs to their spherical-harmonic/radial convention.

Use only the mapped RS-LMTO form. Never add/remove 4pi by numerical experiment.

## Algorithm
Construct the approved discrete form of:
\[
\Gamma U=m_z.
\]

Use a stable solve appropriate to the derived matrix.

If singular/ill-conditioned beyond what the publication supports:
- report condition numbers;
- do not add regularization by choice;
- BLOCK if necessary.

## Independence
This route must not call the direct ALSDA kernel and rescale it.

Construct it independently from:
- chiKS(0);
- m(r);
- the mapped sum-rule equation.

Compare afterward.

## Tests
- analytic discretized Gamma/m problem;
- one-site radial case;
- multisite unequal moments;
- Goldstone null vector;
- comparison to direct radial Bxc/m as a diagnostic only.

## Checklist
- [ ] ASA mapping used
- [ ] 4pi mapped by derivation
- [ ] U solved from Gamma U=m
- [ ] no ad-hoc rescale
- [ ] conditioning reported
- [ ] multisite test
- [ ] Goldstone identity
- [ ] direct ALSDA remains independent

## Commit
`td-dft: add radial Goldstone sum-rule interaction`
