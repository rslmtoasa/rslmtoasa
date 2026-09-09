# LR-05 — Implement the approved radial/angular response basis

## Prerequisites
LR-01 through LR-04 PASS.

## Goal
Implement only the response coordinate system derived in LR-02/03. No chiKS accumulation yet.

## Suggested file
`source/lr_response_basis.f90`

## Representation
Use exactly the approved basis/discretization from `docs/LR_RESPONSE_BASIS_MAPPING.md`.

The data model must naturally support a future super-index
\[
\lambda=(i,s,L,\rho)
\]
without implementing non-collinear physics yet.

## Required capabilities
- deterministic index mapping;
- site/sphere metadata;
- angular L;
- radial index/weights;
- overlap metric if nonorthogonal;
- exact orthonormalization only if required by the approved mapping;
- projector/reconstruction helpers with documented normalization.

Do not imitate a Löwdin step from literature unless the chosen RS-LMTO basis actually requires the same metric treatment.

## Tests
- radial quadrature;
- spherical-harmonic normalization;
- metric/orthogonality;
- index round trip;
- projection/reconstruction of analytic radial-angular functions;
- multisite separation.

## Forbidden
- no chiKS;
- no Kxc;
- no Dyson;
- no site compression;
- no unapproved l truncation.

## Checklist
- [ ] matches LR-02 mapping
- [ ] radial weights correct
- [ ] angular normalization correct
- [ ] metric correct
- [ ] multisite indexing correct
- [ ] future density-channel axis possible
- [ ] no new physics approximation

## Commit
`linear-response: add radial angular response basis`
