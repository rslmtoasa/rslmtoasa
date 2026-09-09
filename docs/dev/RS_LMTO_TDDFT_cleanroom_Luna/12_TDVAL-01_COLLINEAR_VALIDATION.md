# TDVAL-01 — Literature-locked collinear validation on bcc Fe and fcc Ni

## Prerequisites
All preceding collinear tasks PASS.

## Goal
Validate the rebuilt response before non-collinear/longitudinal expansion. This task changes no physics equations.

## Reference states
Use clean, converged crystalline restarts.

Record:
- XC functional;
- lattice and ASA-sphere convention;
- radial ground-state data;
- local moment;
- Fermi level;
- exact Git commit;
- clean working tree;
- binary hash.

## Convergence hierarchy

1. **Radial grid/basis:** static chiKS, raw Goldstone residual, low-energy response.
2. **Angular content:** only the range approved in LR-02; demonstrate convergence where a hierarchy exists.
3. **k mesh:** at least three meshes.
4. **eta:** at least three broadenings.
5. **q:** first mesh-commensurate q as numerical oracle, then arbitrary q after convergence.

## Required routes
For both Fe and Ni compare separately:
1. raw direct radial ALSDA;
2. sum-rule interaction;
3. published eigenvalue correction if explicitly selected.

No hybrid.

## Invariants
- q=0 rigid-rotation null vector;
- q->-q covariance;
- global spin reversal;
- static/dynamic bridge;
- no odd-in-q contamination in centrosymmetric SOC-off systems after convergence;
- quadratic long-wavelength ferromagnetic branch.

## Literature comparison
Compare converged dispersion/loss behavior with published Fe/Ni TDDFT results from the reference literature.

Do not tune parameters to match.

## Pass rule
A route is validated only if algebraic invariants, numerical convergence and physical energy scale all pass.

Failure of one route does not authorize modifying another route.

## Report
`docs/TDDFT_COLLINEAR_REVALIDATION.md`

## Checklist
- [ ] clean Fe restart
- [ ] clean Ni restart
- [ ] radial convergence
- [ ] angular convergence/status
- [ ] k-mesh ladder
- [ ] eta ladder
- [ ] commensurate-q oracle
- [ ] arbitrary-q convergence
- [ ] q parity
- [ ] spin reversal
- [ ] static/dynamic bridge
- [ ] direct ALSDA result
- [ ] sum-rule result
- [ ] corrected result if selected
- [ ] literature comparison
- [ ] no empirical tuning

## Commit
`tests: validate radial collinear TDDFT response`
