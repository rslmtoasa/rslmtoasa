# TDVAL-01 — Literature-locked collinear validation on bcc Fe and fcc Ni

## Prerequisites

Required:

- LR-03 PASS;
- LR-04 PASS;
- LR-05 PASS;
- LR-06 PASS;
- KXC-01 PASS;
- TDDY-01 PASS.

Strongly preferred before final sign-off:

- LR-GF-02 PASS;
- GSR-01 PASS.

GCR-01 is optional and must remain explicitly selected.

## Nature

VALIDATION CAMPAIGN.

Do not change physics equations to improve agreement.

If a failure is found, document it and return to the responsible implementation
task.

## Goal

Validate the rebuilt collinear Pauli/no-SOC response on canonical itinerant
ferromagnets before any noncollinear/longitudinal expansion.

Use at least:

- bcc Fe;
- fcc Ni.

## Provenance

For every reference state record:

- exact Git commit;
- clean/dirty worktree;
- input deck;
- structure/lattice constant;
- ASA sphere setup;
- XC functional and provenance;
- converged moment;
- Fermi level;
- electronic temperature;
- reciprocal mode;
- SR→Pauli closure metric from the foundation campaign;
- binary/compiler information.

Do not mix results from different ground states in one convergence ladder.

## Validation hierarchy

### Stage A — static invariants

Before plotting a magnon dispersion establish:

- LR-03 covariance;
- q=0 raw Ward/Goldstone residual;
- rigid-rotation overlap;
- comparison of direct ALSDA and sum-rule interaction;
- optional corrected route;
- static spectral/GF backend agreement if LR-GF-02 is available.

### Stage B — radial/angular convergence

The direct production radial mesh is the reference.

Test any allowed radial decimation only if explicitly implemented as a numerical
optimization.

For angular content, distinguish:

- complete product-space cutoff;
- converged reduced cutoff.

Never call a reduced cutoff exact.

### Stage C — k mesh

Use at least three systematically increasing k meshes.

Track:

- static response norm;
- Goldstone residual;
- representative finite-q spectra;
- low-energy peak position.

### Stage D — eta

Use at least three broadenings.

Do not extrapolate eta→0 unless the procedure is explicitly documented.

Separate true Landau damping from artificial eta broadening.

### Stage E — q

Start with:

1. q=0;
2. mesh-commensurate finite q;
3. arbitrary q;
4. long-wavelength sequence.

Verify BZ-folding/gauge invariance in production calculations.

## Interaction routes

Keep separate plots/tables for:

1. direct radial ALSDA;
2. Lounis sum-rule interaction;
3. optional BES-corrected ALSDA.

No hybrid curves.

## Required invariants

For the certified centrosymmetric no-SOC systems assess:

- q ↔ -q;
- spin reversal;
- static/dynamic bridge;
- absence of converged odd-in-q contamination;
- q=0 Goldstone behavior for routes expected to enforce it.

## Long-wavelength dispersion

Only after convergence has been established, identify the collective branch from
the loss matrix.

Demonstrate the expected small-q behavior

\[
\omega(q)=Dq^2+O(q^4)
\]

for the ferromagnetic acoustic mode.

Do not force a quadratic fit over q values outside the asymptotic regime.

Report fit window sensitivity.

## Damping

Separate:

- explicit eta;
- intrinsic width associated with the Stoner continuum.

Do not subtract eta by an undocumented formula.

## Literature comparison

Compare with published Fe/Ni TD-DFT results using the closest compatible physical
setup.

Document differences in:

- lattice;
- XC;
- basis;
- ASA/full potential;
- SOC;
- temperature/broadening.

Do not tune eta, kernel scale, or moments to match literature.

## Cross-backend validation

If LR-GF-02 exists, compare the reciprocal spectral and GF bare response at
representative q/omega points before Dyson.

Agreement of the two shared-vertex routes does not independently validate the
radial vertex; state that limitation.

## Failure triage

If the acoustic mode is not quadratic or Goldstone behavior is wrong:

1. do not patch TDVAL;
2. identify whether failure originates in:
   - transition vertex;
   - occupations;
   - chiKS;
   - kernel;
   - metric/Dyson;
   - correction route;
3. return to that task.

## Report

Create:

`docs/TDDFT_COLLINEAR_REVALIDATION.md`

Use explicit evidence categories:

- algebraic;
- numerical convergence;
- physical material validation;
- literature comparison.

## PASS rule

A route is validated only when all applicable categories pass.

A passing corrected route does not rescue a failing raw direct route.

A passing sum-rule route does not validate ALSDA.

## Completion record

- [x] branch, starting HEAD, worktree state, compiler, and binary provenance recorded;
- [x] LR-03 through TDDY-01 implementation-level prerequisites re-audited;
- [x] LR-01/LR-02N radial and Pauli closure gates rerun;
- [x] historical pre-purge Fe/Ni TD-DFT response evidence excluded;
- [x] production-driver availability checked and the current blocker documented;
- [x] `docs/TDDFT_COLLINEAR_REVALIDATION.md` created;
- [ ] current rebuilt bcc-Fe response reference available;
- [ ] current rebuilt fcc-Ni response reference available;
- [ ] Fe/Ni static invariants and interaction-route comparison completed;
- [ ] radial/angular, k-mesh, eta, and q convergence ladders completed;
- [ ] converged collective branch, damping, or stiffness promoted;
- [ ] physical material validation PASS claimed.

The unchecked material items are intentional.  The current branch contains the
certified response modules and finite fixtures, but no production caller that
connects an SCF/eigenpair material state to those modules.  TDVAL-01 therefore
records a BLOCKED material-validation campaign rather than promoting legacy
pre-purge spectra or synthetic fixture results.

## Commit

`tests: validate rebuilt collinear TDDFT response`
