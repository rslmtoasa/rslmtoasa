# TDDFT-09 — three-backend equivalence campaign and invariant harness

## Agent mission

Build a formal cross-backend verification harness so eigenpairs, K-space GF and native R-space GF are treated as independent numerical realizations of the same chi_KS rather than three loosely comparable plotting modes.

## Context

The strongest benefit of supporting all three backends is diagnostic independence. Pointwise matrix comparisons should precede spectral peak comparisons.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Define common deterministic fixtures and convergence settings.
2. Compare complex matrix elements of chi0 at selected q,w, including w=0 and points inside/outside the Stoner continuum.
3. Compare full matrices using absolute/relative norms and eigenvalues.
4. Compare Ward residuals and spectral sum-rule residuals backend by backend.
5. Establish separate convergence axes: k mesh, Rmax, energy/contour resolution, eta/gamma.
6. Produce a small machine-readable evidence table committed under `results/` or the project’s established evidence location.
7. Add negative controls that deliberately flip a sign/factor in a test-only path and verify the harness catches it.
8. Do not set overly loose tolerances just to make all backends pass; explain differences and converge them.

## Acceptance checklist

- [x] All three backends compared pointwise.
- [x] Static and dynamic points included.
- [x] Matrix-norm/eigenvalue comparisons included.
- [x] Ward and sum-rule diagnostics compared.
- [x] Negative controls demonstrate test sensitivity.
- [x] Tolerances are justified by convergence evidence.
- [x] Evidence is machine-readable and documented.

## Completion record

Completed with `UnitTddftBackendEquivalence` and the committed evidence file
[`results/tddft_backend_equivalence.json`](../../../../results/tddft_backend_equivalence.json).
The deterministic periodic fixture covers `q=(0,0.25,0.5)`, a frequency
inside and outside its Stoner support, exact static eigenpair/K-GF agreement,
finite-eta `omega=0` Ward diagnostics for all three routes, and independent
k-mesh, R-cutoff, energy, contour, and eta/gamma convergence axes.  The native
sampled-real-axis route explicitly advertises exact static as unsupported;
that capability boundary is part of the check.  The summary and scope limits
are documented in [`docs/TDDFT_BACKEND_EQUIVALENCE.md`](../../../TDDFT_BACKEND_EQUIVALENCE.md).

## Required evidence / deliverables

Add `docs/TDDFT_BACKEND_EQUIVALENCE.md` summarizing what is now proven and what remains backend-specific.

## One-line commit message

`tddft: cross-validate chi0 backends`
