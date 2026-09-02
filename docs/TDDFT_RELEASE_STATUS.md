# TDDFT release status

**Status date:** 2026-09-02
**Decision:** implementation-level fixture validation is complete; the Fe/Ni
material release gate is **not passed**.

The current checkout may describe the transverse subsystem as a validated
framework on deterministic fixtures. It must not describe Fe/Ni transverse
spectra, or the three-backend material path as a whole, as
production-validated. The open material gates are recorded in
[TDDFT_FE_NI_VALIDATION.md](TDDFT_FE_NI_VALIDATION.md).

## What is validated

| Area | Status | Evidence |
| --- | --- | --- |
| Pauli, circular, source/measurement, frequency, and loss conventions | Validated on analytic fixtures | `UnitResponseConventions`, `UnitResponseVertices`, `UnitTddftWardConventions` |
| Explicit eigenpair `chi0` reference | Validated on deterministic fixtures | `UnitTddftChiKS`, `UnitTddftTransitionWorkspace`, TDDFT-03 report |
| K-space Lehmann GF and direct resolvent | Validated on deterministic fixtures | `UnitTddftGreenChiKS`, `UnitLehmannChain`, `KSPACE_GF_VALIDATION.md` |
| Mixed contour vs direct K-GF integration | Validated on deterministic fixtures | `UnitTddftContour`, `TDDFT_GF_INTEGRATION.md` |
| Native real-space `G(R)` response, exact static contour, and q-batched Fourier transform | Validated on deterministic fixtures | `UnitTddftRealspaceGF`, `UnitTddftBackendEquivalence`, `TDDFT_REALSPACE_GF.md` |
| Common backend interface and fixture equivalence | Validated within recorded fixture envelopes | `UnitTddftBackend`, `UnitTddftBackendEquivalence`, `TDDFT_BACKEND_EQUIVALENCE.md` |
| Ward/Goldstone diagnostics and controlled correction policies | Validated as algebraic/software behavior | `UnitTddftWard`, `UnitTddftGoldstone`, `TDDFT_GOLDSTONE_WARD.md` |
| Loss-matrix modes, Stoner classification, and linewidth gates | Validated as post-processing behavior | `UnitTddftDysonModes`, `TDDFT_SPECTRAL_ANALYSIS.md` |
| AF/ferri ordered circular channels | Validated on a two-sublattice analytic fixture | `UnitTddftAfFerriChirality`, `TDDFT_AF_FERRI_CHIRALITY.md` |
| MPI work decomposition and native R-space q reuse | Validated by the recorded performance campaign | `UnitTddftPerformance_mpi`, `TDDFT-14_PERFORMANCE_MPI_REPORT.md` |

These are software and controlled-fixture claims. They establish contracts,
signs, factors, data ownership, numerical controls, and failure behavior;
they do not replace a converged material campaign.

## Current material evidence

The reproducible Fe/Ni decks and raw outputs are documented in
[TDDFT_USER_GUIDE.md](TDDFT_USER_GUIDE.md) and
[TDDFT_FE_NI_VALIDATION.md](TDDFT_FE_NI_VALIDATION.md). The bounded TDDFT-11
probe used `nsp=1`, SOC off, `ham_only`, site projection, a `4x4x4` response
mesh, three q points, and three frequencies.

| Material | Three-backend bare-`chi0` gate | Largest raw K-space Ward residual | Material conclusion |
| --- | ---: | ---: | --- |
| bcc Fe | FAIL; max difference `46.36215447` | `0.2570132937` | Not production validated |
| fcc Ni | FAIL; max difference `61.13543750` | `0.0321499941` | Not production validated; retained restart has no connected fcc hopping neighbor |

The native real-space runs exercised the intended native path, but their
bounded clusters did not include an R-tail convergence sweep. The current
material probe also has no demonstrated quadratic small-q enhanced mode,
independent same-ground-state LKAG/frozen-magnon/GBT stiffness comparison,
Ni path/backfolding audit, or eta/k-mesh/frequency convergence envelope.
Those omissions are reasons to keep the release gate open, not reasons to
apply a spectral shift or empirical kernel rescaling.

## Production support boundary

The supported initial boundary is collinear, SOC-free, first-order,
orthogonal `ham_only` response with site projection. The three backend
choices and their current limitations are summarized in
[TDDFT_USER_GUIDE.md](TDDFT_USER_GUIDE.md). Generalized overlap/HOH, SOC and
noncollinear response, Hubbard U/V, GBT/explicit textures, CCOR, and
site-orbital projection are rejected explicitly. Native RS provides direct
bare `chi0` and the validated exact static input for the transverse Ward
diagnostic; its mixed-contour and Dyson paths remain experimental/deferred.

The longitudinal and full charge-spin routes are capability-gated and
fixture-tested, but are not part of this transverse production decision.
Sternheimer/Liouville-Lanczos and relativistic/noncollinear response remain
design tracks.

## Verification commands

The closeout verification is:

```sh
cmake --build build --target rslmto.x UnitTddftBackend UnitTddftGreenChiKS \
  UnitTddftContour UnitTddftRealspaceGF UnitTddftWard UnitTddftDysonModes \
  UnitTddftConfig UnitTddftBackendEquivalence UnitTddftAfFerriChirality -j2
ctest --test-dir build --output-on-failure -L tddft
python3 tests/unit/test_tddft_dispatch.py
python3 tests/regression/tddft_validation/test_validation.py
```

The material collector is read-only with respect to the raw inputs and can be
rerun with:

```sh
python3 tests/regression/tddft_validation/materials/collect_tddft11_evidence.py
```

Every production response file includes a provenance block with the
requested/canonical backend, implementation, q and k grids, endpoint or
native-R description, temperature/Fermi/electron-count data, band window,
integration and broadening controls, R-space tail controls, source/kernel
provenance, Goldstone policy, output switches, and MPI ownership. Runtime
`TDDFT_PERF_PLAN` output records the decomposition strategy. Completion of a
run is not reported as convergence.

## Re-opening the material release gate

Before calling the transverse path production-validated, rerun Fe and Ni from
consistent converged ground states and require all of the following:

- all three `chi0` backends agree within a documented converged uncertainty;
- raw SOC-free q=0 Ward/Goldstone residual decreases under k, energy, eta,
  and response-basis refinement;
- the acoustic branch is gapless and quadratic at small q;
- stiffness agrees with independent same-ground-state LKAG and
  frozen-magnon/GBT references;
- Stoner continuum and linewidth behavior are stable under the convergence
  ladders; and
- the Ni neighbor connectivity and reciprocal path/backfolding are audited.

Until then, use `goldstone_mode = 'diagnose'`, retain raw products, and treat
`sum_rule`/`projected` or `goldstone_mode = 'correct'` as explicitly audited
numerical cleanup only. Empirical global lambda factors and arbitrary energy
shifts remain forbidden.
