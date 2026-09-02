# RS-LMTO TD-DFT baseline status

**Status date:** 2026-09-02
**Scope:** TDDFT-R2-01 through TDDFT-R2-06
**Release decision:** implementation and controlled-fixture validation are
complete; the bcc Fe/fcc Ni material gate remains open.

## Executive summary

RS-LMTO-ASA now has three implemented transverse bare-response routes:

- explicit eigenpair transitions;
- an eigenpair-backed K-space Lehmann Green-function bubble; and
- the native real-space `G(R,z) -> chi0(R,omega) -> chi0(q,omega)` route.

The response conventions, exact static operations, backend interface, native
R-space MPI reduction, contour/truncation controls, diagnostics, and failure
guards are covered by deterministic tests. The complete serial TDDFT label
suite passed 27/27 tests, and the native R-space MPI comparison passed at 1,
2, and 4 ranks.

This does not constitute a production material release. The stored bcc Fe
probe is **PARTIAL** and the stored fcc Ni probe is **FAIL**. The three bare
`chi0` matrices do not agree at the required material tolerance, the raw
SOC-free Ward diagnostics are not converged to a gapless result, and the
required small-q, material-convergence, and independent stiffness gates are
not all satisfied. No Goldstone correction, empirical frequency shift, or
response rescaling is used to conceal those results.

In this document, **production-supported** means approved for material
results, not merely accepted by the input parser or executable. The current
software boundary is therefore narrower than the set of implemented and
fixture-tested components.

## Support matrix

| Physics or capability | Implemented | Unit-tested | Material-validated | Production-supported | Notes |
| --- | --- | --- | --- | --- | --- |
| Collinear FM, SOC=0 transverse bare response | Yes | Yes | Fe **PARTIAL**; Ni **FAIL** | **No** | Software boundary is available, but the material release gate is open. |
| Explicit eigenpair backend | Yes | Yes | No | **No** | Reference transition route; material three-backend gate fails. |
| K-space Lehmann GF backend | Yes | Yes | No | **No** | Eigenpair-backed resolvent/bubble; not the native reciprocal-GF storage path. |
| Native real-space GF backend | Yes | Yes | No | **No** | Direct bare `chi0` and exact static Ward input are available; material equivalence fails. |
| Exact static Ward operation | Yes | Yes | No | **No** | Genuine divided-difference/static-contour paths; material raw residuals remain release evidence. |
| Native R-space MPI | Yes | Yes (1/2/4 ranks) | No material MPI evidence | **No** | Deterministic complex-`chi0` layout comparison passes; material gate is independent. |
| Native R-space mixed contour | Yes, provider/API level | Yes, controlled fixture | Not assessed | **No** | Complex-energy provider and contour route are tested; the current production material source is sampled real-axis and rejects this mode. |
| R-space full-tail and production truncation | Yes | Yes | No material R-tail convergence | **No** | `R_source` and `R_max` are distinct; production truncation requires prior convergence evidence. |
| AF/ferri transverse chirality | Limited | Yes, analytic two-sublattice fixture | No | **No** | Material validation remains future work. |
| Longitudinal charge-spin response | Limited/prototype | Yes, module/capability fixtures | No | **No** | Capability-gated; no FeSe/second-sound validation. |
| Full charge-spin response | Limited/experimental | Yes, algebraic/capability fixtures | No | **No** | Selected-XC derivatives are guarded and are not a relativistic release. |
| SOC response | No response path | Guard tested | No | **No** | Unsupported before response tensors are built. |
| Noncollinear response | No response path | Guard tested | No | **No** | Torque/full spinor kernel is not implemented. |

The supported input boundary for the executable is consequently:
collinear, SOC-free, first-order, orthogonal `reciprocal_mode='ham_only'`,
site-projected transverse response. Generalized overlap/HOH, Hubbard U/V,
GBT/explicit textures, CCOR, site-orbital projection, SOC, noncollinear,
longitudinal production response, and full relativistic response remain
explicitly rejected or capability-gated.

## Backend audit

### Explicit eigenpairs

The eigenpair route is the transparent transition-sum reference. It retains
separate `k` and exact folded `k+q` endpoints, finite-temperature occupation
differences, ordered circular vertices, and a genuine zero-frequency divided
difference including the Fermi-derivative limit. It is implemented and
fixture-tested, but the Fe/Ni material response is not accepted as a
production validation because the three-backend comparison and raw Ward
gates fail.

### K-space Green function

The K-space route evaluates the retarded Green-function bubble using an
eigenpair-backed Lehmann resolvent. Direct real-axis integration is the
reference; the mixed-contour path is tested with a genuine complex-energy
source. This route is algorithmically distinct from the explicit transition
denominator but is not a claim that K-space and native R-space storage are
the same implementation. The material probe does not pass the bare-response
equivalence gate.

### Native real-space Green function

The native route consumes R-space Green blocks, forms `chi0(R,omega)`, and
Fourier transforms the response itself over a q batch. Pair ownership is
disjoint across MPI ranks; local q contributions are reduced once at the
complex response level, with no division by MPI size and no second reduction.
The exact static operation uses the independently implemented retarded/
advanced contour identity rather than a finite-eta dynamic sample.

The arbitrary-complex-energy provider and mixed-contour machinery exist and
are covered on controlled fixtures. The current production material source,
however, attaches sampled real-axis native blocks, so its supported route is
direct integration; the native production driver guards mixed contour,
longitudinal/full response, and Dyson enhancement where the required native
source/kernel is unavailable.

## Raw Goldstone/Ward status

The raw diagnostic is retained before any optional correction. The material
evidence records the following maximum K-space/native comparisons and raw
residuals:

| Material | Largest bare-backend matrix difference | Largest raw Ward residual | Decision |
| --- | ---: | ---: | --- |
| bcc Fe | `46.36215447` | `0.2570132937` | **PARTIAL** |
| fcc Ni | `61.13543750` | `1.0000000000` | **FAIL** |

The native static Ward residuals are separately retained in the raw files;
they are not replaced by a corrected value. No converged interacting acoustic
gap or quadratic small-q branch is claimed for either material.

## R2 evidence review

| Milestone | Closeout finding | Evidence |
| --- | --- | --- |
| R2-01 MPI correctness | Native R-space complex `chi0` agrees with the rank-one reference at 1/2/4 ranks; local work is disjoint and normalization is not rank-divided. | [MPI closeout](dev/RS_LMTO_TDDFT_R2_PROMPT_PACK_2026-09-02/TDDFT-R2-01_MPI_CORRECTNESS.md), [MPI unit](../tests/unit/test_tddft_realspace_mpi.f90) |
| R2-02 exact static RGF | Independent static RGF path is enabled and checked against the K-space/eigenpair fixture; material static equivalence is not passed. | [static closeout](dev/RS_LMTO_TDDFT_R2_PROMPT_PACK_2026-09-02/TDDFT-R2-02_STATIC_RGF_WARD.md), [backend evidence](../results/tddft_backend_equivalence.json) |
| R2-03 contour RGF | Complex-energy provider, direct reference, contour convergence, q reuse, and timing are covered on controlled fixtures; material production contour support is not claimed. | [contour closeout](dev/RS_LMTO_TDDFT_R2_PROMPT_PACK_2026-09-02/TDDFT-R2-03_COMPLEX_CONTOUR_RGF.md), [R-space guide](TDDFT_REALSPACE_GF.md) |
| R2-04 truncation/performance | Full-tail validation and early production truncation are distinct; source-radius insufficiency is reported; deterministic speedup and MPI behavior pass. | [truncation closeout](dev/RS_LMTO_TDDFT_R2_PROMPT_PACK_2026-09-02/TDDFT-R2-04_RSPACE_TRUNCATION_PERFORMANCE.md), [performance report](dev/TDDFT-14_PERFORMANCE_MPI_REPORT.md) |
| R2-05 Fe/Ni validation | Fe is **PARTIAL** and Ni is **FAIL**; failed gates remain visible and no correction is substituted for raw evidence. | [material report](validation/TDDFT_BCCFE_FCCNI_VALIDATION.md), [machine-readable evidence](../results/validation/TDDFT-11_FE_NI/evidence.json), [material decks](../tests/regression/tddft_validation/materials) |

## Verification performed for this closeout

The current build produced the following results:

```text
ctest --test-dir build --output-on-failure -L tddft             PASS 27/27
ctest --test-dir build-r2-mpi-openmpi --output-on-failure \
  -R 'UnitTddftRealspaceMPICorrectness_mpi_(1|2|4)'             PASS 3/3
python3 tests/regression/tddft_validation/test_validation.py   RESULT: PASS
```

All six bcc Fe/fcc Ni material decks were rerun in isolated temporary copies
with `OMP_NUM_THREADS=1`: explicit eigenpairs, K-space Lehmann, and native
R-space GF for each material. The regenerated response, Goldstone, and
manifest artifacts match the stored evidence; only profiling timing metadata
changes between runs. The machine-readable material collector reproduces the
stored backend-difference and raw-Ward values when given the stored evidence
revision label.

## Production boundary and remaining risks

The current release may claim a validated TDDFT software framework and a
guarded initial transverse input boundary. It may not claim production-
validated Fe/Ni spectra, converged three-backend material equivalence, a
gapless material acoustic branch, or a converged material stiffness. Raw
Goldstone diagnostics must remain visible, and `goldstone_mode='diagnose'`
is the appropriate material evidence mode until the release gate is reopened.

Future work is deliberately out of scope for this closeout:

- AF/ferri material validation;
- longitudinal FeSe/second-sound validation;
- Sternheimer/Savrasov response;
- SOC and noncollinear response.

The material gate can be reopened only after consistent converged ground
states, all three bare-response comparisons, systematic raw Ward refinement,
gapless quadratic small-q behavior, independent same-ground-state stiffness
references, stable Stoner/linewidth convergence, and the Ni connectivity/path
audit are complete.

## Primary documentation

- [TDDFT release status](TDDFT_RELEASE_STATUS.md)
- [TDDFT user guide](TDDFT_USER_GUIDE.md)
- [TDDFT backend contract](TDDFT_BACKENDS.md)
- [TDDFT conventions](TDDFT_CONVENTIONS.md)
- [TDDFT Goldstone/Ward diagnostics](TDDFT_GOLDSTONE_WARD.md)
- [TDDFT Fe/Ni validation report](validation/TDDFT_BCCFE_FCCNI_VALIDATION.md)
