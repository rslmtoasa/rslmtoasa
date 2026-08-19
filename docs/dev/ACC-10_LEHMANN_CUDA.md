# ACC-10 — GPU Lehmann Green-function contractions

ACC-10 moves the validated strict-Lehmann contraction behind the existing
reciprocal execution-backend seam.  CPU Fourier assembly and the host
eigenpair arrays remain authoritative.  A CUDA request transfers one batch of
eigenvalues, eigenvectors, k-points, complex energies, and directed pair
geometry, evaluates all requested `(pair, energy, row, column)` outputs on the
GPU, and copies the result back into the existing `green%gij/gji` arrays.

The normal contour and the established 64-point Fermi-eta contour are submitted
together.  Torque-resolved arrays continue to be derived by the existing
`pauli_decompose_block` routine, and no GPU Green object or device-residency
framework was introduced.

## Correctness evidence

The isolated CUDA contract is registered as `ReciprocalCudaLehmann`.  It
compares a complete three-pair request, including onsite and intersite phase
factors, pair offsets, all requested energies, and complex output values against
the CPU reference calculation.  The NVIDIA validation host reported:

```text
PASS: CUDA Lehmann contraction max_error=5.64785e-17
```

The production VAL-05 campaign was run with both current CPU and CUDA builds
on the fixed bcc-Fe potential.  All k-mesh, eta, and energy-window cases
passed.  The selected onsite `G_ii` and first-neighbour intersite `G_ij`
outputs agree to the printed precision, while the established Dyson(Sigma=0)
invariant remains below `1e-10` in the CUDA campaign.

The downstream CUDA triads also passed:

| Consumer | CUDA Lehmann | Dyson(Sigma=0) | Recursion |
|---|---:|---:|---:|
| bcc-Fe Jij, pair 1-335 | 0.2547380660120 | 0.2547380660129 | 0.5078779010703 |
| bcc-Fe Jij, pair 1-336 | 0.3132323415068 | 0.3132323415068 | 0.3861927945335 |
| bcc-Fe conductivity `sigma_xx` | 3.235054 | 3.235054 | 3.239956 |
| bcc-Fe damping `alpha` | 0.002527619 | 0.002527619 | 0.001341155 |

The route-specific recursion differences are the existing documented
broadening/mesh envelopes; CUDA Lehmann and Dyson remain the tight pair.

## Performance evidence

For the representative bcc-Fe `12^3` Green request (GNU Fortran 13.3,
one OpenMP thread, CUDA 13.3, NVIDIA validation GPU, full-bz orthonormal
scope), the production timing records were:

| Backend | Total contraction request | H2D | GPU contraction | D2H |
|---|---:|---:|---:|---:|
| LAPACK CPU | 16.770 s | 0 | 16.770 s | 0 |
| CUDA | 2.813 s | 0.001641 s | 2.804782 s | 0.001679 s |

This is a measured end-to-end contraction speedup of approximately `5.96x`
for this fixture, including transfers and canonical output materialization.
The isolated microfixture is intentionally small and remains transfer/launch
dominated; its timings are retained by the unit test but are not used as a
production speedup claim.

## Scope boundary

ACC-10 supports standard orthonormal Lehmann eigenpairs and the existing
full-BZ reciprocal Green route.  Generalized-overlap problems, GPU H(k)
assembly, device-resident eigensystem handoff, and a separate GPU Green object
remain outside this task.

**Checks:** `Acc10LehmannSource`, `UnitLehmannChain`,
`UnitDysonEquivalence`, `ReciprocalCudaLehmann`, production VAL-05 CPU/CUDA,
and the CUDA Jij/conductivity/damping triads.
