# Optional BES Goldstone eigenvalue correction

## Status and claim level

The optional Buczek–Ernst–Sandratskii (BES) numerical correction is implemented
in [`source/lr_goldstone_correction.f90`](../source/lr_goldstone_correction.f90)
and registered as the standalone `UnitLrGoldstoneCorrection` test.  The route is
disabled by default and can only run when the caller explicitly sets
`request%enabled = .true.` and supplies all selection gates.

The focused tests establish algebraic consistency and independent finite-matrix
numerical evidence.  They do not establish a converged Fe/Ni spectrum,
material Goldstone accuracy, or literature-spectrum agreement.

The published contract is BES-04, Buczek, Ernst, and Sandratskii, Phys. Rev. B
84, 174418 (2011), Eqs. (38–39): diagonalize the static denominator, replace
only the Goldstone eigenvalue by zero, reconstruct the denominator, and obtain
the corrected kernel from the same calculated static susceptibility.

## Canonical LR-04 representation

The evaluator consumes the canonical right-weighted LR-04 matrices
`static_susceptibility` and `xc_kernel`.  It forms

\[
D = I - \chi_{KS}(0)K_{xc}
\]

with `response_compose_operators`; no raw pointwise matrix is multiplied and no
local quadrature weight is inserted by GCR-01.  The raw full-space `D` and
`Kxc` are retained in the result.

The response mesh contains the exact LR-04 null-measure origin.  The eigensystem
and correction use every coordinate with `metric_weights > 0`, preserving the
LR-04 order.  The corrected full-space objects use the identity/zero extension
on the null-origin coordinate; no inverse origin weight or epsilon floor is
introduced.

## Actual eigensolver selection

The implementation checks the active representation rather than assuming that
the physical mode is real.  It reports normalized ordinary-Hermitian and
metric-Hermitian defects.

| representation | condition | solver and reconstruction |
|---|---|---|
| Hermitian | `D = D^H` within `1e-10` | `zheev`; inverse of the returned eigenvector matrix |
| metric-Hermitian | `D^H W = W D` within `1e-10` | `zheev` on `W^(1/2) D W^(-1/2)`, transformed back; metric left vectors retained |
| complex non-Hermitian | neither condition | `zgeev`; right and left vectors retained, with an explicit inverse of the right-vector matrix |

For all cases the correction is reconstructed as

\[
D_{corr}=V\,\operatorname{diag}(\lambda_1,\ldots,0,\ldots,\lambda_n)V^{-1}.
\]

The returned `eigenvector_inverse` is used for the non-Hermitian case and is
also retained for diagnostics in the other cases.  This preserves every
untouched eigenvalue and eigenvector; the only changed spectral entry is the
identified Goldstone eigenvalue.

## Goldstone identification and refusal policy

The candidate is the eigenvalue of smallest magnitude.  It is accepted only if
all of the following hold:

* its magnitude is at most `0.25` times the nearest competing eigenvalue;
* its right eigenvector has normalized LR-04 metric overlap at least `0.80`
  with the supplied LR-03 rigid-rotation vector;
* at least two consistency samples are supplied, including an eta or k-mesh
  change, and every sample retains the same isolation and rigid character;
* the caller marks LR-06 convergence as passed;
* the static point is explicitly marked q=0, omega=0, with eta and k-mesh
  provenance supplied; and
* SOC, external, and constraining fields are all absent.

The direct KXC-01 contract is also explicit: units must be `Ry bohr^3`, and the
caller must verify provenance as `KXC-01 direct ALSDA LR-03`.  Missing or
mismatched provenance is a blocker.  A blocked request retains no promoted
corrected route and reports the reason in `result%status`.

These gates refuse multiple nearby small modes, weak rigid overlap, poor
LR-06 convergence, physically gapped zero modes, and defects that are not
isolated from the rest of the spectrum.  GCR-01 never applies a global scalar
rescaling and never modifies the supplied bare susceptibility.

## BES kernel reconstruction

After the mode passes the gates, only its eigenvalue is changed:

\[
\lambda_G\longmapsto 0.
\]

The corrected kernel is obtained by solving the published relation on the
positive-measure active space:

\[
\boxed{K_{xc}^{corr}=\chi_{KS}(0)^{-1}(I-D_{corr})}.
\]

The implementation uses a linear solve with `chiKS`, not an explicit inverse,
then verifies

\[
I-\chi_{KS}(0)K_{xc}^{corr}=D_{corr}.
\]

If the active static susceptibility is singular or the reconstruction residual
exceeds `1e-10`, correction is refused.  Both `raw_denominator` /
`raw_kernel` and `corrected_denominator` / `corrected_kernel` remain available
for diagnostics.  A later Dyson implementation must consume the corrected
kernel only when this route is explicitly selected; GCR-01 does not wire it in
as a default or compute a dynamical susceptibility.

## Public request/result contract

```fortran
type(lr_goldstone_correction_request) :: request
type(lr_goldstone_correction_result) :: result

request%response_space => space
request%static_susceptibility = chi_static_canonical
request%xc_kernel = kxc_canonical
request%rigid_rotation_vector = rigid_lr03
request%enabled = .true.
request%lr06_converged = .true.
request%static_q_zero = .true.
request%static_frequency = 0.0_rp
request%eta = eta
request%k_mesh_points = nk
request%kernel_units = 'Ry bohr^3'
request%kernel_provenance = 'KXC-01 direct ALSDA LR-03'
request%kernel_provenance_verified = .true.
request%consistency_samples = samples
call evaluate_lr_goldstone_correction(request, result)
```

The result reports the representation class, selected eigensolver, raw and
corrected eigenvalues, active right/left eigenvectors and inverse, mode index,
metric overlap, isolation ratio, raw/corrected rigid residuals, untouched
eigenvalue change, and kernel reconstruction residual.

## Test evidence

`UnitLrGoldstoneCorrection` uses independent canonical finite-matrix fixtures.
It covers:

1. the disabled default no-op;
2. an exact metric-Hermitian Goldstone mode;
3. a controlled perturbation of only that mode;
4. a wrong small eigenvalue with weak rigid overlap;
5. two nearby small eigenvalues with ambiguous assignment;
6. a complex non-Hermitian denominator using `zgeev` and a right-eigenvector
   inverse;
7. unchanged non-Goldstone eigenvalues and exact corrected rigid-mode action;
8. BES kernel reconstruction; and
9. KXC provenance mismatch refusal.

Focused verification:

```text
cmake --build build --target UnitLrGoldstoneCorrection -j2
ctest --test-dir build --output-on-failure -R '^UnitLrGoldstoneCorrection$'
```

Result: `1/1` CTest passed.  The executable reported zero corrected rigid-mode
residual and zero untouched-eigenvalue change in the exact and perturbed
metric-Hermitian fixtures, zero kernel-reconstruction residual in both the
metric and non-Hermitian fixtures, and `BLOCKED` for the weak-overlap,
two-small-mode, and provenance-mismatch controls.

This is algebraic/finite-fixture evidence only.  The route remains downstream
of converged material validation and separate from the Lounis sum-rule
interaction.

## Scope boundary

GCR-01 provides only the optional static BES correction.  It does not change
the LR-06 susceptibility, implement Dyson/loss matrices, infer convergence
from a single calculation, correct SOC or field-gapped modes, use the GSR-01
interaction, or claim Fe/Ni validation.
