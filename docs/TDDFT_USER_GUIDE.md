# RS-LMTO TDDFT user guide

This is the user-facing entry point for the current linear-response TDDFT
implementation. It describes the three `chi0` backends, their numerical
controls, supported input boundary, output provenance, and reproducible Fe/Ni
smoke decks.

The implementation is validated on deterministic fixtures, but the current
Fe/Ni material release gate is still open. See
[TDDFT_RELEASE_STATUS.md](TDDFT_RELEASE_STATUS.md) before treating a material
spectrum as production validated.

## Quick start

TDDFT is an opt-in post-processing route selected with
`post_processing = 'susceptibility'` and an `&tddft` group. A minimal
transverse run is:

```fortran
&calculation
post_processing = 'susceptibility'
/
&tddft
enabled = .true.
channel = 'transverse'
chi0_backend = 'eigenpairs'
q_mode = 'list'
q_coordinates = 'direct'
q_file = 'q_points.dat'
nomega = 201
omega_min = 0.0
omega_max = 0.5
eta = 0.01
response_projection = 'site'
goldstone_mode = 'diagnose'
output_chi0 = .true.
output_stoner = .true.
/
```

All energies and frequencies are in Rydberg. `q_file` contains a positive
point count followed by three fractional reciprocal-lattice coordinates per
line. The response mesh is generated without symmetry reduction so that the
`k` and exact folded `k+q` endpoints are retained.

The route currently requires a collinear, SOC-free, first-order,
orthogonal (`reciprocal_mode = 'ham_only'`) ground-state Hamiltonian with
`nsp = 1`, `magnetic_representation = 'periodic_nc'`, and a standard supported
DFT XC path. Unsupported branches fail before response tensors are built.

## Choosing a `chi0` backend

The blueprint names the routes `eigenpairs`, `gf_k_lehmann`, and `gf_r_native`.
The accepted input names and their canonical output labels are:

| User choice | Input value | Canonical output | Best use |
| --- | --- | --- | --- |
| Explicit transitions | `eigenpairs` | `eigenpairs` | Transparent reference, static Ward diagnostics, and pair-potential Xi |
| K-space GF | `kspace_lehmann` (aliases: `green`, `kspace_gf`) | `kspace_lehmann` | Periodic calculations and independent GF-bubble validation |
| Native real-space GF | `realspace_gf` (aliases: `realspace`, `rs_gf`) | `realspace_gf` | Dense q batches, surfaces/films, and native RS workflows |

### `eigenpairs`

This is the transparent explicit transition sum. It retains the separate
`k` and `k+q` endpoint eigenpairs, supports the exact static divided
difference used by the Ward diagnostic, and is the only backend currently
allowed to build the direct LMTO pair-potential Xi route. Use it first when
debugging signs, factors, q mapping, or Goldstone behavior.

### `kspace_lehmann`

This backend evaluates the K-space Green-function bubble using the validated
orthogonal Lehmann resolvent. It is algorithmically independent of the
explicit transition denominator while sharing the same response vertices.
`gf_integration = 'direct'` is the real-axis reference. The validated
`'mixed_contour'` option deforms the same-half-plane terms and retains the
near-Fermi mixed term; it requires a genuine complex-energy source.

The current K-space source is an eigenpair-backed resolvent. It is a GF
response validation route, not a claim that the K-space and native RS storage
paths are the same implementation.

### `realspace_gf`

This is the native path

```text
G_ab(R,z), G_ba(-R,z) -> chi0_ab(R,omega) -> chi0_ab(q,omega)
```

It consumes native RS-LMTO Green blocks, builds the real-space susceptibility
once per frequency batch, and Fourier transforms `chi0(R,omega)` to the
requested q batch. It never requires an intermediate `G(k,z)`.

The current native source accepts sampled real-axis blocks, so it supports
`gf_integration = 'direct'` and bare `chi0` output only. Mixed contour, exact
static Ward input, Xi/Dyson enhancement, and longitudinal/full response are
rejected until a complex-energy native source and matching static kernel are
available.

## Support matrix

The matrix describes the current production driver boundary, not a planned
feature. “Rejected” means the run stops with an explicit capability error;
it does not silently use the collinear formula.

| Capability / input | `eigenpairs` | K-space GF | Native RS GF | Notes |
| --- | --- | --- | --- | --- |
| Collinear `nsp=1`, SOC off | Supported | Supported | Supported | Initial transverse boundary |
| SOC or noncollinear response | Rejected | Rejected | Rejected | Relativistic torque and full spinor kernel are not implemented |
| `reciprocal_mode='ham_only'` | Required | Required | Required | Orthogonal response metric |
| Generalized overlap / HOH / second order | Rejected | Rejected | Rejected | No `delta S`/metric response derivative |
| Hubbard U/V response | Rejected | Rejected | Rejected | Hubbard response kernel is not implemented |
| GBT / explicit spin texture | Rejected | Rejected | Rejected | Response tangent is not derived for this route |
| CCOR-modified Hamiltonian | Rejected | Rejected | Rejected | Response derivative is not implemented |
| Site projection | Supported | Supported | Supported | Validated production projection |
| Site-orbital projection | Rejected | Rejected | Rejected | Requires orbital-resolved ground-state kernel/moment data |
| Transverse bare `chi0` | Supported | Supported | Supported | Native RS uses direct real-axis integration |
| Transverse Xi/Dyson | Supported | Supported with legacy site kernel | Rejected | Pair-potential Xi is currently eigenpair-only |
| Exact static Ward diagnostic | Supported | Supported | Rejected | Native source is dynamic-first |
| Mixed contour | Not applicable | Supported | Rejected | Native complex-energy source is not available |
| Longitudinal `(0,z)` response | Limited | Limited | Rejected | Selected-XC/Hartree path is capability-gated; not part of transverse release |
| Full `(0,x,y,z)` response | Limited | Limited | Rejected | Selected-XC derivative capability is required and remains experimental |

The `channel='longitudinal'` and `channel='full'` routes are exposed for
validated implementation-level fixtures and selected capability checks. They
must not be described as a production longitudinal or relativistic TDDFT
release. In particular, no longitudinal Goldstone condition or LLB parameter
is inferred.

## Convergence controls

Convergence must be demonstrated for the response, not only for the ground
state or one-electron Green function. Record the q path, frequency grid,
response basis, electronic temperature, Fermi level, band window, and every
backend-specific control in the output metadata.

### Common controls

- `eta` is a numerical response broadening. It is not a physical linewidth and
  is never subtracted automatically from a mode width.
- The static Ward calculation uses the real zero-frequency Fermi divided
  difference and no dynamic `eta`.
- The response Fermi level is resolved on the complete response mesh at the
  ground-state target electron count. A mismatch is fatal.
- Increase q and frequency resolution together when fitting a small-q mode;
  do not infer a quadratic dispersion from a single finite q point.
- Compare the bare Stoner spectrum and enhanced loss spectrum. A peak without
  a well-conditioned Xi unity crossing remains a Stoner/continuum feature.

### Explicit eigenpairs

Converge `nk1,nk2,nk3`, `band_first/band_last`, frequency spacing, and `eta`.
Keep `occupation_tolerance = 0` for the transparent reference; a positive
value is an explicit transition-pruning approximation. The exact endpoint
mapping is part of the result, so do not replace the complete response mesh
with a symmetry-reduced mesh.

The batched transition accumulator is the normal implementation. The scalar
path remains a numerical oracle for small fixtures and performance checks.
For a static Ward sweep, vary the k mesh and response basis while keeping
`eta` out of the static calculation.

### K-space Green functions

Converge these controls independently:

```fortran
green_eta          = 0.0       ! uses eta/2
green_energy_min  = ...
green_energy_max  = ...
green_energy_points = 2001
gf_integration     = 'direct'  ! or 'mixed_contour'
contour_points     = 64
contour_subdivisions = 8
near_fermi_points  = 128
contour_height     = 0.0       ! safe default when unset
```

For direct integration, increase the real-axis energy points and verify that
the energy window contains the source spectrum and thermal/frequency tails.
For mixed contour, vary contour nodes, horizontal subdivisions,
near-Fermi nodes, and contour height; reject unstable or noncausal results.
The effective one-particle half-width is `green_eta`, or `eta/2` when the
input is zero, so the bubble has response broadening `eta`.

### Native real-space Green functions

Converge the real-space response with `realspace_rmax` or an explicit audited
pair list. The provider reports retained and omitted pair points, shell
count, outer-shell norm, omitted-tail norm, tail ratio, and
`real_space_tail_assessed`. A bounded local zone is a performance choice, not
a convergence claim.

For bulk use all three Fourier axes. For films, set
`realspace_representation = 'film'` and list the two in-plane axes in
`realspace_fourier_axes`; the remaining coordinate stays real-space. The
native route amortizes the real-space build over a q batch, so dense q paths
are attractive after the R-space tail is converged.

The current native source is direct-only. A requested mixed contour is an
explicit error, not a real-axis interpolation disguised as analytic
continuation.

## Fe/Ni reproducible examples

The bounded TDDFT-11 decks are under
[`tests/regression/tddft_validation/materials`](../tests/regression/tddft_validation/materials).
They use the same current-branch restart inputs, direct q coordinates, and
three backend choices for each material. From either material directory run:

```sh
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_eigenpairs.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_kspace_lehmann.nml
OMP_NUM_THREADS=1 ../../../../../build/bin/rslmto.x input_realspace_gf.nml
```

The committed raw q-resolved files and manifests make these examples
reproducible. To regenerate the machine-readable comparison from the checked
outputs:

```sh
python3 tests/regression/tddft_validation/materials/collect_tddft11_evidence.py
```

This is a bounded smoke campaign, not a release pass. The current Fe/Ni
numbers show backend disagreement, nonzero raw K-space Ward residuals, no
native R-tail sweep, and an Ni restart with no connected fcc hopping neighbor.
The complete interpretation is in
[TDDFT_FE_NI_VALIDATION.md](TDDFT_FE_NI_VALIDATION.md).

The deterministic implementation gates can be run with:

```sh
ctest --test-dir build --output-on-failure -L tddft
```

The CI fixture checker is separate and does not make a material claim:

```sh
python3 tests/regression/tddft_validation/test_validation.py
```

## Output and provenance

The response writer emits one self-describing file per q point. The headers
identify:

- requested and canonical backend, implementation, and compatibility alias;
- response channel, circular ordering, Pauli/source convention, units, and
  spectral-loss sign;
- q coordinates, k mesh, frequency grid, Fermi level, temperature, electron
  counts, band window, and occupation pruning;
- `eta`, its numerical-only role, effective `green_eta`, energy integration,
  real-axis/contour controls, and physical-linewidth policy;
- native R-space cutoff, representation, Fourier axes, tail diagnostics, and
  q-batch reuse when applicable;
- source Hamiltonian, measurement/source/kernel separation, Xi/Dyson choice,
  Goldstone policy, correction status, and unsupported-feature policy;
- MPI ownership/decomposition provenance and output switches.

The `TDDFT_PERF_PLAN` line printed at runtime records the complete work
decomposition. Raw Goldstone and raw pair-potential products remain alongside
any explicitly corrected product. No empirical global scale, energy shift, or
finite-eta static inverse is permitted.

## Developer and physics references

- [TDDFT_CONVENTIONS.md](TDDFT_CONVENTIONS.md) — Pauli, circular, source, and
  loss conventions.
- [TDDFT_BACKENDS.md](TDDFT_BACKENDS.md) — common backend interface and
  ownership contract.
- [TDDFT_GOLDSTONE_WARD.md](TDDFT_GOLDSTONE_WARD.md) — Ward identity and
  controlled diagnostic policies.
- [TDDFT_GF_INTEGRATION.md](TDDFT_GF_INTEGRATION.md) and
  [TDDFT_REALSPACE_GF.md](TDDFT_REALSPACE_GF.md) — backend-specific numerical
  details.
- [TDDFT_SPECTRAL_ANALYSIS.md](TDDFT_SPECTRAL_ANALYSIS.md) — loss matrix,
  Stoner continuum, branch tracking, and linewidth gates.
- [TDDFT-14 performance report](dev/TDDFT-14_PERFORMANCE_MPI_REPORT.md) —
  measured backend regimes and MPI reuse strategy.
- [TDDFT_RELEASE_STATUS.md](TDDFT_RELEASE_STATUS.md) — current release
  decision and open gates.
