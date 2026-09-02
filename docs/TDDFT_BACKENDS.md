# TDDFT chi0 backends

For the user-facing backend choice, convergence controls, support matrix, and
reproducible material decks, start with
[TDDFT_USER_GUIDE.md](TDDFT_USER_GUIDE.md). The current release decision is in
[TDDFT_RELEASE_STATUS.md](TDDFT_RELEASE_STATUS.md).

The TD-DFT response layer consumes a common `tddft_chi0_batch_result`.  Each
entry in `q_response(:)` is the compatibility-preserving
`tddft_chi0_result` used by the Dyson, loss-matrix, and text-output code.  A
batch result also carries the requested `(3,nq)` q grid, the frequency batch,
and aggregate provenance metadata.

The public Fortran contract is in `source/tddft_backend.f90`:

```fortran
type(tddft_chi0_request)
   real(rp), allocatable :: q_points(:, :)  ! (3,nq)
   integer, allocatable :: q_indices(:)     ! optional endpoint-data indices
   real(rp), allocatable :: omega(:)
end type

type(tddft_chi0_batch_result)
   type(tddft_chi0_result), allocatable :: q_response(:)
   real(rp), allocatable :: q_points(:, :)
   integer, allocatable :: q_indices(:)     ! returned endpoint/order provenance
   real(rp), allocatable :: omega(:)
end type

type, abstract :: tddft_chi0_backend
contains
   procedure :: evaluate(request, batch_result)
   procedure :: evaluate_one(q, omega, result)
   procedure :: evaluate_frequency_batch(q, omega(:), result)
   procedure :: evaluate_q_batch(q_points(:, :), omega, batch_result)
   procedure :: evaluate_grid(q_points(:, :), omega(:), batch_result)
   procedure :: evaluate_static_grid(q_points(:, :), batch_result)
end type
```

`evaluate` is the one deferred operation.  The convenience operations build a
request and call it, so every backend has the same scalar and batch semantics.
The q-batch result is a list of per-q response matrices rather than a new
four-dimensional replacement for `tddft_chi0_result`; this leaves the existing
Dyson and output interfaces unchanged.

## Backend names

| Requested name | Canonical name | Contract implementation |
| --- | --- | --- |
| `eigenpairs` | `eigenpairs` | Explicit transition reference backend |
| `kspace_lehmann` | `kspace_lehmann` | K-space Lehmann GF bubble |
| `green`, `kspace_gf` | `kspace_lehmann` | Compatibility aliases |
| `realspace_gf` | `realspace_gf` | Native real-space provider adapter |
| `realspace`, `rs_gf` | `realspace_gf` | Compatibility aliases |
| `mock` | `mock` | Deterministic unit-test fixture |

`make_tddft_chi0_backend` allocates the selected contract.  Initialization is
separate because the calculation owns the response mesh and endpoint data.
The current production calculation attaches the periodic eigenpair endpoint
adapter for `eigenpairs`, the K-space Lehmann adapter for `green` and
`kspace_lehmann`, and the native RS-LMTO provider for `realspace_gf`. The
native route cannot fall back to a k-space calculation. Its current source is
direct real-axis and supports both dynamic bare-`chi0` and the separate exact
static Ward operation; mixed contour and Dyson enhancement remain explicitly
rejected.

## Provenance and ownership

Every returned result records the canonical backend, implementation,
integration convention, broadening, k-mesh or real-space controls, q/omega
batch sizes, convergence status, and whether a Goldstone correction was
applied.  “Not assessed” is reported when a backend cannot establish numerical
convergence; a provider must not claim convergence merely because evaluation
completed.

Endpoint eigenpair arrays are copied into the eigenpair adapters during
initialization.  Provider-backed real-space objects use polymorphic allocatable
ownership and are copied with the existing Fortran allocate-from-source
pattern.  No backend returns borrowed inner-loop work arrays.

## Evaluation sequence

```text
caller
  |
  | make_tddft_chi0_backend(name)
  v
abstract chi0 backend
  |
  | evaluate(request: q_points, q_indices, omega)
  +------------------+--------------------+---------------------+
  |                  |                    |
eigenpairs     K-space Lehmann       native RS-GF provider
  |                  |                    |
transition     G(k,z) bubble          chi0(R,omega), then
sum             (same response        requested q/mixed/real
                basis)                 representation
  +------------------+--------------------+---------------------+
                         |
              tddft_chi0_batch_result
                         |
                 common Dyson/loss layer
```

The native real-space branch constructs and reuses `chi0(R,omega)` inside its
provider.  The common interface does not request, require, or hide a
`G(R,z)->G(k,z)` transformation.
