# GBT WP5 / Gate G5 report

Date: 2026-08-03.

## Outcome

WP5 implementation cleanup is complete. Both solvers consume the same `ee`
operator built from the linked primitive `S` blocks, and reciprocal assembly
uses only `fourier_transform_array`.

**Gate G5: PASS (accepted/open with an operational bound).** The shared-bond
and pre-eigensolver Hermiticity requirements pass. The strict unconstrained
monotonic RS/k-space sequence was not demonstrated, but the project owner has
accepted `lld <= 28` as the meaningful recursion range for the approximately
1000-atom PBC cluster; larger depths are not used as gate evidence. WP6 may
proceed with this limitation recorded.

## Deletion inventory

Deleted from `source/reciprocal_fourier.f90`:

- the `fourier_transform_gbt` implementation;
- the `fourier_transform_gbt_array` implementation, including its `h0`/`bz`,
  cone-angle, and `exp(+-i q.R)` reconstruction;
- all GBT-specific transforms of `ee`, `eeo`, and `eecc` in first- and
  second-order reciprocal assembly;
- every reciprocal `gbt_kspace` conditional.

Deleted from `source/reciprocal.f90`:

- both type-bound procedure bindings;
- both module interfaces.

Deleted outside reciprocal assembly:

- all three `gbt_kspace` physics conditionals in `hamiltonian_ccor.f90`;
- redundant frozen-magnon assignments to `hamiltonian%gbt_kspace`.

`gbt_kspace` remains only in the Hamiltonian input declaration, stored
compatibility field, default initialization, and input-boundary warning. If
true, it is warned about, ignored, and immediately cleared; users must select
`magnetic_representation='gbt_single_q'` explicitly. It cannot change the
magnetic representation, reciprocal route, or CCOR physics.

`tests/unit/test_gbt_wp5_source_contract.py` prevents restoration of either
legacy transform or use of `gbt_kspace` outside that input boundary.

## Cache invalidation contract

`hamiltonian%operator_generation` advances at the common `build_bulkham`
boundary before rebuilding the real-space operator. This single generation
covers every consumed q, cone angle, sublattice reference axis, and potential
parameter without a collision-prone floating-point fingerprint.

`reciprocal%cached_operator_generation` records the generation used for H(k).
On mismatch, `invalidate_if_operator_changed` discards:

- `hk_bulk`, `hk_so`, `hk_total`, and `sk_overlap`;
- mesh and path eigenvalues/eigenvectors;
- total/projected DOS and integrated DOS;
- band moments and Cartesian spin-density projections;
- tetrahedron-derived storage and canonical N/EBAND validity.

The k points and weights remain reusable because they do not depend on q,
cone/reference axes, or potential parameters. `build_kspace_hamiltonian`
invalidates unconditionally before rebuilding and records the new generation;
the eigensolver, DOS, and reciprocal density entry points also check the
generation. Real-space block/Lanczos objects are rebuilt per recursion call;
the persistent fast Chebyshev cache is explicitly reset at every Chebyshev
recursion entry.

`UnitKspaceOccupations` independently changes q, cone, reference-axis, and
potential rebuild generations. For every change it verifies that H(k), both
eigensystems, DOS, moments, and density projections are deallocated while the
k weights survive. All four cases pass.

## Hermiticity evidence

`reciprocal%diagonalize_hamiltonian` retains fatal assertions immediately before
`zheev`/`zhegv` for both H(k) and O(k), with tolerance
`1e-10 * max(1,max|matrix|)`. Across the finite-q meshes and the independent
Green-function comparison, the largest observed value was:

```text
max |H-H^H| = 2.7767e-16
```

No generalized-overlap GBT case was enabled; that combination remains guarded
for WP6.

## Identical-bond and convergence comparison

The fixtures in `tests/gbt_wp5/cases.json` use the same bcc-Fe restart,
`q=(0,0,+/-0.05)`, `theta=20 degrees`, first-order S-only Hamiltonian, and full
chemical BZ. Only solver resolution changes. All final `gbt_bonds.out` files
for k meshes 4^3, 6^3, 8^3, 12^3, 16^3 and block-recursion depths 12, 20, 28,
36, 44 are byte-identical:

```text
SHA-256 7044ac05ae1cb3d648d2cde6b5055b8d1660b8165bb160a0aaea5312bc607474
```

The raw fixed-potential band energies expose the practical recursion-depth
limit. `Delta E(q)` is `EBAND(q=+0.05)-EBAND(q=0)` in Ry.

| Route / resolution | EBAND(q=0) (Ry) | Delta E(+0.05) (Ry) |
| --- | ---: | ---: |
| k 4^3 | -1.99686650 | +2.031e-5 |
| k 6^3 | -1.99993633 | +2.8338e-4 |
| k 8^3 | -1.99532044 | -1.02336e-3 |
| k 12^3 | -1.99801565 | -3.99e-6 |
| k 16^3 | -1.99732106 | -4.0563e-4 |
| RS block lld=12 | -1.99622940 | -1.843e-5 |
| RS block lld=20 | -1.99762396 | +2.5549e-4 |
| RS block lld=28 | -1.99829850 | +3.6348e-4 |
| RS block lld=36 | -1.94632551 | -4.084591e-2 |
| RS block lld=44 | -1.98738532 | +1.455091e-2 |

The k sequence oscillates for this metallic Fermi surface, and the block
recursion sequence loses stability beyond lld=28 under the present terminator,
energy-integration settings, and approximately 1000-atom PBC cluster. Pairing
all increasingly fine entries therefore does not give monotonic route
differences. The project owner accepts `lld <= 28` as the useful range and has
left G5 open. The remaining oscillation is numerical convergence debt, not an
operator-routing discrepancy: the byte-identical bonds and roundoff
Hermiticity isolate it downstream of assembly.

An additional same-run `kspace_green` comparison at finite q, lld=20, 8^3,
and eta=0.02 Ry also consumed the shared operator. It passed Hermiticity
(`2.2377e-16`) and Lehmann/Dyson equivalence (`1.89929e-13`), but the recursion
and Lehmann DOS were not converged to each other (`rms=2.582917`, maximum
`16.41406`), consistent with the documented unequal broadening/terminator
behavior.

## Commands and tests

```bash
cmake --build build_13 -j2
ctest --test-dir build_13 --output-on-failure -L unit
build_13/bin/UnitKspaceOccupations
python3 tests/unit/test_gbt_wp5_source_contract.py

/Users/andersb/envs/p311/bin/python tests/run_test.py \
  --binary build_13/bin/rslmto.x \
  --cases-json tests/gbt_wp5/cases.json \
  --case-name <wp5 case> \
  --scratch-root /private/tmp/rslmto_wp5_g5_fixed --force-serial
```

Results:

- build: PASS;
- unit suite: 17/17 PASS;
- WP5 source contract: PASS;
- q/cone/reference-axis/potential cache invalidation: PASS;
- all ten frozen-magnon resolution cases: PASS as executions;
- finite-q same-run Green comparison: PASS as execution;
- identical bond dump: PASS;
- pre-solver Hermiticity: PASS;
- unconstrained monotonic RS/k-space convergence: not demonstrated;
- operational recursion range: `lld <= 28` for this cluster;
- **G5: PASS (accepted/open by project-owner decision)**.
