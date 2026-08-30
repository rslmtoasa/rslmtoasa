# GBT re-audit WP02: gauge algebra, Hermiticity, and shifted k

## Status

PASS for the scalar-relativistic fixed-potential WP02 operator gate.
The current primitive structure-constant gauge satisfies link unitarity, bond
reversal, closed-loop purity, contracted real-space Hermiticity, reciprocal
Hermiticity, and the exact one-sublattice shifted-momentum identity.  No
production GBT, SCF, or energy-bookkeeping source was changed.

This report establishes the tested algebraic scope only.  It does not establish
HOH/CCOR covariance, commensurate-supercell folding, rotating-frame density or
constraint semantics, frozen-magnon energy semantics, LKAG consistency, or
fully self-consistent GBT.

## Dependency and implementation

WP01 passes at the audited base commit `0fb13b6c` (`test: enforce exact q-zero
reduction of GBT`).  Its direct unit gate was rerun before this WP02 test was
added.

The WP02 addition is:

- `tests/unit/test_gbt_wp02_gauge_shifted_k.f90`: deterministic fixed-potential
  unit/functional gate;
- CTest target `UnitGbtWp02GaugeShiftedK`.

The fixture has one site, two scalar orbitals, an onsite block and two opposite
pairs of complex directed hopping blocks.  The ordinary reference is built
independently as two scalar spin-sector Fourier sums.  The GBT candidate is
built with the production `gbt_endpoint_link`, `gbt_lift_orbital_block`, and
`gbt_contract_collinear` path.  The fixture is fixed-potential: all orbital
blocks and local potential factors are deterministic constants, with no SCF or
energy integration.

## Convention derivation

The repository's conventions are fixed by the implementation, rather than by
assigning a sign from an external formula:

1. `gbt_bond_phase(q,d,alat)` returns
   `alpha = 2*pi*q.d/alat`.
2. For common collinear endpoint frames,
   `gbt_endpoint_link` is
   `diag(exp(-i*alpha/2), exp(+i*alpha/2))`, in the repository's
   `(up, down)` block ordering.
3. Reciprocal assembly in `reciprocal_fourier.f90` uses
   `exp(+i*2*pi*k.R)`.
4. For a lattice translation `d = R*alat`, the spin blocks therefore carry

   ```text
   up:   exp(+i*2*pi*k.R) exp(-i*pi*q.R) = exp(+i*2*pi*(k-q/2).R)
   down: exp(+i*2*pi*k.R) exp(+i*pi*q.R) = exp(+i*2*pi*(k+q/2).R)
   ```

Consequently, the exact WP02 oracle is

```text
H_GBT(k,q,theta=0) = H_up(k-q/2) direct-sum H_down(k+q/2).
```

The test compares both complete matrices and sorted Hermitian eigenvalue
multisets at three generic `(k,q)` pairs.  The independent ordinary side never
calls `gbt_endpoint_link`.

## Numerical results

The threshold was selected before running the test: `1e-12` for absolute,
relative, and eigenvalue residuals; negative controls must exceed `1e-3`.
Compiler and build details are recorded in §4.

| check | max absolute | max relative | max eigenvalue |
|---|---:|---:|---:|
| endpoint-link unitarity, one/multi-sublattice frames | `1.1102e-16` | `1.1102e-16` | `0.0000e+00` |
| endpoint-link independent SU(2) oracle | `0.0000e+00` | `0.0000e+00` | `0.0000e+00` |
| bond reversal, one/multi-sublattice frames | `1.5701e-16` | `1.5701e-16` | `0.0000e+00` |
| one-sublattice translated closed loop | `0.0000e+00` | `0.0000e+00` | `0.0000e+00` |
| multi-sublattice translated closed loop | `2.3066e-16` | `2.3066e-16` | `0.0000e+00` |
| contracted real-space bond reversal | `5.5943e-17` | `5.5943e-17` | `0.0000e+00` |
| reciprocal `H(k)` Hermiticity, three generic points | `1.7554e-16` | `1.7554e-16` | `0.0000e+00` |
| shifted-k matrix identity, three generic points | `3.3307e-16` | `1.5895e-16` | `2.2204e-15` |
| centrosymmetric `spec H(k,q) = spec H(-k,-q)` | — | — | `0.0000e+00` |

The test-only negative controls are sensitive by large margins:

| deliberate test perturbation | residual |
|---|---:|
| missing half phase | `9.7685e-01` |
| q sign flip | `1.7336e+00` |
| reversed endpoint order | `5.3062e-01` |

The q-reversal assertion is intentionally the appropriate centrosymmetric
pointwise relation, `spec H(k,q) = spec H(-k,-q)`.  At a generic fixed `k`,
`spec H(k,q) = spec H(k,-q)` is not the required identity for a spin-split
collinear reference.

## External cross-check

The result is consistent with the two relevant external formulations:

- The public [ELK `genwfsvp_sp.f90` source](https://elk.sourceforge.io/doxygen/html/genwfsvp__sp_8f90_source.html)
  explicitly adds `+q/2` to the first variational spin component and `-q/2`
  to the second when `spinsprl` is enabled.  ELK's [magnetism material](https://elk.sourceforge.io/CECAM/Nordstrom-magnetism.pdf)
  writes the equivalent generalized Bloch spinors with the two components
  shifted by opposite half wavevectors.  Which component is called up, and the
  sign assigned to `q`, is convention-dependent; the repository's derivation
  above fixes its assignment as `(up: k-q/2, down: k+q/2)`.  Reversing q maps
  the two assignments without changing the mathematical content.
- In the LMTO-ASA formulation of [Shallcross, Nordström, and Sharma, Phys.
  Rev. B 76, 054444 (2007)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.76.054444),
  the local-frame potential is spin diagonal while the spiral structure matrix
  is formed from `S(k+q/2)` and `S(k-q/2)`.  At `theta=0` its spin-mixing terms
  vanish and the same two shifted spin sectors remain.  WP02 tests precisely
  this theta-zero limit using the repository's real-space LMTO link and its
  ordinary downstream Fourier phase.

## Reproduction

Environment and scope:

- base commit tested: `0fb13b6c`;
- compiler: GNU Fortran 13.3.0 (`/bin/gfortran`);
- build type: Release;
- MPI: OFF; OpenMP: ON; libXC: OFF; CUDA: OFF;
- feature set: scalar-relativistic fixed-potential fixture; SOC, Hubbard,
  HOH/CCOR, local-axis, SCF, and total-energy paths are not exercised;
- q and k: dimensionless fractional reciprocal coordinates in the repository's
  Cartesian `2*pi/alat` convention; the synthetic bond vectors are in lattice
  translations and `alat = 3.4`;
- cone angle: `theta=0` for the shifted-k and q-reversal tests; arbitrary fixed
  endpoint reference frames are used for link/loop/reversal algebra;
- k mesh: three deterministic generic points, not a Brillouin-zone quadrature.

Commands:

```text
cmake -S . -B build-ci-reference-serial -DRUN_UNIT_TESTS=ON -DRUN_REG_TESTS=ON
cmake --build build-ci-reference-serial --target UnitGbtWp02GaugeShiftedK -j2
ctest --test-dir build-ci-reference-serial \
  -R '^(UnitGbtWp01Q0|UnitGbtWp02GaugeShiftedK)$' --output-on-failure
```

Observed result: both focused tests passed.  The WP02 executable reported an
overall maximum absolute residual of `3.3307e-16`, relative residual of
`2.3066e-16`, and eigenvalue residual of `2.2204e-15`.

## Completion checklist

- [x] Link unitarity tested.
- [x] Bond reversal tested for one- and multi-sublattice endpoint frames.
- [x] Closed-loop pure-gauge identity tested with nontrivial translations.
- [x] Contracted real-space bond Hermiticity tested.
- [x] Reciprocal Hamiltonian Hermiticity tested.
- [x] Fourier/q sign convention derived explicitly.
- [x] `k-q/2` / `k+q/2` identity passes at generic points.
- [x] q↔-q spectral symmetry checked in the appropriate inversion-related form.
- [x] Meaningful half-phase, sign, and endpoint-order negative controls included.
- [x] ELK and LMTO correspondence documented.
- [x] No SCF or energy-fitting changes made.
