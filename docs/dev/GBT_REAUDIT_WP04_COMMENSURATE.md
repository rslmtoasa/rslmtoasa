# GBT re-audit WP04: matched commensurate-supercell oracle

## Scope and verdict

WP04 adds an exact fixed-potential matched-operator oracle for the scalar-
relativistic, one-sublattice GBT core. Period-2 and period-3 commensurate
cases compare the complete supercell eigenvalue multiset with the correctly
folded primitive GBT spectrum. The test is intentionally synthetic and uses a
small LMTO-like orbital block table; it does not claim production bcc-Fe or
SCF equivalence.

The implementation is split into the reusable constructor
`tests/unit/gbt_wp04_matched_operator.py` and the executable regression
`tests/unit/test_gbt_wp04_commensurate_supercell.py`. Both representations
receive the same immutable five-record directed table: one onsite block and
two explicitly stored reverse-bond pairs. No neighbor search, cutoff, or
empirical energy shift is used on the supercell side.

## Construction and folding derivation

The repository convention is q in Cartesian units of `2*pi/alat`, with
directed primitive cell displacement `d` and ordinary Fourier factor
`exp(+2*pi*i*k.d)`. For each table row the primitive side constructs

```text
G_ab(d) = U_a^H exp(-i*pi*q.d*sigma_z) U_b
H_GBT,ab(d) = W_a [S_ab(d) tensor G_ab(d)] W_b.
```

The explicit side makes the same finite operator in the lab frame. With

```text
P_na = exp(-i*pi*q.n*sigma_z) U_a
W_na = P_na W_a P_na^H,
```

each translated directed row contributes

```text
W_na [S_ab(d) tensor I] W_mb,
```

where `m=n+d` before reduction into the N-cell basis. Onsite local terms are
rotated by `P_na` and added once per site. A supercell boundary contributes
the ordinary phase `exp(+2*pi*i*k_super.d)`.

For `q_x=1/N`, the spinor frame obeys `P_(n+N,a)=-P_(n,a)`, even though the
physical rotated potential is N-periodic. Consequently the unitary map from
folded primitive GBT blocks to the explicit supercell uses

```text
k_l = k_super + (l + 1/2)/N * xhat,  l=0,...,N-1,
Q_(n,l) = P_n exp(+2*pi*i*(l+1/2)*n/N) / sqrt(N).
```

The test checks the direct matrix identity

```text
H_SC Q = Q [ H_GBT(k_0) direct-sum ... direct-sum H_GBT(k_(N-1)) ],
```

as well as `Q^H Q=I`, sorted complete eigenvalue multisets, trace, Frobenius
norm, and pre-diagonalization Hermiticity. The wrong fold without the
spinor half-period shift is a negative control.

The utility retains three-component cell and physical displacements and
validates `physical_displacement = alat * cell_displacement` for this
one-cell fixture. A multi-sublattice extension is supported by the
`source`/`target` fields and per-sublattice frames. For a basis with nonzero
`tau_b-tau_a`, the production physical-displacement convention requires that
the basis phase be included consistently in the endpoint frame or in the
translation coordinate; that extension is documented but intentionally not
accepted by this one-dimensional q-period fixture.

## Evidence

The threshold was declared before evaluating results: exact matrix, unitary,
invariant, Hermiticity, and spectrum residuals must be below `2e-12`; the
wrong-fold negative control must exceed `1e-3`. The test was run against
commit `e51ab219d26f95e4ac627b0cb7fb85b08a64fe62` with the WP04 working-tree
changes applied, using Python `3.12.3` and NumPy `2.5.1`. The isolated CMake
validation build used GNU Fortran `13.3.0`, Release mode, and
`ENABLE_MPI=OFF`, `ENABLE_OPENMP=OFF`, `ENABLE_CUDA_PLUGIN=OFF`,
`ENABLE_CUDA_RECIPROCAL=OFF`, `ENABLE_MKL_KERNELS=OFF`, and
`ENABLE_SPGLIB=OFF`; the GBT feature fixture is scalar-relativistic and
fixed-potential with SOC, SCF, and generalized-overlap paths absent.

The fixture uses `alat=3.7`, a generic reference cone
`theta=0.71`, `phi=0.37`, and primitive-reciprocal supercell k points
`k_super=(0.019,0,0)` for period 2 and `(-0.021,0,0)` for period 3. The
periods use `q=(0.5,0,0)` and `q=(1/3,0,0)`, respectively. The exact residuals
are:

| case | matrix | spectrum | trace | Frobenius | unitary | Hermiticity | wrong fold |
|---|---:|---:|---:|---:|---:|---:|---:|
| period 2 | `4.635e-16` | `5.551e-16` | `2.224e-16` | `2.220e-16` | `1.346e-16` | `1.178e-16` | `5.134e-01` |
| period 3 | `5.241e-16` | `8.882e-16` | `4.523e-16` | `4.441e-16` | `4.448e-16` | `3.942e-16` | `5.029e-01` |

All exact residuals pass the declared threshold and both negative controls
are detected. The failure class exercised by this oracle is
representation/folding-level; no SCF, quadrature, or energy error is involved.
The registered CTest entry and focused existing GBT targets
`UnitGbtStructure`, `UnitGbtWp01Q0`, and `UnitGbtWp02GaugeShiftedK` also passed
in that isolated build.

## Reproduction

```text
python3 tests/unit/test_gbt_wp04_commensurate_supercell.py
ctest --test-dir <configured-build> --output-on-failure \
  -R UnitGbtWp04CommensurateSupercell
```

The test is a lean Python/NumPy unit fixture and does not require a production
SCF input or compiler-dependent executable. The production VAL-16 test remains
a separate diagnostic: its historically generated finite clusters do not use
this matched operator and therefore are not the WP04 acceptance oracle.

## Evidence boundary

This WP04 result establishes exact fixed-potential spectral equivalence for
the finite synthetic primitive operator at periods 2 and 3, including a
generic nonzero cone angle. It does not establish real LMTO structure-constant
production equivalence, HOH/CCOR beyond WP03's separate covariance scope,
rotating-frame density semantics, constraint fields, energy bookkeeping,
corrected frozen magnons, or self-consistent GBT.

## Completion checklist

- [x] One canonical primitive directed operator feeds both representations.
- [x] Period-2 commensurate case is automated.
- [x] Period-3 commensurate case is automated.
- [x] Folding map and spinor half-period are documented.
- [x] Complete eigenvalue multisets are compared.
- [x] Trace, Frobenius, Hermiticity, and unitary diagnostics are reported.
- [x] No empirical energy shift is used.
- [x] No independent neighbor or cutoff generation is used.
- [x] A generic cone angle is tested.
- [x] Multi-sublattice extension is considered and documented.
- [ ] Higher-level production fixture: deferred; the existing VAL-16 inputs
      use unmatched finite supports and are not a valid WP04 oracle.
