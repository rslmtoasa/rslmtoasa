# GBT re-audit WP01: exact q=0 reduction

## Status

PASS for the WP01 scalar-relativistic fixed-potential q=0 operator gate.
No production Hamiltonian source change was required.  The only implementation
correction made while writing the test was in the new test's reciprocal
assembly helper: it now selects the directed bond's site subblock before
Fourier accumulation.

This result establishes the requested q=0 reduction for the tested scope.  It
does not establish production readiness for nonzero q, SOC, HOH, CCOR,
Hubbard terms, local-axis mode, constraints, SCF, LKAG, or the full validation
ladder in the master blueprint.

## Dependency and implementation

WP00 is available as
`docs/dev/GBT_REAUDIT_WP00_CURRENT_STATE.md` at the audited baseline
`accd1706` (`docs: map current GBT implementation and validation evidence`).
The current architecture identified there is used unchanged: raw directed
structure constants are linked at the endpoints, lifted to spinor blocks, and
contracted through the ordinary LMTO potential path.

The WP01 additions are:

- `tests/unit/test_gbt_wp01_q0.f90`: independent deterministic operator gate;
- `tests/gbt_wp01/cases/bccFe`: committed fixed-potential bcc-Fe fixture;
- `tests/validation/val18_gbt_wp01_q0.py`: production `kspace_green` check
  comparing ordinary `periodic_nc` and `gbt_single_q` at q=0;
- CTest registrations `UnitGbtWp01Q0` and `Val18GbtWp01Q0BccFe`.

The ordinary reference in the unit test is built with
`lmto_bond_value`/`lmto_hhmag_to_spinor`; GBT uses the production primitive
link/lift/contract helpers.  The expected frame transformations are formed by
an explicit test-side SU(2) construction, not by calling the GBT endpoint-link
helper on the reference side.

## Fixture and conventions

The synthetic operator fixture has two sites, two scalar orbitals per site,
four directed bonds (two onsite and a forward/reverse pair), complex
Hermitian-compatible raw blocks, distinct potential factors, and three generic
k points.  It checks the ordinary primitive-cell operator before any
diagonalization.

The realistic fixture is scalar-relativistic bcc Fe with `nsp = 3`, `hoh =
.false.`, `strux_backend = 'strux_lib'`, `strux_want_sdot = .false.`, `nstep =
0`, and the committed Fe fixed potential.  The production runner uses
`q_ss = (0, 0, 0)` in Cartesian `2*pi/alat` units, `theta_ss = 0`, no symmetry
reduction, no time reversal, `green_eta = 0.02 Ry`, and an `8 x 8 x 8`
reciprocal mesh.  The fixture's `kspace_green` route emits 5,890 numeric
records across `kspace_green_c1.dat` and `kspace_green_gf.dat`.

The unit-only rotation cases are:

- global SU(2): `theta = 37 degrees`, `phi = 0.41`;
- two-sublattice: `(theta, phi) = (0.63, 0.27)` and `(1.11, -0.58)`;
- negative control: source endpoint of one directed bond perturbed by `0.17`
  radians.

## Exact operator results

All values are maxima over the tested directed bonds or three generic k
points.  The acceptance threshold is `1e-12` for both absolute and relative
operator residuals; eigenvalue residuals are reported separately.

| check | max absolute | max relative | max eigenvalue |
|---|---:|---:|---:|
| q=0 collinear directed blocks | 2.2204e-16 | 1.1484e-16 | 0.0000e+00 |
| q=0 collinear onsite blocks | 2.2204e-16 | 1.6489e-16 | 0.0000e+00 |
| q=0 collinear reciprocal `H(k)` | 2.2204e-16 | 1.1794e-16 | 1.7764e-15 |
| q=0 global SU(2) transformed blocks | 2.2248e-16 | 2.4594e-16 | 0.0000e+00 |
| q=0 ordinary global rotation | 1.1189e-16 | 1.2735e-16 | 0.0000e+00 |
| q=0 global SU(2) reciprocal `H(k)` | 2.2888e-16 | 1.6837e-16 | 1.3323e-15 |
| q=0 two-sublattice transformed blocks | 2.2253e-16 | 2.5131e-16 | 0.0000e+00 |
| q=0 two-sublattice onsite blocks | 2.2253e-16 | 2.5131e-16 | 0.0000e+00 |
| q=0 two-sublattice reciprocal `H(k)` | 2.7756e-16 | 2.3457e-16 | 1.7764e-15 |

The deliberate negative-control residual is `8.8014e-02`, well above its
`1e-3` detection threshold.

## Production bcc-Fe result

The ordinary and GBT fixed-potential runs agree on all 5,890 compared numeric
records:

| comparison | max absolute | max relative |
|---|---:|---:|
| `periodic_nc` vs `gbt_single_q`, q=0 | 1.110223e-14 | 5.551115e-15 |
| GBT endpoint-link vs identity | 0.000000e+00 | — |

The production result is an application-level Green-function comparison from
the same q=0-built operator.  The direct pre-diagonalization block and `H(k)`
acceptance metrics are owned by `UnitGbtWp01Q0`.

## Reproduction

Environment used:

- commit baseline: `accd17061fb06536a731b3994e5658186a29539a`;
- compiler: GNU Fortran 13.3.0 (`/bin/gfortran`);
- build type: Release;
- MPI: OFF; OpenMP: ON; libXC: OFF; CUDA: OFF;
- Fortran flags: `-ffree-line-length-0 -cpp -mtune=generic` plus the
  configured Release/OpenMP flags.

Commands:

```text
cmake -S . -B build-ci-reference-serial -DRUN_UNIT_TESTS=ON -DRUN_REG_TESTS=ON
cmake --build build-ci-reference-serial --target UnitGbtWp01Q0 -j2
ctest --test-dir build-ci-reference-serial \
  -R '^(UnitGbtWp01Q0|Val18GbtWp01Q0BccFe)$' --output-on-failure
```

Observed result: both CTest entries passed; total CTest time was 18.40 s on
the stated serial build.
