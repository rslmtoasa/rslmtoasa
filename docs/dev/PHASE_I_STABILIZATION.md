# Phase-I stabilization: architecture and status

This is the developer-facing record for Phase I — structural stabilization of
RS-LMTO-ASA.  It describes the code as it is organized at the current STAB
head.  It is intentionally not a feature roadmap: use the code and the CTest
configuration as the source of truth when this document becomes stale.

For detailed workflow entry points and kernels, see
[`docs/DEVELOPER_MAP.md`](../DEVELOPER_MAP.md).

## Core architecture

RS-LMTO-ASA is primarily a real-space LMTO / Green-function magnetic
electronic-structure code.  The common physical stack is:

```text
lattice / charge / potential
          |
        H(R)
       /    \
      RS   reciprocal
      |        |
  recursion  eigensystem
      |       /       \
      GF  Lehmann GF  response
      |
magnetic observables / transport
```

`H(R)` is the real-space LMTO Hamiltonian.  The real-space branch uses
recursion and Green functions for SCF and magnetic observables.  The
reciprocal branch Fourier-transforms the same real-space Hamiltonian and uses
the resulting eigensystem for k-space SCF, bands, DOS, spectral quantities,
and the reciprocal Green-function/response work.  It complements the
real-space/Green-function architecture; it does not replace or supersede it.

## Current module and submodule map

The parent modules below hold the public derived types, their type-bound
procedure declarations, and the contracts used by callers.  The named
submodules contain the corresponding implementation partitions.  A file not
listed as a submodule here is either a parent contract file or a standalone
support/kernel module.

| Area | Public contract / parent | Implementation submodules and principal helpers |
|---|---|---|
| Hamiltonian | `source/hamiltonian.f90` (`hamiltonian_mod`) | `hamiltonian_build.f90`, `hamiltonian_ccor.f90`, `hamiltonian_hubbard.f90`, `hamiltonian_paoflow_io.f90` |
| Recursion | `source/recursion.f90` (`recursion_mod`) | `recursion_core.f90`, `recursion_haydock.f90`, `recursion_chebyshev.f90`, `recursion_transport.f90`; hot kernels remain in `haydock_fast.f90` and `chebyshev_fast.f90` |
| Calculation / workflow dispatch | `source/calculation.f90` (`calculation_mod`) | `calculation_preprocessing.f90`, `calculation_reciprocal.f90` |
| Green functions | `source/green.f90` (`green_mod`) | `green_lifecycle.f90`, `green_lanczos.f90`, `green_block.f90`, `green_chebyshev.f90` |
| SCF and atomic stack | `source/self.f90` (`self_mod`) | `self_reciprocal.f90`, `self_xc_response.f90`; supporting physical state is in `potential.f90`, `charge.f90`, and related lattice/atomic files |
| Charge / electrostatics | `source/charge.f90` (`charge_mod`) | `charge_interface.f90`, `charge_madelung_2d.f90` |
| Exchange | `source/exchange.f90` (`exchange_mod`) | `exchange_dynamics.f90` (damping and inertia); `conductivity.f90` remains its own transport module |
| Reciprocal space | `source/reciprocal.f90` (`reciprocal_mod`) | `reciprocal_lifecycle.f90`, `reciprocal_backend.f90`, `reciprocal_fourier.f90`, `reciprocal_bands.f90`, `reciprocal_dos.f90`, `reciprocal_projection.f90`, `reciprocal_occupations.f90`, `reciprocal_spin_density.f90`, `reciprocal_moments.f90`, `reciprocal_green.f90`, `reciprocal_bsf.f90` |

When extending an area, add the public type-bound declaration and interface
contract to its parent module, and put a sizeable implementation partition in
the appropriate submodule.  Do not expose a new internal helper merely to
avoid a submodule boundary.

## Protected atomic LMTO core

The mature atomic, radial, and potential-building routines in `self.f90` and
the closely related atomic/potential files are protected scientific
infrastructure.  Their numerical contracts are established by accumulated
physics use, not by their age or formatting.

Do not casually modify these routines for cleanup, modernization, abstraction,
or test convenience.  **Refactor around them before refactoring through
them.**  A change that truly must enter this region needs a stated numerical
contract, a pre-change baseline, and the smallest possible patch.

## Test taxonomy

CTest labels are configuration-dependent.  The following labels exist in the
Phase-I gate configuration and are useful entry points:

```sh
ctest -L unit --output-on-failure
ctest -L quick --output-on-failure
ctest -L rs --output-on-failure
ctest -L kspace --output-on-failure
ctest -L green --output-on-failure
ctest -L exchange --output-on-failure
ctest -L conductivity --output-on-failure
ctest -L functional --output-on-failure
ctest -L tooling --output-on-failure
ctest -L performance --output-on-failure
ctest -L validation --output-on-failure
```

- `unit` tests small numerical, lifecycle, and source-contract invariants.
- `functional` runs production workflows with committed reference checks;
  `scf`, `postproc`, `rs`, `kspace`, `exchange`, and `conductivity` narrow it
  to a route or observable. `quick` is the selected lean functional subset.
- `tooling` checks manifest, dispatch, and campaign tooling. It is not a
  substitute for a numerical workflow test.
- `performance` contains profiling/performance checks. Passing it is not a
  correctness or maturity claim.
- `validation` is reserved for the configured validation campaign and should
  be read with that campaign's scope and tolerances.

There are currently no `surface`, `interface`, or `vacuum` CTest labels in
this configuration.  Select those fixtures by their actual test names (for
example `ctest -R 'Example_(surface|interface)_' --output-on-failure` and
`ctest -R 'UnitVacuum' --output-on-failure`) or add a label deliberately in a
test-maintenance task.  Do not document nonexistent label commands.

## Fixture vocabulary

Fixtures are compact coverage instruments, not interchangeable materials
benchmarks:

| Fixture family | What it is intended to cover |
|---|---|
| Si / `sp` | lean nonmagnetic, Chebyshev, real-space/k-space equivalence, and basis coverage |
| bcc Fe / `spd` | magnetic Block/Lanczos recursion and magnetic reciprocal paths |
| Metallic multisublattice | multiple site and species behaviour, including magnetic multisublattice stacks |
| Cu | surface, impurity, interface, and vacuum-region wiring |
| Pt | SOC and transport/exchange-conductivity paths where the fixture enables them |

Passing a fixture establishes only the contract it checks.  In particular, do
not describe every fixture as a converged materials benchmark.

## Structure-constant policy

`strux_lib` is the normal structure-constant backend for the test suite.
Explicit `legacy_strux` tests remain to preserve and compare the legacy
backend.  This test-fixture policy does not change the application's
production default.

## Feature maturity

The labels below are deliberately conservative.  “Validated” means that the
current suite has direct, relevant checks; it is not a claim of broad material
benchmarking.  Compilation or an unasserted smoke run alone is insufficient.

| Capability | Status | Evidence at this head |
|---|---|---|
| Atomic/radial/potential LMTO core | Production | Protected established SCF infrastructure; retained behind its existing numerical contracts |
| Real-space SCF and recursion (Block, Lanczos, Chebyshev) | Validated | Regression/functional fixtures and the `quick` suite exercise representative magnetic and nonmagnetic paths |
| Real-space Green-function lifecycle and recursion GF variants | Validated | Green lifecycle/algorithm unit coverage plus functional real-space workflows |
| Reciprocal SCF, bands, and DOS | Validated | `kspace` functional group covers Si and magnetic bcc Fe, including tetrahedron, CCOR, and HOH cases |
| Reciprocal Green function / Lehmann route | Experimental | Unit-level chain coverage exists, but the Phase-I gate has no broad material-level validation for all downstream consumers |
| Exchange and conductivity | Validated | `quick` exchange coverage and the `conductivity` functional group, including Pt and PAOFLOW routes |
| Interfaces and vacuum leads | Experimental | Focused unit/wiring and Cu fixtures exist; unresolved interface/vacuum limitations remain below |
| Spin dynamics | Experimental | The bulk production-step smoke remains non-physical; VAL-11 separately validates the deterministic abspinlib Depondt solver only for its scoped prescribed-field problems |
| TDDFT / response | Development | Extensive conventions and unit tests exist, but active response work and limited end-to-end validation remain |
| GPU reciprocal acceleration | Development | Not part of this CPU Debug gate; no production-maturity evidence is recorded here |

## Known issues retained after STAB

Resolved entries already recorded in `tests/KNOWN_ISSUES.md` remain historical
evidence; they are not reopened here.  The unresolved entries are retained,
including the multi-sublattice k-space frozen-magnon Goldstone check, the
Cu(111) interface residual, the interface row-count alignment behaviour,
multi-layer `A | vacuum` triage, and the no-op magnetic constraining field.
Temporary guards and wiring fixes are not treated as general resolutions.

## Phase-I gate record

Configuration: `build-rf-debug`, CMake Debug, GNU Fortran 13.3.0
(`/bin/gfortran`); `BUILD_TESTING=ON`, `ENABLE_MPI=OFF`, `ENABLE_OPENMP=OFF`,
`ENABLE_CUDA_PLUGIN=OFF`, `ENABLE_LIBXC=OFF`, `ENABLE_SPGLIB=ON`,
`ENABLE_FUSED_RECIPROCAL=ON`, and `ENABLE_MARCH_NATIVE=ON`.

The CMake build directory was reconfigured and rebuilt immediately before the
recorded gate.  This matters: an earlier run used an out-of-date unit-test
executable and falsely reported `UnitMadl2dExhGuard` as failing.  The rebuilt
test passes and confirms the moved 2-D Madelung implementation correctly.

| Gate | Result | Record |
|---|---|---|
| `ctest -L unit --output-on-failure` | Pass | 44/44 passed |
| `ctest -L quick --output-on-failure` | Pass | 14/14 passed |
| `ctest -L kspace --output-on-failure` | Pass | 10/10 passed |
| `ctest -L conductivity --output-on-failure` | Pass | 4/4 passed |
| Focused exchange coverage | Pass | included in `quick`: 1/1 `exchange` fixture passed |

The lean gate intentionally excludes the wider validation campaign and the
performance/profile cases; MPI, CUDA, and MKL-specific variants are also
outside this serial CPU Debug configuration.  Those exclusions are scope
limits, not passes.

### Closure checklist

- [x] Architecture documentation reflects the current STAB organization.
- [x] Module/submodule map is documented.
- [x] Protected atomic LMTO core is documented.
- [x] Test taxonomy is documented using existing labels only.
- [x] Fixture vocabulary is documented.
- [x] `strux_lib` test policy is documented.
- [x] Feature maturity is stated conservatively.
- [x] Resolved known-issue history is retained and unresolved issues are not deleted.
- [x] Unit gate passes.
- [x] Quick gate passes.
- [x] Focused reciprocal, exchange, and conductivity functional groups pass.
- [x] Phase-I status is recorded.

**Phase I is complete at this head.**  The rebuilt lean and focused functional
gates pass with no observed STAB-introduced regression.
