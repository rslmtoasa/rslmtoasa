# Phase 0 — environment and unmodified baseline

## Provenance and workspace

- Repository: `rslmtoasa/rslmtoasa`; branch: `fable_v4b_u` (user-selected).
- Audited commit: `7be98f7a8b5a0e7d1d3227821891403152b95608` is an ancestor of HEAD `56429aa18c1e9190d7612e96b759ac7b07fe2e7f` (2026-10-13).
- A direct diff of the audited LDA+U source and test files between those commits is empty. The current HEAD therefore contains the audited LDA+U implementation unchanged.
- Initial `git status --short`: only `?? cleanup/`, containing the user-supplied audit/spec files. Those files were preserved. No physics source or test files were modified in this phase.
- The audit asks for a new branch/worktree. The user's explicit instruction selected the existing `fable_v4b_u` branch, so evidence was produced in-place without switching branches or overwriting user changes.

## Toolchain and build

Host tools: macOS, CMake 4.4.4, GNU Fortran 16.2.0, Python 3.15.0 by default, Accelerate BLAS/LAPACK. The selected serial configuration used Debug, MPI/OpenMP/CUDA/LibXC/spglib/fused reciprocal disabled, and unit/example/regression targets enabled. CMake warned that GNU Fortran 16 on macOS skips `-fcheck=all` because of allocatable-descriptor runtime corruption. It also emitted existing conversion/pointer-lifetime warnings; there were no build errors.

Commands and exits:

| Command | Exit/result | Evidence |
|---|---:|---|
| `cmake -S . -B /tmp/rslmto-ldau-phase0-build -DCMAKE_BUILD_TYPE=Debug -DRUN_UNIT_TESTS=ON -DRUN_EXAMPLE_TESTS=ON -DRUN_REG_TESTS=ON -DENABLE_MPI=OFF -DENABLE_OPENMP=OFF -DENABLE_CUDA_PLUGIN=OFF -DENABLE_CUDA_RECIPROCAL=OFF -DENABLE_SPGLIB=OFF -DENABLE_FUSED_RECIPROCAL=OFF` | 0 | `/tmp/rslmto-ldau-phase0-cmake.log` |
| `cmake --build /tmp/rslmto-ldau-phase0-build -j 4` | 0, all targets built | `/tmp/rslmto-ldau-phase0-build.log` |
| `ctest --test-dir /tmp/rslmto-ldau-phase0-build -N` | 0, 200 registered tests discovered | `/tmp/rslmto-ldau-phase0-ctest-list.log` |

The initial example runner used Python 3.15, which lacks `f90nml`. CTest failed before starting either example with `ModuleNotFoundError`; no repository files were changed and no network install succeeded because package-index DNS is unavailable. Python 3.11.17 already has `f90nml` 1.4.4, so the build was reconfigured with `-DPython3_EXECUTABLE=/opt/homebrew/bin/python3.11 -DEXAMPLE_PYTHON_EXECUTABLE=/opt/homebrew/bin/python3.11`, then rebuilt successfully. All reported example/VAL-04 runs use this existing environment. The test wrapper dependency is therefore available; no package was installed.

## Fixed-U dependency map

The inspected executable path is:

1. `source/hamiltonian_build.f90::build_from_file` reads the `&hamiltonian` U/J inputs and converts eV to Ry (roughly lines 103–280). The general-U path skips an orbital channel if either its U or J entry remains at the sentinel; this confirms the audited U-only input issue.
2. `source/potential.f90::potential` stores `ldm` and `ldm_flatten` as real arrays (declaration around line 116; allocate/flatten/expand/namelist around lines 357–467, 782–831, 1012, 1543–1569). `source/mix.f90::save_ldm_to` mixes/stores flattened occupations (562–607). Restart and output use the same real flattened representation.
3. Reciprocal occupations enter through `source/reciprocal_projection.f90::calculate_ldm_from_eigenvectors` (838–972), then `sync_lattice_ldm` (237–274) handles distributed synchronization. The accumulation's complex matrix is reduced to real occupation storage. The build configured here had MPI off; MPI behavior was mapped statically but not exercised.
4. Real-space moments/occupations enter through `source/bands.f90::calculate_moments` (863–1127). For the Hubbard matrix it integrates `-Im G_ij/pi`, clips diagonal entries and applies a trace cap, diagnoses anti-Hermiticity, flattens the real array, and explicitly reduces it with MPI when enabled (around lines 975–1084).
5. `source/hamiltonian_hubbard.f90::calculate_hubbard_u_potential_general` (11–350) contracts the occupation with `Coulomb_mat` and constructs the real spin-resolved correction. `source/math.f90::Coulomb_mat`, `a_k`, and `realgaunt` (1780–1806; real Gaunt implementation at 484–551) construct tensor elements in real harmonics. Meanwhile `source/math.f90::hcpx` (998–1102), called during Hamiltonian assembly, transforms the spd Hamiltonian basis. The audited mismatch is present in the executable statements.
6. `source/hamiltonian_build.f90::build_bulkham` inserts the Hubbard matrix (around 1239–1276); `build_locham` inserts it for local/recursion Hamiltonians (1676–1709). The same potential object supplies the reciprocal onsite block.
7. `source/self.f90::is_converged` (2224–2238) checks the charge/energy delta and optional constraints; it does not require fixed-U occupation/potential convergence. `NEWRHO` energy bookkeeping sets `SUMEV`, `EKIN`, and `ETOT` (around 2645–2663) without an explicit scalar Hubbard interaction-minus-double-counting term. This is an accounting gap to derive in Phase 3, not evidence that a term should simply be appended.

This map is a navigation summary of current code and reported structural limitations; numerical claims below are from the unmodified executable. No Phase 1 tensor or orbital reference was attempted.

## Baseline tests and campaign

| Check | Measured | Expected/tolerance | Command and raw record |
|---|---|---|---|
| Standalone unit suite | 110/110 pass, total 80.51 s | all registered unit tests pass | `ctest --test-dir /tmp/rslmto-ldau-phase0-build -L '^unit$' --output-on-failure -j 4`; `/tmp/rslmto-ldau-phase0-unit.log` |
| No-U reciprocal bcc-Fe quick case | Runs deterministically but 9/18 compared values fail; Fermi level −0.0782523237 vs reference −0.0463779468 Ry; canonical band energy −1.92743828 vs −1.89229630 Ry; moment 2.000240 vs 2.002214 | reference abs/rel tolerance 1e−4/1e−5 | CTest and repeat direct invocation, both recorded below |
| Repeat of no-U reciprocal case | Same compared values and same 9/18 failures | same reference tolerance | `python3.11 tests/run_test.py --binary /tmp/rslmto-ldau-phase0-build/bin/rslmto.x --cases-json tests/scf/cases.json --case-name Example_k_space_scf_bccFe --scratch-root /tmp/rslmto-ldau-phase0-kspace-run2 --compare-ref tests/scf/references --abs-tol 1e-4 --rel-tol 1e-5 --force-serial`; `/tmp/rslmto-ldau-phase0-kspace-run2.log`, raw run under `/tmp/rslmto-ldau-phase0-kspace-run2/Example_k_space_scf_bccFe` |
| Reciprocal bcc-Fe magnetic quick case | Runs but 9/21 compared values fail; Fermi level −0.0785617377 vs −0.0461367419 Ry; moment 2.038397 vs 2.001427 | reference abs/rel tolerance 1e−4/1e−5 | `ctest --test-dir /tmp/rslmto-ldau-phase0-build -R '^Example_k_space_scf_bccFe_spd_magnetic$' --output-on-failure`; `/tmp/rslmto-ldau-phase0-fe-magnetic-smoke.log`; raw case under `/tmp/rslmto-ldau-phase0-build/Testing/scf/Example_k_space_scf_bccFe_spd_magnetic` |
| Real-space bcc-Fe Hubbard fixture smoke (U=2 eV, J=0, one iteration) | Pass | committed reference | `ctest --test-dir /tmp/rslmto-ldau-phase0-build -R '^(Example_k_space_scf_bccFe|Example_bulk_bccFe_nsp2_hubbard)$' --output-on-failure -j 2` (the Hubbard case passed); `/tmp/rslmto-ldau-phase0-fe-smoke-py311.log`; raw case under `/tmp/rslmto-ldau-phase0-build/Testing/scf/Example_bulk_bccFe_nsp2_hubbard` |
| VAL-04 U sweep | Script reports PASS. No-U and explicit U=0 `Fe_out.nml` are byte-identical. 2/4 eV conversion, d moment/charge and onsite splitting respond; all configured finite-step checks pass. | script's stated checks; U=0 vs no-U <1e−9; stored LDM Hermiticity <1e−7 and runtime anti-Hermitian residual <1e−6 | `python3.11 tests/validation/val04_lda_u.py --binary /tmp/rslmto-ldau-phase0-build/bin/rslmto.x --scratch-root /tmp/rslmto-ldau-phase0-val04`; `/tmp/rslmto-ldau-phase0-val04.log`; five preserved case directories under `/tmp/rslmto-ldau-phase0-val04/` |

VAL-04's outputs are finite-iteration records, not converged ground states. The four-step cases (baseline/U=0/U=2/U=4) and eight-step U=2 case all end without SCF convergence; the U=2 eight-step charge/energy residual sequence has first value `2.300011e-2` and minimum `4.693231e-3`. Across the U=2/U=4 cases, final logged anti-Hermiticity diagnostics are below `1e-6`; largest seen is `9.414e-8`. The script checks U=2/U=4 and the eight-step run, not baseline/U=0 occupation quality.

The repeated reciprocal no-U mismatches are categorized as an environment/reference baseline discrepancy: both failing examples run with no Hubbard input, the audited LDA+U files are byte-identical to the pinned source, and the repeated first case returns identical measured differences. The current macOS/GNU Fortran 16 debug build uses Accelerate and the generic committed references have no Mac-specific variant for these cases. This evidence does not identify which platform/compiler numerical pathway causes the offset; it is not attributed to LDA+U. The no-U RS Hubbard fixture, unit suite, and VAL-04 pass. No reference values or tolerances were edited.

## VAL-04 measured observables

All entries below are the finite-step outputs from the raw case directories. d traces are read from the stored `ldm_flatten`; U and J inputs are in eV and internal U is in Ry. U=0 has no retained Hubbard LDM on the inactive path, so its d trace is reported as unavailable by the campaign. `n_up/down` here denote the d-shell traces.

| Input U (eV), J=0 | Internal U (Ry) | d trace up/down | d moment | d charge | final avg potential split (Ry) | ETOT (Ry) | SUMEV (Ry) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| no Hubbard | 0 | unavailable | 2.199743 | 6.499257 | 0 | −2541.981942267 | −1.998357396 |
| explicit 0 | 0 | unavailable | 2.199743 | 6.499257 | 0 | −2541.981942267 | −1.998357396 |
| 2 | 0.1469971714 | 4.530831 / 2.031707 | 2.473729 | 6.555574 | −0.073223 | −2542.202638383 | −2.555956234 |
| 4 | 0.2939943429 | 4.711353 / 1.933212 | 2.635084 | 6.597957 | −0.163705 | −2542.379739744 | −3.030712458 |

The VAL-04 command, exact case inputs, full logs, and `Fe_out.nml` files remain in `/tmp/rslmto-ldau-phase0-val04`. `baseline/Fe_out.nml` and `u0/Fe_out.nml` have the same SHA-256 (`6ee007e1…bb3f06e`), as expected from the script's exact U=0 recovery check.

## Reproduction and artifact integrity

Build directory: `/tmp/rslmto-ldau-phase0-build`. Raw case/suite outputs are outside the repository. Selected SHA-256 values:

| Artifact | SHA-256 |
|---|---|
| `bin/rslmto.x` | `597480eaeecbae3292bf8e215673b0cd4189ccf4d07df055b2689fe910b56d2f` |
| build log | `75043e3847a974fe2ccc5c29ca7fbb16073251fed8bfb377956c0e2a784469b8` |
| unit log | `5a5056fa37364f0f51cc2af225364de033c338db0c37e1cb4623fefe4c710238` |
| VAL-04 console log | `e1dba9e864e846ca034baf4ca2e2c73b38b842e72c339427d470d10baaf8caee` |
| reciprocal no-U repeat comparison log | `9bd385ea3fa986f6de563eb1d9b9e59b2d4a9149d33dcb7d5dcca13518a85dd6` |
| reciprocal no-U first-run raw `testrun.log` | `1bbac07fff4e116d05e13bbd4fc7e56f07da03155192d86bed9ed87e81151d0f` |
| reciprocal no-U repeat raw `testrun.log` | `a33a609d7e3a87337acfe2a8750c43da2ae83d5e77256924fc3de4fb15bf3d1a` |

The VAL-04 raw output fingerprints are:

| Case | `Fe_out.nml` SHA-256 | `testrun.log` SHA-256 |
|---|---|---|
| no-U baseline | `6ee007e1d156ac890937096ec4d159343c8e239191bbb76e3f6f92827bb3f06e` | `b7746f610e767506b0f7eca86e3c1f8ceaa7b1225203e067aac0c1e61f1763c3` |
| explicit U=0 | `6ee007e1d156ac890937096ec4d159343c8e239191bbb76e3f6f92827bb3f06e` | `66442c68b8feec50153a6e0cad89dd52561d3601dd3124706ded2d2fbf59274d` |
| U=2 eV | `86bce7034152837806d09fb487759709834501f223c5675cbcf47c0dd1118d3c` | `965a4c7fdd700d4992c765b77fd262abe7f667cfa0ae4209c4cd45b84a425aa2` |
| U=4 eV | `ea5ed42494773f791a314281b44588c13ab5da17ef7a4bc5b18de8be3a5e4f59` | `f779e2c3f117178b3eeec4eba4fc4c3cda5ae8a4ff82c747112b495af2b0eaf8` |
| U=2 eV, eight steps | `bb29f165e31df306dcff2f3f427f12297f65956f01245a403d87e288ce5654ca` | `3ed19598960a8229444e1c6528a636c7e794234656aa3e9299b12bd9ef7685fa` |

Command logs `/tmp/rslmto-ldau-phase0-*.log` and raw case directories are intentionally not committed. Source tree remains unmodified; only this evidence directory is new.

## Gate 0 disposition

Build and focused existing tests ran. Unit tests, the no-U RS bcc-Fe Hubbard smoke, and the existing VAL-04 campaign pass; the required pre-fix no-U/U=0/U=2/U=4 records are preserved. The two reciprocal reference discrepancies were repeated and identified as no-U environment/reference baseline failures, without altering references or tolerances. Phase 0 evidence is ready for root review; the next agent should take Phase 1 only after acknowledging the documented no-U reference baseline.
