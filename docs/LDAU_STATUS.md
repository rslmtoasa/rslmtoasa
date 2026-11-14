# LDA+U status

Dated 2026-11-14. Measured values and open items only; status labels are assigned by the developer. Tree described: `main_rc_ldau_hoh` with the
branches `feat/ldau-surface-nsp4-test`, `feat/local-energy-clusters`, `feat/local-energy-hoh`, `fix/arbitrary-k-test-setup`,
`feat/local-radial-default` and `feat/ldau-status` merged in that order (tip of the last: the commit that adds this file).
Measurements below were made with gfortran 13, OpenBLAS, Release, no OpenMP, no spglib, 4 cores, one thread per run, unless stated otherwise.
Details and untriaged findings are in `tests/KNOWN_ISSUES.md`; the test index is `tests/ldau/README.md`; the theory and keyword text is
`docs/source/theory/lda_u_workflow.rst`.

## 1. What the code does

- **Hubbard form.** Fixed Liechtenstein U/J on the d shell with `lmax = 2`. Input through `hubbard_u_general`/`hubbard_j_general`, the legacy
  `hubbard_u`/`hubbard_j`, or `hubbard_u_impurity`/`hubbard_j_impurity` (impurity and surface clusters). Self-consistent U (SC-U) and intersite V
  are separate options that `local` does not implement.
- **Spin representations.** `nsp = 1` (scalar), `nsp = 2` (collinear spinor with SOC), `nsp = 3` (noncollinear, SOC off), `nsp = 4` (noncollinear
  with SOC). Antiferromagnet convention: at `nsp = 1` the second site's file has up and down swapped; at `nsp >= 3` both sites use the same file with
  `mom` +z and -z. At `nsp >= 3` set `density_policy = 'constrained_spiral'` for collinear-in-z work (section 5).
- **`hubbard_radial_moments`.** `local` takes the d-shell radial moments (c, q2) from the LDA part of the Hamiltonian with V_U's own contribution
  removed; `spectral` takes them from the Green function of the full H. `auto` (the default) is `local` when a nonzero U or J is present and `local`
  is implemented for the input, otherwise `spectral` (always `spectral` without U, where the mode makes no difference). The inputs `local` does not
  implement (`lmax` other than 2, a non-Liechtenstein form, SC-U, intersite V, SOC at `nsp = 1` or 3, hoh on a cluster outside the spinor routes)
  take `spectral` and the `HUBBARD_RUN` log line gives the reason; an explicit `local` stops with that message. An explicit setting is never changed.
- **Moving between modes.** A potential or frozen host written in the other mode has its density matrix reset, a warning is logged, and the run
  continues. A host with no recorded mode is accepted as is, with a warning. A frozen host reset this way is not self-consistent in the run mode.
- **Energy.** In `local` mode `physical_total_energy` includes `E_Hub` of the P_out of the bulk sites, or of the cluster sites for impurity and
  surface runs (host types excluded), at first order and with hoh. In `spectral` mode the Hubbard energy is not routed.
- **Routes with U.** Bulk real space (block and Chebyshev recursion), bulk k-space (first and second order), impurity clusters (calctype I) and
  surfaces (calctype S), with and without hoh and `ccor_2c` where the tests below say so.

## 2. Tests and what they agree to

Differences are `max |difference|` over the quantities the test lists (Diff, E_F, spin moments, d density-matrix entries, `E_Hub_P_OUT`,
`physical_total_energy`). The nsp = 3 against nsp = 1 and nsp = 4 against nsp = 3 rows compare two spin representations of the same collinear state
through the same code: they are consistency checks, not references. The `engaged` column is the measured size of the effect the test must see
(a skipped correction or routing would make the parities hold trivially).

| Test | Compares | Measured | Tolerance | `engaged` |
|---|---|---|---|---|
| `LdaUKspaceLocalNsp3` | k-space nsp 3 against 1 (AF, FM, AF with hoh), local against spectral | all differences 0.0 | 2e-6 | 1.3e-2 Ry (E_F) |
| `LdaULocalNsp4` | real-space nsp 4 SOC off against 3; FM with SOC against nsp 2; z against x | all differences 0.0 | 2e-6, 1e-5, 1e-5 | 6.3e-3 muB |
| `LdaULocalNsp4Workflows` | k-space and impurity nsp 4 (first order, hoh) | all differences 0.0 | 2e-6, 1e-7 | 1.3e-2 Ry (E_F) |
| `LdaUImpurityNsp3` | B2FeCo impurity nsp 3 against 1, first order and hoh | 1.5e-9 (hoh 0) | 1e-7 | 1.6 (hoh 2.5) muB |
| `LdaUSurfaceNsp3` | Ni(001) surface nsp 3 against 1, first order and hoh | 1.51e-8, 1.35e-8 | 1e-7 | 2.5e-2 muB |
| `LdaUSurfaceNsp4` | Ni(001) surface nsp 4 (SOC off) against nsp 3 | 0.0, 0.0 | 1e-7 | 2.5e-2 muB |
| `LdaULocalEnergy` | real-space energy: nsp 3 against 1, nsp 4 against 3, FM nsp 4 against 2 | `E_Hub` 8.6e-14, 0, 5.2e-13 Ry; total 3.1e-10, 0, 3.2e-10 Ry | 1e-8, 1e-6 Ry | 0.272, 0.246 Ry |
| `LdaUKspaceEnergy` | k-space energy, same pairs | `E_Hub` 5.9e-13, 0, 7.0e-12 Ry; total 4.4e-9, 0, 1.1e-9 Ry | 1e-8, 1e-6 Ry | 1.1e-1 Ry |
| `LdaUClusterEnergy` | impurity and surface energy, nsp 3 against 1 | impurity `E_Hub` 4.8e-10, total 1.8e-9 Ry; surface 5.4e-10, 1.5e-8 Ry | 1e-8, 1e-6 Ry | 1.0e-1 Ry |
| `LdaUHohEnergy` | energy with hoh: k-space AF, real-space AF, impurity, surface | `E_Hub` 9.2e-13, 1.7e-13, 7.4e-12, 3.5e-10 Ry; total 7.4e-10, 3.0e-10, 4.5e-10, 1.4e-8 Ry | 1e-8, 1e-6 Ry | 0.246, 0.253, 0.0975, 0.200 Ry |
| `LdaUDefaultRadial` | default against explicit local and spectral; no-U deck; restarts and a frozen host in the other mode | default equals explicit local bit for bit; against spectral 1.2e-3 muB | equality, > 1e-4 | n/a |
| `LdaURefSurfaceNsp2Soc0` | stored row: nsp 2 SOC off against nsp 1 on the surface | 1.51e-8 | 1e-7 | n/a |
| `LdaURefRsImpurityNsp2` (validation tier) | stored row, same for the impurity | 1.19e-9 | 1e-7 | n/a |

Other LDA+U tests (guards, GBT, hoh Fermi shift, noncollinear rotations, missing density-matrix file, NaN guard, impurity equals bulk, stored
reference rows) are indexed in `tests/ldau/README.md`. The independent oracles are the numpy/SymPy generators in `tests/ldau/oracles/` behind
`UnitHubbardSpinor` and `UnitLdaULocalMoments`.

## 3. Finite difference of the energy against U

J = 0, so `E_Hub` is linear in U and `dE_Hub/dU` at fixed P is `E_Hub/U`. Central differences of the converged `physical_total_energy`
(`conv_thr` 1e-9, 45 to 100 steps per run), reference `E_Hub_P_OUT/U` at U = 2. Criterion set before the bulk runs: 5e-3. The documented ASA base
residual (`docs/validation/fixed_ldau/ENERGY_ACCOUNTING.md`) is 1.4e-3 to 1.7e-3.

| System | FD (Ry/eV) | Reference (Ry/eV) | Relative difference |
|---|---|---|---|
| bcc Fe, nsp 1, RS block, first order (U = 1.8/2.2) | 0.05424312 | 0.05442266 | -3.3e-3 |
| same, hoh | 0.05342474 | 0.05355919 | -2.5e-3 |
| Fe impurity in an Fe host (host E_F, first order) | 0.05383704 | 0.05442260 | -1.08e-2 |
| B2FeCo impurity, host E_F, first order | 0.05151809 | 0.04780857 | +7.76e-2 |
| B2FeCo impurity, deck E_F 43 mRy above the host's | 0.05102763 | 0.04650623 | +9.7e-2 |

The bulk values are 2 to 3 times the documented residual. The cause of the B2FeCo difference is not known (`tests/KNOWN_ISSUES.md`, 2026-11-13).
The first-order finite difference of a cluster whose E_F differs from its host's is further off: nothing in the code checks the two E_F values.

## 4. Comparison with Quantum ESPRESSO

MnO AFM-II, PBE, U = 2 eV, J = 0, shift of the occupied d-up centroid against the occupied O p centroid (U = 2 minus U = 0), from
`tests/KNOWN_ISSUES.md` (2026-11-08): `local` -1.11 eV, `spectral` -1.72 to -2.00 eV, QE -1.24 / -1.29 eV. NiO and FeO comparisons were reported
by the developer and their data are not stored in this repository; the FeO comparison is not used as a target because of methodological differences
with QE. The derivation of `local` is in `docs/validation/fixed_ldau/LOCAL_MOMENTS.md`; no literature source for the ASA radial-moment step was found.

## 5. Open items

Not tested or not resolved, in the order I would take them:

1. **B2FeCo energy.** The finite difference differs from `E_Hub/U` by +7.76e-2 where Fe in Fe gives -1.08e-2. Cause not known. Untested guesses:
   an energy term between the two inequivalent cluster types that is not stationary; the Co site's response.
2. **SOC-on clusters and surfaces at nsp 3 and 4.** Not tested here (a separate effort). `LdaUSurfaceNsp4` holds SOC at zero and is bit-identical to
   nsp 3 because the SOC matrix is zero and the correction code is shared for `nsp >= 3`; it does not test the correction itself. nsp 4 cluster
   energies take the same route and are not tested.
3. **`relaxed_reference` axis update at `nsp >= 3`** (not specific to U). With the default `density_policy`, a near-zero moment can flip its axis:
   the Ni(001) surface at nsp 3 differs from nsp 1 by 2.6e-2, and with `magbeta = 0` (axes held) it agrees to 1.5e-8. Antiferromagnets in the
   `nsp = 1` file convention turn ferromagnet-like at `nsp >= 3` under the default policy. `constrained_spiral` avoids both.
4. **Stored references and old fixtures.** The stored impurity reference rows and `test_impurity_hoh_u.py` still use a deck E_F of -0.041385 Ry with a
   4-step host (about 43 mRy off the converged host). `reference_rows.py` rows are not re-baselined. `host_and_impurity` in
   `test_impurity_nsp3_u.py` converges the host and copies its E_F and window.
5. **PBE.** `atomsc` diverges to NaN for the spectral radial moments under PBE (`tests/KNOWN_ISSUES.md`, OPEN, partly diagnosed); the proposed fixes
   need edits in `self.f90` and are untested.
6. **Other entries** in `tests/KNOWN_ISSUES.md` that concern LDA+U and that I did not re-check in this round: RS rc=100 non-finite runs, Co bulk at
   U=0, the Fe U>0 embedding offset, SOC+U bulk untested corners, SOC+U impurity at nsp 2, the constraint controller in a KS nsp 2 smoke, and the
   nsp 3 block-recursion dependence on the quantization axes for the AFM bcc Fe cell.
7. **SC-U and ACBN0 wording.** The sections "Self-consistent U update (SC-U)" and "Optional ACBN0-like potential path" of `lda_u_workflow.rst`
   were not corrected. The developer's explanation is that SC-U computes a U that is inert in H here, and that the ACBN0-like routines are not ACBN0;
   no fatal guard was added for SC-U by decision.
8. **Not tested at all:** rotated rhombohedral FeO with J != 0, `local_axis` combined with U (no guard), nsp 3 finite difference, J != 0 finite
   difference, hoh with `ccor_2c` energy, Chebyshev hoh finite difference.

## 6. Platform notes

On the container used for these measurements (gfortran 13, OpenBLAS), without any LDA+U change: `UnitLdaULocalMoments` fails (a bitwise comparison
at 2.2e-16, platform-dependent, `tests/KNOWN_ISSUES.md` 2026-11-05); `LdaUGuards` fails the scenarios `surface_nsp1_host`, `surface_u_less_host`,
`tilted_axis_block` and `mixed_axes`; the stored rows `LdaURefKsBulkNsp2SocUTilted` (`ks_soc_u2_x`, max difference 2.21e-02) and
`LdaURefNoUKeyOutputs` (`noU_rs`, 1.90e-09 against 1e-10) fail, and so does `SpinResponseBaselineRepro`. `UnitArbitraryKEigenpairs` segfaulted until the test setup was fixed (the hand-built
lattices lacked the potential parameters that the NaN guard reads). The tier results of the final tree are in the next section.

## 7. Full standard tier of the final tree

Build with `-DRUN_UNIT_TESTS=ON` only (example, regression and validation tests not enabled), `ctest -LE 'validation|performance' -j3`:
83 tests, 78 passed, 5 failed, 2402 s on 4 cores. All five failures also occur on the unmodified `main_rc_ldau_hoh` tip (f016c6d) in the same
container:

| Failing test | Measured | Tolerance |
|---|---|---|
| `UnitLdaULocalMoments` | bitwise sub-case `spinor_frame_kernel` 2.22e-16 | 0 |
| `LdaUGuards` | 4 of 6 scenarios fail: `surface_nsp1_host`, `surface_u_less_host`, `tilted_axis_block`, `mixed_axes` | n/a |
| `LdaURefKsBulkNsp2SocUTilted` | `ks_soc_u2_x` max difference 2.21e-02 | 1e-10 |
| `LdaURefNoUKeyOutputs` | `noU_rs` max difference 1.90e-09 | 1e-10 |
| `SpinResponseBaselineRepro` (not LDA+U) | `spin_response_q001.dat` 1.318, `q002` 0.581, `q003` 0.309, `state` 9.7e-3 | 1e-8 |

The `LdaUGuards` scenarios and the two stored-row numbers were not examined further; `UnitArbitraryKEigenpairs` and the other quick-tier tests pass.
The developer's own full run (144 tests, with the example and regression tests) is not reproduced here.

## 8. How to run

`ctest -L '^quick$'` (quick tier), `ctest -L lda_u` (LDA+U tests), `ctest -LE 'validation|performance'` (full standard tier). Build with
`-DRUN_UNIT_TESTS=ON`; ctest does not rebuild. A single test: `python3 tests/ldau/<test>.py --binary <rslmto.x> --scratch-root <dir>`. Typical
runtimes on 4 cores with one thread per run: `LdaUKspaceEnergy` 24 s, `LdaUDefaultRadial` 206 s, `LdaUClusterEnergy` 394 s, `LdaUImpurityNsp3` 401 s,
`LdaUSurfaceNsp3` and `LdaUSurfaceNsp4` about 800 to 1300 s, `LdaUHohEnergy` 878 s under load.
