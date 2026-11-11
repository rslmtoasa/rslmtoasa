# LDA+U tests

| Path | Content |
|---|---|
| `test_guards.py` | `LdaUGuards`: every LDA+U guard scenario in one test (nsp=1 host at nsp=2, U-less host, tilted axis with block recursion, mixed axes, host-axis mismatch, k-space local at nsp=4); prints PASS/FAIL per scenario |
| `test_gbt_hoh_u.py` | `LdaUGbtHoh`: real-space `gbt_single_q` at q = 0 against `periodic_nc` with U, `hoh`, `ccor_2c` (ferromagnet) |
| `test_kspace_gbt_u.py` | `LdaUKspaceGbt`: k-space `gbt_single_q` at q = 0 against `periodic_nc` with U at first and second order, ferromagnet and antiferromagnet |
| `test_kspace_local_nsp3.py` | `LdaUKspaceLocalNsp3`: `hubbard_radial_moments='local'` in k-space at `nsp=3` against `nsp=1`, antiferromagnet and ferromagnet, and `local` against `spectral` at `nsp=3` |
| `test_impurity_hoh_u.py` | `LdaUImpurityHoh`: `hubbard_radial_moments='local'` with `hoh` on the B2FeCo impurity, `nsp=2` (SOC held at zero) against `nsp=1`, and `local` against `spectral` |
| `test_impurity_nsp3_u.py` | `LdaUImpurityNsp3`: `hubbard_radial_moments='local'` at `nsp=3` on the B2FeCo impurity (first order and hoh) in an `nsp=3` k-space host against `nsp=1`, and `local` against `spectral` |
| `test_local_nsp4_workflows_u.py` | `LdaULocalNsp4Workflows`: `hubbard_radial_moments='local'` at `nsp=4` in k-space (AF and FM, first order and hoh) and on the B2FeCo impurity: SOC off against `nsp=3`, SOC on against `nsp=2` |
| `test_local_noncollinear_u.py` | `LdaULocalNoncollinear`: `hubbard_radial_moments='local'` at `nsp=3` with noncollinear site axes (real space): rotated, exchanged and generally rotated moment pairs against the same pair, parallel against `nsp=1`, `local` against `spectral` |
| `test_surface_hoh_u.py` | `LdaUSurfaceHoh`: `hubbard_radial_moments='local'` with `hoh` on the reduced Ni(001) surface (calctype S), `nsp=2` (SOC held at zero) against `nsp=1`, and `local` against `spectral` |
| `test_surface_nsp3_u.py` | `LdaUSurfaceNsp3`: `hubbard_radial_moments='local'` at `nsp=3` on the reduced Ni(001) surface (first order and hoh, `density_policy='constrained_spiral'`) against `nsp=1`, and `local` against `spectral` |
| `test_surface_nsp4_u.py` | `LdaUSurfaceNsp4`: `hubbard_radial_moments='local'` at `nsp=4` on the reduced Ni(001) surface with SOC held at zero (first order and hoh, `density_policy='constrained_spiral'`) against `nsp=3`, and `local` against `spectral` |
| `test_local_nsp4_u.py` | `LdaULocalNsp4`: `hubbard_radial_moments='local'` at `nsp=4` (real space): SOC held at zero against `nsp=3`, a ferromagnet with SOC against `nsp=2` and z against x, and `local` against `spectral` |
| `test_local_energy_u.py` | `LdaULocalEnergy`: the Hubbard energy of `hubbard_radial_moments='local'` at `nsp=3` and `4` (real space): AF `nsp=3` against `nsp=1`, SOC-off `nsp=4` against `nsp=3`, FM with SOC `nsp=4` against `nsp=2` |
| `test_missing_dm.py` | `LdaUMissingDm`: a fixed-U restart whose potential names a missing density-matrix file warns and continues (the density matrix then starts from the first Green function) |
| `test_impurity_equals_bulk.py` | `LdaUImpurityEqualsBulk`: Fe impurity in an Fe host (hoh, local, U = 2) reproduces the bulk spin moment and potential parameters when the impurity deck carries the bulk Fermi energy; the same run with a stale Fermi energy must miss it |
| `test_nonfinite_potential_guard.py` | `LdaUNonfinitePotentialGuard`: a k-space run with NaN potential parameters stops with the B1 guard message naming atom and field (before the eigensolver) |
| `reference_rows.py` | stored-reference rows: `reference_rows.py --list`, `reference_rows.py ROW --binary rslmto.x` |
| `references/` | JSON scalars of every row; five rows are registered with ctest, the others are archived data |
| `data/` | fixtures read by the Fortran oracles (`UnitHubbardSpinor`, `UnitLdaULocalMoments`) |
| `oracles/` | the independent numpy/SymPy generators of those fixtures (not run by ctest) |
| `decks/ni001_soc_u3/` | reduced Ni(001) surface deck and frozen hosts used by the surface row |
| `rotate_spin_frame.py` | rotates a host and seeds to a tilted axis (see `docs/validation/fixed_ldau/TESTING_POLICY.md`) |

References are rewritten only with `reference_rows.py ROW --binary BIN --update "reason"`; see `docs/validation/fixed_ldau/TESTING_POLICY.md`
and `tests/README.md` (reference policy). Rows with a `<row>.<platform>.json` use that file on that platform.
