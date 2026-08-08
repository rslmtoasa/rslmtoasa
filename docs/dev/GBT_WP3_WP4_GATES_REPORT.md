# GBT WP3/WP4 operator-gate report

Date: 2026-08-03.

## Scope and architecture

This change completes only the documented WP3 representation split and WP4
first-order primitive-`S` link. WP5 deletion and higher-order feature
enablement were not started.

The magnetic representation is now explicit and independent of solver choice:

- `periodic_nc` retains `hmfind -> chbar_nc -> ham0m_nc` and its ordinary
  per-type potential moments;
- `explicit_texture` uses a site-indexed moment provider on the unchanged NC
  pair pipeline and rejects per-type bulk reuse if equal type labels carry
  unequal moments (true magnetic-supercell reuse remains valid);
- `gbt_single_q` bypasses `ham0m_nc`, requires `strux_backend='strux_lib'`,
  and constructs `ee` from raw primitive orbital `S` through a temporary
  spinor block.

`gbt_kspace` is accepted only as a deprecated input alias and is translated
once to `gbt_single_q`; it is then cleared. The retained reciprocal GBT
reconstruction is therefore inactive on the new production path, but remains
in the source for the WP5 deletion gate.

The GBT path uses one helper contract:

\[
\alpha=2\pi\,\mathbf q\cdot\mathbf d/a_{\rm lat},\qquad
G_{ab}=(U_a^0)^\dagger R_z(\alpha)U_b^0.
\]

The complete directed Cartesian neighbor vector supplied by the strux-lib
mapping is used. Raw persistent `lattice%sbar`/`sdot` remain orbital-sized and
are not mutated. `S` is lifted to a temporary `nb x nb` spinor block and
contracted with unchanged diagonal collinear `wx0/wx1` factors. Onsite
`cx0/cx1` terms remain collinear and unphased. Each production build writes
`gbt_bonds.out`, containing the directed geometry, dense link, raw `S`, and
gauged spin blocks.

No Sdot consumer was enabled: HOH/`eeo`, CCOR, overlap modes, Hubbard-V, SOC,
and unaudited GBT `local_axis` use remain fatal. `strux_want_sdot` is therefore
not required by the first-order path.

## Changed files

- `source/magnetic_representation.f90`: representation constants,
  normalization, site endpoint provider, and per-type-reuse validation.
- `source/gbt_structure.f90`: pure phase, endpoint-frame/link, orbital lift,
  and collinear contraction helpers.
- `source/hamiltonian.f90`, `source/hamiltonian_build.f90`, and the Hamiltonian
  namelist include: representation state/routing, explicit-texture provider,
  GBT builder, guards, bounds checks, and bond dumps.
- `source/calculation.f90`: frozen-magnon drivers select `gbt_single_q` for
  both reciprocal and recursion consumers.
- `source/reciprocal_lifecycle.f90`: representation-based full-BZ and feature
  validation.
- `source/reciprocal_bands.f90`: reported pre-eigensolver GBT Hermiticity
  maximum.
- CMake and unit-test files: WP3/WP4 tests plus production finite-q fixtures.

The untracked user MX/MY/MZ and LocalAxis fixture trees were read/copied only;
they were not edited, deleted, or regenerated.

## Validation evidence

### Algebra and representation

Commands:

```bash
cmake --build build_13 -j2
ctest --test-dir build_13 --output-on-failure -L unit
build_13/bin/UnitGbtStructure
```

All 16 unit tests pass. Relevant maxima/budgets:

| Check | Observed | Budget |
| --- | ---: | ---: |
| Existing G1 non-cubic/dense/reverse/Stoner oracle | `1.780e-15` | `<1e-12` |
| Production-helper dense link/lift/contraction suite | `1.2413e-16` | `<1e-12` |
| q=0 common-frame `max|G-I|` | `1.1102e-16` | `<1e-12` |
| q=0 pre-solver `max|H-H^H|` | `2.2206e-16` | `<1e-10*scale` |
| finite-q pre-solver `max|H-H^H|` | `1.6711e-16` | `<1e-10*scale` |

`UnitMagneticRepresentation` verifies two sites of the same chemical type with
different endpoint moments, a five-site skyrmion-like texture, rejection of
invalid per-type reuse, acceptance of true magnetic-supercell reuse, and that
`gbt_single_q` cannot enter the ordinary NC moment provider.

`UnitGbtStructure` uses an independent dense `2x2` construction and checks
phase units, unequal endpoint frames/species factors, reverse-link covariance,
raw-S preservation, q=0 common-frame identity, a q=0 relative-sublattice
non-identity, a finite-q non-no-op link, and the one-orbital Stoner
`k+/-q/2` identity.

### Ordinary NC and local-axis preservation

Six one-step production runs were made from the user's MX/MY/MZ fixtures, with
`local_axis` both false and true. Every run reports:

```text
Band energy of system: -1.9855801812 Ry
```

This reproduces the handoff baseline exactly at its reported precision.

### G2O production probes

Commands:

```bash
/Users/andersb/envs/p311/bin/python tests/run_test.py \
  --binary build_13/bin/rslmto.x --cases-json tests/gbt_wp2/cases.json \
  --case-name wp2_q0_fixed_potential_rotation \
  --scratch-root /private/tmp/rslmto_wp4_final --force-serial

/Users/andersb/envs/p311/bin/python tests/run_test.py \
  --binary build_13/bin/rslmto.x --cases-json tests/gbt_wp2/cases.json \
  --case-name wp4_finite_q_non_noop \
  --scratch-root /private/tmp/rslmto_wp4_final --force-serial

/Users/andersb/envs/p311/bin/python tests/run_test.py \
  --binary build_13/bin/rslmto.x --cases-json tests/gbt_wp2/cases.json \
  --case-name wp4_finite_q_realspace_consumer \
  --scratch-root /private/tmp/rslmto_wp4_final --force-serial
```

The fixed-potential q=0 common tilt reports the same canonical values for the
reference and rotated probes:

```text
EF     = -6.36108857e-2 Ry
N      = 8.00000000
EBAND  = -1.99686650 Ry
```

The resulting q=0 branch energy is `3.70650519e-15 Ry`, below the declared
`1e-10 Ry` acceptance budget. The largest electron-count residual is
`9.3863e-12`, below the WP2 `2e-10` budget.

The full-BZ finite-q reciprocal probe gives:

| q | canonical EBAND (Ry) | omega (Ry) |
| ---: | ---: | ---: |
| `0` | `-1.99686650` | `0` |
| `+0.05` | `-1.99684619` | `3.49543146e-4` |
| `-0.05` | `-1.99684619` | `3.49543146e-4` |

Thus the finite-q path differs from q=0 by `2.031e-5 Ry` and cannot pass as a
disabled probe. q and -q agree at the printed `1e-8 Ry` energy precision, and
their pre-eigensolver Hermiticity maxima are both `1.6711e-16`. The maximum
finite-q electron residual is `1.0481e-12`.

The recursion-consumer smoke test also completes on the same built `ee` path.
At the deliberately short recursion settings its q and -q EBAND values differ
by about `2e-8 Ry`; this is recorded as a convergence limitation, not treated
as the later G5 real-space/reciprocal agreement gate.

The narrow/wide DOS fixtures both retain canonical
`EBAND=-1.99686650 Ry`; the direct unit oracle reports zero canonical movement
under DOS window/grid changes.

The final two-rank MPI occupation test passes:

```bash
ctest --test-dir build_wp2_mpi --output-on-failure \
  -R UnitKspaceOccupations_mpi
```

The legacy-backend GBT launch stops with the required fatal diagnostic:

```text
gbt_single_q requires strux_backend='strux_lib'; the legacy backend is unsupported.
```

The ordinary legacy-backend MX/MY/MZ regressions above remain valid.

## Gates and remaining risks

- **G3: PASS.** The three modes are explicit, ordinary NC/local-axis behavior
  is preserved, site identity is tested, and unsafe per-type texture reuse is
  rejected.
- **G4: PASS for first-order S.** Dense/reverse/Stoner/q-symmetry/q=0/non-no-op
  and production Hermiticity evidence meet their budgets. Persistent S/Sdot
  storage is unchanged; Sdot-dependent terms remain guarded.
- **G2O: PASS.** The q=0 rotating-frame reduction and finite-q non-no-op
  production checks both pass.
- **Final G2: PASS.** G2E was already passing; G2O now closes the operator slice.

Unresolved work is intentionally deferred: WP5 must remove the now-inactive
reciprocal GBT reconstruction and perform converged RS/k-space bond/operator
comparisons. WP6 must audit Sdot, HOH/overlap, CCOR, Hubbard-V, SOC policy, and
other derived operators before any guard is relaxed. No physical bcc-Fe
dispersion interpretation or golden-data update was made.
