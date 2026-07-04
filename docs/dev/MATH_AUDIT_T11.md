# T11 Math Module Audit

Scope searched: `source/`, `tests/`, and `example/` for `*.f90`, `*.F90`,
`*.cu`, and `*.cpp`. Counts below are conservative text-level live-use counts:
definition/end lines and comment-only lines are excluded, but same-name
procedures in other modules can still collide because Fortran is
case-insensitive. No routines were deleted in T11.

## P9 update (Phase 2): deletions executed

All 26 routines in the "Uncalled Candidates" table below were removed from
`source/math.f90` (719 lines). Before deleting, re-verified every entry
against the full source tree (not just the T11 scope) and specifically
checked for generic-interface/`module procedure` references, since a name
absent from direct call sites can still be live through an interface block —
`math.f90` has no `interface`/`generic` blocks at all, so none of these were
reachable that way. The two entries the audit itself flagged as needing a
closer look before trusting the raw grep count were re-checked by hand:

- `sph2cart`: all 8 text hits are the string literal `'sph2cart'` passed to
  `hcpx(..., 'sph2cart')` in `hamiltonian_paoflow_io.f90` — a mode argument,
  not a call to the `sph2cart` function. Confirmed dead, deleted.
- `trace`: all 3 text hits are inside comments/adjacent words ("orbital
  trace", "trace cap") in `bands.f90`/`exchange.f90`, not calls. Confirmed
  dead, deleted.

The three "Additional likely-unused routines with ambiguous text hits"
(`inverse`, `cartesian_to_direct`, `direct_to_cartesian`) were **not**
deleted — the audit's own hedged language for those ("apparent hits are
comments... only the function result declaration was found") reflects lower
confidence than the 0-count "Uncalled Candidates" table, and P9 only
approved the latter. They remain candidates for a future, separately
verified pass.

**Transitive orphans found and also removed:** the "Duplicate/Overlap Notes"
section below already noted that `erf_qe`/`erfc_qe` were only called by the
now-deleted QE smearing family (`theta_function`, `delta_function`,
`integrated_delta_function`). Re-grepped after deleting that family:
zero remaining external references to either. Removed both (99 lines) as a
direct, mechanical consequence of the approved deletions — not a new
independent audit finding, just following through on what this document
already predicted would happen.

Verified after all deletions: clean rebuild, full regression + example
matrix green (bit-identical results expected and observed, since these were
confirmed-dead code paths, not behavior changes).

## Uncalled Candidates (all removed, see above)

These had no live call/use sites in the search scope:

| Routine | Line | Count | Notes |
| --- | ---: | ---: | --- |
| `zeros1` | 372 | 0 | Array constructor helper; replaceable by allocation/assignment at call sites if ever used. |
| `zeros2` | 380 | 0 | Same as above. |
| `ones1` | 388 | 0 | Same as above. |
| `ones2` | 396 | 0 | Same as above. |
| `diag1` | 404 | 0 | Small diagonal-matrix helper; no live calls. |
| `diag2` | 418 | 0 | Small diagonal extraction helper; no live calls. |
| `kron` | 537 | 0 | Kronecker delta helper; no live calls. |
| `anticommutator` | 759 | 0 | Real matrix helper; no live calls. |
| `canticommutator` | 773 | 0 | Complex matrix helper; no live calls. |
| `ceye` | 785 | 0 | Identity-matrix helper; no live calls. |
| `eye` | 799 | 0 | Identity-matrix helper; no live calls. |
| `is_hermitian_r4` | 812 | 0 | Diagnostic helper; no live calls. |
| `unique_int` | 998 | 0 | Duplicates the generic `array_mod::unique`/`funique` machinery. |
| `sph2cart` | 1074 | 0 | No direct calls; `hcpx` has a string mode named `'sph2cart'`, which is separate. |
| `nm2rho` | 1083 | 0 | Density-matrix conversion helper; no live calls. |
| `fermi_function` | 1116 | 0 | Duplicates current `fermifun` family conceptually; no live calls. |
| `fermi_function_derivative` | 1153 | 0 | Duplicates `dfermifun` conceptually; no live calls. |
| `indexx` | 1171 | 0 | Local sort-index helper; no live calls. |
| `theta_function` | 1314 | 0 | QE smearing helper; no live calls. |
| `delta_function` | 1393 | 0 | QE smearing helper; no live calls. |
| `integrated_delta_function` | 1464 | 0 | QE smearing helper; no live calls. |
| `trace` | 2013 | 0 | Complex trace helper; `imtrace`, `rtrace`, and `trace9` are used instead. |
| `rotmag` | 2120 | 0 | `rotmag_loc` is the live rotation helper. |
| `sph2car` | 2374 | 0 | Legacy coordinate transform; no live calls. |
| `setup_rotation_matrix` | 2485 | 0 | No live calls. |
| `Coulomb_mat_direct_gaunt` | 2584 | 0 | Alternative Coulomb implementation; no live calls. |

Additional likely-unused routines with ambiguous text hits:

| Routine | Line | Apparent Count | Note |
| --- | ---: | ---: | --- |
| `inverse` | 926 | 1 | Apparent hits are comments/CUDA text, not live Fortran calls. |
| `cartesian_to_direct` | 2090 | 1 | Only live-looking hits are commented PAOFLOW lines or the function result declaration. |
| `direct_to_cartesian` | 2108 | 1 | Only the function result declaration was found. |

## Duplicate/Overlap Notes

- `unique_int` overlaps with the generic `array_mod::unique` and
  `array_mod::funique` helpers in `source/include_codes/array/`.
- `QSort` in `math.f90` textually collides with `array_mod`'s `qsort`
  procedures. The apparent count of 9 includes those same-name array-module
  calls and recursive calls inside `math.f90`; it should not be treated as a
  clean external use count.
- `zeros*`, `ones*`, `eye`, `ceye`, and `diag*` duplicate simple allocation,
  assignment, and diagonal loops. If deleted later, prefer direct code at the
  small number of call sites rather than a new abstraction.
- `fermi_function`/`fermi_function_derivative` overlap conceptually with the
  live `fermifun`/`dfermifun` pair, but their signatures differ.
- The QE smearing helpers (`theta_function`, `delta_function`,
  `integrated_delta_function`, `erf_qe`, `erfc_qe`) were internally related:
  the top-level smearing helpers were unused, and `erf_qe`/`erfc_qe` were
  only called by that unused family — so once P9 removed the top-level
  three, `erf_qe`/`erfc_qe` became dead too and were removed as well (see
  "P9 update" above).

## Live Low-Count Routines

Keep these for now despite low counts because they are part of live call
chains or physics-specific paths:

| Routine | Line | Count | Representative use |
| --- | ---: | ---: | --- |
| `factorial2` | 451 | 1 | `symbolic_atom.f90` crystal-field path. |
| `real_spharm` | 733 | 1 | `symbolic_atom.f90` crystal-field path. |
| `rho2nm` | 1094 | 1 | Internal warning/use in the conversion routine. |
| `dfermifun` | 1141 | 1 | Used by `simpson_f`. |
| `lorentz_kernel` | 1839 | 1 | Conductivity kernel damping. |
| `chebyshev_gauss_quadrature` | 1921 | 1 | Band moment integration. |
| `trace9` | 2077 | 1 | Exchange path. |
| `updatrotmom_single` | 2397 | 1 | Bands moment update. |

## Full Count Snapshot

`init_math_operators` 2; `fact` 13; `wigner3j` 2; `gaunt` 13;
`realgaunt_ord` 6; `realgaunt` 2; `genew` 1; `associated_lp` 1;
`complex_spharm` 3; `cross_product` 19; `normalize` 5; `erodrigues` 3;
`determinant` 2; `distance` 11; `angle` 9; `pos` 32;
`second_derivative` 4; `inverse_3x3` 3; `cart2sph` 1; `fermifun` 2;
`erf_qe` 2; `erfc_qe` 2; `hcpx` 40; `simpson_m` 15; `simpson_f` 40;
`jackson_kernel` 7; `t_polynomial` 6; `gauss_legendre` 5; `imtrace` 8;
`rtrace` 3; `scalar_multiply` 12; `imtrace9` 28; `rtrace9` 6;
`rotmag_loc` 11; `ROTMAT` 2; `DSs` 14; `binom` 1; `factint` 2;
`car2sph` 2; `tabulated_slater_integrals` 4; `a_k` 4; `Coulomb_mat` 9.
