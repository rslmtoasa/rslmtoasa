# nio_af: antiferromagnetic NiO (AFM-II), inputs only

Copied unchanged from `example/bulk/fccNiO/AF` (inputs only), plus the B1 variant of `Nid` and the script that makes it.

| item | value |
| --- | --- |
| cell | 4 sites, rhombohedral AFM-II cell, `lattice.nml` (vectors in units of alat) |
| alat | 4.1684002458, `wav = 1.2929345694`, `ct = 3.3` for all four types |
| sites (units of alat, order of the `label` list in `input.nml`) | 1 Niu (0,0,0), 2 Nid (1,1,1), 3 O1 (1/2,1/2,1/2), 4 O2 (3/2,3/2,3/2). The comment block inside `lattice.nml` lists the sites in another order; it is unchanged and not used. |
| Hamiltonian | `hoh = .true.`, `ccor_2c = .true.`, `screening = 'fitted'`, `strux_want_sdot = .true.` |
| exchange-correlation | `txc = 106` (libXC); the tests set `txc = 1` (native Barth-Hedin) because the CI build has no libXC |
| Ni moments (muB per Ni, collinear SCF, `txc = 1`, hoh + CCOR) | 0.881 on the deck's 12x12x12 mesh (35 steps, SCF difference 2e-7); 0.795 on the 6x6x6 mesh of the I2 test (converged) |

## Two conventions for the down-spin Ni

- **Mirrored** (the deck as it stands: `Nid.nml`): `mom = (0, 0, +1)`; the spin-1 channel of every spin-indexed field holds the *minority* values of Nid, so `Nid` is the spin-swapped `Niu` and the AFM order is carried by the potential data.
- **B1** (`Nid_B1.nml`): `mom = (0, 0, -1)`; the spin-1 channel holds the *majority* values, as for any site with `mom(3) < 0` in the k-space GBT assembler (the frame carries the orientation). Made by `make_nid_b1.py`: every spin-indexed field of `&par` (`center_band`, `width_band`, `gravity_center`, `shifted_band`, `obar`, `c`, `enu`, `ppar`, `qpar`, `srdel`, `vl`, `pl`, `pnu`, `qi`, `dele`, `ql` (third index), `ldm_flatten` (second index), `xi_p`, `xi_d`, `rac`, `hyper_field`) is swapped, `vrmax(2)` changes sign, `mom = (0, 0, -1)`. The script stops on a key it cannot classify.

`AppGbtKspaceNio_r2` shows which fields the fixed-potential Hamiltonian reads: leaving `center_band`, `width_band`, `shifted_band`, `obar` or `mom` unswapped changes the eigenvalues by 6e-3 to 4e-2 Ry; the others (`c`, `enu`, `ppar`, `qpar`, `srdel`, `vl`, `pl`, `ql`, `gravity_center`, `xi_*`, `rac`) change nothing at fixed potential and are swapped by construction only. `vrmax(2)` enters only with `ccor_vmt_mode = 'pair_surface'`, where the ordinary route is not invariant under the two encodings (`tests/KNOWN_ISSUES.md`) and the k-space GBT glue stops; the tests use `surface_scalar`.
