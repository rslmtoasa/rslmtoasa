# Spinor d-shell occupation storage and extraction (S2 inventory)

Index map as in `SPINOR_FUNCTIONAL_CONTRACT.md`: compact spinor index `i=(sigma-1)*5+(m+3)`, which is
Hamiltonian row `(sigma-1)*spin_off + l*l + (m+l+1)`; for `l=2`, rows 5..9 (up) and `5+spin_off`..`9+spin_off`
(down). `P(a,b)=sum_n f_n psi_n(a) conj(psi_n(b))`.

## Inventory (read from source at `15a278b0`; no runs)

| item | finding |
|---|---|
| Basis | `basis_init` gives `norb=9`, `nb=2*norb=18`, `spin_off=norb=9` for `lmax=2` in every `nsp`; the log of an `nsp=1` run already reads `nb=18 spin_off=9`. `nsp=1` has two decoupled spin blocks; `nsp>=2` couples them (`build_lsham`, spin-off-diagonal blocks). |
| RS Green function | `green%g0(nb, nb, channels_ldos+10, atoms_per_process)` (`green_lifecycle.f90`). All four spin blocks live in the one 18x18 array: up-up rows/cols 1..9, down-down 10..18, spin-off-diagonal elements `g0(i, j+spin_off)` and `g0(i+spin_off, j)`. The Chebyshev driver fills every `(l,m)` of `mu_n(nb,nb,...)` and `g0(l,m,...)` (`green_chebyshev.f90`); `bands.f90` already reads `g0(j, j+spin_off)` for the transverse spin DOS. Block recursion filling of the off-diagonal blocks is not verified here. |
| RS occupations | `bands%calculate_moments` (`bands.f90`) builds `A_ab=i(G_ab-G_ba^*)/(2 pi)` (`hubbard_spectral_density_element`) only for `iorb+soff, jorb+soff` of the same spin and integrates with `simpson_m` (weight power 0 for P, 1 for M1) up to `en%fermi` on the mesh `en%ene`, `nv1`, `e1`. M1 goes to a local/optional `m1_target(nrec, l+1, spin, m, m')`. |
| KS eigenvectors | `reciprocal%eigenvectors(nb, nb, nk_local)` indexed (site offset + spin-major orbital, band), `eigenvalues(nb, nk_local)`; `calculate_ldm_from_eigenvectors` (`reciprocal_projection.f90`) already accumulates the full site matrices `DM` and `DM_m1` (all spin blocks) with Fermi-Dirac weights and k weights, then copies out only the `(ispin,ispin)` l-blocks into `ldm` and `m1_out`; MPI reduction covers real and imaginary parts of `m1_out`. |
| KS and `nsp=2` | the reciprocal path carries `n_spin_components` and a spinor-site-moment routine for the `nsp=2`/SOC path (`self_reciprocal.f90`), and the canonical fixture `Example_k_space_scf_bccFe_spd_magnetic` runs `nsp=2`. Not exercised by S2 yet. |
| ldm storage | `potential%ldm(l+1, spin, m, m')` complex, `(4,2,7,7)`, spin-diagonal only; `ldm_flatten`/`ldm_flatten_imag(l,spin,:)` for legacy inline restart. |
| dm sidecar | `<label>_dm.nml`, `dm_re/dm_im(l, spin_row, spin_col, m1, :)`; the namelist arrays are `(lmax+1,2,2,...)` but `potential.f90` writes and reads only `spin_row=spin_col` blocks. |
| M1 storage | not kept in `potential`; computed on demand (local mode, trace). |

## Implementation (S2)

- Key `hubbard_spinor_dm_diagnostic` (`&hamiltonian`, default `.false.`). With `nsp>=2` and the key on (or fixed U active, unreachable while the nsp>=2 guard stands) the 10x10 d-shell `P` and `M1` are stored in `potential%ldm_spinor` and `potential%m1_spinor`; both are allocated only then. The Hamiltonian does not read them. M1 is recomputed each iteration and not persisted.
- RS: `bands%calculate_spinor_dm`, called after `calculate_pl` in `calculate_moments`. It forms `A_ab=i(G_ab-G_ba^*)/(2 pi)` for the 10x10 block of `g0`, integrated with `simpson_m` on `en%ene` with `nv1`, `e1` and the Fermi endpoint of the collinear P (weight power 0 for P, 1 for M1); real and imaginary parts are summed over MPI ranks.
- KS: `reciprocal%calculate_spinor_dm`, called from the KS SCF after the density build. `P=sum f psi psi^dagger`, `M1=sum f eps psi psi^dagger` with the k weights and Fermi-Dirac occupations of the collinear producer; MPI reduction of real and imaginary parts.
- dm sidecar: with `ldm_spinor` allocated the file gains `spinor_d = 1` and the `dm_re/dm_im(3, s, t, :, :)` blocks for `s /= t`; the spin-diagonal blocks stay those of `ldm`. Files without the array are unchanged and carry no spinor content.
- Log line `HUBBARD_SPINOR_DM backend=rs|ks`: antiherm, trace, trace minus `ql_up+ql_dn`, trace of the up block minus `ql_up`, eigenvalue bounds, `max_up_dn`.

## Measured (S2)

Frozen-H probe `tests/validation/fixed_ldau/projector_spinor_frozen_h.f90` (18x18 H: crystal field with complex orbital mixing, exchange splitting, xi L.S with xi=0.35, up-down coupling; exact eigenpairs from `zheev`, E_F in the widest gap; two k points, 6000 channels, eta=0.004):

| check | measured | tolerance |
|---|---|---|
| KS P, M1 vs exact eigenpair sums | 1.4e-17, 2.3e-16 | 5e-12 |
| RS P vs finite-window Lorentzian | 4.2e-5 | 1e-4 |
| RS M1 vs finite-window Lorentzian [Ry] | 6.1e-7 | 3e-4 |
| RS vs independent Simpson at 6000 channels | 1.1e-16 | 5e-12 |
| mesh refinement 750..6000 decreases (P, M1) | 0 violations | 0 |
| Hermiticity (KS 5e-17; RS by construction) | 5.4e-17, 0 | 1e-14 |
| eigenvalue excursion outside [0,1] | 1.6e-15 | 1e-8 |
| Tr P_ss vs band-moment q0, both spins (KS) | 4.4e-16 | 5e-12 |
| Tr P_ss vs ql(1,d,s) (RS; same Green function and quadrature, by construction) | 3.1e-15 | 1e-10 |
| xi=0, no coupling: spin-off-diagonal blocks, KS and RS | 0, 0 | 1e-12 |
| xi=0: diagonal blocks vs the nsp=1 extraction, KS and RS | 1.4e-17, 0 | 1e-12 |

Controls that must exceed their tolerance: conjugated KS eigenvectors (P 0.31, M1 0.19 Ry), transposed P (0.31), swapped spin blocks (0.85 KS, 0.77 RS), RS Green function of H^T (0.31, 0.19 Ry), xi=0.35 without the coupling (0.44 in the off-diagonal block). The radial moments are resolved along the density-policy axis; the probe uses `constrained_spiral` (global z), as the production default does.

dm round trip (`UnitLdaUHamiltonian`): the four spin blocks of `ldm_spinor` survive write and read bitwise; a collinear file carries no `spinor_d` and no off-diagonal entries. A reader that zeroes the off-diagonal blocks makes the test fail.

Production smoke, U=0, flag on, 4 steps, single-threaded (RECORDS, not gates):

| deck | max antiherm | trace minus ql | eig_min, eig_max | max abs P_up,dn |
|---|---|---|---|---|
| RS `bulk/bccFe`, nsp=2, block recursion | 0 (by construction) | 7e-15 | 0.311, 0.948 | 1.245e-2 |
| KS `k_space_scf/bccFe`, nsp=2, lmax=2, 2x2x2 mesh (documented fixture namelists) | 1.1e-17 | 1.8e-15 | 0.217, 0.996 | 0.9e-2 to 1.8e-2 |

`control%has_soc()` is true for nsp=2 and nsp=4; the nonzero spin-off-diagonal block in a collinear z-magnetized Fe run is the effect of the l.s term. The RS block-recursion run populates the spin-off-diagonal blocks, so no Chebyshev run was needed.

## Inventory for inserting a spinor V (S3 T0; read from source at `9dcbecee`, no runs)

| item | finding |
|---|---|
| `l.s` operator | `hamiltonian%build_lsham` (`hamiltonian_build.f90`) fills `lsham(nb,nb,ntype)`: `prefac*Lz` (up-up), `prefac*(Lx-i Ly)` (up-down), `prefac*(Lx+i Ly)` (down-up), `-prefac*Lz` (down-down) with `prefac=0.5*xi` (p: `sqrt(xi_p(1) xi_p(2))`, d: `sqrt(xi_d(1) xi_d(2))`). Called at `nsp==2` or `4` from `calculation.f90:1849` and `calculation_reciprocal.f90`. |
| RS bulk | `locham = ee(:,:,m,type) + lsham(:,:,type)` (`hamiltonian_build.f90:1281`); the Chebyshev recursion adds `lsham` as a separate `zgemm` term (`recursion_chebyshev.f90`) and in the scaled cache (`chebyshev_fast.f90`). `build_bulkham` adds `hubbard_u_pot(:,:,type)` to `ee(:,:,1,type)` as a full `nb x nb` matrix (`hamiltonian_build.f90:1665`), as does the surface path (`:1904`). |
| RS cluster/impurity | `build_locham` adds `hubbard_u_pot(:,:,u_type)` to `hall(:,:,1,nlim)` as a full `nb x nb` matrix (`:2112`). `hall` is spin-resolved and already carries spin-off-diagonal blocks (`hmag` x,y). |
| KS | the k-space Hamiltonian is assembled from the same `ee` and `lsham` (`hamiltonian_build.f90:2556, 2598`), so a `hubbard_u_pot` entry in `ee(:,:,1,type)` reaches KS. |
| Consequence | because insertion is of the whole `hubbard_u_pot(nb,nb,ntype)` matrix, a spinor V written into its d rows and columns `(sigma-1)*spin_off + 4 + i` (including the up-down blocks) is added by the existing code, once per build. Only `calculate_hubbard_u_potential_general` (collinear, same-spin blocks) has to change. |
| Current guard | `hamiltonian_build.f90:756`: fixed U with SOC or noncollinear spinors is fatal. |
| SOC scaling | `soc_scale` in `&self` (default 1.0) multiplies `xi_p` and `xi_d` of every symbolic atom (`self.f90:1787`); there is no knob in `&hamiltonian` or `&control`. |
| Density policies | `control%density_policy` is `constrained_spiral` (default) or `relaxed_reference`; any other value is fatal (`control.f90:793`). Both resolve a per-site axis from `constraint_target` if nonzero, else `potential%mom` (z if that is zero) via `resolve_site_axis`/`set_axis`, and every later radial spin projection uses that axis. Under `constrained_spiral` the reference stays fixed and the unabsorbed transverse density is reported as residual and torque; under `relaxed_reference` the reference `potential%mom` follows the density. The axis is therefore z only when `mom` is along z and no constraint target is set. |

## Spinor V insertion at nsp=2 (S3; experimental, not merged)

Key `hubbard_soc_experimental` (default `.false.`). With fixed U, `nsp=2` and the key on: the nsp>=2 guard is lifted for `nsp=2` only; `calctype` other than `B`, a `density_policy` other than `constrained_spiral`, a `constraint_target`, or a site reference axis with |x| or |y| above 1e-8 is fatal; `hubbard_radial_moments` `local` and `auto` become `spectral` with a warning. `calculate_hubbard_u_potential_general` fills the d block of `hubbard_u_pot` from `hubbard_spinor_functional` and the spinor P (`ldm_spinor`, or the collinear `ldm` blocks before the first extraction), rows and columns `(sigma-1)*spin_off+4+i`, including the up-down blocks; the existing build routines add it once. Energy routing is unchanged (off).

Operator test `UnitLdaUSocInsertion` (bcc Fe nsp=2, U=2 eV, J=0.5 eV, a frozen Hermitian P with |V_up,dn| = 1.9e-2): inserted block minus the no-U build against the S1 V, RS `ee` 4.3e-19 and KS H(k) 1.4e-17 and 2.8e-17 (limit 1e-14); controls that must exceed 1e-14: 2V (count) 1.9e-2, transposed V 3.4e-2, spin-diagonal blocks only 1.9e-2. Collinear limit (zero coherences) against the existing path 1.9e-16 (limit 1e-14); control with coherences 2.7e-2. Guard cases (all stop with the stated message): `Example_impurity_soc_u_guard`, `Example_bulk_bccFe_soc_u_policy_guard`, `Example_bulk_bccFe_soc_u_tilted_axis_guard`; `Example_bulk_bccFe_nsp2_hubbard` (flag off) still stops at the SOC guard.

SOC -> 0 in production (bcc Fe, RS block, nsp=2 with `soc_scale=0` against nsp=1, 4 steps, fresh deck copies, single-threaded; differences nsp=2 minus nsp=1):

| quantity | (c0) U=0 baseline | (c1) U=2/J=0 spectral, flag on | bound max(1e-10, 10 x c0) | (c1) soc_scale=1 control |
|---|---|---|---|---|
| d occupation up (`ql`) | 8.0e-5 | 3.8e-5 | 8.0e-4 | 3.6e-5 |
| d occupation down (`ql`) | 3.6e-5 | 2.0e-5 | 3.6e-4 | 1.6e-5 |
| moment from `ql` | 1.7e-5 | 5.8e-5 | 1.7e-4 | 5.2e-5 |
| E_F (6 printed decimals) | 7.2e-5 | 3.0e-4 | 7.2e-4 | 3.0e-4 |
| logged "Spin moment of atom 1" | 8.8e-5 | **2.83e-3** | 8.8e-4 | 2.89e-3 |

The logged spin moment exceeds the bound. Per-iteration values of the two (c1) runs are identical at iteration 1 (spin moment 2.101202, E_F -0.003690, Diff 2.30e-2 in both) and differ from iteration 2 (1.8e-3, 2.5e-3, 2.8e-3). The U=0 pair differs at iteration 2 by 2.0e-5. The +U runs are far from converged (Diff 2.6e-2 after 4 steps against 9e-6 for U=0), so the same path difference is amplified; this reading was not tested further. Only the iteration-1 spin moment separates soc_scale=1 from soc_scale=0 (8.1e-4); later iterations are dominated by the path difference, so the control does not discriminate at 4 steps (E_F and the `ql` quantities are below the bound for the control). Not run: the 100-step SCF smoke (T2(d)) and the merge.

### S3.1: which P builds V at iteration >= 2 (hypothesis: raw P_spinor against mixed `ldm`)

Code: the collinear `ldm` is saved as old and new around the DOS step and blended by `mix%mix_ldm_linear` with `ldm_beta`, whose default is 1.0 (no blending); the result is written back to `potential%ldm` before the next Hamiltonian build. `potential%ldm_spinor` is overwritten by the extraction every iteration. With the default `ldm_beta` both are the latest P_out. The trace machinery (`write_fixed_u_trace`) is limited to nsp=1, so a temporary log line (not committed) was added to the spinor branch: at every Hamiltonian build of the nsp=2, soc_scale=0, U=2 eV, J=0 run the diagonal blocks of `ldm_spinor` minus `ldm` were 0.0E+00 (first build: spinor array not yet allocated, collinear blocks used). V diagonal averages (up, down), nsp=2 spinor path against the nsp=1 `HUBBARD_VDIAG` line: build 1 0.073499, 0.073499 in both; build 2 -0.051767, 0.012573 in both; build 3 -0.058575, 0.013292 against -0.058611, 0.013310; build 4 -0.059511, 0.013827 against -0.059560, 0.013853. Both paths therefore build V from the same P; V is identical through iteration 2 and the difference first appears in the P_out of iteration 2. The hypothesis is not confirmed.

Spin-moment definitions: the logged `Spin moment of atom` is `potential%mtot = |m|` of the projected-DOS spin integrals (all orbitals, up to E_F, that iteration's unmixed output); the `ql`-derived moment is `sum_l [ql(1,l,up) - ql(1,l,dn)]` from the potential file, the mixed band-moment charges after the last iteration.

### S3.1: converged comparison (replaces T2(c))

bcc Fe, RS block recursion, `recur='block'`, beta 0.01 Broyden, fresh deck copies from the fixture seed, 200 steps with `conv_thr=1e-30`, single-threaded; nsp=2 runs use `soc_scale=0`, the +U runs U=2 eV, J=0, `spectral`, flag on for nsp=2. Acceptance on the last 20 steps: U=0 nsp=2 Diff 1.5e-8 (window maximum 4.0e-8), nsp=1 Diff 7.6e-9 (3.3e-8), spreads at the 1e-6 print resolution; nsp=2 +U Diff 1.8e-8 (2.2e-8). The nsp=1 +U run did not meet the smoke criteria after 200 steps (Diff 2.7e-6, 20-step spread of the logged moment 1.2e-5) nor after the ladder step (200 more steps, beta 0.005: Diff 1.9e-6, spread 1.2e-5); its moment moved by 1e-6 between the two, so it is PROVISIONAL and its state is not in doubt at the size of the differences below.

| quantity (nsp=2 minus nsp=1) | U=0 baseline | U=2, J=0 | ratio | bound max(1e-10, 10 x baseline) | |
|---|---|---|---|---|---|
| d occupation up (`ql`) | +2.7e-4 | -2.99e-3 | 11.2 | 2.7e-3 | FAIL |
| d occupation down (`ql`) | -1.8e-4 | +2.48e-3 | 13.9 | 1.8e-3 | FAIL |
| moment from `ql` | +2.34e-4 | -5.53e-3 | 23.6 | 2.3e-3 | FAIL |
| logged spin moment | +2.34e-4 | -5.52e-3 | 23.6 | 2.3e-3 | FAIL |
| E_F (6 printed decimals) | -1.1e-4 | -2.9e-5 | 0.3 | 1.1e-3 | PASS |

States: U=0 nsp=2 / nsp=1 moment 2.104064 / 2.103830; U=2 nsp=2 / nsp=1 2.380491 / 2.386010. The nsp=2 minus nsp=1 difference changes sign between U=0 and U=2, so it is not a simple amplification of the U=0 offset. The cause is not identified. Established: the V insertion matches the collinear V at the first two builds (S3.1 trace), the spinor block fill equals the collinear path for zero coherences (1.9e-16), and the nsp=2 and nsp=1 converged states already differ by 2.3e-4 at U=0.
