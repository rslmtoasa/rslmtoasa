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

### S3.2: the SOC-to-zero comparisons need `fix_soc`

Cause of the nsp=2 versus nsp=1 offsets recorded above (S3 c0/c1 and S3.1): `soc_scale` multiplies `xi_p` and `xi_d` once, when the restart potential is read (`self.f90:1787`); after every atomic SFC `xi_p`, `xi_d` and `rac` are recomputed from the radial potential unless `fix_soc = .true.` (`self.f90:5711`). With `soc_scale=0` alone SOC is off in iteration 1 only. That matches iteration-1 equality and divergence from iteration 2, the same offset in RS and KS, and `soc_scale=0` giving the same converged state as `soc_scale=1`. Those earlier comparisons therefore measured the real SOC effect, not a spinor-path difference. The valid SOC-to-zero test is `soc_scale=0` with `fix_soc=.true.`.

KS gate (bcc Fe, `k_space_scf/bccFe` namelists with `lmax=2`, 12x12x12 mesh, fresh deck copies, all runs converged to Diff < 1e-8; nsp=2 with `soc_scale=0`, `fix_soc` against nsp=1; U=2 eV, J=0, spectral, flag on for nsp=2):

| quantity (nsp=2 minus nsp=1) | A0, U=0 | A1, U=2 | bound max(1e-7, 10 x A0) | A1 control: `soc_scale=1` |
|---|---|---|---|---|
| d occupation up (`ql`) | 0 | 9.4e-10 | 1e-7 | 1.27e-3 (exceeds) |
| d occupation down | 0 | 5.2e-10 | 1e-7 | 1.00e-3 (exceeds) |
| moment from `ql` | 0 | 1.4e-9 | 1e-7 | 2.40e-3 (exceeds) |
| logged spin moment (6 decimals) | 0 | 0 | 1e-7 | 2.39e-3 (exceeds) |
| E_F (6 decimals) | 0 | 0 | 1e-7 | 2.7e-4 (exceeds) |
| P diagonal blocks (dm sidecar) | n/a | 6.7e-11 | 1e-7 | 2.06e-2 (exceeds) |

GATE passes. The same runs without `fix_soc` (soc_scale only) differ from nsp=1 by 2.4e-4 (U=0) and 1.2e-3 to 2.4e-3 (U=2) in the d occupations and moments, and equal the `soc_scale=1` runs to 1e-8.

RS (bcc Fe, RS block, 200 steps, `conv_thr=1e-30`, same decks as S3.1 plus `fix_soc`): nsp=2 with `soc_scale=0`, `fix_soc` reproduces nsp=1 at every iteration: the per-iteration logged spin moment, E_F and Diff differ by 0.0 over 200 iterations at U=0 and at U=2 (J=0, spectral, flag on), and the final `ql` occupations agree in every printed digit; the potential files differ only in `xi_p`, `xi_d`, `rac` and the sidecar's `spinor_d` blocks (zeros). The nsp=1 +U trajectory itself stalls at Diff 2.7e-6 (beta 0.01), so that pair is compared at equal step count, not at convergence. The RS recursion-terminator diagnosis (S3.2 B) was not needed: there is no RS offset to explain.

Real SOC effect (converged, nsp=2 full SOC minus nsp=1; moment from `ql`): RS U=0 +2.3e-4 (2.104064 against 2.103830), RS U=2 -5.5e-3 (2.380491 against 2.386010); KS U=0 -2.8e-4 (2.086904 against 2.087188), KS U=2 -2.35e-3 (2.353949 against 2.356302).

Smoke (T2(d)): bcc Fe, nsp=2, RS block, U=2 eV, J=0, flag on, spectral, `soc_scale=1` (default), fixture seed, 100 steps, `conv_thr=1e-30`, `hubbard_spinor_dm_diagnostic` on. GATE (final Diff <= 1e-7, 20-step spreads <= 1e-5, no guard or NaN): Diff 1.0e-8 (window maximum 4.8e-8), 20-step spread of the logged moment 1.0e-6 and of E_F 0, no fatal or NaN: passes. Records: spin moment 2.380490 muB, orbital moment 0.0979 muB along z (0.0417 muB at U=0 converged, same code and seed); max |P_up,dn| 1.26e-2 at the first step and 1.53e-2 at the last (1.245e-2 at U=0, S2); P eigenvalues in [0.270, 0.964]; Hermiticity 0 and |trace - ql| 1.6e-14 (both by construction in RS).

## Spinor local radial moments (S4; experimental, flag `hubbard_soc_experimental`)

Kernel `correct_d_shell_moments_spinor` (`math.f90`), per global-z spin channel s with Pi_s the spin-s projector on the 10-component d block and V the d block of `hubbard_u_pot` (spin-off-diagonal blocks included): `q0_s = Tr(Pi_s P)` (passed in), `t_s = Re Tr(Pi_s V P)`, `c'_s = c_s - t_s/q0_s`, `q2'_s = q2_s - 2 Re Tr(Pi_s V M1c_s) + Tr(V Pi_s V P) - t_s^2/q0_s`, `M1c_s = M1 - c_s P`, with `q2_s = sum_n f_n (eps_n - c_s)^2 a_n^dagger Pi_s a_n`. Derived independently from `b'_s = Pi_s[(eps - c'_s) a - V a]` by expanding `sum_n f_n |b'_s|^2` with `Pi_s^2 = Pi_s` and Hermitian V; it agrees with the formulas given for S4. Only V is subtracted.

Wiring: with the flag and nsp=2, `local` and `auto` stay local (provenance `local`, energy routing off); RS (`bands.f90`) and KS (`self.f90`) extract the spinor P and M1 first and apply the spinor correction to both spin channels. R05 for the flagged nsp=2 run uses the full 10x10 `P_in`, `P_out`, `V_in` (as inserted) and `V_out(P_out)` through `fixed_u_residual_metrics` (now accepting one block per site), with the usual backend targets (RS P 5e-7, V 1e-7 Ry). `RSLMTO_SPINOR_TRACE=<file>` dumps P, M1, V_in and the stored q0, c, q2 as they enter `calculate_pl`.

Unit tests (all with controls that must exceed their tolerance):

| check | measured | tolerance | control |
|---|---|---|---|
| (a) block-diagonal spinor kernel vs `correct_d_shell_moments` | 0.0 | 1e-14 | V spin blocks swapped 0.80; spin-off-diagonal V dropped 4.3e-2 |
| (b) finite H (xi L.S + spin-mixing coupling + spinor V_U, exact `zheev`): c' vs the direct b' definition | 1.1e-16 | 1e-12 | off-diagonal V dropped 0.17 |
| (b) q2' vs the direct b' definition | 1.4e-17 | 1e-12 | off-diagonal V dropped 2.5e-2; Pi_s V V in place of V Pi_s V 9.8e-4 |
| (c) q2' >= 0 on 40 random physical states | minimum 0 (none negative) | 1e-12 | Tr(V Pi_s V P) dropped: 1.21 negative |
| R05: only off-diagonal blocks change over 6 steps | never converged | | diagonal-block (collinear) layout of the same sequence converges |

Production checks (bcc Fe, U=2 eV, J=0, local mode, flag on for nsp=2; every `soc_scale=0` run also `fix_soc=.true.`):

- T3(b) cross-run identity (RS, one step each from the converged S3 smoke state, spectral and local; the local run logs `restart seed in spectral, run uses local; seed will be rewritten` at info level, no frozen-host guard in a bulk run): P, M1 and V_in are bit-identical between the runs. Stored local c' minus the formula on the spectral c, q2: 0.0 and 0.0 (limit 1e-12); q2': -1.4e-17 and 0.0 (limit 1e-10). Control without the spin-off-diagonal V blocks: c' misses by 1.7e-5 and 3.8e-5, q2' by 2.0e-5 and 2.9e-5, so the coherence contribution is small but resolved.
- T3(a) KS, converged (Diff 3.5e-9, 40 steps each), nsp=2 minus nsp=1 local: d occupation up 7.8e-10, down 7.3e-12, moment from `ql` 1.7e-10, logged spin moment 0, E_F 0 (6 decimals), P diagonal blocks 3.7e-11; bound 1e-7 (the U=0 difference is 0): passes. State: moment 2.51078 muB, E_F -0.08295 Ry.
- T3(a) RS, first 20 iterations, nsp=2 against nsp=1 (local, `fix_soc`): iterations 1-11 equal in every printed digit (Diff to 10 digits, spin moment, E_F, R05 residual); from iteration 12 they differ (Diff 6.5e-9 at 12, 6.6e-8 at 15, 6.2e-7 at 20; spin moment 1e-6 to 1.2e-5), so the gate "identical or within 1e-10" FAILS as written from iteration 12. Control: an nsp=1 run with its seed changed by a relative 1e-13 differs from the unperturbed nsp=1 run by 7.3e-9 in Diff already at iteration 3 and by 1.1e-6 at iteration 20 (spin moment up to 7e-6), as much as the spinor run. The U=2 local transient has large steps (E_F -0.103 to -0.071 at iteration 7), which amplify any tiny difference. The nsp=2 iteration-20 P diagonal blocks were not recorded because the runs are not identical.
- T3(c) smoke (RS block, nsp=2, SOC on, U=2/J=0, local, `hubbard_spinor_dm_diagnostic` on, fixture seed, 150 steps, `conv_thr=1e-30`): final Diff 9.2e-9 (window maximum 1.9e-8), 20-step spreads 0 (logged moment) and 0 (E_F), no fatal or NaN: passes. Records against the spectral smoke of S3: spin moment 2.551792 (2.380490), orbital moment 0.0688 muB (0.0979), E_F -0.083471 (-0.134657); max |P_up,dn| 1.262e-2 at step 1, 1.294e-2 at step 150 (1.26e-2, 1.53e-2); P eigenvalues in [0.259, 0.966]. R05 on the full spinor matrices: rp_abs 2.14, 0.119, 1.5e-3 (iteration 5), 1.1e-2 (10), 9.1e-4 (20), 4.4e-7 (50), 9.1e-8 (100); rv_abs 0.314, 1.7e-2, 2.2e-4, 1.6e-3, 1.3e-4, 6.4e-8 (50); the three-update streak is first reached at iteration 49.

### S4: the RS gate as redefined (replaces the iteration-by-iteration criterion)

(i) nsp=2 (`soc_scale=0`, `fix_soc`) and nsp=1 RS local trajectories (U=2 eV, J=0, fixture seed, beta 0.01) are identical through iteration 11 in every printed digit: PASS. (ii) The iteration-by-iteration criterion over 20 iterations is retired for local-mode +U transients: an nsp=1 run with its seed changed by a relative 1e-13 differs from the unperturbed nsp=1 run by 7.3e-9 in Diff at iteration 3 and 1.1e-6 at iteration 20 (spin moment up to 7e-6), so no equivalent implementation can meet it. Converged-state check:

RS block recursion, fixture seed `tests/scf/cases/bulk/bccFe/Fe.nml` (sha256 `6ded298a…`), U=2 eV, J=0, local, beta 0.005 (the RB2 continuation setting; the RB2 seed, sha `64f34f28…`, and its run directories are not available, so the recipe could not be matched and both partners were run), 150 steps, `conv_thr=1e-30`, per-iteration potential snapshots. Both reach the RS acceptance window: nsp=2 final Diff 1.56e-8 (window maximum 2.7e-8), spreads of d up, d down and moment over the last 20 iterations 1.1e-8, 1.6e-8, 1.0e-8; nsp=1 Diff 9.9e-9 (1.7e-8), spreads 3.8e-9, 1.2e-9, 9.3e-9.

| quantity | nsp=2 | nsp=1 | difference | bound |
|---|---|---|---|---|
| d occupation up | 4.615565472 | 4.615565620 | 1.48e-7 | 1e-6 |
| d occupation down | 1.895463860 | 1.895463740 | 1.21e-7 | 1e-6 |
| moment from `ql` (muB) | 2.558679698 | 2.558679966 | 2.68e-7 | 2e-6 |
| E_F (Ry, 6 decimals) | -0.083263 | -0.083263 | equal | equal |
| P diagonal blocks (dm sidecar) | | | 1.17e-7 | 1e-6 |

GATE passes. For reference the RB2 accepted nsp=1 state (other seed and recipe) has d occupations 4.615579, 1.895454 and moment 2.558681 (`C1_RB_RESULTS.md`). The iteration-20 P diagonal blocks of both runs (and the iteration-150 ones) are in `tests/validation/fixed_ldau/data/fe_bcc_u2_local_p_diagonal_blocks.txt`; at iteration 20 the two runs differ by 3.95e-6, at iteration 150 by 1.17e-7.

## nsp=2 energy candidates and Hellmann-Feynman checks (S5; experimental, flag `hubbard_soc_experimental`, energy not routed)

Code: with the flag, `HUBBARD_MOMENT_ENERGY_SITE/TOTAL` take `E_int`, `E_DC`, `E_Hub` from `hubbard_spinor_functional` on the full 10x10 `ldm_spinor` (P_out) (`spinor_site_hubbard_energy`, `self.f90`); `candidate_b = E_phys + E_Hub` (local), `candidate_c = E_phys + E_Hub` (spectral), `_F` versions subtract T*S (KS). `physical_total_energy` is unchanged (`E_Hub_P_OUT` stays 0 at nsp=2). The W1 test operator is accepted at nsp=2 with the flag (`hamiltonian_build.f90`). Driver and analysis: `tests/validation/fixed_ldau/soc_u_hf.py`, `soc_u_fixed_p_derivative.f90`.

T1 SOC to 0 (KS, bcc Fe, 12x12x12, U=2 eV, J=0, local, both converged; nsp=2 with `soc_scale=0`, `fix_soc`; the nsp=1 partner is the converged S4 KS run `k_nsp1`, same code path, run with the S4 binary): candidate b -2541.862184127762 against -2541.862184127180 Ry, |delta E| 5.8e-10 (gate 1e-7); E_Hub 0.110479378744488 against 0.1104793787442344, |delta E_Hub| 2.5e-13 (gate 1e-9); candidate b_F differs by 5.8e-10. nsp=2 run: 40 steps, Diff 1.36e-8. Control (SOC on): candidate b_F of the SOC-on KS runs at U=1.95 and 2.05 eV averaged (-2541.862686138 Ry, midpoint of two runs, not a run at 2.0) differs from the `soc_scale=0` value (-2541.862322177) by -3.6e-4 Ry, far above the 1e-7 gate.

T2(a) W1, U=J=0, SOC on, local, lambda = +-0.005 Ry, identity on the d block; F = candidate b (KS: F with -TS) + lambda Tr(W P), Tr(W P) from the d occupations `ql` (the logged `e_external` is 0 at U=J=0, see KNOWN_ISSUES). dF/dlambda against Tr(W P): KS 6.489010 against 6.498365 (-0.144%), RS 6.488623 against 6.499235 (-0.163%); gate 1%. Spin-selective runs (spin 1 only) and `R^T g = B`: g_s (Ry/electron) KS 0.00636, 0.00919; RS 0.00747, 0.00896.

T2(b) U derivative at U0=2 eV, J=0, central differences. Target dE_Hub/dU at fixed P (spinor functional at U +- eps, P of the two runs averaged) and measured dF_b/dU, Ry/eV:

| | step (eV) | measured D_b | target | b - target | A_U | b within A_U |
|---|---|---|---|---|---|---|
| KS | 0.05, 0.10 averaged | 0.05508311 | 0.05522711 | -1.44e-4 | 2.65e-3 (|D(h1)-D(h2)| 1.8e-6) | yes |
| RS (150 steps, beta 0.005) | 0.10 only | 0.05427951 | 0.05452495 | -2.45e-4 | 4.80e-3 (no step-difference term) | yes |

Candidate c (spectral, 0.10 eV only): KS D_c -0.08237415, miss -1.401e-1 = 52.8 A_U; RS D_c -0.08599122, miss -1.434e-1 = 29.9 A_U (gate: more than 10 A_U). Last-20-step spread of candidate b in the RS runs 1.3e-8 to 6.1e-8 Ry; all KS runs converged with Diff at most 2.5e-8 (W1) and 1.4e-8 (U, J).

T2(c) J derivative, KS, U=2 eV, J0=0.5 eV, local, step 0.10 eV only: D_b -0.06388335, target -0.06348434, miss -3.99e-4, A_J 3.10e-3 (no step-difference term). Not run for lack of KS budget (14 runs used): the 0.05 eV step, and candidate c at J. The gate on the average of two step sizes and the c miss test for J are therefore uncovered, and nothing is routed.

The driver's E_Hub at the sidecar P equals the logged E_Hub to 1.3e-15 (KS) and 4e-17 (RS): by construction (same functional, same P); it checks the sidecar orientation only.
