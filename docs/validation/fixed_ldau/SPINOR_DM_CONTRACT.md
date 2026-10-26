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
