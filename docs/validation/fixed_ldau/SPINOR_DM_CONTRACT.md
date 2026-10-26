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
