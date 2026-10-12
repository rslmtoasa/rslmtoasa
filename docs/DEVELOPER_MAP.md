# RS-LMTO-ASA Developer Map

A cheat-sheet for humans and AI assistants: for each workflow, the entry
point, the classes involved, the hot kernels, and where new functionality
plugs in. Every pointer below was checked against the code at the time of
writing — if something has moved, trust the code and fix this file, not the
other way around.

For the stabilized core architecture, module/submodule ownership, protected
atomic LMTO region, test taxonomy, fixture vocabulary, and the Phase-I gate
record, start with [the decision log](DECISIONS.md).
This map then provides the detailed workflow-level entry points and kernels.

For B1 frozen magnons, VAL-17 ([GBT harmonic and Goldstone report](validation/VAL-17_GBT_HARMONIC_GOLDSTONE.md))
records the current reciprocal FeCo route: a gapless, in-phase Gamma candidate
on the current deck, the reciprocal same-q collinear gauge reference used for
finite-q probes, and the open small-q stiffness mesh-convergence limitation.

---

## 1. Execution skeleton

```text
main.f90
  -> calculation(args%input)      ! constructor: restore_to_default + build_from_file
  -> calculation_obj%process()    ! source/calculation.f90
       pre_processing select-case  -> one pre_processing_* subroutine
       processing     select-case  -> processing_sd (only non-'none' value)
       post_processing select-case -> one post_processing_* subroutine
```

All three select-cases run in sequence and are independent: a case can set
any combination of `pre_processing`, `processing`, and `post_processing` in
its `&calculation` namelist (default `'none'` for all three). Each
`pre_processing_*`/`processing_*`/`post_processing_*` subroutine is a
**self-contained pipeline** — it constructs its own `control`/`lattice`/
`charge`/`hamiltonian`/... objects from `this%fname` rather than reusing
state from a sibling subroutine. Where one route needs the output of
another (e.g. `post_processing_orbital_modern` reading a converged
`Fe_out.nml`), that happens through the filesystem (`atoms%database`), not
shared Fortran state.

| `&calculation` key | Value | Routine | Notes |
|---|---|---|---|
| `pre_processing` | `'bravais'` | `pre_processing_bravais` | Bulk SCF (`calctype='B'`). Runs `self%run()`. |
| `pre_processing` | `'buildsurf'` | `pre_processing_buildsurf` | Surface SCF (`calctype='S'`). `build_surf_full()`, no `newclu()`. |
| `pre_processing` | `'newclubulk'` | `pre_processing_newclubulk` | Impurity-in-bulk (`calctype='I'`, bulk host). `newclu()`, no `build_surf_full()`. |
| `pre_processing` | `'newclusurf'` | `pre_processing_newclusurf` | Impurity-in-surface. `build_surf_full()` + `newclu()`. **No example/test input anywhere uses this route** (see `tests/KNOWN_ISSUES.md`). |
| `pre_processing` | `'buildinterface'` | `pre_processing_buildinterface` | Two-sided layered/interface SCF (`calctype='L'`, B7.5). `build_interface_full()` (region A \| active \| region B), then `surfmat()` (kernel reused unchanged) with its one-sided registry overwritten by `charge%build_interface_registry()`. Per-iteration Madelung update is `charge%interfacepot`, not `surfpot` (`self.f90` dispatch). `buildsurf` itself is untouched and remains the permanent one-sided regression oracle. |
| `processing` | `'sd'` | `processing_sd` | Spin dynamics. Rebuilds its consumer stack from the selected pre-processing route; the Depondt predictor/corrector now performs a predictor electronic refresh, corrected update, and post-correction electronic refresh. `Val13AbInitioSpinDynamics` validates the deterministic one-site bcc-Fe zero-torque loop; `Example_impurity_B2FeCo_sd_smoke` covers the production impurity output path. Broader dynamics remain out of scope. |
| `post_processing` | `'exchange'` | `post_processing_exchange` | Real-space intersite J_ij/D_ij; optional `do_damping=T` evaluates the route-agnostic Gilbert tensor and optional `do_inertia=T` emits the experimental raw magnetic-inertia diagnostic. Both consume the canonical Green functions filled by `gf_route`. |
| `post_processing` | `'exchange_q'` | — | **Removed in Stage 0b; the value stops at startup.** Campaign-era reciprocal q-resolved exchange; it `use`d `linear_response_mod`, so it was removed with it (tag `lr-campaign-archive-2026-10`). Not an LKAG oracle: use `'exchange'`. |
| `post_processing` | `'exchange_p2rs'` | `post_processing_exchange_p2rs` | Same, Hamiltonian sourced from a PAOFLOW-format import instead of `build_bulkham()`. |
| `post_processing` | `'conductivity'` | `post_processing_conductivity` | Real-space conductivity tensor. |
| `post_processing` | `'conductivity_p2rs'` | `post_processing_conductivity_p2rs` | Same, PAOFLOW-imported Hamiltonian. |
| `post_processing` | `'paoflow2rs'` | `post_processing_paoflow2rs` | Base PAOFLOW import path: builds the lattice normally but replaces `hamiltonian%build_bulkham()` with `hamiltonian%build_from_paoflow_opt()` (reads `paoham.dat`). |
| `post_processing` | `'orbital_modern'` | `post_processing_orbital_modern` | Orbital moments via `recursion%chebyshev_orbital_mod()` — always Chebyshev, regardless of `control%recur`. Loops over *every* atom in the cluster (`this%lattice%kk`), so runtime scales with cluster size, not recursion depth. |
| `post_processing` | `'band_structure'` / `'density_of_states'` | `post_processing_band_structure` / `_density_of_states` | k-space route via `reciprocal_mod`, not the real-space recursion machinery at all. |
| `post_processing` | `'fermi_surface'` | `post_processing_fermi_surface` | Writes the dense k-space eigensystem for Fermi-surface analysis (`reciprocal_export.f90`) on the `&reciprocal` `fs_nk1..3` mesh, always the full BZ. Its frozen-state setup (database potential → `hamiltonian` → `reciprocal`, `fermi_level` from `&energy`) is the template for new k-space post-processing routes. |
| `post_processing` | `'spin_response'` | `post_processing_spin_response` → `spin_response%prepare/interaction/run/write_state` (`spin_response.f90`) | Transverse spin response chi0(q, omega), U and Dyson for Juelich-d amplitudes on the frozen k-space state, q list from `&spin_response` (`q_list` or `q_file`), several sites per cell. See "Transverse spin response (B11)" below. |
| `post_processing` | `'bsf'` | `post_processing_bsf` → `reciprocal%calculate_bsf` (`reciprocal_bsf.f90`) | Bloch spectral function A(k,E) = −1/π Im Tr G(k,E+iη) along the canonical spglib k-path (milestone B3). Consumes the B2 engine's `dyson_kspace_inverse` per (k,E) (Σ=0 ⇒ backend E; Σ-ready for CPA/DMFT). η = `&reciprocal` green_eta, E grid = n_energy_points/dos_energy_min,max, path = `&kpath` nk_per_segment. Writes `bsf.dat` (total/up/down) + `bsf_bands.dat` overlay. Partial-trace convention in `bsf_kernel.f90` (`bsf_spectral_trace`). |
| `post_processing` | `'kspace_green'` | `post_processing_kspace_green` | B2 validation driver: fills `green%gij` via recursion then via the k-space engine (`reciprocal%fill_green`, backend E + D≡E check) and cross-checks on-site DOS / m_z. Report-only. |
| `post_processing` | `'frozen_magnon'` | `post_processing_frozen_magnon` | Sweeps `hamiltonian%q_ss` over a `&frozen_magnon` q-list, preferably from `q_file` (`q_coordinates='cartesian'` for `2*pi/alat` Cartesian components or `'direct'` for reciprocal-lattice coordinates), writing total energy, band energy, per-sublattice moment magnitude, and `omega(q)` to `frozen_magnon.dat`. `mode='mft'` (default) converges SCF once at the reference point, reuses that potential for a single-iteration band-energy pass at every other q, and computes `omega` from band-energy differences; `mode='scf'` re-converges at every q and computes `omega` from total-energy differences. `branch_mode='auto'` builds multi-sublattice magnon branches in `frozen_magnon_branches.dat`/`frozen_magnon_modes.dat` via the direct GBT frozen-magnon method (second derivatives of the force-theorem band-energy surface w.r.t. sublattice cone angles; Essenberger PRB 84, 174425 Eq. 26). **Single-sublattice is validated; for the multi-sublattice acoustic branch see VAL-17 above, `tests/KNOWN_ISSUES.md`, and B1 in `docs/ROADMAP.md`.** See `docs/DECISIONS.md` for the archived campaign record. |
| `post_processing` | `'linear_response'` | — | **Removed in Stage 0b; the value stops at startup.** Removed campaign (tag `lr-campaign-archive-2026-10`; see `docs/DECISIONS.md`). Transverse spin response is rebuilt as B11 (`'spin_response'`, `docs/ROADMAP.md`). |
| `post_processing` | `'pauli_projection'` | — | Campaign-era (LR-02N); ran inside the bravais SCF handoff. **Removed in Stage 0b; the value stops at startup.** |

Routines live in `calculation.f90` except the `bravais`, `buildsurf`,
`newclubulk` and `newclusurf` pre-processing routines
(`calculation_preprocessing.f90`) and the k-space post-processing routines
`band_structure`, `density_of_states`, `bsf`, `fermi_surface`,
`kspace_green` and `spin_response` (`calculation_reciprocal.f90`).

`exchange_p2rs`/`conductivity_p2rs`/`paoflow2rs` all funnel through the
shared helper `prepare_post_processing_stack(this, use_paoflow, ...)` in
`calculation.f90`, which branches on `control%calctype` for the non-paoflow
path and always uses the `'B'`-style bulk sequence when `use_paoflow=.true.`
(PAOFLOW import is bulk-only — `'S'`/`'I'` call `g_logger%fatal`).

---

## 2. Per-workflow chains (entry → orchestration → kernel)

### Real-space SCF
```text
calculation%pre_processing_bravais / _buildsurf / _newclubulk / _newclusurf
  -> self%run()                         [self.f90 — legacy SCF loop core, off-limits per ground rule 10]
       -> hamiltonian_build.f90         (build_bulkham / build_locham, +U/+V via hamiltonian_hubbard.f90)
       -> recursion(hamiltonian_obj, energy_obj, sparse_obj)   [recursion.f90]
            select case (control%recur)
              case ('lanczos')   -> recursion%recur()    -> crecal()
              case ('block')     -> recursion%recur_b()  -> crecal_b() / hop_b[_hoh]()
              case ('chebyshev') -> recursion%chebyshev_recur()
                   -> gpu_plugin_ready()? -> gpu_backend%chebyshev_moments()  [CUDA, bypasses cheb_backend]
                   -> else: cheb_moments_cpu() dispatcher (recursion.f90)
                        -> chebyshev_fast.f90: cheb_moments_fast /
                           _batched / _mkl_batch / _mkl_sparse, or legacy path
       -> green.f90               (sgreen / block_green / chebyshev_dos_dispatch per recur)
       -> charge.f90 / mix.f90    (charge transfer, mixing, convergence check)
```

### GBT SCF frame contract (WP12)

`self%run_dos()` has two density producers: the real-space Green-function
route and the optional reciprocal eigenvector route. Both feed the shared
rotating-frame density contract before radial up/down projection. At the
magnetic-state boundary, `density_policy='constrained_spiral'` keeps
`potential%mom` as a fixed frame marker; only charge and longitudinal radial
channels are mixed, while the measured transverse vector remains a residual
for the constraint controller. `relaxed_reference` retains full Cartesian
mixing for ordinary SCF but is explicitly unsupported with `gbt_single_q`
until a moving-reference GBT Hamiltonian contract is audited.

Set `gbt_scf_diagnostics = .true.` in `&control` to emit
`gbt_scf_diagnostics.dat`. Each row joins incoming frame/moment/radial
channels, the solver-frame density matrix and moment, target/field,
magnetic-mixer state, next radial channels and marker, exchange parameters,
physical energy, and separated constraint diagnostics. The trace is opt-in
and does not alter the equations.

### k-space (reciprocal)
```text
post_processing_{band_structure,density_of_states,bsf,fermi_surface,kspace_green,spin_response}
  (calculation_reciprocal.f90)
  -> reciprocal_mod (reciprocal.f90 + submodules reciprocal_*.f90)
       -> reciprocal_fourier.f90: k-space H(k) built via Fourier transform of
          real-space hoppings; branches on hamiltonian%ccor_2c and
          hamiltonian%hoh (both bypass the "plain" Fourier sum).
          calculate_eigenpairs_at_kpoints diagonalizes at arbitrary k.
       -> reciprocal_bands.f90 / reciprocal_dos.f90: diagonalization,
          tetrahedron or gaussian DOS integration (&reciprocal dos_method)
       -> reciprocal_projection.f90: orbital-projected DOS, band moments
       -> reciprocal_occupations.f90: Fermi level, Fermi-Dirac occupations,
          electron count, band energy
       -> reciprocal_green.f90: B2 Green-function engine (fill_green,
          backends 'lehmann' (E) and 'dyson' (D))
       -> reciprocal_bsf.f90, reciprocal_moments.f90 (KPM conductivity
          moments), reciprocal_spin_density.f90, reciprocal_backend.f90
          (execution backends), reciprocal_export.f90 (eigensystem export)
  -> kpoint_workset.f90: owned fractional k-points and weights; shifted(q)
     builds the k+q set (weights valid for full-BZ meshes only)
```
Note the k-space route does not go through `recursion`/`green`/`chebyshev_*`
at all — it is a parallel diagonalization-based path, not a recursion-based
one. `etot`/`ws_r`/`vmad` on these cases usually come from a *preserved*
restart potential (`preserve_outputs: true` cases), not from a fresh SCF
loop — see `Example_density_of_states_bccFe_ccor_2c` in
`tests/postproc/cases.json` for why its checks target `dos_kspace.dat`
instead.

### Exchange (J_ij), conductivity, orbital moments — intersite chain
```text
calculation.f90: run_intersite_moments(control_obj, recursion_obj)
  -> recursion%recur_b_ij() / recursion%chebyshev_recur_ij()   [recursion.f90]
  -> green%calculate_intersite_gf_core() / green%calculate_intersite_gf_eta()
       -> green%bgreen()   [green.f90 — per-pair, per-energy Green's function terminator]
  -> exchange.f90 (calculate_exchange / _twoindex / _gauss_legendre)
     or conductivity.f90 (calculate_gamma_nm / calculate_conductivity_tensor)
```
The energy integrals in `exchange.f90` (`simpson_f`) use the electronic
temperature `&energy temperature` (default 0 K, where the behaviour is
bit-identical to before the key existed). At finite T, `e_mesh` in `energy.f90`
warns unless the grid reaches 10 kT above E_F with a spacing of at most kT/2.
Real-space LKAG J_ij are robust for near shells; beyond about 3a they are not
converged in cluster size, recursion depth or cutoff (`docs/DECISIONS.md`, B11
Stage 3b), which limits any small-q stiffness taken from them.

### Spin dynamics
```text
calculation%processing_sd()
  -> prepare_post_processing_stack(..., preprocessing_route=this%pre_processing)
       rebuilds the selected route's consumer stack from the saved state
  -> self(bands_obj, mix_obj)
  -> spin_dynamics(self_obj)          [spin_dynamics.f90]
       -> sd_run(): per step, self%run() -> bands%calculate_magnetic_torques()
          -> asd_pred_euler / asd_corr_euler (Landau-Lifshitz-Gilbert integrator)
       -> include_codes/abspinlib/abSpinlib.f90   (lower-level ASD numerics)
```
`Example_bulk_bccFe_sd_smoke` checks that a one-step production trajectory is
emitted. The deterministic one-site zero-torque loop is validated by
[`VAL-13`](validation/VAL-13_AB_INITIO_SPIN_DYNAMICS.md), using the LMTO field/torque
seam from [`VAL-12`](validation/VAL-12_LMTO_MAGNETIC_FIELDS_TORQUES.md) and the
Depondt evidence from [`VAL-11`](validation/VAL-11_ABSPINLIB_INTEGRATORS.md). The
impurity output path has a production smoke, but broader impurity physics and
general trajectories remain outside the scoped claim.

### PAOFLOW import
```text
hamiltonian%build_from_paoflow_opt()   [hamiltonian_paoflow_io.f90]
  reads paoham.dat, format: "idxi idxj idxk global_orb_i global_orb_j Re[eV] Im[eV]"
  (the "legacy7" record shape — see write_rs_tb_records)
```
The *export* side (`hamiltonian%export_rs_tb_all` / `rs2pao`, selected via
`hamiltonian%export` in `calculation.f90`'s `pre_processing_bravais`,
select-case near line 641) writes files in the same format
`build_from_paoflow_opt` reads — `export='python'` produces
`<base>_paoham.dat` in exactly the legacy7 layout. This is how
`tests/postproc/cases/paoflow/*/paoham.dat` were generated (see
`tests/postproc/README.md` for provenance).

---

## 3. Kernel inventory

- **`ham_apply_t`** (`chebyshev_fast.f90`): the operator-apply callback
  interface — `subroutine ham_apply_t(x1, x0, y, alpha, beta, ld, nb)`
  computing `y = alpha*H*x1 + beta*x0`. Every fast Chebyshev kernel
  (`cheb_moments_fast`, `_batched`, `_mkl_batch`, `_mkl_sparse`) is written
  against this one interface via `cheb_three_term_moments`, so the recursion
  itself (three-term Chebyshev recurrence) is written once and the backend
  only supplies `apply_h`.
- **`cheb_cache_t`** (`chebyshev_fast.f90`): holds the scaled/BSR/orthogonalized
  Hamiltonian, velocity operators, and MKL pointer arrays, each behind its
  own `cheb_cache_fingerprint_t` — a cheap hash of the inputs that last built
  that cache slot. A cache slot is rebuilt only when its fingerprint changes,
  so repeated calls with the same Hamiltonian (e.g. across energy points)
  don't redo the scale/BSR-pack/orthogonalize work. Single-rank assumption:
  the cache is a module-level singleton, not distributed across MPI ranks.
- **`hsweep_sp` / `happly_sp`** (`chebyshev_fast.f90`): single-precision (`sp`)
  fused OpenMP matvec sweep over all cluster atoms — the "no per-atom scalar
  loop, no lightcone masking" trick described in the module's header comment.
  `happly_sp` layers the hoh two-sweep (bond + onsite + hop) on top of
  `hsweep_sp` when `do_hoh` is set.
- **Block-Lanczos fast kernels** (`haydock_fast.f90`): `block_lanczos_fast`
  dispatching to `block_lanczos_sp`/`_dp` — the `recur_b`/`crecal_b` hot path
  for `control%recur='block'`.
- **k-space pure kernels** (`lehmann_kernel.f90`, `dyson_kernel.f90`,
  `bsf_kernel.f90`, `moment_kernel.f90`): plain-array procedures with no
  object state, called by the reciprocal submodules. New k-space numerics
  should follow this shape so analytic tests can call them directly.
- **CUDA plugin surface**: `source/rsrec_cuda_plugin.f90` (Fortran
  `iso_c_binding` wrapper, type `rsrec_cuda_backend`) ↔
  `source/cuda/rsrec_gpu.cu` (device kernels) ↔ `source/cuda/rsrec_cuda.h`
  (C API). `control%gpu_plugin=.true.` makes `chebyshev_recur()` (and the
  intersite/transport variants in `recursion_transport.f90`) dispatch
  straight to `gpu_backend%chebyshev_moments()` — this **bypasses
  `cheb_backend` entirely**, regardless of what CPU backend the case also
  specifies. Layout convention: site-major, and the plugin's internal
  arithmetic is fp32/fp64 mixed the same way `chebyshev_fast.f90`'s CPU
  kernels are (see that module's header). Cache/fingerprint lifecycle
  mirrors `cheb_cache_t`: `ensure_context`/`upload_bsr` are only re-run when
  the uploaded Hamiltonian's fingerprint changes.
- **CUDA reciprocal backend** (`source/cuda/reciprocal_cuda.cpp`,
  `reciprocal_lehmann.cu`, built with `ENABLE_CUDA_RECIPROCAL`, off by default):
  a batched Hermitian eigensolver and a strict-Lehmann contraction behind the
  reciprocal execution-backend interface, with the GPU tests
  `ReciprocalCuda{ArbitraryK,Lifecycle,Eigensolver,Lehmann}`. Not yet audited on
  `fable_v4b`, and not used by the spin response (B4).
- **Spin-response kernels** (`spin_response.f90`): plain-array procedures
  (amplitudes, `accumulate_chi0`, `u_juelich`, `solve_dyson`, ...) that the
  analytic tests call directly; see the B11 recipe below.

---

## 4. "Where to add X" recipes

### Porting `strux_lib` to GPU
Structure-constant generation lives in `include_codes/strux_lib/` (a dozen
files: `strux_api.f90`, `strux_tb.f90`, `strux_bessel.f90`,
`strux_solver_assembly.f90`, ...), called from `lattice_strux.f90` when
`control%strux_backend == 'strux_lib'` (required whenever `lmax > 2`; the
legacy backend only supports up to `spd`). The boundary is: `strux_lib`
consumes lattice geometry (positions, `lmax`) and produces structure-constant
matrices consumed by `hamiltonian_build.f90`. Mirror the existing plugin
pattern for a GPU port: a C API + `iso_c_binding` wrapper module + CPU
fallback, exactly as `rsrec_cuda_plugin.f90` does for the recursion kernels
— do not invent a second GPU convention.

### Transverse spin response (B11)

The `linear_response` code was removed from `fable_v4b` (tag
`lr-campaign-archive-2026-10`); its dispatch stops at startup, and its docs,
numbers and tests are not references.

**In the tree:** `source/spin_response.f90` (`spin_response_mod`, about 1,100
lines) with its `&spin_response` namelist in
`source/include_codes/namelists/spin_response.f90`, run by `post_processing =
'spin_response'` (`post_processing_spin_response` in `calculation_reciprocal.f90`).
It works on the frozen database state, never an SCF. The model is collinear. The response space is indexed by kappa = (site, channel),
channel = s, p, d, f from `&spin_response channels` (default `'d'`; must contain d, letters
beyond the basis lmax are fatal), with one static U per site acting on the d channel only.

- `prepare`: full-BZ mesh (symmetry reduction and time reversal switched off,
  and logged), E_F policy and electron-count check, channel amplitudes (coefficient
  "Mills" and Juelich, `enu` and `ppar` of each l), moments per (site, channel) in
  `occ_mills`/`occ_juelich` and for d in `d_occ_*`, spin relabelling when the d
  moment is negative.
- `interaction`: U_Juelich from the Goldstone condition on the static chi0 at
  eta = 0; U_Mills from the d band centres (`center_band`), with the Mills
  residual and implied gap; both printed. It also builds `kmat`, the matrix over
  kappa that `run` uses: U on the d rows (`kernel = 'goldstone'`, default), or
  U = -K from the ALSDA matrix (`'alsda'`; `'alsda_scaled'` times the factor
  lambda, least squares along M of chi0(0,0) lambda U M = -M). With an ALSDA
  kernel `_state.dat` prints lambda, the q = 0 Rayleigh eigenvalue of
  1 + chi0 U along M for the unscaled kernel, and K. LDA only (`txc = 1`).
  The default `kernel = 'goldstone'` is the production kernel; `'alsda'` and
  `'alsda_scaled'` are experimental (known gaps in `tests/KNOWN_ISSUES.md`).
- `run`: chi0 for all q in one pass over the k-chunks (`bare_response_list`: the
  k endpoint once per chunk, then each q's k+q endpoint), then per q the Dyson
  solve with K over kappa (U on the d rows, zero elsewhere), tr chi0, tr chi, tr L
  and the smallest |eig(I + chi0 K)|; peak for one site, crossing for one site
  and one channel (`n/a` otherwise). Output: `<prefix>_q<NNN>.dat`,
  `<prefix>_dispersion.dat` (last column: smallest real part of the eigenvalues of
  I + chi0 U at the first omega, the signed delta_q for `n_omega = 1`, `omega = 0`,
  `eta = 0`) and `<prefix>_state.dat`. With `n_omega = 1`,
  `omega = 0` and `eta = 0` it gives the static chi0(q, 0).
- Adiabatic frequencies, `<prefix>_adiabatic.dat`: eigenvalues of
  Omega = -N (chi0(q, 0)^-1 + K), N = diag of the signed moment of each row kappa
  (`adiabatic_frequencies`; static chi0 from the same eigenpairs through the
  optional `chi0_static` of `bare_response_list`). Global spin frame: an
  antiparallel site has m < 0 and a negative eigenvalue; no per-site frame is
  used (a local-frame vertex on a down site makes chi0_ij zero between up and
  down sites). Rows without interaction (K = 0, the s and p channels) enter with
  their own moment and have no physical meaning as frequencies.
- `accumulate_chi0` lists the band pairs with a nonzero vertex (exact zeros of the
  spin-block amplitudes) and f_n - f_m /= 0, then does one ZGEMM per block of 32
  omega points, blocks in parallel (OpenMP), so the output is byte-identical for
  any thread count when `MKL_NUM_THREADS = 1`. Static runs are serial; for many
  q, split the q-list over single-thread processes.
- Both endpoints come from `calculate_eigenpairs_at_kpoints` through
  `eigenpairs_chunk`, which applies the site-dependent Bloch phase
  exp(-2 pi i G.tau_i) for a k+q folded by G, needed for cells with several
  sites (`tests/KNOWN_ISSUES.md`, resolved entry). Occupations and constants come
  from `reciprocal_mod` (`fermi_dirac_occupation`, `kB_Ry_per_K`,
  `find_fermi_level_from_eigenvalues`); `enu`, `ppar` and `center_band` from the
  potential.
- Kernels are public procedures on plain arrays: `d_amplitudes_coefficient`,
  `d_amplitudes_juelich`, `accumulate_chi0`, `u_juelich`, `u_mills`,
  `mills_residual`, `interaction_values`, `radial_kernel_matrix`, `solve_dyson`,
  `spectral_trace`, `pole_estimates`.
- ALSDA kernel, K_ab = (1/4 pi) int dr P_a P_b K(r)/r^2 per site (spherical
  vertex, L = 0), K = [V_xc,up - V_xc,dn]/m with the core in m, P = SR overlap of
  the up and down radial functions. A response run has no radial data, so
  `source/alsda_radial.f90` (`alsda_radial_mod`, outside `spin_response*`)
  re-runs `self%lmtst` on a copy of each site (the lattice potential parameters
  are untouched), builds phi at E_nu (coefficient amplitudes) or the frozen
  projector (Juelich amplitudes), and `alsda_kernel` hands the blocks to
  `spin_response%set_alsda_kernel`; the driver calls it after `prepare`.
- q is given in direct coordinates of the code's primitive basis for
  `crystal_sym`; check the Cartesian columns of `_dispersion.dat` when reusing a
  q-path from another code.

**References and decks:** developer-owned, under `tests/spin_response/oracles/`:
the analytic models A1-A4 and the radial kernel integrals R1, `tolerances.nml` (whose keys are declared once in
`tests/unit/spin_response_tolerances.f90`), and the LKAG + LSWT magnon curves
`fe_lswt.dat` and `ni_lswt.dat` from the real-space `exchange.f90` (see
[`EXCHANGE_VALIDATION_MAP.md`](validation/EXCHANGE_VALIDATION_MAP.md)); their
small-q end is not converged in the J(r) cutoff. The decks are
`tests/spin_response/fe_bcc` and `tests/spin_response/ni_fcc`; the Stage 2 and
3a records are at tag `spin-response-stage3a`.

**Where to add:** a 2 x 2 spin structure per channel enlarges kappa and the vertex builder
in `accumulate_chi0`. New kernels follow the
plain-array pattern above, so analytic tests can call them. The radial module
`radial_ground_state` is read only through `alsda_radial.f90`;
`lmto_radial_augmentation` is not used.

### Lehmann-representation Green's functions
Implemented (B2): `reciprocal_green.f90` provides `reciprocal%fill_green`
with backend `'lehmann'` (E, Σ = 0) or `'dyson'` (D), using the pure kernels
in `lehmann_kernel.f90` and `dyson_kernel.f90`. Results land in the same
`green` containers the real-space routes fill, so exchange, conductivity and
damping consume either route. New Green-function consumers on the k-space
side belong in the reciprocal family, not in `green.f90`.

---

## 5. Testing map

| Workflow | Suite / label | Case-file location |
|---|---|---|
| Regression backend matrix (Chebyshev CPU backends × hoh × ccor_2c, block sp/dp) | `regression;backend` | `tests/regression/cases.json` |
| Bulk/surface/impurity SCF, all recursion modes | `example;scf` | `tests/scf/cases.json` |
| Exchange, conductivity, bands, DOS, orbital moments, PAOFLOW import | `example;postproc` | `tests/postproc/cases.json` |
| Serial-vs-MPI consistency (representative subset) | `launch_serial` / `launch_mpi` | `"launch_modes": ["serial","mpi"]` field, either `cases.json` |
| Quick tier (PR-fast subset; `ctest -L '^quick$'` before each commit) | `quick` | `"quick": true` field in either example `cases.json`; `quick` in the `labels` of a `tests/regression/cases.json` case; `LABELS` in `CMakeLists.txt` for unit and tooling tests |
| MKL kernels (`mkl_batch`/`mkl_sparse`) | gated by `requires_cmake_option: ENABLE_MKL_KERNELS` | regression + `Example_bulk_bccFe_nsp2_chebyshev_mkl_batch` in `scf/cases.json` |
| CUDA plugin | compile-only in CI (`cuda_compile` job); real-GPU consistency via `tests/run_gpu_matrix.sh` | n/a (manual, off-CI) |
| Standalone Fortran unit tests | `unit` (+ topic labels), built with `-DRUN_UNIT_TESTS=ON` | `tests/unit/test_*.f90`, registered with `add_fortran_unit_test` in `CMakeLists.txt` |
| Spin response: analytic kernels (A1-A3, E1, E3, E4), bcc-Fe components C1-C4 (incl. C4.3) and S2 | `unit;spin-response;quick` | `tests/unit/test_spin_response_*.f90`; the Fe programs run in the scratch deck made by `SpinResponseFeDriver` (`tests/spin_response/fe_scratch.py`) |
| Spin response: channels s, p, d (`UnitSpinResponseChannelsKernel` synthetic; `UnitSpinResponseChannels` on the Fe scratch deck, needs `channel_d_rel` in `tolerances.nml`) | `unit;spin-response;quick` | `tests/unit/test_spin_response_channels*.f90`; the spd driver runs are made by `fe_scratch.py` |
| Spin response: S3, three-atom bcc Fe cell against the one-atom cell (`SpinResponseS3`) | `unit;spin-response;quick` | `tests/spin_response/s3_supercell.py` |
| Spin response: 1 against 4 OpenMP threads, byte-identical (`SpinResponseOmpIdentity`) | `unit;spin-response;quick` | `tests/spin_response/omp_identity.py` |
| Spin response: 24^3 reproduction against frozen outputs, Fe and Ni (`SpinResponseBaselineRepro`, `SpinResponseNiBaselineRepro`) | `regression;spin-response` | `tests/spin_response/baseline_repro.py`; frozen files in `tests/spin_response/baseline/` and `baseline_ni/` (`--freeze` rewrites them) |
| Exchange energy integration at finite temperature (`UnitSimpsonFTemperature`) | `unit;quick` | `tests/unit/test_simpson_f_temperature.f90` |
| Validation scripts | `validation` | `tests/validation/*.py` (spin response: `tests/validation/spin_response_ni/`), records in `docs/validation/` |

**Adding a case:** see `tests/scf/README.md` / `tests/postproc/README.md` for
the full case-file format. Short version: add an `input.nml` (+ any
supporting data files) under `cases/<workflow>/<system>/`, add an entry to
the suite's `cases.json`, generate its reference with
`tests/generate_references.py --case <TestName>`, commit both.

**Regenerating references:** see `tests/README.md` — run
`generate_references.py` against a trusted build, review the diff of
reference files like any other code change, never regenerate just to "fix"
an unexplained failure.

---

## 6. Conventions box

- **KISS.** No argument-validation ladders, wrapper layers, or "just in
  case" guards on internal routines — checks belong at true boundaries
  (namelist parsing, file I/O, plugin/library availability).
- **Class-based architecture stays.** Derived types with constructors,
  `destructor`, `restore_to_default`, and type-bound procedures. New code
  follows the same pattern; submodules split *implementation*, never
  flatten types into bare module procedures.
- **Legacy files are off-limits.** `self.f90` and `symbolic_atom.f90` contain
  physics-dense legacy code scheduled for a separate audit — files may be
  edited around the edges (a `use` statement, a comment), but their bodies
  are not refactoring targets at this stage.
- **One task, one commit.** Single-line commit messages, no co-author
  trailer (`CLAUDE.md`, rule 2).
- **Expected values never come from the code under test.** See
  "Verification integrity" in `CLAUDE.md`.
- **Bit-level behavior is the contract.** Run the quick tier
  (`ctest -L '^quick$'`) before each commit, and
  `ctest -LE 'validation|performance'` in a fresh build at the end of each stage
  (`CLAUDE.md`, rule 1).

---

## 7. Python module outlook (P11 — design awareness, not implementation)

A future goal is exposing this code as a Python module for notebooks,
high-throughput screening, and tutorials. Nothing below is implemented; it
is a standing lens applied while doing Phase 2 work, recorded here so a
future binding effort starts with a blocker list instead of rediscovering
these from scratch.

**Recommended binding route:** given the class-based design, `f90wrap` is
the natural fit — it wraps derived types with type-bound procedures into
Python classes directly. Fallback: a thin explicit `iso_c_binding` API layer
(the `rsrec_cuda_plugin` C-API pattern, in reverse) wrapped with
ctypes/cffi — more work, but zero toolchain magic. A pybind11 route would
only make sense via that same C API.

**Python-binding blockers found so far:**

1. **Hard-coded input/output filenames.** Every `pre_processing_*` /
   `processing_*` / `post_processing_*` subroutine reads `this%fname`
   (an `input.nml` path) and writes fixed filenames into the current working
   directory (`Fe_out.nml`, `totaldos.out`, `jij.out`, ...) rather than
   returning data or accepting an in-memory config object. A Python session
   driving multiple calculations would need a working-directory-per-call
   convention (already how the test suites use it) or a refactor to accept
   config objects directly — out of scope for Phase 2.
2. **`g_logger%fatal` calls `stop`/`error stop`.** Every subroutine here
   calls `g_logger%fatal(...)` on bad input (see `check_pre_processing`,
   `check_processing`, `check_post_processing`, the PAOFLOW file-not-found
   checks, the `ENABLE_MKL_KERNELS` guard in `chebyshev_fast.f90`, ...) —
   this terminates the whole process. A Python session must not be killed
   by one bad call; the fatal path would need to become an exception at the
   binding boundary (or `g_logger%fatal` itself would need a non-fatal mode
   for library use — a bigger design decision, not attempted here).
3. **Each pipeline subroutine rebuilds everything from scratch.** There is
   no "give me a `lattice`/`charge`/`hamiltonian` object and let me call
   `post_processing_exchange`-equivalent logic on it" entry point — the
   `pre_processing_*`/`post_processing_*` subroutines are monoliths that
   construct control→lattice→charge→hamiltonian→recursion→green→bands→self
   internally. A binding would want these decomposed into reusable
   construction + a thin dispatch layer, which is exactly the kind of
   structural change Phase 1/2 ground rule 10 and the "no further
   structural refactoring this phase" rule put out of scope for now.
4. **`mpi_mod`'s MPI-range partition is bare module-level mutable state,
   read directly (not through `this`) by most of the physics core.**
   `start_atom`, `end_atom`, `atoms_per_process`, and `g2l_map` are set once
   by `get_mpi_range` and then read as free module variables in
   `bands.f90`, `electrostatics_multipole.f90`, `exchange.f90`, `green.f90`,
   `reciprocal.f90`, `reciprocal_green.f90`, `reciprocal_lifecycle.f90`,
   `recursion_haydock.f90`, `recursion_chebyshev.f90`,
   `recursion_transport.f90`, and `self.f90` — found while writing the
   `green.f90` P7 docstrings (2026-08-09), where nearly every routine's
   atom loop is bounded by `start_atom`/`end_atom` and indexed through
   `g2l_map`. A single-process Python session could call `get_mpi_range`
   with a 1-rank partition to make this safe, but multiple independent
   `calculation` objects in the same interpreter would still corrupt each
   other's partition — this is a real blocker, not just a caches-shaped one
   like `cheb_cache_t`/the CUDA context singleton (those are already
   fingerprinted and per-instance-safe). Not fixed here: moving the
   partition onto `this` touches every one of those 11 files' call sites,
   which is the structural change Phase 2 defers.

**Results as data:** the physically meaningful outputs per workflow are
documented on each object by the now-complete P7 docstring pass — e.g.
`self`/`bands`'s moment arrays, `green`'s `g0`/`gij`/`gji` (+ the
`Gi{x,y,z}`/`Gj{x,y,z}` spin-decomposed and `G00ij`/`G01ij`/... two-index
families), `exchange`'s J_ij/D_ij, `conductivity`'s tensor,
`reciprocal_dos`'s DOS arrays. That documentation is the binding surface
specification for a future Python module.
