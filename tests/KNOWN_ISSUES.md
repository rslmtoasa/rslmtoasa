# Known issues

Bugs found during other work and recorded here rather than fixed in place (see "Scope
discipline" in `CLAUDE.md`). This file keeps open items only. Each entry is a candidate
for a future bug-fix task and separates what was verified from what is suspected.
Entries that were fixed, resolved or superseded were dropped; they are recoverable from
the archive tags `archive/fable_v4b-2026-11-05`, `archive/fable_v4b_u-2026-11-05`,
`archive/fable_v4b_gbt-2026-11-05` and `archive/main-2026-11-05`. Campaign evidence for
fixed-U LDA+U is summarised in `docs/validation/fixed_ldau/VALIDATION_SUMMARY.md`.
Issues in the linear-response code removed in Stage 0b (`lr-campaign-archive-2026-10`)
are not tracked here; see `docs/DECISIONS.md`.

## Reference and tolerance questions

### Bitwise sub-case of `UnitLdaULocalMoments` fails with `-march=native` on FMA hardware — 2026-11-05
- **Verified:** sub-case `spinor_frame_kernel` (`test_hubbard_spinor_moments.f90`, check "(f) n = z vs the S4 z-block loop (bitwise)") gives
  2.2e-16 against 0 on x86_64 Linux, gfortran 13, OpenBLAS, `-DENABLE_MARCH_NATIVE=ON`. It gives 0.0 with `-DENABLE_MARCH_NATIVE=OFF` (the CI
  configuration, `.github/workflows/tests.yml`). The failure is the same on the unmodified `fable_v4c` (`9a87c09`) as test `UnitHubbardSpinorMoments`.
- **Suspected, not diagnosed:** fused multiply-add contraction in two production routines (`correct_d_shell_moments_spinor` and the S4 z-block
  loop) that sum in the same order but contract differently. The check is not loosened.

### LDA+U stored-reference rows: tolerances were set on the developer platform — 2026-11-05
- **Verified (x86_64 Linux, gfortran 13, OpenBLAS, 4 steps, single thread):** `imp_rs_soc0` (RS-host impurity, nsp=2 with SOC = 0 against nsp=1) differs by
  1.0e-6 in the 4-step signature against the tolerance 1e-7 (iteration 1: 0.0 native, 1.6e-8 CI configuration against 1e-9); the KS-host row
  `imp_ks_soc0` gives 1.7e-7 (native) and 5.1e-8 (CI configuration). `LdaURefRsImpurityNsp2` is therefore tier `validation`. The surface row
  `surf_soc0` agrees to 0.0.
- **Verified:** the committed SOC-on snapshot (developer platform) is matched to 6e-10 (KS) and 2e-6 (RS `rs_soc_u2`) against the tolerance 1e-10, so
  `ks_soc_u2_x` has a `linux-x86_64` variant (`tests/ldau/references/ks_soc_u2_x.linux-x86_64.json`, written from the CI configuration with a recorded reason).
- **Verified:** the stored references (`tests/ldau/references/`) are bound to the CI configuration (`-DENABLE_MARCH_NATIVE=OFF`; reference build `90f57c3` built the
  same way). A `-march=native` build on this FMA CPU deviates from them: `noU_rs` 1.2e-6, `U_rs` 7.6e-9, `ks_soc_u2_x` 6.0e-10, against a tolerance of 1e-10;
  `surf_soc0` agrees to 0.0. With the CI configuration the candidate matches the reference build to 0.0 on the no-U, `U_rs`, `U_imp` and `U_ks` rows.
- **Not diagnosed:** whether the 1e-6 impurity difference is only the documented local-mode transient amplification (see "RS +U transient") or a platform
  effect of the BLAS used. Owner decision needed: platform tolerance, or a different row.

### Long examples time out under `ctest -j4` — 2026-10-03
- **Observed:** `Example_orbital_modern_bccFe` (about 480 s serial), `Example_exchange_bccFe_hoh` and `Example_exchange_conductivity_fccPt_hoh` hit
  their 600-900 s limits under `ctest -j4`; serially they finish within the limits. Not investigated.

### Example-suite etot tolerance is 5.8e-4 to 0.37 Ry — 2026-10-04
- **Verified:** `_check_value` in `tests/run_test.py` fails a value only when both `abs_tol` and `rel_tol` are exceeded, so for etot of order 1e3 Ry
  the effective tolerance is `rel_tol * max(abs(ref), 1)`. Of 50 example cases in `tests/scf/cases.json` and `tests/postproc/cases.json`, 34 compare
  etot and all 34 are looser than 1e-5 Ry (diamond Si 5.8e-4, Fe cases 2.5e-3 to 2.5e-2, `Pt2MnGa` 0.37).
- **Consequence:** etot shifts of 5.8e-6 to 1.9e-5 Ry, which the `Regression_*` matrix detects (absolute 1e-6), are invisible to every example case. A
  reference shifted by +2e-6 passes.
- **Proposal, not done:** per-quantity absolute etot tolerances chosen from measured run-to-run and platform spread (or compare with "or" instead of
  "and"). The `Regression_*` and Triad tests use their own rules and are not affected.

### Triad J_ij and the recursion noise floor — 2026-10-04
- **Verified:** the recursion-route J_ij of the Triad deck (bcc Fe) reproduces only to about 5e-8 relative: a `lattice%alat` perturbation of 1e-12
  gives |dJ/J| = 5.1e-8 (J[1_335]), `-O0` against default differs by 5.0e-8 and 2.3e-8; lehmann and dyson agree at 1e-15. It does not scale with the
  perturbation, so it is rounding noise. Its source in the recursion is not diagnosed. A check on recursion J tighter than about 1e-6 relative would
  be below this floor.
- **Open accuracy question (not measured):** after the `simpson_f` Fermi-level fix (`source/exchange.f90` now passes the grid point `en%ene(enpt +
  1)`), with `enpt` = 141 the E_F point is the centre of a Simpson panel, with f = 0.5 there. Whether the integral up to E_F is accurate, and how it
  depends on the parity of `enpt` or on `channels_ldos`, was not measured. The fixed lehmann J[1_335] moved +0.4% (0.257612 to 0.255696); neither
  value is known to be the correct J.
- **Suspected, not tested:** the original residual dependence came from `energy_min + edel*i` (`source/energy.f90:214`) contracted to FMA at `-O2`;
  the regression tests were not rebuilt with `-ffp-contract=off`. Tests: `tests/regression/triad_bccFe_exchange` (`Triad_triad_bccFe_jij`,
  `golden_rtol` 1e-2).

### Lanczos etot shift at `2cae269`, three screening-alpha tables — 2026-10-05
- **Verified:** the legacy Lanczos deck etot moves +6.693e-6 Ry at `2cae269` ("Repair DRESP-03TG screening representation consistency"); Block and
  Chebyshev move 3e-7 and 4e-7; strux_lib is bitwise unchanged. Located in the `qm` read of `predls` (`source/symbolic_atom.f90`): `micha` constants
  are unchanged, but `potential%screening_alpha` is now filled on the legacy path, so `qm(1)` goes from 0.348485 to 0.3485. Why Lanczos moves ten
  times more than Block/Chebyshev is not diagnosed. `Regression_bccFe_lanczos` pins the current value (reference regenerated, no provenance recorded);
  cross-platform agreement at 1e-6 is not known.
- **Open decision, three alpha tables:** `legacy_micha_alpha = [0.3485, 0.05303, 0.010714, 0.00337]` (`source/lattice.f90:52`),
  `default_screening_alpha_values = [0.3485, 0.0530, 0.0107, 0.00674]` (`lattice.f90:48`, used by strux_lib; 0.00674 = 2 x 0.00337 has no source in
  the repository history before `1355a50`; the vendored `strux_tb.f90` has 0.00535), and `qm_canonical = [.348485, .053030, .010714]`
  (`source/math.f90:110`, used when `screening_alpha` is unallocated). `qm_lmto47` (`math.f90:112`) has no reader. Which table each backend should
  use, and the l=3 value, is not decided.
- **Suspected, from code reading only:** `ccor_2c` also reads `potential%screening_alpha` (`source/hamiltonian_ccor.f90:496-499`), so a legacy-backend
  `ccor_2c` deck changed at `2cae269`.
- **Related, not diagnosed:** the k-space SCF reference `Example_k_space_scf_bccFe/ref.json` is 7.047e-6 Ry above the current etot; the step localizes
  to the commit "Reconcile legacy LDA XC kernels against fixed-density references" (mechanism inferred from the title only).

### Four validation scripts fail in the CI full run — 2026-10-04
- **Verified (last full run, 4 of 179 failing):** `Val04LdaUPhysics` (Hermiticity residual of the stored occupation matrix 1.08e-7 against 1e-7 for
  the u4 deck; later also stops at the Fixed-U guard, see below), `Val12LmtoFieldsTorques` (SOC-free global rotation changes the torque invariant by
  8.281e-03 T against 5.0e-3), `Val16GbtCommensurateSupercells` (the `nstep = 100` reference supercell run does not reach `Converged!`; script
  retired, archive tag) and `Val17GbtHarmonicGoldstone` (`tests/validation/val_gbt_harmonic.py`: canonical occupation N = 8.0069 against `|N-8| <
  1e-8`, and sin^2(theta) spread 10.39% against 5%). A later cone-angle remeasurement gives a 1.49% spread, so the Val17 status is not reconciled.
- **Not diagnosed.** Scripts in `tests/validation/` (`val04_lda_u.py`, `val12_lmto_fields_torques.py`).

## Fixed-potential and frozen-magnon (GBT) findings

### Small-q GBT stiffness is mesh sensitive; multi-sublattice k-space not re-tested — 2026-08-16
- **Verified:** the fine-mesh small-q frozen-magnon energies are quadratic in q at 12^3 and 16^3, but the inferred stiffness changes by about 55%
  between those meshes. The material stiffness and the multi-sublattice maturity claim stay open until denser or shifted meshes and an independent
  finite-q reference exist.
- **Verified:** the reciprocal FeCo acoustic Gamma value is 6.48e-15 Ry (real-space route -4.16e-17 Ry, in-phase eigenvectors; two-sublattice bcc
  FeCo, `branch_mode = 'auto'`). The original ~0.28 Ry violation was measured on the pre-WP1 k-space route and was not re-checked for every system;
  one system, one theta_probe.
- **Test impact:** `Example_frozen_magnon_bccFe_auto` and `_auto_scf` are single-sublattice smoke cases only; they pin no multi-sublattice values.
  `GBT+SOC`, GBT local-cluster/impurity and GBT intersite Hubbard-V are unsupported combinations (`docs/DECISIONS.md`).

### k-space EBAND mesh refinement is mildly non-monotonic at some q — 2026-08
- **Verified (bcc Fe, `q_ss = (0,0,0.5)`, nk = 8/12/16):** `eband` = -1.98668388, -1.98714329, -1.98785909 Ry; the step grows (4.59e-4 then 7.16e-4
  Ry). At H the same sweep shrinks. Not diagnosed; plausibly tetrahedron discretization noise.

### GBT: suspected double sign flip for sites with `mom(3) < 0` (F1) — 2026-10-17
- **Read in the source:** `gbt_endpoint_angles` (`source/hamiltonian_build.f90`) sets `axis_sign = -1` and returns `theta = pi - theta_ss`, `phi = pi`
  when `potential%mom(3) < 0`; the two flips may cancel or be one intended flip, depending on the stored `wx1`/`cx1` convention. Not traced.
- **Measured (AFM bcc Cr, `nsp = 3`, q = 0, one iteration):** with `mom(3) < 0` on a site, `gbt_single_q` bands are 4.2e-2 Ry from `periodic_nc` and
  the moment is about 4.3 times smaller (0.207 against 0.884 muB); the same antiferromagnet with swapped channels and `mom(3) = +1` agrees to 2.7e-15
  Ry. Real-space GBT is frozen for the restart.
- **Measured again (2026-11-10, real space, two-atom Fe cell, `nsp = 3`, block recursion, q = 0, U = 0, no `hoh`, one step, hand-built deck):** with `mom = +z, -z` on the two sites `gbt_single_q` gives E_F -0.069256 Ry and |m| 2.104043 muB against `periodic_nc` -0.043287 Ry and 2.179464 muB; with U = 2 eV and `hoh` the antiferromagnet under `gbt_single_q` reproduces the ferromagnet to all printed digits. Not traced.

### `cmplx` without `kind` rounds double-precision arguments to single in four places — 2026-10-18
- **Read, not run:** `source/bands.f90:1366` and `:1444` (`cmplx(0.0_rp, res)`, Gauss-Legendre nodes), `source/green_block.f90:214`,
  `source/hamiltonian_paoflow_io.f90:236` (PAOFLOW hoppings). Expected size about 1e-8 relative in the real-space Green-function quadrature. Not
  measured; fixing changes reference values.

### Real-space LKAG exchange does not follow the full-ASA reference with Hamiltonian order — 2026-10-22
- **Measured (bcc Fe, `fe_bcc` potential, zero SCF steps), J(0)-J(q) at q1 / q2 [mRy]:** rigid rotation oracle 3.82 / 11.01; recursion `hoh = .false.`
  2.26 / 6.82; recursion `hoh = .true.` 4.54 / 15.64; k-space GBT first order 3.35 / 9.64, second order 3.60 / 10.20. J_01 from recursion 0.755 (`hoh`
  off) and 2.175 mRy (`hoh` on). Not matched between routes: temperature (0 against 300 K), energy window, terminator, 5a cutoff (about 1%).

### `tests/spin_response/oracles/fe_lswt.dat` is not reproducible from the repository — 2026-10-22
- **Verified:** `make_fe_lswt.py` only converts UppASD output; the Hamiltonian order, temperature, window and `nsp` of the run behind the file are not
  recorded.

### k-space GBT constrained spiral: controller diverges, S9 relation does not hold — 2026-10-20
- **Measured (bcc Fe, 16^3, 300 K, first order, q = (0,0.25,0), theta_ss = 30 deg):** with the `constrain.f90` controller (`i_cons = 3`, lambda = 10),
  m_perp goes 3.0e-4 -> 0.035 -> 2.8 -> 8.0 muB in three iterations and the eigensolver fails; the other sign overshoots equally.
- **Verified:** dE/dtheta = -M B^x does not hold in this scheme (SCF dE/dtheta 2.318e-3 Ry/rad at 30 deg against -M B^x = -7.026e-4; the force theorem
  gives 2.279e-3). The transverse moment in `write_gbt_scf_diagnostics` (`source/self.f90`) is a square root of a difference of squares and carries
  cancellation error near convergence (4.03e-8 against 4.17e-8 muB).

### Explicit-supercell SCF diverges at default mixing; Cr flat spirals collapse — 2026-10-20
- **Measured:** collinear `nsp = 1` SCF of the 8-atom bcc Fe supercell (`tests/gbt/fixtures`, `beta = 0.1`, Broyden) drops to Diff 1e-5 after about 10
  iterations, then grows. For `tests/spin_response/cr_bcc` (2.95 A, AFM) the flat spiral at q = (0,0,1/4) on 8x8x16 converges to |m| = 0 on both sites
  (GBT and supercell); on 4x4x8 even the collinear AFM state collapses. Not diagnosed.

### `ccor_vmt_mode = 'pair_surface'` is not invariant under the spin-swapped encoding — 2026-10-22
- **Measured (AFM-II NiO, `periodic_nc`, `nsp = 3`, 4x4x4, `ccor_2c`, fixed potentials):** max |E(mirrored) - E(B1)| is 1.42e-3 Ry (first order) and
  1.43e-3 Ry (second order) for `pair_surface`, against 5e-15 and 1e-14 Ry for `surface_scalar`. The k-space GBT route with the B1 file matches
  `periodic_nc` to 1e-14 Ry, so the difference is in the ordinary route.
- **Suspected, from reading:** `ccor_vmt_pair_surface` (`source/hamiltonian_ccor.f90:598-603`) averages `vrmax` index by index before the moment
  orientation is applied (`:661-673`); the same average is in `build_ccor_pair_block_gbt` (`:210-215`), not run. The k-space GBT glue refuses
  `pair_surface`. Which encoding is intended is a developer decision. Test retired; archive tag.

### I3 with `ccor_2c` needs a tighter SCF criterion — 2026-10-22
- **Measured:** with the default `conv_thr` (5e-9) the two collinear SCFs of the CsCl cell and the 8-atom supercell differ by |dE/atom| = 1.39e-7 Ry
  (limit 1e-7) at second order although both report `Converged!`; `conv_thr = 5e-10` gives 3.3e-9. First order needed `conv_thr = 1e-9` and
  `OMP_NUM_THREADS=1` (1.84e-7 Ry before). The effect of `conv_thr` and of the thread count was not separated; a stall at Diff about 6e-10 in a Debug
  `-O0 -fcheck=all` build was not reproduced. Test: `test_kspace_scspiral.py` family; the I3 script is retired (archive tag).

### k-space `'orthogonal'` SCF of AFM NiO without CCOR collapses, then LAPACK fails — 2026-10-22
- **Measured (NiO AFM fixture, `nsp = 1`, 6x6x6, 300 K, no `hoh`, no `ccor_2c`):** with `beta = 0.1` the SCF difference runs 0.23, 0.50, 1.2, 3.4 and
  `LAPACK eigensolver failed` follows in iteration 3; with `beta = 0.02` Broyden the Ni moment falls 1.53 -> 0.004 muB in eight iterations, then the
  same failure. Untriaged: physics of the full-ASA Hamiltonian without CCOR, the start, or the mixing. Fe converges. An NiO `'orthogonal'` SCF waits
  for CCOR with `'orthogonal'`.

### `nsp = 1` AFM with `mom = -z` in a potential file collapses the moment — 2026-10-21
- **Verified (12^3, `tests/spin_response/cr_bcc` as an `nsp = 1` deck):** `Cr1.nml` with `mom = -z` gives site moments 0, 0 (E = -4195.4943 Ry)
  against 0.919776 for the shipped `+z`; the same collapse at first order and `'orthogonal'`.
- **Suspected, not checked:** other code reads `mom(3)` in the `nsp = 1` path (`source/hamiltonian_build.f90`, obarm/enim lines;
  `source/calculation.f90`), swapping the channels of a mirrored file. Convention: `nsp = 1` AFM is encoded by mirrored potentials with `mom = +z`.

## Recursion, Green-function and GPU findings

### Fast block Lanczos with `hoh` appears to count `lsham` twice — 2026-10-20
- **Read, not run:** `block_lanczos_sp` folds `lsham` into the shell-1 block of `he`/`ha`, so the eeo sweep acts on it and the hoh combine adds `lsham
  psi` again (`ons = enim + lsham`); legacy `hop_b_hoh` applies `lsham` once and outside the sweep. The operators would differ by `lsham psi - eeo
  (lsham psi)`.
- **Measured (block exchange, `tests/postproc/cases/exchange/bccFe`, `hoh`, `lld = 20`):** J of the first two shells at `nsp = 2` is 1.919232,
  0.371056 (legacy) against 1.908953, 0.369334 (fast, fp32 and `fast_dp` agree, so not precision); at `nsp = 1` legacy and `fast_dp` agree (1.933942).
  Cause untested: `nsp = 2` differs from `nsp = 1` in more than `lsham`.
- **Proposal, not started:** one fast-versus-legacy case with SOC on and `hoh = .true.`; if it differs, remove the shell-1 fold from the hoh branch of
  the two kernels (about 10 lines).

### RS block recursion at `nsp = 2` is not invariant under spin-axis rotation — 2026-11-01
- **Verified (bcc Fe, U=0, `soc_scale=0` with `fix_soc`, `mom` along z and along x, iteration 1):** block recursion (lld 21) differs, dE_F 3e-4 Ry,
  dmoment 1.8e-4 muB; Chebyshev about 1e-7; the KS solver agrees within 1e-9. `sym_term = .true.` does not resolve it (converged E_F -0.068999 z,
  -0.068975 x; moment 2.099999 and 2.10006). With U=2/J=0 local mode the converged moment differs by 1.1e-2 muB. Cause open; a termination in the
  global frame is a suspect.
- **Handling:** `nsp = 2` U runs with a tilted spin axis (RS, `recur` other than `chebyshev`) stop with a fatal. Related: global-frame recursion
  (`local_axis = F`) is not rotation invariant even without SOC (band energy +2.3e-4 Ry at lld 13, -1.9e-5 at lld 25); `local_axis = T` reproduces z
  to 6e-6 Ry. Block recursion at `nsp = 4` breaks down at large lld (|m| 2.104, 2.087, 1.776 muB at lld 21, 31, 41), consistent with the lld <= 28
  bound in `docs/DECISIONS.md`; cause not diagnosed.

### CPU Chebyshev DOS reconstruction is single precision for every `cheb_backend` — 2026-10-21
- **Read:** `chebyshev_green_core` (`source/green_chebyshev.f90`) reconstructs the DOS on the CPU with `cheb_green_fast` (single precision) in every
  branch without `eta`; `cheb_backend` is not consulted. So `legacy` is not an fp64 oracle for the SCF energy.
- **Measured (one SCF step):** |etot(GPU fp64) - etot(CPU legacy)| is 9.5e-10 (diamond Si), 5.8e-9 (`bccFe_chebyshev_fast_hoh`) but 2.6e-8
  (`bccFe_chebyshev_legacy_hoh`) and 2.3e-8 (`_fast_hoh_ccor_2c`). That the 2e-8 to 5e-8 comes from the CPU reconstruction was not verified. The fp64
  rows of `tests/regression/run_gpu_matrix.py` therefore use two loose criteria and cannot detect an fp64 kernel error below about 1e-7.
- **Update:** `cheb_backend` is consulted since the default became `fast_dp`: `fast_dp` and `legacy` reconstruct with `cheb_green_dp` (double precision, 4e-16 against a literal sum in `UnitChebyshevGreenPrecision`); `fast` and the other fast backends keep `cheb_green_fast` (single precision, 2e-7). The 2e-8 to 5e-8 attribution to the reconstruction above was not re-measured.

### `SpinResponseOmpIdentity` fails in the CUDA/MKL build — 2026-10-21
- **Measured:** `tests/spin_response/omp_identity.py` on `fe_bcc` reports differing `spin_response_q*.dat` files between 1 and 4 threads in the
  `build_gpu` binary. Not diagnosed; suspects: MKL threading of the CUDA build, and a library mismatch in that setup (binary linked to MKL 2026.2,
  `LD_LIBRARY_PATH` at 2026.1).

## Layered (`calctype = 'L'`) and surface findings

### Raising `&charge nlay_a/nlay_b` breaks the alignment fixed point — standing trap
- **Verified (`example/interface/fccCu111_AA`, 5 iterations):** `&lattice nlay_a/nlay_b` 4/4 keeps `V(B) = 0.000000` Ry; `&charge nlay_a/nlay_b` 4/4
  gives `V(B) = -0.4498` Ry. `&charge` counts rows of the synthetic 2D Madelung stack (fixed size `max(49, ...)`, `source/lattice_cluster.f90:709`),
  so raising it relabels charged interior rows as frozen (`source/charge.f90:1711`). Widen the buffer with `&lattice` and keep `&charge` at 1.
- **Open:** the `align_regions` "widen the active zone" drift warning also fires in runs where the identity holds; whether the deep-probe drift
  threshold needs recalibration is open.

### `A | vacuum` and `A | vacuum-gap | B` — 2026-08-16
- **Open:** `A | vacuum-gap | B` needs a four-region layout and a rework of the type-block arithmetic in `build_interface_full`, which hardcodes three
  regions; unavailable.
- **Residual, `A | vacuum` at `nlay = 3`:** the undamped alignment/vacuum feedback can still oscillate or grow from a bad initial offset; `vmix` is
  exposed in `&charge` (0.2 used in the multilayer test); the vacuum-onset warning is intentional. Example: `example/interface/fccCu111_Avac`.

### `processing = 'sd'` workflow scope — 2026-08
- **Limitation:** `Example_bulk_bccFe_sd_smoke` is an execution smoke test only. Broader impurity and multi-site spin dynamics are untested; MPI
  launch of the B2FeCo deck was not reproduced.

### `write_strux_block` / `buildsurf` with `strux_lib`: multi-species surfaces unchecked — 2026-10-25
- **Not checked:** the local branch of `structb_strux` (`source/lattice_strux.f90`) reads the species of each cluster atom from `this%izp(1:kk)`,
  which `build_surf_full` does not fill per cluster atom (`this%iz` holds the types). Single-species surface decks run; multi-element surfaces and
  interfaces were not run.

### `constraints_bfield` is pre-sized like the former `constraints_mom_ref` — 2026-11-01
- **Not tested:** `source/control.f90` pre-sized `constraints_mom_ref` (1,1), which made indexed input fail; that array was fixed.
  `constraints_bfield` has the same pattern. A zero field is valid there, so the trailing-column cut used for `constraints_mom_ref` does not apply.

## spin_response

### kB defined ten times — 2026-10-05
- **Verified:** `kB_Ry_per_K = 6.3336814e-6_rp` is declared locally nine times (`source/reciprocal.f90`, `reciprocal_dos.f90`,
  `reciprocal_occupations.f90`, `reciprocal_spin_density.f90`, `self_reciprocal.f90`), and `math.f90` holds `0.633362019d-5` (`kB_simpson_f`), a
  relative 9.7e-6 different; not unified because `simpson_f` is on the bit-identical exchange path. Proposal: one definition (about 10 lines; a 1e-5
  relative change in kT at finite T).

### `<prefix>_dispersion.dat` rows run `n/a` into the previous number — 2026-10-07
- **Verified:** with no peak or crossing, rows read `2.19599654E-01n/a` (format `(7es16.8,2a24)` with left-justified texts); seen for every row with
  `n_omega = 1`. Proposal: a leading `1x` before the two `a24` items (one line in `run`); not applied.

### Radial data at r = 0 — 2026-10-11
- **Verified (bcc Fe):** `radial_ground_state` holds `vxc_up(1) = vxc_down(1) = 0` at r = 0 while `n_up(1) - n_down(1) = -0.8069` bohr^-3
  (`source/radial_ground_state.f90:257`), so dV/m is 0 there, not its limit; `radial_kernel_matrix` starts every sum at the second point.

### Experimental ALSDA kernel: gaps in the checks — 2026-10-11
- **Verified (bcc Fe, 12^3, `alsda`, Juelich, lambda 1.195):** the projector arithmetic of `site_radial_data` has no independent check: x from the
  opposite spin moves lambda to 1.362 (+14%), flipping the sign of x to 1.675 (+40%), core left out of m to 1.1897; all pass. A transposed or
  conjugated vertex in `accumulate_chi0` passes A4.
- **Measured (24^3 and 60^3, scripts not committed):** the q = 0 mode stays gapped with spd (smallest eigenvalue of 1 + chi0 U -3.3e-3 Juelich,
  -2.1e-2 Mills at 24^3); Fe at 60^3 has delta_q < 0 for `alsda_scaled` spd at several xi (-0.022 to -0.112 for Mills). Per-channel scaling not
  pursued.

### `source/spin_response.f90` is past the 1,000-line guideline — 2026-10-12
- **Verified:** 1187 lines now (1146 when recorded). The developer decided not to split now. Proposal: move the pure kernels (`accumulate_chi0`,
  `u_juelich`, `solve_dyson`, `adiabatic_frequencies`) to a submodule (about 600 lines).

### bcc Cr antiferromagnet: response frequencies 1.6 to 2.4 times LSWT — 2026-10-13
- **Measured (32^3, `tests/validation/spin_response_cr`, 94 q):** the adiabatic frequency over the LSWT `omega_plus` of `oracles/cr_lswt.dat` is 1.626
  to 2.374 (median 1.727; 64^3 the same); the force-theorem route is within 0.968 to 1.087 of LSWT for q >= 1/8, and Omega_RPA / Omega_force equals
  1/sqrt(1 - lambda_u) = 1.9657 on Gamma-X. So the factor is the ratio of two routes built from one chi0. Cause of the route ratio not diagnosed; the
  dynamic tr L peak (1.6 to 1.8 times Omega_RPA) and its maximum at the window edge for 68 of 94 q (eta = 2e-3 Ry, window +-0.06 Ry) were not
  examined; the Fe-style pole weight is not defined here.
- **Verified:** `spin_swapped` follows the sign of the total d moment, +5e-14 to +2.5e-13 muB here; a flip on another compiler would change
  `_state.dat` and fail `SpinResponseCrBaselineRepro`.
- **Deferred:** site-resolved spectral output (L per row) is not written. Proposal: a dynamic run with a wider window and several `eta` (about 20
  lines).

### SOC gap in cubic bcc Fe carries an O(xi^2) bias (NC-2, Gate 1) — 2026-10-14
- **Measured (`tests/validation/spin_response_fe_soc`, `nsp = 2`, q = 0, omega = 0):** omega_adiabatic -8.0e-6 (12^3) to -2.3e-5 Ry (48^3), scaling as
  xi^2 (log-slope 1.995, 1.998); symmetry expects fourth order (K_1 about 1e-10 Ry) and a gap below 1e-8 Ry. SOC runs print a warning.
- **Suspects:** (1) the Juelich factor multiplies up and down components with different (enu, p) (a one-site toy with per-spin vertex scaling
  reproduces the sign and xi^2 slope); (2) the `K = U m0/m` prescription; (3) on-site L.S of the first-order k-space H; (4) deck handling. Only the xi
  = 0 identity (`SpinResponseSocZeroIdentity`) checks the setup.

### Symmetry log cannot check physical inversion — 2026-10-14
- **Measured:** two Fe atoms with different type labels, second atom at (1/2 + 0.05, 1/2, 1/2): the log prints P4mm (#99) (no inversion) while spglib
  with element labels gives P4/nmm (#129) with inversion; omega(q) = omega(-q) to 1e-16 Ry. Likely cause (read): `source/symmetry.f90` passes
  `lattice%iz(i)` as the species. Check inversion with spglib and element labels.

## Fixed-U LDA+U (RS and KS)

Evidence and run history: `docs/validation/fixed_ldau/VALIDATION_SUMMARY.md`; older records in
the archive tags.

### k-space `second` SCF of AFM NiO with `hoh` and `ccor_2c` breaks the Ni1/Ni2 symmetry at U=0 — 2026-11-07
- **Verified (user-supplied deck `NiO_AF_v2`: 4-atom AFM-II cell, `nsp=1`, `hoh=.true.`, `ccor_2c=.true.`, `use_kspace=.true.`, `kspace_ham_order='second'`, nk=8, beta 0.1, 20 steps, pristine potentials, no U):**
  Ni1/Ni2 spin moments (muB) 0.774/0.774 (step 1), 0.753/0.713, 0.593/0.001 (step 3), 0.568/0.034 (step 10), 0.509/0.021 (step 20); Diff 3.0e-2 -> 4.9e-4 (looks like convergence). beta 0.05: 0.508/0.019 at step 20. beta 0.2: erratic, E_F between -0.158 and -0.010 Ry.
  The same deck in real space (Chebyshev, lld 200) keeps Ni1/Ni2 = 0.7117/0.7109 -> 0.7142/0.7143 over 4 steps (0.74 by step 9).
- **Not verified:** the cause (candidates not tested: k-mesh and symmetry handling for the spin-polarised cell, E_F in the gap at 300 K). The deck is hand-built and not reproduced with a second NiO deck.
- **Possibly the same cause as a fixed k-space `nsp=1` AF bug** (a mirrored-channel site had `potential%mom` flipped to -z each iteration and its density projected on the measured axis; fixed for the two-atom Fe cell and MnO, `ks_af_mirror_nsp1`): this deck was not re-run.
- **Related:** real-space +U on the same deck: Ni moment falls with the SCF steps at U=6 and 8 (0.39 to 0.52 and 0.14 to 0.19 muB after 6 to 8 steps) while E_F drifts to -0.4 Ry; not diagnosed. k-space with U not run.

### AFM bcc Fe cell with U=2 eV: block recursion at `nsp=3` depends on the quantization axes (mostly the terminator) — 2026-11-08
- **Verified (`af_deck` in `tests/ldau/test_nsp34_u.py`, bcc Fe 2-atom AF cell, U=2 eV, J=0, spectral radial moments):**
  Block recursion, `nsp=3`, 2 steps, moments (z, x) against the same state rotated to (x, -z): E_F identical (+0.01929, -0.04356 Ry), site |m| step 1 2.33066/2.32922 against 2.32922/2.33066 and step 2 2.50357/2.50950 against 2.50950/2.50357 muB: the two sites swap their values, i.e. the result is not invariant under a global spin rotation. With `sym_term = .true.` (orbital-independent terminator, &control) the same pair gives E_F +0.018950, -0.043334 Ry on both and site |m| step 1 2.333923/2.333990 against 2.333990/2.333923, step 2 2.508328/2.508447 against 2.508447/2.508328 muB: the swap shrinks from 1.4e-3 / 5.9e-3 to 6.7e-5 / 1.2e-4 muB but is still above the log precision. The default (`sym_term = .false.`, orbital-dependent terminator) is what the numbers above were run with. The same rotated pair with Chebyshev agrees to the log precision (the `af_angle` check).
  Block recursion for the collinear AF (m2 = -z) at `nsp=3` and at `nsp=1` (swapped files) agree with each other (|m| 2.29663, 2.59689 on both sites).
- **Not verified:** causes; the `sym_term` comparison suggests the default orbital-dependent terminator is the main part, the 1e-4 residual is not explained. Nothing else was triaged: the block-recursion spinor path for non-collinear starting axes and the k-space `nsp=1` AF assembly are untouched suspects. The k-space `nsp=1` behaviour looks like the NiO entry above (same symmetry breaking, there at U=0). A hand-built deck is a live suspect for the block case until a second deck reproduces it.
- **Coverage:** `LdaUNsp34Noncollinear` uses Chebyshev for the rotation checks and k-space `nsp=3` for the mirror check; block recursion at `nsp>=3` with non-collinear axes is not tested.
- **Proposal (not started, owner decision):** one fatal or warning for block recursion with `nsp>=3` and U; or a triage of the block path (size not estimated).

### `Example_bulk_bccFe_nsp2_hubbard` stops at the Fixed-U guard — 2026-10-18
- **Verified:** the run stops with "Fixed-U LDA+U currently requires collinear scalar-spin input" (`source/hamiltonian_build.f90:768`);
  `Val04LdaUPhysics` stops at the same guard (its script rewrites `nsp = 1` to `nsp = 2`). In a scratch copy with `nsp = 1` the guard passes and the
  run aborts in `update_fixed_u_residuals` (`source/self.f90:1843`), "Index '4' of dimension 2 of array `this%mix%_data%ldm_new` above upper bound 3".
  Not checked: Release behaviour, Val04 at `nsp = 1`.
- Both tests were disabled (`DISABLED TRUE`, 2026-10-26). Either they move to `nsp = 1` after the bound is triaged, or the guard is extended once a
  spinor-basis oracle exists; developer decision.

### Co bulk U=0 non-finite — OPEN
- **Verified:** RS records of the fixed Co bulk deck become non-finite after about 19-40 iterations at U=0; reproduces on `8b8f7ef`. A fixed Co block
  deck run to its `nstep = 25` cap stayed finite and ended at Diff = 0.0054412342 Ry with 19 valid checkpoints (not accepted). Cause not diagnosed.

### Fe U>0 embedding offset — OPEN, low priority
- **Verified:** an Fe impurity in an Fe host at U=3.8 eV, J=1.1 eV differs from RS bulk by -6.6e-3 muB, about 735 times the U=0 control, with
  identical Hubbard operators, Fermi level and Green's-function evaluation. Ni does not show it.
- **Probable cause (not tested):** cluster mismatch, rc=80 for bulk against rc=100 for the impurity; the matched test is blocked by the rc=100 entry
  below.

- **Update (impurity-in-host test decks, 2026-11): the dominant offset is the fixed Fermi energy of the impurity deck, not the Hamiltonian.** `calctype I`
  keeps E_F at its input value; the bulk step determines E_F from the electron count and it moves with hoh and U. bcc Fe, nsp = 1, one step from the same
  state, rc = 40 (impurity 'Fe-imp' copy of the host, identical potentials, same cluster size of 748 atoms in the bulk and impurity `clust` files).
  Verified (bulk against impurity spin moment, muB): U = 0 first order 2.0837 against 2.2215 with the deck's E_F = -0.042265 Ry, 2.0797 with the bulk step's
  E_F (-0.069773); U = 0 hoh 2.273531 against 2.273516 (E_F -0.084057, max parameter difference 8e-6); U = 2 hoh local 2.521310 against 2.521297 (E_F
  -0.110565, 9e-6); U = 2 first order local 2.349640 against 2.349847. The earlier offsets (hoh 2.52 against 2.35) were the stale E_F: the bulk E_F moves
  by 0.04 to 0.07 Ry with hoh, U. Registered as `LdaUImpurityEqualsBulk` (rc = 20, hoh, local, U = 2: 2.471169 against 2.471392; with the stale E_F the
  impurity gives 2.3775, so the test can fail). So build_locham and build_bulkham agree for the hoh local route (the hoh blocks, enim/obarm at the cluster
  edge, V_U insertion): no builder bug found.
- **Residuals, open (not diagnosed):** first order (non-hoh) at U = 0: 4e-3 muB; first order U = 2: 2e-4; hoh spectral U = 2: 3.5e-3 (2.514966 against
  2.518487, parameters 7e-4). The impurity run reports a nonzero charge transfer at the central atom (-2.9e-3 e first order U = 0, -9.8e-3 e hoh
  spectral; zero to 4e-6 where the moments agree to 1e-5), i.e. the cluster electron count at the bulk E_F is not the bulk count. Moving E_F by +1e-4 and
  +5e-4 Ry changes the first-order impurity moment by -5.4e-3 and +2.4e-3 muB (non-monotonic, 2.0797, 2.0743, 2.0821), so these residuals are within the
  sensitivity to E_F (the bulk E_F is printed with six decimals and the Green's-function energy integration is noisy); SPECULATIVE that this is all there is.
  Excluded by runs: cluster size (same atom count), different E_F precision is not excluded. Not checked: recursion start vectors/terminator differences, Madelung.

### Spectral radial moments under PBE: NaN potential, stop in the k-space eigensolver — OPEN (partly diagnosed)
- **Verified (instrumented run, NiO AFM-II 4-atom cell, k-space nsp = 1, hoh + ccor_2c, txc = 108, spectral, U = 2, J = 0, restart from the converged PBE U = 0
  state, `ldm_beta` = beta = 0.1, one thread; stops at step 33 in every run):** the cause is **not** `calculate_pl`/POTPAR. At the failing step PL(s,p,d) of Niu
  are 4.607, 4.317, 3.786 (the previous step: 4.607, 4.317, 3.787, all POTPAR solves finite, p eigenvalue -0.55 Ry); the atomic self-consistency loop of that step
  (`atomsc`, fixed PL and Q) does not converge: its radial-density residual DRHO starts at 3.4e3 in every step (as in the earlier ones, 57-59 iterations to the
  tolerance), but at steps 24, 26, 27 and 32 it hits the 80-iteration cap and at step 33 it grows to 1e2 and then, at iteration 60, becomes NaN. During these
  iterations the spin-2 radial density goes negative (min -0.9 to -0.1 e/bohr), the spin-2 valence eigenvalues oscillate (s: 0.31, -0.54, 1.25, 2.1 Ry), the
  radial eigenvalue search finds no solution (`RSEQ NIT GT 400 AND BAD NODES`, E = 10 to 50), the density and then V become NaN (495 of 495 radial points), and
  POTPAR returns NaN C, Q, srdel and enu(d) with the sentinel 10 for enu(s,p). The NaN reaches the k-space Hamiltonian and the run stops in
  `reciprocal_backend.f90:810`. With the Fermi energy of the base deck (-0.093245 instead of -0.263911) the run did not crash within 80 steps (Diff 4e-3), and with
  `ldm_beta` 0.3 / radial beta 0.03 it still crashed in earlier runs, so the trigger is the state of the atomic loop and not the outer mixing.
- **Changed:** `reciprocal_fourier.f90` (`build_kspace_hamiltonian`) stops at the first non-finite C, Q or srdel with `B1 non-finite guard: atom=Niu site=1 channel=s l=0
  spin=1 ... field=potential parameter C (atomic radial solve failed)`; the run stops at the same step as before, only with a clear message (test
  `LdaUNonfinitePotentialGuard`, fails without the guard). The Diff of the LDA NiO local U = 2 run (4 steps) is identical to the stored log to ten digits.
- **Not fixed (needs `self.f90`, off-limits):** the atomic loop in `atomsc` (about line 3915; `rho_in = BETA1*RHO + (1-BETA1)*rho_in` with `BETA1` 0.3/0.5 from
  `atomsc_mixing_beta`, no guard on the sign of `rho_in`, GGA gradient taken from the negative density in `VXC0SP`) is the place to act. SPECULATIVE proposals,
  untested: clamp `rho_in` to non-negative before `VXC0SP` (PBE only), lower the atomic BETA when DRHO grows between iterations, and stop with a message naming
  the atom and iteration when `rho` is not finite. It is not known why PBE does so and LDA (txc = 106) does not at this U; not investigated whether the negative
  density comes from the spectral correction terms (Q2, hoh) or from the GGA potential.

### RS rc=100 runs non-finite or non-converging — OPEN
- **Verified:** the Fe rc=100 impurity-style deck becomes non-finite after 102 iterations at U>0 and 34 at U=0; the rc=100 U=0 bulk is not accepted
  after 600 steps; reproduces on `8b8f7ef`. In a Debug trajectory the d-spin occupation exceeded its five-state capacity at iteration 90 (`ql =
  16.0232`). Collinear `nsp = 1` SCF now stops when q0 exceeds `2l+1` per spin by more than 1e-6 and rejects non-finite mixed values. A Chebyshev
  rc=100 run reached Diff 1.40e-8 Ry. Fe rc=140 decks are invalid setups and not part of this entry.
- **Note:** RS/KS results are not bitwise reproducible across thread counts (about 1e-15 relative); use one thread for identity checks.

### Ni impurity W1 lambda-zero SIGSEGV in optimized builds — 2026-10-21
- **Verified:** the frozen Ni lambda=0 W1 deck exits SIGSEGV in the first recursion in Release and RelWithDebInfo; the Debug bounds run passed 36
  loops. One cause was found: an unassociated `pl_bands%control` read in fixed-U replay (`source/self.f90` replay setup, `source/bands.f90:2254`).
  After binding the pointer, one complete iteration ran in Release without a crash.
- **Open:** longer-run behaviour and impurity energy (RD4) are untested; impurity energy routing stays off.

### KS-local host to RS impurity equivalence for Ni — 2026-10-22
- **Not established:** `auto` selects local impurity moments only when every U-carrying frozen host records local mode, for RS-local and KS-local
  hosts. No accepted Ni KS host exists, so the Ni KS-host-to-RS-impurity comparison is untested. Impurity Hubbard energy is excluded.

### Surface LDA+U (experimental): Ni baseline and local-density physicality — 2026-10-22
- **Verified:** Ni(001) U=3/J=0.7 probes with local mode stop on NaN projected DOS / `g0` or on a negative rotating-frame density eigenvalue (site 2,
  l=0, `lambda_min = -3.16e-6`); widening the window from 0.5 to 1.5 Ry does not change it. Two-step Diffs on a two-layer fixture are 0.0163 and 0.387
  (U=3) against 0.0162 and 0.384 (U=0), so the large Diff is not isolated to LDA+U; it arises in the surface charge response from the bulk-potential
  seed and is not fixed by lowering beta.
- **Update 2026-11-04:** collinear Ni(001) U=3/J=0.7 (types 1, 3, 4, U=3 bulk host) ran 200 steps to Diff 3.5e-6 (first below 5e-6 at step 179;
  spreads 1e-6 to 3e-5); no accepted no-U surface reference exists. Surface Hubbard-energy routing stays off.
- **nsp=2 surface (S6c), short runs only:** the SOC->0 surface equals the `nsp = 1` surface at iteration 1; with SOC on, the orbital moment keeps
  rising (Ni-S 0.066 -> 0.109 muB in 5 steps) and a 20-step Tier 2 run shows a Diff excursion near iteration 15-18 in three runs (unconverged host,
  converged host, J=0); a mixing-history cause is a hypothesis, not tested. No `nsp = 2` surface was run to its convergence floor.

### `calctype='B'` with an impurity cluster stops in the local Hubbard correction at U=0 — 2026-10-26
- **Verified:** the Ni impurity deck (`pre_processing='newclubulk'`, two types, Chebyshev, `hubbard_radial_moments='local'`, all U/J zero) with
  `calctype='B'` stops with a fatal in the local Hubbard correction. **Not verified:** the cause; the `B` path uses `itype_local = ib(na_glob)`, which
  can exceed the number of Hubbard types.
- **Proposal (not done):** a fatal at input time, or skip the correction when no atom carries U/J (about 10 lines).

### `soc_scale` rescales only the restart SOC strength — 2026-10-26
- **Verified:** `soc_scale` multiplies `xi_p` and `xi_d` once at `&self` read (`source/self.f90:1787`); after each atomic SFC they are recomputed from
  the radial potential unless `fix_soc = .true.`. So `soc_scale=0` needs `fix_soc`. Proposal: document it, or apply `soc_scale` after the
  recomputation (about 3 lines).

### RS local-mode +U transient amplifies 1e-13 seed perturbations — 2026-10-26
- **Verified (bcc Fe, RS block, `nsp = 1`, U=2 eV, J=0, local, beta 0.01 Broyden, 20 steps):** a 1e-13 relative change in one seed occupation changes
  the logged Diff substantially, so two implementations that differ at round-off cannot be required to agree to 1e-10 over the first 20 iterations.
  Converged `nsp = 1` +U RS runs stall at Diff about 2e-6 with beta 0.01.
- **Not diagnosed:** which step amplifies (suspect the discrete Fermi-level grid search).

### SOC+U bulk adoption: untested corners — 2026-10-27
- **Not tested:** RS J derivative at `nsp = 2`; k-space meshes other than 12^3 and materials other than bcc Fe; multi-site cells; `nsp = 2` with U on
  more than one shell. The RS U=2.05 eV Hellmann-Feynman run ended at Diff 1.57e-6 after 150 steps (the transient above).
- **Unchanged limits:** `calctype S` with U at `nsp = 2`, and `nsp = 3/4` with U, stop with a fatal.

### SOC+U impurity at `nsp = 2`: untriaged findings — 2026-10-30
- **Verified:** Ni (RS Chebyshev, local) U=3/J=0.7 impurity runs at beta 0.05 and `ldm_beta=1` oscillate (Diff 1e-4 to 7e-2 over 200 steps) at `nsp =
  1`, SOC->0 and SOC on; they need damping. The KS-host smoke (bcc Fe `nsp = 2`, U=2/J=0, `auto` resolves to local) ended at Diff 1.3e-7
  (`ldm_beta=1`) and 2.5e-6 (`ldm_beta=0.3`) after 100 steps; the 20-step spread criterion was not met. The impurity local-mode transient amplifies
  round-off between the `nsp = 2` (10x10) and `nsp = 1` (5x5) paths over 4 steps (7.15e-9 measured against a 1e-7 row tolerance).
- **Note:** the element `symbol` in the impurity seed names the output files; relabel a copied seed.

### Constraint controller diverges in a KS `nsp = 2` bcc Fe smoke (untriaged) — 2026-11-01
- **Verified:** KS bcc Fe (`k_space_scf/bccFe`, 12^3, `nsp = 2`, SOC on) with `&constraints`, `constraints_enable`, `constraints_mom_ref(:,1) =
  1,0,0`, seed `mom` along x: the target is read and resolved, then the controller diverges.
- **Not verified:** whether the hand-built deck is wrong; the only constraint deck in the repository
  (`example/bulk/bccFe/input_frozen_magnon_constrained.nml`) is an `nsp = 4` frozen-magnon run. So there is no converged constraint-target result for
  the S6b path.


### Radial moments with U: `spectral` moves the Mn d band against O p about twice as far as `U/2` and as Quantum ESPRESSO; `local` agrees with QE — 2026-11-08
- **Measured (MnO AFM-II, ideal rocksalt a = 4.446 A, PBE `txc = 108`, U = 2 eV, J = 0; shifts U=2 minus U=0 of the occupied d-up centroid against the occupied O p centroid):** HOH + ccor_2c `nsp = 3` spectral -2.00 eV (RS and k-space agree to 0.005 eV,
  start from the no-U state gives -1.99); first order `nsp = 1` spectral -1.72 to -1.75 eV (not fully settled, `Diff` ~2e-5 to 1.8e-4); first order `nsp = 1` `local` -1.11 eV; Quantum ESPRESSO PBE (GBRV Mn, PAW O, U_eff = 2, `ortho-atomic` / `atomic`) -1.24 / -1.29 eV.
  `V` alone gives about -0.9 eV. The Mn charge change is +0.05 e (spectral, HOH) and -0.014 e (`local`); QE -0.011 e.
- **Measured, isolating the part (first order, `nsp = 1`, 30 steps, temporary switch, not committed):** applying only the centre correction `c' = c - Tr(VP)/q0` reproduces `local` (-1.02 eV, bare Mn-O -0.20/-0.33 eV); applying only the `q2` correction stays near `spectral` (-1.58 eV, -0.89/-0.81 eV).
  First order against HOH + ccor_2c at `nsp = 3`: the extra shift of HOH is about 0.2 to 0.25 eV; `nsp = 1` and `nsp = 3` first order give identical `Diff` sequences to seven digits.
- **Measured with `local` under HOH (k-space `second`, HOH + ccor_2c, `nsp = 1`, 40 steps, `Diff` 3e-5 to 5e-5; the real-space twin of the same deck agrees on the Mn moment to 0.1% over its first two steps):** d-up -1.12 eV, d-down -0.25 eV, split -0.87 eV against the same deck's no-U state; Mn charge -0.019 e; Mn moment +0.127 muB
  (HOH `spectral`: -1.93, -0.21, -1.72 eV, +0.046 e, +0.104 muB; QE PBE: -1.24 / -1.29, -0.45 / -0.34, -0.80 / -0.95 eV, -0.011 e, +0.11 / +0.12 muB). The committed HOH `local` test (`af_fe_local_hoh`) compares `nsp = 3` with `nsp = 1`; there is no independent reference for the HOH correction itself except this comparison.
- **Read in the source:** `gravity_center` (the linearization energy, through `calculate_pl`) includes `V` in `spectral` mode and excludes it in `local` (`correct_d_shell_moments`, `math.f90`).
- **Not verified:** the mechanism beyond the above (a derivation from `psi = phi u + phidot (h u)` makes `local` the h-moment construction, algebra checked against `correct_d_shell_moments`; no literature source for the ASA radial-moment step was found); one material; QE at 6^3 k and one cutoff.

### Ni(001) surface at `nsp = 3` differs from `nsp = 1` from step 2 (untriaged, hand-built deck) — 2026-11-10
- **Verified:** `tests/ldau/decks/ni001_soc_u3` (calctype S, U = 3, J = 0.7, 4 steps, Chebyshev, window -3 to 3 Ry) with the `host_soc0` host relabelled `hubbard_nsp = 3` (the nsp = 2 spinor dm kept) and the seeds as shipped, moments along +z: step 1 agrees with `nsp = 1` (all atom moments identical to 6 digits), from step 2 the spin moments differ (ES 0.0051 against 0.0310 muB, Ni-S 0.7858 against 0.7999, Ni-S-1 agrees to 3e-4; max over 112 signature quantities 2.6e-2). The same with `hubbard_radial_moments = 'spectral'` (2.6e-2) and with hoh and `local` (1.1e-2): so it does not come from the local correction or hoh.
- **Not verified:** the cause. Suspects, in order: the relabelled host/seed files (the seeds carry zero spin coherences and no `spinor_d` seed for the surface layers; a spinor-P mixing start different from the scalar one), then the nsp = 3 surface route itself. No control at U = 0 or with a regenerated nsp = 3 host was run, so a surface nsp = 3 parity test for `local` is not available.
