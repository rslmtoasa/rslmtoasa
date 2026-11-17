# Fixed-U (LDA+U) campaign: measured results

Summary of the evidence directories `docs/validation/fixed_ldau/{evidence,rs}` and the top-level
phase and contract documents, written before their deletion. Format per row: quantity, measured value,
tolerance or comparison, source file. Unit conventions: energies in Ry unless stated; U, J in eV.
Scope of all results: fixed Liechtenstein U/J, d shell, bcc Fe / fcc Ni / B2 FeCo fixtures.

## 1. Bulk KS and RS, nsp=1 and nsp=2

How: independent SymPy/Wigner references for the functional, finite-H operator probes with direct
eigenpair or Lorentzian references, no-U byte-identity runs against a pinned baseline, and finite-difference
derivative checks. Controls that must miss the tolerance were run next to each oracle.

Functional, tensor, occupation, Hamiltonian derivative (nsp=1)
- Coulomb tensor, 2,118 p/d elements vs SymPy: max residual 2.09e-17, tolerance 2e-12 (PHASE1.md).
- d-shell energy, U=1.4 Ry, J=0.22 Ry: Eint 6.530967462948, EDC 5.334218750000, EHub 1.196748712948,
  Tr(VP) 0.402247425896 asserted against the SymPy reference. Tr(V dP) vs central differences, steps 1e-3, 1e-4, 1e-5 (PHASE3.md).
- Complex Hermitian P: RS and KS one-rank vs two-rank MPI difference 1.11e-16; scalar RS U=2 `Fe_out.nml` byte-identical at 1 and 2 ranks (PHASE2.md).
- VAL-04 sweep, U=0/2/4 eV: eight-step U=2 residual falls from 2.300011e-2 to a minimum of 4.693231e-3 (finite-iteration check, not converged) (PHASE1.md).
- Baseline: 112/112 unit tests, 181.79 s; Debug MPI 116/116. Reciprocal reference replay: 10 of 18 reference values differ at the pinned base and at current source
  (EF -0.0782523237 vs reference -0.0463779468 Ry; moment 2.000240 vs 2.002214). Not attributed to the U changes (REMEDIATION_GATES.md).

No-U identity
- Single-thread RS no-U output 27/27 artifacts and KS 24/24 artifacts match the registered baselines, repeated at several tips (C6_2_RESULTS.md, C6_4_RESULTS.md, C6_6_T1_RESULTS.md, C6_10_RESULTS.md).
- In C6_6_T1_RESULTS.md the RS comparison against the original C5.2 reference differs in 12/27 artifacts; it matches 27/27 against a refreshed oracle. Cause of the 12/27 difference: unresolved.
- Explicit U=J=0 vs no-U: byte-identical hashes for one-step bcc Fe RS and KS (RA_NO_U_IDENTITY.md); two-step Ni(001) surface, Diff 7.71e-8 in both (C6_7_ZERO_U_COMPARISON/RESULTS.md).

Local-moment (radial) kernels, nsp=1
- Finite-H extraction oracles, one site: KS M1 error 7.76e-18 (limit 5e-12); RS vs analytic Lorentzian 2.89e-7 (limit 3e-4);
  RS vs 6000-channel Simpson 1.96e-17; trace(M1)-q0*c, KS 1.11e-16 and RS 8.88e-16 (limit 1e-10); wrong-operator controls 7.2e-2 (limit above 3e-3) (m0_1_finite_operator.json).
- Local correction vs direct eigenvectors, two-site RS: 2.78e-16; host-site swap control at least 3.84e-4 (C5_2_PART_D_RESULTS.md).
- KS carrier test, accepted U=2, J=0 state (iteration 32): |dSUMEV - sum t| = 1.585921e-3 against allowance A = 1.702238e-3; q0 up/down 4.5276/2.0841 (a2_m03_report.md).
- RS carrier test: |dSUMEV - sum t| = 1.861597e-3 against A_RS = 1.917943e-3; min q2 local/spectral 1.0757e-2/9.745e-3 Ry^2 (LOCAL_MOMENTS.md, supplemental RS M1 verdict).
- Scalar shift scaling (A1.1, KS): fitted zero-shift relative slope intercept a0 = -7.865e-6 against limit 1e-6, so that test failed; q0 = 6.498480 (a11_scaling_summary.md).
- RS C4 RC1: W1 derivative residual -0.171244% (h=0.005 Ry, accepted endpoints), -0.158988% (h=0.0025 Ry, provisional) (C4_RC1_RESULTS.md).

Fixed U/J on bulk, converged states (RB records, Fe, d up / d down / moment muB)
- RS local, U=2/J=0: 4.615579 / 1.895454 / 2.558681; true KS: 4.612752 / 1.892685 / 2.546925 (C1_RB_RESULTS.md).
- U=4/J=0: RS 4.746883 / 1.743735 / 2.785370; KS 4.753320 / 1.728713 / 2.808594. U=2/J=0.9: RS 2.198745 vs KS 2.366575 muB moment.
- U=0 KS-minus-RS moment difference -0.100339 muB; the U2/J0.9 difference-in-differences is +0.268170 muB.
- Two KS cases hit the 183-step cap without a Converged marker; the U=0 KS continuation kept beta 0.05 (deviation recorded).

nsp=2 bulk (SOC insertion and spinor kernels) (SPINOR_DM_CONTRACT.md)
- Spinor P, M1 vs direct eigenpair sums: KS 1.4e-17 and 2.3e-16 (tol 5e-12); RS P vs Lorentzian 4.2e-5 (tol 1e-4).
- Spinor V insertion (`UnitLdaUSocInsertion`, bcc Fe U=2, J=0.5): RS ee 4.3e-19, KS H(k) 1.4e-17 and 2.8e-17 (limit 1e-14); controls 1.9e-2 to 3.4e-2.
- Spinor local kernel vs direct b' definition: c' 1.1e-16, q2' 1.4e-17 (tol 1e-12); q2' >= 0 on 40 random states, minimum 0.
- SOC to 0, KS bcc Fe U=2 (`soc_scale=0` plus `fix_soc`), nsp=2 minus nsp=1: d up 9.4e-10, d down 5.2e-10, moment 1.4e-9 (bound 1e-7); SOC-on control 1.27e-3 to 2.40e-3.
- Without `fix_soc` the same comparison differs by 2.4e-4 (U=0) and 1.2e-3 to 2.4e-3 (U=2); this was the real SOC effect, not a spinor-path error.
- RS block recursion, 200 iterations, `fix_soc`: per-iteration moment, EF, Diff differ by 0.0 at U=0 and U=2.
- RS local mode U=2, J=0, converged: nsp=2 minus nsp=1 d up 1.48e-7, d down 1.21e-7, moment 2.68e-7 muB (bounds 1e-6, 1e-6, 2e-6).
- Real SOC effect, moment from ql, nsp=2 minus nsp=1: RS U=0 +2.3e-4, RS U=2 -5.5e-3, KS U=0 -2.8e-4, KS U=2 -2.35e-3 muB.
- RS smoke, SOC on, U=2/J=0, local: orbital moment 0.0688 muB (0.0979 in spectral mode); max |P_up,dn| 1.294e-2; the three-update R05 streak is first met at iteration 49.
- Iteration-by-iteration RS agreement stops at iteration 12 (Diff 6.5e-9); a 1e-13 relative seed change gives 7.3e-9 at iteration 3, so that criterion was retired.

## 2. Impurities

- Energy diagnostics, B2 FeCo `calctype='I'`: E_int, E_DC, E_Hub vs separate evaluator, max difference 8.891e-13 Ry, limit 1e-12 (RA1_RESULT.md).
- Provenance guards: 13/13 fixture scenarios met (RA2 guards); nsp=2 impurity decks stop at the SOC guard and, with local mode, at the `requires nsp=1` guard (RA2_RESULT.md, RA3_RESULT.md).
- Frozen host and embedding, C3 Part A (Fe): host-site on-site blocks equal bulk to 0; hopping difference 0; bulk EE references agree to 1.17e-15 Ry (C3_PART_A_RESULTS.md).
  The rc/lld scaling run was an invalid setup (bulk non-finite at iteration 184), so no finite-cluster trend was obtained.
- RA4 impurity-minus-bulk records (|M| difference, muB): Fe U=0 9.05e-6; Fe U=3.8/J=1.1 -6.64748e-3 (ratio to U=0: 734.9, target below 10 not met);
  fcc Ni U=3/J=0.7 +3.12e-8 (ratio 0.17 to 0.25 of U=0). Fe multistep discrepancy: unresolved (C1_RA4_RESULTS.md).
- RD2 frozen Ni replay: carrier ratio 0.99142992 per spin, q2' non-negative (C5_4_FINAL_STATUS.md). The original C5.2 replay failed; the C5.4 offline fit changed that verdict.
- RD3(a) Ni impurity vs host ratios: d up 0.64463, d down 0.59821, moment 0.80695. RD3(b) unevaluated: no accepted Ni KS endpoint (C5_4_FINAL_STATUS.md).
- Fe in KS-local host, RS impurity smoke (C5.12 run 2, iterations 81-100): max Diff 5.71e-8 (limit 3e-8 for acceptance, 1e-7 for smoke); d spreads 1.95e-7 and 8.27e-7;
  moment spread 7.66e-7 muB. It missed the C5.1 acceptance limits and met the C5.13 smoke limits (C5_12_AUDIT.md, C5_13_RESULTS.md).
  Impurity minus host: d up +0.019297, d down -0.030948, moment +0.201405 muB.
- nsp=2 impurity oracle (`UnitLdaUSpinorImpurityFrozenH`, frozen host + ES + two U layers): c', q2' vs spinor formula 0 and 0 (1e-12, 1e-10);
  vs direct eigenpairs 1.5e-5 and 2.1e-5 (3e-4); controls 4.3e-2 to 2.4e-1 (SPINOR_DM_CONTRACT.md).
- Impurity energy (RD4): not measured; the W1 run crashed (SIGSEGV) before SCF. Impurity `auto` reproduction at a KS host: unevaluated (C5_4_FINAL_STATUS.md).
- C5.7 KS Fe U=2 `auto` runs: both exited 139 (segmentation fault), cause not found in sources (C5_7_RESULTS.md).

## 3. Surfaces

- Fixture: fcc Cu(001) `calctype='S'`, `nlay=2`; Ni(001) decks derived from it (alat 3.5234, fermi -0.061141 Ry) (C6_T0_RESULTS.md, C6_2_RESULTS.md).
- Surface auto-local resolver extended to `S` at nsp=1: `UnitAutoImpurityLocalMode` 1/1; surface Hubbard energy routing stays off (C6_T1_RESULTS.md).
- Ni(001) two-step RS probe, d occupation up/down at the surface: U=0 4.231290/4.231290; U=3/J=0.7 local 4.265438/4.265438. Surface minus bulk d change with U: -0.006188 e (C6_6_SURFACE_PROBE_RESULTS.md).
- Two-step Diff: U=0 control 0.0162421, 0.384026; U=3/J=0.7 0.0163260, 0.386918, against the requested 5e-8, so the response criterion failed (same file).
- U=4/J=1 on the surface layer only (C6_10_RESULTS.md): surface Ni moment 1.081094 then 0.112752 muB; U=0 control 0.689004 muB (stable).
  Diff 0.1149764639 then 0.0389620957; the finite-U moments were still changing, so these are response values only.
- Some U>0 surface runs hit NaN projected DOS and `Invalid local d-shell Hubbard moment correction; status=5` (C6_4_RESULTS.md, C6_7, C6_8). The wide-window U=4/J=1 case gave n_up-n_down 0.003501 at the surface after step 1.
- nsp=2 surface with SOC to 0 vs nsp=1, same host numbers (SPINOR_DM_CONTRACT.md, S6c): reduced surface 4-step regression, 0.0 over 112 quantities (tol 1e-7);
  full-size iteration 1, ql 6.9e-18, Diff 3.1737e-6 in both (tol 1e-9); against an independently converged nsp=1 host, ql 2.5e-10 (tol 1e-9).
- Surface operator test `UnitLdaUSocSurfaceInsertion`: on-site ee difference 1.4e-17, 2.2e-19, 1.9e-34 (limit 1e-14); off-site and lsham differences 0; controls 2.0e-2 to 5.0e-2.
- nsp=1 host guard: an nsp=2 surface with a frozen U host written at nsp=1 stops with "was written at nsp=1, run nsp=2. Regenerate the host at nsp=2." (`Example_surface_soc_u_nsp1_host_guard`).
- SOC-on short runs from nsp=1 seeds: Ni-S orbital moment 0.066 to 0.109 muB in 5 steps and 0.118 to 0.175 muB in a 20-step run; Diff 1.2e-4 to 2.7e-4; J=0 diagnostic did not remove the drift (surface SCF instability, undiagnosed).
- Zero-U Ni(001) one-step: spin moments atoms 1-3 0.007355, 0.689005, 0.640933 muB, identical for no-U and explicit-zero-U.
- Surface stage state at the end: paused; no accepted converged no-U Ni(001) or Fe(001) nsp=1 reference was available (C6_5_RESULTS.md).

## 4. Tilted axes and the RS floors

Tilted and mixed spin axes
- Source inventory only in the docs: axis resolver, frame guards (`Example_bulk_bccFe_soc_u_tilted_axis_guard`, `Example_bulk_bccFe_soc_u_tilted_block_guard`), per-site axes. Tilted-axis block recursion is fatal (SPINOR_DM_CONTRACT.md).
- Measured tilted-axis number: iteration 1 only, tolerance 1e-9 stated in TESTING_POLICY.md; value not found in sources (it sits in the deck README `decks/ni001_soc_u3/README`, outside this tree).
- Mixed-axis (different axes on different sites) results: not found in sources. The doc notes nothing enforces a common axis at nsp=2.
- The spinor kernel uses global-z spin projectors; for an axis of -z the channels are swapped (read from code, not run). The test seed for tilted axes is a rotated state (rotate_spin_frame.py); a collinear z seed along x relaxed slowly (TESTING_POLICY.md).
- Surface rotation floor at a converged state: open (SPINOR_DM_CONTRACT.md).

RS floors (bcc Fe, U=0, RS block, Broyden)
- Convergence plateau: final-quarter (steps 181-240) max Diff 2.67e-8, mean 9.06e-9, OLS slope 7.19e-11 per step; requested conv_thr 1e-9 was not reached; RS-only threshold 3e-8 (a2_rs_floor_3e8.json).
- Carrier floor D0_RS = 5.0565e-6 Ry (|Tr M1 - SUMEV|). Uncertainty u_RS = 6.6267e-8 Ry (channel refinement 2500 to 6000: D -5.0565e-6 to -4.9902e-6; repeat 4.44e-9) (a3_rs_floor_uncertainty_results.json).
- KS floor D0 value: not found in sources.
- RS fixed-U states: accepted window Diff 1.56e-8 (max 2.7e-8) for nsp=2; nsp=1 9.9e-9 (max 1.7e-8) (SPINOR_DM_CONTRACT.md, S4).
- RS nsp=1 U=2 trajectories are seed-sensitive: 1e-13 seed change gives Diff differences 7.3e-9 (iteration 3) and 1.1e-6 (iteration 20).
- RS no-U reference at the time of the plan did not converge (PROGRESS.md; later floor runs above).

## 5. Energy

- Candidate b (local base free energy plus production E_Hub on direct P_OUT), native KS, Ry/eV, against the fixed-P target (a5_m22_u2_derivative_result.md, a5_m23_j_derivative_result.md):
  - U=2, J=0: D_b +0.0550977, target +0.0552397, difference -1.42e-4, allowance interval [1.02e-3, 1.55e-3].
  - U=2, J=0.5: difference -3.73e-4, interval [1.19e-3, 1.75e-3]. U=2, J=0.9: difference +1.25e-4, interval [1.18e-3, 2.58e-3].
- Candidate c (spectral) misses: -0.14180 (U=2), +0.12617 (J=0.5), +0.10496 (J=0.9); about 114, 30 and 92 times the upper allowance. Candidate a is report-only.
- RS clean bulk, candidate b mean error vs mean allowance A (C4_RC2_RESULTS.md): J=0: 2.31e-4 vs 4.86e-3; J=0.5: 3.20e-4 vs 8.32e-3; J=0.9: 1.45e-4 vs 5.69e-3.
  Candidate c errors 0.1404, 0.1171, 0.1060, each above 10 A.
- Earlier RC attempt (C3_RC_RESULTS.md): at J=0.5, h=0.05 the b error exceeded A by 4.45e-7 Ry/eV, which stopped that campaign; the C4 allowance formula was later amended.
- Base-energy W1 non-stationarity: -0.14% KS, -0.16% to -0.17% RS (not reduced by the new route). Spin-selective response matrix determinant 2.58 nominal, allowed range -1.47 to +6.63 (a5_spin_w1_results.md).
- nsp=2 energy (SPINOR_DM_CONTRACT.md, S5): SOC to 0 candidate b difference 5.8e-10 Ry (gate 1e-7), E_Hub 2.5e-13; W1 dF/dlambda -0.144% KS, -0.163% RS (gate 1%);
  U derivative D_b 0.05508311 vs target 0.05522711 (KS), c miss 52.8 A_U; J derivative KS D_b -0.06389993 vs -0.06351369; RS J not tested.
- Routed energy: RS nsp=2 base -2541.981867214201 plus E_Hub 0.1144260991844432 gives -2541.867441115017 Ry (SPINOR_DM_CONTRACT.md). KS nsp=1 U=2/J=0 candidate b -2541.8621841266 Ry (a4_3_like_for_like_observables.md).
- Core algebra (PHASE3.md): the radial path response is nonlinear, so the SUMEV double-count coefficient was not identified in the pre-A5 phase; the A5 comparison above is the scoped numerical resolution.
- Matched serial/MPI reproducibility and restart checks for the routed energy: indexed in evidence/a4_4_control_evidence_index.json; values not found in sources read here.

## Status and limits

- Excluded from the energy route and not tested: SC-U, intersite V, nsp=3/4 with U beyond the symmetry checks of section 6, impurity U, impurity/surface hosts for energy, nonzero test fields, constrained SCF.
- Energy routing is off for surfaces and impurities (RD4 not measured).
- No energy-anisotropy E(x)-E(z); no converged SOC-on surface; no converged nsp=2 surface vs converged nsp=1 surface.
- R06 (converged scalar Fe U/J sweep), R07 (noncubic J), R09 SOC ladder, R10, R11, R12, R13, R15: NOT_RUN in REMEDIATION_GATES.md; R02 aggregate NOT_RUN; R14 aggregate NOT_RUN.
- Energy results hold to the ASA base-functional accuracy of about 0.14% to 0.17% stationarity residual; no a-versus-b preference was claimed; no general radial energy identity was derived.
- hcp Co impurity: the pre-existing non-finite failure was not resolved. Fe U>0 impurity vs bulk ratio target below 10 not met.
- Surface Ni(001) response gate (Diff <= 5e-8 in two steps) failed; the threshold was later declared unsuitable by the owner and no replacement was set.

Raw evidence is at git tag `archive/fable_v4b_u-2026-11-05` under docs/validation/fixed_ldau/{evidence,rs}.

## 6. HOH and CCOR_2c (2026-11-07, branch main_rc_ldau_hoh)

- V_U enters as an on-site term in `enim` with `hoh`, in `ee`/`hall` without. Fixed U=2/J=0 eV, bcc Fe (`tests/scf/cases/bulk/bccFe_hubbard_nsp1` seed), one step from the repo deck,
  E_F(U=2) - E_F(U=0), Ry: `hoh` block recursion +0.0697 (lld 40), Chebyshev +0.0682 (window +-4 Ry, 600 moments) and +0.0679 (+-3 Ry, 300 moments), k-space second order +0.0677;
  no `hoh` Chebyshev +0.0665. Spin moment after step 1, U=0 -> U=2 (muB): `hoh` block 2.182 -> 2.147, Chebyshev 2.193 -> 2.166.
- NiO AFM-II, 4-atom cell, `hoh` + `ccor_2c`, U=2/J=1 eV on both Ni, one step from the same potentials: E_F -0.130883 -> -0.103250 Ry (real-space Chebyshev), -0.13393 -> -0.10789 Ry (k-space second order, 8^3);
  Ni moment 0.7117 -> 0.7834 muB (real space), 0.7741 -> 0.8355 muB (k-space). Five real-space steps: Ni moment 0.78, 1.07, 1.19, 1.21, 1.11 muB, Diff 4.8e-3, 1.1e-2, 1.3e-2, 1.3e-2, 8.3e-3; E_F -0.103 to -0.183 Ry, not converged. U=0 baseline, 6 steps: Ni moment 0.711 to 0.721 muB, Diff about 3.0e-4.
- The Fe seed has no `vrmax`, so `ccor_2c` on it changes no digit of the step-1 result: those decks do not test CCOR. The default kernel, `cheb_backend='fast_dp'` and `'legacy'` agree to the printed digits for `ccor_2c` + U (U=2: E_F -0.001913 Ry, moment 2.185995 muB).
- The code rewrites `fermi` in `input.nml` after a run and step 1 depends strongly on it; a deck copied from a run directory starts from a different state.
- Not run: energy accounting with `hoh`, spectral-vs-local comparison with `hoh`, impurity or surface with `hoh`, no-`hoh` `ccor_2c` with U on a seed that has `vrmax`, an external-code comparison.
- Derivative test of the spectral energy candidate (2026-11-07, KS bcc Fe `tests/scf/cases/k_space_scf/bccFe`, nsp=1, 12^3 k, 300 K, U=1.8 to 2.2 eV, J=0, every run converged, final Diff <= 1.4e-8, built from the pristine deck).
  Central difference of the free energy against the fixed-P Hubbard derivative E_Hub(P_OUT(U0))/U0, Ry/eV, U0 = 2 eV. `local` candidate b: D = +0.055098, target +0.055240, miss -1.4e-4 (documented: +0.0550977, +0.0552397).
  `spectral` candidate a (F - Tr(V_in P_OUT) + E_Hub): D = +0.064771, target +0.057828, miss +6.9e-3 (step 0.2 eV); +6.6e-3 (step 0.4 eV); at 16^3 k: +6.0e-3, with `local` b -1.8e-4. `spectral` candidate c (without the trace term): miss -0.141.
  The documented allowance at U=2, J=0 is [1.0e-3, 1.6e-3]. At U=2 eV the spectral free energy F_a lies 8.7 mRy above the local F_b (7.7 mRy at 1.9 eV, 9.7 mRy at 2.1 eV). Spin moments: spectral 2.356302 muB (matches the earlier record), local 2.510777 muB.
  The spectral state is therefore not a stationary point of the functional that the local state satisfies to 1.4e-4 Ry/eV; spectral energy routing is not adopted. Not tested: RS, J != 0, `hoh`, other materials.
- `hoh` with `nsp=2` and SOC (spinor U, `tests/scf/cases/bulk/bccFe_hubbard`, Chebyshev, one step, U=2 J=0): E_F(U=2) - E_F(U=0) +0.0680 Ry with `hoh`, +0.0668 Ry without; orbital moment 0.0425 muB, max |P_up,dn| 1.3e-2.
  Only the spin-diagonal placement of V_U is tested; an earlier guard (`hoh` requires `nsp=1`) was removed. Not tested: the spin-off-diagonal part of V_U with `hoh`, `ccor_2c` at `nsp=2`, tilted axes.
- Symmetry checks of `nsp=3`/`4` with U (`tests/ldau/test_nsp34_noncollinear.py`; no external reference; bcc Fe, U=2 eV, J=0, Chebyshev unless stated; tolerance 2e-6 Ry or muB): ferromagnet with the moment along z against (0.612, 0.612, 0.5): E_F and |m| differences 0 to the log precision, spinor-P eigenvalues differ by 8.8e-12, |P_up,dn| in the global frame 0.297 for the tilted run;
  two-atom AF with moments (z, x) against (x, -z) and with the sites exchanged: E_F and site |m| equal over 2 steps; two parallel moments at `nsp=3` against the `nsp=1` ferromagnet: equal; `nsp=4` with SOC, moment along z against x: 1.0e-6 (E_F, |m|); SOC off at `nsp=4` against `nsp=3`: equal;
  k-space second order `nsp=3` AF, |m1| - |m2| over 3 steps: 0 (|m1| 2.385568, 2.632155, 2.652151 muB); restart from `_out.nml` + `_dm.nml` after 2 steps against the continuous step 3: E_F -0.052606 Ry and |m| 2.526930 muB on both.
  With bugs injected locally (not committed): spin-off-diagonal V_U blocks dropped: rotation E_F difference 5.0e-3 Ry, |m| 1.6e-1 muB, AF angle 6.4e-1, SOC cubic 2.2e-1; U skipped on site 2: k-space mirror 7.9e-1 muB; spin coherences not read from the DM file: restart E_F -0.051938 against -0.051902 Ry.
  Not covered: block recursion at `nsp>=3` with noncollinear axes (not rotation invariant, KNOWN_ISSUES), k-space `nsp=1` AF (breaks the site symmetry, KNOWN_ISSUES), `local_axis`, constraints, impurity/surface.
