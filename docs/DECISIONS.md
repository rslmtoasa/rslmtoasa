# RS-LMTO-ASA decisions

One line per closed campaign or decision: what it concluded, and where the
evidence lives. Campaign reports are not kept in the tree.

Archive tags (fetch them with `git fetch --tags`):

- `lr-campaign-archive` (`70eb73b`): the campaign reports. A path written
  `lr-campaign-archive:<path>` resolves there, and at any commit up to `592ebc4^`.
- `lr-campaign-archive-2026-10` (`e95a17c`, the last `fable_v4` commit): the
  linear-response source, tests and `docs/linear_response/` removed from
  `fable_v4b`.
- `spin-response-baseline-v1`: the verified B11 baseline (Stage 2).
- `spin-response-stage3a`: the B11 Stage 2 and 3a records and scripts
  (`docs/validation/VAL-18_*`, `VAL-19_*`, `tests/validation/spin_response_stage2/`,
  `spin_response_stage3a/`).

Until 2026-10-11 this file held one row per campaign report; that version is in
the history of `docs/DECISIONS.md` at `2477858`.

## Campaigns before 2026-08

- 2026-06 to 2026-08 — ACC (accelerators): CPU baselines recorded; the
  real-space CUDA recursion plugin passed its route comparisons; a reciprocal
  CUDA backend with a strict-Lehmann contraction passed its comparison; H(k)
  assembly kept on the host; KPM transport on CUDA complete. Reports at
  `lr-campaign-archive:docs/dev/ACC-*`, `RS_LMTO_ASA_PHASE_III_A_*`,
  `RS_LMTO_ASA_ACC_PERFORMANCE_RESCUE.md`.
- 2026-07-16 — B2 / RF05 (k-space Green functions): B2 completed, intersite
  normalization resolved, both gates signed; typed reciprocal execution
  backends with a LAPACK factory. Reports at
  `lr-campaign-archive:docs/dev/B2_RECIPROCAL_GREEN_HANDOVER.md`,
  `RF05_RECIPROCAL_EXECUTION_BACKENDS.md`.
- 2026-08-03 to 2026-08-08 — GBT (generalized Bloch theorem): classified as
  prototype evidence. G2E, G5 (recursion bound `lld ≤ 28`), G7 and G8 passed;
  the production operator slice G2O failed; G6 stayed open; G9 failed
  (cone-angle scaling, band-energy residuals), so G10 could not close.
  Reports at `lr-campaign-archive:GBT_RS_LMTO_completion_blueprint.md` and
  `docs/dev/GBT_WP*`.
- 2026-08-27 — GBT re-audit: fixed-potential operator gates passed for the
  tested scope (q = 0, shifted k, commensurate supercells); composite terms
  passed except CCOR, SOC, Hubbard-V, local axes and constraints; frame and
  constraint contracts defined; corrected constrained MFT for the acoustic
  branch; converged small-q curvature not established; q = 0 SCF closure passed,
  finite-q constrained SCF unqualified. Reports at
  `lr-campaign-archive:docs/dev/GBT_REAUDIT_WP*`.
- 2026-08-20 — KPM transport: CPU and GPU benchmark rows passed their
  correctness checks; the CUDA reconstruction was optimized and the benchmark
  harness consolidated. Reports at `lr-campaign-archive:docs/dev/RS_LMTO_ASA_KPM_*`.
- 2026-08 — SCF benchmarks (B0C, B1, B1R, B1R2): a shared CPU/GPU SCF benchmark
  harness; the lean desktop tier completed 18 of 18 cases; real-space against
  k-space accuracy remained inconclusive. Reports at
  `lr-campaign-archive:docs/dev/RS_LMTO_ASA_SCF_*`.
- 2026 — Phases I–III: Phase I stabilization complete; Phase II validation kept
  scoped, with no broad GBT or TDDFT promotion; Phase 3 status recorded.
  Reports at `lr-campaign-archive:docs/dev/PHASE*`.
- 2026-07-03 — Refactoring: no-physics-change rules; an independent clean build
  and regression matrix passed 8 of 8 runnable cases; uncalled math routines
  removed. Reports at `lr-campaign-archive:docs/dev/REFACTORING_PLAN.md`,
  `REVIEW_NOTES_T1-T6.md`, `MATH_AUDIT_T11.md`, `REFACTORING_PHASE2.md`.
- 2026-08-30 — XC: legacy and libXC scope accepted, unsupported requests fail
  explicitly; radial GGA reconciled; the LDA Barth–Hedin derivative and
  zero-spin-channel handling repaired (commit `80390cc`, which shifts bcc-Fe
  total energies by 6e-6 to 1.5e-5 Ry); a converged high-spin fcc-Fe magnetic
  SCF not established. Reports at `lr-campaign-archive:docs/XC_*`.

## Linear-response (LR) campaigns, 2026-08-09 to 2026-10-03

Many PASS or "certified" verdicts in these campaigns rest on checks that held by
construction (one-site cells, q = 0 identities, residuals aliasing each other).
None of them is validated physics, and none is a reference for new work.
Reports at `lr-campaign-archive:docs/DRESP_*`, `TG_FZ_*`, `LR*`, `TDDFT_*`,
`RSGF_*`, `GOLDSTONE_EIGENVALUE_CORRECTION.md`; the condensed formulation
documents at `lr-campaign-archive-2026-10:docs/linear_response/`.

- 2026-10-03 — LR-METHOD-05R: SR ALSDA kernel repaired; blocked, the
  response-representation closure not established (Fe raw Ward relative
  residual ≈ 0.243).
- 2026-10-03 — LR campaign closed (`fable_v4`, about 175 source commits): no
  validated magnon dispersion. External audits of the earlier `tddft_*` code
  found self-referential gates, a moment-sign normalization error and
  Fermi-level contract violations. Code archived at
  `lr-campaign-archive-2026-10`; B11 restarted on `fable_v4b`.

## Spin response restart (B11), from 2026-10-03

- 2026-10-06 — Baseline (Stage 2, tag `spin-response-baseline-v1`): bcc-Fe
  Juelich-d k-space transverse response, verified by the analytic suite
  (A1–A3), contract tests (C1–C4), S2, S3 (three-atom cell, max difference
  5.7e-13) and a 24³ reproduction snapshot; mutations 1–19 each fail at least
  one test. A folded k+q needs a site-dependent Bloch phase for several sites
  (`tests/KNOWN_ISSUES.md`, resolved). 60³ scans against LKAG + LSWT: peak
  deviations −0.284 to +0.134 for ξ ≤ 0.4333 along Γ–H; projector sensitivity
  −16.4 %. Records at tag `spin-response-stage3a`.
- 2026-10-08 — Stage 3a (stiffness): the Juelich χ⁰ zeroth moment is 2.349
  against the Juelich moment M = 2.068, and the q = 0 pole weight is W = 2.545;
  with ω = U·W·δ_q, the raw stiffness at 160³ is D(ξ = 0.025) = 315.9 / 298.6 /
  278.9 meV Å² at 100 / 300 / 600 K and D(ξ = 0.1) = 271.9 / 268.9 / 260.5. The
  small-q rise of D depends on kT; mesh steps at 100 K are −3.4 % (ξ = 0.025)
  and −2.5 % (ξ = 0.05) for 80³ → 160³, −0.16 % at ξ = 0.05 for 120³ → 160³.
  The tagged VAL-19 table's zeroth moment 0.00223 used a ±0.05 Ry window that
  excludes the transitions. Records at tag `spin-response-stage3a`.
- 2026-10-09 — Stage 3b: `&energy temperature` added (default 0 K, the T = 0
  path bit-identical; Sommerfeld unit test). Feasibility study closed: the
  real-space LKAG far shells (beyond about 3a) are not converged in cluster size
  (rc 80 → 144), recursion depth (lld 31 → 61) or cutoff, so the small-q
  stiffness of the reference is undetermined at affordable cost; the near
  shells and the mid-zone dispersion are robust; finite T = 600 K reduces the
  far shells by about 17 %. No further reference work.
- 2026-10-10 — bcc Fe physical fidelity (120³, η = 1e-3, Juelich-d): converged
  100³ → 120³ to under 1.5 % at ξ = 0.1–0.4; the tr L ridge follows LKAG + LSWT
  over most of the full path (Γ–H ratio 0.90 / 0.79 / 0.96 / 1.03 at
  ξ = 0.1–0.4) and matches the static pole ω_s = U·W·δ_q below ξ ≈ 0.4; steps
  along Γ–H coincide with continuum bands in −Im χ⁰ (ξ ≈ 0.2–0.27, 0.4–0.44);
  Landau-damped beyond ξ ≈ 0.45; a damped Fermi-surface softening near
  ξ ≈ 0.75 (δ_q 0.030 → 0.016) whose depth depends on the projector (the
  coefficient projector gives δ_q < 0, an unphysical instability). δ_q > 0 at all
  90 non-Γ points of the full path; small-q δ_q/|q|² isotropic within 4 %.
- 2026-10-11 — fcc Ni physical fidelity (128³, η = 1e-3, Juelich-d): moment
  0.659 μ_B; static 128³ and 120³ agree to 4e-4; δ_q > 0 at all q ≠ 0 of the
  path; small-q δ_q/|q|² isotropic within 0.7 %; D(W) ≈ 905 meV Å², about 3.2 ×
  Fe, W/M = 1.38 (Fe 1.23); the ridge follows ω_s near Γ and dissolves from
  about a quarter of Γ–X (|q| ≈ 0.45 Å⁻¹). Against an LKAG reference at the
  response's Fermi level (32³), ω_s/LSWT = 1.06–1.15 away from small q, rising
  towards Γ where the reference is unreliable. The two do not share one
  potential: these response runs used the 16³-SCF state (total spin moment
  0.6586 μ_B), the J(r) of the reference comes from the 32³-SCF state (0.6831
  μ_B), and its AMS moment is M = 0.658589 from `momfile` (the 16³ value;
  recomputing ω from `jfile` reproduces the file to 0.9 % or better with 0.658589
  and to 2.7–4.5 % with 0.6831 at 8 of 11 q tried, the three smallest |q| being
  off by 16–40 % with either). The `ni_fcc` deck was moved to the 32³ state
  afterwards (`2ce62be`); the 128³ maps have not been rerun on it.
- 2026-10-11 — Next steps, in order: 3c–3e (s, p, d, f channels; ALSDA kernel
  matrix; KL adiabatic output), then non-collinear and SOC (NC-0a/0b feasibility
  and cross-route comparison with real-space recursion; NC-1 to NC-3), then
  GPU (B4), with LDA+U on its own branch.
- 2026-10-13 — NC-1 closed on bcc Cr (cr_bcc; 32³ and 64³ static, 32³
  dynamic). Structural: p = 0.992 (64³); smallest eig(1+χ0U) at q ≠ 0
  +1.75e-3 (64³), none negative; |ω₊+ω₋|/ω₊ ≤ 4.0e-12; U₁, U₂ equal to
  5.1e-12; spin_swapped T. The adiabatic ω_RPA is 1.65–1.94 × the LKAG + LSWT
  magnons: this is the force-theorem error, not a normalization.
  ω_force = N U (1 + χ0 U) from the same χ0 matches LSWT to 0.968–1.087
  (84 q, max component ≥ 1/8), and ω_RPA/ω_force = 1.9637 at q = 1/32
  against 1/√(1−λ_u) = 1.9657 (λ_u = 0.741); A5 antiferromagnet 1.4092
  against 1.4091. The broad dynamic tr L peak (1.62–1.83 ω_RPA) is not
  examined; site-resolved output deferred.
- 2026-10-14 — NC-2 (SOC, anisotropy gap) stopped as a negative result.
  In the tree: the consistent kernel K = U m0/m(ξ) (fixed exchange field; a
  frozen U would itself create a spurious gap), its ξ = 0 identity test, and
  a warning on every SOC run. Gate 1, cubic bcc Fe null test: the q = 0 gap
  is −2.3e-5 Ry and scales as ξ² (slopes 1.99–2.00, 12³ to 48³), where
  symmetry requires fourth order and far below 1e-8 Ry. Cause confirmed in
  a one-site toy: the projected vertex breaks the Goldstone sum rule through
  its spin asymmetry; the Goldstone-fitted U absorbs that at ξ = 0, but its
  O(ξ²) change appears as a gap (complete vertex −2.6e-17 Ry; gap
  proportional to the asymmetry; with the Fe projector factors −7.1e-5 Ry,
  same sign as Gate 1). Consequence: SOC gaps below about 1e-5 Ry cannot be
  measured with the projected vertex; FePt (1–3e-4 Ry) was not attempted.
