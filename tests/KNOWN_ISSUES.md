# Known issues found via Phase 2 test coverage

Bugs surfaced while closing test-matrix coverage gaps (Phase 2, P1). Recorded
here rather than fixed in place, per the "no further structural refactoring"
rule for this phase — each entry is a candidate for a future bug-fix task.

## `recur = 'lanczos'` + `nsp = 2` produces NaN DOS and `lmom`

- **Symptom:** with `control%recur = 'lanczos'` and `control%nsp = 2`, every
  row of `totaldos.out` is NaN in the DOS column, and `lmom` in the output
  namelist is NaN. `etot`, `ws_r`, and `mom` are unaffected (physically
  sane, match the `nsp=1` lanczos and `nsp=2` block/chebyshev runs).
- **Scope:** reproduces on both `hamiltonian%hoh = .true.` and `.false.`.
  `nsp=1` lanczos is unaffected (`totaldos.out` and `lmom` both sane).
  `nsp=2` with `recur = 'block'` or `'chebyshev'` is unaffected.
- **Likely area:** `green%sgreen()` (the lanczos-path DOS/Green's-function
  routine, dispatched from `calculation.f90`'s `recur` select-case) or the
  orbital-moment calculation in `bands.f90` — both paths presumably assume
  a code path shared with `nsp=1` that doesn't generalize to two spin
  channels.
- **Test impact:** `tests/scf/cases.json` entries `Example_bulk_bccFe_nsp2_lanczos`
  and `Example_bulk_bccFe_nsp2_lanczos_hoh` check only `etot`/`ws_r`/`mom`;
  the `totaldos.out` and `lmom` checks are intentionally omitted until this
  is fixed, so the reference doesn't pin down an all-NaN result.
- **Found:** Phase 2, P1 (adding lanczos coverage to the SCF example suite —
  this combination had zero prior test coverage under any phase).

## `frozen_magnon` `branch_mode = 'auto'`: multi-sublattice acoustic magnon not gapless at Γ

- **Symptom:** the multi-sublattice magnon branches from
  `post_processing_frozen_magnon_auto` (`&frozen_magnon branch_mode = 'auto'`,
  `calculation.f90`) do not reproduce the Goldstone theorem for systems with
  **inequivalent** magnetic sublattices: the acoustic branch has a finite
  `omega(Γ)` (~0.28 Ry on a two-sublattice bcc FeCo k-space test) instead of
  going to zero, and the dispersion is nearly flat.
- **Scope:** the **single-sublattice** limit is correct — `omega(Γ) ≈ 0`
  (naturally, not enforced) with a clean quadratic dispersion, on the
  real-space recursion path. The failure is specific to ≥2 inequivalent
  sublattices; it appears on the k-space path (the only one exercised so far
  for multi-sublattice).
- **Method (for reference):** the auto-branch implements the direct GBT
  frozen-magnon method (Essenberger et al., PRB 84, 174425 (2011), Eq. 26;
  Sandratskii, Carva & Silkin, PRB 111, 184436 (2025)): the magnon matrix is
  the second derivative of the frozen-magnon energy surface w.r.t. sublattice
  cone angles, evaluated with the magnetic force theorem (band energy at the
  fixed reference potential), and the magnon energies are the eigenvalues of
  the real symmetric matrix `√(M_μM_ν)·Re[J̃_μν^q]`. This construction gives a
  gapless acoustic mode at Γ **iff** the band energy is invariant under a
  global (uniform) spin rotation.
- **Likely root cause:** the reciprocal band-energy evaluation
  (`reciprocal%calculate_band_energy_from_moments`, via
  `build_kspace_hamiltonian`/`diagonalize_hamiltonian`) is not exactly
  invariant under a uniform rotation of all sublattice moments — suspected
  contributors are the per-probe Fermi-level re-determination
  (`auto_find_fermi = .true.`) shifting the band-energy zero, or a
  moment-projection term in the band-energy sum. Diagnostic not yet run:
  compare `E_ref` against the uniform-tilt pair energy `E_{12}(Γ)` (should be
  equal), and run the same case on the real-space recursion path to isolate
  k-space vs. general.
- **Test impact:** `tests/scf/cases.json` `Example_frozen_magnon_bccFe_auto`
  and `_auto_scf` are single-sublattice smoke cases only (no committed
  reference); they exercise the code path but do not pin multi-sublattice
  values. The plain acoustic `Example_frozen_magnon_bccFe` (single-branch
  flat-spiral sweep) is unaffected and is the validated `frozen_magnon`
  deliverable.
- **Status:** deferred deliberately. The multi-branch magnon spectrum is the
  validation target for blueprint **B11** (linear-response TDDFT / transverse
  magnons); this is the natural place to resolve the global-rotation-invariance
  question. See `docs/dev/B1_GBT_SPIN_SPIRAL_PLAN.md` (T5) and commit
  `d86fe42`.

## `processing = 'sd'` (spin dynamics) has no working pre-processing path

- **Symptom:** `calculation%processing_sd()` (`calculation.f90`, dispatched
  when `&calculation` sets `processing = 'sd'`) always runs the sequence
  `build_data → bravais → build_surf_full → newclu → structb`, regardless of
  the case's own `pre_processing` setting. That sequence is byte-for-byte
  identical to `pre_processing_newclusurf`'s (impurity embedded in a
  surface) — it is not the `bravais` (plain bulk) or `buildsurf` (plain
  surface) sequence those `pre_processing` values normally trigger elsewhere
  in `calculation.f90`.
- **Reproduction:** patching a plain bulk case's `&calculation` with
  `processing = 'sd'` and adding an `&sd` block crashes inside
  `build_surf_full` (`lattice_cluster.f90:154`, "Bad integer for item 1 in
  list input") — `build_surf_full` expects surface-specific namelist keys
  (`&surface`-style: `nlay`, `surftype`, ...) that a bulk input doesn't
  have.
- **Blocker:** no example or test input anywhere in the repository uses
  `pre_processing = 'newclusurf'` (impurity-in-surface) — that whole
  pipeline has zero precedent, so there is no known-good namelist shape to
  build a spin-dynamics test case from without first constructing (and
  validating) a `newclusurf` input from scratch.
- **Likely root cause:** `processing_sd`'s pre-processing block looks like
  it was copy-pasted from `pre_processing_newclusurf` and never generalized
  — spin dynamics conceptually should run on top of whichever
  `pre_processing` route (`bravais`/`newclubulk`/`buildsurf`/`newclusurf`)
  the case already uses, matching the same sequence
  `pre_processing_bravais`/etc. use, not a fixed one.
- **Test impact:** P1's requested `Example_bulk_*_sd` smoke case is not
  added — there is currently no reachable code path to exercise. Revisit
  once `processing_sd`'s pre-processing is fixed to reuse the calling
  case's own route (likely sharing code with the corresponding
  `pre_processing_*` subroutine rather than hardcoding one).
- **Found:** Phase 2, P1 (adding spin-dynamics coverage — `processing_sd`
  had zero prior test coverage under any phase).
