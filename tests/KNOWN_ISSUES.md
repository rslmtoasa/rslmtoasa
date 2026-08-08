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

## [RESOLVED 2026-07-23, commit 8b42928] Exchange `J_ij` NaN — `simpson_f` out-of-bounds read

- **Was filed as:** "Lehmann first-order exchange: latent uninitialized-local NaN
  in `J_ij` (layout-sensitive)", found during B5.3 (2026-07-16). **The
  uninitialized-local diagnosis was wrong** — see the actual root cause below.
- **Symptom (as observed):** `post_processing='exchange'` + `gf_route='lehmann'`
  (or `'dyson'`) could produce `J_ij = NaN` under `-O3`, presenting as a
  heisenbug: clean under `-O0`/`-finit-real=snan`/`-fcheck`, and clean under `-O3`
  until an unrelated `calculation`-type layout change (the B5.3 `do_damping`
  `logical`) flipped it to NaN.
- **Actual root cause:** a **one-element out-of-bounds heap read in
  `math.f90::simpson_f`**, not an uninitialized local, and **not**
  Lehmann-specific — it was latent in the recursion route too, masked by heap
  layout. The Fermi/dFermi branches declared their arrays `dimension(NPTS+10)`
  and looped `I = 2, NPTS+9, 2`, reading index `NPTS+10`. Callers pass arrays of
  length `en%channels_ldos+10`, but `NPTS = en%nv1 = channels_ldos+1` (every real
  input has even `channels_ldos`), so the true extent is `NPTS+9` and the loop
  read one element past the end of every integrand and of `en%ene`. The garbage
  byte read as a NaN bit-pattern under some heap layouts and as a finite/near-zero
  value under others — hence the layout sensitivity and the false "uninitialized"
  signature. `-finit-real=snan`/`-fcheck=bounds` miss it because the read is legal
  against the callee's *declared* `NPTS+10` dummy; only past the *actual*
  allocation. Valgrind on Linux `-O3` pinned it: *"Invalid read of size 8 at
  math.f90:1128 … 0 bytes after a block of size 2480"* (= 310 doubles =
  `channels_ldos+10`), origins `energy.f90::e_mesh` and the `exchange` integrand
  allocations.
- **Fix:** declare the `simpson_f` dummies `dimension(NPTS+9)` (the true extent)
  and cap the Fermi/dFermi loops at `NPTS+8`, so the last index read is `NPTS+9`
  (in bounds). Drops one Simpson triple whose Fermi weight is ~0; integrals shift
  only ~1e-6 (recursion `J_ij` 1_335: 0.5078779 → 0.5078765; lehmann/dyson move at
  machine epsilon; the B5.2 σ triad is unchanged to 6 dp). The `triad_bccFe_jij`
  golden was regenerated on the fixed build; `triad_bccFe_sigma` was unchanged.
- **Confirmed:** on the Linux `-O3` + `pad_dummy` trigger layout (where it
  originally NaN'd), `valgrind --track-origins=yes` went from 1509 errors / 37
  contexts to **1 error / 1 context** — the sole remainder being the unrelated
  pre-existing `clusba` read at `lattice_strux.f90:889` (a separate issue). On
  gfortran-13 the pre-fix and post-fix recursion values are identical, confirming
  the result no longer depends on the out-of-bounds byte.
- **Unblocks:** the deferred B5.3 `do_damping` wiring + Gilbert-damping α triad
  (`docs/dev/B5.3_gilbert_damping_audit.md`) can now re-land — the layout
  perturbation that re-triggered the NaN no longer does.

## `frozen_magnon` `branch_mode = 'auto'`: multi-sublattice acoustic magnon not gapless at Γ — RE-MEASURED CLEAN on the real-space route (WP9, 2026-08-07), k-space route not re-checked

- **WP9 update (2026-08-07), real-space route:** re-measured on the current
  `fable_v2_gbt_v2` architecture (post WP1–WP8), two-sublattice bcc FeCo
  (`example/bulk/bccFeCo`, `nsp=3`, `recur='block'`, `lld=21`,
  `strux_backend='strux_lib'`), `branch_mode='auto'`, `mode='mft'`,
  `theta_probe=20°`, real-space recursion (`&self` does not set
  `use_kspace`, so this is the RS route, not the k-space route the original
  finding below was measured on). Reference SCF: two inequivalent moments
  1.976 μB (Fe) / 1.657 μB (Co), `Total RMS Diff = 7.8e-4` (not fully
  converged to `conv_thr`, but stable and non-oscillating once seeded from a
  seperately-converged starting potential — see
  `tests/regression/wp9_validation/multisublattice_goldstone/`). Measured
  `omega(Γ)`, acoustic branch: **`-2.15e-22 Ry`** — zero to numerical noise,
  not the ~0.28 Ry violation below. The eigenvector at Γ has both
  sublattices in phase (amplitudes 0.7375/0.6753, both phase −π), exactly
  the uniform-rotation pattern a true acoustic mode requires; the optical
  branch (finite gap, `omega(Γ)=7.04e-3 Ry`) has them out of phase, as
  expected. Small-q growth (`omega` at q3=0.02/0.05: `1.18e-5`/`7.87e-5` Ry)
  is smooth and roughly quadratic, not flat. This is real, measured evidence
  (ran the binary, read `frozen_magnon_branches.dat`/`_modes.dat` directly),
  not inferred — but it is **one system, one route, one theta_probe**, and
  does **not** by itself confirm the k-space route (where the original
  ~0.28 Ry number below was measured) is also fixed; that was not re-tested
  in this task. Do not treat this as a full resolution of the entry below
  without re-checking k-space. Plausible explanation for the improvement:
  WP3's explicit `theta_ss_sublattice`/`phi_ss_sublattice` plumbing (see
  `source/calculation.f90::post_processing_frozen_magnon_auto`,
  `set_reference_sublattice_angles`) looks materially different from, and
  more careful than, what finding K5 (referenced below) describes for the
  pre-WP1 code this entry was originally written against — read directly,
  not inferred from the old prose.
- **Original symptom (pre-WP1, k-space route, NOT re-verified against
  current HEAD):** the multi-sublattice magnon branches from
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
- **Status:** partially re-measured (see the WP9 update above — real-space
  route now looks clean on one system), not closed. Still open: re-run the
  same check on the k-space route (where the original violation was
  measured) before declaring this fixed; the multi-branch magnon spectrum
  remains the validation target for blueprint **B11** (linear-response
  TDDFT / transverse magnons) for the general (non-collinear-reference)
  case. See `docs/dev/B1_GBT_SPIN_SPIRAL_PLAN.md` (T5) and commit `d86fe42`.

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

## `calctype = 'L'` (111) site DOS deviates ~2e-3 from bulk; (001) is exact

- **Symptom:** in the cross-calctype fcc Cu oracle (see
  `tests/scf/README.md`, "Cross-calctype oracle"), a single Cu layer treated
  as an interface between Cu regions must reproduce bulk fcc Cu, because every
  region is the same material starting from the same parameter set with
  `vmad ~ 0`. `surftype = '0 0 1'` does so **exactly** — zero difference at
  every sampled DOS row. `surftype = '1 1 1'` deviates by **2.05e-3** at
  row 1200 (E = 0.686, near the d-band peak), against a peak DOS of 48.3.
  `Example_interface_fccCu001_chebyshev` vs
  `Example_interface_fccCu111_chebyshev`.
- **Not the cause:** the TB-LMTO Hamiltonian. Instrumenting `build_bulkham`
  and `build_locham` with a geometry-keyed dump (per-neighbour displacement
  vector plus the hopping block's Frobenius norm, matched across calctypes by
  vector rather than by neighbour index) shows the on-site block and all 19
  fcc neighbour hoppings **bit-identical** across `B`/`I`/`L`. `etot` also
  agrees to ~1e-8 Ry and `ws_r` exactly, for both orientations.
- **Therefore:** the residual is downstream of the Hamiltonian, and is
  specific to the 111 surface normal. Candidates not yet investigated: the
  layer ladder / `zstep` determination in `build_interface_full` for a
  non-cubic normal (`dx,dy,dz` are transformed for hcp/111 cases before the
  layer scan), the resulting selected-cluster boundary, or the
  representative-site choice interacting with 111's different in-plane
  periodicity.
- **Deliberately captured, not tolerated away:** the committed reference pins
  the current 111 values, so any change to this residual shows up as a test
  failure rather than passing silently. Do not widen `abs_tol`/`rel_tol` to
  make the two orientations agree.
- **Found:** B7.5 (`calctype='L'` wiring, commit `97f1e0e`). The earlier
  validation used 001 only and reported agreement at print precision; the 111
  deviation surfaced when both orientations were added to the suite.

## `calctype = 'L'`: raising `&charge nlay_a/nlay_b` breaks the alignment fixed point

- **Symptom:** in the `A | A` identity geometry (`example/interface/fccCu111_AA`,
  two frozen Cu regions around one active Cu layer, all from the same converged
  parameter set), the alignment solver must converge to `V(B) = 0`. Raising the
  **`&charge`** row counts breaks that; raising the **`&lattice`** layer counts
  does not. Single-variable runs, 5 iterations each:

  | `&lattice nlay_a/nlay_b` | `&charge nlay_a/nlay_b` | converged `V(B)` | identity |
  | --- | --- | --- | --- |
  | 1 / 1 | 1 / 1 | `0.000000` Ry | holds |
  | 4 / 4 | 1 / 1 | `0.000000` Ry | holds |
  | 1 / 1 | 4 / 4 | `-0.4498` Ry | broken |

- **The two knobs are different, and only one is the buffer width.**
  `&lattice nlay_a/nlay_b` count **atomic layers** and are the correct way to
  widen the frozen buffer — that is what `build_interface_full` bins sites into
  (`lattice_cluster.f90:576-624`), and widening it is harmless.
  `&charge nlay_a/nlay_b` count **rows of the synthetic 2D Madelung stack**,
  whose size is a fixed constant (`this%nbas = max(49, ...)`,
  `lattice_cluster.f90:709`) deliberately decoupled from the physical layer
  count so `set2d`'s NLAMA/NLAMB split stays balanced.
  `build_interface_registry` computes `nlay_active = nbas - nlay_a - nlay_b`
  over those rows (`charge.f90:1711`), so raising them relabels *interior*
  Madelung rows — rows carrying real deviation charge — as frozen boundary and
  moves the deep probe onto a charged row.
- **Consequence for users:** widen the frozen buffer with `&lattice`, and leave
  `&charge nlay_a/nlay_b` at 1 unless the Madelung row partition is genuinely
  what you mean to change. The source comments at `lattice.f90:320-331` already
  state that the two are deliberately separate (LAYERS vs SITES); the trap is
  that the names are identical.
- **Still open:** the `align_regions` "widen the active zone" drift warning
  fires even in the `&lattice`-widened runs where the identity holds exactly, so
  it is noisier than its wording implies. Whether the deep-probe drift threshold
  should be recalibrated is a B7.7 question (G-B7-3 revisits the tolerances).
- **Found:** B7.6 (examples and documentation). Documented for users in
  `docs/source/user_guide/examples/interface_fcccu111.rst`.

## `calctype = 'L'`: no vacuum region, and `vacuum_lead` has no caller

**RESOLVED for `A | vacuum`** (B7.6 wiring). `&lattice region_b_kind = 'vacuum'`
now makes region B a vacuum region: `build_from_interface` takes a `kind_b`
argument, `build_interface_registry` passes `region_kind_vacuum`, and
`refresh_vacuum_region` generates the frozen parameters per run and regenerates
them each iteration at the solved vacuum level. Example:
`example/interface/fccCu111_Avac`.

**Still open:** `A | vacuum-gap | B`. That needs a genuine four-region layout
(`lead_a | active-vacuum | active | lead_b`) and a rework of the type-block
arithmetic in `build_interface_full`, which still hardcodes three regions.

- **Original symptom:** two of the four geometries B7 §1.2 scopes for the
  layered path — `A | vacuum` and `A | vacuum-gap | B` — could not be
  constructed.
- **Original cause (no vacuum region):** `region_registry_build_from_interface`
  hardcoded kinds `lead_a`, `active`, `lead_b`; nothing on the `buildinterface`
  route could assign `region_kind_vacuum`.
- **Original cause (generator unwired):** `source/vacuum_lead.f90` was a tested
  component with no consumer — its only mention outside its own file and
  `CMakeLists.txt` was a comment in `self.f90:263`.
- **Found:** B7.6 (examples and documentation).

## `calctype = 'L'`: interface electrostatics were identically zero — **FIXED**

**RESOLVED.** Three separate index-space bugs made `Q`, `P`, `step` and every
deep probe come out exactly zero for every layered run. All three are fixed;
this entry is kept because the failure mode is instructive and because the
committed `tests/scf` interface references were generated while it was active.

1. **Active charge written to a frozen row.** `interfacepot` conflated two
   index spaces: `atomrec` (1..nrec, the active TYPE counter — indexes `dq`,
   `chargetrf_type`, `symbolic_atom(nbulk+.)`) and the Madelung ROW (1..nbas,
   which is what `tdq`/`tq10`/`vm` and every registry array are indexed by).
   The active zone starts at row `nlay_a+1`, so writing to `tdq(atomrec)` put
   the charge on region A's frozen boundary row. `compensation_sites` then
   returned `ideep_lo` = that same row and subtracted the whole residual from
   it — exact cancellation, `tdq ≡ 0`, `vm ≡ 0`. Fixed with an explicit
   `irow = nlay_a + atomrec` in both the charge loop and the write-back.
2. **`boundary_nef` used the same wrong mapping** (`nbulk + isite` with
   `isite` a Madelung row), which is out of range for a boundary row. It
   returned 0 for *both* sides, so the N(E_F) compensation weighting collapsed
   to the 50/50 fallback and **vacuum silently received half the compensation
   charge** — exactly what B7 §1.5 warns about ("compensation placed there
   does not perturb the work function, it SETS it"). Now reads the registry's
   per-site reference type.
3. **`reference_type` was filled by cycling active-type values across all
   rows**, so frozen boundary rows carried no valid type — which is what made
   (2) silent. Now assigned per region: region A's types on its rows, region
   B's (or the single vacuum type) on its rows, `chargetrf_type` on active rows.

**Also added: the active zone is now centred in the Madelung stack by default.**
The active charge occupies rows `nlay_a+1 .. nlay_a+nlay`, and the deep probes
are the extreme frozen rows either side. An off-centre split puts one probe
adjacent to the charge and the other ~`nbas` rows away, so the solver correctly
reports a nonzero `V_B` for a physically symmetric cell. Measured on
fccCu111_AA (nbas=49, one active layer): `&charge` 1/1 → `V(B) = -0.0109` Ry;
centred 24/24 → exactly 0. Leaving `&charge nlay_a/nlay_b` unset now derives
`(nbas - nlay)/2` either side and logs the choice. Explicit values still win.

**Why it hid:** the shipped oracle is the `A | A` identity, whose *correct*
answer is exactly zero for all five reported quantities — a spuriously zero
result is indistinguishable from a right one. B7 §5.3 anticipated precisely
this and specifies the real oracle as A-against-A-with-a-rigid-offset; the
shipped example did not honour that.

**Verification after the fix:**

| case | Q | P | step | V(B / vacuum) |
| --- | --- | --- | --- | --- |
| `A \| A` identity | 0 | 0 | 0 | 0 |
| `A \| vacuum` | 0 | -1.12e-2 | -0.0949 Ry | +0.0949 Ry |

The identity still holds exactly; the surface now produces a real dipole
barrier where it previously produced none. **Buffer-width convergence** (`&lattice nlay_a`/`nlay_b`, fcc Cu(111)):

| buffer | step |
| --- | --- |
| 1 / 1 | -0.0949 Ry |
| 2 / 2 | -0.0977 Ry |
| 4 / 4 | -0.0977 Ry |

Converged by width 2 and stable to five digits at width 4 — the barrier is a
real converged quantity, not a buffer artefact.

Cross-check against the independent one-sided `buildsurf` route on the same
system: `buildsurf` gives `vmad1 = vm1 - vbulk = -0.1236` Ry against the
interface route's converged `-0.0977` Ry — same sign and magnitude, ~21%
apart. The two probe different points of the same profile (`buildsurf` uses
rows 1 and `nbas`; the interface route uses the frozen-region extremes), so
exact agreement is not expected, but closing that gap is a B7.7 item.

**Consequence for the committed references:** `tests/scf` interface cases pass
unchanged because they pin DOS and moments, not electrostatic quantities. Any
future reference regeneration will capture genuinely different `vmad` values.

- **Found:** B7.6 (vacuum-lead wiring), while checking why the generated vacuum
  lead produced no measurable dipole barrier.


## `calctype = 'L'`: `A | vacuum` diverges with more than one active layer — UNTRIAGED

- **Symptom:** with `&lattice nlay = 3` (three active layers) on the
  `A | vacuum` geometry, the reported potential step and the active atoms'
  `vmad` come out physically impossible — hundreds of Rydberg:

  | case, `nlay = 3` | step | `vmad`, active layers 1 / 2 / 3 |
  | --- | --- | --- |
  | `A \| A` metallic | `0.0002` Ry | `-1.6e-4`, `-1.8e-3`, `-4e-15` |
  | `A \| vacuum` | **`-334` Ry** | **`167`, `56`, `0.38`** |

  The metallic three-layer case is sane, so whatever this is, it involves the
  vacuum region specifically. At `nlay = 1` both geometries are well behaved
  (`A | vacuum` gives a converged `step = -0.0977` Ry).
- **NOT YET DIAGNOSED, and the test deck is a live suspect.** The three-layer
  decks above were hand-built by editing `ntype`, `ct(:)` and the `&atoms`
  `label(:)` list of the one-layer examples and copying `Ac1.nml` to `Ac2.nml`
  / `Ac3.nml`. That is exactly the kind of setup that can be silently wrong —
  in particular the active-layer type block and `chargetrf_type` assignment
  (`build_interface_full`) were not independently verified for `nlay > 1` on
  the vacuum path. **Do not treat this as a confirmed code bug until a
  known-good multi-layer deck reproduces it.**
- **If it is real, the likely area** is the compensation weighting across
  multiple active rows: with one active row the compensation sites sit
  symmetrically either side of it, and with three they do not. That is
  precisely the machinery G-B7-3 exists to sign off, so a genuine finding here
  belongs to B7.7 rather than to a spot fix.
- **First triage steps** (in order): build a multi-layer deck from
  `buildsurf`'s own working three-layer surface case rather than by hand;
  confirm `chargetrf_type` and the layer→type map for `nlay > 1`; then dump
  `tdq` per row to see whether the divergence enters through the charge, the
  compensation, or the kernel.
- **Found:** B7.6 (vacuum-lead wiring), while probing whether the new
  electrostatics path had multi-layer coverage. Not investigated further —
  out of scope for the wiring task.

## Magnetic constraining field (`&control constraints_enable`) is a no-op

- **Symptom:** none visible to the user — the mechanism silently does
  nothing. `constraints_enable = .true.` (with `constraints_i_cons`,
  `constraints_mom_ref`, `constraints_bfield`, etc.) runs without error and
  produces ordinary unconstrained SCF output.
- **Root cause, verified by reading the code:** both call sites of
  `cfd::constrain` (`source/self.f90:1082-1102`, inside `run_dos`'s real-space
  SCF branch, and `source/exchange.f90:1475-1497`, inside
  `calculate_exchange`) allocate local `mom_in`/`mom_ref`/`bfield`, call
  `constrain(mom_in, mom_ref, bfield, nrec)` (`source/include_codes/abspinlib/constrain.f90:60-234`),
  then immediately `deallocate` all three without ever writing `bfield` back
  into `potential%mom`, any Hamiltonian array (`ee`/`hall`), or
  `this%control%constraints_bfield`. The constraint energy `etcon` computed
  inside `constrain` (e.g. `constrain.f90:89,96,149`) is a plain local
  `real(dblprec)`, never returned, never logged, and never added to the total
  energy. A user-supplied `constraints_bfield` seed is also discarded: both
  call sites zero-initialize `bfield(:, ia_loc) = 0.0_dblprec` immediately
  before calling `constrain`, so any namelist value never reaches it.
- **Scope:** this is independent of the Generalized Bloch Theorem work —
  the mechanism is equally inert under `periodic_nc`, `explicit_texture`, and
  `gbt_single_q`. There is currently no Hamiltonian/energy term for any
  representation mode to gauge or gate.
- **Also noted while reading `constrain.f90`, unrelated to the above:** for
  `constraints_i_cons` 2 and 3, the local accumulator `etcon` is read
  (`etcon = etcon + ...`) inside the `do na` loop (`constrain.f90:89,96`)
  with no `etcon = 0.0_dblprec` initialization beforehand — a latent
  uninitialized-read, though moot for now since `etcon` is discarded anyway.
- **Found:** WP6c (GBT audit of Hubbard/constraints/velocity/torque/SOC
  terms), while checking what frame the constraining field uses. Not fixed
  here — wiring it up is a real feature-completion task (decide the target
  frame, thread `bfield` into the potential/Hamiltonian, add the energy term,
  fix the `constraints_bfield` seed discard, fix the `etcon` initialization),
  estimated at least one dedicated task, not a WP6c-sized fix.

## Five `tests/scf/cases.json` GBT fixtures fatal on `strux_backend='legacy'` (pre-existing, not a WP6c regression)

- **Symptom:** `Example_bulk_bccFe_nsp4_block_spiral_qplus`,
  `Example_bulk_bccFe_nsp4_block_spiral_qminus`, `Example_frozen_magnon_bccFe`,
  `Example_frozen_magnon_bccFe_auto`, and `Example_frozen_magnon_bccFe_auto_scf`
  all abort with `gbt_single_q requires strux_backend='strux_lib'; the legacy
  backend is unsupported.` (`hamiltonian_build.f90::build_gbt_bulkham`)
  instead of producing output to compare against their reference.
- **Root cause, verified:** `calculation.f90:1836` and `:1977`
  unconditionally set `hamiltonian_obj%magnetic_representation = gbt_single_q`
  for the bulk-spiral (`q_ss`/`theta_ss` active) and `frozen_magnon`
  post-processing workflows. `lattice%strux_backend` defaults to `'legacy'`
  (`lattice_lifecycle.f90:404,654,920`) unless a case's namelist explicitly
  sets `strux_backend='strux_lib'`. These five `tests/scf/cases.json` entries
  do not set it, so every run through these workflows now hits the WP4
  legacy-backend guard (`hamiltonian_build.f90`) before producing any
  Hamiltonian.
- **Confirmed pre-existing, not introduced by WP6c:** reproduces identically,
  same fatal, same line, on the pre-WP6c commit `0188a6a` with no source
  changes (checked via `git stash` + rebuild before diagnosing this). WP6a
  and WP6b's own gate evidence only reports `ctest -L unit` pass counts
  (18/18, 19/19) — the SCF example suite (`ctest -L scf`) does not appear to
  have been run as part of closing those gates, so this breakage was not
  caught then either. It predates WP6c and is unrelated to the
  Hubbard-U/V/constraints/velocity/torque/SOC terms WP6c actually audited.
- **Two more `tests/scf/cases.json` failures are separate and older still:**
  `Example_bulk_bccFe_nsp2_block`/`_hoh` fail a single `totaldos.out` row by
  ~4-6e-5 (rel ~2e-6) — this matches the pre-existing gfortran-13 DOS
  tolerance delta already recorded in `docs/dev/GBT_WP0_G0_REPORT.md`, not a
  new issue.
- **Proposed fix (not applied — regression-reference-touching, needs a
  green light):** add `"lattice": {"strux_backend": "strux_lib"}` (and
  whatever `screening`/`strux_want_sdot` the GBT path additionally requires)
  to the five affected `tests/scf/cases.json` entries, then regenerate their
  reference values on the corrected `strux_lib` GBT path — the numbers will
  differ from the old `legacy`-backend GBT output, so this is a reference
  regeneration, not a one-line config patch. Estimated size: small
  (config + one `tests/generate_references.py` run + review), but touches
  golden values and so needs explicit sign-off before doing it.
- **Found:** WP6c, while running the full `ctest -L scf` suite as required by
  the project's rule 1 regression gate before/after every task.

## `nsp4_block_spiral_qplus`/`_qminus` never enable GBT (diagnosed)

- **Symptom:** both cases fail on `Fe_out.nml:mom[3]`, `run = 1.000000e+00`
  against `ref = 6.856242e-09`. The converged moment is `+z` where the
  reference has it in-plane, as a `theta_ss = 90` cone should be.
- **Cause (diagnosed 2026-08-06, WP7):** the two `tests/scf/cases.json`
  entries set `hamiltonian.q_ss = (0,0,±0.05)` and `theta_ss = 90.0` but
  **never set `magnetic_representation`**. Since WP3/WP5 made the
  representation explicit, it defaults to `periodic_nc`, under which `q_ss` is
  stored and never read — so both cases have been converging an ordinary
  collinear ferromagnet along `+z`. `mom[3] = 1.0` is exactly that. There is
  no spiral in these runs and there has not been one since the representation
  split landed. They also set `nsp = 4`, i.e. SOC on, which GBT rejects
  outright.
- **How it surfaced:** the WP7 input guard
  (`validate_spiral_keys_are_consumed`, `hamiltonian_build.f90`) now makes this
  fatal, so the failure mode changed from a silent wrong value to an explicit
  rejection naming the missing key. The committed reference `mom[3] ≈ 6.9e-9`
  presumably dates from when the spiral was enabled implicitly (via
  `gbt_kspace` or the absolute-position branch), both since deleted.
- **Proposed fix (not applied — touches golden references, needs a green
  light):** add `"magnetic_representation": "gbt_single_q"` to both cases'
  `hamiltonian` namelist, change `nsp` from `4` to `3` (SOC is fatal under
  GBT), add `"lattice": {"strux_backend": "strux_lib"}` (the WP4 requirement),
  then regenerate both references. The values *will* move — these cases have
  not computed a spiral in a long time — so this is a reference regeneration
  with a physics justification, not a config patch. It would also restore the
  `E(q) = E(-q)` pin the pair was written to provide, which is currently
  vacuous (both sides compute the same ferromagnet).
- **Found:** WP7, answering "does GBT functionally work?".

## [RESOLVED 2026-08-07, test-fixture bug] `magnetic_representation = 'gbt_single_q'` decks must set `mom` to the collinear-rotating-frame convention `(0,0,±1)`, never to the physical cone direction — WP9's commensurate-supercell decks got this backwards

- **Root cause, fully traced (2026-08-07):** in `gbt_single_q` (for a
  `hoh=.false.` deck — the WP9 decks below never call `build_obarm`/
  `build_enim` at all, confirmed empirically, see below), the bond Hamiltonian
  (`gbt_contract_collinear`, driven by `gbt_endpoint_angles`) and the on-site
  correction are built **entirely** from `theta_ss`/`phi_ss_sublattice` plus
  the *sign* of `potential%mom(3)` (a binary up/down-sublattice selector,
  `source/hamiltonian_build.f90:1619`) — `potential%mom(1)`/`mom(2)` and the
  magnitude of `mom(3)` never enter the Hamiltonian construction at all.
  *But* `ql`'s up/down decomposition — what every moment readout (including
  the mixing/SCF-convergence machinery) is built on — separately **projects
  the density onto `potential%mom` as its quantization axis**. So for
  `gbt_single_q`, `mom` must always be set in the code's internal
  collinear-rotating-frame convention (nominally `(0,0,1)`, i.e. "+z always
  means the rotating-frame reference," exactly what the already-validated
  `tests/regression/wp8_littlegroup/base/Fe.nml` does even at
  `theta_ss=90`) — **not** pre-rotated to the physical lab-frame cone
  direction `m₀=(sinθ,0,cosθ)`. The WP9 decks below (ported from the older,
  pre-representation-split `3fd21c0` fixture, where that convention may have
  been different) set `mom=(1,0,0)` for `theta_ss=90`, i.e. physically
  pre-rotated — wrong for the current architecture.
- **Direct proof:** with the wrong `mom=(1,0,0)`, the physical density is
  unaffected (`source/bands.f90`'s `DENSITY_POLICY` diagnostic, which is
  independent of `mom`, reads `m_long=0.000000 |m_transverse|=2.395328` —
  the moment is present, magnitude ≈2.395 μB, matching the supercell almost
  exactly) but `ql`'s `mom`-projected up/down split reads that as ≈0 (a
  z-pointing vector projected onto x). Editing only `mom` from `(1,0,0)` to
  `(0,0,1)` on an otherwise byte-identical deck reproduces the **exact same**
  `DENSITY_POLICY` physics (confirming the Hamiltonian truly doesn't depend
  on `mom`'s value) but now `ql` correctly reads ≈2.395 μB, matching
  `DENSITY_POLICY`'s `m_long`.
- **Fix applied:** `mom(:) = 0.0d0, 0.0d0, 1.0d0` in all six `gbt_single_q`
  decks (`gbt_supercell/{q050,q033}/gbt{,_scf,_constrained}/Fe.nml`). Effect,
  re-running the registered `WP9CommensurateSupercell_*` ctests:

  | case | moment gap before | moment gap after | eband gap (unaffected — frame-invariant) |
  | --- | --- | --- | --- |
  | q050 MFT | 2.396 μB | **5.78e-4 μB** (~4100x better) | 7.19e-4 Ry (unchanged) |
  | q033 MFT | ~2.33 μB | **6.38e-3 μB** (~360x better) | 5.98e-3 Ry (unchanged) |
  | q050 SCF | already ~1.3e-3 μB | unchanged (already correctly axis-aligned once mixing converges) | 2.73e-4 Ry (unchanged, **now the sole remaining failure**) |

  All four `WP9CommensurateSupercell_*` cases still **FAIL** their derived
  tolerances — but now on the band-energy residual alone, not the moment
  direction. That residual is real, smaller, and separate; see below.
- **What remains open (the actual physics question, not a fixture bug):**
  a band-energy gap of 7.2e-4–6.0e-3 Ry/atom (MFT leg, larger for the
  6-atom q033 case than the 4-atom q050 case) and 2.7-2.8e-4 Ry/atom (SCF
  leg, q050 only) between the `gbt_single_q` single-atom result and the
  explicit lab-frame supercell, exceeding tolerances derived with 6x+
  headroom over the measured noise floor. Leading hypothesis, not confirmed:
  the MFT leg's frozen starting potential is carried over from the older
  `3fd21c0` kernel and may not be self-consistent under the current one
  (the SCF leg, which lets the same starting potential relax, shows a
  smaller gap, consistent with this). Not traced to a specific source line.
- **Time-reversal checked and ruled out as a contributor to this specific
  finding** (WP9 integrator, 2026-08-07, prompted by a question about
  local/global axis handling): `force_full_bz_for_nonzero_q_gbt`
  (`source/reciprocal_lifecycle.f90:526`) unconditionally forces both
  `use_symmetry_reduction` and `use_time_reversal` to `.false.` for any
  nonzero-q `gbt_single_q` k-space build, regardless of the namelist;
  confirmed by log inspection (`"nonzero-q GBT rebuilding the full chemical
  BZ mesh"`, full unreduced mesh point count). Not relevant to this
  real-space-route battery in the first place, but checked because the same
  question applies to Battery B (see the Gamma-H sweep entry below).
- **What I verified vs. guessed:** the `mom`/`ql`/`DENSITY_POLICY` mechanism
  above is directly verified by reading `gbt_endpoint_angles`,
  `gbt_contract_collinear`'s call site in `build_gbt_bulkham`, and the `hoh`
  default (`source/hamiltonian_build.f90:559`, confirms `build_obarm`/
  `build_enim` are dead code for these `hoh`-unset decks), plus the
  before/after `mom` experiment run directly. The remaining band-energy
  residual's cause (stale starting potential) is still a guess, not traced
  further.
- **Original symptom, verified directly (build_13 binary, current `fable_v2_gbt_v2`
  HEAD, bcc Fe, `alat=2.8612`, `nsp=3`, `recur='block'`, `lld=16`,
  `strux_backend='strux_lib'`, single-atom cell, `magnetic_representation=
  'gbt_single_q'`, starting potential carried over unmodified from the WP9
  `_b1r5_reference` decks (commit `3fd21c0`, a materially different, older
  GBT kernel)):
  - **Frozen-potential evaluation (`nstep=1`, `beta=1.0`, `magbeta=0`, the
    force-theorem convention):** at `q_ss=(0,0,0.5)`, `theta_ss=90°`, the
    output `Fe_out.nml` has `ql(1,:,1)` within ~1e-8 of `ql(1,:,2)` (up ≈
    down), i.e. the moment reads as ≈0 rather than the ≈2.4 μB the paired
    explicit supercell converges to (see the WP9 battery report for the
    supercell numbers). Band energy differs by only ~7e-4 Ry/atom, i.e. the
    scalar/energy channel is close but the magnetic channel reads as absent.
  - **At `q_ss=0, theta_ss=0`** (a literal collinear FM, which per
    `docs/dev/plans/B1_gbt_frozen_magnons_v2.md` §2.4 "must remain
    bit-identical to today's collinear/noncollinear-FM output"), the same
    frozen-potential, `nstep=1` evaluation gives `ql(1,:,1)`
    **bit-identical** to `ql(1,:,2)` under `gbt_single_q`, versus a real
    ≈1.209 μB moment under `magnetic_representation='periodic_nc'` on the
    otherwise-identical deck (band energy differs by only ~1.3e-4 Ry:
    `-2.3429851382` gbt vs `-2.3428537567` periodic_nc).
  - **Full SCF from the same starting potential** (`nstep=25`,
    `mixtype='broyden'`, `beta=0.15`, `magbeta` at its default of 1.0) at
    `q_ss=(0,0,0.5)`, `theta_ss=90°` converges smoothly (monotone decreasing
    mixing residual, no charge-conservation warnings) to a **real** moment of
    ≈2.253 μB and `eband=-1.9926859349` Ry — close to, but not matching, the
    supercell's ≈2.4 μB. This is a much smaller gap than the frozen-potential
    evaluation shows, which points at the *stale, foreign-kernel starting
    potential* rather than the Hamiltonian construction itself as (at least
    part of) the explanation for the `nstep=1` near-zero reading.
  - **The same full-SCF recipe at `q_ss=0, theta_ss=0`, by contrast,
    diverges**: `mix.f90` logs "too much charge in the external atom!",
    charge transfer of `-8.0`, and the run ends at a degenerate state
    (`Band energy of system: 0.0000000000`, moment "considered induced").
    This was **not** cross-checked against `periodic_nc` under the same
    full-SCF recipe, so whether this specific divergence is `gbt_single_q`-
    specific or a starting-point/mixing-schedule artifact common to both
    representations is unknown.
- **What I verified vs. what I am guessing:** all six numbers above are
  measured (ran the binary, read `report.out`/`Fe_out.nml`/the run log), not
  inferred. What I am guessing: that the frozen `3fd21c0` potential is
  simply not self-consistent under the current kernel's contraction (a
  starting-point mismatch) rather than the moment-zeroing being a permanent
  property of `gbt_single_q` regardless of input — the SCF-recovers-most-
  of-it observation supports this reading but I did not trace it to a
  specific line in `source/hamiltonian_build.f90`/`source/gbt_structure.f90`
  (`build_gbt_bulkham`, `gbt_contract_collinear`, `gbt_endpoint_angles`) or
  rule out `self.f90`'s potential-parameter orthogonalization step
  (off-limits per repo rule 5) as a contributor. I also did not diagnose the
  `q=0` SCF divergence at all — it surfaced from an extra side-experiment
  beyond the assigned task and is reported here only because it is a real,
  reproducible symptom, not because I understand it.
- **Impact (superseded by the fix above, kept for the historical trail):**
  before the `mom` fix, the MFT leg looked like it was destroying the
  moment entirely; it was not — the fixture was reading the density's
  projection onto the wrong axis. The `q=0, theta_ss=0` full-SCF divergence
  noted above was **not reproduced** by a later, independent q=0 control
  (see the follow-up note two paragraphs up in git history / superseded
  text) using a fresh (non-carried-over) starting potential, which converged
  cleanly and matched `periodic_nc` — pointing at the foreign starting
  potential, not `gbt_single_q` itself, consistent with the current
  understanding above.
- **Remaining follow-up (not applied, unrelated to the fix above):**
  diagnose the band-energy residual described in the "what remains open"
  paragraph above — likely start by regenerating the MFT leg's frozen
  starting potential from an SCF convergence under the *current* kernel
  rather than reusing the `3fd21c0` state, to test the leading hypothesis
  directly. Estimated size: small, mostly rerunning with a modified fixture.
- **Found:** WP9 Battery A (commensurate-supercell known-answer test), while
  porting the `3fd21c0` supercell decks to the current architecture; root
  cause traced by the WP9 integrator the same day after a targeted question
  about local-vs-global axis handling.

## `frozen_magnon` k-space cone-angle (theta_ss) scaling is not theta-independent, by a large margin

- **Symptom, verified directly** (bcc Fe, WP9 Battery B,
  `tests/regression/wp9_validation/gammaH_sweep/base_cone/`, k-space route,
  `nk=12`, fixed small `q_ss=(0,0,0.05)`, `theta_ss` swept 5/10/15/20
  degrees, `omega(q) = 4[E(q)-E(0)]/(M sin^2 theta)`): `omega` measures
  `-6.417e-3, -1.293e-3, -3.426e-4, -8.270e-6` Ry respectively — a **~776x**
  spread across the window, decreasing monotonically as theta grows, rather
  than the flat line the harmonic-regime self-diagnostic
  (`docs/dev/plans/B1_gbt_frozen_magnons_v2.md` §2.10) predicts.
- **Time-reversal/BZ-reduction checked and ruled out:** confirmed by direct
  log inspection that `force_full_bz_for_nonzero_q_gbt` correctly forces the
  full, unreduced k-mesh (no time-reversal, no spatial symmetry reduction)
  for every one of these nonzero-q builds, regardless of the deck's
  `use_symmetry_reduction`/`use_time_reversal` namelist settings — so an
  under-reduced or over-reduced BZ integral is not the explanation.
- **Not diagnosed further:** a naive noise-amplification argument (fixed
  absolute noise in the band-energy difference, divided by `sin^2(theta)`)
  predicts only a ~15x range across 5-20 degrees, an order of magnitude
  short of the measured 776x — so simple additive Fermi-search/DOS noise
  does not by itself explain the spread. Whether the noise itself grows at
  small theta (a smaller induced-moment signal proportionally noisier
  against the same DOS/mesh discretization) or something else is at play is
  unknown.
- **Impact:** the code cannot currently be trusted to report a stable
  small-cone spin-stiffness estimate — the answer depends heavily, and
  unphysically, on which small test angle is chosen.
- **Found:** WP9 Battery B (bcc-Fe Gamma-H frozen-magnon sweep).

## `frozen_magnon` k-space EBAND mesh refinement is mildly non-monotonic at some q

- **Symptom, verified directly** (same deck as above, `q_ss=(0,0,0.5)`,
  `nk=8/12/16`): `eband` = `-1.98668388, -1.98714329, -1.98785909` Ry — the
  refinement step **grows** from `nk=8->12` (4.59e-4 Ry) to `nk=12->16`
  (7.16e-4 Ry) instead of shrinking. At `q_ss=(0,0,1.0)` (H) the same sweep
  is well-behaved (step shrinks `5.72e-4 -> 5.25e-4` Ry). Both steps are the
  same order of magnitude and well inside the overall spread
  (1.175e-3 Ry) — not a divergence, just non-monotonic.
- **Not diagnosed:** plausibly ordinary tetrahedron-integration
  discretization noise at that particular q-point rather than anything
  systematic, but not investigated.
- **Found:** WP9 Battery B, while reworking the battery to a k-space-only
  design (RS route dropped per Anders: long-wavelength spirals can have real
  numerical problems from real-space cluster truncation/PBC artifacts).

## `strux_backend='strux_lib'` badly breaks a `crystal_sym='file'` custom-lattice supercell that `legacy` handles cleanly

- **Symptom, verified directly** (same bcc-Fe 4-atom `q_ss=(0,0,0.5)`
  explicit supercell deck as the entry above, `periodic_nc`, `crystal_sym=
  'file'`, hand-built `lattice.nml`, `nstep=1` force theorem): with the
  default `strux_backend='legacy'`, all four (translationally-equivalent-
  by-construction) sites converge to charge `8.0000000-8.0000010` e and
  moment `2.3959047-2.3959059` mu_B (agreeing to ~1e-6, as expected for
  equivalent sites). Adding `strux_backend='strux_lib'` to the same deck
  (otherwise byte-identical) gives charge `7.972040` e (site 1) vs
  `8.347495` e (site 2) — charge not even close to conserved per site,
  breaking translational equivalence outright — and moment `2.112684` vs
  `2.452723` mu_B, a ~15% site-to-site spread. Band energy also moves by a
  large amount (`eband_total = -7.1317212193` vs `-7.9268686408` legacy,
  ~0.20 Ry/atom).
- **What I verified vs. guessed:** ran both decks directly and read
  `report.out`/`Fe*_out.nml`; the numbers are measured. I did not
  investigate why `strux_lib` fails specifically on a `crystal_sym='file'`
  hand-built lattice with 4 distinct types (vs. `crystal_sym='bcc'` +
  `ntype=1`, where `strux_lib` is required and works fine for the
  `gbt_single_q` side of the same battery) — that is unexplored.
- **Impact:** confirms `legacy` (the default) is the right backend choice
  for the WP9 commensurate-supercell battery's `super/`-side `periodic_nc`
  decks; `strux_lib` was tried per this task's instruction to check
  empirically and is unusable there. Flagging per this task's own
  instruction to report backend disagreement beyond noise rather than
  paper over it — this is *far* beyond noise.
- **Found:** WP9 Battery A, checking empirically (as instructed) whether
  the `super/` side benefits from `strux_backend='strux_lib'` for a clean
  backend-matched comparison against the `gbt/` side.
