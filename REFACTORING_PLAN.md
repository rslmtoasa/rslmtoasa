# RS-LMTO-ASA (`fable_v1`) — Refactoring Instructions

**Audience:** coding assistant (Claude Code / Codex) working on this repo.
**Goal:** reduce structural redundancy accumulated during recent AI-assisted feature work, so that new features can be added *into* the architecture instead of *beside* it. Physics and numerics must not change.
**Scope:** structure, modularity, dead code, and performance hygiene. No new features.

---

## 0. Ground rules (read first, apply to every task)

1. **Bit-level behavior is the contract.** Every task below must leave
   `tests/regression` (bccFe_block, bccFe_chebyshev, bccFe_lanczos) and
   `tests/scf` passing at the same tolerances as before the change. Run the
   suite **before** starting (record baseline) and **after every task**.
2. **One task = one commit/PR.** Do not batch tasks. Each commit message
   references the task ID below (e.g. `refactor(T3): unify green eta variants`).
3. **Before deleting any routine**, verify it is uncalled by (a) `grep` across
   `*.f90 *.F90 *.cu *.cpp`, including type-bound `procedure ::` bindings and
   generic interfaces, and (b) a clean build with
   `-Wall -Wextra -Wunused` (gfortran) showing no new warnings. Move deleted
   code to the commit history, not to comment blocks.
4. **Do not change public namelist keys** (`cheb_backend`, `gpu_backend`, etc.).
   Backend names in the input file remain valid; only their internal routing
   changes.
5. **Do not touch numerics inside kernels** (loop order, precision of
   accumulators, fp32/fp64 boundaries, cgemm call shapes) unless the task
   explicitly says so. Cosmetic hoisting is fine; reassociation is not.
6. When a task says "extract shared helper", the helper goes in the module the
   task names — do not create new single-routine modules ad hoc.
7. If a task turns out to conflict with something you discover in the code,
   **stop and report** rather than improvising a different design.
8. **KISS.** This is a research code maintained by computational scientists,
   not a software-engineering showcase. Prefer lean, readable, obviously
   correct code over defensive machinery. Concretely:
   - Do not add argument-validation ladders, status/error return codes,
     wrapper layers, getters/setters, or "just in case" guards to internal
     routines. A physics kernel given a wrong input should fail loudly during
     testing, not carry permanent runtime checks.
   - Where recently added (AI-generated) routines already carry excessive
     validation, logging, or paranoid `if (allocated(...))` scaffolding on
     hot paths, **trim it** as part of whichever task touches that code —
     keep checks only at true boundaries (namelist parsing, file I/O,
     plugin/library availability).
   - One clear implementation beats a configurable abstraction. Only
     introduce an abstraction where this plan explicitly calls for it (e.g.
     the operator-apply interface in T4).
   - Fewer lines doing the same thing correctly is always the preferred
     direction.
9. **Preserve the class-based architecture.** The modernized code deliberately
   uses derived types with constructors, `destructor` finalization,
   `restore_to_default`, and type-bound procedures. This design stays: new
   helpers and refactored code follow the same pattern, and the submodule
   split (T9) must not flatten types into bare module procedures. Splitting
   the *implementation* across submodules while the type definition,
   constructor interface, and procedure bindings remain in the parent module
   is exactly what Fortran submodules are for.
10. **Legacy LMTO-ASA routines in `self.f90` and `symbolic_atom.f90` are
    off-limits at this stage.** They contain physics-dense legacy code
    scheduled for a separate, later audit. The *files* may be edited
    (updating a `use` statement, removing commented-out code around them,
    call-site adjustments required by other tasks), but the bodies and
    numerics of the legacy routines themselves must not be refactored,
    "cleaned", or restructured — unless a change is trivially safe, fully
    covered by existing tests, and explicitly reported as such.

---

## Phase 1 — Safety net and dead weight (low risk, do first)

### T1. Extend regression coverage to every backend being merged later
Currently the regression tests exercise only default paths. Before touching
kernels, add regression cases (smallest possible bccFe cell, few recursion
steps) that pin down:
- `cheb_backend = fast | batched | mkl_batch | mkl_sparse | legacy`
  (mkl variants conditionally, behind `ENABLE_MKL_*` detection in the runner)
- `hoh = .true.` with `cheb_backend = fast` and `legacy`
- `ccor_2c = .true.` path
- block Lanczos sp vs dp (`fast` vs `fast_dp`)

Reuse `tests/regression/test_regression.py` machinery; reference data comes
from the current `fable_v1` build (it is the trusted state). This task creates
the safety net that authorizes everything below.

**Also expand the fast example-based suite** behind the existing
`RUN_EXAMPLE_TESTS` CMake option so that, together with the regression cases,
every user-facing feature has at least one covered exercise path: the
`example/` families (bulk, surface, impurity, band_structure,
density_of_states, exchange, conductivity, k_space_scf) plus the recursion
backends and hoh/ccor variants above. Use minimal cells and few
recursion/energy points — these are smoke-and-compare tests, not production
runs; total added suite runtime should stay in minutes.

**Important — where testing lives (KISS):** all of this coverage belongs in
`tests/` and `example/`-driven CI, *not* inside the physics sources. Do not
add assertion scaffolding, self-check branches, or verification code paths to
the Fortran modules as part of this task. External tests pin behavior; the
sources stay lean.

### T2. Dead code removal
- `source/recursion.f90`: delete `chebyshev_recur_full` and
  `chebyshev_recur_s_ll` (already commented out of the type-bound procedure
  list as "obsolete") and their bodies.
- **Keep `bpopt` and `emami`** — they are live parts of the recursion
  terminator estimator chain (Green's function routines → `get_terminf` →
  `get_cinf` → `bpopt` → `emami`), invoked via type-bound calls. Add a brief
  header comment on each documenting this call chain so a future cleanup
  pass does not misidentify them as dead (a plain grep for `call bpopt`
  misses the `this%bpopt(...)` invocations).
- Remove commented-out executable code: `recursion.f90` has ~220 such lines
  (e.g. the `!!! call this%crecal_b()` blocks), `lattice.f90` and
  `hamiltonian.f90` ~50 each, `self.f90` ~30. Keep explanatory comments;
  delete only commented-out *code*.
- `source/cuda/kpm_reference.py` and `rsrec_validate.py`: if these are
  validation scripts worth keeping, move them to `tests/`; otherwise delete.
- Decide `cheb_moments_legacy`'s fate explicitly: it is still the mandatory
  fallback for `hoh` + BSR backends and `ccor_2c+hoh`, so **keep it**, but add
  a header comment stating exactly which feature combinations require it, so
  it is not "cleaned up" accidentally later.

### T3. Repo hygiene
- Fix the stray tab-indented lines in `hamiltonian.f90`
  (`build_ccor_pair_surface_block`, `ccor_lambda_components`,
  `ccor_spin_product`, `ccor_apply_spin_spiral`) — normalize to the 3-space
  indentation used everywhere else. Add/verify an `.editorconfig` or fprettify
  config so future tool-generated code conforms.
- `source/chebyshev_mkl_batch.cpp` (47 lines) — if it only wraps
  `cblas_cgemm_batch`, fold the wrapper into the existing MKL interface used
  by `chebyshev_fast.f90` via `iso_c_binding` and delete the file, adjusting
  `CMakeLists.txt`.

---

## Phase 2 — Kill the kernel clone families (main event)

### T4. Unify the Chebyshev CPU kernel family in `chebyshev_fast.f90`
Current state: `cheb_moments_fast`, `_fast_batched`, `_fast_mkl_batch`,
`_fast_mkl_sparse`, `_stochastic_fast`, `_orbital_fast` each carry private
nested copies of `apply_step*`, `cherk_full*`, `swapm*`, and
`hsweep/happly` — ~5 near-identical clones of each.

Target design:
1. Create one internal abstraction: an **operator-apply interface**
   ```fortran
   abstract interface
      subroutine ham_apply_t(x1, x0, y, alpha, beta)
   ```
   with four concrete implementations (dense-neighbor `apply_step`, BSR
   `apply_step_bsr`, MKL-batch, MKL-sparse) selected once per call by
   `cheb_backend`. The hoh two-sweep (`apply_step_hoh`) remains a separate
   implementation of the same interface (it already only exists for `fast`).
2. `cherk_full`, `swapm`, and the Chebyshev three-term loop become **single**
   shared routines used by all backends. The only per-backend code left is the
   operator apply and its buffer preparation.
3. `cheb_moments_stochastic_fast` and `cheb_moments_orbital_fast` reuse the
   same operator-apply objects instead of their private `hsweep/happly`
   copies; their outer moment-accumulation loops stay as they are.
4. Preserve the exact OpenMP scheduling clauses of each apply variant
   (`schedule(dynamic, 32)` vs `static`) — they are performance-tuned.

Validation: T1 regression cases for all four backends, plus a timing check
(each backend within noise of its pre-refactor wall time on the bccFe cases).

### T5. Replace the module-global cache state with a cache handle type
`chebyshev_fast.f90` lines ~35–60 hold ~25 `save`'d arrays with hand-rolled
validity flags (`scaled_cache_valid`, `bsr_cache_valid`, `vel_cache_*`, ...).
Create a derived type `cheb_cache_t` grouping: scaled Hamiltonian blocks
(hee/hha), ortho blocks, BSR arrays, velocity blocks, work buffers, MKL
pointer arrays, plus one `fingerprint` component (nb, nnmax, ntype, nmax, kk,
hoh, precision) that replaces the scattered `*_cache_nb` integers.
- One module-level instance of `cheb_cache_t` preserves current single-context
  behavior; `cheb_fast_reset_cache()` becomes a type-bound `reset`.
- The `ensure_*` routines become type-bound procedures on the cache and gain a
  single shared "fingerprint changed → rebuild" check instead of five
  duplicated comparison ladders.
This is pure restructuring — allocation sizes and rebuild triggers must be
identical.

### T6. Resolve operator selection once in the dispatcher
`recursion.f90 :: cheb_moments_cpu` repeats the same
`if (ccor_2c) then <ee_ccor_work,...> else <hamiltonian%ee,...>` if/else
inside **every** backend case. Refactor:
```fortran
ee  => merge_ptr(this%ee_ccor_work,  this%hamiltonian%ee,  use_ccor)
hall=> merge_ptr(this%hall_ccor_work, this%hamiltonian%hall, use_ccor)
```
(pointer or associate-based; arrays must not be copied), then a single call
per backend case. The select-case shrinks from ~100 lines to ~30. Apply the
same pattern to the equivalent triplicated `block_lanczos_fast` call sites at
recursion.f90 ~1845–1859 and ~2080–2094.

### T7. Deduplicate `green.f90`'s ×4 clone families
Three families each exist in up to four variants:
`{block_green, chebyshev_green, calculate_intersite_gf} × {plain, _eta} × {cpu, _gpu}`.
Target: **one routine per family** with signature
`(this, eta, fermi_point, g_ef)` where `eta` and `fermi_point/g_ef` are
`optional` — the plain variant is the call without them (this is exactly the
relationship between e.g. `block_green_ij` and `block_green_ij_eta` today).
GPU routing happens inside via the single early
`if (this%control%gpu_plugin) ...` branch that already exists in
`block_green`. Keep thin wrappers with the old names only where they are
type-bound procedures referenced from other modules, and mark them
deprecated in comments.
Constraint: the `_eta` variants integrate to the Fermi point with a finite
imaginary part — the merged routine must reproduce both numerical paths
exactly; do not "simplify" the eta=0 limit.

### T8. One GPU context layer
Two parallel Fortran wrappers exist over the same `librsrec` C API:
- `rsrec_cuda_plugin.f90` (`rsrec_cuda_backend`) — used by `recursion.f90`
- `cuda/recursion_gpu_mod.f90` (`rsgpu`) — used by `green.f90`

Additionally `green.f90` holds **five** separate `type(rsgpu), save ::
gpu_backend` instances (lines ~409, 628, 832, 1223, 1387), each initializing
its own context and uploading its own device state.

Target:
1. Pick `rsrec_cuda_plugin.f90` as the survivor (it carries the backend-format
   enum and stub handling). Port any capabilities that only exist in `rsgpu`
   (`block_dos`, `chebyshev_gf_eta` reconstruction entry points, precision
   selector) into it.
2. Introduce a **single shared context accessor**
   (`get_gpu_context(lattice_fingerprint)`) owned by the plugin module;
   `recursion` and `green` both use it. The five `save`'d instances in
   `green.f90` are removed. Context invalidation happens when the lattice /
   Hamiltonian fingerprint changes (same fingerprint idea as T5).
3. Delete `cuda/recursion_gpu_mod.f90` once nothing references it.
4. Do not modify `cuda/rsrec_gpu.cu` in this task beyond what the header
   consolidation strictly requires.

Validation: GPU regression run if hardware is available; otherwise build with
`ENABLE_CUDA_PLUGIN=ON` against the stub (`rsrec_stub.c`) and confirm the CPU
fallback paths are untouched.

---

## Phase 3 — Structure of the monster modules

### T9. Split the four largest modules using Fortran submodules
`reciprocal.f90` (222 KB), `recursion.f90` (200 KB), `lattice.f90` (195 KB),
`hamiltonian.f90` (164 KB). Use **submodules** so module interfaces (and all
`use` statements elsewhere) stay untouched:
- `recursion.f90` → parent module with the type + interfaces; submodules:
  `recursion_haydock` (recur, recur_b, crecal*, hop*, terminators zsqr/
  get_cinf/get_terminf), `recursion_chebyshev` (cheb_* drivers, dispatcher,
  moments), `recursion_transport` (velocity matmuls, recur_b_ij,
  chebyshev_recur_ij, stochastic/orbital moments).
- `hamiltonian.f90` → submodules: `hamiltonian_build` (build_* core),
  `hamiltonian_ccor` (all `*ccor*` routines — they are already a coherent
  cluster), `hamiltonian_paoflow_io` (rs2pao, build_from_paoflow*,
  export_rs_*), `hamiltonian_hubbard` (hubbard U/V).
- `lattice.f90`, `reciprocal.f90`: same pattern; propose the split boundaries
  from the routine clustering before executing, and get sign-off (report the
  proposed grouping first).
Rule: moving code only — zero logic edits in this task, so the regression
diff is exactly "none". Per ground rule 9, the derived type, its constructor
interface, `destructor`, `restore_to_default`, and the full
`procedure ::` binding list remain in the parent module; only implementations
move to submodules (`module subroutine` / `module procedure` bodies). Note:
`self.f90` and `symbolic_atom.f90` are **not** in scope for this split (ground
rule 10) even though they are mid-sized — leave them intact for the later
legacy audit.

### T10. Deduplicate the `calculation.f90` p2rs orchestration pairs
`post_processing_exchange` vs `post_processing_exchange_p2rs` and
`post_processing_conductivity` vs `post_processing_conductivity_p2rs` share
their recursion/moment scaffolding (both call `recur_b_ij` +
`chebyshev_recur_ij`). Extract the shared scaffold into a private helper
(`run_intersite_moments(...)`) parameterized by the Hamiltonian source
(native vs paoflow); the four public entry points become thin.

### T11. Audit `math.f90` (88 KB, 72 procedures)
Report — do not yet delete — procedures that are (a) uncalled, (b) duplicates
of intrinsics or BLAS/LAPACK wrappers that exist elsewhere in the code, or
(c) duplicated between `math.f90` and `include_codes/array/*` or
`abspinlib`. Produce the list with call counts; deletion happens in a
follow-up commit after human review.

---

## Phase 4 — Build system and options

### T12. Simplify the kernel option matrix
`ENABLE_MKL_BATCH` and `ENABLE_MKL_SPARSE` gate two of the CPU backends.
After T4 the backends share one driver, so collapse to a single
`ENABLE_MKL_KERNELS` option that compiles both (they link the same MKL);
keep runtime selection purely via the `cheb_backend` namelist key. Preserve
graceful degradation: if built without MKL, requesting an MKL backend must
produce the same fatal-with-hint logger message as today.

### T13. CMake target hygiene
Review `CMakeLists.txt` for the same modernization you applied in UppASD:
per-feature `INTERFACE` targets instead of global flags, correct usage
requirements on the CUDA plugin target, and make `RUN_REG_TESTS` register the
new T1 cases. No behavior change for default builds.

---

## Suggested execution order and effort

| Order | Task | Risk | Est. size |
|-------|------|------|-----------|
| 1 | T1 tests | none | small |
| 2 | T2, T3 dead code/hygiene | low | small |
| 3 | T6 dispatcher | low | small |
| 4 | T5 cache handle | medium | medium |
| 5 | T4 kernel unification | **high** | large |
| 6 | T7 green dedup | medium | medium |
| 7 | T8 GPU layer | medium | medium |
| 8 | T9 submodule split | low (mechanical) | large |
| 9 | T10–T13 | low | small–medium |

T4 is the highest-value and highest-risk task; do not start it before T1, T5
and T6 have landed, since they shrink the surface it has to touch.

## Definition of done (whole effort)
- All regression + SCF tests pass at baseline tolerances, including the new
  backend matrix from T1.
- `chebyshev_fast.f90` contains exactly one copy each of the apply/cherk/swap
  logic per mathematical variant (dense, BSR, MKL-batch, MKL-sparse, hoh).
- Exactly one Fortran GPU wrapper module; at most one live GPU context per
  rank.
- No source file above ~80 KB except mechanically split parents and the
  deliberately untouched `self.f90`.
- No commented-out executable code in `source/`.
- Net line count of `source/` is **lower** than at the start (tests excluded);
  no new runtime validation/defensive code on hot paths.
- Class-based structure (constructors, destructors, type-bound procedures)
  intact across all touched modules.
- Legacy routine bodies in `self.f90` / `symbolic_atom.f90` byte-identical
  except where a change was explicitly reported and approved.

## Progress checklist

Baseline:
- [x] Existing `tests/regression` baseline recorded: `bccFe_lanczos`,
  `bccFe_block`, and `bccFe_chebyshev` passed with `build/bin/rslmto.x`.

Tasks:
- [x] T1. Extend regression coverage to every backend being merged later.
- [x] T2. Dead code removal.
- [x] T3. Repo hygiene.
- [x] T4. Unify the Chebyshev CPU kernel family in `chebyshev_fast.f90`.
- [x] T5. Replace the module-global cache state with a cache handle type.
- [x] T6. Resolve operator selection once in the dispatcher.
- [x] T7. Deduplicate `green.f90`'s clone families.
- [x] T8. One GPU context layer.
- [x] T9. Split the four largest modules using Fortran submodules.
- [x] T10. Deduplicate the `calculation.f90` p2rs orchestration pairs.
- [ ] T11. Audit `math.f90`.
- [ ] T12. Simplify the kernel option matrix.
- [ ] T13. CMake target hygiene.
