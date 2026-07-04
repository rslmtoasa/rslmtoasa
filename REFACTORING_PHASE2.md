# RS-LMTO-ASA — Phase 2: Test Matrix, CI, Documentation, Developer Map

**Audience:** coding assistant working on branch `fable_v2` (or its successor).
**Precondition:** Phase 1 (T1–T13 + R1–R3) is complete and independently
verified: clean build, regression backend matrix and example suites pass
against pre-refactor references.
**Goal of this phase:** no further structural refactoring. This phase hardens
what exists: complete the test matrix (P1–P3), make CI actually run it
(P4–P6), document the refactored code (P7), and produce a developer map
(P8). P9–P10 are small cleanups.

All Phase 1 ground rules remain in force, in particular: KISS (rule 8),
class-based architecture (rule 9), legacy `self.f90`/`symbolic_atom.f90`
routines untouched (rule 10), one task per commit, regression suite green
after every task.

---

## Progress checklist

- [x] P4 — CI runs the regression backend matrix (`-DRUN_REG_TESTS=ON`),
      split into backend/example steps; per-case timeouts raised to 300s.
- [x] P2 — `launch_modes: ["serial", "mpi"]` schema; applied to one bulk
      (Pt2MnGa block), the surface (fccCu001), the impurity (B2FeCo), one
      chebyshev (Pt2MnGa chebyshev).
- [x] P1 — Close workflow coverage gaps: lanczos (nsp2 NaN bug found, see
      tests/KNOWN_ISSUES.md), spin dynamics (blocked, see KNOWN_ISSUES.md,
      no test added), paoflow (all 3 routes), orbital moments, hubbard,
      k-space ccor_2c/reciprocal-hoh, export paths (rs2pao/python).
- [x] P3 — GPU/MKL coverage strategy: MKL CI job (found and fixed a stale
      T12 option-name bug along the way — mkl_batch/mkl_sparse regression
      cases were silently skipping unconditionally since T12 collapsed the
      option names), CUDA compile-only CI job, `tests/run_gpu_matrix.sh`
      validated end-to-end on real GPU hardware.
- [x] P5 — CI trigger/matrix hygiene: `pull_request` trigger (self-build
      path, no upstream artifact for unmerged PRs), concurrency groups,
      "quick" label + nightly-full/PR-quick split, best-effort Windows
      test job (`continue-on-error`, unverified — no Windows runner
      available here; see tests/README.md for the documented fallback).
- [ ] P7 — Docstring blocks for `recursion.f90`, `hamiltonian_*`,
      `lattice_*`, `green.f90`/`calculation.f90` gaps, `reciprocal_*`.
- [x] P8 — `docs/DEVELOPER_MAP.md` (all 7 sections, every pointer verified
      against the code) + repo-root `CLAUDE.md`.
- [x] P6 — Reference-update workflow documentation in `tests/README.md`.
- [x] P9 — Executed the T11 math-audit deletions: 26 confirmed-uncalled
      routines removed from `math.f90`, plus `erf_qe`/`erfc_qe` (became
      dead as a direct consequence). 818 lines removed, clean rebuild,
      full matrix green.
- [x] P10 — Housekeeping: R2 executable-stack note expanded to explain the
      gfortran-trampoline cause (`chebyshev_fast.f90`'s internal-procedure
      callbacks); planning docs moved to `docs/dev/`; MKL fatal-with-hint
      behavior verified manually (both `mkl_batch` and `mkl_sparse` exit 1
      with a clear "requires CMake option ENABLE_MKL_KERNELS=ON" message
      against a fresh `-DENABLE_MKL_KERNELS=OFF` build) rather than adding
      an automated negative test — that would need a new "expect-fatal"
      case type plus an inverse `requires_cmake_option`-style gate (skip
      when the option *is* on), more than the "if trivial" bar the doc set.
- [ ] P11 — Python-binding awareness lens (standing, applied throughout;
      final blocker list recorded in the developer map).

---

## Part A — Completing the test matrix

Current coverage (verified): regression backend matrix (10 cases: 5 CPU
Chebyshev backends × hoh × ccor_2c + block sp/dp), SCF examples (19 cases:
bccFe nsp2/3/4 × block/chebyshev × hoh, fccCu001 surface, B2FeCo impurity,
Pt2MnGa Heusler, k-space SCF), postproc examples (6 cases: exchange,
exchange+conductivity, bands, DOS), and 6 MPI cases at 2–4 ranks.

### P1. Close the workflow coverage gaps
Add example-suite cases (same `cases.json` + reference mechanism) for the
workflows that currently have **zero** coverage:

1. **`recur = lanczos`** in the modern matrix. Only the legacy `oneliner.sh`
   baseline exercises Lanczos; add `Example_bulk_bccFe_nsp2_lanczos` (and a
   hoh variant) to the SCF suite so all three recursion modes live in one
   framework.
2. **Spin dynamics** (`processing_sd`, calculation.f90): one short-trajectory
   case (few time steps, small cell). Compare final moments/energies to a
   reference with a relaxed tolerance if the integrator is
   accumulation-order sensitive — determine empirically and document the
   chosen tolerance in the case entry.
3. **PAOFLOW import path** (`post_processing_paoflow2rs`, `exchange_p2rs`,
   `conductivity_p2rs`): T10 refactored this scaffolding but nothing tests
   it. Add one small paoflow-input case per route. If no small paoflow input
   file exists in `example/`, generate a minimal one once, commit it as test
   data, and note its provenance in the case README.
4. **Orbital moments** (`post_processing_orbital_modern`).
5. **Hubbard corrections** (`hamiltonian_hubbard.f90`): one SCF case with
   `hubbard` on (bcc Fe with a small U is fine — this is a smoke-and-compare
   test, not physics).
6. **k-space post-processing with ccor_2c / hoh** (the recent tetrahedron
   work): one bands or DOS case with `ccor_2c = .true.`, one with
   reciprocal-mode hoh.
7. **Export paths** (`hamiltonian` export select-case, calculation.f90:641):
   one case per export format that simply checks the exported file exists
   and parses.

Keep every case minimal (small cell, few recursion steps / energy points,
target < 2 min serial). Regenerate references with
`tests/generate_references.py` from the current trusted build and commit
them. Do NOT add in-source verification code (ground rule 8 / T1 note).

### P2. Serial-vs-MPI consistency cases
The 6 MPI cases run *only* under MPI, so a rank-decomposition bug that
changes results relative to serial would not be caught. Add a
`launch_modes: ["serial", "mpi"]` mechanism to the case schema and
`run_test.py`/CMake registration: a case listing both modes registers two
CTest entries (`<name>` and `<name>_mpi`) compared against the **same**
reference file. Apply it to a representative subset (one bulk, the surface,
the impurity, one chebyshev) rather than everything — KISS, CI minutes are
finite. Register MPI tests with `PROCESSORS` property so ctest scheduling
accounts for ranks.

### P3. GPU and MKL coverage strategy
GitHub-hosted runners have no GPU, and MKL is not installed by default.
Concretely:

1. **MKL (runnable in CI):** add an Ubuntu CI job that installs Intel oneAPI
   MKL (apt repo `intel-oneapi-mkl-devel`), configures with
   `-DENABLE_MKL_KERNELS=ON`, builds, and runs *only* the two MKL regression
   cases plus one MKL example case. This finally exercises
   `mkl_batch`/`mkl_sparse` in CI.
2. **CUDA (compile-only in CI):** add a job installing `nvidia-cuda-toolkit`
   (no driver needed to compile), building with `-DENABLE_CUDA_PLUGIN=ON`,
   and running the *CPU* regression subset against that binary (plugin
   present but unused). This catches interface drift between
   `rsrec_cuda_plugin.f90`, `rsrec_cuda.h`, and `rsrec_gpu.cu` at every push
   — the most common GPU breakage mode — without needing hardware.
3. **GPU (runnable, off-CI):** add `tests/run_gpu_matrix.sh` that runs the
   regression matrix plus the `_gpu` green/DOS routes on a local CUDA
   machine and prints a PASS/FAIL table. Document in the tests README that
   this is the manual pre-release step on a GPU workstation. Optionally
   support a self-hosted runner label (`runs-on: [self-hosted, gpu]`) behind
   `workflow_dispatch` so it can be triggered when such a runner is
   registered — but do not make any required check depend on it.

---

## Part B — CI workflows

### P4. Make CI run the full matrix, not just examples
`tests.yml` currently runs `ctest --label-regex example` only — the
`Regression_*` backend matrix (labels `regression;backend`) never runs in CI.
Fix:
- Configure step: add `-DRUN_REG_TESTS=ON`.
- Split the ctest invocation into two steps with separate summary sections:
  `--label-regex backend` (regression matrix) and `--label-regex example`.
  Keep `--parallel 1`; these are OpenMP jobs.
- Raise per-case `timeout` values in the three `cases.json` files from 120 s
  to 300 s. Independent verification showed load-sensitive timeouts at 120 s
  on slower machines; the timeout exists to catch hangs, not to enforce
  performance.

### P5. Trigger and matrix hygiene
- Add a `pull_request` trigger path: the current chain
  (binaries.yml → workflow_run → tests.yml) only covers pushes. Simplest
  KISS-compliant fix: give tests.yml a `pull_request` trigger that builds
  from source itself (one `cmake --build` step, reusing the same steps) when
  not invoked via workflow_run. Alternatively have binaries.yml also run on
  `pull_request` — pick whichever keeps the YAML smaller, and add
  `concurrency` groups with `cancel-in-progress: true` to both workflows.
- Windows is built by binaries.yml but never tested. Add a Windows test job
  running only the fast SCF subset under msys2 python (install `f90nml` via
  pip). If path/namelist issues make this disproportionate work, timebox it:
  report the blocker and mark Windows as build-only in the README instead.
- Add the P3 MKL and CUDA-compile jobs to the matrix.
- Nightly `schedule:` run of the **full** matrix (all labels, all cases) on
  Linux + macOS, while `pull_request` runs a `quick` labelled subset (tag
  ~10 representative cases with an extra `quick` label). This keeps PR
  turnaround under ~15 min while nothing rots.

### P6. Reference-update workflow documentation
Add a short `tests/README.md` section: how to intentionally regenerate
references when physics changes (run `generate_references.py` against a
trusted build, review the diff of reference files like any code change,
never regenerate to "fix" an unexplained failure). This is process
documentation, not tooling.

---

## Part C — Documentation

### P7. Docstring blocks for the refactored modules
The codebase already has a Doxygen-ish convention (`!> @brief`, `@author`,
`!> ...` blocks — see `element.f90`, `chebyshev_fast.f90` which is at ~90%
coverage). Bring the terse files up to the same standard. Measured current
state: `recursion.f90` ~0/35 routines documented, `hamiltonian_build.f90`
~0/22, `green.f90` 18/27, `calculation.f90` 13/23; the other new submodule
files (`lattice_*`, `hamiltonian_*`, `reciprocal_*`) are similarly bare.

For every public and type-bound procedure in these files add a header block:

```fortran
!> @brief One-line purpose.
!> @details 2–6 lines: the math (which operator/equation, notation as in the
!>          docs/ccor_2c.md style), the algorithm choice if non-obvious, and
!>          which workflow(s) call it.
!> @param[in]  x   meaning, shape convention (e.g. site-major (ld,nb))
!> @param[out] y   meaning
!> @note side effects: which caches it reads/populates, MPI collectives, I/O.
```

Rules (KISS applies to prose too):
- Document **what and why**, never restate the code line-by-line. No inline
  comment noise inside loops except where the physics is genuinely opaque
  (e.g. the hoh two-sweep, the KPM doubling formulas — those two deserve a
  short block each if not already present).
- Priority order: `recursion.f90` → the four `hamiltonian_*` submodule files
  → `lattice_*` submodule files → `green.f90` gaps → `calculation.f90` gaps
  → `reciprocal_*`. Legacy `self.f90`/`symbolic_atom.f90` bodies are still
  off-limits; a file-level header describing what the legacy module does is
  allowed and welcome.
- Also document the two module-level contracts introduced in Phase 1: the
  `ham_apply_t` callback interface (semantics of `alpha/beta`, y = α·H·x1 +
  β·x0) and `cheb_cache_t` (fingerprint invalidation semantics, single-rank
  assumption).
- One commit per file or file-group; no code changes mixed in, so the diff
  is pure comments and trivially reviewable.

---

## Part D — Developer map

### P8. Create `docs/DEVELOPER_MAP.md`
A cheat-sheet for humans and AI assistants: for each workflow, the entry
point, the classes involved, the hot kernels, and where new functionality
plugs in. Structure to produce (verify every pointer against the code while
writing; the skeleton below reflects the post-T9 layout):

1. **Execution skeleton:** `main.f90` → `calculation` type
   (`calculation.f90`): `pre_processing` (bravais / newclubulk / newclusurf /
   buildsurf) → `processing` (SCF drivers per `calctype` B/S/I) →
   `post_processing` (exchange, conductivity, bands, DOS, orbital, paoflow,
   sd). Table: namelist key → routine.
2. **Per-workflow chains**, one short section each, in the form
   *entry → orchestration → kernel*, e.g.:
   - Real-space SCF: `processing` → `self` (SCF loop, legacy core) →
     `hamiltonian_build` → `recursion` (`recur`: lanczos `crecal`/block
     `recur_b`/chebyshev `cheb_moments_cpu` dispatcher) →
     `chebyshev_fast.f90` kernels or GPU plugin → `green` → `charge`/`mix`.
   - k-space: `reciprocal*` route incl. tetrahedron integration and where
     ccor/hoh enter.
   - Exchange (J_ij), conductivity, orbital moments: intersite chain
     `run_intersite_moments` → `recur_b_ij`/`chebyshev_recur_ij` →
     `green%calculate_intersite_gf_core` / `bgreen` → `exchange.f90` /
     `conductivity.f90`.
   - Spin dynamics: `processing_sd` → `spin_dynamics.f90` →
     `include_codes/abspinlib`.
3. **Kernel inventory:** the `ham_apply_t` implementations (dense fused,
   BSR, MKL batch, MKL sparse, hoh two-sweep), the shared
   `cheb_three_term_moments` driver, `hsweep_sp`/`happly_sp`, block-Lanczos
   fast kernels, and the CUDA plugin surface (`rsrec_cuda_plugin.f90` ↔
   `cuda/rsrec_gpu.cu`), with one line each on layout conventions
   (site-major, fp32) and the cache/fingerprint lifecycle.
4. **"Where to add X" recipes** — concrete insertion points with file names
   for at least these three (requested) plus any that fall out naturally:
   - **Porting `strux_lib` to GPU:** structure-constant generation lives in
     `include_codes/strux_lib` called from `lattice_strux.f90`; describe the
     data in/out at that boundary, and recommend mirroring the existing
     plugin pattern (C API + `iso_c_binding` wrapper module + CPU fallback,
     as in `rsrec_cuda_plugin.f90`) rather than inventing a second GPU
     convention.
   - **Linear-response TDDFT:** the χ₀ building blocks are the retarded
     Green's functions on the energy mesh (`green.f90` cores,
     `energy_mesh`/`en` handling in `bands.f90`) and the intersite GF
     machinery; describe where a new `post_processing_susceptibility` would
     register in `calculation.f90` (`check_post_processing`,
     `prepare_post_processing_stack`) and which existing pieces it reuses.
   - **Lehmann-representation Green's functions:** entry point belongs in the
     **reciprocal family**, not `green.f90` — the required eigenpairs only
     exist on the k-space diagonalization route. Plan for a new
     `reciprocal_green.f90` submodule of `reciprocal_mod` (post-T9 pattern),
     consuming eigenvalues/eigenvectors from the existing `reciprocal*`/
     `bands` machinery. Design constraint to record in the map: the routine
     should return the Green's function in the same container/layout the
     real-space routes produce (the `green` type's per-site/intersite
     arrays), so that downstream consumers — exchange, transport/
     conductivity, damping, and future linear-response — can use either
     spectral route interchangeably. Implementation details beyond this
     interface contract are deliberately deferred to the feature itself; the
     map should state the insertion point and the contract, nothing more.
5. **Testing map:** which suite/label covers which workflow, how to add a
   case, how to regenerate references (link P6).
6. **Conventions box:** KISS rules, class-based architecture requirement,
   legacy-file policy, commit conventions — one screen, lifted from the
   Phase 1 ground rules so a fresh AI session inherits them by reading this
   one file.

Also add a repo-root `CLAUDE.md` (or extend it if present) pointing at
`docs/DEVELOPER_MAP.md`, `REFACTORING_PLAN.md`, and the tests README, and
restating the five most important ground rules in ~10 lines.

---

## Part E — Small follow-ups

### P9. Execute the T11 math audit deletions
`MATH_AUDIT_T11.md` lists ~16 uncalled helpers (`zeros1/2`, `ones1/2`,
`diag1/2`, `kron`, `eye`/`ceye`, (c)`anticommutator`, `is_hermitian_r4`,
`unique_int`, `sph2cart`, `nm2rho`, `fermi_function`, …). The report-only
phase is done and the deletions are hereby approved: remove them, keeping
any that turn out to be referenced via generic interfaces (re-check
interface blocks, not just call sites). Update the audit file to record what
was removed.

### P10. Housekeeping
- Verify the executable-stack note from R2 made it into the user-facing
  build documentation (a note was added to `docs/source/getting_started.rst`
  — confirm it explains the gfortran/trampoline cause and the
  hardened-system symptom).
- Move `REFACTORING_PLAN.md`, `REVIEW_NOTES_T1-T6.md`, `MATH_AUDIT_T11.md`
  into `docs/dev/` to declutter the repo root, updating any links.
- Confirm `ENABLE_MKL_KERNELS` produces the same fatal-with-hint logger
  message when an MKL backend is requested in a non-MKL build (T12
  acceptance criterion; add a tiny negative test if trivial — a case that
  expects the fatal — otherwise verify manually and note it).

### P11. Keep the door open for a Python module (design awareness, not implementation)
A future goal is exposing the code as a Python module for notebooks,
high-throughput screening, and tutorials. **Do not implement anything in
Phase 2**, but while executing P1–P10, observe and enforce the properties
that will make it cheap later:

1. **Library/driver separation.** Everything a Python binding would call
   must be reachable without `main.f90`: construct `control`/`lattice`/
   `charge`/... objects, run SCF or a post-processing step, read results
   from the objects. Where a task touches code that only works via the full
   `calculation%process()` chain plus file I/O in the working directory,
   note it in `DEVELOPER_MAP.md` under a "Python-binding blockers" list
   (do not fix it yet). Known suspects to check while documenting:
   hard-coded input/output filenames, `stop`/`error stop` in library code
   paths (a Python session must not be killed by a fatal — the `g_logger%
   fatal` pattern is the thing to audit), and any remaining bare module-level
   mutable state (the Phase-1 cache/context handles are already the right
   shape).
2. **Results as data, not files.** When writing P7 docstrings, explicitly
   document which arrays on which objects hold the physically meaningful
   results (moments, J_ij, DOS, bands, conductivity tensors) — that
   documentation doubles as the binding surface specification.
3. **Binding route (record as recommendation in the map, decide later):**
   given the class-based design, `f90wrap` is the natural fit (it wraps
   derived types with type-bound procedures into Python classes); the
   fallback is a thin explicit `iso_c_binding` API layer (the
   `rsrec_cuda_plugin` C-API pattern, in reverse) wrapped with ctypes/cffi
   — more work but zero toolchain magic. A pybind11 route would only make
   sense via the C API. Add a short "Python module outlook" section to
   `DEVELOPER_MAP.md` stating this and linking the blocker list.
4. **CMake:** when doing P4/P5, keep the object-library structure such that
   a `SHARED` library target of the physics core (everything except
   `main.f90`) could be added with ~5 lines. If a trivial
   `-DBUILD_SHARED_CORE=ON` option falls out naturally from T13's target
   hygiene, adding it (build-only, untested, undocumented to users) is in
   scope; anything more is not.

---

## Suggested order

| Order | Task | Why first |
|-------|------|-----------|
| 1 | P4 | CI must run what already exists before adding more |
| 2 | P2 | small schema change, unblocks P1 case design |
| 3 | P1 | the big coverage push |
| 4 | P3 + P5 | new CI jobs once cases exist |
| 5 | P7 | docs on a now-stable code base |
| 6 | P8 | map written last = most accurate |
| 7 | P6, P9, P10 | cleanups |

P11 has no slot of its own — it is a standing lens applied while doing P1,
P4/P5, P7, and P8.

## Definition of done
- Every `processing`/`post_processing` select-case branch in
  `calculation.f90` has at least one CI-run test case (or a documented,
  justified exemption).
- CI runs regression + example matrices on Linux and macOS for every PR
  (quick subset) and nightly (full), plus MKL-run and CUDA-compile jobs.
- A serial-vs-MPI consistency pair exists for at least four representative
  cases.
- Every public routine in the Phase-1-refactored files carries a `!>`
  docstring block; no line-noise comments added.
- `docs/DEVELOPER_MAP.md` exists, and every file/routine pointer in it is
  verified against the code.
- Repo root contains a `CLAUDE.md` handing future sessions the ground rules.
- `DEVELOPER_MAP.md` contains the "Python module outlook" section with a
  concrete binding-blocker list gathered during Phase 2.
