# RS-LMTO-ASA Developer Map

A cheat-sheet for humans and AI assistants: for each workflow, the entry
point, the classes involved, the hot kernels, and where new functionality
plugs in. Every pointer below was checked against the code at the time of
writing (Phase 2, P8) — if something has moved, trust the code and fix this
file, not the other way around.

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

| `&calculation` key | Value | Routine (`calculation.f90`) | Notes |
|---|---|---|---|
| `pre_processing` | `'bravais'` | `pre_processing_bravais` | Bulk SCF (`calctype='B'`). Runs `self%run()`. |
| `pre_processing` | `'buildsurf'` | `pre_processing_buildsurf` | Surface SCF (`calctype='S'`). `build_surf_full()`, no `newclu()`. |
| `pre_processing` | `'newclubulk'` | `pre_processing_newclubulk` | Impurity-in-bulk (`calctype='I'`, bulk host). `newclu()`, no `build_surf_full()`. |
| `pre_processing` | `'newclusurf'` | `pre_processing_newclusurf` | Impurity-in-surface. `build_surf_full()` + `newclu()`. **No example/test input anywhere uses this route** (see `tests/KNOWN_ISSUES.md`). |
| `pre_processing` | `'buildinterface'` | `pre_processing_buildinterface` | Two-sided layered/interface SCF (`calctype='L'`, B7.5). `build_interface_full()` (region A \| active \| region B), then `surfmat()` (kernel reused unchanged) with its one-sided registry overwritten by `charge%build_interface_registry()`. Per-iteration Madelung update is `charge%interfacepot`, not `surfpot` (`self.f90` dispatch). `buildsurf` itself is untouched and remains the permanent one-sided regression oracle. |
| `processing` | `'sd'` | `processing_sd` | Spin dynamics. **Known bug:** its pre-processing block is hardcoded to the `newclusurf` sequence regardless of the case's actual route — see `tests/KNOWN_ISSUES.md`. |
| `post_processing` | `'exchange'` | `post_processing_exchange` | Real-space intersite J_ij/D_ij. |
| `post_processing` | `'exchange_p2rs'` | `post_processing_exchange_p2rs` | Same, Hamiltonian sourced from a PAOFLOW-format import instead of `build_bulkham()`. |
| `post_processing` | `'conductivity'` | `post_processing_conductivity` | Real-space conductivity tensor. |
| `post_processing` | `'conductivity_p2rs'` | `post_processing_conductivity_p2rs` | Same, PAOFLOW-imported Hamiltonian. |
| `post_processing` | `'paoflow2rs'` | `post_processing_paoflow2rs` | Base PAOFLOW import path: builds the lattice normally but replaces `hamiltonian%build_bulkham()` with `hamiltonian%build_from_paoflow_opt()` (reads `paoham.dat`). |
| `post_processing` | `'orbital_modern'` | `post_processing_orbital_modern` | Orbital moments via `recursion%chebyshev_orbital_mod()` — always Chebyshev, regardless of `control%recur`. Loops over *every* atom in the cluster (`this%lattice%kk`), so runtime scales with cluster size, not recursion depth. |
| `post_processing` | `'band_structure'` / `'density_of_states'` | (see `prepare_post_processing_stack` callers) | k-space route via `reciprocal_mod`, not the real-space recursion machinery at all. |
| `post_processing` | `'bsf'` | `post_processing_bsf` → `reciprocal%calculate_bsf` (`reciprocal_bsf.f90`) | Bloch spectral function A(k,E) = −1/π Im Tr G(k,E+iη) along the canonical spglib k-path (milestone B3). Consumes the B2 engine's `dyson_kspace_inverse` per (k,E) (Σ=0 ⇒ backend E; Σ-ready for CPA/DMFT). η = `&reciprocal` green_eta, E grid = n_energy_points/dos_energy_min,max, path = `&kpath` nk_per_segment. Writes `bsf.dat` (total/up/down) + `bsf_bands.dat` overlay. Partial-trace convention in `bsf_kernel.f90` (`bsf_spectral_trace`). |
| `post_processing` | `'kspace_green'` | `post_processing_kspace_green` | B2 validation driver: fills `green%gij` via recursion then via the k-space engine (`reciprocal%fill_green`, backend E + D≡E check) and cross-checks on-site DOS / m_z. Report-only. |
| `post_processing` | `'frozen_magnon'` | `post_processing_frozen_magnon` | Sweeps `hamiltonian%q_ss` over a `&frozen_magnon` q-list, preferably from `q_file` (`q_coordinates='cartesian'` for `2*pi/alat` Cartesian components or `'direct'` for reciprocal-lattice coordinates), writing total energy, band energy, per-sublattice moment magnitude, and `omega(q)` to `frozen_magnon.dat`. `mode='mft'` (default) converges SCF once at the reference point, reuses that potential for a single-iteration band-energy pass at every other q, and computes `omega` from band-energy differences; `mode='scf'` re-converges at every q and computes `omega` from total-energy differences. `branch_mode='auto'` builds multi-sublattice magnon branches in `frozen_magnon_branches.dat`/`frozen_magnon_modes.dat` via the direct GBT frozen-magnon method (second derivatives of the force-theorem band-energy surface w.r.t. sublattice cone angles; Essenberger PRB 84, 174425 Eq. 26). **Single-sublattice is validated; the multi-sublattice acoustic branch is not yet gapless at Γ — see `tests/KNOWN_ISSUES.md`, deferred to B11.** See `docs/dev/plans/B1_gbt_frozen_magnons.md` §1.5. |

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

### k-space (reciprocal)
```text
post_processing_band_structure / _density_of_states
  -> reciprocal_mod (reciprocal.f90 + reciprocal_{lifecycle,fourier,bands,dos,projection}.f90)
       -> reciprocal_fourier.f90: k-space H(k) built via Fourier transform of
          real-space hoppings; branches on hamiltonian%ccor_2c and
          hamiltonian%hoh (both bypass the "plain" Fourier sum — the recent
          tetrahedron-integration work reads these same two flags)
       -> reciprocal_bands.f90 / reciprocal_dos.f90: diagonalization,
          tetrahedron or gaussian DOS integration (control via &reciprocal
          dos_method)
       -> reciprocal_projection.f90: orbital-projected DOS, band moments
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

### Spin dynamics
```text
calculation%processing_sd()
  -> self(bands_obj, mix_obj)         ! constructed but NOT run() here —
                                       ! spin_dynamics%sd_run() calls self%run() itself
  -> spin_dynamics(self_obj)          [spin_dynamics.f90]
       -> sd_run(): per step, self%run() -> bands%calculate_magnetic_torques()
          -> asd_pred_euler / asd_corr_euler (Landau-Lifshitz-Gilbert integrator)
       -> include_codes/abspinlib/abSpinlib.f90   (lower-level ASD numerics)
```
See `tests/KNOWN_ISSUES.md` for why this route currently has no test case
(its pre-processing is hardcoded to the never-exercised `newclusurf` shape).

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

### Linear-response TDDFT
The χ₀ building blocks are the retarded Green's functions on the energy mesh
(`green.f90`'s per-site/intersite cores) and the `energy`/`en` mesh handling
in `bands.f90`. A new `post_processing_susceptibility` would register in
`calculation.f90` the same way `exchange`/`conductivity` do: add the string
to `check_post_processing`'s allowed list (line ~1195) and add a case to the
`post_processing` select-case in `process()` (line ~203), then reuse
`prepare_post_processing_stack` for the common pre-processing/Hamiltonian
setup rather than duplicating it.

### Lehmann-representation Green's functions
Entry point belongs in the **reciprocal family**, not `green.f90` — the
required eigenpairs only exist on the k-space diagonalization route
(`reciprocal_bands.f90`). Plan: a new `reciprocal_green.f90` submodule of
`reciprocal_mod` (does not exist yet, following the post-T9 submodule
pattern — see `reciprocal_{lifecycle,fourier,bands,dos,projection}.f90` for
the shape), consuming eigenvalues/eigenvectors already computed by
`reciprocal_bands.f90`. Design constraint: return the Green's function in
the same container/layout the real-space routes use (the `green` type's
per-site/intersite arrays), so downstream consumers (exchange, conductivity,
damping, future linear-response) can use either spectral route
interchangeably. Implementation details beyond this interface contract are
deliberately out of scope here.

---

## 5. Testing map

| Workflow | Suite / label | Case-file location |
|---|---|---|
| Regression backend matrix (Chebyshev CPU backends × hoh × ccor_2c, block sp/dp) | `regression;backend` | `tests/regression/cases.json` |
| Bulk/surface/impurity SCF, all recursion modes | `example;scf` | `tests/scf/cases.json` |
| Exchange, conductivity, bands, DOS, orbital moments, PAOFLOW import | `example;postproc` | `tests/postproc/cases.json` |
| Serial-vs-MPI consistency (representative subset) | `launch_serial` / `launch_mpi` | `"launch_modes": ["serial","mpi"]` field, either `cases.json` |
| PR-fast subset | `quick` | `"quick": true` field, either `cases.json` |
| MKL kernels (`mkl_batch`/`mkl_sparse`) | gated by `requires_cmake_option: ENABLE_MKL_KERNELS` | regression + `Example_bulk_bccFe_nsp2_chebyshev_mkl_batch` in `scf/cases.json` |
| CUDA plugin | compile-only in CI (`cuda_compile` job); real-GPU consistency via `tests/run_gpu_matrix.sh` | n/a (manual, off-CI) |

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
- **One task, one commit.** Commit messages reference the task/phase ID
  (e.g. `test(P1.3): ...`, `ci(P5.1): ...`).
- **Bit-level behavior is the contract.** The regression + example suites
  must pass at the same tolerances after every task — run them before
  starting (baseline) and after every change.

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
4. **No remaining bare module-level mutable state found beyond the
   Phase-1-introduced caches** (`cheb_cache_t`, the CUDA plugin's context
   singleton) — those are already the right shape (fingerprinted, explicit
   invalidation) for a future binding to reason about, so no new blocker
   there.

**Results as data:** the physically meaningful outputs per workflow are the
arrays this map's P7 docstring pass will document on each object — e.g.
`self`/`bands`'s moment arrays, `exchange`'s J_ij/D_ij, `conductivity`'s
tensor, `reciprocal_dos`'s DOS arrays. That documentation doubles as the
binding surface specification once P7 is done.
