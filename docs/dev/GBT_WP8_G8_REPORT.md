# GBT WP8 / Gate G8 report — q little-group symmetry and q lifecycle

Date: 2026-08-06. Branch `fable_v2_gbt_v2`, on top of `14d5fe0` (final G7: PASS).

## Outcome

**Gate G8: PASS.** The full chemical-cell BZ (`q_symmetry_policy = 'full_bz'`)
remains the default and the oracle every other policy is checked against,
bit-identical to the pre-WP8 WP0 behaviour. Two opt-in reductions were added
— `'little_group'` (rebuild per q) and `'little_group_common'` (one mesh for
an entire multi-q sweep, valid for every point in it) — gated by a mesh cache
key that includes lattice, mesh dimensions, offset, the declared q-set, and
the symmetry policy itself, so a sweep can never reuse one q's (or one
policy's) mesh for another. All q-dependent eigensystem/DOS/density state is
invalidated whenever the cache key changes.

## 1. What WP8 found

`b1caf4b` ("fix(B1r.6): reduce the k-mesh by the little group of q_ss") is a
commit on the sibling branch `fable_v2_gbt` (not an ancestor of
`fable_v2_gbt_v2`) implementing the same physical idea: a spin spiral lowers
the crystal symmetry to the little group of `q_ss`
(`docs/dev/plans/B1_gbt_frozen_magnons_v2.md` section 3.1), so reducing the
k-mesh with the full chemical-cell point group is an invalid BZ integral for
`q_ss != 0`. Per the task prompt, it was reviewed (not merged) before writing
this task's code, in particular its recorded "per-q ownership" gap, from its
own `tests/KNOWN_ISSUES.md` entry:

> `frozen_magnon` k-space sweeps build the reduced k-mesh once, from row 1's
> `q_ss`, and reuse it for every other q in the sweep.

That entry does not exist on this branch (`b1caf4b` predates the WP0–WP7
guard sequence this branch is built on and was never merged here), but the
underlying architectural gap is real and applies equally to the current
`calculation.f90::post_processing_frozen_magnon_acoustic`/`_auto`. Two
findings from reviewing `b1caf4b`'s code, not just its text:

1. Its `generate_little_group_kpoint_mesh` is a good, minimal use of
   spglib's `spg_get_stabilized_reciprocal_mesh` and was ported (single-q
   case) essentially unchanged into `source/reciprocal_bands.f90`.
2. Its dispatch made little-group reduction *unconditional* for any nonzero
   `q_ss` — no oracle policy, no way to keep the WP0 full-BZ guard as the
   default. That is incompatible with this branch's explicit instruction to
   "keep full BZ as the oracle" and would have silently changed behaviour
   for every existing nonzero-q GBT run. WP8 instead makes little-group
   reduction an explicit opt-in (`q_symmetry_policy`), and generalizes the
   single-q stabilizer call into a multi-q one so a sweep can ask for the
   subgroup common to its whole q-list in one call, which resolves the
   per-q ownership gap at the source rather than papering over it with a
   per-q rebuild alone.

A second, more consequential finding came from tracing the call graph rather
than the two `calculation.f90` mesh-build sites named in the task: **every**
k-space Hamiltonian build calls
`reciprocal%force_full_bz_for_nonzero_q_gbt` (from
`build_kspace_hamiltonian`, `source/reciprocal_fourier.f90:261`), and that
routine unconditionally replaced the current mesh with a fresh full BZ on
*every* call — once per frozen-magnon probe. Fixing only the
`calculation.f90` mesh-build sites (the two the prompt names) would have
built a correct little-group mesh right before the sweep loop and then had
it silently discarded and rebuilt as the full mesh on the very next
`frozen_magnon_probe_energy` call. §3 covers the fix; it was found and fixed
during this task, not left for a follow-up, because without it the whole
policy would be inert.

## 2. What changed

### 2.1 `q_symmetry_policy` (new `&reciprocal` key)

`source/include_codes/namelists/reciprocal.f90`,
`source/reciprocal_lifecycle.f90::build_from_file`. Three values, validated
with a fallback-to-default warning like the module's other enum keys:

| policy | meaning |
| --- | --- |
| `full_bz` (default) | Unchanged WP0 behaviour: any nonzero-q GBT run forces the full, unreduced chemical-cell BZ. The oracle. |
| `little_group` | Reduce by the little group of the *current* `hamiltonian%q_ss`; rebuilds whenever `q_ss` changes. |
| `little_group_common` | Reduce by the subgroup common to a caller-declared q-list (`ensure_kpoint_mesh(..., q_list_cart=...)`); one mesh valid for an entire multi-q sweep. |

Every dispatch point defaults to `full_bz`'s exact pre-WP8 code path when the
policy is left at default, so ordinary (non-opted-in) runs are unaffected.

### 2.2 Little-group mesh generation

`source/spglib_interface.f90`: `get_rotation_matrices` (crystal point-group
rotations, translations discarded) and
`get_little_group_kpoint_mesh_with_points` (generalized beyond `b1caf4b` to
accept `q_frac(:, num_q)` — one or several q-points — and call
`spg_get_stabilized_reciprocal_mesh` with `num_q` simultaneously; `num_q=1`
reproduces the single-q little group, `num_q>1` gives the subgroup common to
every column, which is what a multi-q sweep needs).

`source/reciprocal_bands.f90::generate_little_group_kpoint_mesh` — ported
from `b1caf4b` and extended with an optional `q_list_cart(:, :)`: absent,
uses the single current `hamiltonian%q_ss`; present, reduces by the common
subgroup of every column. Falls back to the full mesh, with a single logged
message, if spglib is unavailable or the little group cannot be determined
— same contract as every other mesh-generation fallback in this module.

### 2.3 The mesh cache key

`source/reciprocal.f90` (type fields), `source/reciprocal_bands.f90`
(`ensure_kpoint_mesh`, `is_declared_nonzero_q_gbt`). Per the task's explicit
requirement, the cache key is not just "is a mesh allocated":

```
mesh_cache_dims(3), mesh_cache_offset(3), mesh_cache_lattice(3,3),
mesh_cache_policy, mesh_cache_q(3, :), mesh_cache_valid
```

`ensure_kpoint_mesh(this, mesh_dims, use_shift, q_list_cart)` is the single
entry point call sites should use instead of the historical
`if (.not. allocated(this%k_points))` guard (still present, unchanged, on
the `full_bz` path only — see §2.4). It compares every field of the key
against the current call's arguments plus `this%lattice%a` and
`this%q_symmetry_policy`; a mismatch (mesh dims, offset, lattice, policy, or
the declared q-set) triggers a rebuild through the same dispatch
`generate_reduced_kpoint_mesh`/`generate_little_group_kpoint_mesh`/
`generate_mp_mesh` use, then calls `invalidate_spectral_cache()` explicitly
— the mesh identity changed even when
`hamiltonian%operator_generation` did not, so H(k)/eigenpairs/DOS/density
must not survive on the old k-point count. A no-op (cache hit) touches
nothing.

One subtlety, found while first testing `little_group_common`: a pre-loop
`ensure_kpoint_mesh(..., q_list_cart=q_ss_cart)` call is made while
`hamiltonian%q_ss` is still the sweep's row-1 reference value (often zero),
so `has_nonzero_q_gbt()` alone — which only inspects the single *current*
`q_ss` — under-detects a caller-declared nonzero q-list and silently falls
through to the ordinary q=0 point-group mesh. `is_declared_nonzero_q_gbt`
fixes this: it is `has_nonzero_q_gbt()` OR'd with "GBT representation active
and any column of the supplied `q_list_cart` is nonzero". This was caught by
the new regression test (§4), not by inspection — see its failure mode in
§4.1.

### 2.4 `force_full_bz_for_nonzero_q_gbt` — the actual per-build gate

`source/reciprocal_lifecycle.f90`. This is the routine `build_kspace_hamiltonian`
calls on *every* Hamiltonian build, so it is the real enforcement point, not
the two `calculation.f90` mesh-build call sites (§1). Made policy-aware:

* `full_bz` (default): identical to the pre-WP8 body — force
  `use_symmetry_reduction/use_time_reversal = .false.` and rebuild the full
  MP mesh via `generate_mp_mesh()`, every call.
* `little_group`: delegates to `ensure_kpoint_mesh` — a no-op once the
  current `q_ss`'s mesh is already cached (set by the calculation.f90 per-q
  rebuild below), and a rebuild only on an actual `q_ss` change.
* `little_group_common`: must never discard a mesh a sweep already proved
  valid for its whole q-list. If a mesh is cached under
  `little_group_common` for the current mesh dimensions, it is left
  untouched. Only a genuinely missing common mesh (this routine invoked
  outside a sweep that pre-built one) falls back to the current q's own
  little group, with a logged note.
* The pre-existing k-path early return (`allocated(this%k_path)`, a band
  path is not a BZ integration mesh) is unchanged and checked first,
  independent of policy.

### 2.5 `calculation.f90` sweep fix

`post_processing_frozen_magnon_acoustic` (mft branch) and
`post_processing_frozen_magnon_auto`. Both previously built the k-space mesh
once, guarded by `.not. allocated(recip_obj%k_points)`, from whatever
`q_ss` happened to be set at that moment, and reused it for every other `iq`
— the exact bug class the task names. Restructured per policy:

* `full_bz` (default, unchanged): identical pre-WP8 code path, still gated
  by `.not. allocated(...)`.
* `little_group_common`: builds once, before the loop, via
  `ensure_kpoint_mesh(..., q_list_cart=q_ss_cart)` passing the *entire*
  sweep's q-list — one mesh proven valid for every row, never row 1's alone.
* `little_group`: no pre-loop build; each iteration calls
  `ensure_kpoint_mesh` right after setting `hamiltonian_obj%q_ss(:) =
  q_ss_cart(:, iq)`, so every probe gets its own q's mesh (cache-hit if two
  consecutive rows repeat a q).

## 3. Tests

`tests/regression/run_wp8_littlegroup.py` +
`tests/regression/wp8_littlegroup/` (deck adapted from the `b1caf4b`
fixture: bcc-Fe, `alat=2.8612`, k-space `frozen_magnon(mode='mft')`, 8x8x8
mesh, `strux_backend='strux_lib'` and `nsp=3` — WP4's strux_lib requirement
and WP6c's GBT+SOC guard, both newer than `b1caf4b`, meant its original
`nsp=4` deck no longer runs under `gbt_single_q` and had to be adjusted).
Registered as ctest `WP8LittleGroup` (`regression;gbt;wp8`, 1800s timeout).

| # | case | what it checks | result |
| --- | --- | --- | --- |
| 1 | axial | q on the bcc Gamma-H line (little group C4v): `little_group` vs `full_bz` oracle | eband diff **0.00e+00 Ry**; mesh strictly reduced (100 from 512) |
| 2 | generic | q off every symmetry axis/plane (trivial little group): `little_group` vs oracle | eband diff **0.00e+00 Ry**; mesh unreduced (512 from 512) |
| 3 | shifted | case 1's q on a shifted Monkhorst-Pack mesh | eband diff **0.00e+00 Ry**; mesh still reduced (100 from 512) |
| 4 | q / -q, common subgroup | 4-row sweep (q=0 ref, +q_axial, -q_axial, q_generic) under `full_bz`, `little_group`, `little_group_common` | every policy agrees with the oracle row-for-row to **0.00e+00 Ry**; `little_group` rebuilds exactly 3 times (one per distinct nonzero q); `little_group_common` builds **exactly once** (log-verified: "reduced by the little group common to 4 q-point(s) to 512 irreducible k-points from 512 total points (1.00x reduction)" — correctly the full mesh, since the set's common subgroup is trivial once the generic point is included) |
| 5 | reordered q-list | same 4 q-points, rows 2–4 shuffled (row 1 = q=0 reference unchanged) | per-q results (matched by q-vector, not row index) identical to case 4 to **0.00e+00 Ry**, for both `little_group` and `little_group_common` |

### 3.1 A real bug the test caught

The first run of case 4's `little_group_common` failed: all three nonzero-q
rows returned the *same* wrong energy (~3e-3 Ry off the oracle), and the
"little group" log line count was 0. Root cause was exactly the
`has_nonzero_q_gbt()` under-detection described in §2.3 — the pre-loop call
silently fell through to the ordinary point-group mesh at q=0 (29
irreducible k-points, matching the reference SCF's own count) instead of the
declared 4-point common-subgroup mesh, and that stale q=0 mesh was then used
for every probe. `is_declared_nonzero_q_gbt` fixed it; re-run gave the
0.00e+00 Ry agreement in the table above. This is the concrete demonstration
of why the task's "cache keys must include ... q ... symmetry policy"
requirement matters: without a correct key, a sweep with an all-zero
reference row is exactly the case where a naive "is q nonzero" check fails
silently.

## 4. Regression evidence

| suite | result |
| --- | --- |
| `ctest -L unit` | **22/22** (unchanged from the pre-WP8 baseline) |
| `WP8LittleGroup` | **PASS** (5 scenarios, 24 sub-checks, all 0.00e+00 Ry against the oracle) |
| `Example_bulk_bccFe_nsp4_block[_hoh]`, `Example_bulk_bccFe_nsp4_chebyshev[_hoh]`, `Example_bulk_bccFe_nsp4_block_spiral_q{plus,minus}`, `Example_frozen_magnon_bccFe[_auto][_auto_scf]`, `Example_k_space_scf_bccFe`, `Example_band_structure_bccFe` | see below |

The full `regression`/`triad`/`backend` ctest suite (Lanczos/Chebyshev/MKL
recursion regressions, exchange triads) was not re-run to completion this
session: those suites exercise real-space recursion and exchange machinery
this task did not touch (the diff is confined to
`source/reciprocal*.f90`, `source/spglib_interface.f90`, and the
`frozen_magnon` sweep functions in `source/calculation.f90`), an initial
attempt hit the same environment-contention timeouts already documented in
the WP7 report (batched long serial runs on this machine), and re-running a
~2.3-hour suite unrelated to the changed code is outside this task's "no
marathons" budget. `ctest -L unit` (fast, deterministic, exercises the
shared reciprocal/GBT/spin-density unit oracles) was re-run clean instead
and is unchanged.

The `scf`-suite examples that *do* exercise the changed code — every
`nsp=4`/k-space/GBT/frozen-magnon case in `tests/scf/cases.json` — were run
directly:

| case | result |
| --- | --- |
| `Example_bulk_bccFe_nsp4_block[_hoh]` | **PASS** (unchanged; `has_nonzero_q_gbt()` false, WP8 code paths not exercised) |
| `Example_bulk_bccFe_nsp4_chebyshev[_hoh]` | **PASS** (unchanged) |
| `Example_k_space_scf_bccFe` | **PASS** (unchanged) |
| `Example_band_structure_bccFe` | **PASS** (unchanged; k-path branch, explicitly untouched by `force_full_bz_for_nonzero_q_gbt`, §2.4) |
| `Example_bulk_bccFe_nsp4_block_spiral_q{plus,minus}`, `Example_frozen_magnon_bccFe[_auto][_auto_scf]` | **FAIL, pre-existing** — identical error messages (`gbt_single_q requires strux_backend='strux_lib'`, `periodic_nc` with `q_ss` set) to those already recorded in `tests/KNOWN_ISSUES.md` ("Five `tests/scf/cases.json` GBT fixtures fatal on `strux_backend='legacy'`", "`nsp4_block_spiral_qplus`/`_qminus` never enable GBT"). These decks never reach `has_nonzero_q_gbt()`; WP8 cannot have changed their failure, and did not — same file, same line, same message. |

6/6 cases capable of running before WP8 still pass after it; the 5 already
broken for unrelated, documented reasons fail identically.

## 5. Completion checklist

- [x] Full-BZ oracle remains available — `q_symmetry_policy='full_bz'`
      (default), bit-identical to the pre-WP8 code path at every dispatch
      point touched.
- [x] Little-group cache includes q and all mesh state — lattice, mesh
      dims, offset, declared q-set, and symmetry policy (§2.3).
- [x] Multi-q sweeps rebuild or use a common subgroup —
      `little_group`/`little_group_common` in
      `post_processing_frozen_magnon_acoustic`/`_auto` (§2.5).
- [x] All q-dependent solver state is invalidated —
      `ensure_kpoint_mesh` calls `invalidate_spectral_cache()` on every
      rebuild; `force_full_bz_for_nonzero_q_gbt` and
      `build_kspace_hamiltonian` already invalidate on every operator
      generation change (pre-existing WP5/WP7 machinery, confirmed still
      wired through the new dispatch paths).
- [x] Axial/generic/shifted-mesh tests pass — §3, cases 1–3.
- [x] q-list reordering leaves results unchanged — §3, case 5.
- [x] Full/reduced comparisons and G8 PASS/FAIL are reported — §3–4,
      **PASS**.

## 6. Files changed

| file | change |
| --- | --- |
| `source/spglib_interface.f90` | `get_rotation_matrices`, `get_little_group_kpoint_mesh_with_points` (multi-q), `spg_get_stabilized_reciprocal_mesh` C interface |
| `source/reciprocal.f90` | `q_symmetry_policy` + mesh-cache-key type fields, `generate_little_group_kpoint_mesh`/`ensure_kpoint_mesh` interfaces |
| `source/reciprocal_bands.f90` | `generate_little_group_kpoint_mesh`, `ensure_kpoint_mesh`, `is_declared_nonzero_q_gbt`; `generate_reduced_kpoint_mesh`'s nonzero-q dispatch made policy-aware |
| `source/reciprocal_lifecycle.f90` | `q_symmetry_policy` defaults/namelist read/validation; `force_full_bz_for_nonzero_q_gbt` made policy-aware (the real per-build gate, §2.4) |
| `source/include_codes/namelists/reciprocal.f90` | `q_symmetry_policy` namelist key |
| `source/calculation.f90` | frozen-magnon acoustic/auto sweeps: per-q rebuild or common-subgroup mesh instead of row-1-only reuse |
| `CMakeLists.txt` | `WP8LittleGroup` ctest registration |
| `tests/regression/run_wp8_littlegroup.py`, `tests/regression/wp8_littlegroup/` | new — WP8 known-answer tests |

## 7. Unresolved risks / next task

* `little_group_common`'s "keep the cached mesh untouched" check in
  `force_full_bz_for_nonzero_q_gbt` (§2.4) only compares mesh dimensions,
  not the full cache key, when deciding whether a pre-built sweep mesh is
  still valid to reuse at the per-build gate. In every sweep call site this
  branch is reached from, the mesh was just built by the same
  `ensure_kpoint_mesh(..., q_list_cart=...)` call one line above, so the
  full key cannot yet have drifted; a hypothetical future caller that
  changes lattice/offset/policy between the sweep's pre-loop build and its
  first probe would not be caught here (it would still be caught the next
  time anything calls `ensure_kpoint_mesh` with the new key, just not by
  this specific gate). Low risk, not exercised by any current call site;
  noted for whoever adds the next `little_group_common` consumer.
- The full `regression`/`triad`/`backend` ctest suite was not re-run to
  completion this session (§4) — this task's diff cannot plausibly touch
  real-space recursion/exchange code, but that is an argument from code
  boundaries, not a measurement. Recommended before the next task that
  touches `source/reciprocal*.f90` or `source/calculation.f90`: run it once,
  idle, outside a contended window.
* WP7's recorded `Example_bulk_bccFe_nsp4_block_spiral_qplus` axis failure
  and the 5 stale-`strux_backend` `scf`-suite cases are pre-existing and
  untouched by WP8; not re-triaged here.

Next allowed task: **WP9** (physical validation).
