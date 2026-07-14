# B2 `reciprocal_green` — session handover (2026-07-13)

Branch: `fable_v2`. Start a fresh session from here.

## TL;DR

Milestone **B2** (k-space Green's-function engine, flagship) is underway. **B2.1
is done; gate G-B2-1 is signed. B2.2 backend-E kernel + wiring landed** with the
pure-math pins at machine precision. **The LMTO-integration C1/C3 now landed too**
as a production `post_processing='kspace_green'` validation driver (DOS-level
cross-check, maintainer's choice) + a dedicated example
(`example/exchange/bccFe_kspace_green`). It runs BOTH routes on the same
converged bcc-Fe potential and reports the on-site DOS
`rho(E) = -1/pi Im Tr G_ii(E)` from the same `green%gij` on-site block. **Result:
both routes satisfy the C4 sum rule (weight -> nb=18) and agree in shape;** the
residual is k-mesh/broadening ripple (B2.6 convergence territory). Report-only —
the acceptance tolerance is a maintainer gate.

**Bug found + fixed this session (via the driver):** the chebyshev intersite
recursion `chebyshev_recur_ij` was **missing the on-site (i==j) special case** the
block path `recur_b_ij` has, so an on-site self-pair started from `1/sqrt2|i>`
(norm^2=1/2) and — moments being quadratic in psi0 — delivered `G_ii` at **exactly
half amplitude**. Fixed in both CPU and GPU branches
(`source/recursion_transport.f90`) by mirroring the block guard. Latent for normal
J_ij (i!=j) runs; no existing golden uses chebyshev intersite/on-site pairs, so
regression stays bit-identical. **FOLLOW-UP (maintainer note): re-check this
1/sqrt2 normalization factor once we proceed from G_ii to the true intersite
G_ij** (i!=j) — the 4-phase combination there must reproduce the correct amplitude.

**B2.3 landed this session** (eta ladder + torque families; C2 resolved). The C2
"local frame" premise turned out to be wrong — the RS intersite `gij` is stored
in the GLOBAL spin frame, so backend E fills it directly (no rotation), validated
on a genuine NC background (Mn3Sn, m_z agrees to 4.2e-4). Regression stays
**10/10 bit-identical**; the collinear `kspace_green` C1 output is byte-identical.
The next task is **B2.4** (backend D + the Σ=0 ≡ backend E invariant). See the
B2.3 section below.

### B2.3 landed this session

- `source/lehmann_kernel.f90` — new public `pauli_decompose_block`: the
  dependency-free (charge,x,y,z) spin-block decomposition
  `gnmag=½(G_uu+G_dd)`, `gz=½(G_uu−G_dd)`, `gy=½(iG_ud−iG_du)`, `gx=½(G_ud+G_du)`,
  reused for BOTH the energy-grid blocks and the Fermi-eta ladder. Pinned by a new
  `UnitLehmannChain` Test 4 (transcription guard on signs + i-factor, 1.1e-16).
- `reciprocal_green.f90::fill_green_lehmann` now also fills:
  - **`gij_eta`/`gji_eta`** (the 64-point Fermi ladder) as strict-Lehmann blocks
    at `z = ene(fermi_point) + i·(1−x)/x` — the same `bgreen` eta contour
    (`z=e(ei)+eta`) and `gauss_legendre(64,0,1)` roots the recursion route uses.
    Stored eta-index-leading `(64,nb,nb,pair)`; only when the arrays are
    allocated. `wgl` weights belong to the consumer, not the fill.
  - **Torque families** `ginmag`/`gi{x,y,z}` (+`gj*` and the `*_eta` partners) via
    `pauli_decompose_block`. `store_eta_torque` transposes the eta axis to front.
  - **C2 spin frame — RESOLVED, and the plan's premise was WRONG.** The
    intersite recursion `recur_b_ij` (recursion_transport.f90) **never** rotates
    to local axes — `local_axis` there only gates GPU eligibility; ONLY the
    on-site DOS recursion `recursion_haydock` rotates. So the RS `gij`/`gji`
    arrays (the ones the contract fills) are stored in the **GLOBAL** spin frame.
    Backend E fills the global-frame block directly (no rotation). Verified on a
    genuine NC background (new example `example/exchange/Mn3Sn`, 120° AFM): the
    driver now also reports the z-projected spin DOS `m_z`, which agrees between
    routes to **4.2e-4** on an off-axis Mn site (moment (-0.5,0.866,0), both
    m_z≈0 in-plane). An early version that rotated ONLY the k-space block into
    the local frame broke this (m_z diff → ~20). A pointer comment in
    `reciprocal_green.f90` records the future option (rotate BOTH routes' GF for
    a local-frame comparison); the intersite i≠j frame is an open question, out
    of scope. Driver also calls `rotate_from_local_axis` to restore global `ee`
    before the k-space build (safety net).

**Still open from B2.2:** re-verify the `1/√2` normalization when extending
on-site → true intersite `G_ij` (i≠j).

### B2.2 landed this session

- `source/lehmann_kernel.f90` — module `lehmann_kernel_mod`, public
  `lehmann_pair_block`: the dependency-free numerical core (only `precision_mod`
  + `math_mod`). Computes one nb×nb intersite block
  `G_ij(z)=(1/N_k)Σ_{k,n} e^{i2π k·dR_ij} ψ_i ψ_j† /(z−ε_nk)`. The module header
  carries the full sign derivation (inverse Bloch transform → `dR_ij = R_i−R_j`,
  `+` phase, reusing the B1 `ham_vec_type_direct` table — C3 verbatim).
- `reciprocal_green.f90::fill_green_lehmann` — the backend-E filler. Ensures
  eigenpairs (`build_kspace_hamiltonian`+`diagonalize_hamiltonian`), builds
  fractional `kfrac`, loops `lattice%ijpair` exactly as
  `calculate_intersite_gf_core` (start_atom..end_atom, g2l_map), maps cluster
  atom → unit-cell site via `iz`, forms `dR = a_cart_inv·(cr_i−cr_j)·alat`, and
  fills `gij`/`gji` (on-site = the i==j pair). `fill_green` now dispatches
  `'lehmann'` here.
- `tests/unit/test_lehmann_chain.f90` — ctest `UnitLehmannChain`
  (`-DRUN_UNIT_TESTS=ON`). 1-band chain `H(k)=−2t cos k`: on-site vs closed form
  `8.7e-16`; intersite m=2 phase vs residue closed form `4.2e-16`; 1/N_k
  normalization `z·G→1` `2.0e-12`. New CMake option + `add_fortran_unit_test`
  helper.

**B2.2 preconditions / deferred (documented in `fill_green_lehmann`):** backend E
currently needs the full unreduced BZ mesh with no MPI k-distribution (each rank
holds all k) and the standard orthonormal eigenproblem (Lehmann completeness).
Symmetry-star averaging + k-parallel reduction attach with B2.5/B4. The
torque-resolved families (`ginmag`/`gi{x,y,z}`) and the `gij_eta` ladder are
**B2.3** — `fill_green_lehmann` fills only `gij`/`gji` so far.

### Original B2.2 brief (for reference)

Plan file (authoritative, read first): `docs/dev/plans/B2_reciprocal_green.md`.
Design context: `BLUEPRINTS.md` §B2. Gate record: `docs/dev/B2_GATE_G-B2-1.md`.

## The B2 contract (do not violate — it is the whole point)

One filler, `reciprocal%fill_green(green_obj, sigma)`, populates **the same
arrays the recursion route fills** on the `green` object, so consumers
(`exchange` via `auxiliary_gij`, damping, future χ) run **unchanged**. Filling
the canonical arrays IS the route-agnosticism — no adapter layer, no
`green_source` abstraction. Two backends behind the one filler:

- **Backend E** (eigenpair / strict Lehmann, Σ=0):
  `G_ij(z) = (1/N_k) Σ_{k,n} e^{ik·ΔR_ij} ψ_{i,n}(k) ψ†_{j,n}(k) / (z − ε_{nk})`
- **Backend D** (direct Dyson, Σ≠0): `G(k,z) = [zS − H(k) − Σ(z)]⁻¹`, Σ via the
  `sigma_provider`. Permanent CI invariant: **D with Σ=0 ≡ E** (that is B2.4).

## Commits made this session (on `fable_v2`)

- `5e78e9c` feat(reciprocal_green): B2.1 skeleton, `sigma_provider`, contour
  adoption.
- `4312133` fix(tests): gave the three `Example_frozen_magnon_bccFe*` cases a
  `control` namelist so `RUN_EXAMPLE_TESTS=ON` configures (was a pre-existing B1
  test-data bug that blocked cmake). They now register as ctest #32–34.
- `5809022` feat(reciprocal_green): exposed `green_eta`/`green_backend` as
  `&reciprocal` namelist inputs; recorded gate G-B2-1 resolution.

## What exists now (B2.1)

- `source/sigma_provider.f90` — module `sigma_provider_mod`: abstract
  `sigma_provider` (deferred `get_sigma(z, isite, sigma_block)`) + concrete
  `sigma_zero`. Its own small module so B8/B10 providers depend only on the
  interface, never on `reciprocal_mod` (no dependency cycle).
- `source/reciprocal_green.f90` — submodule `(reciprocal_mod) reciprocal_green`:
  - `build_green_contour(this, en, z_contour)` — **concrete**, the gate-signed
    contour: `z(ie) = en%ene(ie) + i·green_eta` over the full `en%ene` grid.
  - `fill_green(this, green_obj, sigma)` — dispatcher; backend cores currently
    raise "not implemented" (E→B2.2, D→B2.4).
- `reciprocal` type gained `green_eta` (real, default 0.01 Ry) and
  `green_backend` (`'lehmann'`|`'dyson'`, default `'lehmann'`); defaults in
  `reciprocal_lifecycle.f90::restore_to_default`; namelist read in the same
  file's `build_from_file` (validated).
- Both new files registered in `source/CMakeLists.txt`.

## Conventions PINNED — do not re-derive (gate G-B2-1, signed)

1. **Energy grid:** `en%ene(1:size(en%ene))`, real-axis mesh from
   `energy%e_mesh` (`energy_min .. above fermi`, `en%channels_ldos + 10` pts).
   Call `en%e_mesh()` before filling (fill_green already does).
2. **Contour:** retarded, `z = en%ene + i·green_eta`. Matches `bgreen`
   (`z = e(ei) + eta`, `eta = i·η`), `green.f90:1555–1587`.
3. **Fermi:** physical `en%fermi`, **never** the `chebfermi`-scaled variable.
4. **Representation:** screened/auxiliary LMTO blocks, **pre-`auxiliary_gij`**
   (consumers call `green%auxiliary_gij` themselves after the fill).
5. **On-site C1 reference (maintainer decision):** `Gii` and `Gij` are **both**
   `2·norb × 2·norb` (= `nb × nb`) **full complex resolvent** blocks from block
   or Chebyshev recursion. The C1 test compares the Lehmann on-site block
   against that (real + imag) at large broadening. The `−iπρ` form stored in
   `sgreen` (`g0`, `green.f90:892`) is the DOS-derived on-site GF for
   moments/charge — **NOT** the C1 reference.

## Next task: B2.2 [OPUS] — backend E core + C1/C3/C4

Fill `gij/gji` and the on-site blocks by strict Lehmann. Acceptance tests (write
first, never loosen a tolerance):

- **1-band chain** `H(k) = −2t·cos k` vs closed-form `G(z)` to 1e-10 (pure math,
  no LMTO machinery — the cleanest unit test).
- **C1**: on-site block vs RS route at large broadening, bcc Fe, elementwise
  (against the full-resolvent block — see convention 5).
- **C3** (phase/bond): reuse B1's `e^{ik·ΔR}` bond convention **verbatim** — the
  spiral phase in `reciprocal_fourier.f90::fourier_transform_gbt_array` uses
  `2π·q·R_direct` with `ham_vec_type_direct`; the Lehmann `e^{ik·ΔR}` must use
  the same neighbor-vector table and sign.
- **C4** (normalization): 1/N_k, spin factor per `nsp`, physical `fermi`; sum
  rule ∫DOS dE = electron count.

**Normative RS references to read (do not guess index order):**
- `green.f90::calculate_intersite_gf_core` (lines ~490–604): the 4-phase `gij`
  combination `gij = 0.5·[g0(1) − g0(2) + (1/i)(g0(3) − g0(4))]` and the torque
  decomposition (lines ~586–601). Backend E fills the true intersite block
  directly (no 4-phase machinery), but the **torque components** (`ginmag`,
  `gi{x,y,z}`) are derived from `gij` by exactly that spin-block algebra — that
  part is B2.3, but read it now so B2.2's `gij` layout matches.
- `green.f90::bgreen` (from line 1461): the resolvent contour and block layout.
- Eigenpairs per k: reuse `reciprocal_bands.f90` (`zheev` machinery) — B4 later
  swaps in a batched GPU eigensolver transparently.

## Session kit for B2.2 (per the token-lean protocol)

`docs/dev/plans/B2_reciprocal_green.md` §1.2–1.3; repo `CLAUDE.md`;
`source/reciprocal_green.f90` (B2.1); `source/reciprocal_bands.f90` lines 1–240;
`green.f90` restore_to_default (arrays, lines ~195–320) and
`calculate_intersite_gf_core`; the B1 phase convention in
`reciprocal_fourier.f90`. Do **not** read the whole repo; prefer targeted grep.

## Build & test

- Build: `cd build && cmake . && make -j4`. Config now succeeds with
  `RUN_EXAMPLE_TESTS=ON` (cache has it ON). ENABLE_CUDA_PLUGIN + ENABLE_MKL_KERNELS
  are ON in this build.
- Regression (must stay 10/10, abs-tol 1e-6, feature-off bit-identical):
  `python3 tests/regression/run_matrix.py --binary build/bin/rslmto.x
  --cases-json tests/regression/cases.json --case-name <case>
  --scratch-root <scratch> --references tests/regression/references
  --abs-tol 1e-6`, looped over the 10 cases in `cases.json`.
- Run **before starting and before committing** (Phase-1 ground rule 1).

## Gotchas / conventions

- **Fenced files (do not edit):** `source/self.f90`, `source/symbolic_atom.f90`.
- `nb = 2·norb` (spd: norb=9, nb=18); BdG doubles it. `spin_off` offsets the
  down block. `i_unit`, `pi` from `math_mod`; `rp` from `precision_mod`.
- `nsp`: 1=collinear scalar-rel, 2=+SOC, 3=non-collinear scalar-rel, 4=+SOC.
- `reciprocal_mod` now `use`s `green_mod`, `energy_mod`, `sigma_provider_mod`.
  `green_mod` does **not** use `reciprocal_mod` — keep it that way (no cycle).
- One task, one session, one commit. The physics spec in the plan is
  authoritative — if it looks wrong, stop and escalate, don't "fix" it (that is
  how B1 failed the first time).

## Memory

`[[b2-reciprocal-green-state]]` and `[[gbt-frozen-magnon-state]]` in the session
memory index carry the short-form state.
