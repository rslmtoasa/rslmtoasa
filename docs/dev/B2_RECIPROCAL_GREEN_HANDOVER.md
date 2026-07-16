# B2 `reciprocal_green` — session handover (2026-07-15)

Branch: `fable_v2`. Start a fresh session from here.

## LATEST (B2.6 landed 2026-07-16) — gate G-B2-2 awaits Anders

**B2.6 is done (code); gate G-B2-2 is Anders' to sign.** The reported "crash"
(`post_processing='exchange'` + `gf_route='lehmann'` crashing in
`calculate_exchange`) **does not reproduce** on the current tree: exchange runs
**unchanged** on the k-space-filled arrays (`gij/gji` + torque families
`Ginmag`/`Gi{x,y,z}`) and returns a physical `J_ij` (isotropic J, DMI ≈ 1e-16 on
collinear bcc Fe). No source fix to the exchange path was needed — only the stale
B2.5 comment on the `lehmann`/`dyson` branch was updated.

**The open `1/√2` intersite normalization is RESOLVED — it is correct.** The
factor-≈2 one sees comparing a fixed-broadening Lehmann `J` against the
(near-sharp) recursion `J` is a **broadening / metallic-Fermi-surface k-convergence
artifact, not a normalization error**: it is **shell-dependent**
(`J_rec/J_leh(η=0.02)` = 1.98 on 1st-NN, 1.31 on 2nd-NN) and **η-dependent**
(converged `J_leh` swings +0.41 → −0.17 as η goes 0.005 → 0.08, and climbs back
toward the recursion 0.51 as η → 0). No global factor does that. The kernel is
independently pinned at machine precision by `UnitGammaSupercell` (intersite block
≡ direct resolvent, `<1e-12`) and the recursion 4-phase algebra reduces to `G_ij`
exactly (`gij = ½[(g0₁−g0₂)+(1/i)(g0₃−g0₄)] = G_ij`, unit-norm starts). Full J vs
N_k / η study + recommended default meshes (16³ @ η=0.02 for ~1%):
`docs/dev/reciprocal_green_convergence.md`. Regression **10/10 bit-identical**;
`UnitGammaSupercell`/`UnitLehmannChain`/`UnitDysonEquivalence` all pass.
**Remaining:** Anders signs default meshes/accuracy (gate G-B2-2); the damping
eta-route cross-check on the same filled `gij_eta` is a natural follow-up.

## PRIOR (B2.5 landed 2026-07-15)

**B2.5 is done.** The `gf_route = recursion | lehmann | dyson` dispatch key was
added to the `&calculation` namelist (`calculation` type field +
`restore_to_default`='recursion' + `check_gf_route` validation) and wired into
`post_processing_exchange`: default `recursion` is byte-identical (regression
10/10), `lehmann`/`dyson` fill the SAME canonical arrays via
`reciprocal%fill_green`, reading mesh/eta from `&reciprocal` and overriding its
`green_backend`. `reciprocal_obj` is kept at **subroutine scope** — a first cut
used a nested helper and the reciprocal (spglib-owning) finalizer **double-freed
on helper return**; the `post_processing_kspace_green` pattern (subroutine-scope
`reciprocal_obj`) is the proven-safe one.

The Γ-only supercell identity + a two-route DOS cross-check landed as a new
dependency-free unit test `tests/unit/test_gamma_supercell.f90` (ctest
`UnitGammaSupercell`): periodic 1-band ring, (1) Γ-only supercell Lehmann ≡
cluster resolvent `[zI−H]⁻¹` 1.3e-15, (2) Γ-only supercell ≡ primitive full-BZ
Lehmann (folding) 1.7e-15, (3) two-route DOS `−1/π Im Tr G` (Lehmann vs Dyson)
6.2e-15 + weight→Nsc (C4). The **real-material** bcc-Fe two-route DOS ships as
the `kspace_green` driver (report-only, weight→18); its elementwise tolerance
stays under gate **G-B2-2**.

**B2.6 note (the B2.5-era "live bug" is closed):** the B2.5 handover flagged
`post_processing='exchange'` + `gf_route='lehmann'` as crashing in
`calculate_exchange` on the intersite torque families. It **does not reproduce** on
the current tree — exchange runs unchanged and returns physical J_ij (see the
LATEST section above and `docs/dev/reciprocal_green_convergence.md`). The
`1/√2` intersite normalization is verified correct (broadening/k-convergence, not a
factor). Regression stayed 10/10 bit-identical; `UnitGammaSupercell`,
`UnitLehmannChain`, `UnitDysonEquivalence` all pass; `kspace_green` driver still
EXIT=0.

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

**B2.4 landed this session** (backend D / direct Dyson + the permanent Σ=0 ≡ E
invariant). New `source/dyson_kernel.f90` (`dyson_kspace_inverse`, dependency-free
LAPACK inverse), `fill_green_dyson` in `reciprocal_green.f90` dispatched from the
formerly-stub `case('dyson')`, shared machinery factored out (`pair_geometry`,
`build_fermi_eta_contour`) so the Lehmann filler stays byte-identical. **S=I is
pinned** (backend E's `zheev` is orthonormal ⇒ D≡E holds only for S=I). Unit
invariant `UnitDysonEquivalence` passes at ~1e-16; end-to-end D≡E on real H(k)
(wired into the `kspace_green` driver) is 4.7e-11 (bcc-Fe) / 3.8e-13 (Mn3Sn NC).
Regression **10/10 bit-identical**. The next task is **B2.5** (finish the
`gf_route` namelist dispatch + the two-route DOS regression case + Γ-only
supercell identity). See the B2.4 section below.

### B2.4 landed this session

- `source/dyson_kernel.f90` — module `dyson_kernel_mod`, public
  `dyson_kspace_inverse(hk, z, sigma_full, gk)`: forms `A = z·I − H(k) − Σ` and
  returns `A⁻¹` via `zgetrf`/`zgetri` (same factor/inverse pair the RS route uses
  in `green.f90`). Dependency-free (only `precision_mod` + LAPACK) so the
  equivalence test needs no LMTO machinery. **S=I is pinned in the module header:**
  backend E diagonalizes the STANDARD Hermitian eigenproblem (`diagonalize_hamiltonian`
  → `zheev`, not `zhegv`), i.e. it assumes an orthonormal basis; the D≡E invariant
  therefore holds only when D also uses S=I. Generalized S(k)
  (`reciprocal_mode='generalized_overlap_proxy'`, `sk_overlap`) is deliberately
  **out of scope** — adding it here would break the invariant unless backend E is
  first re-cast as a generalized eigenproblem. Do not add S(k) without also
  generalizing E.
- `reciprocal_green.f90::fill_green_dyson` — the backend-D filler, dispatched from
  `fill_green`'s `case('dyson')` (previously a not-implemented stub). **Streams
  one nb·nsite × nb·nsite inversion per (k,z)** and distributes the sub-blocks to
  **every** pair (design §1.4 — never re-invert per pair, never materialize all
  (k,z)). Reuses the **exact** pair→site map, bond vector, `e^{i2π k·dR_ij}`
  inverse-Bloch phase, 1/N_k factor, and `pauli_decompose_block` torque step as
  backend E. Σ enters block-diagonal via the provider (`build_sigma_full` places
  each site's `get_sigma(z,isite,·)` block on the diagonal); `sigma_zero` ⇒ Σ=0
  ⇒ backend E. Fills gij/gji + the eta ladder + torque families exactly as the
  Lehmann filler.
- **Shared machinery factored out** (handover ask): `pair_geometry`
  (pair→ioff/joff/dR) and `build_fermi_eta_contour` (the 64-pt Gauss-Legendre
  Fermi ladder). `fill_green_lehmann` now calls both — pure code motion, its
  numeric output is **byte-identical** (both example C1 data tables unchanged).
- **Acceptance test (permanent CI):** `tests/unit/test_dyson_equivalence.f90`
  (ctest `UnitDysonEquivalence`, `RUN_UNIT_TESTS=ON`). Four pins, all ~1e-16:
  (1) 1-band chain on-site Dyson vs closed form; (2) intersite m=2 phase vs
  closed form; (3) small **multiband** H(k) (2 sites × 2 orb, nmat=4) intersite
  block with a nonzero bond phase, **Dyson route ≡ Lehmann route** (eigenpairs
  from `zheev`); (4) a **Σ=s₀·I sign pin** — Dyson at `z` must equal Lehmann at
  `z−s₀` (`A = z·I − H − Σ`), so a sign error in how Σ enters is caught.
- **End-to-end D≡E on real H(k)** wired into the `post_processing='kspace_green'`
  driver (`calculation.f90`): after the Lehmann fill it re-fills with
  `green_backend='dyson'` and reports `max|gij_dyson − gij_lehmann|` over the FULL
  gij (all pairs → intersite phase exercised). bcc-Fe = **4.7e-11**, Mn3Sn NC =
  **3.8e-13**. Report-only (`kspace_green_c1.dat`, `out.kgf`; both untracked).

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

- `feat(reciprocal_green): B2.4 backend D (direct Dyson) + Σ=0 ≡ E invariant` —
  `dyson_kernel.f90`, `fill_green_dyson` + shared `pair_geometry`/
  `build_fermi_eta_contour`, `UnitDysonEquivalence`, driver end-to-end D≡E check.

Prior-session B2.3 commits: `fe870d3` (eta ladder + torque families + the
initial, later-reverted C2 rotation), `7e9ac2c` (C2 — RS intersite `gij` is
global-frame, drop the rotation; Mn3Sn NC example + `m_z` metric).

Earlier B2.1/B2.2 commits (prior sessions): `5e78e9c` (B2.1 skeleton +
`sigma_provider` + contour), `5809022` (`green_eta`/`green_backend` namelist +
G-B2-1), `def7d8b` (B2.2 Lehmann kernel + wiring), `dff1369` (chebyshev on-site
×0.5 fix), `a61748d` (B2.5 kspace_green C1/C3 driver + bccFe example).

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

## Next task: B2.5 [SONNET] — finish namelist dispatch + regression cases

Backends **E (Lehmann) and D (Dyson) are both done**, plus the C1/C3/C2 pins,
the eta ladder + torque families, and the permanent D≡E invariant. What remains
of B2.5 (the driver + example already landed; see the checklist `[~]`):

- **`gf_route = recursion | lehmann | dyson`** control key in the post-processing
  dispatch (`calculation.f90`) so a production SCF/postproc run can select the
  k-space filler instead of the recursion route (today only the `kspace_green`
  validation driver calls `fill_green`). `green_backend='lehmann'|'dyson'` already
  exists on the `reciprocal` type + namelist; the missing piece is routing a real
  consumer (exchange/DOS) through `fill_green` behind the new key.
- **Two-route DOS regression case:** bcc-Fe k-space DOS vs recursion DOS as a
  permanent cross-validation (the maintainer-gated tolerance from B2.6 applies).
- **Γ-only supercell identity:** Γ-only Lehmann on an N-site supercell ≡ the
  real-space cluster result (blueprint validation item 2 — a self-contained
  identity worth having as a test).

Then **B2.6** [SONNET/OPUS]: run `exchange.f90` + damping on the Lehmann/Dyson-
filled arrays with **zero consumer changes** (that is the acceptance), plus the
J_ij(R) vs N_k convergence study; **gate G-B2-2** (Anders signs default meshes +
documented accuracy). Doc: `docs/dev/reciprocal_green_convergence.md`.

**Open thread from B2.2 (still not closed):** re-verify the `1/√2` normalization
when extending on-site → true intersite `G_ij` (i≠j) on the recursion side — see
the `chebyshev_recur_ij` note below. Backend E/D fill i≠j directly (no 4-phase),
so this is a recursion-route concern, not a k-space one.

**Overlap note carried forward:** backend D pins **S=I** (matches backend E's
orthonormal `zheev`). If a future task wants generalized S(k), backend E must be
re-cast as a generalized eigenproblem FIRST, then both re-validated — do not add
`sk_overlap` to `dyson_kspace_inverse` in isolation (it would silently break D≡E).

## Session kit for B2.5 (per the token-lean protocol)

`docs/dev/plans/B2_reciprocal_green.md` §2 (task B2.5); repo `CLAUDE.md`;
`source/calculation.f90` (post-processing dispatch + the `kspace_green` driver as
the existing `fill_green` caller); `source/reciprocal_green.f90` (the finished
`fill_green` / `fill_green_lehmann` / `fill_green_dyson`); an existing `&namelist`
as a template. Do **not** read the whole repo; prefer targeted grep.

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
