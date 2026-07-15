# B2 — `reciprocal_green`: k-space Green's functions, two backends (flagship)

**Effort:** M–L (3–6 weeks). **Lead:** OPUS (kernels + conventions), SONNET
(plumbing). **Depends on:** B1 (phase/bond conventions — B1 gate G-B1-1 must
be signed first). **Unlocks:** B3, B5, B8, B10, B11.

---

## 1. Design (maintainer-approved; supersedes nothing, refines BLUEPRINTS §B2)

### 1.1 The contract that makes everything else cheap

New submodule `source/reciprocal_green.f90` of the `reciprocal_*` family
(post-T9 submodule pattern). Public surface ≈ one type-bound routine:

```
call reciprocal%fill_green(green_obj, mesh, pairs, sigma_provider)
```

which populates **the same arrays the recursion route populates** on the
existing `green` object (`source/green.f90` lines ~62–64):
`gij/gji(nb, nb, nE, njij)` on the `en%ene` grid, the `gij_eta(64, nb, nb,
·)` Fermi-point ladder, and the torque-resolved components (`ginmag`,
`gix/giy/giz` and j-side partners) — using the same `njij` pair bookkeeping
from `lattice`. Consumers (`exchange.f90` via `auxiliary_gij`, damping,
future χ) then run **unchanged, by construction**. There is no adapter
layer, no `green_source` abstraction: filling the canonical arrays IS the
route-agnosticism. (This is binding; it is recorded in
`DEVELOPER_MAP.md`.)

### 1.2 Two backends behind the one filler

**Backend E — eigenpair / strict Lehmann (Σ = 0).**
```
G_ij(z) = (1/N_k) Σ_{k,n} e^{ik·ΔR_ij} ψ_{i,n}(k) ψ†_{j,n}(k) / (z − ε_{nk})
```
One `zheev` per k (reuse `reciprocal_bands.f90` machinery; B4 later swaps
in batched GPU zheevd transparently), then all energies and all pairs
amortized over the eigenpairs. Broadening: explicit iδ, with tetrahedron
weights as a later option (host-side, out of B4's scope).

**Backend D — direct Dyson inversion (Σ(z) ≠ 0).**
```
G(k,z) = [ z·S − H(k) − Σ(z) ]⁻¹        per (k, z)
```
batched nb×nb `zgetrf/zgetrs`; same accumulation into the same arrays.
Σ enters via the `sigma_provider` abstract interface defined here
(first real implementation: B8's `sigma_cpa`; `sigma_zero` default now;
`sigma_static` from DFT+U as the cheap plumbing test — see B10).

**Permanent CI invariant:** backend D with Σ = 0 ≡ backend E to solver
tolerance, every run of the unit suite.

### 1.3 Convention checklist — this is the actual work

Each item is a task with its own known-answer test. Conventions, not math,
are where this milestone can silently fail.

| # | Convention | Pin |
|---|---|---|
| C1 | **LMTO representation:** deliver G in the same (screened/auxiliary) representation the RS route delivers *pre-*`auxiliary_gij`; downstream transforms then apply identically | one on-site block vs RS route at large broadening, elementwise |
| C2 | **Local spin frames:** rotate global-frame eigenvector products into the per-site local frames used by the noncollinear RS blocks (`local_axis` machinery in `hamiltonian_build`) | noncollinear test case, G_ii spin structure vs RS |
| C3 | **Phase/bond convention:** e^{ik·ΔR} with the true bond vector — adopt **verbatim** the convention pinned in B1.3's doc-test | reuse B1.3's asserted α factor; same signs |
| C4 | **Normalization:** 1/N_k, spin-degeneracy per `nsp`, energy zero (physical `fermi`, never the `chebfermi`-scaled variable) | DOS sum rule ∫DOS dE = electron count |

### 1.4 Efficiency spec

- Never materialize G(k,z) for all (k,z): backend E streams over k
  accumulating into the (pair, z) targets; backend D streams (k, z) tiles.
- Diagonalize once per k; reuse eigenpairs for all energies, all pairs,
  and the `gij_eta` ladder in the same pass.
- MPI: distribute k (backend E) or (k,z) tiles (backend D); reduce into
  the green arrays with the same pattern the recursion route uses for
  `atoms_per_process` (see `green.f90` allocation, lines ~200–290).
- fp64 (`rp`) throughout; contour/η owned by `energy.f90` — never
  hard-coded.

## 2. Tasks

### B2.1 [OPUS] — Module skeleton + `sigma_provider` + contour adoption
Type, constructor, `restore_to_default`, abstract `sigma_provider` +
`sigma_zero`, energy-grid adoption from the `en` object. Doc comments carry
§1 equations. **Gate G-B2-1:** Anders signs that the adopted contour/grid
convention matches what `exchange.f90`'s contour loop assumes, before
B2.2. *Kit:* this file; `green.f90` lines 1–320; `energy.f90` grid
accessors; `reciprocal.f90` type region.

### B2.2 [OPUS] — Backend E core + C1/C3/C4 pins
Lehmann fill of `gij/gji` + on-site blocks. Unit test: hand-coded 1-band
chain H(k) = −2t cos k vs closed-form G(z) to 1e-10; C1 test (on-site block
vs RS at large broadening, bcc Fe); C4 sum rule. *Kit:* B2.1;
`reciprocal_bands.f90` lines 1–240; B1.3 convention test file; this file
§1.2–1.3.

### B2.3 [OPUS] — C2 local frames + `gij_eta` ladder + torque components
Fill the eta ladder and `ginmag`/`gi{x,y,z}` families exactly as the RS
route defines them (read the recursion-route filler as the normative
reference — locate via `DEVELOPER_MAP.md`, do not guess index order).
Noncollinear C2 test. *Kit:* B2.2; the RS gij_eta filler routine;
`hamiltonian_build.f90` local_axis region.

### B2.4 [OPUS] — Backend D + equivalence invariant
Per-(k,z) `zgetrf/zgetrs` with `sigma_provider`; the Σ=0 ≡ backend E test
(chain fixture + one bcc-Fe H(k), 1e-9 elementwise) wired as a permanent
unit test. *Kit:* B2.1–B2.3; this file §1.2.

### B2.5 [SONNET] — Namelist + dispatch + regression cases
`gf_route = recursion | lehmann | dyson` key in post-processing control;
`&reciprocal_green` namelist (mesh, η, backend); regression case: bcc Fe
k-space DOS vs recursion DOS (the permanent two-route cross-validation);
Γ-only Lehmann on an N-site supercell ≡ real-space cluster result
(blueprint validation item 2 — a beautiful self-contained identity). *Kit:*
B2.1–B2.4 interfaces; `calculation.f90` post-processing dispatch; an
existing namelist as template.

### B2.6 [SONNET] — J_ij + damping through the filled arrays + convergence
Run `exchange.f90` and the damping eta-route on Lehmann-filled arrays with
**zero consumer changes** (that's the acceptance); J_ij(R) vs N_k
convergence study for bcc Fe (the known Lehmann-route disappointment mode —
Friedel tails at metallic E_F). **Gate G-B2-2:** Anders signs default
meshes + documented accuracy. Doc:
`docs/dev/reciprocal_green_convergence.md`. *Kit:* B2.5; `exchange.f90`
lines 260–300 (consumer shape only).

## 3. Progress checklist
- [~] B2.1 skeleton + provider — **code landed; G-B2-1 OPEN** (see
      `docs/dev/B2_GATE_G-B2-1.md`). `sigma_provider.f90`,
      `reciprocal_green.f90` (contour + dispatcher), type fields + defaults,
      CMake registration; regression 10/10 bit-identical.
- [~] B2.2 backend E + C1/C3/C4 — **kernel + wiring landed; pure-math pins
      PASS to machine precision.** `source/lehmann_kernel.f90`
      (`lehmann_pair_block`, dependency-free numerical core);
      `reciprocal_green.f90::fill_green_lehmann` wires eigenpairs + pair→site
      map + bond vectors and fills `gij/gji` (on-site = i==j pair).
      `tests/unit/test_lehmann_chain.f90` (ctest `UnitLehmannChain`, gate
      `RUN_UNIT_TESTS=ON`) pins the 1-band chain on-site closed form (8.7e-16),
      the intersite `e^{ik·ΔR}` phase/C3 (4.2e-16), and the 1/N_k
      normalization/C4 (2e-12). Regression 10/10 bit-identical (fill path not
      production-wired). **LMTO-integration C1/C3 DONE** via a production
      `post_processing='kspace_green'` driver (`calculation.f90`) + example
      `example/exchange/bccFe_kspace_green`: runs BOTH routes on the same
      converged bcc-Fe potential and reports the on-site DOS
      `-1/pi Im Tr G_ii(E)` from the same `green%gij` on-site block (DOS-level
      cross-check, maintainer's choice). Both routes -> C4 sum rule weight ≈ nb=18
      and agree in shape; report-only (`kspace_green_c1.dat`), tolerance is a
      maintainer gate. **Side bug fixed:** `chebyshev_recur_ij` lacked the
      block-path on-site (i==j) special case -> `G_ii` came out ×0.5 (start
      `1/sqrt2|i>`, norm²=½, moments quadratic). Fixed in both CPU+GPU branches;
      latent for i≠j J_ij, no golden affected. **FOLLOW-UP:** re-verify this
      normalization when extending on-site → true intersite `G_ij` (i≠j).
- [x] B2.3 eta ladder + torque + C2 (spin frame) — **DONE, validated on a
      genuine NC background (Mn3Sn, 120 deg).**
      `fill_green_lehmann` now also fills the 64-point `gij_eta`/`gji_eta` Fermi
      ladder (Lehmann blocks at `z = ene(fermi_point) + i*(1-x)/x`, the same
      `bgreen` eta contour + `gauss_legendre(64,0,1)` roots the recursion route
      uses) and the torque-resolved families `ginmag`/`gi{x,y,z}` (+ `gj*` and
      the `*_eta` partners). The Pauli (charge,x,y,z) spin-block decomposition
      was extracted to `lehmann_kernel_mod::pauli_decompose_block` (dependency-
      free, reused for both energy-grid and eta blocks) and pinned by a new
      `UnitLehmannChain` Test 4 (transcription guard, 1.1e-16).
      **C2 RESOLVED — the plan's premise was wrong.** The intersite recursion
      `recur_b_ij` (recursion_transport.f90) **never** rotates to local axes
      (`local_axis` there only gates GPU eligibility; only the on-site DOS
      recursion `recursion_haydock` rotates), so the RS `gij`/`gji` arrays are
      stored in the **GLOBAL** spin frame. An early implementation that rotated
      only the k-space block into the local frame BROKE the match on Mn3Sn
      (m_z diff 4e-4 -> ~20). Backend E therefore fills the global-frame block
      directly (no rotation); the driver's on-site z-projected spin DOS m_z
      agrees between routes to **4.2e-4** on an off-axis Mn site (moment
      (-0.5,0.866,0), both m_z ~ 0 as expected in-plane). Validation via a new
      NC example `example/exchange/Mn3Sn` (8³ mesh, local_axis off) + the extended
      `post_processing='kspace_green'` driver (now also reports m_z and the
      on-site block diff; `rotate_from_local_axis` restores global ee as a safety
      net). A future LOCAL-frame comparison would have to rotate BOTH routes'
      GF (pointer comment left in `reciprocal_green.f90`); the intersite i/=j
      frame (i, j moments differ) is an open question, deliberately out of scope.
      Regression 10/10 bit-identical; collinear `kspace_green` C1 byte-identical.
- [x] B2.4 backend D (Dyson) + Σ=0 invariant — **DONE.**
      `source/dyson_kernel.f90` (`dyson_kspace_inverse`, dependency-free LAPACK
      `zgetrf`/`zgetri` core; **S=I pinned** — backend E's `zheev` is orthonormal,
      so D≡E holds only for S=I; generalized S(k) deferred). `reciprocal_green.f90`
      gained `fill_green_dyson`, dispatched from `fill_green`'s `case('dyson')`
      (was a stub): streams **one nmat×nmat inversion per (k,z)** and distributes
      the sub-blocks to **every** pair with the same 1/N_k inverse-Bloch phase +
      pair→site map + `pauli_decompose_block` torque step as backend E (§1.4 — no
      re-inversion per pair, no full-(k,z) materialization). Σ enters block-diagonal
      via the provider (`build_sigma_full`); `sigma_zero` ⇒ backend E. Shared
      machinery factored out (`pair_geometry`, `build_fermi_eta_contour`) — the
      Lehmann filler now calls them, byte-identical output (both example data
      tables unchanged). **Acceptance (permanent CI):** `UnitDysonEquivalence`
      (ctest, `RUN_UNIT_TESTS=ON`): 1-band chain Dyson vs closed form (on-site
      7.8e-16, intersite 1.8e-16), a small multiband H(k) D≡E block with a nonzero
      bond phase (6.7e-16), and a Σ=s₀·I sign pin (9.9e-16). **End-to-end** D≡E on
      real H(k) wired into the `kspace_green` driver: bcc-Fe `max|gij_dyson −
      gij_lehmann|` = 4.7e-11 (4096 k, all pairs), Mn3Sn NC = 3.8e-13. Regression
      10/10 bit-identical.
- [x] B2.5 dispatch + DOS regression + Γ-only identity — **DONE.**
      **`gf_route = recursion | lehmann | dyson`** key added to the `&calculation`
      namelist (`calculation` type field + `restore_to_default`='recursion' +
      `check_gf_route` validation). Wired into `post_processing_exchange`: the
      default `recursion` reproduces the legacy path byte-for-byte (regression
      10/10 bit-identical), while `lehmann`/`dyson` fill the SAME canonical arrays
      from `reciprocal%fill_green` (backend E/D), reading mesh/eta from the
      `&reciprocal` namelist and overriding its `green_backend`. `reciprocal_obj`
      is kept at subroutine scope so its spglib-owning finalizer runs exactly as
      in `post_processing_kspace_green` (a nested helper double-freed on return).
      **Γ-only supercell identity + two-route DOS** landed as a permanent,
      dependency-free unit test `tests/unit/test_gamma_supercell.f90` (ctest
      `UnitGammaSupercell`, `RUN_UNIT_TESTS=ON`): a periodic 1-band ring proves
      (1) Γ-only supercell Lehmann ≡ cluster resolvent `[zI−H_sc]⁻¹`
      (`dyson_kspace_inverse`) 1.3e-15; (2) Γ-only supercell ≡ primitive-cell
      full-BZ Lehmann (the folding identity) 1.7e-15; (3) two-route DOS
      `−1/π Im Tr G` (Lehmann site-summed vs Dyson trace) agree 6.2e-15 with the
      integrated weight → Nsc (C4). The **real-material bcc-Fe** two-route DOS
      cross-validation already ships as the `kspace_green` driver (report-only,
      weight → nb=18); its *elementwise* acceptance tolerance stays under gate
      **G-B2-2** (B2.6). **B2.6 entry point:** running `post_processing='exchange'`
      with `gf_route='lehmann'` reaches `fill_green` cleanly, then
      `calculate_exchange` crashes reading the intersite torque families
      (Ginmag/Gj{x,y,z}) — arrays the DOS driver never exercised. Producing
      correct J_ij on the k-space arrays (and re-verifying the 1/√2 intersite
      normalization) is B2.6's acceptance, not this dispatch key's.
- [ ] B2.6 J_ij/damping zero-change run + convergence (G-B2-2 signed) — **START
      HERE:** `gf_route='lehmann'`/`'dyson'` on `post_processing='exchange'`
      currently crashes in `calculate_exchange` on the k-space-filled torque
      families (see B2.5 entry). Fix the consumer/fill so exchange runs unchanged,
      then the J_ij(R) vs N_k convergence study + gate G-B2-2.
