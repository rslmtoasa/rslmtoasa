# RS-LMTO-ASA — Feature Blueprints (Phase 3+)

**Audience:** developers and coding assistants planning major feature work on
the post-refactor codebase (`fable_v2` layout and later).
**Status:** design document. Nothing here is implementation instruction yet;
each blueprint is promoted to a task plan (Phase-1 style) when its turn comes.
**Standing rules:** all Phase 1/2 ground rules apply (KISS, class-based
architecture, legacy `self.f90`/`symbolic_atom.f90` fence, known-answer tests
before an estimator ships, route-agnostic consumers). Strategic identity: this
is *the* real-space code for defected/disordered/interfacial magnetism, with
k-space machinery as a cross-validating and enabling companion — not another
k-space total-energy code.

**Parked (decided, do not revisit without new arguments):** Full-charge
density (Vitos). Rationale: very large implementation whose payoff overlaps
territory where FCD-EMTO and plane-wave codes are unbeatable; the combined
correction is this code's pragmatic accuracy upgrade.

---

## Dependency DAG

```
B1 GBT fix ──────────────┐
                          ├─→ B11 LR-TDDFT (magnons)
B2 reciprocal_green ──┬──┤
   (two backends)     │   ├─→ B3 Bloch spectral functions
                      │   ├─→ B5 route-agnostic G(E) post-processing
B4 GPU k-space ───────┘   │      (Jij / damping / transport via Lehmann)
                          ├─→ B8 CPA + DLM (k-space)
                          └─→ B10 DMFT (Σ-provider API) ─→ B12 couplings
B6 surface electrostatics ─→ B7 boundary self-energies (leads + vacuum)
B7 ──→ B9 RS-CPA / RS-DLM (shares effective-medium machinery with B8)
```

Rough effort scale used below: **S** ≈ days, **M** ≈ 1–3 weeks, **L** ≈
1–3 months, **XL** ≈ 3+ months of focused effort (AI-assisted pace).

---

## B1. Generalized Bloch theorem spin spirals (bug fix + frozen magnons)

**Motivation.** The existing spiral branch is known-broken; fixing it is the
cheapest route to magnon physics (frozen-magnon J(q), adiabatic spectra) and
establishes the phase/bond-vector conventions that B2 must reuse.

**Diagnosis (already made).** `ham0m_nc` uses absolute site positions instead
of the bond vector, and applies what appears to be a full-angle rather than
half-angle SU(2) rotation. Correct GBT: for spiral vector q, the inter-site
block transforms with half-angle bond rotations,
H_ij → U(q·ΔR_ij/2) H_ij U†(q·ΔR_ij/2), with U the SU(2) z-rotation, and the
cone angle handled in the local frames.

**Insertion points.** `hamiltonian_build.f90` (spiral branch of the
noncollinear builder); no consumer changes — SCF and recursion are agnostic.

**Validation (known-answer).**
1. q = 0 reproduces the collinear/noncollinear reference exactly.
2. Flat spiral in bcc Fe at commensurate q vs. explicit supercell total
   energy: agreement to SCF tolerance.
3. E(q) symmetry: E(q) = E(−q) without SOC; smoothness through zone
   boundaries.
4. Small-q expansion vs. spin-wave stiffness from the J_ij sum rule
   (independent machinery — a true cross-check).

**Deliverables beyond the fix.** A `frozen_magnon` post-processing mode:
E(q) sweep and adiabatic magnon dispersion ω(q) from E(q) and moments. This
later becomes the validation target for B11.

**Effort/Risk:** S–M / low. **Dependencies:** none. **Do first.**

---

## B2. `reciprocal_green` — Lehmann + Dyson Green's functions (flagship)

**Motivation.** One k-space GF engine that (a) feeds the *existing*
real-space consumers untouched, (b) is the substrate for BSF, CPA, DMFT,
LR-TDDFT, and all Σ(E) physics.

**Design: one filler API, two backends.**
New submodule `reciprocal_green.f90` of `reciprocal_mod` (post-T9 pattern),
public surface ≈ one type-bound routine:
`call reciprocal%fill_green(green_obj, mesh, pairs, sigma_provider)`
which populates the **same arrays the recursion route populates** on the
`green` object: `gij/gji(nb, nb, nE, njij)` on `en%ene`, the
`gij_eta(64, ...)` Fermi-point ladder, and the torque-resolved components
(`Ginmag`, `Gi{x,y,z}`, ...), using the same pair bookkeeping (`njij`
ordering from `lattice`).

- **Backend E (eigenpair / strict Lehmann):** Σ = 0.
  G_ij(z) = (1/N_k) Σ_{k,n} e^{ik·ΔR_ij} ψ_{i,n}(k) ψ†_{j,n}(k) / (z − ε_{nk}).
  One diagonalization per k, all energies and all pairs amortized. GPU-batched
  via B4. Broadening: explicit iδ or tetrahedron weights (option).
- **Backend D (direct Dyson inversion):** Σ(z) ≠ 0.
  G(k,z) = [z S − H_k − Σ(z)]⁻¹ per (k, z); batched 18×18 (or nb×nb)
  `getrf/getrs` — trivially GPU-batched. Same accumulation into the same
  arrays. Σ enters through the provider interface (B10/B12).

**Interface contract (route-agnostic consumers).** Downstream code (exchange,
damping, transport, future χ) must not know which route or which Σ produced
the arrays. This contract is already recorded in `DEVELOPER_MAP.md`; it is
binding.

**Convention checklist — the actual work.** Each item gets its own
known-answer test before the estimator ships:
1. **LMTO representation:** deliver G in the same (screened/auxiliary)
   representation the RS route delivers pre-`auxiliary_gij`; the downstream
   transforms then apply identically. Establish by comparing a single on-site
   block against the RS route at large broadening.
2. **Local spin frames:** rotate global-frame eigenvector products into the
   per-site local frames used by the noncollinear RS blocks.
3. **Phase/bond convention:** e^{ik·ΔR} with the true bond vector — adopt
   verbatim the convention fixed in B1.
4. **Normalization:** 1/N_k, spin-degeneracy factors per `nsp`, and the
   energy zero (`fermi` vs `chebfermi` scaled variables — use physical).

**Validation (known-answer).**
1. On-site G_ii and DOS: Lehmann vs recursion on bcc Fe, elementwise at
   matched broadening; sum rule ∫DOS dE = electron count.
2. Γ-only Lehmann on an N-site supercell ≡ real-space cluster result.
3. J_ij(R) from the unchanged `exchange` module: Lehmann vs RS route,
   converged in k-mesh — **include the convergence study J_ij vs N_k**;
   intersite G at metallic Fermi surfaces converges slowly in k and this is
   the known disappointment mode of Lehmann-route exchange.
4. Damping ladder: gij_eta from backend E vs RS eta route.
5. Backend D with Σ = 0 ≡ backend E to solver tolerance (the two-backend
   consistency test, permanent in CI).

**Effort/Risk:** M–L / medium (conventions, not math).
**Dependencies:** B1 (conventions), B4 helpful not required.

---

## B3. Bloch spectral functions

**Motivation.** A(k,E) = −(1/π) Im Tr G(k, E+iδ): band-structure-like output
for perfect crystals now, and *the* headline observable once CPA (B8) and
Σ(E) physics (B10/B12) exist — disorder/correlation-broadened bands.

**Design.** Thin consumer of B2: a `post_processing_bsf` registered in
`calculation.f90` (`check_post_processing` + `prepare_post_processing_stack`),
taking the band-path machinery from `post_processing_band_structure` and the
per-(k,z) resolvent from B2 backend E or D. Orbital-, spin-, and
site-projected variants are trace weights, not new machinery.

**Validation.** Sharp-limit test: δ → 0 peaks coincide with eigenvalues from
the band-structure route on the same path (automatable: peak positions vs
`bands` output within δ).

**Effort/Risk:** S / low. **Dependencies:** B2.

---

## B4. GPU port of k-space functionality

**Motivation.** Batched dense linear algebra over k-points is the ideal GPU
workload; accelerates B2/B3/B8/B10/B11 across the board.

**Design.** Extend the proven plugin pattern (C API + `iso_c_binding`
wrapper + CPU fallback, as in `rsrec_cuda_plugin.f90`) with a
`rsksp_gpu` plugin exposing three batched primitives:
1. batched Hermitian eigensolve (`zheevd`-batch / cuSOLVER `syevjBatched`),
2. batched linear solve `getrf/getrs` (for B2 backend D and CPA),
3. batched GEMM for eigenvector-product accumulation (moment/GF assembly).
Fortran side decides batching over k (and energy for backend D); device
memory reuse follows the fingerprint/context lifecycle from Phase 1 (one
shared context, `ensure_context`). HIP variant mirrors CUDA as in `librsrec`.

**Explicitly out of scope:** GPU-porting the tetrahedron bookkeeping or
symmetry reduction — host-side, cheap.

**Validation.** Bitwise-tolerance comparison of eigenvalues/G blocks vs CPU
path for every batched primitive (unit-level, in `tests/cuda`), plus the B2
suite run through the GPU path on hardware (manual matrix per Phase 2 P3).

**Effort/Risk:** M / low-medium (well-trodden).
**Dependencies:** none hard; do alongside or after B2.

**Related decision — GPU strux (from the original wishlist):** gate on a
profile. Structure-constant setup normally amortizes against recursion time;
port only if profiling a large interface/surface workload shows it as a
bottleneck, and then via the same plugin pattern at the
`lattice_strux.f90` ↔ `include_codes/strux_lib` boundary.

---

## B5. Route-agnostic G(E) post-processing (Jij, damping, transport)

**Motivation.** "All G(E) post-processing accessible from Lehmann without
duplication." Mostly falls out of B2 by construction; this blueprint covers
the two items that don't.

**Design.**
1. **Exchange & Gilbert damping:** no code changes. They consume
   `gij/gji/gij_eta` + torque components; B2 fills them. Add a dispatch key
   (`gf_route = recursion | lehmann | dyson`) in the post-processing control
   so the user selects the producer; consumers untouched.
2. **Conductivity (moment-native):** the KPM pipeline consumes Chebyshev
   moments, not gij. Add an **exact moment generator**: in the eigenbasis
   T_n(H̃_k) = U T_n(ε̃) U†, so μ_n are computed exactly from eigenpairs and
   written into the *same* moment arrays consumed by
   `finish_conductivity_moments`. The whole downstream Kubo machinery is
   reused unchanged; the Lehmann route becomes an exact moment producer.
   (Bonus: exact-vs-recursion moments on the same crystal is a direct KPM
   error bound — make it a regression test.)
3. **Gilbert damping unification (Kamberský torque-correlation):** the RS
   implementation exists in the exchange/damping path via the eta ladder;
   after (1) it runs on Lehmann input as-is. Add SOC-derivative torque
   operators only if the current implementation lacks them (audit first —
   report, don't assume).

**Validation.** Same-system triads: J_ij, α (damping), σ (conductivity)
computed via recursion vs Lehmann vs Dyson(Σ=0); agreement within documented
broadening/k-mesh envelopes, pinned as regression cases.

**Effort/Risk:** S–M / low. **Dependencies:** B2 (+B4 for speed).

---

## B6. Dipole charge moments — surface/interface electrostatics
(Skriver–Rosengaard; Vitos treatment)

**Motivation.** Correct work functions, surface dipoles, band alignment at
interfaces. **Prerequisite for B7 being physically trustworthy.**

**Design.** Extend the ASA Madelung problem from monopoles (q_R) to include
l = 1 (and optionally l = 2) multipole moments of the sphere charges:
1. Compute Q_{RL} moments from the density inside spheres (integrals over
   the same radial meshes the charge machinery owns — `charge.f90` is a
   consumer here, but the *new* code lives outside the legacy fence as a
   separate module, e.g. `electrostatics_multipole.f90`).
2. Generalized Madelung matrix M_{RL,R'L'} — structure-constant-like lattice
   sums; slab/2D Ewald handling for surface geometry (reuse/extend the 2D
   conventions from the preprocessing toolkit; L_z remains derived, never
   input).
3. Feed corrected potential shifts v_R (and l=1 components where the
   representation allows) into the SCF potential update.

**Fence note.** The potential-update touchpoint sits close to legacy
`self.f90`; implement as: new module computes shifts, legacy code receives
them through the *existing* variable it already reads (monopole Madelung
shift), extended additively — the minimal-diff strategy inside the fence.

**Validation.**
1. Bulk limit: all dipole moments → 0, total energy and potential shifts
   identical to current code (bit-level regression).
2. fcc(100)/(111) work functions vs published LMTO-ASA-dipole values
   (Skriver–Rosengaard tables) within the accuracy claimed there.
3. Slab thickness convergence of the surface dipole.

**Effort/Risk:** M / medium (Ewald/2D sums are the risk).
**Dependencies:** none. Can run parallel to B2.

---

## B7. Boundary self-energies — interface leads + vacuum GF
(one abstraction, two implementations)

**⚠ Superseded — see `docs/dev/plans/B7_interfaces_and_vacuum_leads.md` for
what shipped.** The design below (López-Sancho/Sancho-Rubio decimation,
principal-layer partition, an energy-dependent boundary self-energy Σ(E))
was the original plan; the implemented B7 took a different, embedding-based
approach — frozen reference regions (bulk/vacuum/lead), a per-site region
registry, a vacuum-lead parameter generator, and a two-sided
deviation-variable electrostatics/alignment solver, with no new
energy-dependent machinery added to the Hamiltonian at all. The text below
is kept as a historical record of the original motivation and is not what a
reader should use to understand the current code.

**Motivation.** Interfacial extension of surface mode (L/R leads bulk or
vacuum) and Skriver–Rosengaard/Turek vacuum GF are the same operation:
terminate a semi-infinite region with an embedding self-energy on the
boundary principal layers.

**Design.** One `boundary_self_energy` type with a common interface
(`build(E) → Σ_L, Σ_R` on boundary-layer orbital blocks), two
implementations:
1. **Bulk lead:** surface GF of the semi-infinite lead by layer decimation /
   López-Sancho iteration on principal-layer H_00, H_01 (built from the
   existing bulk Hamiltonian machinery for the lead material).
2. **Vacuum:** vacuum/free-space boundary GF in the
   Skriver–Rosengaard/Turek construction, replacing today's empty-sphere
   stacks-of-vacuum-layers pragmatics with a true semi-infinite vacuum.

**Where it enters.** Real-space side: as modified terminators/edge
self-energies for recursion clusters whose surface layer touches a lead —
i.e., an energy-dependent Σ on boundary orbitals, which means those
boundary-touched recursions run per energy (accept this cost; interior sites
keep the fast energy-independent machinery). k-space-of-the-2D-cell side:
standard interface GF assembly G = [E − H_C − Σ_L − Σ_R]⁻¹ per (k∥, E) —
another consumer of B2 backend D and B4 primitive 2.

**Validation.**
1. Lead = same material as center ≡ bulk (embedding identity test).
2. Vacuum implementation vs current empty-sphere surface for a thick slab:
   layer-resolved DOS/moments converge to each other away from the surface.
3. Decimation convergence and analyticity checks (Im Σ ≤ 0 on E + iδ).

**Effort/Risk:** L / medium-high (the RS per-energy boundary recursion is the
novel part). **Dependencies:** B6 strongly recommended first; B2/B4 for the
k∥ route.

---

## B8. CPA + DLM in k-space

**Motivation.** Substitutional disorder and finite-temperature spin disorder
(DLM ⇒ paramagnetic state, and with existing J_ij machinery, Tc estimates —
arguably the bigger scientific payoff).

**Design.** Standard TB-LMTO-CPA (Turek et al. formulation):
single-site coherent potential Σ_CPA(E) from the self-consistency
Σ = Σ_c x_c t_c[Σ] on each disordered sublattice; coherent GF via B2 backend
D with Σ_CPA as (the first shipped) Σ provider — CPA is thereby also the
integration test of the B10 provider interface. DLM = CPA with up/down
(or Lebedev-mesh) orientations of the same atom at 50/50 (or Weiss-field
weighted) concentrations.

**Observables staging (honesty about vertex corrections):** ship
CPA-DOS, BSF (B3 handles Σ automatically), total energies and moments first.
Disorder-averaged J_ij (Lichtenstein-CPA) next. Transport with CPA **without
vertex corrections is wrong by construction** — stage vertex-corrected
conductivity (Velický ladder) as its own later work item; until then the
code must refuse or loudly warn on CPA+conductivity.

**Validation.**
1. x → 0/1 limits reproduce pure crystals.
2. Dilute limit vs single-impurity real-space calculation (this code's own
   impurity mode — a rare self-contained cross-check most codes can't do).
3. Published CuNi / FeV CPA-DOS benchmarks.
4. DLM bcc Fe: local moment vs literature; vanishing net magnetization.

**Effort/Risk:** L / medium. **Dependencies:** B2 (backend D), B3 for BSF
payoff, B10 interface co-designed.

---

## B9. CPA / DLM for the real-space part

**Motivation.** Disordered *hosts* for the impurity machinery and
inhomogeneous disorder (concentration profiles at interfaces) — exactly
where the real-space identity shines and k-space CPA can't follow.

**Design.** Effective-medium recursion: sites carry Σ_CPA(E) from a
single-site condition evaluated with the *local* GF from recursion on the
effective medium. Energy-dependent on-site Σ means per-energy recursion for
the medium construction — the same cost pattern as B7's boundary recursion;
share the "per-energy recursion driver" infrastructure between B7 and B9
(design it once, in B7). Impurity mode then embeds a real defect in the CPA
host via the existing Dyson impurity machinery. Inhomogeneous CPA
(layer-resolved concentrations at an interface) composes with B7.

**Validation.** Homogeneous-bulk RS-CPA vs B8 k-space CPA on the same alloy
(the two-route consistency test again); dilute limit vs impurity mode.

**Effort/Risk:** L–XL / high (this is research-adjacent; blueprint promotes
only after B7+B8 exist). **Dependencies:** B7 (shared driver), B8 (reference
results, shared single-site solver).

---

## B10. DMFT via a Σ-provider API (solvers external)

**Motivation.** Colleagues provide the full solver stack (Hubbard-I, FLEX,
CT-QMC); this code provides G_loc, the downfolding/upfolding, and the DFT+DMFT
charge self-consistency. Minimal basis + B4 makes the lattice side cheap.

**Design.**
1. **Σ-provider interface** (co-designed with B8): object providing
   Σ(z) blocks per correlated site/orbital subset on an arbitrary complex
   mesh. Implementations: `sigma_cpa` (B8), `sigma_static` (existing DFT+U,
   as a trivial provider — good plumbing test), `sigma_file/socket`
   (external solver exchange).
2. **DMFT loop:** G_loc from B2 backend D → hybridization Δ(z) extraction on
   the correlated subspace (projector = orbital block in the minimal basis:
   this is where LMTO-ASA is unusually clean) → export (Δ, ε_imp, U-matrix;
   reuse the Slater-integral machinery from the DFT+U code) → external solver
   → import Σ → double-counting (reuse the audited DFT+U DC expressions) →
   iterate. Matsubara vs real-axis meshes both supported by backend D by
   construction (z is just complex).
3. **API alignment with the colleagues' stack is a deliverable in itself:**
   agree on the exchange format (recommend: TRIQS-compatible HDF5 or their
   native format — decide with them, document in the blueprint promotion).
   Hubbard-I remains the in-house fallback solver so CI can run a full loop
   without external dependencies.

**Validation.**
1. Σ = 0 loop ≡ LDA (fixed point identity).
2. Static-Σ provider ≡ existing DFT+U results (reuses the audited
   implementation as reference).
3. Hubbard-I on an isolated-atom limit vs analytic atomic GF.
4. One published DFT+DMFT benchmark with the external solver (e.g. γ-Fe or
   NiO with Hubbard-I) at the semi-quantitative level appropriate to ASA.

**Effort/Risk:** L (lattice side; solvers excluded) / medium.
**Dependencies:** B2 backend D (hard), B4 (strongly), B8 (shared interface).

---

## B11. Linear-response TDDFT — transverse magnons

**Motivation.** True dynamical magnons ω(q) with Landau damping, beyond the
adiabatic J(q) picture; the flagship spectroscopy deliverable.

**Design.** Transverse susceptibility route:
χ₀^{+−}(q, ω) from eigenpair pairs (B2 backend E products across k and
k+q; GPU-batched), then RPA/ALDA Dyson χ = χ₀ [1 − f_xc χ₀]⁻¹ with the
ASA-local xc kernel. Registered as `post_processing_susceptibility`
(insertion points already recorded in `DEVELOPER_MAP.md`). Enforce the
Goldstone condition by the standard kernel-rescaling (sum-rule) correction —
document it, don't hide it.

**Validation (this is where B1 pays off).**
1. Goldstone: ω(q→0) → 0 after sum-rule enforcement; report the bare
   violation as a quality metric.
2. Small-q slope vs spin-wave stiffness from B1 frozen magnons and from
   J_ij sums — three-way consistency.
3. bcc Fe / fcc Ni dispersion vs published LR-TDDFT and experiment
   (semi-quantitative).

**Effort/Risk:** L / medium-high (kernel + Goldstone hygiene).
**Dependencies:** B1 (validation target + conventions), B2, B4 (χ₀ cost).

---

## B12. Couplings roadmap — e–magnon, e–phonon, phonon–magnon

**Motivation.** Long-horizon: interacting quasiparticle physics on top of the
GF frameworks. Recorded here as architecture guidance so B2/B10 decisions
don't foreclose it.

**Design stance (not implementation).**
1. **Electron–magnon:** Σ_em(z) from χ (B11) convolved with the electron GF
   — enters as another Σ provider (B10 interface). First target: spin-polaron
   renormalization of BSF (B3 displays it for free).
2. **Electron–phonon:** phonons/deformation potentials are not in-code;
   design the Σ provider to accept externally computed λ/α²F or g(k,q)
   (from a plane-wave code) rather than computing phonons internally —
   consistent with the FCD decision (don't rebuild what other codes do
   better). Eliashberg-level Σ_ep(iω) on Matsubara via backend D.
3. **Phonon–magnon:** different character (boson–boson); the in-code
   contribution is *adiabatic ingredients*: dJ_ij/du from finite-displacement
   J_ij calculations (the exchange machinery + displaced clusters — a
   real-space specialty). Dynamical coupling itself lives downstream in
   UppASD-style simulations; define the export format (dJ/du tensors) with
   that consumer in mind.

**Effort/Risk:** XL cumulative / research. **Dependencies:** B2, B10, B11.
Promote piecewise; never as one task.

---

## Suggested campaign order

| Wave | Items | Rationale |
|------|-------|-----------|
| 1 | B1, then B2 (+B4 in parallel) | keystone + conventions; B4 is independent mechanical work |
| 2 | B3, B5 | near-free payoffs of wave 1; establish the three-route regression triads |
| 3 | B6, then B7 | electrostatics before boundaries; B7 builds the per-energy RS driver |
| 4 | B8 (+ Σ-provider interface with B10 in mind) | CPA/DLM payoff; provider plumbing proven |
| 5 | B10, B11 | DMFT loop with external solvers; magnons |
| 6 | B9, B12, vertex-corrected CPA transport | research-adjacent tail |

## Promotion protocol
A blueprint is promoted by writing a Phase-1-style task plan (numbered tasks,
per-task validation, ground rules restated) in `docs/dev/`, reviewed before
implementation starts. Every promoted feature must ship its known-answer
tests in the same campaign — an estimator without its tests does not merge
(standing rule).
