# B10 — DMFT via a Σ-provider API (solvers external, Hubbard-I in-house)

**Effort:** L for the lattice side (solvers excluded) / medium risk.
**Lead:** OPUS (API + downfolding + double counting), SONNET (IO/formats).
**Depends on:** B2 backend D (hard), B4 (strongly), B8 (shared provider
interface — co-design, B8 is the first shipped provider).
**Division of labor:** colleagues provide the solver stack (Hubbard-I,
FLEX, CT-QMC); this code provides G_loc, downfolding/upfolding, and the
DFT+DMFT charge self-consistency. Minimal basis + B4 makes the lattice
side cheap.

## 1. Design spec

### 1.1 The `sigma_provider` interface (co-owned with B2/B8)
Object providing Σ(z) blocks per correlated site/orbital subset on an
arbitrary complex mesh (Matsubara and real-axis both — backend D takes
any complex z by construction). Implementations, in build order:
1. `sigma_zero` (B2.1 — exists first, keeps everything testable),
2. `sigma_static` — the existing audited DFT+U potential recast as a
   trivial (z-independent) provider. **This is the cheapest end-to-end
   plumbing test of the whole architecture:** loop with `sigma_static`
   must reproduce existing DFT+U results identically.
3. `sigma_cpa` (B8),
4. `sigma_file` / `sigma_socket` — external-solver exchange.

### 1.2 The DMFT loop
```
G_loc(z) = projector · [B2 backend D fill] · projector†     (correlated block)
Δ(z)     = z − ε_imp − G_loc⁻¹(z) − Σ(z)                     (hybridization)
export (Δ, ε_imp, U-matrix)  →  external solver  →  import Σ(z)
double counting (reuse audited DFT+U DC expressions)  →  mix  →  iterate
```
- **Projector:** the correlated subspace is an orbital block in the
  minimal LMTO basis — LMTO-ASA is unusually clean here; no Wannier
  construction. The projector convention (screened vs orthogonal
  representation for the correlated block) is a pin with a unit test:
  static-limit consistency with how DFT+U applies its potential.
- **U-matrix:** reuse the Slater-integral machinery from the DFT+U code
  (audited in the DFT+U/+J/+V campaign — 5 bugs fixed; treat that code as
  the reference, do not duplicate).
- **Double counting:** FLL/AMF expressions from the audited DFT+U module,
  exposed per correlated site; scheme name recorded in output headers.
- **Charge self-consistency:** stage 2 — v1 is one-shot (fixed DFT
  potential); the charge-SC loop (Σ-dependent density → new potential) is
  its own task after the one-shot loop validates.

### 1.3 Exchange format (Gate G-B10-1)
Agree the Σ/Δ exchange format **with the collaborators** before coding
`sigma_file`: recommend TRIQS-compatible HDF5, accept their native format
if that is faster politically. Anders owns this gate; the decision is
recorded in `docs/dev/dmft_exchange_format.md` with a versioned schema.
Hubbard-I remains the in-house fallback solver so CI runs a full loop with
no external dependencies.

## 2. Validation (known-answer)
1. **Σ = 0 fixed-point identity:** the full loop with `sigma_zero`
   reproduces the LDA charge/DOS exactly (no drift after N iterations).
2. **`sigma_static` ≡ DFT+U:** existing audited results as reference,
   elementwise on the occupation matrices and total energy.
3. **Hubbard-I atomic limit:** isolated-atom (hopping → 0) GF vs the
   analytic atomic Green's function (closed form for given N, U, J).
4. One published DFT+DMFT benchmark at ASA-appropriate semi-quantitative
   level (γ-Fe or NiO with Hubbard-I; cite the comparison figure in the
   test doc).

## 3. Tasks
- **B10.1 [OPUS]** Provider interface finalization (with B2.1/B8.1 already
  consuming it) + `sigma_static` + validation 2.2. *Kit:* this file; B2.1
  provider stub; the DFT+U potential-application site in
  `hamiltonian_hubbard.f90`.
- **B10.2 [OPUS]** Projector + G_loc + Δ extraction + static-limit pin
  test. *Kit:* B10.1; B2 module interface; this file §1.2.
- **B10.3 [OPUS]** In-house Hubbard-I solver + atomic-limit test 2.3.
  Self-contained module; the Slater/U machinery pasted from
  `hamiltonian_hubbard.f90` accessors. *Kit:* this file; the
  Slater-integral routines (±60 lines).
- **B10.4 [SONNET]** Loop driver + mixing + DC wiring + namelist + Σ=0
  fixed-point test 2.1. *Kit:* B10.1–B10.3 interfaces; B8.3 loop driver as
  the structural template.
- **B10.5 [SONNET]** `sigma_file` per the G-B10-1 schema + round-trip
  tests; CI full-loop case with Hubbard-I. *Kit:* B10.4; the signed schema
  doc.
- **B10.6 [OPUS]** Published benchmark case (2.4) + charge-SC scoping note
  (stage-2 decision memo for Anders).

## 4. Checklist
- [ ] B10.1 providers + DFT+U identity
- [ ] B10.2 projector + Δ
- [ ] B10.3 Hubbard-I + atomic limit
- [ ] B10.4 loop + Σ=0 fixed point
- [ ] B10.5 file exchange (G-B10-1 signed)
- [ ] B10.6 benchmark + stage-2 memo
