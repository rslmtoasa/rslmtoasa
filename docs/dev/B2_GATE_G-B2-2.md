# Gate G-B2-2 — default meshes + documented accuracy for k-space exchange

**Status:** OPEN — awaiting Anders. **Task:** B2.6 (J_ij/damping through the
Lehmann/Dyson-filled arrays + convergence). **Branch:** `fable_v2`.
**Backing evidence:** `docs/dev/reciprocal_green_convergence.md` (the J vs N_k / η
study), commit `c83e1ca`.

## What you are signing (three things)

1. **The normalization question is closed.** Confirm you accept that the k-space
   intersite Green's function is correctly normalized — the factor-≈2 seen when
   comparing a fixed-η Lehmann `J` against the near-sharp recursion `J` is a
   broadening / metallic-Fermi-surface convergence effect, **not** a `1/√2` bug.
   (Evidence: `UnitGammaSupercell` pins the intersite block to the direct resolvent
   `<1e-12`; the ratio is shell- and η-dependent, which no global factor produces;
   the recursion 4-phase algebra reduces to `G_ij` exactly.) → **no code change
   rides on this; it is an acceptance of the physics interpretation.**

2. **The production defaults** for `gf_route='lehmann'|'dyson'` exchange (the two
   numbers below). These are user-set namelist keys today; your sign-off blesses
   the recommended values and, if you want, promotes them to type defaults.

3. **The documented accuracy envelope** — i.e. that the Lehmann-route `J_ij` is a
   *broadening-defined* estimator (report η alongside J), and that matching the
   recursion (sharp-limit) `J` requires the joint `η → 0`, `N_k → ∞` limit. This is
   the honest characterization the gate certifies as "documented accuracy."

## Recommended defaults (for your approval / edit)

| Key | Current | Recommended | Basis |
|---|---|---|---|
| `green_eta` (`&reciprocal`) | 0.01 Ry (type default from G-B2-1); examples use 0.02 | **0.02 Ry** as the documented working point | J k-converges cleanly (~1% at 16³); a defensible physical smearing |
| `nk1=nk2=nk3` (`&reciprocal`) | user-set | **16** (16³ = 4096) for ~1%; **≥20** to quote <1%; raise substantially if η < 0.01 | Table A in the study: 8³ is ~20% high, 16³→32³ flat within ~1% |

## The accuracy statement being certified

- **At fixed η, k-convergence is well-behaved.** 1st-NN J at η=0.02: 8³ = 0.314,
  16³ = 0.259, 32³ = 0.256 — converged to ~1% by 16³. (The metallic "disappointment
  mode" only bites as η → 0, where far finer meshes are needed.)
- **The Lehmann J is a broadening-renormalized quantity.** At the 0.02 Ry working
  point it is *not* equal to the recursion J (1st-NN: 0.256 vs 0.508). Reducing η
  raises it back toward recursion (0.005 → 0.406, and → ~0.51 as η → 0 with k → ∞).
  This is physics, not error: quote `J_ij(η, N_k)` with η stated.
- **Backend D (Σ=0) ≡ E**, so `gf_route='dyson'` gives the same J here; with a real
  Σ (CPA/DMFT, B8/B10) the finite η is physical and this ambiguity disappears.

## Open decisions only you can make

1. **Which operating point is "production"?** Option A (recommended): ship
   η = 0.02 Ry as a *documented broadening-defined* estimator (cheap, ~1% k-converged,
   J is η-renormalized — fine for trends/DMFT/CPA where η is physical). Option B:
   define the deliverable as *recursion-matching* J, which needs the joint
   η → 0 / N_k → ∞ limit (expensive, and currently gated by the backend-E
   full-unreduced-undistributed-BZ requirement — no symmetry reduction / k-parallel
   yet; those attach with B4). Pick A or B (or "A now, B when B4 lands").

2. **Promote the recommended `green_eta`/`nk` to type defaults**, or leave them
   user-set with the study as guidance? (I recommend leaving them user-set and
   pointing the docs at this gate + the study.)

3. **Damping eta-route:** the same `gij_eta` + torque ladder that exchange uses is
   filled identically by both backends, so `calculate_gilbert_damping` should run on
   the k-space arrays with zero changes too. I did **not** run a damping cross-check
   in B2.6 (exchange was the acceptance target). Do you want that as a follow-up
   before the gate signs, or is the exchange acceptance + the shared-fill argument
   sufficient?

## How to sign

Edit this file: set **Status: SIGNED (Anders, <date>)**, record your choices for the
three open decisions and any change to the recommended defaults. I fold the
decisions into code/docs (mirroring how G-B2-1 was closed) and mark B2.6 `[x]` in
the plan checklist.
