# B2.6 — Lehmann/Dyson exchange: route agreement + k-mesh/broadening convergence

**Milestone:** B2 (`reciprocal_green`), task B2.6. **Branch:** `fable_v2`.
**Status:** study complete; **gate G-B2-2** (default meshes + documented accuracy)
awaits Anders' signature. This document is the accuracy record the gate signs off.

## 1. What B2.6 had to establish

Two things, per the plan (`docs/dev/plans/B2_reciprocal_green.md` §2 B2.6) and the
blueprint (`docs/dev/plans/BLUEPRINTS.md` §B2, validation item 3):

1. **Zero-consumer-change exchange.** `post_processing='exchange'` with
   `gf_route='lehmann'` (or `'dyson'`) must run `exchange.f90` **unchanged** on the
   k-space-filled `green` arrays and produce a physical `J_ij`.
2. **The open `1/√2` intersite normalization question** (carried since B2.2) must
   be closed: is the k-space intersite block correctly normalized against the
   recursion route the exchange consumer was written for?

Both are now settled. The exchange consumer runs unchanged (see §2), and the
intersite normalization is **correct** — the factor-≈2 one sees in a naive
fixed-broadening `J` comparison is a **broadening / metallic-Fermi-surface
k-convergence artifact, not a normalization error** (see §3–§4).

## 2. Zero-consumer-change acceptance

`post_processing='exchange'` + `gf_route='lehmann'|'dyson'` fills the SAME
canonical `green` arrays (`gij/gji`, the `gij_eta` ladder, and the torque families
`Ginmag`/`Gi{x,y,z}` + j-partners) that the recursion route fills, then calls
`calculate_exchange` / `calculate_exchange_twoindex` **with no consumer changes**.
On collinear bcc Fe the k-space route returns the correct tensor structure:
isotropic `J` on the diagonal, `D_ij ≈ 1e-16` (zero DMI, as required by the
collinear symmetry), `A_ij` diagonal. The dispatch is the B2.5 key; B2.6 is its
acceptance.

Reproduction (reuses the tracked `example/exchange/bccFe_kspace_green` potential;
change only the `&calculation` and pair inputs):

```
&calculation
  post_processing = 'exchange'
  gf_route        = 'lehmann'     ! or 'dyson' (Σ=0 ≡ lehmann) or 'recursion'
/
&lattice
  ...pbc bcc cell, n1=n2=n3=10...
  njij = 1
  ijpair(1, :) = 1, 335           ! 1st-NN pair: atom 1 (0,0,0), atom 335 (-½,-½,-½)
/
&reciprocal
  nk1 = 20  nk2 = 20  nk3 = 20
  use_symmetry_reduction = .false.  ! backend E needs the full unreduced BZ
  green_backend = 'lehmann'
  green_eta = 0.02
/
```

## 3. Why a fixed-broadening `J` comparison is misleading

The two routes broaden the Green's function **differently**:

- **Recursion** fills the energy-grid `gij` with the analytic continued-fraction
  terminator (`block_green_ij`/`chebyshev_green_ij` pass `eta = (0,0)` to `bgreen`).
  Effective broadening is very small — the on-site DOS is essentially a sharp
  spectrum (the `kspace_green` driver shows recursion `dos_rs ≈ 1e-5` at the band
  edge where the η=0.02 Lehmann DOS is ≈ 0.086).
- **Lehmann/Dyson** broaden by the explicit retarded `green_eta` (`z = E + iη`).

`J_ij` is a Fermi-surface-weighted **integral** `∝ ∫^{E_F} Im Tr[d G_ij d G_ji] dE`.
Unlike the DOS **sum rule** (broadening-independent — both routes integrate to the
electron count, weight → nb = 18), this integral is strongly broadening-dependent
for a metal. So comparing `J` at a single finite η against the near-sharp recursion
value does **not** test normalization; it tests broadening.

## 4. The normalization is correct — three independent lines of evidence

**(a) Machine-precision kernel identity (the decisive one).** The Γ-only supercell
unit test `tests/unit/test_gamma_supercell.f90` (ctest `UnitGammaSupercell`) pins
the strict-Lehmann **intersite** block `lehmann_pair_block(...,0,m,...)` against the
direct cluster resolvent element `[zI − H_sc]^{-1}_{0,m}` for **every** `m`, to
`< 1e-12`. `UnitLehmannChain` independently pins the intersite `e^{ik·ΔR}` phase
against a residue closed form (`4.2e-16`). The k-space intersite kernel therefore
carries **no** `1/√2` (or any) normalization error.

**(b) The recursion 4-phase algebra reduces to `G_ij` exactly.** `recur_b_ij`
builds four unit-norm start vectors `(|i>±|j>)/√2`, `(|i>±i|j>)/√2` (norm² = 1 for
`i≠j`), so `g0_n = <ψ_n|G|ψ_n> = ½(G_ii+G_jj) ± …`. The combiner in
`calculate_intersite_gf_core`, `gij = ½[(g0₁−g0₂) + (1/i)(g0₃−g0₄)]`, gives
`½[(G_ij+G_ji) + (G_ij−G_ji)] = G_ij` exactly. Both routes therefore target the
identical object; and `pauli_decompose_block` reproduces the recursion route's
(charge, x, y, z) decomposition (green.f90:593–596) line for line. The on-site
`i==j` special case (unit-norm start) was the B2.2 `×0.5` fix; the `i≠j` case never
needed it.

**(c) Shell- and η-dependence rule out a global factor.** A normalization factor
would be a single constant on every bond at every broadening. The measured ratio
`J_recursion / J_lehmann(η=0.02)` is **1.98 on the 1st-NN bond but 1.31 on the
2nd-NN bond**, and the converged `J_lehmann` **swings with η from +0.41 to −0.17**
(Table B). No constant does that. As `η → 0` the Lehmann `J` climbs back toward the
recursion (near-sharp) value.

### Table A — `J_ij` vs N_k, 1st-NN bond, η = 0.02 Ry (mRy, ×10³/4π)

| N_k (nk³) | J_ij (Lehmann) |
|-----------|----------------|
| 8³  = 512    | 0.3142 |
| 16³ = 4096   | 0.2593 |
| 20³ = 8000   | 0.2557 |
| 24³ = 13824  | 0.2564 |
| 28³ = 21952  | 0.2555 |
| 32³ = 32768  | 0.2557 |

k-converged to ~1% by **16³**; 8³ is ~20% high. This is the well-behaved,
broadening-regularized convergence (finite η tames the metallic Fermi-surface
sampling — the "disappointment mode" only bites as η → 0, where far finer meshes
are needed).

### Table B — `J_ij` vs broadening η, 1st-NN bond, N_k = 20³ (mRy)

| η (Ry) | J_ij (Lehmann) |
|--------|----------------|
| 0.005  |  0.406 |
| 0.01   |  0.400 |
| 0.02   |  0.256 |
| 0.04   |  0.023 |
| 0.08   | −0.167 |
| recursion (η ≈ 0) | 0.508 |

### Table C — `J_ij(R)` two shells, N_k = 20³, η = 0.02 Ry (mRy)

| shell | pair | recursion (η≈0) | Lehmann (η=0.02) | ratio |
|-------|------|-----------------|-------------------|-------|
| 1st-NN | (1, 335), `a/2·(1,1,1)` | 0.508 | 0.256 | 1.98 |
| 2nd-NN | (1, 336), `a·(0,0,1)`   | 0.386 | 0.296 | 1.31 |

## 5. Recommended defaults (for gate G-B2-2)

- **k-mesh:** `nk1=nk2=nk3 = 16` (16³ = 4096) is the recommended production default
  for bcc-Fe-like metals at `green_eta = 0.02 Ry` — J converged to ~1%. Use `≥ 20³`
  when quoting to better than 1%, and expect to raise the mesh substantially if η is
  pushed below ~0.01 Ry.
- **Broadening:** `green_eta = 0.02 Ry` is a reasonable default. It is a genuine
  physical smearing of `J`, **not** a numerical nuisance to be minimized: quoting a
  Lehmann-route `J_ij` requires quoting its η (and, for η ≲ 0.01, demonstrating
  k-convergence). To reproduce the near-sharp recursion `J`, take `η → 0` **with**
  `N_k → ∞` jointly.
- **Backend D ≡ E:** with Σ = 0 the Dyson route reproduces Lehmann to solver
  tolerance (permanent invariant `UnitDysonEquivalence`; end-to-end 4.7e-11 on the
  full `gij`), so `gf_route='dyson'` gives the same `J` as `'lehmann'` here.

## 6. Known-answer tests backing this milestone

Normalization/correctness is pinned at machine precision by the existing unit
suite — no new machine-precision test is possible for the real-material `J`
(convergence/broadening-limited by physics):

- `UnitGammaSupercell` — Γ-only supercell Lehmann ≡ direct resolvent, **intersite**
  blocks, `< 1e-12`; two-route (Lehmann vs Dyson) DOS `6.2e-15`.
- `UnitLehmannChain` — on-site `8.7e-16`, intersite phase `4.2e-16`, 1/N_k `2e-12`,
  Pauli-decomposition transcription guard `1.1e-16`.
- `UnitDysonEquivalence` — D(Σ=0) ≡ E, chain + multiband bond-phase + Σ sign pin,
  all `~1e-16`.

The real-material artifact is this convergence study on the tracked
`example/exchange/bccFe_kspace_green` potential (report-only; the acceptance is the
maintainer gate G-B2-2).
