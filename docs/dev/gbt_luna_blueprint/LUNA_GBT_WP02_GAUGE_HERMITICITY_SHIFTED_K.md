# Luna Prompt — GBT-WP02: Gauge Algebra, Hermiticity, and Exact k±q/2 Identity

## Mission

Prove that the current single-q real-space GBT link implements the generalized Bloch theorem itself, independent of SCF and total-energy integration.

This is the strongest low-level certification stage.

## Dependency

WP01 must pass.

## Scientific oracle

For a one-sublattice collinear scalar-relativistic reference, the generalized Bloch construction must be equivalent to ordinary spin sectors evaluated at shifted momenta `k +/- q/2`, modulo the repository's Fourier/sign convention.

Do not assume the sign. Derive it from:

- the definition of `gbt_bond_phase`;
- the Fourier exponent used by the reciprocal Hamiltonian builder;
- the spinor ordering.

Then verify the derived identity numerically.

## Required algebraic tests

### A. Endpoint-link unitarity

For generic q, displacement, and endpoint reference frames:

\[
G^\dagger G=I.
\]

### B. Bond reversal

For each tested directed bond:

\[
G_{ba}(-\mathbf d)=G_{ab}^{\dagger}(\mathbf d).
\]

Test both one- and multi-sublattice reference frames.

### C. Closed-loop gauge identity

For deterministic closed loops in the lattice graph:

\[
G_{12}G_{23}\cdots G_{n1}=I.
\]

Use loops that include nontrivial lattice translations and, if available, multi-sublattice endpoints.

### D. Real-space Hamiltonian bond reversal

After LMTO contraction in the minimal feature set:

\[
H_{ji}(-R)=H_{ij}^{\dagger}(R).
\]

### E. Reciprocal Hermiticity

At generic k and q:

\[
H(k)=H^{\dagger}(k).
\]

Use strict roundoff-scale thresholds.

## Exact shifted-k test

Build a one-sublattice collinear fixed-potential fixture.

For several generic deterministic `(k,q)` pairs:

1. build `H_GBT(k,q,theta=0)` through the production GBT path;
2. independently build ordinary spin-resolved Hamiltonians at `k+q/2` and `k-q/2` through the ordinary path;
3. determine the spin-block ordering and q signs from the code convention;
4. compare either matrices in the appropriate block basis or complete sorted eigenvalue multisets.

Prefer direct matrices where feasible, followed by eigenvalues as an integration check.

Acceptance: roundoff-level agreement.

## q ↔ -q symmetry

For SOC-free centrosymmetric fixtures, confirm the appropriate spectral symmetry under q reversal. Keep this distinct from the stronger shifted-k identity.

## Negative controls

At least one test-only negative control should prove sensitivity to one of:

- missing factor `1/2` in spinor phase;
- q sign flip;
- reversed endpoint order;
- wrong full displacement.

Do not commit intentionally broken production code; implement the negative control in isolated test logic if possible.

## External cross-check

Document the correspondence to:

- ELK's explicit `k +/- q/2` spin-spiral construction;
- an LMTO spin-spiral structure-matrix formulation.

Do not copy external code. Use it to verify mathematical convention.

## Deliverables

- automated algebraic unit tests;
- an exact shifted-k functional test;
- `docs/dev/GBT_REAUDIT_WP02_GAUGE_SHIFTEDK.md` with the sign derivation and numerical residuals.

## Completion checklist

- [ ] Link unitarity tested.
- [ ] Bond reversal tested.
- [ ] Closed-loop pure-gauge identity tested.
- [ ] Contracted real-space bond Hermiticity tested.
- [ ] Reciprocal Hamiltonian Hermiticity tested.
- [ ] Fourier/q sign convention derived explicitly.
- [ ] `k+q/2` / `k-q/2` identity passes at generic points.
- [ ] q↔-q spectral symmetry checked.
- [ ] At least one meaningful negative control included.
- [ ] ELK and LMTO correspondence documented.
- [ ] No SCF or energy-fitting changes made.

## Commit

`test: prove GBT gauge and shifted-k identities`
