# LR-05 — Implement the certified Pauli transition-vector primitive

## Prerequisites

- LR-03 PASS.
- LR-04 PASS.
- LR-02R formal Pauli mapping available.
- The SR→Pauli numerical closure consumed by LR-03 is available.
- If LR-03 did not identify a supported initial Pauli/no-SOC capability, STOP.

## Nature

ONE-ELECTRON RESPONSE-VERTEX IMPLEMENTATION.

No susceptibility accumulator.
No XC kernel.
No Dyson.

## Goal

Implement the production primitive

\[
T^\mu_{nm;I}(\mathbf k,\mathbf q),
\qquad I=(a,L,M,i,\mu),
\]

for the exact baseline certified by LR-02R/LR-03.

This converts two electronic eigenstates into a vector in the production response
space.

## Supported initial capability

Require explicitly:

- `ham_only`;
- orthogonal second-order/HOH;
- collinear;
- no SOC;
- no Hubbard/additive operator;
- `sp` or `spd`;
- Pauli/no-SOC response;
- production direct radial mesh.

Fail closed outside this tuple.

Capability checks must be called by the public evaluator itself. Do not leave them
as an optional caller responsibility.

## Physics

Use the LR-02R definition, not a newly invented one:

\[
\Psi^P_{n\mathbf k,a}(\mathbf r)
=
\frac1r
\sum_{lm\sigma}
c^{n\mathbf k}_{alm\sigma}
U^{n\mathbf k}_{al\sigma}(r)
Y_{lm}(\hat r)|\sigma\rangle,
\]

and

\[
T^\mu_{nm;aLM}(r)
=
\int d\Omega\,
Y^*_{LM}
\Psi_{n\mathbf k}^{P\dagger}
\Gamma^\mu
\Psi^P_{m,\mathbf k+\mathbf q}.
\]

Use the existing:

- LMTO radial augmentation;
- code harmonic convention;
- Gaunt utilities;
- finite-q site gauge;
- LR-03 Pauli/circular operator definitions.

Do not use the scalar-relativistic lower component in this Pauli primitive.

Do not insert `GFAC`.

## API design

The core evaluator should accept an **explicit operator matrix** or a tightly typed
channel whose matrices are owned by the LR-03 convention layer.

Avoid magic integers for `plus_minus`.

A useful separation is:

- electronic endpoint state;
- explicit operator;
- response-space destination vector.

The evaluator should not know about frequency, occupations, or denominators.

## k and k+q endpoints

Use the production arbitrary-k eigensystem.

Do not manually rephase eigenvectors unless required by the LR-02R gauge contract.

The transition vector must be covariant/invariant under reciprocal folding exactly
as derived.

## Degenerate states

Do not attach physical meaning to individual vectors inside an exact degenerate
subspace.

The primitive may evaluate a chosen eigenvector, but tests of physically summed
quantities must be invariant under unitary rotations in degenerate subspaces.

Document this distinction.

## Tests

### 1. Direct angular quadrature oracle

For synthetic coefficient vectors and real production radial arrays, compare the
Gaunt-based production `T` against brute-force independent sphere quadrature.

Cover:

- s-s;
- p-p off diagonal;
- p-d;
- d-d;
- charge;
- sigma-z;
- one transverse Pauli operator.

### 2. Hermitian endpoint relation

For Hermitian `Gamma` and equal endpoints, verify the expected conjugation
relation between `T_nm` and `T_mn`.

### 3. BZ crossing

Use a two-site basis with nontrivial basis positions and a `k+q` point that folds
across the Brillouin-zone boundary.

Verify the complete **production transition vector**, not only the coefficient
gauge.

### 4. Ground-state occupied diagonal closure

Using the real converged fixture consumed by LR-02N/LR-03, sum diagonal occupied
charge and sigma-z transition vectors and reproduce the documented Pauli-projected
ground-state radial densities.

This is a cross-layer numerical oracle.

Do not compare to scalar-relativistic `nSR` and call the difference a failure; use
the already documented `nP`.

### 5. Capability guards

Verify hard rejection of:

- SOC;
- generalized overlap;
- noncollinear;
- `spdf`;
- unsupported additive operators.

## Performance

Correctness first.

It is acceptable initially to evaluate band-pair transition vectors directly.

Avoid allocating the full `Nr × Norb²` product repeatedly if a simple reusable
workspace suffices, but do not introduce caching that changes gauge/provenance
semantics.

## Deliverable

Create:

`docs/LR_PAULI_TRANSITION_VERTEX.md`

Include:

- public API;
- supported capability tuple;
- exact equation;
- endpoint/gauge convention;
- radial and angular normalization;
- tests and claim levels;
- known SR→Pauli approximation inherited from LR-02N/LR-03.

## Forbidden

- no occupation difference;
- no frequency denominator;
- no `chiKS`;
- no Kxc;
- no Ward enforcement;
- no site-only vertex;
- no empirical normalization factor.

## Checklist

- [x] public evaluator fail-closed;
- [x] explicit LR-03 operator convention used;
- [x] production LMTO augmentation reused;
- [x] production Gaunt convention reused;
- [x] full response vector emitted;
- [x] independent angular oracle passes;
- [x] folded k+q end-to-end test passes;
- [x] occupied Pauli ground-state closure passes;
- [x] unsupported modes rejected;
- [x] no susceptibility code added.

## Implementation outcome

`source/lr_pauli_transition_vertex.f90` provides the public evaluator and
typed capability/endpoint contracts. `UnitLrPauliTransitionVertex` and its
five expected-failure capability variants pass. The implementation claim is
the LR03/LR04 algebraic one: production Pauli/no-SOC transition vectors exist
in the certified direct response space. It does not claim exact
scalar-relativistic response, susceptibility correctness, XC correctness, or
material/literature validation.

## Commit

`linear-response: add Pauli radial transition vertices`
