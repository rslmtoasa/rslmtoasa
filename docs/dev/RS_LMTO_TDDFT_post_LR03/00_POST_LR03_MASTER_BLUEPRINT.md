# Post-LR-03 master blueprint

## Starting point

The clean-room work has already established the one-electron GF, radial LMTO
augmentation, converged radial magnetic ground state, radial/angular Pauli response
mapping, SR→Pauli approximation boundary, and the canonical transverse-response
conventions.

The remaining problem is now implementation, not discovery of missing basic
objects.

Before every downstream task:

1. confirm branch `fable_v4`;
2. record exact starting HEAD;
3. require a clean worktree unless the task explicitly continues an uncommitted
   blocker audit;
4. read the current evidence documents rather than old blueprints;
5. stop if a prerequisite document itself declares BLOCKED for the capability
   required by the task.

## Architectural spine

The implementation should have four physically distinct layers.

### Layer A — response coordinate space

Owns:

- site / response-harmonic / radial-point / channel indexing;
- radial integration metric;
- operator application rules;
- conversion between raw pointwise kernels and any weighted matrix
  representation.

It owns **no electronic structure**.

### Layer B — one-electron transition vertices

Owns:

\[
T^\mu_{nm;I}(\mathbf k,\mathbf q)
\]

for the certified Pauli/no-SOC baseline.

It consumes electronic structure and LMTO augmentation, but owns no susceptibility
denominators and no XC physics.

### Layer C — bare response

Owns:

\[
\chi^{KS}
\]

from transitions and occupations.

The first backend is the explicit spectral/Lehmann band sum.

A second reciprocal GF backend is independent and must not call the spectral
accumulator internally.

### Layer D — interaction and enhanced response

Separate components:

- direct ALSDA radial kernel;
- Lounis sum-rule interaction;
- optional BES eigenvalue correction;
- Dyson solver;
- loss matrix.

No route may silently rescale another.

## Discrete response-space warning

The direct radial representation is not an orthonormal finite vector space by
default. A physical contraction contains the established radial quadrature.

If raw point values are stored, an operator relation schematically has the form

\[
y_I=\sum_J A_{IJ}W_Jx_J.
\]

Therefore ordinary dense multiplication of raw kernels is generally wrong unless
the weight placement has first been transformed away by a rigorously defined
representation.

LR-04 exists specifically to fix this once.

Every later prompt must use the LR-04 algebra and must not insert or remove
quadrature weights locally.

## Claims gates

### Gate 1 — after LR-05

May claim:

> production transition vectors exist in the certified Pauli response space.

May not claim susceptibility correctness.

### Gate 2 — after LR-06

May claim:

> the bare spectral KS susceptibility is implemented and passes its algebraic and
> numerical oracles.

May not claim magnons or Goldstone correctness.

### Gate 3 — after KXC/GSR/GCR/TDDY

May claim interacting response implementation only.

### Gate 4 — after TDVAL-01

Only then may claim validated Fe/Ni collinear TD-DFT within the stated
Pauli/no-SOC/scalar-relativistic-ground-state approximation.

## Retired old tasks

The old LR-04 radial-interface prompt is retired: LR-01 already established and
persisted the converged radial ground-state contract.

The old LR-05 response-basis prompt is superseded by LR-04 + LR-05 in this pack:
LR-02R already supplied angular, gauge, radial-measure and super-index primitives,
while the production operator algebra and transition evaluator remain to be
implemented.
