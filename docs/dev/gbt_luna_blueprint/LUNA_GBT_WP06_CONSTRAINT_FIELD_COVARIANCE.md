# Luna Prompt — GBT-WP06: Constraining-Field Covariance and Real Nonzero-Field Fixture

## Mission

Validate the existing constraining-field machinery as actual constrained magnetism under GBT, separately from any density-projection/mixing policy.

A successful test must require and converge a **nonzero** constraint field.

## Dependency

WP05 frame contract must be explicit and tested.

## Concepts that must remain separate

1. **Reference/density policy:** whether the target magnetic direction is held fixed during SCF/mixing.
2. **Constraining field:** onsite Pauli field adjusted to cancel transverse torque and enforce the target direction.
3. **Controller penalty/merit function:** any numerical object used to converge the field.

Do not merge these concepts in names, logging, or energy semantics.

## Audit current implementation

Trace:

- source and storage of target directions;
- call into `constrain.f90`;
- update equation/controller;
- units of `mag_cfield`;
- frame of `mag_cfield`;
- sign with which it enters the Hamiltonian;
- insertion into ordinary and GBT onsite magnetic blocks;
- convergence criterion.

Verify that GBT does not apply any translational phase to the onsite field.

## Define the GBT constraint-frame convention

Recommended:

- store the primitive-cell constraint target in the GBT rotating frame;
- store/update `mag_cfield` in the same rotating frame;
- reconstruct lab-frame site fields only for output and explicit-supercell comparison.

For cell n and sublattice a:

\[
\mathbf B^C_{na,lab}=\mathcal R[U_{na}]\,\bar{\mathbf B}^C_a.
\]

## Required nonzero-field fixture

Do **not** validate only a torque-free Fe spiral.

Choose at least one state known/constructed to be nonstationary under the frozen ordinary potential, with a target direction displaced enough that the unconstrained result develops a measurable transverse residual.

Recommended hierarchy:

1. simple deterministic toy/Stoner-like fixture if existing infrastructure supports it;
2. fcc Ni single-sublattice finite-q finite-cone state;
3. B2 FeCo multi-sublattice state.

The unconstrained control must demonstrate a nonzero orientation error/torque. The constrained run must then reduce it while producing nonzero `B^C`.

## Diagnostics to implement

Per sublattice and iteration, expose machine-readable values for:

- target unit vector;
- actual integrated moment unit vector;
- moment magnitude;
- angular error;
- transverse moment residual;
- `B^C_x`, `B^C_y`, `B^C_z` in the declared frame;
- `|B^C|`;
- `B^C dot e_target`;
- torque proxy / cross product;
- controller convergence metric.

Avoid excessive default logging; put detailed iteration traces behind a diagnostics switch or dedicated evidence fixture.

## Explicit-supercell covariance check

For a commensurate spiral, compare the converged primitive GBT constraint field with the explicit-supercell site fields after frame rotation.

Require:

- matching magnitudes;
- correct lab-frame orientation cell by cell;
- matching target/actual moment angles.

This is a stronger check than energy alone.

## Transversality

Determine from the actual constraint formalism whether the physical control field should be purely transverse to the target direction. If yes, test

\[
\mathbf B^C\cdot \mathbf e^{target}\approx0.
\]

If the controller intentionally contains a longitudinal component, document why and what part is physical versus numerical.

## Deliverables

- nonzero constraint regression fixture;
- GBT/supercell constraint covariance test;
- clear diagnostics;
- `docs/dev/GBT_REAUDIT_WP06_CONSTRAINT_FIELD.md`.

## Completion checklist

- [ ] Constraint target frame identified.
- [ ] `mag_cfield` frame identified and made explicit.
- [ ] Hamiltonian coupling sign/units verified.
- [ ] Onsite/no-phase property verified.
- [ ] Unconstrained control shows measurable transverse error.
- [ ] Constrained run converges a nonzero field.
- [ ] Target-angle residual converges to declared tolerance.
- [ ] GBT/supercell field covariance passes.
- [ ] Transverse/longitudinal field content is understood.
- [ ] Detailed diagnostics are available without polluting normal output.

## Commit

`test: validate GBT constraining-field covariance`
