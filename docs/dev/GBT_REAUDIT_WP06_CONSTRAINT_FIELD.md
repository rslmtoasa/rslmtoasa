# GBT re-audit WP06: constraint-field covariance and nonzero-field fixture

## Status and scope

This report records WP06 on source baseline commit
`46d5a726c1fe57e89c85267e9fa098b1a824ee77` with the WP06 changes present in
the working tree.  WP05's explicit frame contract is a prerequisite and is
consumed here through `gbt_frame_for_cell` and
`gbt_rotating_to_lab_vector`.

The verification build is CMake Debug mode with GNU Fortran 13.3.0, OpenMP
enabled, MPI/libXC/CUDA disabled, `RUN_UNIT_TESTS=ON`.  The focused test set
was run through CTest after rebuilding the affected targets.

## Audited contract

The three concepts remain separate:

1. `constraint_reference`/`constraint_target` is the fixed target direction;
2. `mag_cfield` is the onsite Pauli coefficient consumed by the Hamiltonian;
3. `constraint_metric` is the numerical SCF convergence merit.

The returned penalty `constraint_energy` is logged/reported as a separate
controller quantity; it is not added to the ordinary DFT total-energy signal
used by the SCF energy convergence or benchmark profile.

For ordinary real-space/reciprocal calculations, targets and fields are in the
active Cartesian spin frame.  For `gbt_single_q`, the primitive target and
`mag_cfield` are in the cell-independent rotating frame.  A GBT explicit-site
comparison reconstructs only the lab field:

```text
B_con(n,a,lab) = R[U(n,a)] B_con(a,rotating)
```

No cell phase is applied to the onsite field during the GBT Hamiltonian build.
The GBT onsite path adds the same Pauli block as the ordinary path, before the
Cartesian-to-spherical conversion, at bond slot `m=1`.  The target and field
therefore do not acquire an endpoint/bond phase.

`mag_cfield` is stored and updated in Rydberg energy units.  The Hamiltonian
uses

```text
H_spin = H_0 sigma_0 + B_x sigma_x + B_y sigma_y + B_z sigma_z
```

so a positive `Bz` raises the up block and favours negative electronic spin.
For method 3, a moment canted in `+x` therefore produces a coefficient in
`-x`; the physics test checks this sign and the finite-difference penalty
derivative.

## Constraint state and GBT reference boundary

The initial `constraints_mom_ref` and `constraints_bfield` arrays are now
checked for exact shape `(3,nrec)` rather than silently discarded.  A zero
target is rejected.  Each active site stores the fixed target in
`symbolic_atom%constraint_target`.

The measured integrated moment is retained in
`symbolic_atom%constraint_actual`.  This is needed because the GBT operator
boundary deliberately exposes `potential%mom = (0,0,+/-1)` as a collinear
rotating-frame reference marker: its z sign selects the radial channel, while
the physical cone remains in the explicit GBT angle controls.  The marker
prevents a transverse lab/rotating moment from being silently accepted as a
GBT reference.  The controller and magnetic mixer use the retained measured
vector, not the marker.

The radial density axes in both the real-space and reciprocal producers use
`constraint_target` when present.  Thus a fixed target cannot be replaced by
the temporary GBT marker during density projection.

## Controller and transversality

`constraint_diagnostics` is returned by `constrain()` for every update.  It
contains, per site, target and actual unit vectors, moment magnitude, angular
error, transverse residual, field components and magnitude, field dot target,
and the target/actual cross-product torque.  The aggregate metric is RMS
angular error in radians.  The optional detailed file is
`constraint_diagnostics.dat`; it is enabled by `constraints_diagnostics` and
has a declared active-frame/Ry/radian header.  Legacy verbose controller
stdout is behind the same switch.

The constraint formalism distinguishes the controller variants:

- method 2 retains the historical non-orthogonal Lagrange update and may have
  a longitudinal component;
- methods 3, 4, and 5 are transverse controllers.  Their output field is
  projected so `B_con dot e_target` is zero, including a longitudinal restart
  seed and the low-moment PID path.

When constraints are enabled, SCF convergence requires both the ordinary energy
criterion and `constraint_metric <= constraints_tolerance`.  The default
tolerance is `1e-6` rad and detailed diagnostics do not pollute normal output.

## Required nonzero-field fixture

`tests/gbt_wp6_constraint/README.md` describes the fixture and
`tests/unit/test_gbt_constraint_covariance.f90` implements it.  It is a
deterministic frozen-response oracle: a zero-field snapshot supplies a canted
unit moment with a measurable orientation error, then eleven successive
prescribed snapshots reduce the angle by a factor of four below `1e-6` rad.
Every snapshot recomputes a nonzero method-3 field and checks the reported
magnitude, angle, transverse residual, and torque.  This isolates the constraint controller from a
torque-free material fixture, which would not prove that a field is actually
needed.

The same executable also checks that the PID and Ma--Dudarev variants remove a
longitudinal restart component.  Method 2 remains intentionally non-
orthogonal, as documented above.

## Explicit period-four covariance

The fixture uses `q=(1/4,0,0)` in the repository's `2*pi/alat` convention and
translations `n=(0,1,2,3)`.  For each site it rotates the local target and
actual moment through the authoritative WP05 frame, runs the constraint update
in that lab-frame site basis, and compares the result with the rotated
primitive field.  It requires:

- field-vector equality cell by cell;
- field-magnitude equality;
- equal target/actual angle;
- zero field component along the site target;
- phases `0, pi/2, pi, 3*pi/2`.

This is the explicit-supercell covariance check requested by WP06, rather than
an energy-only comparison.

## Verification

The following tests pass in the stated build:

```text
UnitConstrainingField
UnitGbtConstraintCovariance
UnitConstrainingFieldSource
UnitGbtWp06Source
```

The existing physics test covers the Hamiltonian sign, penalty derivative,
global spin rotation covariance, and transverse restart-field rejection.  The
new executable covers the nonzero-field residual/field fixture and the
period-four GBT/supercell covariance.  The source-contract test covers the
single insertion path, fixed-target storage, GBT marker/actual split, and
diagnostic/convergence plumbing.

## Evidence boundary

Established by WP06:

- target and field storage frames are explicit;
- GBT onsite fields are unphased and inserted once;
- the Hamiltonian sign and Ry units are documented and tested;
- a nonzero field is required and remains nonzero in the converged fixture;
- angular residual, torque, transversality, and controller diagnostics are
  machine-readable;
- primitive GBT and period-four explicit-site fields are covariant.

Not established by this focused WP06 gate:

- convergence of a full material SCF run with constraints enabled;
- constrained frozen-magnon or LKAG closure;
- GBT support for SOC, local-axis, or relaxed-reference combinations rejected
  by the existing operator boundary.
