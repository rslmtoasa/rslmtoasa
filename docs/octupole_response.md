# Orbital octupole response

These are local orbital-angular-momentum multipoles, not collective magnetic
multipoles of a multi-site spin configuration. L matrices use hbar=1.

## Components and normalization

For a != b, the new pure rank-three component is

    Oabb = {La,Lb^2}/2 - {La,L^2}/10 - La/10
         = {La,Lb^2}/2 - (L^2 + I/2) La/5.

Oxyz is the six-permutation symmetrization of Lx Ly Lz divided by six.
The physical PDF convention (cubic angular momentum divided by hbar^2)
is obtained by multiplying these dimensionless matrices by hbar.

| Preferred pol | Harmonic-label alias | Old input (now rejected) | Explicit mixed observable |
|---|---|---|---|
| Oxyy | Ox5y2r2 | Ox3y2r2 | Mx3y2r2 |
| Oyzz | Oy5z2r2 | Oy3z2r2 | My3z2r2 |
| Ozxx | Oz5x2r2 | Oz3x2r2 | Mz3x2r2 |
| Oxzz | Ox5z2r2 | Ox3z2r2 | Mx3z2r2 |
| Oyxx | Oy5x2r2 | Oy3x2r2 | My3x2r2 |
| Ozyy | Oz5y2r2 | Oz3y2r2 | Mz3y2r2 |
| Oxyz | -- | Oxyz remains valid | -- |

Lowercase aliases are accepted. Bare xyy/yzz/zxx/xzz/yxx/zyy/xyz also work.
Harmonic-label aliases have exactly the same normalization as Oabb: they do
NOT multiply the operator by five. These seven components are independent
but need not be mutually orthogonal or equally normalized.

The explicit M labels reproduce the old formula {La,Lb^2}/2 - {La,L^2}/6.
They mix dipole and octupole and must not be reported as pure octupoles.
Old Oa3b2r2 inputs fail with a migration message rather than silently changing
the physical observable or silently retaining a mixed moment.

Pure components vanish in s and p, have zero dipole projection in every shell,
and obey sum_a [La,[La,O]]=12 O. The distinction applies to both accumulation
and conventional current. Selecting an M label deliberately requests a mixed
observable, including its dipole-current contribution.

## Input slots

| Setting | Role |
|---|---|
| linear_out | Measured response operator, slot a |
| pol_alpha | Multipole component of the output operator |
| v_alpha | Output current direction, ignored by a local accumulation |
| linear_in | Driving operator, slot b |
| pol_beta | Multipole component of the input operator |
| v_beta | Input current direction; electric-field direction for charge input |

Both slots accept `octupole` and its synonym `octupole_current` for J={O,v}/2,
and `octupole_accumulation` for the local O. Operator keywords are lowercase.
Use `charge` input for the electric-field-driven response. Other input
operators require defining the conjugate external perturbation and units.
Directions are normalized; their magnitudes are not field strengths. Both
vectors must be finite and nonzero because setup builds both bare velocities.
pol_alpha/pol_beta have length 10; every documented component fits. Defaults
are 'z', which is not an octupole component: set the relevant pol explicitly.

## Hall-current example

Merge `example/transport/octupole_hall.nml` into the matching sections of a
complete converged material input. It is a response-settings fragment, not
a standalone material fixture. It requests J_y^(Oxyy) driven by E_x,
i.e. the transverse coefficient sigma_yx^(Oxyy). A nonzero value depends on
the material symmetry; this input does not guarantee a nonzero Hall effect.
Do not create duplicate &control or &hamiltonian sections.

For accumulation, replace linear_out='octupole' with
linear_out='octupole_accumulation'. The accumulation example provides this.
For an input octupole current instead, use linear_in='octupole', pol_beta,
and v_beta; set linear_out independently with pol_alpha and v_alpha as needed.

## Current scope and limitations

Both octupole accumulation and current responses require hoh=.false. and
nmax=0. Setup explicitly rejects HOH and local impurity Hamiltonians because
the existing velocity path is incomplete. The new current builder also rejects
HOH directly. A broader HOH velocity repair is deliberately separate.

This change adds and tests the operators, not a validated material conductivity.
Previously identified energy-integration bounds, normalization/unit conversion,
and local-projection issues in the common conductivity code remain outstanding.
Quantitative Hall results must wait for those repairs and end-to-end validation.
The current is the conventional anticommutator; no conserved-current or
torque-dipole correction is implied.

Malformed &control and &hamiltonian namelists now fail immediately rather than
continuing with partially read settings. Obsolete cond_type/js_alpha inputs
must be migrated to linear_in/linear_out and pol_alpha/pol_beta.

## Validation

Run from the repository root:

    python tests/regression/test_quadrupole_operators.py
    python tests/regression/test_octupole_operators.py

The octupole checks evaluate the actual selector expressions and actual hcpx
basis transformation extracted from source. They test all aliases, pure rank,
p-shell vanishing, dipole orthogonality, legacy migration, seven-component
independence and conventional-current matrix algebra for both slots. Source
routing and example namelist checks cover the integration contract. These are
not compiled Fortran or end-to-end Kubo tests.
