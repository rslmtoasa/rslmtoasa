# GBT re-audit WP07: constraining-field energy semantics

## Status and evidence boundary

This report records WP07 on source baseline commit `1c2e8b73` with the WP07
working-tree changes applied.  The focused verification uses the existing CMake
Debug build with GNU Fortran, OpenMP enabled, MPI/libXC/CUDA disabled, and
`RUN_UNIT_TESTS=ON`.

The WP06 deterministic nonzero-field fixture is used for the controller-level
experiment.  Its target is `(0,0,1)`, its canted moments are unit vectors, and
the explicit covariance leg uses `q=(1/4,0,0)` in the repository's Cartesian
`2*pi/alat` convention.  No material constrained-SCF fixture was available in
this work package, so the energy-invariance result below is a fixed-state
controller oracle, not a claim of converged material energetics.

## Code audit

The constraint update is `cfd::constrain` in
`source/include_codes/abspinlib/constrain.f90`.  Its `etcon` output is now
named and exposed as `controller_penalty_energy` in the structured diagnostics.
The legacy `etcon_out` and `self%constraint_energy` interfaces remain for
compatibility, but both are explicitly aliases for this diagnostic quantity.

For a measured moment `m`, target unit vector `e`, and the controller's current
strength `lambda_t`, the implemented branches are:

* method 2: `P2 = lambda_t/2 * |m/|m| - e|^2`;
* method 3: `P3 = lambda_t/2 * |(I - e e^T)(m/|m| - e)|^2`;
* method 4: `P4 = lambda/2 * |m_delta|^2`, where `m_delta` is the PID
  proportional signal, including its first-iteration damping;
* method 5: `P5 = lambda_t * (|m| - m dot e)` after the target is scaled to
  the measured moment magnitude.

These are not one common variational functional.  In particular, method 4
contains integral/derivative controller history, and methods 2, 3, and 5 use
adaptive controller strengths.  Therefore `etcon` is a numerical
penalty/controller merit value.  It is not a contribution to the DFT total
energy and must not be used as a frozen-magnon physical energy shift.

The field is a separate Hamiltonian object.  The code inserts

`V_con = sum_i B_con(i) dot sigma`

once in the local onsite block.  The diagnostics now also report

`field_coupling_energy = sum_i B_con(i) dot m(i) = <V_con>`.

This is the contribution of the auxiliary Pauli operator to an occupied
eigenvalue sum when that operator is present.  It is a field/eigenvalue bookkeeping diagnostic;
it is not added to the physical DFT energy and is not the controller penalty.
The field and moment are evaluated in the active frame, so the same expression
is valid for ordinary calculations and for the primitive rotating GBT frame.

## Physical, auxiliary, and algorithmic quantities

`self::atomsc` computes each atomic physical energy as

`ETOT = EKIN + UTOT + RHOEPS`

from the current density and radial effective potentials.  The ordinary SCF
and report paths now expose the sum of these terms as
`physical_total_energy`.  The constraining field changes the state whose DFT
energy is evaluated, but the field operator itself is not added to `ETOT`.
This is the physical constrained-state energy used by the SCF and frozen-SCF
bookkeeping.

The variational constrained-DFT functional has the schematic form

`L[rho,B] = E_DFT[rho] - sum_i B_i dot (m_i - M_i)`

up to the code's sign convention.  At an exactly constrained stationary state
the residual term vanishes, leaving `E_DFT[rho]`.  This code constrains
directions while allowing magnitudes to relax, and its field is produced by a
controller rather than by an exact multiplier solve.  It therefore does not
represent a separately evaluable Lagrange functional or a Lagrange residual
energy.  `controller_penalty_energy` must not be mislabelled as one.

There is no auxiliary total-energy accumulator in the production SCF path.
The public energy contract is consequently:

| quantity | API/output | meaning | included in physical total? |
| --- | --- | --- | --- |
| physical DFT energy | `self%physical_total_energy`, `Total energy of system` | `sum(potential%etot)` | yes, by definition |
| controller penalty | `self%constraint_penalty_energy`, `constraint_penalty_energy` | current `etcon` merit value | no |
| legacy penalty alias | `self%constraint_energy` | compatibility alias of controller penalty | no |
| field coupling | `self%constraint_field_coupling_energy`, `constraint_field_coupling_energy` | `<V_con>` diagnostic | no |
| raw MFT band energy | `bands%eband` | occupied eigenvalue sum in the selected Hamiltonian | no; it is an MFT observable |

The report retains the historical `Total energy of system` line and the
historical `Constraint energy` line for parser compatibility, while adding
unambiguous physical and diagnostic labels.  New consumers should use the
dedicated names.

### Migration note

The historical `Total energy of system` value changes meaning at the source
level: it is now the physical DFT total and no longer includes the controller
merit value.  The historical `Constraint energy` key remains available as the
controller-penalty alias.  Existing parsers should migrate to
`Physical DFT total energy of system`, `Constraint penalty/controller energy
(diagnostic only)`, and `Constraint-field coupling <V_con> (band diagnostic
only)` so that no diagnostic is mistaken for an energy correction.

## Band-energy / force-theorem conclusion

For a one-shot Hamiltonian containing `V_con`, the raw occupied eigenvalue sum
contains `<V_con>` in addition to the rest of the one-particle sum.  The code
therefore exposes `field_coupling_energy` so a future force-theorem workflow
can make that choice explicit.

The literature distinguishes the observables rather than identifying the
controller penalty with a total energy.  Dederichs-type constrained DFT writes
the constrained functional with a multiplier term and notes that the term
vanishes at the self-consistent constraint, so the physical constrained energy
is the unconstrained DFT functional evaluated at the constrained density.  The
same review gives the alternative path integral of the converged multiplier
with respect to the constrained moment.  See [Constrained DFT for magnetic
systems, Eq. 10--18](https://juser.fz-juelich.de/record/153090/files/FZJ-2014-02765.pdf).

For magnetic force-theorem differences, Bruno identifies the occupied
single-electron sum as the non-self-consistent observable and discusses its
systematic renormalization; see [Bruno, Phys. Rev. Lett. 90, 087205
(2003)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.90.087205).
Jacobsson et al. explicitly define bare MFT as eigenvalue sums without the
constraining field and corrected frozen magnons as fixed-potential one-shot
eigenvalue sums with the converged field added; see [Jacobsson et al.,
Scientific Reports 12, 18987 (2022)](https://www.nature.com/articles/s41598-022-20311-7),
especially Eqs. 4 and the corrected-frozen-magnon procedure following Eq. 6.

Accordingly, WP07 makes the following implementation decision:

1. fully constrained SCF reports `physical_total_energy`, not a raw
   eigenvalue sum and not `physical_total_energy + etcon`;
2. bare MFT continues to use its existing raw band-energy convention;
3. no corrected-MFT field-coupling subtraction or new corrected-MFT mode is
   implemented in WP07; WP08 must choose and test its complete force-theorem
   convention using the exposed `<V_con>` diagnostic and the usual fixed-
   potential double-counting assumptions.

## Numerical invariance experiment

The WP06 canted fixed-state snapshots were evaluated with method 3 after
resetting the controller, first with `lambda_t=5` and then with `lambda_t=50`.
The controller penalty and field coupling changed with the tuning parameter,
as expected, while the fixed-state physical-energy oracle was identical to
machine precision.  The regression test records both facts: tuning-dependent
controller quantities are diagnostics, and no code path adds them to the
physical-energy quantity.

The experiment is a fixed-state controller oracle and deliberately does not
claim that two different finite-step controller trajectories are the same
electronic state.  A full material
invariance run, with two converged constrained-SCF calculations and identical
moments, remains a WP08/WP12 integration task.

## Regression coverage

`UnitGbtWp07ConstraintEnergy` checks:

* zero residual gives zero controller penalty and zero field coupling for all
  supported controller variants;
* the method-3 penalty and `<V_con>` diagnostic obey their direct formulas;
* changing controller strength changes only controller diagnostics for the
  same fixed state;
* the structured diagnostic explicitly reports that no variational Lagrange
  functional is represented.

`GbtWp07ConstraintEnergySource` checks:

* physical SCF/report energy is `sum(potential%etot)`;
* no ordinary or benchmark total-energy path adds the penalty;
* the new fields and labels are wired;
* the existing `constraint_energy` member remains only a compatibility alias.

## Evidence boundary

Established by WP07:

* current `constraint_energy` is classified as a controller penalty/merit term;
* physical DFT energy and auxiliary field coupling are separate named values;
* ordinary and GBT SCF/report paths do not double-add the penalty;
* raw MFT band energy is kept distinct from physical SCF total energy;
* zero-field and fixed-state controller invariants are regression tested.

Not established by WP07:

* a converged material controller-parameter invariance curve;
* the final corrected-MFT double-counting convention;
* a corrected constrained frozen-magnon mode;
* LKAG or fully self-consistent GBT energy closure.
