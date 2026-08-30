# GBT re-audit WP05: rotating-frame / lab-frame contract

## Status and scope

This report records the WP05 implementation on the `fable_v3` source baseline
at `e51ab219` (the working-tree changes were not part of that commit when this
report was written).  WP04 is treated as the preceding fixed-potential gate.
WP05 defines frame semantics and covariance helpers; it does not tune the SCF
loop, radial potentials, or constraining-field controller.

The verification build used CMake Debug mode with GNU Fortran 13.3.0,
OpenMP enabled, MPI/libXC/CUDA disabled, and `RUN_UNIT_TESTS=ON`.  The new
test was run from `build/bin/UnitGbtFrameContract` and through CTest.

The reference problem is scalar-relativistic, SOC-free `gbt_single_q`.  The
unit tests use double precision (`rp`), no Hubbard-V, no SOC, no local-axis
machinery, and no external material fixture.

## Authoritative convention

`gbt_structure_mod` now owns the frame contract through `gbt_frame_t`:

```fortran
type :: gbt_frame_t
   complex(rp) :: U(2,2) ! rotating spinor -> laboratory spinor
   real(rp)    :: R(3,3) ! rotating vector  -> laboratory vector
end type
```

For a sublattice reference direction `(theta, phi)` and cell phase `alpha`,

```text
U = D(alpha) U0,       R = Rz(alpha) R0,
D(alpha) = exp(-i alpha sigma_z / 2),
U0 = D(phi) Ry(theta) D(-phi).
```

The positive phase advances the lab azimuth by `+alpha`.  The phase is
`alpha = q_cart . R_cart`, where `q_cart` is in the repository's Cartesian
`2*pi/alat` convention and `R_cart` is a Cartesian length.  A directed link is
always formed as `U_a^dagger U_b`; `gbt_endpoint_link` now delegates to this
same frame construction.

The GBT primitive potential/reference contract is:

- LMTO potential factors live in the primitive rotating frame.
- The reference axis is cell-independent in that frame.
- `potential%mom` is a collinear rotating-frame reference marker in the GBT
  builder.  Its z sign selects the local up/down channels; its transverse
  components are not a cone-angle input.
- The physical cone and sublattice offsets are supplied by `theta_ss`,
  `theta_ss_sublattice`, and `phi_ss_sublattice`.
- The density accumulated from GBT eigenvectors is rotating-frame density.
  Lab-frame moments are reconstructed only for output or explicit-supercell
  comparison, using the centralized helper.
- In a constrained spiral, transverse density is residual/torque information
  unless a future, explicitly supported relaxed-reference implementation says
  otherwise.

## Central operations

The module provides frame construction for a phase or a physical cell
translation, endpoint links, vector and spinor forward/inverse transforms, and
density-matrix forward/inverse transforms.  Density matrices use
`rho = psi psi^dagger`, so `rho_lab = U rho_rot U^dagger`.  `R` and `U` are
constructed together, rather than allowing the density and Hamiltonian paths to
carry independent rotation conventions.

`spin_density_mod::lab_frame_moment` now calls the central vector transform.
The existing rotating-frame density storage and projection policy are otherwise
unchanged.

## `potential%mom` audit

| Consumer | WP05 interpretation |
| --- | --- |
| `hamiltonian_build::gbt_endpoint_angles` | Reads only the sign of z to select the collinear channel/reference orientation. The new `validate_gbt_reference_moments` boundary rejects transverse values before `build_gbt_bulkham`. |
| `reciprocal_spin_density::accumulate_spin_density_kspace` | Accumulates the complete 2x2 eigenvector density in the primitive rotating spin basis; it does not read `potential%mom`. |
| `self_reciprocal::compute_kspace_spin_moments_spinor` | Returns Cartesian components of the same active eigenvector spin basis. Under GBT these are rotating-frame components; no cell phase is applied. |
| `self::run_dos` / density update | Maps the reciprocal producer's moments into the existing SCF quantities. This downstream full-SCF feedback remains an explicitly deferred WP12 audit; WP05 does not reinterpret it silently. |
| `self::report` | Labels GBT `potential%mom0` output as rotating primitive-frame data; lab phases are not invented for a primitive report. |

An input such as `potential%mom = (mx,my,mz)` with a transverse component above
`1e-10` is rejected by the GBT Hamiltonian boundary with an error explaining
that the physical cone belongs in the theta/phi spiral controls. This makes the
previous silent “read z, ignore x/y” behavior a deliberate contract violation.
The currently available density `relaxed_reference` policy is not treated as a
license to feed a transverse `potential%mom` into the GBT Hamiltonian; its full
operator support remains a later audit item.

## Tests and thresholds

`UnitGbtFrameContract` is registered in the standard CMake unit-test workflow.
It checks, with a predeclared absolute threshold of `2e-12`:

1. SU(2) spinor and SO(3) vector forward/inverse consistency;
2. density-matrix covariance, charge invariance, moment-magnitude invariance,
   and spinor/density agreement;
3. positive phase and sublattice-offset reconstruction for a period-four,
   two-sublattice explicit site table, site by site;
4. frame-link agreement with the legacy public endpoint-link interface and
   reverse-link conjugation;
5. acceptance of collinear rotating references and rejection of transverse
   references.

Observed maximum residual for the new executable: `4.4755e-16`.

The period-four test is the lightweight algebraic supercell oracle: it compares
each reconstructed site against the explicit lab-frame moment

```text
(sin(theta_a) cos(phi_a + q.R_n),
 sin(theta_a) sin(phi_a + q.R_n), cos(theta_a)).
```

The full matched finite-operator spectrum comparison remains the WP04 test;
WP05 adds the moment/frame layer rather than duplicating that operator oracle.

## Evidence boundary

Established by this WP05 change:

- one documented rotating/lab convention;
- mutually paired SU(2)/SO(3) transforms;
- density charge and moment-norm covariance;
- spiral azimuth and sublattice-offset reconstruction;
- explicit `potential%mom` boundary semantics and rejection predicate;
- identification of the reciprocal density producer as rotating-frame data.

Not established here:

- full reciprocal or real-space SCF convergence under GBT;
- constraining-field covariance or nonzero-field convergence;
- physical versus controller constraint-energy bookkeeping;
- corrected constrained frozen magnons or LKAG closure;
- unrestricted `relaxed_reference` GBT operator feedback.
