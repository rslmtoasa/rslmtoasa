# bcc Fe `exchange_q`

Run from this directory with:

```text
../../../build/bin/rslmto.x input.nml
```

The deck consumes the accepted bulk k-space SCF state and writes
`exchange_q.dat`.  The 12x12x12 k mesh is the electronic Brillouin-zone
integration.  The independent 51-point q path is the magnetic perturbation
path; it is not a second k mesh.  The default path is

```text
q = (0, 0, xi)_direct, 0 <= xi <= 0.5
```

For bcc, `direct` means reciprocal-lattice coordinates in the primitive basis
used by the reciprocal Hamiltonian.  The output also reports the corresponding
Cartesian vector in units of `2*pi/alat` and in inverse Angstroms.  The raw
observable is `DeltaJ(q) = J(Gamma) - J(q)`, with separate torque-torque,
contact, and total finite-H components.  `DeltaJ/q^2` is reported only for
nonzero q and is a diagnostic, not an automatic stiffness fit.

Set `native_turek = .true.` to enable the native screened-LMTO/Turek contour
route. `native_crosscheck = .true.` remains a compatibility alias. The native
columns report ordered `J_ud(q)`, `J_du(q)`, their visible historical average
`J_sym(q)`, raw `DeltaJ_sym(q)`, and the expected complex-amplitude curvature
`2*DeltaJ_sym(q)`; no fitted scale is applied. The contour is controlled by
`native_contour_points`, `native_contour_margin`,
`native_contour_height_fraction`, `native_contour_account_fermi_poles`, and
the optional even `native_contour_target_fermi_poles`. The closure validation
deck uses a fixed 256-pole target: the native coefficient spectrum has a
larger complex span than the compact legacy unit fixture, so this target is
the smallest tested even target that leaves all native poles strictly inside
the production ellipse while retaining explicit residue accounting.

The production finite-H formulation is selected by
`finite_h_spectral_mode = 'metallic'` and inherits the reciprocal SCF
temperature and Fermi level.  For exact commensurate paths, use the compact
validation decks `input_commensurate_12.nml` or `input_commensurate_24.nml`;
their headers record `endpoint_mode = mesh_reuse`.  The original `input.nml`
keeps an off-mesh path to exercise the exact
`endpoint_mode = explicit_diagonalization` route.  The diagnostic
`input_legacy_24.nml` selects the retired occupied-only spectral expression for
the same 24³ q points.
