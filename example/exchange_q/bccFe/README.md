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

Set `native_crosscheck = .true.` only after the LKAG q-space oracle has passed
its independent validation and convergence checks; this enables the controlled
metallic native-LKAG columns. For convergence studies, set
`native_green_eta` and a positive `native_energy_points` in `&exchange_q`;
the latter rebuilds the native Simpson energy mesh for that diagnostic only.

The production finite-H formulation is selected by
`finite_h_spectral_mode = 'metallic'` and inherits the reciprocal SCF
temperature and Fermi level.  For exact commensurate paths, use the compact
validation decks `input_commensurate_12.nml` or `input_commensurate_24.nml`;
their headers record `endpoint_mode = mesh_reuse`.  The original `input.nml`
keeps an off-mesh path to exercise the exact
`endpoint_mode = explicit_diagonalization` route.  The diagnostic
`input_legacy_24.nml` selects the retired occupied-only spectral expression for
the same 24³ q points.
