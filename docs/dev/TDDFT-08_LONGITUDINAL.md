# TDDFT-08: on-site longitudinal susceptibility and relaxation

- [x] Reuse the generalized TDDFT vertex, chi_KS, and Dyson infrastructure.
- [x] Select the `MZ` (`sigma_z`) vertex, which is same-spin in a collinear basis.
- [x] Retain full site-site response matrices and report on-site diagonals.
- [x] Provide a symmetric finite-field static-response driver.
- [x] Require matched `+Delta B` and `-Delta B` data and test linearity as `Delta B -> 0`.
- [x] Calibrate `U_parallel = chi_KS(0)^-1 - chi_parallel(0)^-1`.
- [x] Enhance dynamic response using the existing Dyson solver.
- [x] Fit and report `T_parallel` and `Gamma_parallel=1/T_parallel`.
- [x] Check dynamic Gamma/static finite-field agreement at the configured tolerance.
- [ ] Direct coupled charge-longitudinal XC-derivative kernel (future work).

## Input and finite-field data

Use `post_processing = 'susceptibility'` and the following additions to
`&tddft`:

```fortran
channel = 'longitudinal'
longitudinal_static_file = 'longitudinal_static.dat'
longitudinal_pair_tolerance = 1.0e-10
longitudinal_linearity_tolerance = 5.0e-2
longitudinal_static_agreement_tolerance = 5.0e-2
longitudinal_fit_omega_min = 0.0
longitudinal_fit_omega_max = 0.05
```

The static file is a clean handoff from independently converged field-SCF
jobs. Its first non-comment line is `nsite nrecords`. Each subsequent record
is `perturbed_site signed_DeltaB_Ry mz_site_1 ... mz_site_nsite`. Supply at
least one matched `+DeltaB/-DeltaB` pair for every source site; multiple
field magnitudes are required to make a real linearity check meaningful.
The central differences reconstruct every column of the site-site static
susceptibility.

The longitudinal route is currently restricted to collinear no-SOC input and
requires a Gamma q-point. It writes the ordinary
site-site `chi_KS` and enhanced Dyson outputs plus
`<prefix>_qXXXXXX_longitudinal.dat`, containing response-projector `m0`, static on-site
susceptibility, diagonal and full `U_parallel`, `T_parallel`,
`Gamma_parallel`, fit range/residual, the static-agreement error, linearity
error, and the numerical `eta` used.

The code uses the retarded `omega + DeltaE + i eta` convention. Consequently
the fitted absorptive form is `chi(0)/(1 + i omega T_parallel)`, the
convention-translated version of the requested relaxation model. A failed
static-agreement or linearity condition terminates the run rather than
emitting an apparently valid parameter.

## Eta and LLB caveats

`eta` is written in every longitudinal result. Repeat the same finite-field
calibration for a decreasing eta series and compare the reported
`T_parallel`, fit residual, and static-agreement error; eta remains a
numerical broadening and is not an intrinsic rate. `Gamma_parallel` is a
microscopic inverse relaxation time. It is explicitly **not** an LLB
`alpha_parallel`; an LLB mapping requires a separately specified dynamical
equation and normalization.
