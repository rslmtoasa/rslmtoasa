# TDDFT-08: longitudinal susceptibility and relaxation (superseded)

## Current production status

The finite-field calibration prototype described here is no longer the
production route.  TDDFT-13 now exposes `channel = 'longitudinal'` as the
coupled site-major `(charge,m_z)` response, with ground-state-derived local
ALSDA and projected Hartree screening.  See
[`docs/TDDFT_LONGITUDINAL.md`](../TDDFT_LONGITUDINAL.md) for the active route.

The static-field reader, fit types, and unit coverage remain as compatibility
APIs for the TDDFT-08 record.  Their controls are optional and are not read by
the TDDFT-13 production driver.  No output supports an LLB parameter or
`alpha_parallel` claim.

The interface and file format below are retained as a design record for that
future reactivation; they are not a supported input recipe.

## Reserved input and finite-field data

The retained compatibility controls are:

```fortran
channel = 'longitudinal'
longitudinal_static_file = 'longitudinal_static.dat'
longitudinal_pair_tolerance = 1.0e-10
longitudinal_linearity_tolerance = 5.0e-2
longitudinal_static_agreement_tolerance = 5.0e-2
longitudinal_fit_omega_min = 0.0
longitudinal_fit_omega_max = 0.05
```

The reserved static file format is a clean handoff from independently
converged field-SCF jobs. Its first non-comment line is `nsite nrecords`.
Each subsequent record is `perturbed_site signed_DeltaB_Ry mz_site_1 ...
mz_site_nsite`. A future implementation must require at least one matched
`+DeltaB/-DeltaB` pair for every source site and multiple field magnitudes to
make a real linearity check meaningful.

## Historical prototype requirements

Before this route is enabled, it must be restricted to a defined collinear,
no-SOC scope and demonstrate all of the following:

- a real-static, `eta`-independent calibration compatible with WR-04;
- matched finite-field linearity and a documented Gamma q-point convention;
- agreement between static and dynamic limits before writing `U_parallel`,
  `T_parallel`, or `Gamma_parallel`; and
- an explicit mapping and normalization before making any LLB claim.

Eta is a numerical broadening, not an intrinsic rate. Even after
reactivation, a microscopic inverse relaxation time is not automatically an
LLB `alpha_parallel`; that mapping needs a separately specified dynamical
equation and normalization.
