# TDDFT-08: longitudinal susceptibility and relaxation (prototype disabled)

## Current production status

`channel = 'longitudinal'` is deliberately **unavailable** in the production
susceptibility driver. It terminates before reading or using a longitudinal
static file. The earlier prototype must be rebuilt on the WR-04 real
static-limit machinery and validated before it can be exposed as a runnable
calculation. In particular, no output from the present code supports an LLB
parameter or `alpha_parallel` claim.

The interface and file format below are retained as a design record for that
future reactivation; they are not a supported input recipe.

## Reserved input and finite-field data

When this route is revalidated, its proposed `&tddft` controls are:

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

## Reactivation requirements

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
