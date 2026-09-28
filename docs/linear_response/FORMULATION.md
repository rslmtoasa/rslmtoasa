# Linear-response formulation

## Scope and status

This is the living entry point for the `linear_response` post-processing
route. The supported baseline is scalar-relativistic, collinear, no-SOC,
orthogonal RS-LMTO with the accepted reciprocal `ham_only` state. The route
contains three deliberately different formulations:

| formulation | role | status |
|---|---|---|
| `rotation` | local transverse rotation dynamics | production route |
| `tddft` | spatial transverse TD-DFT | strict-ASA `L=0` target; `L>0` remains research scope |
| `projected` | site-space projected response | model and cross-check, not the spatial response |

The three routes share input parsing and conventions but do not share a
claim of physical equivalence. Longitudinal response, SOC, noncollinearity,
generalized-overlap response, additive operators, and `spdf` response are
outside this baseline.

Sources: `lr-campaign-archive:docs/TDDFT_FORMULATION.md`, `lr-campaign-archive:docs/LR_TDDFT_CONVENTIONS.md`

## Rotation dynamics

For a transverse rotation at `q`, the production object is the retarded
kernel

`K^R(q,omega) = C(q) + Pi^R(q,omega)`,

where `C` is the frequency-independent contact term and `Pi^R` is the
retarded bubble of torque vertices evaluated with the accepted second-order
Hamiltonian `H2`. The rotation coordinates are the two transverse components
per site. The static reduction is the symmetrized zero-frequency kernel and
is compared with the force-theorem Hessian. The `q=0` mode uses the exact
uniform-rotation/Berry identity of this formulation; it does not subtract a
mode or apply a Goldstone correction. Finite-frequency damping is still
retained.

The pole scan uses the configured circular channels. Its default controls are:

| control | default |
|---|---:|
| `rotation_eta_ladder` | `1e-4, 2.5e-5, 6.25e-6` Ry |
| `rotation_probe_omega` | `1e-5` Ry |
| `rotation_probe_eta` | `1e-9` Ry |
| `rotation_slope_step` | `1e-5` Ry |
| `rotation_pole_window_floor` | `1e-3` Ry |
| `rotation_pole_window_scale` | `2.5` |
| `rotation_pole_window_max` | `2e-2` Ry |
| `rotation_pole_coarse_points` | `61` |
| `rotation_pole_fine_points` | `41` |
| `rotation_pole_refinement_half_width` | `2.0` |

The scan compares the positive-real-kernel crossing, minimum `|K|`, and loss
peak. A pole is not promoted when these diagnostics disagree or the eta
ladder is unresolved.

`native_turek` (and its compatibility alias `native_crosscheck`) is optional.
When enabled, the Turek path is evaluated as an independent static diagnostic;
it does not select the circular branch, pole window, pole acceptance, static
normalization, or any Goldstone treatment. The pole-window estimate uses only
the production finite-H/H2 static curvature and the explicit window controls.
With the diagnostic disabled, the rotation output marks the Turek field as
missing rather than writing a fabricated zero.

Sources: `lr-campaign-archive:docs/NATIVE_ROTATION_DYNAMICS.md`, `lr-campaign-archive:docs/TDDFT_FORMULATION.md`

## Spatial TD-DFT

The bare transverse response is the retarded Kohn-Sham susceptibility
`chi_KS` in the canonical response space. The accepted response representations
are `radial_points` and `product_compact`; the product representation is a
rank-revealing compact basis for the same endpoint branches. `response_lmax`
is `-1` for inference and is bounded by the current `0..4` implementation;
for an `spd` basis, the inferred complete product cutoff is `4`.

The direct transverse ALSDA interaction is local in the radial response
space. The enhanced response solves the canonical Dyson equation

`(I - chi_KS K_xc) chi = chi_KS`.

The loss matrix is the retarded anti-Hermitian part
`L = -(chi - chi^dagger)/(2 i pi)`; point-space callers use the response
metric adjoint and compact orthonormal callers use the ordinary dagger. The
explicit `lcmm` interaction is an independent Goldstone-sum-rule route; it
is not silently substituted for ALSDA.

`tddft` accepts `L=0` strict-ASA work as the current target. `L>0` requires
the complete non-spherical response space and remains research scope. The
real-space GF route is a provider-backed alternative bare response, not an
implicit change of the response conventions.

Sources: `lr-campaign-archive:docs/LR_KS_SUSCEPTIBILITY.md`, `lr-campaign-archive:docs/KXC_ALSDA_TRANSVERSE_KERNEL.md`, `lr-campaign-archive:docs/GOLDSTONE_SUMRULE_INTERACTION.md`, `lr-campaign-archive:docs/TDDFT_DYSON_AND_LOSS.md`, `lr-campaign-archive:docs/LR_LMTO_PRODUCT_RESPONSE_BASIS.md`

## Projected response

The projected route forms a site susceptibility `chi0` from the reciprocal
spin-transition contract and fits a local site interaction. `stoner_fit`
uses the projected Mills/Stoner scalarization; `lcmm` solves the projected
site Goldstone equation while keeping the site response fully coupled. The
route reports rank, conditioning, residuals, and eta stability. It is a
controlled site model and cross-check, not a replacement for the spatial
`chi_KS` or the rotation kernel.

Sources: `lr-campaign-archive:docs/DRESP_01_PROJECTED_SITE_SPIN_CONTRACT.md`, `lr-campaign-archive:docs/DRESP_02_PROJECTED_RECIPROCAL_CHI0.md`, `lr-campaign-archive:docs/DRESP_04_PROJECTED_MILLS_RPA.md`, `lr-campaign-archive:docs/DRESP_05_PROJECTED_JUELICH_LCMM.md`

## Explicit boundaries

The longitudinal channel is not implemented. `bare_response='kspace_resolvent'`
is reserved and fails explicitly. Real-space rotation dynamics are not
implemented. The real-space GF provider currently supplies the TD-DFT bare
response seam; it does not turn the auxiliary path operator into a physical
susceptibility input.

Sources: `lr-campaign-archive:docs/TDDFT_PRODUCTION_DRIVER.md`, `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md`, `lr-campaign-archive:docs/TG_FZ_R5_PRODUCTION_H_REPRESENTATION_BRIDGE.md`

## `&linear_response` reference

All values below are the defaults loaded by `linear_response_config`. Energies,
frequencies, broadenings, contour margins, and pole controls are in Ry unless
noted. Counts are dimensionless; `rotation_axis` is dimensionless; q values
are direct reciprocal coordinates unless `q_coordinates='cartesian'`; file
values are paths.

| key | default |
|---|---|
| `formulation` | `'rotation'` |
| `representation` | `'radial_points'` |
| `bare_response` | `'lehmann'` |
| `realspace_solver` | `'auto'` |
| `interaction` | `'alsda'` |
| `projection` | `'spd'` |
| `diagnostics` | `'none'` |
| `channel` | `'chi_plus'` |
| `q_coordinates` | `'direct'` |
| `q_file` | `''` |
| `n_q` | `1` |
| `n_q_points` | `0` (rotation input list) |
| `q_list` | empty; maximum 2000 points |
| `n_omega` | `1` |
| `use_omega_grid` | `.false.` |
| `omega_grid` | zero array; used when enabled |
| `omega_min`, `omega_max` | `0.0`, `0.0` |
| `eta` | `0.01` |
| `n_eta` | `1` |
| `eta_grid` | first value `eta` |
| `response_lmax` | `-1` |
| `write_full_matrix` | `.true.` |
| `rotation_axis` | `[1.0, 0.0, 0.0]` |
| `finite_h_spectral_mode` | `'metallic'` |
| `finite_h_response_backend` | `'spectral'` |
| `contour_points` | `32` |
| `contour_shape` | `'ellipse'` |
| `contour_margin` | `0.25` |
| `contour_height_fraction` | `0.35` |
| `contour_account_fermi_poles` | `.true.` |
| `native_turek` | `.false.` |
| `native_crosscheck` | `.false.` |
| `native_green_eta` | `1e-3` |
| `native_energy_points` | `0` |
| `native_contour_points` | `64` |
| `native_contour_margin` | `0.25` |
| `native_contour_height_fraction` | `0.35` |
| `native_contour_account_fermi_poles` | `.true.` |
| `native_contour_target_fermi_poles` | `0` |
| `native_rsgf_provider` | `'auto'` |
| `gf_integration_points` | `2001` |
| `gf_integration_eta` | `0.0` |
| `gf_energy_margin` | `1.0` |
| `output_file` | `'rotation_dynamics.dat'`; non-rotation implicit default becomes `'tddft_response.dat'` |
| `rotation_eta_ladder` | `[1e-4, 2.5e-5, 6.25e-6]` |
| `rotation_probe_omega` | `1e-5` |
| `rotation_probe_eta` | `1e-9` |
| `rotation_slope_step` | `1e-5` |
| `rotation_pole_window_floor` | `1e-3` |
| `rotation_pole_window_scale` | `2.5` |
| `rotation_pole_window_max` | `2e-2` |
| `rotation_pole_coarse_points` | `61` |
| `rotation_pole_fine_points` | `41` |
| `rotation_pole_refinement_half_width` | `2.0` |
| `rotation_n_omega` | `0` |
| `rotation_omega_min` | `0.0` |
| `rotation_omega_max` | `2e-2` |
| `rotation_eta` | `1e-4` |
| `rotation_grid_file` | `'rotation_response_grid.dat'` |

`gf_integration_eta` is the auxiliary real-axis discontinuity broadening and
is distinct from the physical response `eta`; when used, it is kept smaller
than the physical broadening. `rotation_n_omega=0` selects the pole-control
path rather than a fixed rotation grid.

Sources: `lr-campaign-archive:docs/TDDFT_PRODUCTION_DRIVER.md`, `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md`

## Compatibility and validation

The code accepts only these `tddft` rows:

| representation | bare response | interaction |
|---|---|---|
| `radial_points` | `lehmann` | `alsda` |
| `radial_points` | `lehmann` | `lcmm` |
| `radial_points` | `realspace_gf` | `alsda` |
| `radial_points` | `realspace_gf` | `lcmm` |
| `product_compact` | `lehmann` | `alsda` |
| `product_compact` | `lehmann` | `none` |

For `projected`, the required pair is `radial_points` plus `lehmann`, the
interaction is `none`, `stoner_fit`, or `lcmm`, and `projection` is `d`,
`spd`, or `both`. `rotation` uses its own q-path and pole validation. In all
formulations, `channel` is `chi_plus` or `chi_minus`, `diagnostics` is
`none` or `invariants`, `realspace_solver` is `auto`, `block`, or
`chebyshev`, and `response_lmax` is `-1` through `4`.

The parser rejects the legacy `&tddft` namelist and the legacy post-processing
values `tddft` and `susceptibility`; use `post_processing='linear_response'`
with this namelist.

Sources: `lr-campaign-archive:docs/LR_RESPONSE_BASIS_MAPPING.md`, `lr-campaign-archive:docs/LR_RESPONSE_SPACE_ALGEBRA.md`, `lr-campaign-archive:docs/TDDFT_PRODUCTION_DRIVER.md`
