# Linear-response formulation

## Scope and status

This is the living entry point for the `linear_response` post-processing
route. The supported baseline is scalar-relativistic, collinear, no-SOC,
orthogonal RS-LMTO with the accepted reciprocal `ham_only` state. The route
contains three deliberately different formulations:

| formulation | role | status |
|---|---|---|
| `rotation` | local transverse rotation / electronic effective-action dynamics | production route |
| `tddft` | transverse direct ALSDA | strict spherical radial points or complete compact product response |
| `projected` | Juelich-d or Mills-1U site-space response | certified within their distinct projection contracts |

The three routes share input parsing and conventions but do not share a
claim of physical equivalence. Longitudinal response, SOC, noncollinearity,
generalized-overlap response, additive operators, and `spdf` response are
outside this baseline.

Sources: `lr-campaign-archive:docs/TDDFT_FORMULATION.md`, `lr-campaign-archive:docs/LR_TDDFT_CONVENTIONS.md`

## Rotation dynamics

For a transverse rotation at `q`, the production object is the retarded
local-rotation effective-action kernel

`K^R(q,omega) = C(q) + Pi^R(q,omega)`,

where `C` is the frequency-independent contact term and `Pi^R` is the
retarded bubble of torque vertices evaluated with the accepted second-order
Hamiltonian `H2`. The rotation coordinates are the two transverse components
per site. The static reduction is the symmetrized zero-frequency kernel and
is compared with the force-theorem Hessian. The `q=0` mode uses the exact
uniform-rotation/Berry identity of this formulation; it does not subtract a
mode or apply a Goldstone correction. Finite-frequency damping is still
retained.

The dynamic rotation kernel is not identified with the full KL transverse
susceptibility. Its inverse is a rotation-coordinate propagator; physical
spin-susceptibility amplitudes require a separate source-coupling derivation.
Its adiabatic limit maps to KL MFT: `H_theta_theta = 2[J(0)-J(q)]`,
`B_Berry = M_band/2`, and `omega = H_theta_theta/B_Berry = 4[J(0)-J(q)]/M_band`
in the positive-moment convention. See [KL_IDENTITY.md](KL_IDENTITY.md) for
the equation mapping and signed-channel/reporting conventions.

The lower-level response evaluator and `K(q,omega)` kernel are dimensioned
for two transverse coordinates per magnetic site (`2*Nsite`). The automatic
production pole workflow is intentionally narrower: its scalar circular
channel reduction, `finite_h(1,1)` static scale, and single-magnetization
reporting are one-site only. A multisite request fails at this workflow
boundary with an explicit capability error; multisite callers must consume
the lower-level `K(q,omega)` response API rather than infer a scalar pole.

The production capability gate requires `reciprocal_mode='ham_only'`,
`kspace_ham_order='second'`, and `hamiltonian%hoh=.true.` so that the pole
workflow cannot enter with a first-order or inactive-HOH state.

The pole scan uses the configured circular channels. Its default controls are:

| control | default |
|---|---:|
| `rotation_eta_ladder` | `1e-4, 2.5e-5, 6.25e-6` Ry |
| `rotation_probe_omega` | `1e-5` Ry |
| `rotation_probe_eta` | `1e-9` Ry |
| `rotation_slope_step` | deprecated compatibility input; ignored |
| `rotation_pole_window_floor` | `1e-3` Ry |
| `rotation_pole_window_scale` | `2.5` |
| `rotation_pole_window_max` | `2e-2` Ry |
| `rotation_pole_coarse_points` | `61` |
| `rotation_pole_fine_points` | `41` |
| `rotation_pole_refinement_half_width` | `2.0` |

The scan compares the positive-real-kernel crossing, minimum `|K|`, and the peak
in `-Im K^{-1}`. A collective local-rotation pole is not promoted when these
diagnostics disagree or the eta ladder is unresolved.

The q=0 Berry slope uses one canonical finite-difference step,
`h=rotation_probe_omega`, for both `K(+h)`/`K(-h)` sampling and the
`2h` denominator. `rotation_slope_step` remains accepted for input
compatibility but has no independent numerical role; differing values are
reported as deprecated and ignored.

### Signed circular-channel semantics

The implemented transverse convention is

`theta_plus = (theta_x - i theta_y)/sqrt(2)` and
`theta_minus = (theta_x + i theta_y)/sqrt(2)`.

For the accepted q=0 state, the independently checked Berry slopes obey
`b_plus = d K_plus/d omega = -B` and
`b_minus = d K_minus/d omega = +B`, where `B` is the signed Berry
commutator, with `B = M_band/2`. Therefore the linearized roots are

`omega_plus = K_plus(q,0)/B` and
`omega_minus = -K_minus(q,0)/B`.

If the static transverse curvature is a common signed value `kappa`, changing
the sign of `kappa` exchanges which circular channel has the positive-
frequency root; the other root is its negative-frequency partner. A negative
static curvature is consequently not, by itself, a proof of an instability.
The signed inverse rotation kernel spectral weight and the retarded pole
diagnostic must also be checked. With
`rotation_spectral_weight = -Im(K**(-1))`, the current causality indicator is
`(d Re K/d omega) * Im K > 0` at the resolved real-axis crossing. Rotation spectral-weight
signs may be opposite in the two channels because the slopes are opposite; taking
an absolute value would erase this convention information.

The focused `UnitLrRotationProductionAdapter` audit evaluates both channels,
their signed linear roots, a signed real-axis crossing, and the rotation
spectral-weight/causality indicator. It also checks the Cartesian covariance
`K_AB(q,omega) = conj(K_AB(-q,-omega))`. The coarse 1x1x1 q=0.05 fixture
provides the following signed example from the existing response-frequency
grid (Ry units, `eta=1e-5`; actual roots are linear interpolations of the
signed `Re K` crossings):

| channel | near-static `Re K(q,0)` | slope | predicted root | actual root | `-Im K^{-1}` sign | slope*`Im K` |
|---|---:|---:|---:|---:|---:|---:|
| `+` | `-2.4011e-3` | `-2.0057` | `-1.1972e-3` | `-1.2052e-3` | negative | positive |
| `-` | `-4.4150e-4` | `+2.0084` | `+2.1982e-4` | `+2.2010e-4` | positive | positive |

The production positive pole is the `-` channel (`+2.9945 meV` at the
smallest configured eta), while the `+` channel has the negative-frequency
partner. This is a chirality/negative-frequency-partner interpretation, not
an absolute-value repair and not a branch-selection rewrite. The simultaneous
coarse-mesh Cartesian q/-q/-omega covariance residual was `6.25e-17`.

`native_turek` is the optional independent static reference. The removed
`native_crosscheck` alias fails explicitly when enabled.
When enabled, the Turek path is evaluated as an independent static diagnostic;
it does not select the circular branch, pole window, pole acceptance, static
normalization, or any Goldstone treatment. The pole-window estimate uses only
the production finite-H/H2 static curvature and the explicit window controls.
With the diagnostic disabled, the rotation output marks the Turek field as
missing rather than writing a fabricated zero.

The rotation workflow is generic when `diagnostics='none'`: it validates the
accepted state, evaluates the q path and dynamic kernel, performs static
reduction and rotation-mode pole/spectral-weight analysis, and writes the
response output without Fe campaign assumptions. `diagnostics='invariants'` enables the optional
historical Fe validation layer. `native_turek` is orthogonal to that choice:
it enables the independent Turek diagnostic but does not enable the Fe
campaign. The one-site, CCOR-off, 300 K, PASS-A/PASS-B, and small-q fit
checks are not requirements of the generic workflow.

Sources: `lr-campaign-archive:docs/NATIVE_ROTATION_DYNAMICS.md`, `lr-campaign-archive:docs/TDDFT_FORMULATION.md`

## Transverse direct-ALSDA TDDFT

The physical local ALSDA derivative uses the accepted scalar-relativistic
SCF number-spin density `m_SR=n_up^SR-n_down^SR` and the same functional's
Pauli XC coefficient `bxc_SR=(Vxc_up^SR-Vxc_down^SR)/2`:
`Kxc_SR=bxc_SR/m_SR`. The electron magnetic moment carries the usual electron
sign; response variables here are number-spin/Pauli coefficients. The factor
2 associated with halved sigma-plus/minus response vertices remains separate
from Kxc. Small finite densities are divided directly, with diagnostics and
no floor; an active zero density fails closed. A null-measure origin extension
cannot change a contraction.

`m_P` is the separately reconstructed Pauli response magnetization.
**m_P does not redefine Kxc.** Radial canonical metric operators and compact
product operators represent this one physical kernel. Compact projection is
`U^H K_point U` in weighted orthonormal coordinates. The raw diagnostic uses
an independently projected accepted SCF XC field: `R_Ward=chiKS*bxc_SR-m_P`,
with metric L2, relative L2, maximum residual and normalized response overlap.
An eta ladder at Gamma and zero frequency tests the limiting behavior;
finite eta is a conditioning diagnostic, not the exact static identity.
Goldstone correction = OFF; these residuals never alter any production input.


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
spatial production interaction is `alsda` only. The static radial GSR service
is independent Ward/sum-rule validation infrastructure; it is not a production
interaction selector. Dynamic radial LCMM is rejected before dispatch: its
former driver supplied scalar-relativistic radial density where the response
contract requires accepted Pauli magnetization. A repair and certification
would require a separate physics campaign. Static-sum-rule U followed by
dynamic Dyson is legitimate theory; this removal concerns the implementation
contract, not that theory.

The radial-point production route requires explicit `response_lmax=0`.
Its non-spherical `L>0` extension is rejected. Compact ALSDA requires the
complete retained six-branch product space of the accepted sp/spd basis.
Neither route claims general full-Halle spin-charge TDDFT. The existing
real-space GF provider supplies the strict spherical radial ALSDA bare
response; `realspace_solver` is its sole public provider selector.

Sources: `lr-campaign-archive:docs/LR_KS_SUSCEPTIBILITY.md`, `lr-campaign-archive:docs/KXC_ALSDA_TRANSVERSE_KERNEL.md`, `lr-campaign-archive:docs/GOLDSTONE_SUMRULE_INTERACTION.md`, `lr-campaign-archive:docs/TDDFT_DYSON_AND_LOSS.md`, `lr-campaign-archive:docs/LR_LMTO_PRODUCT_RESPONSE_BASIS.md`

## Projected response

The projected route offers `none`, `mills_1u`, and `lcmm`. `mills_1u` uses
the accepted transformed d-shell center splitting and coefficient-space d
moment to define one physical site interaction, then evaluates the local d
spin-flip bubble over the full accepted spd bands. It is classified as
`MILLS-1U — CONTROLLED RS-LMTO MODEL PROJECTION`; see
[`MILLS_IDENTITY.md`](MILLS_IDENTITY.md) for its equations, normalization,
and model-reduction residual. It requires `projection='d'`, an accepted
converged reciprocal state, scalar-relativistic collinear `ham_only` data,
and the full spd basis. It does not alter the Jülich `lcmm` route.

The former `stoner_fit` keyword fails fast as removed. The projected site
response remains a controlled model and is not a replacement for spatial
`chi_KS` or the rotation kernel.

`projected + lcmm + projection='d' + channel='chi_plus'` is the certified
**JUELICH-d — LITERATURE REFERENCE**. It uses the frozen-EF normalized d
projector and a decreasing ladder of at least four positive eta values to
construct the static eta→0 interaction, followed by site-space Dyson.
Gamma is required in the q list. `spd`, `spdf`, and `both` are rejected at
configuration validation because their dynamic generalization is not certified.
See [JUELICH_IDENTITY.md](JUELICH_IDENTITY.md).

`interaction='none'` is a bare-response validation seam: projected d/spd/both
checks validate site operators and bubbles, and compact bare checks validate
transition vertices. Neither executes an interacting production method.

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
| `native_crosscheck` | removed; `.true.` rejected; `.false.` is inert input transition |
| `native_green_eta` | `1e-3` |
| `native_energy_points` | `0` |
| `native_contour_points` | `64` |
| `native_contour_margin` | `0.25` |
| `native_contour_height_fraction` | `0.35` |
| `native_contour_account_fermi_poles` | `.true.` |
| `native_contour_target_fermi_poles` | `0` |
| `native_rsgf_provider` | removed selector; only inert default `'auto'` accepted; use `realspace_solver` |
| `gf_integration_points` | `2001` |
| `gf_integration_eta` | `0.0` |
| `gf_energy_margin` | `1.0` |
| `output_file` | `'rotation_dynamics.dat'`; non-rotation implicit default becomes `'tddft_response.dat'` |
| `rotation_eta_ladder` | `[1e-4, 2.5e-5, 6.25e-6]` |
| `rotation_probe_omega` | `1e-5` |
| `rotation_probe_eta` | `1e-9` |
| `rotation_slope_step` | deprecated compatibility input; ignored |
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
| `radial_points` | `realspace_gf` | `alsda` |
| `product_compact` | `lehmann` | `alsda` |
| `product_compact` | `lehmann` | `none` |

For `projected`, the required pair is `radial_points` plus `lehmann`, the
interaction is `none`, `mills_1u`, or `lcmm`, and `projection` is `d`,
`spd`, or `both` only for the bare `none` seam. Both interacting projected
routes require `projection='d'`; Juelich-d also requires `chi_plus` and its
static eta ladder. `rotation` uses its
own q-path and pole validation. In all
transverse formulations, `channel` is `chi_plus` or `chi_minus` (Juelich-d is plus-only), `diagnostics` is
`none` or `invariants`, `realspace_solver` is `auto`, `block`, or
`chebyshev`, and `response_lmax` is `-1` through `4`.

The parser rejects the legacy `&tddft` namelist and the legacy post-processing
values `tddft` and `susceptibility`; use `post_processing='linear_response'`
with this namelist.

Sources: `lr-campaign-archive:docs/LR_RESPONSE_BASIS_MAPPING.md`, `lr-campaign-archive:docs/LR_RESPONSE_SPACE_ALGEBRA.md`, `lr-campaign-archive:docs/TDDFT_PRODUCTION_DRIVER.md`

## Deferred capabilities

Full KL finite-frequency chi^{+-} is not implemented; rotation is not its
alias. Juelich spd/spdf have a static aggregate sum-rule identity but no
certified dynamic closure. The Mills full Coulomb tensor is the literature
parent and is not implemented. Longitudinal, charge, SOC/noncollinear response,
general full-Halle spatial response, and radial LCMM dynamics await separate
certification. These are roadmap statements, not runnable selectors.

The complete selector and validation inventory is in [PUBLIC_API.md](PUBLIC_API.md).
