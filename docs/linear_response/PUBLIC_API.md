# Linear-response public API

Production uses accepted scalar-relativistic, collinear, no-SOC, orthogonal
second-order `ham_only` data. Unsupported spin structures fail at the state
capability boundary; unsupported methods and combinations fail at input or
prepared-request validation. No selector substitutes a different physics method.

## Canonical production matrix

| Physical identity | formulation | representation | bare_response | interaction | projection / capability |
|---|---|---|---|---|---|
| Rotation effective-action dynamics | rotation | defaults (response axes do not apply) | defaults | defaults | local transverse rotation; one-site automatic pole workflow |
| JUELICH-d — LITERATURE REFERENCE | projected | radial_points | lehmann | lcmm | d, chi_plus; frozen-EF normalized projector; Gamma; ≥4 decreasing positive eta values |
| MILLS-1U — CONTROLLED RS-LMTO MODEL PROJECTION | projected | radial_points | lehmann | mills_1u | d interaction, full spd propagation; native center splitting and coefficient moment; no fit |
| Transverse spherical direct ALSDA | tddft | radial_points | lehmann or realspace_gf | alsda | explicit response_lmax=0 |
| Transverse compact direct ALSDA | tddft | product_compact | lehmann | alsda | complete retained six-branch sp/spd product space |

The two ALSDA representations and the registered spherical GF provider are
numerical realizations of direct transverse ALSDA, not additional literature
methods. Existing provider certification is bounded; native material convergence
is not implied by finite/provider integration tests.

Rotation outputs K=C+Pi and its inverse rotation-coordinate propagator. The
static adiabatic MFT observable maps to KL; finite-frequency spin susceptibility
amplitudes are not identified. Neither full KL nor full-Halle general spatial
spin-charge TDDFT is offered.

## Selector inventory, before and after LR-METHOD-04

The pre-edit inventory included rotation plus nine response combinations: four
radial rows (Lehmann/GF × ALSDA/LCMM), two compact rows (ALSDA/none), and three
projected interactions (none/Mills-1U/LCMM), with projection cross-products.

| Route / selector | Before | After | Reason |
|---|---|---|---|
| Juelich-d | production, deeper gates | production, early gates | closed frozen-EF d identity |
| Juelich spd/both | parser accepted; driver rejected | explicitly rejected before dispatch | dynamics uncertified |
| Juelich spdf | generic projection error | explicit Juelich certification error | dynamics uncertified |
| mills_1u | production, d-only | production, d-only | controlled unfitted interaction |
| stoner_fit | rejected | rejected | removed fitted U; no automatic migration |
| mills_fit/projected_mills/projected_stoner/scalarized_mills | no canonical production route | explicit rejection | no alternate fitted/projected Mills method |
| rotation | production | production | effective-action kernel |
| kl/kl_dynamic/katsnelson/lichtenstein | unsupported formulations | explicit full-KL rejection | not implemented, never aliases for rotation |
| Radial ALSDA | production baseline | strict spherical production baseline | surviving direct ALSDA scope |
| Radial L>0 | exposed research deck | rejected; deck/test deleted | uncertified full-spatial semantics |
| Compact ALSDA | production | production | complete accepted compact product response |
| Radial LCMM + Lehmann | executable dynamics | rejected; dynamic driver branch deleted | SR density violated Pauli response contract |
| Radial LCMM + realspace_gf | executable dynamics | rejected; dynamic driver branch deleted | same contract mismatch |
| Dynamic GSR/corrected ALSDA Dyson labels | accepted by low-level requests | explicitly rejected; route constants/test removed | no unvalidated dynamic interaction |
| Compact none | bare validation | bare validation | transition/bubble oracle, no interacting dynamics |
| Projected none | bare validation | bare validation | d/spd/both site-operator and bubble checks |
| Static radial GSR | independent service | validation-only service | static Ward/sum-rule algebra |
| Compact GSR | independent service | validation-only service | independent product-space static Ward action |
| native_crosscheck=true | alias for native_turek | rejected; use native_turek | one static-oracle selector |
| native_rsgf_provider non-default | overlapped realspace_solver | rejected; use realspace_solver | one GF-provider selector |
| block_recursion input spelling | provider alias | rejected; use block | canonical spelling |
| legacy_occupied | executable alternate curvature | rejected; dispatch deleted | superseded production curvature |
| contour-only rotation backend | alternative curvature | rejected; use both for validation | oracle cannot replace production curvature |
| rotation_slope_step | ignored deprecated input | inert transition input | no independent physics; validates canonical probe-step independence |
| plus/minus channel aliases | accepted by low-level response requests | rejected; use chi_plus/chi_minus | canonical transverse spelling |
| kspace_resolvent | reserved/rejected | rejected | not implemented |
| old &tddft / susceptibility dispatch | rejected | rejected | canonical &linear_response interface |

`formulation` is rotation/tddft/projected. Public `backend` does not exist:
runtime backend labels are derived privately from canonical combinations.
`representation`, `bare_response`, `interaction` and `projection` are constrained
by the matrix above. Bare projected validation accepts d/spd/both only.
`realspace_solver` is auto/block/chebyshev and applies only to realspace_gf.
`channel` is chi_plus/chi_minus, with Juelich-d plus-only; longitudinal/charge
are unsupported. `diagnostics` is none/invariants and cannot choose physics.
`native_turek` is a rotation static-reference flag. Removed native aliases
accept only inert defaults for explanatory input transition, never execution.

For rotation, `finite_h_spectral_mode` is metallic only and
`finite_h_response_backend` is spectral or both; both adds the finite-H contour
oracle while keeping the spectral production curvature. Contour shape is ellipse.
All methods use explicit q/frequency/broadening grids and existing capability
limits. Juelich-d and Mills-1U accept eta ladders; other response rows require
one eta per run. `native_green_eta` and `native_energy_points` are inert old
input fields; native static validation uses the named contour controls only. No numerical tolerances or reference values were changed for this cleanup.

## Retained non-production execution and services

| Path / service | Independent purpose for surviving production |
|---|---|
| native_turek static reference | Turek/LKAG MFT cross-check of rotation curvature; never sets pole branch/window or production Hessian |
| finite-H resolvent contour (both) | validates spectral rotation Hessian against independent resolvent construction |
| zero-temperature finite-H occupied-projector routines in unit fixtures | independent analytic static/finite-difference checks of rotation derivatives; not a runtime production mode |
| Hamiltonian rotation / finite-difference oracles | validate torque/contact algebra and static Hessian |
| diagnostics=invariants for rotation | validates static reduction, Berry/circular/covariance/causality checks and the existing bounded Fe campaign |
| diagnostics=invariants for compact ALSDA | validates static denominator and interacting exact q/-q covariance |
| compact bare none / direct transition band loops | validate compact Pauli transition vertices and bare Lehmann accumulation without Kxc/Dyson |
| projected bare none / direct band-loop tests | validate projected d/spd site operators, moments and bare bubbles, not an alternate Juelich projector |
| dense inverse, eigenpair GF and native block/Chebyshev fixtures | validate directed coefficient GF and endpoint augmentation against independent finite resolvents/Lehmann sums |
| evaluate_lr_goldstone_sumrule / solve_lr_goldstone_equation | static radial Pauli Ward equation, metric-aware action and rank/residual tests; future full-spatial certification infrastructure |
| evaluate_compact_goldstone_sumrule | UnitLrCompactGsrAction compares independent point-space actions, Gamma columns and projected static Ward action |
| manual finite-GF quadrature / mixed-eigenvector audit fixtures | report spectral moments and GF/Lehmann differences across integration ladders; independent finite bubble/normalization validation, not production methods (two existing CTest entries remain disabled) |
| direct ALSDA and small complex Dyson fixtures | validate Kxc normalization, matrix solve, loss and covariance; known-answer fixtures use direct-ALSDA provenance; no dynamic GSR/corrected label remains |
| frozen-EF Juelich projector/direct-loop fixtures | validate Juelich-d normalization, moments, static U and site Dyson independently |
| conditioning, low-m and response-rank reports | measure solve/basis quality without changing any interaction or response |

The existing accepted-real-space-potential adapter builds the reciprocal
reference eigensystem at the accepted EF for spherical ALSDA and bare validation
fixtures. It is shared state preparation required by their existing tests, not a
second interaction or an abandoned response equation. Compact ALSDA, Mills-1U
and Juelich-d require the accepted reciprocal SCF state directly.

## Deferred physics

Full KL finite-frequency chi^{+-}: not implemented. Juelich spd/spdf: static
aggregate identity derived, dynamic closure uncertified. Mills full Coulomb
tensor: literature parent, not implemented. Longitudinal/charge, SOC/noncollinear
transverse, full-Halle spin-charge, non-spherical radial ALSDA and radial LCMM
dynamics require future certification. No half-live runtime selector is retained.

Static-sum-rule U followed by dynamic Dyson is legitimate LCMM theory. The
removed radial implementation was not certified because its SR density input
did not belong to the Pauli susceptibility representation. Replacing that input
alone would not certify a repaired route.
