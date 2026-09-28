# Linear-response conventions

## Scope, units, and representations

The certified baseline is scalar-relativistic, collinear, no-SOC, two-spin
channel, orthogonal RS-LMTO in reciprocal `ham_only` mode. Energies,
frequencies, and physical broadenings are in Ry; lengths are in bohr; q and k
use fractional/direct reciprocal coordinates; response coefficients use the
declared radial and angular basis measures. `RHO` is the stored spherical
density `4*pi*r^2*n(r)` and is integrated with `dr`, not with a second `r^2`.

The direct response index is `I=(site,L,M,radial_point,Pauli_channel)`, with
site first. The response cutoff is `Lmax=response_lmax`, inferred as twice
the orbital cutoff when `response_lmax=-1`; the current `sp` and `spd` complete
product cutoffs are `2` and `4`.

Sources: `lr-campaign-archive:docs/LR_TDDFT_CONVENTIONS.md`, `lr-campaign-archive:docs/LR_RESPONSE_BASIS_MAPPING.md`

## Spin, signs, and Fourier phase

The coefficient spinor is ordered `(all orbitals, up; all orbitals, down)`.
The Pauli matrices are the measurement/source operators. `s_z=n_up-n_down`
is number-spin density; the reported positive code moment is `M_code=mu_B*s`,
whereas the physical electron moment has the electron minus sign.

The live lattice phase is

`exp(+i 2*pi*k_direct dot R_direct)`.

For an endpoint from site `a` to `b`, use
`R + tau_b - tau_a`; do not drop the site endpoint phase. The response bubble
uses state `(n,k)` and `(m,k+q)`, occupations `f_nk-f_m,k+q`, denominator
`omega + epsilon_nk - epsilon_m,k+q + i eta`, and measurement vertex first.
The endpoint is folded by a reciprocal-lattice integer only; it is not
replaced by a nearby mesh point.

Sources: `lr-campaign-archive:docs/LR_TDDFT_CONVENTIONS.md`, `lr-campaign-archive:docs/LR_RESPONSE_BASIS_MAPPING.md`, `lr-campaign-archive:docs/LR_KS_SUSCEPTIBILITY.md`

## Circular channels and eta

The halved ladder matrices are `sigma^+ = (sigma_x+i sigma_y)/2` and
`sigma^- = (sigma_x-i sigma_y)/2`; the density combinations use the unhalved
`Sigma^+ = sigma_x+i sigma_y` and `Sigma^- = sigma_x-i sigma_y`. The
`chi_plus` channel measures `s^+` with `s^-` as source and its positive-
frequency transition is occupied up to empty down. `chi_minus` is the reverse.
With halved vertices the circular factor is `2`; with unhalved vertices it is
`1/2` on the corresponding product.

`eta>0` is retarded and belongs in the response denominator. A real-axis GF
provider may also use `integration_eta` for the spectral discontinuity; that
auxiliary value is separate from the physical response eta and must be
smaller. Eta ladders are convergence diagnostics, not pole shifts or repairs.

Sources: `lr-campaign-archive:docs/LR_TDDFT_CONVENTIONS.md`, `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md`, `lr-campaign-archive:docs/TDDFT_DYSON_AND_LOSS.md`

## Response-space algebra

Point-space response vectors carry their radial metric explicitly. A raw
nonlocal operator is converted to the canonical right-weighted operator
`B=A*W`; the Dyson service then composes canonical operators without inserting
another hidden weight. In compact orthonormal coordinates, if `U` contains
the weighted-orthonormal product modes, use
`c=U^H W^(1/2) x` and `x=W^(-1/2) U c`; compact multiplication is ordinary.

The canonical Dyson solve is

`D chi = chi_KS`, with `D=I-chi_KS K`.

The retarded loss is `L=-(chi-chi^dagger)/(2*i*pi)`. In point space the
dagger is the metric adjoint; in compact orthonormal space it is the ordinary
dagger. No inverse, denominator shift, pseudoinverse, or implicit eta is used.

Sources: `lr-campaign-archive:docs/LR_RESPONSE_SPACE_ALGEBRA.md`, `lr-campaign-archive:docs/TDDFT_DYSON_AND_LOSS.md`

## Product basis and endpoint branches

The compact product basis is built from the LMTO endpoint branches and a
weighted SVD. Its complete second-order branch labels are
`00,10,01,11,20,02`; rank and singular-value diagnostics remain part of the
contract. The direct point basis and compact basis are two representations of
the same accepted response space, not independent physical models.

The Pauli/no-SOC endpoint vertex retains the physical `Phi` and `Phidot`
augmentation maps. A GF endpoint pair keeps all four branches:

```text
Phi G Phi^dagger
Phidot (h G) Phi^dagger
Phi (G h) Phidot^dagger
Phidot (h G h) Phidot^dagger
```

Onsite `-I` terms and the `-h^gamma` contact term are retained where required;
an offsite reverse block is supplied independently and is never inferred by
blind transpose.

Sources: `lr-campaign-archive:docs/LR_LMTO_PRODUCT_RESPONSE_BASIS.md`, `lr-campaign-archive:docs/RSGF_ENDPOINT_AUGMENTATION.md`, `lr-campaign-archive:docs/LR_RESPONSE_BASIS_MAPPING.md`

## Pauli transition vertex

For response component `I=(a,L,M,i,mu)`, the transition amplitude is the
angular projection of the Pauli-projected LMTO state with the operator
`Gamma^mu` inserted between the bra at `(n,k)` and ket at `(m,k+q)`. The
radial state is `G + (epsilon-E_nu_work) Gdot`; the direct pointwise vertex
retains its `1/r^2` factor, which cancels only when a later volume measure is
applied. The exact scalar-relativistic density operator is not silently
replaced by this Pauli projection.

The angular basis is complex-harmonic and ordered by `(l,m)` in the live `sp`
or `spd` list. No Hermitian reduction, real-harmonic relabeling, extra `4*pi`,
or unannounced spin factor is permitted.

Sources: `lr-campaign-archive:docs/LR_PAULI_TRANSITION_VERTEX.md`, `lr-campaign-archive:docs/LR_TDDFT_CONVENTIONS.md`, `lr-campaign-archive:docs/LR_RESPONSE_BASIS_MAPPING.md`

## Kernel sign and interaction routes

The ALSDA transverse kernel is derived from the stored up/down XC arrays. Its
spin coefficient is the energy-valued Pauli splitting
`B_xc^sigma=(V_xc,up-V_xc,down)/2`; it is not a field in tesla. The canonical
local interaction has the response-space measure required by the operator
composition. The projected LCMM interaction is solved from its own site Ward
equation and does not inherit a radial `4*pi` by analogy.

The Dyson metadata must identify the selected route and its provenance:
direct ALSDA, Goldstone sum rule, or an explicitly corrected route. A raw
Goldstone residual is reported as a diagnostic; no route silently repairs the
denominator.

Sources: `lr-campaign-archive:docs/KXC_ALSDA_TRANSVERSE_KERNEL.md`, `lr-campaign-archive:docs/GOLDSTONE_SUMRULE_INTERACTION.md`, `lr-campaign-archive:docs/DRESP_05_PROJECTED_JUELICH_LCMM.md`
