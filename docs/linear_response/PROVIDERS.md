# Linear-response providers

## Provider boundary

The real-space bare-response route is provider-based. A provider supplies a
directed coefficient-space Green-function pair for
`(left_site,right_site,R)` at the requested complex energy:

```text
gij       = G_(left_site,right_site)(z)
gji       = G_(right_site,left_site)(z)
hgamma_ij = h^gamma_(left_site,right_site)
hgamma_ji = h^gamma_(right_site,left_site)
```

The first index is the destination/row and the second is the source/column.
Retarded and advanced arguments are `E+i eta` and `E-i eta` at the provider
boundary. `gji` is an independently supplied directed block; the response
code never manufactures it by transposing `gij`.

Sources: `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md`, `lr-campaign-archive:docs/RSGF_ENDPOINT_AUGMENTATION.md`

## Abstract and concrete providers

The endpoint augmentation seam exposes an abstract effective-Hamiltonian
action provider. Its action applies the production `H_eff` to a localized
source-site seed and returns the destination block. A dense provider is the
small finite-system implementation for this contract. The susceptibility
provider has three distinct coefficient-GF adapters:

| provider | role | reported controls |
|---|---|---|
| dense finite inverse | exact finite-basis oracle | complex energy and directed slices |
| block recursion | native/provider implementation | recursion depth and terminator |
| Chebyshev | native/provider implementation | polynomial order and kernel |

The adapters have separate callbacks and controls. A converged block or
Chebyshev result must approach the reciprocal spectral result while reporting
coefficient-GF error, endpoint-augmentation error, and bubble/integration
error separately.

Sources: `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md`, `lr-campaign-archive:docs/RSGF_ENDPOINT_AUGMENTATION.md`

## Oracle hierarchy

The eigenpair-built Green function is an oracle for a finite Hamiltonian: it
is the same finite spectral resolvent used to check a provider and is not an
independent production route. The dense inverse provider is the separate
finite-matrix implementation used to test directed complex-energy blocks.
Agreement between these two finite constructions does not certify a native
recursion or Chebyshev provider until its controls are converged and its
metadata is recorded.

The production reciprocal Lehmann susceptibility remains the reference for
the provider cross-check. The provider route must use the same state,
endpoint, Pauli projection, response space, q phase, occupation, frequency,
and physical eta as that reference.

Sources: `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md`, `lr-campaign-archive:docs/LR_KS_SUSCEPTIBILITY.md`, `lr-campaign-archive:docs/DRESP_03G_FINITE_H_CONTOUR_GF.md`

## Endpoint augmentation and bubble

Every directed coefficient block enters the physical response only after the
four endpoint branches are assembled with `Phi` and `Phidot`. The onsite and
offsite contact terms from the endpoint contract remain present. The augmented
blocks are then used in the real-axis bubble:

```text
A(E) = i [G^R(E) - G^A(E)] / (2*pi)
G^R(E+omega) with +eta
G^A(E-omega) with -eta
f(E; EF,T)
composite Simpson on [energy_min-margin, energy_max+margin]
```

The periodic pair phase is
`exp(+i 2*pi*q dot (R+tau_right-tau_left))` in direct coordinates. The
assembled raw matrix is converted to canonical `B=chi_raw*W` before Dyson.

Sources: `lr-campaign-archive:docs/RSGF_ENDPOINT_AUGMENTATION.md`, `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md`, `lr-campaign-archive:docs/LR_RESPONSE_BASIS_MAPPING.md`

## H2, H_exact, and the path operator

Rotation dynamics and the accepted linear-response Hamiltonian use
`H2 = E_nu + hbar - hbar*obar*hbar`, the second-order truncation. The exact
LMTO transformed Hamiltonian is
`H_exact = E_nu + hbar*(1+obar*hbar)^(-1)`; its difference from `H2` is a
known truncation distinction and is not silently corrected in the response.
The native Turek/Hamiltonian bridge keeps `H2`, `H_exact`, and any auxiliary
diagnostic separate.

`lmto_path_operator` supplies the auxiliary `g^alpha` path operator for the
Turek exchange trace. It does not reconstruct the physical LMTO coefficient
Green function and is never a `chi0` provider input. A path-operator exchange
trace and a response bubble therefore remain different contracts.

Sources: `lr-campaign-archive:docs/TG_FZ_R5_PRODUCTION_H_REPRESENTATION_BRIDGE.md`, `lr-campaign-archive:docs/TG_FZ_R7_SCREENING_REPRESENTATION_REPAIR.md`, `lr-campaign-archive:docs/TG_FZ_R8R_RESOLVENT_CONTACT_CLOSURE.md`, `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md`

## Requirements for a future F-01 provider

Any provider promoted beyond the finite/provider baseline must:

1. return directed `gij` and `gji` blocks at the same requested complex energy;
2. expose the effective-Hamiltonian provenance needed for `h^gamma`;
3. pass the independent endpoint augmentation oracle, including contact and
   `Phidot` branches;
4. report recursion/Chebyshev controls and convergence separately from the
   response norm;
5. reproduce the reciprocal Lehmann result on the same finite and material
   states at fixed physical eta;
6. preserve the direct-coordinate endpoint phase and all ordered q pairs;
7. carry the accepted no-SOC/collinear/orthogonal capability flags and fail
   closed for unsupported SOC, overlap, noncollinear, additive, or `spdf`
   requests; and
8. keep provider, augmentation, energy-integration, and Dyson residuals in
   output provenance.

These are provider requirements, not a claim that a native production
real-space route is already registered in `calculation.f90`.

Sources: `lr-campaign-archive:docs/TDDFT_RS_GF_BACKEND.md`, `lr-campaign-archive:docs/RSGF_ENDPOINT_AUGMENTATION.md`, `lr-campaign-archive:docs/TG_FZ_R5_PRODUCTION_H_REPRESENTATION_BRIDGE.md`
