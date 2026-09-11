# RSGF-01R — Native real-space GF susceptibility backend

## Scope and status

The implementation is `source/lr_rs_gf_susceptibility.f90`, in
`lr_rs_gf_susceptibility_mod`.  It is an independent real-space coefficient-GF
route to the same Pauli/no-SOC radial/angular `chiKS` object used by LR-06 and
LR-GF-02.

The backend is intentionally provider-based.  A native block-recursion or
Chebyshev implementation supplies directed coefficient blocks through
`lr_rs_gf_provider%get_pair`; the response code does not call
`green%auxiliary_gij` and does not treat an auxiliary/screened block as a
physical radial GF.  The certified path is:

```text
native coefficient G_ij/G_ji(z)
        -> RSGF-00 endpoint augmentation
        -> LR-03 retarded/advanced GF bubble
        -> real-space pair phase and assembly
        -> LR-04 canonical B = chi_raw * W
```

The dense finite inverse provider is included as an exact finite-system oracle.
The block-recursion and Chebyshev provider types are separate adapters, with
independent controls and callback seams; neither provider is required for the
other to be validated.

Production-driver integration is **not claimed here**.  The main calculation
driver remains reciprocal-only for TDDFT, as specified by TDRUN-01.  A
production material validation through `calculation.f90` therefore remains
outside this task.

The subsequent [`RSGF_CAPABILITY_CLOSURE.md`](RSGF_CAPABILITY_CLOSURE.md)
reconciles this scope with LR-REP-00: R0–R2 are closed for the documented
finite/provider baseline, while production-driver registration is R3 and
remains pending TDRUN-02. This document does not promote finite/backend
evidence to production or Fe/Ni material validation.

## Coefficient-GF contract

For a pair `(left_site,right_site,R)`, the provider returns

```text
gij       = G_(left_site,right_site)(z)
gji       = G_(right_site,left_site)(z)
hgamma_ij = h^gamma_(left_site,right_site)
hgamma_ji = h^gamma_(right_site,left_site)
```

The first index is the row/destination endpoint and the second index is the
column/source endpoint.  `gji` is an independently supplied directed block at
the same complex energy.  The backend never makes an advanced/reverse block by
blindly transposing `gij`; retarded and advanced calls use `z=E+/-i eta` at the
provider boundary.

The exact finite oracle is `lr_rs_dense_gf_provider`.  It inverts the explicit
coefficient Hamiltonian at every requested complex energy and returns both
directed slices.  The recursion and Chebyshev adapters report, respectively,
`recursion_depth`/terminator and `polynomial_order`/kernel in the response
metadata.  Their callback is the intended connection to the native two-sweep
or coefficient-moment implementations.

## Endpoint augmentation correctness

Every directed block entering the bubble is passed through
`augment_lr_gf_endpoint_pair` from RSGF-00.  The four branches are retained:

```text
Phi G Phi^dagger
Phidot (h G) Phi^dagger
Phi (G h) Phidot^dagger
Phidot (h G h) Phidot^dagger
```

The onsite contact terms and the offsite `-h^gamma_ab` term are therefore kept
by the existing certified adapter.  The native backend builds its response
vertex with the physical `Phi`/`Phidot` endpoint maps after this seam; it does
not substitute reciprocal energy moments or a site-only scalar.

The independent endpoint-layer oracle is
`tests/unit/test_lr_gf_endpoint_augmentation.f90` and is documented in
`docs/RSGF_ENDPOINT_AUGMENTATION.md`.  It separately tests dense physical
augmentation, effective-action provenance, spectral endpoint reconstruction,
and negative controls for omitted contact/`Phidot` terms.

## Bubble correctness and energy integration

The native route reuses the LR-GF-02 real-axis contract:

```text
A(E) = i [G^R(E) - G^A(E)] / (2*pi)
G^R(E+omega) with +eta
G^A(E-omega) with -eta
f(E; EF,T) from lr_fermi_dirac_occupation
composite Simpson integration on [energy_min-margin, energy_max+margin]
```

The same two Kubo terms and factor of two are accumulated.  `integration_eta`
is only the spectral-discontinuity broadening and must be smaller than the
physical response `eta`.  No independent occupation, Fermi level, or smearing
convention is introduced.

For every frequency and pair, the backend contracts the augmented directed
branches in the same order as the reciprocal GF bubble.  The result is first a
pointwise/raw matrix and is then converted with the existing
`response_raw_to_canonical` routine.  Its representation string is exactly
`LR-04 canonical right-weighted B=chi_raw*W`.

## Periodic q transform

The live reciprocal convention is imported from
`response_basis_mapping_mod` as `response_fourier_phase_sign=+1`.  The native
pair phase is

\[
  e^{+i2\pi q\cdot(R+\tau_{right}-\tau_{left})},
\]

where `q`, `R`, and `tau` are fractional/direct coordinates.  This is the
endpoint displacement used by the reciprocal assembler; no historical RS
response sign is guessed and basis-site positions are not dropped.

The API accepts a pair list so a provider can expose the complete real-space
translation set.  The focused test includes a nonzero-q, translation-sensitive
two-cell fixture.  It uses a translation-independent dense oracle only to
isolate the phase multiplication; a material periodic block provider must
still supply the physically correct translated `gij/gji` blocks.

## Validation evidence

`tests/unit/test_lr_rs_gf_susceptibility.f90` uses a finite diagonal Hamiltonian
with nonzero `Phi` and `Phidot` radial maps and the complete `L=0..2` response
product space.  It checks:

| layer | comparison | status |
|---|---|---|
| coefficient GF | dense complex-energy inverse versus the same finite spectral resolvent through the reciprocal state | PASS through the exact provider oracle |
| endpoint seam | native blocks enter RSGF-00 before the bubble; independent endpoint oracle remains separate | PASS with the RSGF-00 test |
| bubble | native real-space result versus explicit LR-06 spectral susceptibility | PASS; full response matrix |
| reciprocal cross-check | native result versus LR-GF-02 on the same finite problem | PASS; full response matrix |
| radial/angular space | all site, response `L,M`, radial points, and canonical metric columns are retained | PASS; no site projection |
| q phase | nonzero q and direct-translation phase fixture | PASS |
| energy integration | 201-point Simpson run with `integration_eta=0.001` versus the spectral reference | PASS within the focused oracle tolerance |

The separate `UnitLrGfSusceptibility` test continues to certify LR-GF-02
convergence and its spectral cross-check.  Together, the two tests isolate the
coefficient/augmentation seam from the response bubble rather than reporting
only one final scalar.

## Provider convergence reporting

The backend result metadata reports:

```text
provider kind and provider controls
phase convention
integration_points
integration_eta
eta
fermi_level
temperature
```

For block recursion, vary `recursion_depth` while holding the energy mesh and
`eta` fixed.  For Chebyshev, vary `polynomial_order` under the same condition.
The converged provider result must approach the reciprocal LR-GF-02 result;
coefficient-GF error, endpoint-augmentation error, and bubble/integration error
must be reported separately.  A production callback is not allowed to hide an
unconverged native block behind the final susceptibility norm.

## Claim boundary

The finite-basis claim established here is:

> A directed native coefficient-space GF can be independently augmented and
> assembled into the same finite-basis Pauli radial/angular KS susceptibility
> as the reciprocal response, to controlled coefficient-GF and energy-
> integration error.

This does not validate Fe/Ni magnons, an XC kernel, a Dyson solve, SOC,
generalized overlap, noncollinearity, additive operators, or production-driver
material execution.  Those remain separate gates.
