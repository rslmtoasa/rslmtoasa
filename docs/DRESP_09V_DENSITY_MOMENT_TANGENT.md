# DRESP-09V: production density-moment tangent

DRESP-09V closes the density-side representation response using the live
reciprocal LMTO energy-moment contract. It does not introduce a generic
orthogonalization matrix `X`.

## Representation

For every accepted k point, the equilibrium coefficient-space density is
constructed from the accepted eigenvalues, eigenvectors, Fermi level,
temperature, and k weights:

```text
rho = f(H)
M0  = rho
M1  = H rho
M2  = H^2 rho
```

The production `reciprocal_spin_density` path stores these matrices in
`spin_density%rho(:,:,l,site,order)`, with orders 1, 2, and 3 corresponding to
`M0`, `M1`, and `M2`. DRESP-09V uses the same site/orbital grouping (s, p, d,
and f where available) and retains the complete complex 2x2 spin block,
including transverse entries, before any local-axis projection.

## Fixed-H versus complete endpoint derivative

The old DRESP-09S diagnostic remains unchanged. Its density-side derivative
holds the endpoint Hamiltonian powers fixed:

```text
delta rho
    |
    v
H^p delta_rho H^q
```

DRESP-09V differentiates the physical endpoint object itself:

```text
H(theta), rho(theta)
    |
    v
D_pq(theta) = H(theta)^p rho(theta) H(theta)^q
    |
    v
delta D_pq
    |
    v
SR radial density
```

For branch order `(00, 10, 01, 11)`, the complete tangent is

```text
delta D00 = delta rho
delta D10 = delta H rho + H delta rho
delta D01 = delta rho H + rho delta H
delta D11 = delta H rho H + H delta rho H + H rho delta H
```

The routine `sr_l0_density_tangent_from_hamiltonian` contracts these complete
branches through the certified DRESP-09S scalar-relativistic bilinears. The
old `sr_l0_density_from_matrix_hamiltonian` remains the fixed-H comparison
oracle.

## Moment and endpoint oracles

The isolated `lr_lmto_density_moment_tangent` service provides two independent
coefficient-space constructions:

- product-rule tangents of `M0`, `M1`, and `M2`;
- divided-difference Fréchet tangents of `g_p(H)=H^p f(H)`.

Equal and near-degenerate eigenvalues use the analytic derivative limit. A
global rigid spin rotation independently checks
`delta M_p = -i [G, M_p]` for all three moments. At equilibrium the endpoint
branches reduce to `M0`, `M1`, `M1`, and `M2`; on the rotated path their
tangents reproduce the corresponding moment tangents.

The finite-dimensional unit fixture uses noncommuting orbital blocks with
unequal spin splitting, finite temperature, a near-degenerate pair, a
finite-angle oracle, production s/p block accumulation, and a spin-degenerate
zero-response control.

## Fe result

The central integration case is the accepted bcc-Fe 4x4x4, 64-k-point,
300-K, second-order `ham_only`, spd state. The coefficient and production
moment gates close:

```text
[H,rho] residual                         2.04034e-12
M0/M1/M2 product-vs-Frechet              1.14792e-14 / 5.90531e-14 / 3.54560e-13
rigid delta-M0/M1/M2 commutators         0 / 1.02330e-12 / 3.04973e-12
production M0/M1/M2 finite difference   1.04167e-6 / 1.04167e-6 / 1.04167e-6
```

The historical DRESP-09S fixed-H residual is `4.08210e-2`. The complete
endpoint tangent has a nonzero endpoint-H contribution norm of `1.83417e-1`,
and the current complete radial residual is `3.05325e-1`, with per-l residuals
`s=9.40638e-1`, `p=1.36148`, and `d=3.05037e-1`. Thus the representation
moment tangent is closed, while the independently reconstructed radial target
remains open. No fitted correction is applied.

The fixed-H, endpoint-H, and total field/density pairings close independently
to the reported numerical precision; the total signed pairing is
`-1.14781e-1`. Core bookkeeping remains separate: the diagnostic core
integral is `-1.92670e-5` and is not used to reduce the valence residual.

The resulting primary classification is:

```text
DENSITY_MOMENTS_CLOSED_RADIAL_TARGET_OPEN
```

The former DRESP-09T density correction is explicitly retired. ALSDA Ward,
Dyson, spectra, and BES/Halle remain outside this milestone (`NOT RUN` and
`OFF`, respectively).
