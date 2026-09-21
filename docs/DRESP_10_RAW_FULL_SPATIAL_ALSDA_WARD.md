# DRESP-10 RAW full-spatial transverse ALSDA Ward

This milestone adds the accepted-state DRESP-10 orchestration for the bcc-Fe
`ham_only`, second-order, 4×4×4 reciprocal snapshot. It is a diagnostic
material gate; it does not run a dynamic spectrum, Goldstone correction,
empirical rescaling, Dyson solve, or mode fit.

The implementation is split into three explicit services:

- `lr_static_frechet_product_response_mod` evaluates the finite-temperature
  Frechet derivative of `f(H)` with the divided-difference `f'(e)` limit for
  degenerate or near-degenerate energies. Its product action is independent of
  the assembled compact susceptibility.
- `lr_full_spatial_alsda_mod` evaluates the full `(site,L,M,r)` source map for
  all six endpoint branches `00,10,01,11,20,02`, split into upper,
  lower-small, lower-rank-0, lower-rank-2, and total pieces. The point-space
  source decomposition is cached before the k-point loop.
- `lr_dresp10_ward_bridge_mod` consumes the accepted reciprocal eigensystem,
  forms the P3 primary target and pointwise `Kxc=Bxc/P3`, checks exact-zero
  active points, compares the raw source to the native rigid-rotation tangent,
  adds the endpoint augmentation contact response, and audits the independent
  static Frechet action.

The production selector is:

```text
backend = 'raw_full_spatial_alsda_ward'
```

It requires one Gamma point, one static `omega=0` point, `response_lmax=4`,
direct ALSDA, no Goldstone correction, and the accepted k-space SCF handoff.

## Fe artifact

The campaign input is
`tests/integration/tddft_driver_smoke/input_dresp10_fe.nml`; the artifact is
`/tmp/dresp10_fe_4k.dat`.

The current run is deliberately classified:

```text
verdict = BLOCKED
classification = RAW_FULL_SPATIAL_WARD_OPEN
```

The source-completeness preflight passes (`response_lmax=4`, product dimension
348, all six branches and all requested SR pieces represented), and the P3
pointwise identity is exact to `1.39e-17`. The material gate then remains open:

```text
field_insertion_relative          = 7.8706583588040424e-01
legacy_l0_field_insertion_relative = 1.6373377259777522e-01
static_frechet_action_relative    = 8.2827917983344587e-02
native_plus_contact_vs_P3         = 9.5871040305000021e-01
ALSDA_plus_contact_vs_P3          = 9.9257265350651058e-01
```

No correction was applied to reduce these residuals. The raw mixed
`K`/augmentation-contact operator and denominator remain explicitly
`NOT_ASSEMBLED`, and dynamic spectra remain `NOT_RUN`.

The next safe continuation is to reconcile the full-spatial L=0 source
normalization/orientation with the certified DRESP-09S field duality, then
resolve the independent product-space Frechet action convention. Only after
those raw gates close should the arbitrary-(L,M) mixed operator and denominator
be assembled.
