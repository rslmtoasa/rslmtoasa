# TDDFT backend equivalence

TDDFT-09 adds a deterministic cross-backend campaign for the non-interacting
transverse response `chi0`.  The campaign exercises the three response routes
through their public backend interfaces:

- explicit eigenpair transitions;
- the K-space Lehmann Green-function bubble;
- the native real-space `G(R,z) -> chi0(R,omega) -> Fourier transform` path.

The test is registered as `UnitTddftBackendEquivalence`.  It can be run after
the normal build with:

```text
build/bin/UnitTddftBackendEquivalence results/tddft_backend_equivalence.json
```

The machine-readable output is committed at
[`results/tddft_backend_equivalence.json`](../results/tddft_backend_equivalence.json).

## Fixture and comparisons

The fixture is a one-site periodic spin-split model with two diagonal bands,
deterministic cosine dispersion, circular transverse vertices, zero
temperature, `q = (0, 0.25, 0.5)` along the first reciprocal axis, and
frequencies `omega = (0.24, 0.60)` Ry.  The lower frequency samples the Stoner
support of the fixture; the higher frequency is outside it.  Each run compares
the complete complex 2x2 response matrix, its Frobenius-norm error, and its
eigenvalues.

Exact static response is tested separately for the eigenpair and K-space GF
routes using the analytic zero-frequency divided difference.  The sampled
real-axis native R-GF provider does not claim an exact static operation and
reports that capability as unsupported.  Its `omega=0` finite-eta result is
still included in a separate dynamic Ward diagnostic, so it is not confused
with the exact static limit.

The campaign also measures k-mesh, R-cutoff, energy-resolution, contour-node,
and eta/gamma convergence.  The native R-GF source is built once for the full
q batch; the q transform then reuses the same `chi0(R,omega)` result.

## Evidence summary

The committed run used `nk=4`, 8001 real-axis energy points over `[-4.5,4.5]`
Ry, `eta=0.02` Ry, and `gamma=eta/2` for the main pointwise campaign.

| comparison or diagnostic | maximum reported value | acceptance envelope |
|---|---:|---:|
| eigenpair vs K-GF matrix | `4.4099e-2` | `5.0e-2` |
| eigenpair vs native R-GF matrix | `4.4098e-2` | `5.0e-2` |
| K-GF vs native R-GF matrix | `4.2528e-4` | `1.0e-3` |
| eigenvalue comparison, global matrix norm | `3.2266e-2` | `5.0e-2` |
| exact static eigenpair vs K-GF | `0` | `1.0e-11` |
| exact static Ward residuals | `0` for both supported routes | `1.0e-11` |
| positive-frequency spectral sum residual | `6.8195e-2` | `8.0e-2` |
| contour relative error | `3.1231e-7` | `5.0e-3` |

The pointwise envelope is set by the stable `nk=2,4,8` results (about
`4.4e-2`) and the independently checked eigenvalue norm (about `3.2e-2`),
not by a loose absolute comparison.  The tighter `1e-3` K/R envelope is
supported by the R-cutoff and eta/gamma ladders: the full-source K/R error is
`4.25e-4`, and the eta-ladder values are `5.70e-4`, `3.49e-4`, and `3.42e-4`
for `(eta,gamma)=(0.04,0.02)`, `(0.02,0.01)`, and `(0.01,0.005)` Ry.  The
eta ladder uses 4001, 8001, and 16001 energy points respectively so energy
quadrature does not dominate the broadening comparison.

The finite-grid spectral integral reaches values `3.8407`, `3.7272`, and
`3.7319` against the fixture target `4.0`; these residuals are recorded as a
finite-grid diagnostic, not as a claim that the production Fe/Ni sum rule is
closed.  Likewise, the controlled fixture does not constitute the Fe/Ni
production equivalence gate.

## Negative controls

The harness perturbs the valid K-GF result only in test memory.  Flipping its
sign gives a relative matrix error of `1.9960`, and multiplying it by two gives
`0.9921`; both controls are detected.  This guards against a comparison that
would accidentally accept sign or factor mistakes.

## Scope boundary

TDDFT-09 proves backend equivalence on the deterministic periodic fixture and
records the numerical uncertainty needed for subsequent work.  It does not
close the master blueprint's Fe/Ni production gate, the SOC/noncollinear
extension, or a fully validated exact static operation for sampled native
real-space sources.
