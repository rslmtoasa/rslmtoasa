# TDDFT-R2-02 — Exact static native-RGF susceptibility and three-backend Ward test

Work on the current `fable_v4` branch.

This is **TDDFT-R2-02: implement and validate the true static limit of the native \(G(R,z)\) TD-DFT backend**.

This task is required so that all three response backends can participate independently in the Goldstone/Ward test.

## Physics requirement

We need

\[
\chi^{0,+-}(q,0)
\]

as a mathematically controlled static susceptibility.

Do **not** define the static limit by merely evaluating the dynamic routine at

\[
\omega=0,\qquad \eta>0
\]

and calling that static.

The eigenpair backend already has a proper divided-difference treatment, including degenerate/intraband limits. Use that as the numerical reference, not as code to be called internally by the RGF backend.

## Tasks

1. Derive the static real-space GF expression consistent with the existing retarded response convention.

2. Implement it as a genuine native-RGF calculation.

3. Preserve:
   - site/orbital ordering;
   - spin-vertex conventions;
   - \(R\leftrightarrow-R\) orientation;
   - normalization;
   - signed moment conventions.

4. Make `supports_static_limit` true only after this route has been validated.

5. Add backend-equivalence tests:

\[
\chi^0_{\rm eig}(q,0)
\approx
\chi^0_{kGF}(q,0)
\approx
\chi^0_{RGF}(q,0).
\]

6. Test both:
   - \(q=0\);
   - finite \(q\).

7. Then use each backend independently in the raw Ward diagnostic:

\[
r =
[1-\chi^0(0,0)K_{\rm xc}]m_{\rm GS}.
\]

The purpose is **not** to force \(r=0\), but to show that the three implementations produce the same raw consistency diagnostic.

## Important safeguards

- Do not tune \(K_{\rm xc}\) separately for different backends.
- Do not invoke a Goldstone correction to demonstrate equivalence.
- Do not compare only eigenvalues of the Dyson matrix; compare raw complex \(\chi^0\).
- Preserve the distinction between the spin measurement/source vertex and the XC kernel.

## Acceptance checklist

- [x] Static RGF formula documented in comments or developer notes.
- [x] No finite-\(\eta\) dynamic surrogate used as the static definition.
- [x] RGF static path implemented independently.
- [x] `supports_static_limit` enabled only after validation.
- [x] eigenpair vs K-GF static equivalence demonstrated.
- [x] eigenpair vs R-GF static equivalence demonstrated.
- [x] \(q=0\) tested.
- [x] finite \(q\) tested.
- [x] raw Ward residual compared across all three backends.
- [x] no Goldstone correction needed for the equivalence test.

## Completion evidence

The static backend operation is `evaluate_static_grid`; native RGF uses the
retarded/advanced contour identity documented in
[`TDDFT_REALSPACE_GF.md`](../../TDDFT_REALSPACE_GF.md). The equivalence fixture
uses the same unadjusted Ward field for all three backends and compares raw
complex response matrices, not only Dyson eigenvalues.

Validated with:

```text
cmake --build build --target UnitTddftBackendEquivalence -j2
ctest --test-dir build --output-on-failure -R UnitTddftBackendEquivalence
```

The committed evidence in
[`results/tddft_backend_equivalence.json`](../../../results/tddft_backend_equivalence.json)
records native-RGF static matrix errors of `2.917456e-3` at q=0 and
`2.930621e-3` at finite q. The raw static Ward residuals for eigenpair,
K-GF, and native RGF are `0`, `0`, and `2.917456e-3`; no Goldstone correction
is applied.

## Commit

`Add exact static real-space TDDFT susceptibility`
