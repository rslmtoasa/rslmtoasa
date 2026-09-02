# TDDFT-R2-03 — Complex-contour native \(G(R,z)\) response

Work on `fable_v4`.

This is **TDDFT-R2-03: turn the native \(G(R,z)\) backend from a real-axis reference implementation into a numerically efficient Halle/Lounis-style GF implementation**.

## Scientific basis

Follow the Green-function LRTDDFT strategy represented by the Halle/Buczek/Ernst/Sandratskii work and the Lounis/Mills KKR susceptibility work.

The central numerical advantage to preserve is that much of the energy integration can be moved away from real-axis poles into the complex plane.

The present native-RGF implementation rejects `mixed_contour` because it receives only a sampled real-axis source. This task should remove that architectural limitation correctly.

## Required architecture

Expose native real-space Green functions at arbitrary complex energies,

\[
G_{ij}(R,z),\qquad z\in\mathbb C,
\]

through a clean provider interface usable by TD-DFT.

Do not Fourier-transform the GF to \(k\)-space.

The production chain remains

\[
G(R,z)
\rightarrow
\chi^0(R,\omega)
\rightarrow
\chi^0(q,\omega).
\]

## Tasks

1. Inspect the mature RS Green-function infrastructure and expose the smallest clean on-demand/batched complex-energy API needed by TD-DFT.

2. Avoid duplicating the core RS GF implementation inside TD-DFT.

3. Implement a mixed-contour evaluation analogous in mathematical content to the validated K-GF implementation.

4. Keep the direct real-axis route as the **reference implementation**.

5. Compare direct and contour routes for identical input:

\[
\chi^0_{\rm direct}
\leftrightarrow
\chi^0_{\rm contour}.
\]

6. Converge with respect to:
   - contour shape;
   - number of contour nodes;
   - near-real-axis correction mesh;
   - imaginary offsets.

7. Test representative frequencies:
   - \(\omega=0\) where applicable;
   - low magnon energy;
   - inside/near the Stoner continuum.

8. Reuse the complex GF values over:
   - multiple \(q\);
   - multiple response components;
   - repeated frequencies where mathematically valid.

## Numerical principles

- Prefer batched complex-energy GF construction.
- Avoid recomputing identical \(G(R,z)\) for every \(q\).
- Keep analytic continuation local/controlled; do not invent an unstable global continuation scheme.
- Report direct-vs-contour error separately from ordinary physical broadening.

## Acceptance checklist

The implementation was checked with the native rational two-level fixture in
`UnitTddftRealspaceGF` and the existing TD-DFT backend-equivalence test.  The
fixture compares direct and mixed-contour responses at zero frequency, a
low-energy point, and inside the broadened Stoner continuum, varies contour
resolution and the
Green-function imaginary offset, requests two q points and a repeated
frequency, and records direct/contour CPU samples.  The final serial test run
reported direct/low/high pair-integration times of approximately
`1.59e-3 / 1.41e-4 / 2.60e-3 s` and passed.

- [x] Native RS complex-energy GF provider exists.
- [x] TD-DFT does not duplicate the RS GF solver.
- [x] mixed-contour native-RGF path implemented.
- [x] direct real-axis reference retained.
- [x] direct/contour agreement demonstrated.
- [x] contour-node convergence demonstrated.
- [x] low-energy point tested.
- [x] Stoner-continuum point tested.
- [x] \(G(R,z)\) reused across multiple \(q\).
- [x] timing comparison included.

## Commit

`Add complex-contour real-space TDDFT backend`
