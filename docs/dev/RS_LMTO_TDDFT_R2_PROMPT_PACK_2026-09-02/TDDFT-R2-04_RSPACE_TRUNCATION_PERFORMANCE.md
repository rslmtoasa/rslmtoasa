# TDDFT-R2-04 — Real-space truncation convergence and actual performance

Work on current `fable_v4`.

This is **TDDFT-R2-04: separate real-space convergence diagnostics from the production \(R_{\max}\) optimization and demonstrate an actual speed benefit**.

## Background

The present implementation computes the pair response for all supplied \(R\)-pairs and only afterwards decides whether each pair belongs inside \(R_{\max}\).

This is excellent for post-hoc truncation diagnostics, but it means \(R_{\max}\) currently does not avoid most of the expensive response integration.

We need two clearly distinct operating modes.

## Mode A — validation/full-tail

Compute the response for all supplied \(R\)-pairs and evaluate diagnostics such as

\[
\epsilon_{\rm tail}
=
\frac{\sum_{R>R_{\max}}|\chi^0(R,\omega)|}
{\sum_R|\chi^0(R,\omega)|}.
\]

This mode is for establishing convergence.

## Mode B — production/truncated

Do not perform the expensive pair response integration for \(R>R_{\max}\).

Use previously established convergence information rather than pretending to know the omitted infinite tail.

## Tasks

1. Separate these two modes cleanly in API/input semantics.

2. Avoid changing the mathematical result in validation mode.

3. In production mode skip expensive pair integrations outside the requested cutoff as early as possible.

4. Do not claim convergence merely because all *supplied* pairs lie inside \(R_{\max}\).

5. Add an explicit distinction between

\[
R_{\max}
\]

used by TD-DFT and the maximum radius

\[
R_{\rm source}
\]

actually available from the underlying GF calculation.

6. Perform nested source-radius calculations:

\[
R_{\rm source}^{(1)}
<
R_{\rm source}^{(2)}
<
R_{\rm source}^{(3)}.
\]

Verify convergence of

\[
\chi^0(q,\omega)
\]

rather than relying only on a shell norm.

7. Benchmark production mode against full-tail mode.

8. Measure:
   - GF time;
   - pair-response integration time;
   - FT time;
   - total TD-DFT backend time.

## Acceptance checklist

- [x] Full-tail validation mode retained.
- [x] Production truncation mode implemented.
- [x] discarded pairs are skipped before expensive integration in production mode.
- [x] `R_source` and `R_max` distinguished.
- [x] no false “converged” claim when the source radius is insufficient.
- [x] nested source-radius convergence demonstrated.
- [x] \(\chi^0(q,\omega)\) convergence reported.
- [x] wall-time speedup demonstrated.
- [x] serial/MPI behavior remains correct.

## Evidence

- `full_tail` remains the default validation mode and evaluates every supplied pair,
  including the omitted tail beyond `R_max`. `production` skips pairs beyond
  `R_max` before pair-response integration and reports that the omitted tail was
  not assessed.
- `realspace_source_rmax` controls the source-pair radius independently of the
  TD-DFT response cutoff `realspace_rmax`. Full-tail diagnostics mark the result
  as unconverged when the supplied source radius does not extend beyond the
  requested response cutoff.
- `UnitTddftRealspaceTruncation` uses 192 synthetic source pairs, three complex
  frequencies, and two q points. It checks exact agreement of the retained
  complex response between full-tail and production modes, verifies 576 versus
  27 pair-response integrations, and performs nested source-radius runs at
  0.875, 3.875, and 23.875 Angstrom. The smallest source is reported
  insufficient; the two larger sources agree at the requested cutoff, while
  the smallest-source response differs from the largest by 2.52%.
- The deterministic benchmark measured 0.214 s for full-tail and 0.0098 s for
  production in the pair-response/FT phase, a 21.9x wall-time speedup. The
  test reports pair-response, FT, GF, and total backend timing fields; the
  synthetic test supplies precomputed Green functions, so its GF phase is zero.
  The production calculation records the actual native GF construction time in
  the TD-DFT metadata.
- Serial checks: `UnitTddftRealspaceTruncation`, `UnitTddftRealspaceGF`,
  `UnitTddftRealspaceMPICorrectness`, `UnitTddftConfig`, and
  `UnitTddftBackendEquivalence` pass. MPI correctness is covered by the
  corresponding MPI test variants.

## Commit

`Optimize and validate real-space TDDFT truncation`
