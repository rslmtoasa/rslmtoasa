# A5 native reciprocal-KS smoke

A short native KS λ=0 preflight was run after the non-KS setup failure was
diagnosed. It used the common RS-origin restart uniformly,
with `&self use_kspace=.true.`, the 12³ reciprocal mesh, 300 K smearing, and
unchanged `conv_thr=1e-9`. The runtime printed both
`run_recursion: use_kspace=.true., skipping recursion solver stage` and
`run_dos: use_kspace=.true. reciprocal SCF branch enabled`; the final total
reported finite entropy and TS, followed by `Converged!0.0000000006`.
This certifies backend routing and that a common RS-origin restart can reach a
KS endpoint. It does not replace A5.3/A5.5 candidate measurements or change
their decision rules.

- Binary SHA256: `52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`
- Source hashes: the five instrumented source hashes listed in
  [`a5_diagnostics.md`](a5_diagnostics.md)
- Pre-run input SHA256: `ebc750a1674d9dbe9a7b28993814fb6fffea6bebc641c8cf7845c3173fffd483`
- Seed `Fe_out.nml` SHA256: `25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5`
- Runtime `Fe_out.nml`/`Fe_scf.nml` SHA256:
  `7f1222f106b6693407476cb11a064cfeac27fe348233cc585fd8d1277fa888b9`
- `stdout.log` SHA256:
  `334d1279756815e1cdc1c7cd17499e6052cb43127c7f9048025be8966c4e2c21`
- `moments_trace.csv` SHA256:
  `a59fd1474a6f86a01e51d64a5a375f1bf6304f50ad0c6a8888456eeb68305d5a`
- Raw run directory:
  `/private/tmp/a5-spin-w1-nativeks-preflight/base_300K_12cubed_lambda_+0.000`
