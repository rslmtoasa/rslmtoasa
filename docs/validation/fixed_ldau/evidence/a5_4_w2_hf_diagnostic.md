# A5.4 W2 Hellmann–Feynman slope diagnostic

This additive note supplements the frozen-state complex-convention result in
[a5_4_w2_complex_convention.md](a5_4_w2_complex_convention.md) and its
unchanged receipt [a5_4_w2_result.json](a5_4_w2_result.json). It records the
finite-difference slope and central zero-field `Tr(W2P)` target. This is a
state diagnostic only; no slope threshold or production gate is applied.

The hook's local external-energy accumulator is zero in these U=J=0 runs, so
reconstruct the free energy as
`F(lambda)=ETOT+lambda*sum_s Tr(W2 P_s)-TS`, using direct `P_OUT` from the
same captured endpoint and the total over both physical spins. At the two
converged reciprocal-KS endpoints:

| λ (Ry) | ETOT (Ry) | TS (Ry) | Σ spin Tr(W2P) | Reconstructed F (Ry) |
|---:|---:|---:|---:|---:|
| −0.01 | −2541.9815158459055 | 9.688935174454278e−5 | 0.3712457147703938 | −2541.985325192405 |
| +0.01 | −2541.981542131514 | 9.030297797852622e−5 | 0.2206149468278093 | −2541.979426285024 |

The central finite difference with half-width 0.01 Ry is
`[F(+0.01)-F(−0.01)]/0.02 = 0.294945369069 Ry/Ry`. The target at λ=0 is
`Tr(W2P(0))=0.295851187793`, with per-spin contributions 0.269057670357 and
0.026793517436. It comes from the compatible A5.2 native-KS 300 K, 12³,
local-moment, U=J=0, `conv_thr=1e-9` λ=0 run. That run used the same initial
seed SHA256
`25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5` and
reciprocal KS backend; its W1 test field has zero coupling, so its saved direct
P is the unperturbed state. The W2 operator was contracted post hoc against
the saved P matrix.

The central difference is lower than the zero-field direct target by
`9.05819e-4 Ry/Ry` (0.306%). The endpoint mean trace, 0.295930330800, is a
separate descriptive approximation, not the λ=0 target. Printed-value and
trace serialization/accumulation contribute at most `1e-10 Ry/Ry` to the
central difference, with a `2e-15` per-endpoint bound for direct `Tr(WP)`.
These capture precision bounds do not cover endpoint SCF sensitivity or the
unquantified O(h²) central-difference truncation. No smaller-step slope was
run, so this mismatch remains a finite-step/response diagnostic rather than
an uncertainty-qualified agreement or disagreement.

The λ=0 target provenance is input SHA256
`ebc750a1674d9dbe9a7b28993814fb6fffea6bebc641c8cf7845c3173fffd483`, P trace
SHA256 `a59fd1474a6f86a01e51d64a5a375f1bf6304f50ad0c6a8888456eeb68305d5a`,
seed SHA256 above, and final restart SHA256
`7f1222f106b6693407476cb11a064cfeac27fe348233cc585fd8d1277fa888b9`.
The machine-readable values and hashes are in
[a5_4_w2_hf_diagnostic.json](a5_4_w2_hf_diagnostic.json).
