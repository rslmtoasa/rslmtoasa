# C5.12 results

**Disposition: T0 mismatch; corrected T2 did not pass; no merge.** The C5.11 input Fermi level was `-0.085152 Ry` against accepted KS host EF `-0.082949627451487762 Ry`, a difference of `-0.002202372548512238 Ry`. Its energy window was `[-1.0, 0.5] Ry`, width `1.5 Ry`. Run 1 and continuation input decks had identical SHA-256 and identical energy values.

T1a shifted the Fermi value and both window bounds by `+0.002202372548512238 Ry`, preserving width exactly. The strict diff against the C5.11 candidate contains only `fermi`, `energy_min`, and `energy_max`. Both run preflights compared input EF to the accepted host trace (`/private/tmp/c58-t2a-fixed-release/run/fixed_u_trace.csv`) and passed at zero difference within `1e-6 Ry`.

The two corrected T2 runs used the frozen Release binary and explicit local recipe. Run 1 completed 100 iterations (rc 0); run 2 continued from run 1's exact checkpoint and completed 100 iterations (rc 0). No accepted 20-step window occurred in either run. Run 2's iterations 81–100 had maximum outer Diff `5.71e-8` (limit `3e-8`), spin d-occupation spreads `1.9525748e-7` and `8.2702855e-7` (limit `5e-7` each), and sphere-moment spread `7.6629247e-7 μB` (limit `1e-6 μB`). The Diff and spin-2 spread fail. Earlier windows also failed; all raw run 2 traces, state snapshots, manifest, log, input and launch record are retained in [T2_run2_full_capture.tar.gz](C5_12_assets/T2_run2_full_capture.tar.gz).

The C5.4 preflight now supports a host-trace EF comparison and records the host path and value. Its passing and failing tolerance cases are covered by a focused test. Existing C5.11 T3 documentation remains in place. C5.12 T2 failed, so no merge is authorized; do not adopt the uncovered C5.3 outcome.

Passed: T0 comparison, strict T1a input correction, EF preflight, both run preflights, and two successful solver completions. Recorded: corrected T2 trajectory and failed acceptance windows. Provisional: none. Stopped: C5.12 T2 gate; no merge. Next: root batched audit of this record.
