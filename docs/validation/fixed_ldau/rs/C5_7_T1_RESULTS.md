# C5.7 T1 results

**Gate: PASS.** The fixed-U provenance guard now follows the SCF ownership tables. Frozen host modes are checked strictly and their stored provenance is preserved. Updated atoms may convert a restart seed; a missing or different mode logs the requested INFO and the written potential records the selected mode and backend. Legacy mode warnings are limited to frozen U-carrying hosts. Impurity `auto` considers frozen host files when selecting local mode.

The fixture suite passed bulk legacy-seed rewrite, bulk spectral-seed to auto-local rewrite, frozen-local host fatal on forced spectral, legacy frozen-host warning with spectral execution, and surface (`S`) host/recursion ownership checks. Single-thread no-U identity also passed: 27 RS and 24 KS artifacts match the frozen C5.2 outputs after excluding `testrun.log` and normalizing only runtime `cpu_s` fields in `str.out`.

The Release binary is `/private/tmp/c57-t1-build/bin/rslmto.x` (SHA-256 `07a94dd35208ed7ace3b3c729f7d122e5771924265035a57606f7b84a446ea10`). The production source SHA-256 is `33c33aa6d6a49a25dc0015e45a9a30143b83f5e9ebb5b4510485df1dab0240ba`. Every fixture and identity invocation used the C5.4 source-deck/restart preflight, full single-thread caps, pre/post `df` snapshots and a live 1 GiB free-space guard.

Detailed fixture inputs, binary/source bindings, outputs and logs are recorded in [`C5_7_T1_FIXTURE_RESULTS.json`](C5_7_T1_FIXTURE_RESULTS.json); exact identity artifacts and hashes are in [`C5_7_NO_U_IDENTITY.json`](C5_7_NO_U_IDENTITY.json). T1 used zero campaign attempts. T2a is the next authorized card.
