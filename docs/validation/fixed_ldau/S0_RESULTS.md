# S0 janitor results

Work branch: `fix/ldau-s0-janitor`; base: `90f57c3`.

## T1 line counts

| File | Before | After |
|---|---:|---:|
| `source/self.f90` | 7,567 | 5,876 |
| `source/ldau_diagnostics.f90` | absent | 1,817 |

The moved diagnostics retain the TRACE family and W test-operator path. The closed A2 scan and C2 audit environment hooks were removed from production source.

## Verification

- Release build: PASS (`cmake --build build -j2`).
- LDA+U and local-moment unit tests, auto impurity/surface resolver fixture, and the four compiled SCF B1 guard tests: PASS (12/12 CTest cases).
- Identity baseline and final runs: 10 four-step runs total; all were launched by `tools/run_rslmto_checked.py` with one thread. The two additional four-step runs compared a DM-backed impurity host against an old inline-LDM host, bringing the total to the 12-run cap.
- The +U RS/KS and impurity potential files have byte-identical non-LDM lines. DM real and imaginary values match the old inline matrices bitwise, and each file hash matches its `dm_sha256` field.
- The DM-backed impurity-host and inline-host runs produce identical solver outputs after excluding the changed input host file and the new DM sidecar.
- No-U RS outputs match byte for byte. No-U KS outputs match except for the measured `cpu_s` value in `str.out`.
- Python syntax compilation, GBT negative-fixture JSON parsing, and `git diff --check`: PASS.

## Uncovered outcomes

The 12-run ceiling was reached before runtime fatal probes for a missing DM sidecar and a DM hash mismatch, before rerunning the 13 guard fixtures, and before a fresh merged-tree identity run. Those outcomes remain uncovered; the read path contains explicit fatal checks for both missing files and hash mismatches. The existing provenance decision path was not changed. After the identity runs, the final source received log-only edits to remove repeated +U/+V messages and compact the U/J summary formatting; the Release build and unit tests passed afterward. No solver runs remain within the authorized budget.

The initial KS+U fixture preparation appended a duplicate `&self` namelist, which made preflight reject the deck before the solver launched. The deck was corrected to set `use_kspace` in its existing group, then preflight and the checked run passed. This was a tooling/deck-preparation correction, not a gate result.

## Parser inventory

See [LOG_PARSER_INVENTORY.md](LOG_PARSER_INVENTORY.md) for the script-by-script solver log token list and the preflight wrapper's parsed fields.
