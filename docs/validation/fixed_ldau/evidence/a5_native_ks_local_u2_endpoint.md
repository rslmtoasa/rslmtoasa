# Native-KS local U=2 endpoint verification

This is the corrected local U=2, J=0, 12³, 300 K, W1, λ=0 case from the
preregistered M2.2/A4.3 inventory. It was run from the common RS-origin
restart and satisfies the native-KS branch, finite entropy/TS, ordinary SCF
convergence, and active-U R05 physical-streak checks. It is reused as the
local center case; no duplicate local U=2 run is needed. This records endpoint
status and provenance only. Candidate energies were not compared or
interpreted.

The original receipt is preserved unchanged. Its `r05_latest_physical=false`
field was a parser defect: the raw source line records Fortran
`physical=.true.`. The supplemental verification record reads the raw line
and retains the original receipt hash.

- Executable SHA256: `52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`
- Pre-run input SHA256: `7f0615d2ab190fe0c1db6fff46829be5e8944690d30c12653ace61bb6794c8c3`
- Common seed SHA256: `25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5`
- `stdout.log` SHA256: `ab6c91d7570183fcfe200d196994173f60e77474f244e2ddd185ef16e9b527db`
- `records.csv` SHA256: `a4b1d904f39d035fc427b41405e184b4ac0bd492704f31b0b6b19cfa16396fc6`
- Final restart SHA256: `bb647254a5612f53d00e2745026c2cde32511f4363a1c7c3a1331d06573ead6b`
- Original receipt SHA256: `0de110c953e25743a4c976b6dcc7cd5853d9b8ce2ac05c8ba08d7ca0bb2f308b`
- Supplemental raw-log verification SHA256: `e21914e74ae04c46f495ab0340f6bb13e76ee65970cb668054652a5f6a5662d8`
- Raw case directory: `/private/tmp/a5-priority-local-u2-ks/u2_local_standard`

The latest raw residual reports `physical=.true.`, streak `9` of required `3`.
The executable prints `Converged!0.0000000006`, a native reciprocal-KS branch
marker, and finite entropy and TS. The supplemental verification record
contains the source and output hashes and the exact residual token.
