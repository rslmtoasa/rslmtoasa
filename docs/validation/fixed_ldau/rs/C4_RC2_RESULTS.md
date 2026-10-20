# C4 RC2 — RS bulk energy candidate checks

| Parameter | Mean allowance A (Ry/eV) | Candidate b mean error (h=.05/.10) | Candidate c mean error (h=.05/.10) | Verdict |
|---|---:|---:|---:|---|
| U=2 eV, J=0 | 0.004864400 | 0.000231354 | 0.140360848 | PASS |
| U=2 eV, J=0.5 eV | 0.008324456 | 0.000320068 | 0.117123942 | PASS |
| U=2 eV, J=0.9 eV | 0.005689116 | 0.000145018 | 0.106045721 | PASS |

All three RC2 parameter cases pass: mean candidate-b error is within mean A, and mean candidate-c error exceeds 10A. Candidate a remains a record only. At J=.9, candidate b errors at h=.05/.10 are 0.000166078/0.000123958 Ry/eV; candidate c errors are 0.106065459/0.106025982 Ry/eV.

The allowance keeps the spin sum and uses the RC1 total gradient because a spin-split W1 gradient is unavailable: `A_X(h)=2 Σ_s |g_total·dq_s/dX(h)| + |D_X(0.05)-D_X(0.10)|`. The J=.9 raw bundle records seven wrong-lineage steps that were discarded; accepted endpoints restart from the accepted J=.8 checkpoint. The exact C4 campaign CLI was not archived and is identified as such in the sidecar.

Raw paths, SHA256 values, accepted iterations, parsed energy fields, fixed-P evaluator invocation, command template, deviations, and the compact accepted J=.9 trace are recorded in [C4_RC2_PROVENANCE.json](C4_RC2_PROVENANCE.json). An independent audit recomputed the J=.9 derivatives and candidate errors before this update.
