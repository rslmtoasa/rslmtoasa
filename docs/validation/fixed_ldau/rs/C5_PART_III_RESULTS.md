# C5 Part III — Stage RD stopped

**Disposition:** RD1–RD5 were not started because the matched no-U RS standard fixture failed byte identity under C5 §2.5(a). No merge SHA exists. Production remains unchanged; the Co cutoff-array sizing fix is only on the unmerged `fix/ldau-rs-impurity` branch.

| Stage | Status | Decision |
|---|---|---|
| RD1 — recursion-site local moments | `NOT_STARTED_DUE_STOP_2_5_A` | No source or oracle changes. |
| RD2 — frozen Ni replay | `NOT_STARTED_DUE_STOP_2_5_A` | No replay or GATE. |
| RD3 — embedding comparisons | `NOT_STARTED_DUE_STOP_2_5_A` | All comparisons below are `NOT_RUN`. |
| RD4 — impurity energy | `NOT_STARTED_DUE_STOP_2_5_A` | Energy decision `NOT_MEASURED`; no energy adoption. |
| RD5 — adoption | `NOT_STARTED_DUE_STOP_2_5_A` | No adoption or documentation changes. |

| RD3 ID | Comparison | Result |
|---|---|---|
| (a) | RS-bulk Ni → Ni impurity | `NOT_RUN` |
| (a′) | RS-bulk Fe → Fe impurity | `NOT_RUN` |
| (b) | KS-bulk Ni → RS Ni impurity | `NOT_RUN` |
| (b′) | KS-bulk Fe → RS Fe impurity | `NOT_RUN` |
| (c) | KS-local host with spectral impurity mode | `NOT_RUN` |
| (d) | Local versus spectral for B2 FeCo and Ni | `NOT_RUN` |

**R2 partial record:** one Fe rc=100, U=3.8 eV/J=1.1 spectral RS Debug diagnostic from the accepted rc=80 C2 restart exited normally at the 102-step cap but was unconverged. Iteration 101 was wholly finite; the first captured NaN snapshot was iteration 102. This is unaccepted diagnostic evidence, not a completed ladder or accepted R2 result; no root cause is assigned. Raw capture files are linked in [`C5_R2_assets`](C5_R2_assets/); coverage and identity audit details are in [`C5_AUDIT.json`](C5_AUDIT.json).

**Interpretation (3 lines):**

The RS byte-identity stop blocks all downstream RD work. No RD comparisons, impurity-energy measurements, or adoption decisions were made. The Co sizing change remains unmerged and production behavior is unchanged.

Passed: R1 cluster setup and KS no-U identity.
Recorded: Part I issue updates, R1 findings, and the partial R2 diagnostic capture.
Provisional: none; the Fe diagnostic state is unaccepted.
Stopped: RS no-U byte-identity failure under C5 §2.5(a); RD1–RD5 not started.
Next: resolve the byte-identity/determinism discrepancy, then resume C5 RD before RE.
