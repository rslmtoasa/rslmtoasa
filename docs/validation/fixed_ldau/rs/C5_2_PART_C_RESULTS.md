# C5.2 Part C — completed outcomes

C1’s fixed-U Chebyshev source audit is complete: SCF forms deterministic full-matrix g0, P and local-mode M1; no stochastic trace is used.
The original C2 U0 Chebyshev run accepted iteration59/window40–59; the later Debug U0 triage separately accepted iteration76/window57–76. Neither supplies the missing accepted U0 block endpoint.
C2 method comparison is NOT_EVALUATED: the corrected U0 block stopped at the finite B1 capacity guard at iteration36, and U2 block was not run.
Release U2 Chebyshev ended with ENOSPC after Diff56/state57; partial iteration58 lacks P_OUT, POST_SCF_RADIAL and Diff. Debug U2 accepted iteration99/window80–99, then exit2 on final save because Fe_out.nml was read-only.
The clean8b U2 attempt is INVALID_COMPARABILITY because its older restart reader lacks metadata; it does not establish a pre-existing numerical failure.
C3 capture integrity and replay pass from accepted Debug U2: both carrier ratios are within the registered 1% band and both q2′ values are nonnegative.
Co U>0 is NOT_ADMITTED. C5.3 leaves production/defaults unchanged; no fix or adoption is made for the uncovered C2 comparison.
U0 was intentionally SIGTERM’d after endpoint preservation (exit143); owned processes are gone and later live outputs are retained separately as not accepted.

| Item | Recorded result |
|---|---|
| C1 | Deterministic full-matrix source audit complete |
| C2 | NOT_EVALUATED; no accepted U0 block endpoint or U2 block run |
| Accepted states | U0 Cheb59/window40–59; Debug U0 triage76/window57–76; Debug U2 Cheb99/window80–99 |
| C3 | GATE PASS: carrier ratios 0.99325122898615015 / 0.99975940370921512; q2′ 0.013665611639277323 / 0.011285244008227068 Ry² |
| Co U>0 | NOT_ADMITTED |

Passed: C1 source audit; listed accepted windows; C3 capture integrity, carrier and q2′ gates.
Recorded: invalid U0 setups, Release U2 ENOSPC, Debug U2 final-save exit2, clean8b INVALID_COMPARABILITY, and disk-relief/recovery history.
Provisional: none; C2 is NOT_EVALUATED because its accepted block arm is absent.
Stopped: none under §2.5; U0 was intentionally terminated after checkpoint preservation, exit143, with all owned processes gone.
Next: proceed to RD1 Ni/Fe work; omit Co and leave production/defaults unchanged.
