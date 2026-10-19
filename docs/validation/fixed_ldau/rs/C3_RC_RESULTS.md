# C3 Stage RC — RS bulk energy

| Item | Result | Verdict |
|---|---|---|
| RC0: occupations | Zero-temperature step occupations integrated to EF; no Fermi factor or entropy enters RS SCF. Functional: E. | PASS |
| RC1: W1 finite differences | h=.0025 Ry residual −0.158988% (provisional endpoint); h=.005 Ry residual −0.171244% (accepted endpoints). Both measured values are within 1%. | UNESTABLISHED_PROVISIONAL_RECORD |
| RC2: U=2,J=0 derivatives | A=0.000320518270191 Ry/eV. At both h, b is within A and c misses by >10A. | U-family PASS |
| RC2: J=.5 | A=.000340883142186 Ry/eV. h=.05: b error .000341328590308 > A by 4.45448e−7; h=.10 b error .000298807186883 ≤ A. c misses by .11713/.11711 (>10A). | FAIL after one consistency check |
| RC2: J=.9 | 1 accepted local endpoint; remaining cases incomplete when the registered J=.5 failure stopped the campaign. | NOT RUN TO COMPLETION |
| RC3: adoption | Conditional adoption tests not run; no routing change made. | NOT RUN |

**Interpretation**

Part A remains OPEN: the rc140 bulk ended non-finite at step 184, and its impurity record is PROVISIONAL; the 12-run cap is exhausted.
RS SCF uses step-function occupations, so the functional is E rather than E−TS.
RC1’s numeric residuals are below 1%; h=.0025 remains a provisional record and the complete RC1 verdict is unestablished.
The ten U endpoints satisfy the registered b/c derivative checks.
The ten J=.5 endpoints use first-accepted phase0 records and complete fixed-P matrices.
An independent consistency check reproduced the J=.5 h=.05 miss and found no phase, units, matrix, energy, or Hubbard double-count error.
The J=.5 b candidate misses A by 4.45448e−7 Ry/eV; the allowance and formulas were left unchanged.
C1 §5(b) STOP applies after the single honest debugging attempt; no J=.9 conclusion is available.
RC3 was not run because the RC2 gate did not pass.

Passed: RC0 classification and the U-family b/c checks; RC1 measured residual magnitudes are below 1%.
Recorded: all ten U endpoints, all ten J=.5 endpoints, one accepted J=.9 endpoint, and the Part A operator/coefficient evidence.
Provisional: Part A rc140 impurity and RC1 h=.0025; rc140 bulk is FAILED_NONFINITE, and interrupted J=.9 runs are partial only.
Stopped: the J=.5 b gate remains failed after one consistency check; no production energy-routing change was made.
Next: no further RS stage runs under this stopped campaign; Part A finite-cluster scaling and RC3 adoption remain unresolved.
