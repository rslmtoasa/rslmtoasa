# C4 RC2 — RS bulk energy candidate checks (partial)

| Parameter | Allowance A (Ry/eV) | Candidate b error, h=0.05 / 0.10 | Candidate c error, h=0.05 / 0.10 | Verdict |
|---|---:|---:|---:|---|
| U=2 eV, J=0 | 0.004864327 / 0.004864472 | 0.000251820 / 0.000210887 | 0.140345658 / 0.140376037 | PASS |
| U=2 eV, J=0.5 eV | 0.008332837 / 0.008316076 | 0.000341329 / 0.000298807 | 0.117133573 / 0.117114311 | PASS |
| U=2 eV, J=0.9 eV | pending | pending | pending | pending |

For the written spin sum, RC1 has a total W1 gradient but not a spin split, so the common `g_total` is applied separately to each spin's RC2 response: `A_X(h)=2 Σ_s |g_total·dq_s/dX(h)| + |D_X(0.05)-D_X(0.10)|`. The auditor confirmed the selected energies and occupation responses use matching accepted endpoints. A collapsed-total response would yield U allowance 0.000122697 and an apparent failure; that interpretation drops the explicit spin sum and was discarded.

The J=0.9 endpoint series is in progress from its accepted local J=0.8 checkpoint. The overall B2 result and RC3 remain pending.
