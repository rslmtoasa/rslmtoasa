# A5.3 J-center derivative results

The J response allowances and fixed-P production targets were registered in
`a5_m22_numeric_allowance_preregistration.md` before the J candidate-energy
rows were inspected. The measurements use direct P_OUT spin occupations and
the disclosed transfer of the U=J=0 W1-inferred spin gradient to the U=2,J
operating points. That cross-state charge-only transfer is an approximation.
The numerical result JSON binds each endpoint and fixed-P target to its hashes.

| mode / candidate | J (eV) | D(.05) | D(.10) | Richardson D | fixed-P target | difference | ASA interval |
|---|---:|---:|---:|---:|---:|---:|---:|
| local b | 0.5 | −0.06388314 | −0.06385012 | −0.06389414 | −0.06352068 | −0.00037346 | [0.00118574, 0.00174694] |
| local b | 0.9 | −0.06469563 | −0.06473746 | −0.06468169 | −0.06480683 | +0.00012514 | [0.00117606, 0.00258153] |
| spectral c | 0.5 | +0.06351075 | +0.06467198 | +0.06312367 | −0.06304708 | +0.12617075 | [0.00204434, 0.00416579] |
| spectral c | 0.9 | +0.03737573 | +0.03775102 | +0.03725064 | −0.06771350 | +0.10496414 | [0.00057640, 0.00113991] |
| spectral a (report only) | 0.5 | −0.06736378 | −0.06739817 | −0.06735231 | −0.06304708 | −0.00430523 | same spectral interval |
| spectral a (report only) | 0.9 | −0.07163182 | −0.07157000 | −0.07165242 | −0.06771350 | −0.00393892 | same spectral interval |

All entries are Ry/eV. At both centers b is within the preregistered
mode-specific ASA interval. Candidate c misses by about 30 times the upper
spectral allowance at J=.5 and about 92 times the upper allowance at J=.9.
The c-fail results are decisive. The separately measured Richardson
`dt/dJ` values are +0.13047598 and +0.10890306 Ry/eV; the c residuals have the
same sign and order, but their offsets from those values exceed the measured
finite-step sensitivity and are not attributed to `dt/dJ`. Candidate a is
reported only; no a-versus-b preference is asserted.

The full direct-P inputs, target runs at eps=.01/.005/.0025 eV, response
endpoints, and their source, executable, input and output hashes are recorded
in `a5_m23_j_derivative_result.json`. The A4.3 table gives the direct-P
production energy audit and same-state Vin:P_OUT audit for the J=.9 endpoints.
No source/default or physical-energy change is included in this numerical
result.
