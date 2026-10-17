# A5.3 independent U-derivative review

This additive review uses the already preregistered fine-grid A5.3 allowance
in `a5_m22_numeric_allowance_preregistration.md`. It does not change the
original or fine-grid preregistrations. Candidate energy rows were inspected
only after that numeric allowance had been recorded.

## Provenance and reductions

The local U=2 endpoint is a verified native-KS, converged, R05-physical state.
Its original receipt remains unchanged; the additive
`/private/tmp/a5-priority-local-u2-ks/u2_local_standard/endpoint_verification.json`
corrects the parser-only `.true.` false negative. It records return code 0,
literal `Converged!0.0000000006`, finite entropy and TS, R05 `physical=.true.`
with streak 9/3, and hashes for the raw endpoint outputs. Its final P matrix
hash is `a07d70e9e86b442555831ab6feead2457b99329cb4d7a6ee9879d4e54283d5cc`.
The production fixed-P driver output hash is
`c1448066fdc8ff7345ab5f6dceb3143faafc60c9d1fc0077ecaff82cc1acd3d7`.
For U steps .01, .005 and .0025 eV, it gives
`dE_Hub/dU=0.0552396890922 Ry/eV`, varying by less than `3e-13 Ry/eV`.

For candidate b, the final logged free-energy values are −2541.865080911355
Ry at U=1.95 eV and −2541.859571132330 Ry at U=2.05 eV. The centered h=.05
derivative is `0.05509779025 Ry/eV`. The h=.10 result from U=1.9/2.1 eV is
`0.05509807945 Ry/eV`, a measured step sensitivity of `2.89e-7 Ry/eV`.
Across the last ten logged endpoint evaluations, the two candidate-b energy
spans are `4.19e-9 Ry` and `1.01e-8 Ry`; their summed span divided by the .1 eV
stencil width is `1.43e-7 Ry/eV`. These are observed tail/step scales, not
statistical confidence limits. The target and candidate-b derivative differ
by `−1.4190e-4 Ry/eV`, which is well within the preregistered local A5.3
allowance interval `[0.00102115,0.00155303] Ry/eV`, even before adding those
small measured derivative scales. Thus b passes the registered local U
comparison at the stated ASA base accuracy.

For spectral candidate c, the final h=.05 derivative is
`−0.08369872376 Ry/eV`; the h=.10 result is `−0.08287859730 Ry/eV`, a step
sensitivity of `8.20e-4 Ry/eV`. The spectral fixed-P target is
`0.057827695387 Ry/eV`, so c's observed miss is `−0.141526419 Ry/eV`, far
outside the spectral allowance interval `[0.00074110,0.00124622] Ry/eV`. The
predeclared `dt/dU` scale is `−0.14862699 Ry/eV`; c's miss differs from that
scale by `0.00710057 Ry/eV` (about 4.8% of its magnitude). The independent
comparison supports a large c miss while reporting that residual explicitly.

The captured spectral `t_total_actual_Vin` values at U=1.95 and 2.05 eV are
`−0.2472708248996794` and `−0.2621335221498775 Ry`, giving
`dt/dU=−0.14862697250 Ry/eV`. The difference between the spectral candidate-a
and candidate-c derivatives is `+0.14862697250 Ry/eV`, matching the negative
of this directly captured slope to the precision of the logged energy rows.
Candidate a is therefore retained as a reported comparison only; this note
does not claim an a-versus-b preference or an adoption decision.

The candidate b endpoint log hashes are
`e627125e…2055a12` (U=1.95) and `fc0080e2…e3fb686` (U=2.05); the candidate
c endpoint hashes are `7124a608…0ef33ab4` and `2ea8cef3…0fcb156d`. The source
campaign summary is
`/private/tmp/a5-priority-ks-main-grid/summary.json`, SHA256
`9b00ccd36bc11fdecc1ba778c4053ea19b265fbb2e2fdac890c44baecb9a81e6`.

J remains a separate allowance and comparison. The spectral J response
stencils are available, with charge-only allowance intervals of approximately
`[0.00191,0.00400] Ry/eV` at J=.5 eV and `[0.000542,0.001096] Ry/eV` at J=.9
eV under the same independent-box propagation. Local J endpoints remain
unrun after the local U=4 nonconverged case stopped the sequential campaign;
no local J allowance or local J candidate conclusion is made here.
