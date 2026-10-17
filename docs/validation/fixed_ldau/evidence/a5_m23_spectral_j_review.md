# A5.3 spectral J derivative review

Historical checkpoint: the local J campaign had not yet launched when this
spectral-only review was written. The completed paired local/spectral J
results supersede that status; see
[`a5_m23_j_derivative_result.md`](a5_m23_j_derivative_result.md). This
original spectral review remains below as historical evidence.

This review follows the spectral J allowance added to
`a5_m22_numeric_allowance_preregistration.md`; the allowance was recorded
before inspecting the J candidate energy rows. All ten spectral J endpoints
(.4, .45, .5, .55, .6, .8, .85, .9, .95, 1.0 eV) are marked CONVERGED with
return code 0, native KS, finite entropy/TS, and R05 PASS in the campaign
summary. The local J endpoints were not launched after the local U=4
nonconverged case stopped the sequential run.

At J=.5 eV, the spectral fixed-P production target is
`dE_Hub/dJ = −0.063047077414 Ry/eV` (stable over .01/.005/.0025 eV target
steps). The measured h=.05 response-based allowance is nominal `.00295699`
with independent-box interval `[.00204434,.00416579] Ry/eV`. At J=.9 eV the
fixed-P target is `−0.067713496800 Ry/eV`, and the allowance is nominal
`.00081927`, interval `[.00057640,.00113991] Ry/eV`.

For spectral candidate c, centered h=.05 derivatives are `+0.06351074677`
Ry/eV at J=.5 and `+0.03737573425 Ry/eV` at J=.9. The h=.10 derivatives are
`+0.06467197695` and `+0.03775101528 Ry/eV`, respectively. Relative to the
fixed-P targets, c misses by `+0.12655782419` and `+0.10508923105 Ry/eV`, far
larger than the corresponding A5.3 allowance intervals. The directly logged
`t_total_actual_Vin` slopes are `+0.13087452255` and `+0.10900754932 Ry/eV`;
the candidate-c miss differs from those scales by `−0.00431669836` and
`−0.00391831827 Ry/eV`. The first is near the upper edge of its allowance,
before adding the observed finite-step sensitivity; the second is several
times its allowance. These are reported as measured differences, without
adding a new c-slope gate.

For report-only candidate a, h=.05 slopes are `−0.06736377578` at J=.5 and
`−0.07163181507 Ry/eV` at J=.9. Its differences from the fixed-P targets are
`−0.00431669837` and `−0.00391831827 Ry/eV`. Candidate a is not used to claim
an a-versus-b preference; spectral candidate b is unavailable in these rows,
and local J comparisons remain unevaluated.

The candidate a-minus-c slope is `−0.13087452255 Ry/eV` at J=.5 and
`−0.10900754932 Ry/eV` at J=.9, matching the negative of the directly
captured `dt/dJ` slope at the logged precision. Last-ten-row candidate-energy
spans are at most `1.48e-7 Ry` per endpoint; observed h=.05/h=.10 derivative
differences are `3.44e-5 Ry/eV` for a and `1.16e-3 Ry/eV` for c at J=.5, and
`6.18e-5` for a and `3.75e-4 Ry/eV` for c at J=.9. These are observed
finite-step/tail scales, not statistical intervals.

The common case summary is
`/private/tmp/a5-m22m23-remaining-ks-grid/summary.json`, SHA256
`4ce3967ade9b02d27fef40ba245fa212fd608b11bb1d6416ddad55f41f4ebda0`. The
fixed-P target inputs are the direct production P matrices from the J=.5 and
.9 center states, SHA256 `c7dea0dc…cffad06` and `2af9e62e…c58f9a34`; the
driver outputs are `0c614841…920fdbd` and `e27c5e89…ba92037` respectively.
These are spectral-only findings and do not complete the local J comparison.
