# VAL-18: bcc Fe TDDFT physical validation gate

Date: 2026-08-16

Protocol: [TDDFT-07](TDDFT-07_VALIDATION.md)
Result: **FAIL — Fe is not promoted to a material-validated TDDFT result.**

The campaign used the existing transverse susceptibility route with
`post_processing='susceptibility'`, `goldstone_mode='diagnose'`,
`xi_backend='compare'`, site response projection, no SOC, the production Fe
band window, 300 K electronic smearing, and a fixed 401-point 0--0.08 Ry
frequency grid. TDDFT inputs were frozen before either independent comparison
was inspected. High-accuracy outputs are in
[`results/validation/VAL-18_bccFe`](../../results/validation/VAL-18_bccFe).

## Gate evidence

### Raw Γ response and Ward/Goldstone

Raw diagnostics were retained before any correction. The legacy site-scalar
route gives residuals 0.25699 (`4^3`), 0.0240553 (`8^3`), 0.00506030
(`12^3`), and 0.0189444 (`16^3`), all with unit magnetization overlap. The
non-monotonic mesh behavior prevents calling the Ward quality converged. The
pair-potential raw residual is instead approximately `3e-14` at `12^3` and
`16^3` (unit overlap), clearly separating a successful numerical Goldstone
correction/backend diagnostic from the unresolved legacy Ward quality.

The `8^3` fixed-physics η sweep has the same raw legacy residual 0.0240553
and pair-potential residual `6.47e-16`; η does not repair the Ward defect.

### Small-q dispersion

Five genuinely small direct q points were used: 0.01, 0.015, 0.02, 0.03,
and 0.04. The fixed `16^3` legacy mode file reports accepted individual
fits at 0.0034, 0.0034, 0.0034, 0.0032, and 0.0032 Ry, respectively, with
positive Xi-projected weights and individual fit residuals from 0.0128 to
0.0585. A zero-intercept fit gives

    omega = D |q|^2,  D = 2.91453 Ry,
    relative residual = 0.65317.

This fails the TDDFT-07 quadratic-dispersion criterion. The pair-potential
route is not silently substituted: its five candidates do not yield usable
isolated mode fits (four fail peak isolation/half-height criteria and the
largest-q point has no Xi=1 crossing). In the `12^3` and `16^3` campaigns,
additional pair candidates were likewise rejected rather than included.

### Independent stiffnesses after freeze

The frozen TDDFT inputs were not retuned after comparison. The independent
GBT/frozen-magnon run gives a zero-intercept fit over q=0.01--0.04 of
`D = 0.161163` in its reported energy/q convention, with relative residual
0.05046. Its individual q² slopes are not stable enough to promote as a
converged material reference. The Jij run produced only two pair records,
so no Jij-derived stiffness can honestly be reported. This is an independent
route incompleteness, not a reason to loosen the TDDFT gate.

### Stoner versus collective response

At fixed `16^3`, q=0.01, η=0.001 Ry, the bare KS/Stoner trace is present in
`tddft_disp_nk16_q000001_chi0.dat`. The legacy enhanced trace has a loss
maximum of about 2675.9 at 0.0034 Ry, with Xi crossing near 0.00333 Ry and
positive projected weight about 2676. The pair-potential enhanced trace is
about 2689 at the zero-frequency endpoint; its Xi crossing is near
1.85e-5 Ry and its mode fit is rejected as having no isolated local maximum.
Thus the files contain both Stoner and enhanced-response diagnostics, but
the pair-potential response does not establish a resolved collective magnon
branch under the required extraction rules.

### η sweep

At fixed physics, η = 0.0005, 0.001, and 0.002 Ry were run. Grid FWHM
observations are 0.0014, 0.0022, and 0.0028 Ry; the largest-η result is
endpoint-censored. A finite width at one η is therefore not called intrinsic
Landau damping. The three-point linear extrapolation has intercept 0.0011 Ry
and relative residual 0.0551, but this does not rescue the failed mode or
Goldstone gates.

### Convergence record

The artifact manifest records all TDDFT-07 axes: k mesh, band window,
response projection, electronic smearing, η, and frequency grid. The k-mesh
series was explicitly run. The band window and response projection are
recorded as the production-supported fixed choices; no further material
promotion is attempted after the gate failure. This is why the result is a
failed physical gate, not a claim of full convergence.

## Failure classification

The immediate blockers are:

- **Goldstone/Ward:** the legacy site-scalar raw residual is not mesh-stable
  and is above the material criterion at the fixed `16^3` mesh.
- **Mode extraction / vertex path:** the pair-potential route has near-unity
  Xi behavior but no accepted isolated small-q mode fits, while the legacy
  route gives a non-quadratic dispersion.
- **Transition convergence:** the non-monotonic legacy k-mesh residual and
  incomplete independent Jij record leave the comparison immature.

No evidence here isolates the electronic structure as the primary defect.
No TDDFT parameter was loosened or tuned to force agreement, no one-η width
was called intrinsic damping, and no Fe maturity promotion was made.

## Checklist

- [x] TDDFT-07 read/followed
- [x] raw q=0 response assessed
- [x] Goldstone/Ward convergence assessed
- [x] convergence axes recorded
- [x] small-q data generated
- [x] zero-intercept q² fit performed
- [x] fit residual documented; result rejected
- [x] GBT stiffness compared after freeze
- [x] Jij route executed after freeze; stiffness rejected as incomplete
- [x] Stoner continuum inspected
- [x] Xi/loss diagnostics recorded; collective feature not established
- [x] η sweep performed
- [x] no one-η linewidth called intrinsic damping
- [ ] Fe maturity updated — correctly left unchanged because the gate failed

The compact campaign manifest is
[`campaign.json`](../../results/validation/VAL-18_bccFe/campaign.json). The
standard checker was run against it and correctly failed with
`small-q dispersion is not quadratic within the campaign limit`.
