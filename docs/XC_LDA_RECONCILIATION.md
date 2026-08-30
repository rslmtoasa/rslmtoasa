# XCR-02 — Reconcile legacy LDA kernels against the XCR-01 oracle

Date: 2026-08-30

Binary: `UnitXcLdaReconciliation`, GNU Fortran 13.3.0, libXC 5.2.3, serial CPU

## Decision

XCR-01 was reviewed before changing production code. It ruled out units,
spin ordering, energy normalization, and exchange as explanations for the
BH discrepancy. It did identify two demonstrated legacy defects:

1. the default Barth–Hedin potential reconstruction was not the spin
   derivative of its own reported energy at finite polarization; and
2. `XCPOT` discarded every finite-density state with one zero spin channel.

The first defect was repaired by differentiating the implemented BH energy
expression. The second was repaired by guarding only the total-density-zero
case; the zero-channel exchange loops in the legacy PBE/LAG helpers were also
made safe for their exact zero contribution. No libXC call replaced a legacy
functional, and `NEWRHO`, `ATOMSC`, radial solvers, and LMTO potential
parameter routines were not changed.

## Pre-edit XCR-01 decision table

The following is the evidence table before the source correction. Values are
maximum absolute differences over the original regular grid, in Ry per
electron for `epsilon_xc` and the potentials.

| Pair | epsilon_xc | Vxc up | Vxc down | Bxc / Stoner | Classification | Action |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| TXC=1 / TXC=101 BH/VBH | 8.83e-9 | 1.55e-2 | 1.55e-2 | 4.87e-9 | LEGACY_DEFECT_FIXED | Repair the common charge-channel derivative and polarized guard |
| TXC=5 / TXC=105 PW/PW | 6.62e-7 | 8.73e-7 | 1.48e-5 | 7.85e-6 | EXPECTED_VARIANT_DIFFERENCE | No formula replacement; retain historical PBE-LDA path |
| TXC=4 / TXC=106 VWN/VWN | 5.08e-4 | 8.36e-4 | 1.07e-2 | 5.35e-3 | EXPECTED_VARIANT_DIFFERENCE | No formula replacement; retain historical VWN path |

The pre-edit BH energy and `B_xc` already agreed closely, while both spin
potentials acquired almost the same offset as polarization increased. The
XCR-01 finite-difference audit of `rho*EXC` localized that offset to the
legacy `ARS/BRS` reconstruction. The old exact fully polarized rows were all
zero because of the one-channel guard and were therefore not valid XC values.

## Conventions and exchange gate

The oracle uses

\[
n=n_\uparrow+n_\downarrow=\frac{3}{4\pi r_s^3},\qquad
\zeta=\frac{n_\uparrow-n_\downarrow}{n},
\]

with `RHO1 = n_down`, `RHO2 = n_up`, `V1 = V_down`, and `V2 = V_up`.
`EXC` is energy per electron; potentials and energies are reported in Ry.
The production libXC wrapper converts native Hartree results by two.

The independent spin-LDA exchange reference is

\[
\epsilon_x=-\frac{3}{2}\left(\frac{6}{\pi}\right)^{1/3}
\frac{n_\uparrow^{4/3}+n_\downarrow^{4/3}}{n},
\qquad
v_{x,\sigma}=-2\left(\frac{6}{\pi}\right)^{1/3}n_\sigma^{1/3}.
\]

The post-repair exchange validation remains at maximum absolute error
`1.55e-15 Ry` and maximum relative error `1.69e-14`. This validates the
sign, factor of two, density convention, and spin ordering. The small
relative maximum is caused by the near-fully-polarized test point and remains
well below the `2e-12` gate.

## Barth–Hedin equivalence audit

The legacy TXC=1 constructor cites von Barth and Hedin, J. Phys. C 5, 1629
(1972), and uses

| Constant | TXC=1 value |
| --- | ---: |
| `XCCP` | 0.0504 |
| `XCCF` | 0.0254 |
| `XCRP` | 30 |
| `XCRF` | 75 |
| interpolation exponent | `FTH = 4/3` |
| interpolation constants | `AA = 2^(-1/3)`, `BB = 1-AA` |

The legacy exchange term is `EPSXP = -0.91633059/rs`. The source’s
polarization variable is `X = RHO1/RHO = (1-zeta)/2`, so it is the
spin-down fraction. Its energy interpolation is the standard normalized
paramagnetic/ferromagnetic interpolation in that variable.

libXC identifies native ID 17 as `XC_LDA_C_VBH`, also citing von Barth and
Hedin. The native ID is documented in the [libXC functional list](https://libxc.gitlab.io/functionals/)
and in the [libXC 5.2.3 function-ID header](https://sources.debian.org/src/libxc/5.2.3-3.1/src/xc_funcs.h).
The original reference is [von Barth and Hedin](https://doi.org/10.1088/0022-3719/5/13/012).

Thus TXC=1 and TXC=101 are the correct same-family comparison. They are not
bitwise identical: the retained historical constants are rounded legacy
values and the independent implementations evaluate the parametrization
differently. The residual energy difference, about `8.9e-9 Ry`, is therefore
classified as `RECONCILED_EQUIVALENT` after the derivative defect is removed.
TXC=3 and TXC=11 use distinct BH-family constants and remain variants, not
aliases of TXC=1.

## Derivative audit and exact correction

Write `x = RHO1/RHO`, `k = 5.1297628`, and

\[
F(x)=\frac{x^{4/3}+(1-x)^{4/3}-2^{-1/3}}{1-2^{-1/3}}.
\]

The legacy energy expression is

\[
\epsilon=\epsilon_x^0+\epsilon_c^P
 +F(x)\frac{CNY+FTH\,\epsilon_x^0}{k},
\qquad CNY=k(\epsilon_c^F-\epsilon_c^P).
\]

For any radial component define

\[
U_a=\epsilon_a-\frac{r_s}{3}\frac{d\epsilon_a}{dr_s}.
\]

For `q = rs/XCRP` or `rs/XCRF`, the source evaluates the needed
`q F'(q)` as

\[
3q^3\ln(1+1/q)-\frac{1+q^3}{1+q}+\frac q2-2q^2.
\]

The exact derivatives of `rho*EXC` are therefore

\[
V_0=\frac43\epsilon_x^0+U_P
 +F\left(U_F-U_P+\frac{FTH}{k}\frac43\epsilon_x^0\right),
\]

\[
S=F'(x)\frac{CNY+FTH\,\epsilon_x^0}{k},qquad
V_\downarrow=V_0+(1-x)S,\qquad V_\uparrow=V_0-xS.
\]

The source correction implements exactly these expressions in the default
BH-family branch. It includes both the total-density derivative and the
polarization derivative. The exchange and energy-interpolation definitions
are unchanged.

At `(rs,zeta)=(0.5,0.99)`, the before/after values are:

| Quantity | Before | After |
| --- | ---: | ---: |
| `Vxc up` | -3.18584800795 | -3.20130354121 |
| `Vxc down` | -1.02235528135 | -1.03781081356 |

The after values agree with symmetric finite differences of the reported
legacy energy to below `1e-8 Ry`; the old common residual was about
`1.55e-2 Ry`.

## Post-repair fixed-density results

The rerun contains 378 rows: 351 regular points and nine exact fully
polarized rows per pair. The regular grid is
`rs = 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6 bohr` and
`zeta = 0, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.25, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999999`.

| Pair | max |Delta epsilon_xc| | max |Delta V_up| | max |Delta V_down| | max |Delta B_xc| | Classification |
| --- | ---: | ---: | ---: | ---: | --- |
| TXC=1 / TXC=101 BH/VBH | 8.88e-9 | 1.18e-8 | 9.16e-9 | 6.46e-9 | RECONCILED_EQUIVALENT |
| TXC=5 / TXC=105 PW/PW | 7.44e-7 | 8.94e-7 | 1.70e-5 | 8.95e-6 | EXPECTED_VARIANT_DIFFERENCE |
| TXC=4 / TXC=106 VWN/VWN | 5.08e-4 | 8.36e-4 | 1.17e-2 | 5.83e-3 | EXPECTED_VARIANT_DIFFERENCE |

For BH/VBH, the largest regular potential difference is now at the numerical
or parameter-rounding level and no longer grows as a common charge offset.
The PW and VWN differences are retained as expected variant differences; no
historical formula was forced to match libXC.

## Small-polarization/Stoner response

The oracle computes the signed magnetic derivative

\[
I_{xc}=\frac{[V_\downarrow-V_\uparrow](+\delta)
 -[V_\downarrow-V_\uparrow](-\delta)}{2n\delta}
\]

at `delta = 1e-2, 1e-3, 1e-4, 1e-5, 1e-6`. The converged `delta=1e-4`
relative differences are unchanged in character by the charge-channel
repair:

| Pair | representative `rs` | Legacy | libXC | Relative difference |
| --- | ---: | ---: | ---: | ---: |
| BH/VBH | 2 | 9.08014137 | 9.08014131 | 6.36e-9 |
| BH/VBH | 6 | 64.89597433 | 64.89597381 | 8.01e-9 |
| PW/PW | 2 | 9.48136335 | 9.48135618 | 7.57e-7 |
| PW/PW | 6 | 57.18527416 | 57.18522030 | 9.42e-7 |
| VWN/VWN | 2 | 9.45401827 | 9.61937341 | 1.72e-2 |
| VWN/VWN | 6 | 57.73400406 | 60.87529921 | 5.16e-2 |

The BH magnetic derivative was already correct because the old common
charge-channel error cancels in `B_xc`. The repair makes the scalar and spin
potentials consistent with the same energy without altering that magnetic
result.

## Polarized-density limits

The old guard

```fortran
if (RHO1 < TOL .or. RHO2 < TOL) then ...
```

was changed to a total-density guard. For a finite total density, the source
now evaluates `zeta = 0, 0.9, 0.99, 0.999999, 1` directly. At exact `zeta=1`,
the vanishing spin contribution is evaluated through its finite
`n_sigma^(1/3)` limit. No artificial density floor is introduced in the
legacy path.

At `rs=2`, the corrected BH endpoint is
`EXC=-0.6618345150 Ry`, `V_down=-0.3503034631 Ry`, and
`V_up=-0.8623961854 Ry`. The libXC wrapper’s comparison endpoint still uses
its existing `1e-20` positive floor for the zero channel, so its endpoint
potential differs by up to `1.17e-5 Ry`; this is a wrapper evaluation
convention, not a legacy discontinuity. The direct legacy regression confirms
finite values for TXC=1, 3, 4, 5, and 11 at all five requested polarizations.

## Legacy functional status

| Legacy selector | Status | Evidence/action |
| ---: | --- | --- |
| TXC=1 Barth–Hedin | LEGACY_DEFECT_FIXED; RECONCILED_EQUIVALENT to TXC=101 after repair | Exact energy derivative and BH/VBH grid agree; full polarized limit finite |
| TXC=2 Slater X-alpha | RECONCILED_EQUIVALENT for its analytic exchange form | No formula change; common polarized guard now does not discard a zero channel |
| TXC=3 Barth–Hedin–Janak | LEGACY_DEFECT_FIXED for the shared BH-family derivative; EXPECTED_VARIANT_DIFFERENCE from VBH | Same source defect fixed with its own constants; no identity with TXC=1 claimed |
| TXC=4 VWN | EXPECTED_VARIANT_DIFFERENCE to native VWN | Energy-gradient audit is at `4.47e-9 Ry`; no formula replacement; polarized guard fixed |
| TXC=5 PBE 1996 LDA | EXPECTED_VARIANT_DIFFERENCE to native PW | Energy-gradient audit is at `7.42e-10 Ry`; no formula replacement; zero-channel helper path fixed |
| TXC=6 Wigner | OUT_OF_SCOPE_BY_LEGACY_CONTRACT | Spin-polarized construction is intentionally unsupported |
| TXC=7 Perdew–Zunger | OUT_OF_SCOPE_BY_LEGACY_CONTRACT | Spin-polarized construction is intentionally unsupported |
| TXC=11 ASW BH variant | LEGACY_DEFECT_FIXED for the shared BH-family derivative; EXPECTED_VARIANT_DIFFERENCE from VBH | Same source defect fixed with its own constants |

TXC=8 and TXC=9 are GGA paths rather than legacy LDA kernels. No claim is
made about their general GGA behavior by this task.

## Artifacts and tests

The complete machine-readable evidence is in:

* [xc_lda_reconciliation.csv](../tests/xc_reconciliation/results/xc_lda_reconciliation.csv)
* [xc_lda_energy_gradient.csv](../tests/xc_reconciliation/results/xc_lda_energy_gradient.csv)
* [xc_lda_stoner_response.csv](../tests/xc_reconciliation/results/xc_lda_stoner_response.csv)
* [xc_lda_exchange_validation.csv](../tests/xc_reconciliation/results/xc_lda_exchange_validation.csv)
* [summary.json](../tests/xc_reconciliation/results/summary.json)
* [xcr02_before_repair_summary.json](../tests/xc_reconciliation/results/xcr02_before_repair_summary.json)

The reproducible command is:

```text
python3 tests/xc_reconciliation/run_xc_lda_reconciliation.py \
  --binary build-xcr/bin/UnitXcLdaReconciliation \
  --output-dir tests/xc_reconciliation/results
```

Focused checks added or rerun:

* `UnitXcLdaLegacyKernel` — corrected BH derivative and polarized limits for
  TXC=1, 3, 4, 5, and 11;
* `UnitXcLdaReconciliation` — full fixed-density grid, exchange gate, Stoner
  response, energy-gradient audit, and endpoint checks; and
* `UnitLibxcXcBaseline` — existing legacy/libXC production baseline.

The bounded fcc-Fe follow-up from XCR-01 remains separate. Its four tested
SCF cases did not converge within 24 iterations, so their last iterates are
not used as an XC formula oracle. Different genuinely parameterized LDA
functionals may still produce different nonlinear solid-state magnetic
solutions.

## Completion checklist

- [x] XCR-01 evidence reviewed before editing production code
- [x] TXC=1 versus TXC=101 equivalence status established
- [x] exchange contribution independently validated
- [x] spin interpolation analytically audited
- [x] small-polarization magnetic kernel audited
- [x] only demonstrated legacy defects changed
- [x] polarized-density limit repaired
- [x] full XCR-01 oracle rerun after the repair
- [x] focused LDA regression tests pass
- [x] fcc-Fe result kept separate from fixed-density proof
- [x] documentation completed

The Barth–Hedin suspicion is closed as a parametrization/implementation
question: TXC=1 is now reconciled with the TXC=101 reference at fixed density,
and the demonstrated legacy derivative and polarized-guard defects are fixed.
