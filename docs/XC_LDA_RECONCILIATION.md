# XCR-01 — Fixed-density legacy/libXC LDA reconciliation

Date: 2026-08-30

Binary: `UnitXcLdaReconciliation`, GNU Fortran 13.3.0, libXC 5.2.3, serial CPU

## Conclusion first

The legacy Barth–Hedin route (`TXC=1`) is reconciled with libXC
`LDA_X + LDA_C_VBH` (`TXC=101`) for the quantities that drive the local
magnetic response:

* the largest regular-grid energy difference is `8.83e-9 Ry` per electron;
* the largest regular-grid `B_xc` difference is `4.87e-9 Ry`;
* the largest relative difference in the symmetric small-polarization
  `I_xc` response is `7.53e-9`.

It is not, however, a full mathematical reconciliation of the two reported
spin potentials. At large polarization the legacy `V_xc^up` and `V_xc^down`
have a common charge-channel offset, reaching `1.55e-2 Ry` at
`(r_s,zeta)=(0.5,0.99)`. An independent finite-difference audit of the legacy
energy shows the same offset: the reported legacy potential is not the exact
spin derivative of the energy expression at those polarized points. The
magnetic difference cancels this common offset, which explains why the
`B_xc` and Stoner results agree while `V_0` does not.

No production XC expression was changed. The fixed-density evidence gate is
complete; a possible legacy repair must be derived and reviewed separately.

## Scope and driver

The driver is [test_xc_lda_reconciliation.f90](../tests/unit/test_xc_lda_reconciliation.f90).
It calls the production `XCPOT` method for the legacy route and the production
`xcpot_libxc_wrapper` for the libXC route. It does not call band structure,
Fermi-level determination, charge mixing, `NEWRHO`, `ATOMSC`, potential
parameter construction, or Hamiltonian construction.

The sweep contains nine Wigner–Seitz radii
`0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6 bohr`, and the requested polarizations
`0, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.25, 0.5, 0.8, 0.95, 0.99`.
The exactly fully polarized point `zeta=1` is emitted separately for every
radius and pair. It is not included in regular-range maxima because the
legacy routine intentionally returns zero when either input channel is below
`1e-20`; the libXC wrapper instead clamps the channel to `1e-20`.

The reproducible runner is
[run_xc_lda_reconciliation.py](../tests/xc_reconciliation/run_xc_lda_reconciliation.py):

```text
python3 tests/xc_reconciliation/run_xc_lda_reconciliation.py \
  --binary build-xcr/bin/UnitXcLdaReconciliation \
  --output-dir tests/xc_reconciliation/results
```

## Conventions established before comparison

The density parameterization is

\[
n = n_\uparrow+n_\downarrow = \frac{3}{4\pi r_s^3},
\qquad
\zeta = \frac{n_\uparrow-n_\downarrow}{n},
\]

with

\[
n_\uparrow=\frac{n(1+\zeta)}{2},\qquad
n_\downarrow=\frac{n(1-\zeta)}{2}.
\]

The historical argument order is established by the call and the wrapper
mapping in [source/xc.f90](../source/xc.f90):

| Quantity | Legacy `XCPOT` argument/output | libXC array/production mapping |
|---|---|---|
| spin-up density | `RHO2` | `rho_libxc(1)` |
| spin-down density | `RHO1` | `rho_libxc(2)` |
| spin-up potential | `V2` | `vrho_libxc(1)` |
| spin-down potential | `V1` | `vrho_libxc(2)` |

The driver therefore passes `RHO1=rho_down` and `RHO2=rho_up` to both
production paths. The fixed-density calls use zero gradients and a positive
dummy radius of `1 bohr`; the radius is irrelevant for LDA.

| Item | Convention |
|---|---|
| density | electrons/bohr³ |
| `r_s` and radius | bohr |
| `EXC` | XC energy per electron, not energy density |
| potential | Rydberg per electron |
| libXC native output | Hartree per electron and Hartree potential |
| libXC conversion | production wrapper multiplies `exc`, `vrho`, and `vsigma` by 2 to obtain Rydberg |
| scalar component | `V0=(V_up+V_down)/2` |
| magnetic component | `B_xc=(V_down-V_up)/2` |
| splitting | `Delta V_xc=V_up-V_down=-2 B_xc` |

The factor-of-two conversion is visible in the wrapper immediately before its
LDA/GGA output dispatch. The exchange-only calibration below independently
checks the density, spin ordering, sign, and Hartree-to-Rydberg conversion.

## Analytic exchange validation

For spin-LDA exchange in the driver's Rydberg convention, the analytic
reference is

\[
\epsilon_x = -\frac{3}{2}\left(\frac{6}{\pi}\right)^{1/3}
\frac{n_\uparrow^{4/3}+n_\downarrow^{4/3}}{n},
\]

\[
v_{x,\sigma}=-2\left(\frac{6}{\pi}\right)^{1/3}n_\sigma^{1/3}.
\]

The driver evaluates native libXC `XC_LDA_X` through the production wrapper
(`TXC=1001`) against these expressions over the non-endpoint grid. The maximum
errors were:

| quantity | maximum absolute error | maximum relative error |
|---|---:|---:|
| exchange energy/potentials together | `1.55e-15 Ry` | `2.95e-15` |

This passes the analytic exchange gate and rules out a unit or spin-ordering
explanation for the later correlation/potential differences.

## Pairwise fixed-density results

The complete machine-readable table is
[xc_lda_reconciliation.csv](../tests/xc_reconciliation/results/xc_lda_reconciliation.csv).
It contains `324` rows: `297` regular points and `27` separately marked fully
polarized points. The requested fields are present, along with per-quantity
absolute/relative differences, `V0`, `B_xc`, `Delta V_xc`, evaluation status,
and classification.

Regular-grid maxima:

| pair | `max |Delta EXC|` (Ry) | `max |Delta V_up|` (Ry) | `max |Delta V_down|` (Ry) | `max |Delta B_xc|` (Ry) | max reported relative difference | classification |
|---|---:|---:|---:|---:|---:|---|
| `TXC=1` / `TXC=101` BH/VBH | `8.83e-9` | `1.55e-2` | `1.55e-2` | `4.87e-9` | `1.56e-2` | `EXPECTED_PARAMETRIZATION_DIFFERENCE` at numerical-level rows; `LEGACY_KERNEL_DEFECT` once the common potential residual exceeds `1e-7 Ry` |
| `TXC=5` / `TXC=105` PW92/PW | `6.62e-7` | `8.73e-7` | `1.48e-5` | `7.85e-6` | `1.49e-5` | `EXPECTED_PARAMETRIZATION_DIFFERENCE` |
| `TXC=4` / `TXC=106` VWN/VWN | `5.08e-4` | `8.36e-4` | `1.07e-2` | `5.35e-3` | `9.43e-2` | `EXPECTED_PARAMETRIZATION_DIFFERENCE` |

Pair B is a useful control, but the legacy `PBEGGA` LDA branch and libXC's
native `LDA_C_PW` are corresponding implementations, not an established
bitwise identity. Pair C shows that the legacy VWN expression and the native
ID-7 VWN implementation cannot be treated as the same version without a
separate parametrization audit.

## Legacy energy-derivative localization

The additional
[xc_lda_energy_gradient.csv](../tests/xc_reconciliation/results/xc_lda_energy_gradient.csv)
file compares the reported legacy potentials to symmetric finite differences
of the reported legacy energy density `n*EXC`. Its maximum residuals are:

| pair | maximum legacy potential-vs-energy-derivative residual |
|---|---:|
| BH/VBH | `1.55e-2 Ry` |
| PW92/PW | `7.42e-10 Ry` |
| VWN/VWN | `4.47e-9 Ry` |

For legacy BH, the first differing expression is the potential reconstruction
in the default Barth–Hedin branch of `XCPOT`, not the reported energy
interpolation. In [source/xc.f90](../source/xc.f90), that is the `ARS`/`BRS`
construction and the two `TRX` assignments following the `EXC` expression.
Writing `x=RHO1/RHO`, `F(x)` for the legacy spin interpolation, and
`K(n)=(CNY+(4/3)EPSXP)/5.1297628`, the reported energy has the form

\[
\epsilon(n,x)=\epsilon_0(n)+F(x)K(n).
\]

The exact derivatives of the energy density `n epsilon` are

\[
v_\downarrow = \epsilon_0+n\epsilon_0'
 +F(K+nK') +(1-x)F'K,
\]

\[
v_\uparrow = \epsilon_0+n\epsilon_0'
 +F(K+nK') -xF'K.
\]

The measured legacy residual is almost the same in both spin channels, so it
is a common charge-channel term. It cancels from `B_xc` and therefore does not
appear in the small-polarization magnetic derivative. This is a localization,
not a repair: the legacy formula remains unchanged in this diagnostic.

## Small-polarization Stoner response

The response CSV is
[xc_lda_stoner_response.csv](../tests/xc_reconciliation/results/xc_lda_stoner_response.csv).
At fixed `n`, the driver uses symmetric `+delta` and `-delta` polarizations
and evaluates

\[
I_{xc}=\frac{[V_\downarrow-V_\uparrow](+\delta)
              -[V_\downarrow-V_\uparrow](-\delta)}
             {2n\delta}.
\]

This is the signed derivative requested in the task, in `Ry bohr³`. The
response was evaluated at `delta=1e-2, 1e-3, 1e-4, 1e-5, 1e-6`; the first
three values demonstrate convergence before roundoff begins to affect the
smallest step.

Representative `delta=1e-4` values:

| pair | `r_s` | legacy `I_xc` | libXC `I_xc` | relative difference |
|---|---:|---:|---:|---:|
| BH/VBH | 2 | `9.08014136` | `9.08014131` | `5.87e-9` |
| BH/VBH | 6 | `64.89597429` | `64.89597381` | `7.52e-9` |
| PW92/PW | 2 | `9.48136335` | `9.48135618` | `7.57e-7` |
| PW92/PW | 6 | `57.18527416` | `57.18522030` | `9.42e-7` |
| VWN/VWN | 2 | `9.45401827` | `9.61937341` | `1.72e-2` |
| VWN/VWN | 6 | `57.73400406` | `60.87529921` | `5.16e-2` |

Plots are diagnostic only:

* [bxc_vs_zeta.png](../tests/xc_reconciliation/results/bxc_vs_zeta.png)
* [stoner_response.png](../tests/xc_reconciliation/results/stoner_response.png)

## Discrepancy classification

The classifications in the CSV use the requested vocabulary:

| classification | use in this run |
|---|---|
| `EXPECTED_PARAMETRIZATION_DIFFERENCE` | Pair B and Pair C, and BH rows whose residual is at numerical/reference level |
| `LEGACY_KERNEL_DEFECT` | BH polarized rows with the reproducible common potential residual, and all exact fully polarized legacy-guard rows |
| `UNIT_OR_CONVENTION_DIFFERENCE` | not observed after the analytic exchange gate |
| `SPIN_ORDERING_DIFFERENCE` | not observed; the exchange test and zero-polarization symmetry pass |
| `ENERGY_NORMALIZATION_DIFFERENCE` | not observed; `EXC` is per particle and the factor two is reconciled |
| `LIBXC_WRAPPER_DEFECT` | not indicated by the exchange gate or the magnetic response |
| `UNCERTAIN` | none in the emitted classification column |

## SCF follow-up

A bounded optional fcc-Fe sanity check was performed for `TXC=1`, `TXC=101`,
`TXC=5`, and `TXC=105`, using `a=3.683 Å`, the corresponding Wigner–Seitz
radius `1.43930285 Å`, `conv_thr=1e-6`, and `24` maximum iterations. None of
the four cases converged; the final RMS differences were `9.58e-3` (BH),
`9.58e-3` (VBH), `5.63e-3` (legacy PW), and `5.64e-3` (libXC PW). Therefore
these last iterates are not used to attribute a magnetic collapse or to claim
an SCF XC spin splitting.

For traceability, the last-iterate total energies, Fermi energies, moments,
and internal d-channel `pl` potential parameters are recorded in
[fcc_fe_scf_followup.json](../tests/xc_reconciliation/results/fcc_fe_scf_followup.json).
The d-channel `pl` up-minus-down splittings were `0.1725833` (TXC=1),
`0.1726112` (TXC=101), `0.1697256` (TXC=5), and `0.1697628` (TXC=105).
Those are unconverged internal potential parameters, not XC-only spin
splittings. Representative fixed-density `Delta V_xc` values remain in the
main reconciliation CSV.

## Test and change checklist

- [x] legacy and libXC paths callable at identical fixed density
- [x] spin ordering established
- [x] units and energy normalization established
- [x] analytic exchange test passed
- [x] `TXC=1` versus `TXC=101` evaluated
- [x] `TXC=5` versus `TXC=105` evaluated
- [x] optional `TXC=4` versus `TXC=106` evaluated
- [x] small-polarization Stoner kernel extracted with step convergence
- [x] discrepancies classified
- [x] no legacy formula changed before diagnosis
- [x] machine-readable evidence and plots produced
- [x] documentation completed
- [x] fcc-Fe SCF follow-up attempted and recorded as non-converged; excluded from conclusions
