# XCR-04 — Magnetic SCF closure after XC reconciliation

## Scope and conclusion

This audit is limited to the ordinary `nsp=2`, `q=0` collinear path.  No GBT
source or GBT diagnostic was changed.  The production path now has opt-in
feedback diagnostics, explicit fresh/restart provenance, a removable magnetic
seed, and a radial ATOMSC residual audit.

The bounded evidence does not support promoting a converged high-spin fcc-Fe
fixed point: the six-step multi-seed campaign remains outside the outer SCF
convergence threshold, as did the earlier 24-step unseeded fcc-Fe follow-up.
It does show a persistent high-spin basin near 2 μB for most seeds, a
seed-dependent low-spin trajectory for the strongest TXC=1/101 seed, and
near-identical legacy/libXC feedback for the reconciled pairs.  The original
low-moment anomaly is therefore classified as **UNRESOLVED**, with the leading
non-XC hypothesis being **SCF_BRANCH_SELECTION**; the available runs are not
converged enough to claim that classification as settled physics.

## Post-XCR fixed-density baseline

The fixed-density results are inherited from XCR-02/XCR-03 and are not
reopened by the solid-state runs below.

| pair | max energy-density difference | max scalar-potential difference | max magnetic-potential difference | status |
| --- | ---: | ---: | ---: | --- |
| TXC=1 / TXC=101 (BH/VBH) | 8.88e-9 Ry | 1.18e-8 Ry | 6.46e-9 Ry | reconciled equivalent |
| TXC=5 / TXC=105 (PW/PW) | 7.44e-7 Ry | 8.94e-7 Ry | 8.95e-6 Ry | expected variant difference |
| TXC=8 / TXC=108 (PBE/PBE) | <9e-16 in `exc` | 5.06e-6 Ry | 5.30e-6 Ry | radial GGA comparison passes |

The LDA evidence is in [XC_LDA_RECONCILIATION.md](XC_LDA_RECONCILIATION.md)
and the GGA evidence is in
[XC_GGA_RADIAL_RECONCILIATION.md](XC_GGA_RADIAL_RECONCILIATION.md).

## Ordinary q=0 feedback map

For each active site and outer SCF iteration, the ordinary loop is:

```text
stored rho0 / current QL, PL
        |
        v
ATOMSC: POISS0 -> VXC0SP -> NEWRHO -> radial residual
        |                 |
        |                 +--> radial rho_up-rho_down and Vxc_up-Vxc_down
        v
RACSI -> POTPAR -> QL and orthogonal C_l spin channels
        |
        v
mix QL/PL -> Madelung -> run_scf -> Hamiltonian -> DOS/moment
```

`magnetic_scf_diagnostics.dat` records, per outer iteration/site/l, the
spin-channel `QL`, `QL_up-QL_down`, radial signed and absolute spin density,
radial signed and absolute `Vxc_up-Vxc_down`, orthogonal `C_l` and its
splitting, and the Cartesian/site moment.  The trace uses the potential
channel convention (`1=up`, `2=down`); the historical `VXC0SP` local variables
are down/up, and the source converts them back to the up/down output channels.

## ATOMSC residual derivation

`VXC0SP` converts the stored radial quantity as

```fortran
tRHO = RHO / (4*pi*r**2)
```

and the electron normalization is accumulated as `WGT*DRDI*RHO`.  Therefore
stored `RHO` is the spherical radial electron density (`4*pi*r**2*n`), and
the discrete L1 residual must be

```fortran
DRHO = DRHO + WGT*DRDI*ABS(RHO_new-RHO_old)
```

The old unweighted expression is retained as `residual_unweighted` for audit
only.  The stopping threshold remains `1e-6`, now applied to the integrated
electron-number norm rather than a mesh-coordinate-dependent sum.  No mixer
parameter was changed.

Representative fcc-Fe TXC=1 audit values from the first bounded run were:

| point | old unweighted norm | corrected integrated norm |
| --- | ---: | ---: |
| outer 1, inner 1 | 5.0772e3 | 2.4700e1 |
| outer 6, final recorded inner step | 4.9529e-5 | 6.5502e-7 |

The complete history is retained in the per-case `*.residuals.dat` traces in
the regression artifact directory.

## Fresh start, restart, and initialization

The legacy default remains restart-compatible: when `fresh_start=.false.` and
`./Fe_out.nml` exists, it is selected before `./Fe.nml`.  Setting
`control%fresh_start=.true.` suppresses the `_out.nml` choice and selects the
named fresh file (or the global database fallback).  The selected exact path
is logged as either:

```text
SCF input provenance: fresh atomic starting potential from ./Fe.nml
SCF input provenance: restart loaded from ./Fe_out.nml
```

The bounded campaign uses fresh-start mode for all physics rows and includes a
separate bcc-Fe restart provenance check.  It does not silently consume a
prior low- or high-spin output.

The initial radial density is exactly unpolarized: `symbolic_atom%rho0`
constructs one normalized radial channel and assigns `rho0(:,2)=rho0(:,1)`.
The initial magnetic asymmetry instead enters through the loaded potential
parameters (`QL`, `PL`, `C`, and the normalized moment direction) in `Fe.nml`
or `Fe_out.nml`.  The ordinary collinear `B_fsm` contribution is zero unless
the existing FSM/constraint route supplies it.

## Temporary magnetic seed

The new self controls are:

```fortran
magnetic_scf_diagnostics = .true.
magnetic_seed_enable = .true.
magnetic_seed_steps = 4
magnetic_seed_field = 0.005       ! Ry, in the existing B_fsm convention
```

The field is added only inside `ATOMSC` for the selected first outer steps.
It is automatically released before the next iteration and logs both the
activation and release.  It does not change `mom`, `QL`, the mixer, or a
constraint target.  It is explicitly rejected together with persistent
magnetic constraints, so a seeded trajectory cannot be reported as a
fixed-spin calculation.

## fcc-Fe seed sweep

The reproducible campaign is
[`run_magnetic_scf_regression.py`](../tests/magnetic_scf/run_magnetic_scf_regression.py),
with `a=3.683 Å`, `r_WS=1.43930285 Å`, common structure/k-point setup, and
nominal basin labels `0, 0.1, 0.5, 1, 2, 3`.  These labels are deliberately
recorded as *nominal seed labels*, not imposed moments: nonzero labels select
the temporary fields `B_fsm=0.005*label Ry`, and the measured release and
post-release moments are stored separately.

| XC/backend | nominal labels | final-M range after 6 outer steps (μB/Fe) | d-channel ΔC range (Ry) | result |
| --- | --- | ---: | ---: | --- |
| TXC=1 legacy | 0…3 | 0.3965–1.9814 | -0.1466…-0.0176 | bounded, not converged |
| TXC=101 libXC | 0…3 | 0.3912–1.9814 | -0.1466…-0.0178 | bounded, not converged |
| TXC=5 legacy | 0…3 | 2.0354–2.0501 | -0.1520…-0.1490 | bounded, not converged |
| TXC=105 libXC | 0…3 | 2.0354–2.0501 | -0.1520…-0.1490 | bounded, not converged |

The pairwise agreement is strong at this resolution, including the d-channel
splitting.  The high-spin basin is visible, but the strong TXC=1/101 seed
trajectory enters a low-moment bounded trajectory.  Since neither trajectory
passes the outer convergence gate, this is evidence for branch sensitivity,
not proof of a second self-consistent ground state.

The earlier 24-step unseeded follow-up reported final moments of 1.9515 and
1.9526 μB for TXC=1 and TXC=101, and 2.0247 μB for TXC=5 and TXC=105, while
remaining unconverged.  It is retained as traceability evidence in
[`fcc_fe_scf_followup.json`](../tests/xc_reconciliation/results/fcc_fe_scf_followup.json).

## bcc-Fe and nonmagnetic controls

The same bounded campaign gives bcc-Fe final moments of 2.1041 μB unseeded and
2.0655 μB after the nominal 2 μB temporary seed, with no runtime or NaN
failure.  This is a control/smoke pass, not a replacement for the existing
bcc-Fe reference workflow.

The spin-symmetric fcc-Cu control finishes at 9.2e-4 μB after the temporary
seed, demonstrating decay rather than artificial ferromagnetic locking in the
bounded run.

## Regression artifact and reproduction

The compact table is
[`magnetic_scf_regression.csv`](../tests/magnetic_scf/results/magnetic_scf_regression.csv)
and its configuration metadata is
[`metadata.json`](../tests/magnetic_scf/results/metadata.json).  The campaign
was run with:

```bash
python3 tests/magnetic_scf/run_magnetic_scf_regression.py \
  --binary build-xcr/bin/rslmto.x \
  --output-dir tests/magnetic_scf/results \
  --max-steps 6 --keep-traces
```

The six-step bound is intentional: the 24-step matrix is too expensive for
the full seed/backend sweep in the available CPU budget.  The artifact marks
every row `BOUND_NOT_CONVERGED`; no final moment is presented as a converged
fixed point.

## Completion status

- [x] post-XCR fixed-density XC baseline recorded
- [x] ordinary q=0 magnetic feedback loop instrumented
- [x] ATOMSC convergence metric derived
- [x] radial Jacobian incorporated where required
- [x] old/corrected residual histories captured
- [x] fresh-start/restart provenance made explicit
- [x] initial radial asymmetry source identified
- [x] temporary magnetic seed is removable and logged
- [x] fcc-Fe six-label bounded sweep completed
- [ ] converged high-spin q=0 fcc-Fe fixed point established
- [x] TXC=1 versus TXC=101 bounded SCF comparison performed
- [x] TXC=5 versus TXC=105 bounded SCF comparison performed
- [x] bcc-Fe bounded control passed without runtime regression
- [x] nonmagnetic fcc-Cu seed/release control passed
- [ ] MnSi/FeGe optional follow-up performed
- [x] GBT left untouched
- [x] original anomaly classified conservatively as UNRESOLVED

Renewed GBT validation should wait for an independently converged ordinary
q=0 fcc-Fe result; the present substrate is instrumented, but not yet closed.
