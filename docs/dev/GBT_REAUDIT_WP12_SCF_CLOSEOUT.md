# GBT WP12 re-audit: SCF diagnosis and frame-contract repair

Date: 2026-08-27

Status: **q=0 SCF closure is established for the scoped one-sublattice,
collinear bcc-Fe k-space fixture.** The finite-`q` constrained-SCF versus
explicit-supercell gate remains **not qualified**
because the existing production supercell counterpart is unmatched and its
convergence evidence is not closed.

This closeout records a localized SCF/frame repair. It does not reopen the
validated primitive `S x G` operator or alter the radial LMTO routines.

## Prerequisite evidence and diagnosis

WP01-WP11 evidence was reviewed before changing the SCF path. WP01-WP08
provide the exact q=0, gauge, covariance, frame, constraint-field,
constraint-energy, and corrected-MFT slices. WP09-WP11 identify unresolved production
qualification issues in reciprocal k-grid convergence, small-q curvature, and
LKAG interaction-range/curve comparison. Those issues are measurement gates,
not evidence of a primitive GBT operator defect.

The localized SCF divergence was at the magnetic-state boundary in
`self%run_dos()`. Both the real-space and reciprocal density producers
previously reached the full Cartesian magnetic mixer. Under the declared
`constrained_spiral` policy, that allowed a measured transverse component to
participate in the state used as the next GBT frame marker. This conflated:

- the fixed rotating-frame reference axis;
- the physical output moment and its transverse constraint residual; and
- the longitudinal/radial channels that are allowed to update.

The Hamiltonian exchange construction and radial routines were not implicated
by the q=0 trace, so neither was modified.

## Repair

The implementation in `source/self.f90` now has an explicit policy boundary:

- `constrained_spiral` keeps the normalized reference axis in `potential%mom`
  and in the magnetic mixer marker. Charge and longitudinal radial channels
  continue through the existing density/mixing path; the measured transverse
  moment is retained as a constraint residual.
- The measured incoming physical moment is captured separately from the frame
  marker, so diagnostics do not mistake the two quantities for one vector.
- GBT constraint-marker sign selection uses the fixed reference direction.
- `relaxed_reference` retains the legacy full-Cartesian mixer for ordinary SCF,
  but `gbt_single_q` with that policy now fails explicitly until a moving-
  reference GBT Hamiltonian contract is implemented and audited.
- The new opt-in `control%gbt_scf_diagnostics` flag emits one row per site and
  SCF iteration in `gbt_scf_diagnostics.dat`. Each row contains the incoming
  frame/moment/radial values, exchange markers, target and field, rotating-
  frame density-matrix components, output moment, longitudinal/transverse
  decomposition, mixer states, next radial channels and marker, physical
  energy, and separated constraint diagnostics.

The trace is off by default and does not change the equations. The public
keyword and developer-map documentation describe the frame and policy
semantics.

## Controlled q=0 closure

The reproducible driver is
[`tests/validation/val21_gbt_wp12_scf.py`](../../tests/validation/val21_gbt_wp12_scf.py).
It runs the same two-iteration, `use_kspace=.true.` bcc-Fe fixture with
`q=[0,0,0]`, `theta=0`, and `density_policy='constrained_spiral'` for both
`periodic_nc` and `gbt_single_q`. It checks trace completeness, row ordering,
and every numeric trace field.

Observed result:

```text
WP12 q=0 trace rows: 2; iterations: [1, 2]
WP12 q=0 ordinary/GBT maximum trace residual: 4.100009e-11
WP12 q=0 SCF diagnostic/closure fixture: PASS
```

The existing independent WP01 production identity was also rerun:

```text
WP01 q=0 production residual (abs, rel): 1.065814e-14 5.464117e-15
WP01 q=0 GBT q=0 link identity residual: 0.000000e+00
```

This is a scoped collinear one-sublattice closure, not a claim that the
unrun global finite-angle or multi-sublattice q=0 campaigns are complete.

## Post-repair regression matrix

| area | check | result and qualification |
|---|---|---|
| q=0 exact operator | `UnitGbtWp01Q0`, `val18_gbt_wp01_q0.py` | PASS; exact/production identity retained |
| shifted-k gauge | `UnitGbtWp02GaugeShiftedK` | PASS |
| commensurate operator | `UnitGbtWp04CommensurateSupercell` | PASS for the matched operator/source fixture; production SCF comparison remains open |
| frame and constraints | `UnitGbtFrameContract`, `UnitGbtConstraintCovariance`, `UnitGbtWp07ConstraintEnergy`, `UnitSpinDensity` | PASS |
| harmonic k-grid | `GbtWp09HarmonicKgrid` | PASS as the existing fixture/source check; production convergence is not qualified |
| small-q | `GbtWp10SmallQ` | PASS as the existing symmetry/analysis check; curvature convergence is not qualified |
| LKAG `J(q)` | `GbtWp11LkagJq` | PASS as the existing analysis check; production range/curve closure is not qualified |
| WP12 SCF trace/closure | `UnitGbtWp12ScfSourceContract`, `val21_gbt_wp12_scf.py` | PASS for the scoped q=0 fixture; maximum trace residual `4.100009e-11` |
| finite-q constrained SCF vs explicit supercell | existing VAL-16/WP04 production counterpart | **NOT QUALIFIED**: unmatched and not converged; no finite-q SCF claim is made |
| global finite-angle and multi-sublattice q=0 SCF | WP12 campaign extensions | **NOT RUN / NOT QUALIFIED** |
| incommensurate production spirals | deferred | **NOT QUALIFIED** until the finite-q matched gate closes |

The focused post-repair CTest subset passed 12/12, including the source
contracts, frame/constraint tests, commensurate operator fixture, WP09-WP11
analysis checks, and GBT oracles. The full project build also completed
successfully.

## Completion checklist

| WP12 requirement | status |
|---|---|
| prerequisite WP01-WP11 evidence reviewed | complete |
| one SCF iteration instrumented end to end | complete; opt-in per-site trace |
| first frame/semantic divergence identified | complete; magnetic mixer/frame-marker boundary |
| `constrained_spiral` policy explicit | complete |
| q=0 GBT SCF matches ordinary SCF | complete for the scoped fixture |
| finite-q matched explicit-supercell comparison | not qualified; existing gate remains open |
| constraint fields/moments compared in the correct frame | complete in the trace and fixed-axis path; broad material campaign remains open |
| radial routines changed only if evidence required it | complete; no radial changes |
| WP07 physical energy semantics preserved | complete; physical and constraint terms remain separate |
| earlier exact GBT gates rerun | complete; focused subset and WP01 production identity pass |
| production support and exclusions documented | complete |

The requested commit message for this work is:

```text
fix: close GBT self-consistency frame contract
```
