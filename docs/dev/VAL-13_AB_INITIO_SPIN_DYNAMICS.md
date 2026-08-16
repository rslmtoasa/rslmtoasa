# VAL-13 — end-to-end ab-initio spin dynamics

## Scope and status

This validation establishes the deterministic, zero-temperature, zero-damping
one-site bcc-Fe scalar-relativistic collinear limit of the production loop:

```text
LMTO electronic state -> LMTO magnetic field/torque
  -> abspinlib Depondt predictor
  -> predicted electronic refresh
  -> abspinlib corrector
  -> corrected electronic refresh
```

The case starts from the established bcc-Fe `+z` state.  VAL-12 supplies the
field/torque seam and VAL-11 supplies the independent Depondt solver evidence.
The expected limit is a longitudinal field and zero transverse torque, so the
unit spin direction remains `+z`.  No stochastic run is used as an oracle.

The scoped result is **Validated** for this deterministic bulk equilibrium
loop.  General ab-initio spin dynamics remains **Experimental**: this does
not establish finite-temperature noise, damping relaxation, SOC precession,
multi-site exchange dynamics, STT/SHE, induced moments, or MPI equivalence.

## Production changes

`processing='sd'` now uses the existing abspinlib Depondt predictor/corrector
entry points.  The application wrapper supplies the solver's required
rank-3 arrays, initializes its Depondt/RNG workspaces, and converts the input
`dt` from the documented ps units to the solver's seconds.  Each step refreshes
the LMTO state at the predictor direction before applying the corrector and
refreshes it once more after the corrected direction.

The magnetic-moment output defect was traced to the site-output naming layer:
per-site diagnostics assumed that `element%symbol` was populated for every
active impurity/post-processing site.  The output layer now provides a generic
deterministic `atom<N>` stem when that metadata is blank and checks the open
operation.  Magnetic-moment output is retained; it is never skipped.
The production `Example_impurity_B2FeCo_sd_smoke` requires both the trajectory
and `Fe_1_spinene.out` output files.  A current-head serial reproduction of the
historical crash with the ordinary B2FeCo deck did not fail after the earlier
SD stack repair; the blank-site metadata condition is now covered by the
generic output guard, and the normal impurity production smoke passes.

## Numerical evidence

The labelled campaign is:

```sh
ctest --test-dir build-rf-debug -R '^Val13AbInitioSpinDynamics$' --output-on-failure
```

It runs fixed short interval trajectories at `dt = 0.01, 0.005, 0.0025 ps`
and a zero-step electronic baseline.  Current GNU Fortran Debug evidence:

| Check | Result |
|---|---:|
| Final refreshed longitudinal field | `B_z = -5.9354513e4 T` |
| Final refreshed transverse field norm | `4.62e-4 T` |
| Maximum trajectory norm error | `0` at output precision |
| Largest transverse component over timestep runs | `5.15e-9` |
| Largest pairwise final-vector timestep spread | `7.12e-9` |
| Zero-damping electronic energy drift vs zero-step baseline | `5.94e-5 Ry` |

The small transverse field is within the electronic-refresh envelope relative
to the longitudinal field; the production direction remains `+z`.  The field
reported above is taken after the corrected electronic refresh, rather than
from the pre-refresh setup diagnostic.
The energy number is a bounded feedback diagnostic, not an independent energy
functional oracle.

## Checklist

- [x] Deterministic bulk trajectory selected
- [x] Full LMTO/field/abspinlib/electronic-update loop exercised
- [x] Spin norm checked
- [x] Trajectory direction/sign checked (`+z` zero-torque limit)
- [x] Timestep convergence checked over a short interval
- [x] Energy behavior checked as a bounded zero-damping feedback envelope
- [x] Impurity empty-site-output condition traced to the output naming layer
- [x] Impurity output root cause fixed generically
- [x] Production impurity SD smoke added and passing
- [x] No stochastic reference used as primary oracle
- [x] Spin-dynamics maturity updated narrowly

## Changed files and tests

- `source/spin_dynamics.f90`
- `source/bands.f90`
- `tests/validation/val13_abinitio_spin_dynamics.py`
- `tests/scf/cases.json`
- `CMakeLists.txt`
- this report and the Phase-II/developer issue records

Passed:

```sh
ctest --test-dir build-rf-debug -R '^(UnitAbspinlibIntegrators|Example_bulk_bccFe_sd_smoke|Example_impurity_B2FeCo_sd_smoke|Val13AbInitioSpinDynamics)$' --output-on-failure
```

Commit message: `Validate end-to-end ab-initio spin dynamics`
