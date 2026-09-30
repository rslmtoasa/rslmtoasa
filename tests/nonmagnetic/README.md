# Per-atom forced nonmagnetic SCF

Set one logical value per self-consistent atom in `input.nml`:

```fortran
&self
  force_nonmagnetic = .false., .true.
/
```

All entries default to `.false.`. Entries follow SCF order `1:nrec`, mapping to
`symbolic_atoms(nbulk + ia)`, so impurity host types are not accidentally selected.
Partial namelist assignments retain false defaults for unspecified entries.

The constraint projects the **SCF charge/spin variables** onto equal spin channels.
It does not freeze the local charge: charge can redistribute during SCF, but each
projection preserves its input spin sum. It does not remove SOC. Hybridization
with magnetic neighbors can still polarize the *unprojected* Green function/DOS;
this is not a Lagrange-multiplier solver enforcing zero raw band spin density.
`mom0`, `mom1`, `mtot` and the reported spin moment are the constrained (zero)
values; `mom` remains a unit reference axis along z. Energy-resolved spin output
retains the raw band response. Unselected atoms retain their usual SCF equations
and may respond physically to the changed neighboring atom.

Application points:

1. Read the flags on every rank and project restart QL/PL and TB parameters before
   the first Hamiltonian construction.
2. Average all three QL moments and PL after band integration, before MPI exchange
   and storage in mixing history. Project old history inputs and mixer output;
   bypass magnetic-direction mixing for selected atoms.
3. Average the initial radial density, the density entering each XC evaluation,
   and each new/mixed density in atomic SCF. Set `B_fsm=0` and clear the atom's
   magnetic constraining field.
4. Average `center_band`, `width_band`, `shifted_band`, and `obar` after `predls`
   and before `build_pot`. The existing Hamiltonian formulas then give exactly
   zero `cx1`, `wx1`, `cex1`, and `obx1`; no Hamiltonian formula changes are needed.
5. Guard zero-norm report directions and angles. Zero is the output convention
   for undefined directions/angles. Report files are written only by rank zero.

`report.out` includes `Nonmagnetic check atom ...` with columns:
`Q_up-Q_down`, `max|cx1|`, `max|wx1|`, `max|cex1|`, `max|obx1|`, `|mom0|`.
The coefficients are evaluated by `build_pot` after the final atomic update.

## Run regression tests

Requires the usual Fortran/BLAS/LAPACK toolchain and Python `f90nml`.

```sh
cmake -S . -B build -DRUN_NONMAGNETIC_TESTS=ON
cmake --build build -j
ctest --test-dir build -R Nonmagnetic --output-on-failure

cmake -S . -B build_mpi -DENABLE_MPI=ON -DRUN_NONMAGNETIC_TESTS=ON
cmake --build build_mpi -j
ctest --test-dir build_mpi -R Nonmagnetic_SCF_MPI --output-on-failure
```

`projection.f90` checks spin sums of all QL/PL entries and a radial density,
parameter averages, the actual `build_pot` coefficients, field removal, and an
untouched magnetic atom. `run_scf.py` reuses the existing Fe/Co inputs to create a
small B2 bulk calculation. It runs four SCF steps with omitted flags, all-false
flags, and constrained Co. It checks exact equality of the default outputs,
all spin-channel entries, the runtime diagnostics, magnetic Fe, and finite
zero-moment output. The nsp=4 case uses a tilted magnetic axis and includes SOC.
The MPI case distributes the two atoms across two ranks. These are short
regressions, not converged material predictions.

For comparison against a separately built, unmodified executable:

```sh
python tests/nonmagnetic/run_scf.py --binary build/bin/rslmto.x \
  --baseline /path/to/unmodified/rslmto.x --workdir build/nonmagnetic/baseline --nsp 4
```
