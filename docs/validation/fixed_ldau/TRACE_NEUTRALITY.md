# Fixed-U trace observer check

## Result and limits

Two fresh one-step RS runs using the same binary, the same saved input deck,
and the same saved `Fe.nml` seed produced byte-identical scientific outputs
with the trace environment enabled and disabled. A second matched pair gave
the same result. In these controlled replays the opt-in trace did not change
the recorded `SUMEV`, `ETOT`, radial density, SCF parameters, or total DOS.
This is bounded evidence for the tested one-step RS case. It does not prove
trace neutrality for every backend, build, or SCF trajectory.

Read-only source review of `source/self.f90` found that the helpers read
`potential%ldm`, radial arrays/energy outputs, `green%g0`, and reciprocal
Hamiltonian/eigenpair caches, then write only the separate CSV file. They do
not call `RSEQSR`, `POTPAR`, `NEWRHO`, `atomsc`, or another physics routine.
The trace path sets private iteration/site selectors; those are only consumed
by the trace guards. The k-space hook is placed after H(k) assembly and before
diagonalization and only reads the H(k) cache. This source review plus
matched-output parity rules out an obvious extra-physics-call or in-place
matrix-mutation mechanism in the inspected path; it does not replace broader
backend/multistep replay.

An earlier trace-on/off comparison reported a `SUMEV` difference of about
`4.95e-8 Ry`. That one-off difference was not reproduced in the matched runs
below and remains unexplained. Do not attribute it to a trace observer effect
or to different starting Fermi energies: `input.nml` is rewritten during a
run with the updated Fermi level, so the post-run file is not evidence of the
pre-run deck. The controlled source deck was preserved separately before
execution.

## Reproduction

Binary SHA256:

```text
641ab24c021780dac1c5a85275a62da140706ac33b94591ab1d7e3336c03a593
```

The immutable pre-run input snapshots are retained at
`/tmp/rslmto-energy-trace/observer_seed_inputs/`:

```text
input.nml  SHA256 0e66d7f32ca5a577aa240097902ff64ac350a138d6d6a6ca939d296e3a16ae5c
Fe.nml     SHA256 fd5d907628cc4ee41ac3b1147aeb1f9e8ec1fe8bb2f38069afe1cf63384ae233
```

The runs used `OMP_NUM_THREADS=1`, `OMP_THREADS_SERIAL=1`, one SCF step, and
the same executable. For the first pair, fresh working directories received
byte copies of those two snapshots. The only difference was whether the four
`RSLMTO_FIXED_U_*` variables were set. Run each command from the indicated
working directory. The trace-on command was:

```sh
env OMP_NUM_THREADS=1 OMP_THREADS_SERIAL=1 \
  RSLMTO_FIXED_U_TRACE=/tmp/rslmto-energy-trace/observer_seed097052_on/records.csv \
  RSLMTO_FIXED_U_RUN_ID=observer_seed097052_on \
  RSLMTO_FIXED_U_LAMBDA=1 \
  RSLMTO_FIXED_U_TRACE_ITERATION=1 \
  /tmp/rslmto-energy-trace/rslmto-trace.x > testrun.log 2>&1
```

Working directory: `/tmp/rslmto-energy-trace/observer_seed097052_on`.

The matching trace-off command was:

```sh
env -u RSLMTO_FIXED_U_TRACE -u RSLMTO_FIXED_U_RUN_ID \
  -u RSLMTO_FIXED_U_LAMBDA -u RSLMTO_FIXED_U_TRACE_ITERATION \
  OMP_NUM_THREADS=1 OMP_THREADS_SERIAL=1 \
  /tmp/rslmto-energy-trace/rslmto-trace.x > testrun.log 2>&1
```

Working directory: `/tmp/rslmto-energy-trace/observer_seed097052_off`.

Both commands exited 0. Across the matched pairs, `Fe_out.nml` was byte
identical (SHA256
`99591cea36b018ba35f6737d60384a6065b9265cfb9bd98bcb3826d22917d7e1`):

```text
SUMEV = -1.9956221522467070 Ry
ETOT  = -2541.9803776218619 Ry
```

For the first pair, the trace-on and trace-off
`radial_ground_state_Fe_1.dat` hashes both equal
`1a7cc0abebddeab3db049b7ecbeb44216184e52b0b7952ba8ae6863019553adf`;
`Fe_scf.nml` and `totaldos.out` hashes also matched. The trace CSV's total
`SUMEV` row exactly matched the trace-off `Fe_out.nml` value.

Scratch runs and logs are retained under `/tmp/rslmto-energy-trace/` in the
`observer_seed097052_*` and `observer_sameinput_*` directories. These are
runtime artifacts, not source-controlled fixtures.
