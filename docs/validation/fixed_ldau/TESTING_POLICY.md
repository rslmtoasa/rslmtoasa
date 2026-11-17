# Testing policy for U + SOC runs

Convergence of a surface or an impurity is slow and sensitive, and a host converged without SOC moves when SOC is switched on. Convergence is therefore not the measure of correctness here. Correctness is established by operator and kernel oracles, by comparison at iteration 1 and over a few steps against an equivalent run, and by watching how the iteration moves.

## Tier 1 (default): few iterations, trajectories recorded

- 1 to 5 iterations from a seed, reduced settings where the check allows (stored-reference rows: 4 steps, lld 50, 500 channels).
- Record per step, from the log: Diff; d trace (or n_d) per layer; spin moment; orbital moment; max |P_up,dn| of the spinor d density matrix. Report the numbers, not "ran" or "converged".
- Equivalence checks (SOC->0 against nsp=1, tilted axis against z): iteration 1 within the stated tolerance, then the 4-step stored-reference rows. The nsp=1 comparison run takes the same host numbers as the nsp=2 run (derive the nsp=1 host from the nsp=2 host: spinor_d line dropped, hubbard_nsp=1, new dm hash); an independently converged host differs at its own convergence level (about 1e-8) and spoils a 1e-9 check.
- Tilted axes: seed the run with the exactly rotated state (`tests/ldau/rotate_spin_frame.py`: host and surface P rotated, mom and hubbard_axis along the axis, dm hashes updated), not with a collinear z-frame seed along a tilted axis. The collinear seed gave a slowly relaxing x host (KNOWN_ISSUES, S6c).
- No convergence studies of surfaces or impurities, no reference-state acceptance windows (last-20 spreads), no energy-anisotropy runs, and no runs of more than 20 iterations, unless the task asks for them by name. Exception to the iteration limit: states used as FROZEN HOSTS for embedded runs must be converged (they are inputs, not results); bulk hosts converge reliably. State the host's iteration count and final Diff in the deck README.

## Stored references

The reference binary of the earlier regression tool is replaced by scalars stored as JSON in `tests/ldau/references/` (written from the reference build `90f57c3` for the no-U and nsp=1 rows, from the committed SOC-on snapshot for the SOC-on rows). `tests/ldau/reference_rows.py ROW` runs one row with the candidate binary and compares by value within the row's tolerance: 1e-10 for rows against a stored state, 1e-9 at iteration 1 and 1e-9 to 1e-7 over 4 steps for SOC->0 against nsp=1. Five rows are registered with ctest (tier `standard`); the others are archived data and run only by name (`reference_rows.py --list`).

A reference is rewritten only with `reference_rows.py ROW --binary BIN --update "one-line reason"`; the reason, binary hash and date are stored in the JSON. Never rewrite a reference to make a check pass.

## Tier 2 (optional): one longer run

At most one run per task, at most 20 iterations, from the step-5 state of the Tier 1 run, `ldm_beta=0.3`. Record the Diff trend and the orbital-moment trajectory. Flag it only if Diff trends upward or goes non-finite. The result is a trajectory, not a converged state, and is reported as such.

## Reporting

Table of check, measured, tolerance, command. Trajectories go in the README of the deck directory. What was not run is listed as not run.
