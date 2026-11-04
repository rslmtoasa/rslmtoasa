# Testing policy for U + SOC runs

Convergence of a surface or an impurity is slow and sensitive, and a host converged without SOC moves when SOC is switched on. Convergence is therefore not the measure of correctness here. Correctness is established by operator and kernel oracles, by comparison at iteration 1 and over a few steps against an equivalent run, and by watching how the iteration moves.

## Tier 1 (default): few iterations, trajectories recorded

- 1 to 5 iterations from a seed, reduced settings where the check allows (regression rows: 4 steps, lld 50, 500 channels).
- Record per step, from the log: Diff; d trace (or n_d) per layer; spin moment; orbital moment; max |P_up,dn| of the spinor d density matrix. Report the numbers, not "ran" or "converged".
- Equivalence checks (SOC->0 against nsp=1, tilted axis against z): iteration 1 within the stated tolerance, then the 4-step regression rows. The nsp=1 comparison run takes the same host numbers as the nsp=2 run (derive the nsp=1 host from the nsp=2 host: spinor_d line dropped, hubbard_nsp=1, new dm hash); an independently converged host differs at its own convergence level (about 1e-8) and spoils a 1e-9 check.
- Tilted axes: seed the run with the exactly rotated state (`tests/validation/fixed_ldau/rotate_spin_frame.py`: host and surface P rotated, mom and hubbard_axis along the axis, dm hashes updated), not with a collinear z-frame seed along a tilted axis. The collinear seed gave a slowly relaxing x host (KNOWN_ISSUES, S6c).
- Frozen hosts need not be converged for a smoke or trajectory run; state the host's Diff and step count in the deck README.
- No convergence studies, no reference-state acceptance windows (last-20 spreads), no energy-anisotropy runs, unless the task asks for them by name.

## Tier 2 (optional): one longer run

At most one run per task, at most 20 iterations, from the step-5 state of the Tier 1 run, `ldm_beta=0.3`. Record the Diff trend and the orbital-moment trajectory. Flag it only if Diff trends upward or goes non-finite. The result is a trajectory, not a converged state, and is reported as such.

## Reporting

Table of check, measured, tolerance, command. Trajectories go in the README of the deck directory. What was not run is listed as not run.
