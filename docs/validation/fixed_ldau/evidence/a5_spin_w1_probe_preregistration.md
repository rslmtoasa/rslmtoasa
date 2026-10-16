# A5.3 spin-selective W1 response preregistration

Measure the local-mode U=J=0 base-functional gradient needed by the approved
A5.3 charge-only chain rule. Start every case from the same immutable best
available common U=J=0 KS restart snapshot. Its mixer residual was
`5.9e-9`, so the starting state is near-converged, not convergence-accepted;
it is only an initial guess. Each final probe SCF must independently meet its
registered `conv_thr`. If the first `lambda=0` standard probe reaches its
240-step cap without literal convergence, preserve its full input, residual
history, trace, and restart hashes, then continue from that case's latest
`Fe_scf.nml` copied as `Fe_out.nml`. Change only Broyden `beta` from 0.1 to
0.05, retain all physical inputs and `conv_thr`, and allow 240 additional
steps. This is a continuation of the same initial-state chain, not a new
independent seed. The endpoint remains usable only with literal convergence;
otherwise report the probe allowance unresolved. Apply the same d-shell W1 shift separately to spin up and spin down, with
`lambda = ±0.005` and `±0.0025 Ry`; capture one `lambda=0` endpoint. Use the
direct P_OUT matrix to record `q_s=Tr(P_s)` and the applied `Tr(WP_s)` for both
spins. Use the logged `F_local` at each signed shift. The historical A5.1
local W1 residual was about `−9.2e−3` occupation; the current spin split and
response matrix are measured from direct production P. Candidate c's
separate U-derivative prediction remains a miss of about `dt/dU≈0.13 Ry/eV`.

For each probe spin `a`, form
`B_a(h)=[F_a(+h)-F_a(-h)]/(2h)-q_a(0)` and
`R_s,a(h)=[q_s(+h_a)-q_s(-h_a)]/(2h)`. Use `h=0.0025 Ry` as the primary
estimate and report `h=0.005 Ry` beside it; do not require an h² scaling
classification. Solve `R^T g=B`, where `g_s` is the residual energy gradient
per d occupation in Ry/electron. Combine this with the independently measured
`dq_s/dU` from the A5.3 U scan to obtain
`N_ASA=sum_s |g_s dq_s/dU|` in Ry/eV. This uses the user-approved
charge-only assumption; it does not represent unmeasured orbital or coherence
gradients.

Use the measured h-step differences in B and R, the λ=±0.005 tightened
repeat, λ=0 endpoint repeat, and ordinary output rounding as the probe
uncertainty. Propagate them through the 2×2 solve and the per-spin U response
to report uncertainty on `N_ASA`. Report the two-spin reconstructed B sum
beside the historical A5.1 local aggregate residual and its uncertainty as a
diagnostic only. A5.1 used an eigen-reconstructed P that differs from this
direct production-P definition; the comparison is not an agreement criterion.
If the measured response cannot determine g within its propagated
uncertainty, report the resulting allowance as unresolved; do not add a
condition-number cutoff. Register the numeric allowance before examining
candidate a/b/c derivatives. A b pass with a c miss by the predicted `dt/dU`
supports the user-approved A5.3 decision; a b miss or c pass stops the
adoption decision.

For the approved J-center checks, use the same measured local charge gradient
`g` under the disclosed charge-only linearization. At U=2 eV, obtain each
spin occupation from the direct P_OUT trace at J centers 0.5 and 0.9 eV.
For each center `J*`, form the primary physical-domain slope
`D_s(J*;0.05)=[q_s(J*+0.05)-q_s(J*-0.05)]/0.10` and the step-sensitivity
check `D_s(J*;0.10)=[q_s(J*+0.10)-q_s(J*-0.10)]/0.20`, in electron/eV.
Report the two centers separately, with
`N_J(J*)=sum_s |g_s D_s(J*)|` in Ry/eV. Propagate the measured `g`
uncertainty and the observed h=0.05 versus h=0.10 response difference; state
any absence of tightened J-state repeats explicitly. These local allowances
apply at J=0.5 and J=0.9 only and do not extrapolate to J=0.

The probe uses the common 12³, 300 K KS state, fixed Liechtenstein U=J=0,
local radial moments, and W1. Record hashes for each source file, executable,
input, seed, and output; do not create tooling locks or post-run criteria.

Pre-run source hashes: `hamiltonian.f90`
`db26fbc766f94c811ae641c0fe2dd918bb913be3fce92f4488b9af826ce29b2f`,
`hamiltonian_build.f90`
`e896d53d82560a0551cc44560191dda5a61a70418a6b75e5c6247c758cd7ca90`,
`hamiltonian_hubbard.f90`
`8dcab0b0fa467b0e8532877e3e7e00a96cef7e47dce38a4ac085dd26dada43ee`,
`include_codes/namelists/hamiltonian.f90`
`afde6c79db8fe06d8d5bd61fbb90720ee4990c9337b6b4efc530da1297edfe70`, and
`self.f90`
`1c29e173c7ae0061f2a9087cc90fbe44db0a24000dfda183340a914d8e20a7ad`.
The final moments-only Release executable hash is
`52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`; the
common 12³/300 K input hash is
`ca65d9695b049f10beb3ee9461c05572b48f3973885de99ea7f8fee3dfc5699c`; the
immutable shared near-converged `Fe_out.nml` restart hash is
`25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5`.
