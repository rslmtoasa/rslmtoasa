# A5.3 KS U/J energy-derivative experiment preregistration

At each spectral and local KS state `P(U=2 eV,J=0)`, compare the U
derivative of each energy candidate with `dE_Hub/dU` evaluated at that fixed
P by the production Liechtenstein `evaluate_hubbard_energy` routine at U0±
{0.01, 0.005, 0.0025} eV. For J, measure centered KS responses at J=`0.5`
and `0.9 eV` using J=`0.4,0.45,0.5,0.55,0.6` and
`0.8,0.85,0.9,0.95,1.0 eV`, respectively, in each radial mode. Evaluate
the fixed-P production J derivative at each center using the same small
{0.01, 0.005, 0.0025} eV steps. Candidate c is predicted to
miss by approximately `dt/dU` (the previous estimate is `0.13 Ry/eV`, to be
recomputed from the captured state); b is predicted to agree within the A5.3
base-functional allowance; a is report-only unless its separation from b
exceeds that allowance. Use the A4.2 U points `1.6, 1.8, 1.9, 1.95, 2.0,
2.05, 2.1, 2.2, 2.4 eV` in spectral and local modes, plus the J stencils
above at U=`2 eV`, each from the same frozen near-converged fresh U=0 KS
restart snapshot. This starting snapshot is common, not an asserted converged
endpoint; every final KS state is assessed against the unchanged SCF
threshold. Reuse
duplicate points across the finite-difference series. Repeat the applicable
tightened-SCF cases to quantify convergence scatter against the ASA allowance;
this is context, not a separate numerical gate. The finite-difference series
provides the truncation diagnostic; all energies and derivatives are
reported in Ry/eV with the input U/J recorded in eV.

If a standard case reaches its 240-step cap without literal convergence,
preserve its full history and final restart, then make one continuation with
only Broyden `beta` changed from 0.1 to 0.05, the same 1e-9 threshold and
physical settings, and a 240-step cap. If that continuation also fails to
resolve the decision-relevant energy scatter, stop and report the measured
floor before extending more cases.

The approved charge-only allowance is obtained from spin-selective W1 probes.
For spin-resolved occupations `q_s` and a separate W1 shift `lambda_a` applied
to spin `a`, measure `R_s,a=dq_s/dlambda_a` and
`B_a=dF/dlambda_a−q_a`, then solve `R^T g=B`. Here `g_s` is the residual
energy gradient per spin d occupation in Ry/electron. From the U scan measure
`dq_s/dU`; the base-functional allowance is
`sum_s |g_s dq_s/dU|` plus the propagated measured uncertainty, in Ry/eV.
For J, use the same registered spin-gradient and the centered per-spin
response `dq_s/dJ` at each prescribed center, giving the analogous
`sum_s |g_s dq_s/dJ|` bound and ordinary finite-difference uncertainty
propagation. The measured response stencils are the nonnegative J points
listed above, so both J centers have symmetric ±0.05 and ±0.10 eV responses.
This transfers the W1-inferred local charge-only gradient at U=J=0 to the
U=2,J=0.5/0.9 operating points; treat that cross-operating-point transfer as
an explicit approximation and report it as a limitation of the J allowance.
This discloses and uses the user-approved charge-only assumption. If the
response matrix cannot resolve `g`, report the measured uncertainty and leave
the b/c decision unevaluated. The finite-difference series provides the
truncation diagnostic; endpoint residuals are reported separately.

The earlier A5.1 eigen-reconstructed density matrix differs from the new
direct production P by about `6.76e-4` in the comparison already audited.
This is a state/definition diagnostic, not a gate. All new W1 responses,
occupation responses and fixed-P targets use the same direct production P
definition within each measured state.

Before inspecting any U/J candidate comparison, preregister the numeric
allowance from the independent spin-selective W1 response and the measured
per-spin `dq_s/dU`, including response-matrix and finite-difference
uncertainty. The prior W1 scale reference is `g0≈6.49848`, with aggregate
residuals `+0.0185609` spectral and `−0.00921488` local; these imply an
order-`1e-2` occupation bias before the chain-rule conversion, not a
Ry/eV allowance by themselves.

If b agrees and c misses by the predicted `dt/dU` under the chain-rule
allowance, apply the user-approved A5.3 adoption rule. If b fails, c agrees
within the allowance, or the bound cannot be formed, stop without changing
the default or `physical_total_energy`. Do not claim an a-versus-b preference
unless their difference exceeds the allowance.

Hashes will be recorded for the source files, executable, each input and
each output. The A5.3/A5.5 campaign shares one explicit U=J=0 KS restart
snapshot at 12³ and 300 K. The copied restart residual was `5.9e-9` under
the existing `1e-9` threshold, with no convergence marker. This common
starting state is not asserted to be converged; each candidate's final state
is assessed against the unchanged SCF threshold.
