# A5.2 local W1 control preregistration and evidence

This is the concise preregistration for the A5 diagnostics owned here. The
first shared no-U attempt used a full trace, stopped at about iteration 43, and
did not converge; it is partial raw evidence and is not an accepted endpoint.
The authorized common initial restart for A5.2 is explicitly near-converged,
not an accepted endpoint: `/private/tmp/a5-common-seed-near/Fe_out.nml`
(SHA256 `25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5`)
with `/private/tmp/a5-baseline-240/input.nml` (SHA256
`ca65d9695b049f10beb3ee9461c05572b48f3973885de99ea7f8fee3dfc5699c`). Each
case gets its own copy of this same input/seed; every resulting SCF endpoint
must still meet existing convergence criteria. The instrumented Release binary hash is
`52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`.
Modified source hashes are: `hamiltonian.f90`
`db26fbc766f94c811ae641c0fe2dd918bb913be3fce92f4488b9af826ce29b2f`,
`hamiltonian_build.f90`
`e896d53d82560a0551cc44560191dda5a61a70418a6b75e5c6247c758cd7ca90`,
`hamiltonian_hubbard.f90`
`8dcab0b0fa467b0e8532877e3e7e00a96cef7e47dce38a4ac085dd26dada43ee`,
the Hamiltonian namelist include
`afde6c79db8fe06d8d5bd61fbb90720ee4990c9337b6b4efc530da1297edfe70`, and
`self.f90` `1c29e173c7ae0061f2a9087cc90fbe44db0a24000dfda183340a914d8e20a7ad`.
Run input and output hashes are captured in the associated receipts. A4
evidence remains unchanged.

The test-only HOH bypass boundary receipt is
[`a5_hoh_guard_boundaries.json`](a5_hoh_guard_boundaries.json), captured with
the final instrumented binary (SHA256
`52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`). All
seven runtime cases passed: the authorized W1/local/U=J=0 HOH case returned
zero; bypass-off, nonzero U, nonzero J, W2, no-W, and spectral negatives each
returned nonzero before SCF iteration markers. The receipt distinguishes the
pre-run seed hash from runtime `Fe.nml` and `Fe_out.nml` output hashes.

The same binary passed the actual one-step native-KS moments-only smoke,
recorded in [`a5_moments_only_smoke.json`](a5_moments_only_smoke.json). It
captured physical `P_OUT` d-shell density matrices for both spins, q0 and
applied `Tr(WP)` rows, ql/radial moments, and ASA energy fields while excluding
onsite H, radial replay vectors, eigenpair, Green-function, and bulk-bond trace
kinds. The CSV trace was 937,750 bytes; its hash and pre-run input/seed hashes
are in the receipt.

## A5.2: local W1 non-stationarity controls

Use the common near-converged KS no-U restart snapshot, scalar-collinear
`nsp=1`, `lmax=2`, local radial moments, W1, and the common 12³/300 K base
setup. Each control's final SCF state must still meet the unchanged threshold.
For each control,
evaluate the W1 free-energy slope at `lambda = 0, ±0.005 Ry`, and compare the
centered finite-difference derivative with `Tr(WP(0))`. Report the signed
residual in absolute occupation units and relative to `Tr(WP(0))`. The
previous local W1 result gives a scale reference only: `g0≈6.49848`,
`R_fine−g0≈−0.0092149` (about −0.14%); its finite-difference uncertainty was
at most `5e-5`.

The three controls change one factor at a time: (a) enable HOH through a
test-only W1/local/U=J=0 path, (b) compare 150 K and 600 K with the 300 K base,
and (c) compare 20³ with the 12³ base k mesh. The HOH run uses the
source-reviewed bypass that is restricted to W1, local moments, and exactly
U=J=0; it must not weaken any production guard. The established spherical
d-shell basis ordering remains `m=-2..+2`.

At each control endpoint, compare the centered free-energy derivative
`[F(+h)-F(-h)]/(2h)` (`h=0.005 Ry`, `F=test_F_local`) with the direct
`Tr(WP)` captured from the actual `P_OUT` density at λ=0. The trace output
contains each physical spin block in the established `m=-2..+2` ordering, its
charge trace, and its applied per-spin `W:P` contraction (zero for an unperturbed
spin). The finite-difference derivative is evaluated with W applied to both spins;
the separate spin-selective chain-rule probes are preregistered by the campaign
owner. This compares the derivative to the same endpoint state, without
reconstructing P from a paired-λ average.

The A5.2 comparisons use the direct production `P_OUT` density matrix in every
case. The historical A4 eigen-reconstructed density differs from this direct
matrix by about `6.76e-4`; that is retained only as a provenance diagnostic and
is not an A5.2 gate or an input to its derivative estimate.

Prediction and expected scale: if the first-order Hamiltonian causes the
residual, HOH reduces its magnitude by about 10×, from order `1e-2` to order
`1e-3` occupation. If smearing dominates, the low-temperature residual is
about 1/4 of the base and the high-temperature residual about 4× the base,
consistent with `T²` scaling. If k integration dominates, 20³ changes the
residual relative to 12³; if k integration accounts for most of the observed
non-stationarity, expect a change of order `1e-2` occupation, while a smaller
change indicates a partial contribution. No direction is preregistered. A control effect
smaller than the combined numerical resolution is unresolved, not evidence
for that mechanism.

Error budget and decision: resolve effects against same-state repeat scatter,
finite-difference truncation, direct trace roundoff, and the SCF endpoint
residual. The existing A4 local W1 bounds put finite-difference uncertainty at
`≤5e-5` against a residual of at least `9e-3`; thus FD alone is sub-percent of
the observed residual. Also report 16-digit log rounding and SCF convergence.
The direct λ=0 trace avoids a separate state-reconstruction uncertainty. The expected HOH signal (`~9e-3`),
temperature shift (`~7e-3` low or `~3e-2` high under pure `T²` scaling), and
mesh signal are interpreted only when they exceed the combined uncertainty.
These are diagnostic outcomes: they identify a likely source of the base
functional non-stationarity and do not change production defaults or physics.

Convergence contingency preregistration: if the first base `lambda=-0.005 Ry`
case reaches its 240-step cap without a `Converged!` marker, preserve its full
stdout, trace, all runtime restart files, and hashes, then perform one bounded
input-only continuation from its latest `Fe_scf.nml`. Change only the mixing
`beta` from 0.1 to 0.05; keep the same `conv_thr=1e-9`, Hamiltonian, mesh,
temperature, λ, and physics, and allow at most 240 continuation steps. Record
the continuation input and restart pre-run hashes and the mutated/restart
output hashes separately. If that continuation does not converge, stop the
remaining local-W1 matrix, report the observed residual/noise floor and the
unassessed controls, and do not repeat the 240-step continuation on every
case. If it converges, use the same mixing setting and convergence criterion
for the remaining A5.2 cases, each from an independent copy of the original
common initial seed.
