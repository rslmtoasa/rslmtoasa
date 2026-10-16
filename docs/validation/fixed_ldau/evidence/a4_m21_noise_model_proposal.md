# A4.1 premeasurement proposal: finite-difference noise and identity decision

Status: reviewer operationalization, approved by the user for A4.1 identity
responses through the root coordinator. No A4.1 response outcomes have been
measured. This note does not amend the user-owned protocol or release any
lambda response; exact inputs and thresholds still require a pre-response
registration. The A4.2 transfer remains unresolved and is not approved here.

## Decision being protected

For each KS mode and each fixed operator `W`, compare the finite-difference
estimate of `dF/dlambda` at zero with the independent matrix contraction
`g0 = Tr(W P(0))`. Both quantities are dimensionless: `F` and `lambda` are Ry,
and `P` is the dimensionless occupation matrix. The A4.1 identity passes only
when the Richardson-extrapolated derivative agrees with `g0` within a bound
obtained independently from SCF repeat differences, finite-difference
truncation, and floating-point evaluation. The measured identity residual is
never included in that bound.

## Proposed measurements and error model

For each mode and each `W`, freeze `P(0)` and `g0` before the nonzero-lambda
responses. At each prescribed lambda (`0`, `+/-0.0025`, `+/-0.005`, `+/-0.01`
Ry), run the ordinary SCF deck and one otherwise-identical deck from the same
starting state with `self%conv_thr` divided by ten and, when fixed-U residual
monitoring is active, both R05 tolerances divided by ten. The production R05
values in `source/mix.f90` are `fixed_u_density_tolerance = 1e-7` and
`fixed_u_potential_tolerance = 1e-7`; the tightened values are each `1e-8`.
Require the existing three consecutive stable updates in both runs. Because
the R05 values are constants rather than namelist inputs, the pre-response
registration must bind the exact test-only mechanism that applies the 10x
scale and prove the standard run retains the unmodified values. All other
deck, binary, mesh, runtime, and seed inputs remain identical. If the test hook
does not activate fixed-U residual monitoring, record that fact; in that case
the ordinary `conv_thr` pair is the applicable tightening.

Let `F_s(lambda)` and `F_t(lambda)` be the same preregistered free-energy
candidate from the standard and tightened runs: in local mode use
`ETOT_local + lambda*Tr(WP) - TS`; in spectral mode use
`ETOT_spectral - TS`. Use the observed paired difference as a conservative
SCF sensitivity proxy, not as a statistical confidence interval:

```
uF(lambda) = abs(F_t(lambda) - F_s(lambda)) + 1e-8 Ry
uD(h) = [uF(+h) + uF(-h)] / (2 h)                 (dimensionless)
D(h) = [F_s(+h) - F_s(-h)] / (2 h)              (dimensionless)
```

The `1e-8 Ry` term is the existing replay energy-comparison guard, applied once
to each energy difference. For the independent target, compare the standard
and tightened zero-lambda matrices using

```
g0 = Tr(W P_s(0))
ug0 = abs(Tr(W [P_t(0)-P_s(0)])) + gamma_256 * sum_ij abs(W_ij P_s(0)_ji)
gamma_n = n*epsilon(rp) / (1 - n*epsilon(rp))
```

The second term is a conservative 5x5 complex trace arithmetic bound, with
`n=256` real floating-point operations. It has units of the derivative,
dimensionless. The matrix difference term is an observed standard/tightened
state sensitivity and remains separate in the report.

Use only the three prescribed central-difference spacings. Define the nested
Richardson estimates

```
R_fine   = (4 D(0.0025) - D(0.005)) / 3
R_coarse = (4 D(0.005)  - D(0.01))  / 3
uR_fine   = (4 uD(0.0025) + uD(0.005)) / 3
uR_coarse = (4 uD(0.005)  + uD(0.01))  / 3
B_trunc   = abs(R_fine - R_coarse)/15 + (uR_fine + uR_coarse)/15
N_unc     = uR_fine + B_trunc + ug0
```

The factor 15 follows from the leading `h^4` error of two adjacent
second-order Richardson estimates; no unprescribed `0.00125 Ry` response is
used. `N_unc` is dimensionless and is fixed without using the identity
residual. `1e-8 Ry` per energy is the only preexisting energy guard in `uF`;
the explicit matrix trace bound covers the arithmetic in `g0`.

To demonstrate the requested `h^2` central-difference scaling, let
`X=D(0.01)-D(0.005)` and `Y=D(0.005)-D(0.0025)`. Leading `h^2` behavior
predicts `X=4Y`. Propagate uncertainty as `uX=uD(0.01)+uD(0.005)` and
`uY=uD(0.005)+uD(0.0025)`. Scaling is resolved only if both differences are
larger in magnitude than their respective propagated bounds and
`abs(X-4Y) <= uX+4*uY`. If this condition is not resolved, classify the
identity result as inconclusive and stop; do not report a pass from a
noise-dominated curvature test.

The identity bias is reported separately as `B_W = R_fine - g0`. The A4.1
identity passes only if scaling is resolved and `abs(B_W) <= N_unc`; otherwise
it fails or is inconclusive and the protocol requires stopping. A nonzero
smooth bias is therefore allowed only when the independent repeat/truncation
budget bounds it. Report `B_W` separately; do not silently absorb it into the
A4.1 allowance.

The proposed independent A4.1 allowance is a derivative error in
`Ry/Ry` (dimensionless). A4.2 differentiates an `Ry` energy with respect to
`eV`-valued U or J, so its derivative/error is `Ry/eV`. The amendment does not
specify a transfer rule from the lambda-W error/bias to the U/J systematic
term, or whether the transfer should instead be made at the energy-response
level after converting the perturbation to Ry. Therefore this proposal does
not authorize using `N_unc + abs(B_W)` directly as A4.2's `N_W`. User
clarification is required to state the transfer rule, unit conversion, and
physical assumption that makes the W-derived systematic applicable to the
U/J derivative. No A4.2 measurement or decision may rely on that transfer
until it is separately preregistered.

## Candidate A4.2 norm-matched transfer for user review

At the frozen U/J reference state `P0`, define the Hermitian potential
derivative for each varied parameter `X` as

```
G_X = dV_U(P0; U,J)/dX             (Ry/eV)
alpha_X = ||G_X||_2 = max(abs(eig(G_X)))  (Ry/eV)
W_X = G_X / alpha_X                (dimensionless, max(abs(eig))=1)
lambda_scale = alpha_X * abs(delta_X)       (Ry)
```

Obtain `G_X` from an analytic derivative of the registered fixed-P
Liechtenstein operator. If an analytic derivative is unavailable, freeze a
fixed-P symmetric/Richardson evaluation and its arithmetic bound before any
U/J response; do not infer `G_X` from the response being tested. The potential
change over `delta_X` has spectral norm `lambda_scale`. If the A4.1 bound is a
dimensionless derivative allowance `N_W` measured in the exact normalized
direction `W_X`, the corresponding U/J derivative allowance is

```
N_X = alpha_X * N_W          (Ry/eV)
```

Equivalently, the energy uncertainty at that perturbation is bounded by
`N_W * lambda_scale` Ry, then divided by `abs(delta_X)` eV. Compare this
quantity to the A4.2 Richardson derivative and include it as the smooth
systematic term alongside separately measured A4.2 truncation and SCF-scatter
terms. Use the `alpha_X` and `lambda_scale` calculated from the locked `P0`
for every prescribed `delta_X`; bind the values and source formula before
responses. If `lambda_scale` falls outside the A4.1 tested lambda interval,
the transfer is unvalidated and must stop or receive a separately approved
extension; no extrapolation is implied.

The prescribed A4.1 operators `W1` and `W2` are not in general equal to the
normalized `W_U` or `W_J`. Thus `alpha_X*N_W` is justified only under an
additional direction-independent per-unit-operator-norm bound. Two tested
directions do not mathematically establish that universal bound. A more
direct alternative is to preregister and measure the A4.1 control for each
exact normalized direction `W_X`, using its own independent `N_W_X`; this
would avoid the direction-transfer assumption but adds controls beyond the
literal W1/W2 list. User clarification is required to choose between this
explicit norm-bound assumption and exact U/J-direction controls before any
A4.2 responses. A unit conversion by `1/13.605703976` alone is insufficient:
it changes parameter units but does not account for the operator norm or
direction of `dV_U/dX`.

## W2 definition and controls

Use the following fixed d-block matrix in the onsite spherical `m=-2..+2`
ordering (the same basis as the production d-block insertion):

```
W2 = [[ 0, +i, 0,    0,     0],
      [-i,  0, 0,    0,     0],
      [ 0,  0, 0.5,  0,     0],
      [ 0,  0, 0,   -0.25,  0],
      [ 0,  0, 0,    0,     0]]
```

It is Hermitian, anisotropic, has a nonzero imaginary off-diagonal, and has
eigenvalues `(+1,-1,+0.5,-0.25,0)`, so `max(abs(eig(W2)))=1` exactly. The
operator, basis convention, and lambda signs must be hash-bound before the
lambda responses. After the `+/-0.01 Ry` W2 states are captured, but before
their frozen dual-control replays, preregister the transpose and conjugate
control predictions from those immutable states and apply A3.1's existing
eligibility rule. If either control is ineligible, record it as such and use
the protocol's prescribed fallback/stop behavior; do not infer activation from
its measured response.

## Why user clarification is requested

A4.1 requires an extracted `N_W` but does not specify its estimator. The user
approved the non-circular identity gate described here: independent SCF,
Richardson, trace, and floating-point allowances; mismatch excluded from the
allowance; stop if mismatch exceeds the bound or `h^2` scaling is unresolved.
The registration must implement the approved 10x tightening for all applicable
SCF and R05 thresholds. A4.2 still calls for a smooth systematic and uses a
derivative with different parameter units; its transfer rule is not approved.
These are operational definitions, not measurements, and no lambda response
may run until its complete registration is frozen.
