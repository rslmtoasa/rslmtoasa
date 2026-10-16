# Proposal: resolve the leading (h^2) term while accounting for (h^4)

Status: **PROPOSED, NOT APPROVED**. This note does not change the frozen A4.1 rule, authorize inspection of response data, or release any nonzero-λ run. The current strict rule remains operative until the user approves a clarification.

## Why clarification is needed

The current h-scaling check requires the two adjacent differences to obey a pure 4:1 ratio within propagated measurement noise. That rejects a smooth response whenever a nonzero (h^4) term is resolved, even where the prescribed Richardson construction measures and removes its leading derivative bias.

For example, let

\[
D(h)=g+10h^2+32h^4,
\]

on the already prescribed grid (h=0.0025, 0.005, 0.01). Set (A=D(0.005)-D(0.0025)), (B=D(0.01)-D(0.005)). Then (A=1.8751875\times10^{-4}), (B=7.503\times10^{-4}), and \(|B-4A|=2.25\times10^{-7}). The latter is a resolved failure under a small roundoff-only bound, although the h⁴ contribution at the coarsest point is only 0.032% of the h² contribution. The same data give (L=6.25\times10^{-5}), (H=1.25\times10^{-9}), and (R_{fine}-g=-4H=-5\times10^{-9}), which is exactly the h⁴ bias the locked Richardson pair estimates.

## Proposed pre-response rule

Assume on the prescribed grid

\[
D(h)=g+b h^2+c h^4+O(h^6).
\]

Let (h=0.0025), (D_1=D(h)), (D_2=D(2h)), and (D_3=D(4h)). Define (A=D_2-D_1), (B=D_3-D_2). The leading and next-order contributions at the fine step are

\[
L=b h^2=\frac{16A-B}{36}, \qquad H=c h^4=\frac{B-4A}{180}.
\]

Propagate the already frozen per-(D_i) uncertainty bounds (u_i) conservatively, without assuming independent errors:

\[
u_A=u_1+u_2+u_{A,fp},\quad
u_B=u_2+u_3+u_{B,fp},
\]

\[
u_L=\frac{16u_A+u_B}{36}+u_{L,fp},\quad
u_H=\frac{u_B+4u_A}{180}+u_{H,fp}.
\]

The (u_{*,fp}) terms bound the actual floating-point operations used to form each difference and linear combination, from their operand magnitudes and operation counts. They must be locked before outcomes; they are not an absolute tolerance.

Require a resolved leading term and dominance of that term over the estimated h⁴ term at the coarsest prescribed point:

\[
|L|>u_L,
\qquad
|L|-u_L>16(|H|+u_H).
\]

The factor 16 follows from the h² and h⁴ scaling at (4h): their magnitudes there are (16|L|) and (256|H|). There is no selected percentage cutoff. The identity decision remains unchanged: require the original ( |R_{fine}-g_0|\le N_{unc}), with (N_{unc}) excluding measured identity bias. Failure of either the leading-term resolution/dominance check, the existing identity decision, or a separate bound that covers the neglected (O(h^6)) remainder is a stop/inconclusive result; no response-selected fit or threshold change is allowed.

This interpretation treats a resolved but sufficiently small h⁴ contribution as truncation measured by the existing Richardson pair. It does not presume that the h⁶ remainder is zero: the pre-response error model must justify its negligibility on this grid or classify the result as inconclusive. The three (D_i) values alone do not independently bound an arbitrary h⁶ term.

## User decision requested

Approve or reject this replacement for the strict pure-h² 4:1 check, keeping the existing Richardson derivative estimate, the original ( |R_{fine}-g_0|\le N_{unc}) identity gate, all frozen response values and grids, and the stop on unresolved/unsupported higher-order behavior. If not approved, the current strict rule remains in force and a resolved nonzero h⁴ mismatch remains a stop even when Richardson estimates its bias.
