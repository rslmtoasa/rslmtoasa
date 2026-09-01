# TDDFT ground-state Ward identity and Goldstone diagnostics

This document records the TDDFT-02 kernel contract implemented in the
current branch.  It is deliberately limited to the transverse, collinear,
SOC-free static identity.  A nonzero residual in a calculation with SOC,
external fields, constraints, noncollinear magnetism, or an unsupported
Hamiltonian representation is not silently interpreted as a Goldstone error.

## Physical identity

The ground-state exchange-correlation magnetic field and the static
transverse kernel are represented in the same active response basis:

\[
 B_{xc}=K_{xc}m, \qquad
 W=\chi_{KS}(q=0,\omega=0)B_{xc}-m,
\]

so that the equivalent Goldstone condition is

\[
 Dm = [I-\chi_{KS}K_{xc}]m=0.
\]

The code reports both normalized residuals,

\[
 r_W=\frac{\|W\|_2}{\|m\|_2},\qquad
 r_D=\frac{\|Dm\|_2}{\|m\|_2},
\]

before any repair.  It also writes the magnetization norm, the `B_xc` norm,
the response basis, and the provenance of both `B_xc` and `K_xc`.

The circular convention is the one shared by `tddft_conventions_mod`:
`m+/- = mx +/- i my`, with unhalved `sigma_x +/- i sigma_y` response
operators.  Consequently the source and measurement factors are kept
separate (`1/2` source factor and `2` circular measurement factor); the Ward
residual does not introduce another ladder-operator factor.

## Ground-state kernel provenance

The production physical transverse operator is obtained from the same
ground-state representation used for the response:

* the XC provider is refreshed through the normal `VXC0SP`/SCF path;
* the site magnetization is reconstructed from occupied spinor eigenvectors
  as the signed `P_site sigma_z` population;
* the direct pair-potential route uses the analytic transverse rotation of
  the ordinary LMTO `ham_only` operator and builds Xi directly;
* the site-projected ALSDA scalar is retained as a documented baseline and
  diagnostic, not inferred from `hxc`, `cx1`, or a finite-eta inverse of
  `chi_KS`.

The direct pair-potential static diagnostic is labelled
`direct LMTO ham_only transverse pair-potential Xi` and uses the same signed
site moments.  Generalized-overlap, HOH, GBT, CCOR, Hubbard, SOC, and
noncollinear representations remain capability-gated where the corresponding
transverse tangent has not been derived.

## Explicit repair policies

The default is `goldstone_policy='diagnose'`.  It records the raw identity and
does not modify the kernel.  The repair routines are separate calls in
`tddft_ward_mod` and retain their raw diagnostic record.

`sum_rule` calls `reconstruct_lounis_kernel`.  In the site-diagonal active
basis it solves the real static system

\[
 \operatorname{Re}(\chi_{KS})\operatorname{diag}(m)k=m,
 \qquad K_{xc}=\operatorname{diag}(k),
\]

using a rank-revealing SVD.  It rejects material static imaginary content,
small moments, rank deficiency, and excessive conditioning.  The resulting
kernel is not an empirical global scale; the change from the physical kernel
is reported explicitly.  A warning or failure threshold can be supplied by
the caller.  There is no hard-coded universal percentage threshold.

`projected` calls `project_goldstone_eigenvalue`.  It selects the closest-to-
one eigenvalue only when its right eigenvector is magnetization-like and its
distance from one is small.  With right/left eigenvectors `v/u`, it applies
the single biorthogonal rank-one update

\[
 \Xi' = \Xi + (1-\lambda)\frac{v u^\dagger}{u^\dagger v}.
\]

The other eigenmodes are not globally rescaled.  The operator change and the
raw/corrected Ward residuals are written separately.  The corrected static Xi
is converted back to the adiabatic kernel by a linear solve when the explicit
policy is used by the legacy active-basis driver.

## Deterministic regression evidence

These are small algebraic fixtures, not Fe/Ni material claims.  Values below
are normalized with the Euclidean norm of the signed response magnetization.

| fixture | raw residual | repaired residual | result |
| --- | ---: | ---: | --- |
| exact two-site `chi_KS B_xc=m` | 0 | — | pass |
| site Xi, `chi=diag(0.3,0.2)`, `K=diag(3,4)`, `m=(2,1)` | 0.1264911064 | — | raw value retained |
| Lounis reconstruction from physical `K=(1.8,4)` | 0.0894427191 | below test tolerance | pass; warning field exercised |
| Halle fixture `Xi=diag(0.9,0.2)`, `m=(1,0)` | 0.1 | below test tolerance | pass; second eigenvalue stays 0.2 |

The executable `UnitTddftWard` covers these cases.  `UnitTddftGoldstone`,
`UnitTddftWardConventions`, `UnitResponseConventions`, and
`UnitTddftConfig` cover compatibility, circular normalization, provider
kernel construction, signed moments, and explicit policy input.

## Production convergence status

The production diagnostic is emitted in the Goldstone output before any
optional correction and includes `ward_residual`, `dm_residual`, basis, and
provenance.  The response basis decision is currently:

* `site` is the validated production basis for the ALSDA site projection;
* the direct LMTO pair-potential route is the richer operator-level physical
  route for the supported `ham_only` representation;
* a site-orbital kernel requires independent orbital-resolved ground-state
  `B_xc` and moment data.  The response channel type can name such a basis,
  but the implementation does not claim site-orbital convergence from a
  scalar site projection.

Fe and Ni k-mesh/energy-window convergence tables are intentionally not
claimed by this implementation pass.  The existing material gates in
`docs/dev/TDDFT-07_VALIDATION.md` remain open/failed and must be rerun after
the corresponding fixtures and production settings are available.  The
deterministic unit evidence above is therefore the complete new numerical
evidence for TDDFT-02 in this branch.

## Output and audit requirements

Raw records are never replaced by corrected records.  Goldstone output uses
the identity header
`chi_KS(0,0) B_xc - m` and retains the `D m` residual.  Explicit repairs write
their policy, decision, conditioning, correction magnitude, and corrected
residual.  No empirical global `lambda`, energy shift, or finite-eta
inverse-kernel substitution is used.
