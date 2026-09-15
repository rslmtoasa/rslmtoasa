# GSR-01 independent Goldstone sum-rule interaction

## Status and claim level

The independent Lounis/Costa/Muniz/Mills (LCMM) static sum-rule route is
implemented in [`source/lr_goldstone_sumrule.f90`](../source/lr_goldstone_sumrule.f90)
and registered in the `rslmto` library.  It consumes a canonical LR-04 static
`chiKS` and an explicitly labelled Pauli/no-SOC magnetization.  It does not
call or rescale the direct ALSDA implementation.

The focused test establishes algebraic consistency, independent numerical
fixtures, and the metric-aware static Goldstone identity.  It does not claim
converged Fe/Ni response, a dynamical Dyson result, or literature-spectrum
agreement.

## LR-03 mapping and discrete equation

The LCMM quantities are the energy-valued down-minus-up splitting and the
up-minus-down number-spin density:

\[
 B_{\rm eff}=V_\downarrow-V_\uparrow,
 \qquad
 U_{\rm LCMM}(r)=\frac{B_{\rm eff}(r)}{4\pi m_z(r)}.
\]

The normalized LR response uses the `L=0` coefficients

\[
 m_{00}(r_i)=\sqrt{4\pi}\,m_z(r_i),
 \qquad
 B_{{\rm eff},00}(r_i)=\sqrt{4\pi}\,B_{\rm eff}(r_i).
\]

Consequently the exact LR-03 conversion is

\[
 B_{{\rm eff},00}(r_i)=4\pi U_{{\rm LCMM}}(r_i)m_{00}(r_i).
\]

LR-04 stores a pointwise susceptibility as the canonical right-weighted
operator `chi_can=chi_raw W`.  Its action already performs the radial source
contraction:

\[
 y_I=\sum_J(\chi_{\rm can})_{IJ}x_J.
\]

For the published spherical approximation, the unknown is therefore a
site/radial scalar `U_LCMM(a,r)`, restricted to the positive-measure `L=0`
coordinates.  Let `p(i)` be the corresponding LR-04 flat index.  Substitution
of `B_eff,00=4*pi U_LCMM m_00` into the static identity gives the active
production equation

\[
 \boxed{
 \Gamma_{ij}U_j=m_{00,i},
 \qquad
 \Gamma_{ij}=4\pi\,
 (\chi_{\rm can})_{p(i),p(j)}m_{00,j}.
 }
\]

This is the actual `Gamma U=m` equation solved by the implementation.  The
factor `4*pi` is retained because LCMM's `U` contains the published
`1/(4*pi)` definition; it is not inserted or removed by fitting a numerical
Goldstone error.  No extra `W_j`, `r_j^2`, or angular quadrature factor is
inserted: `chi_can` already contains the LR-04 source metric and `m_00`
already contains the normalized `Y_00` conversion.

The interaction passed to a later LR-04 Dyson layer is the canonical local
operator

\[
 K_{\rm eff}=4\pi U_{\rm LCMM},
 \qquad
 B_{{\rm eff},00}=K_{\rm eff}m_{00}.
\]

It is stored as `canonical_interaction`; `u_lcmm` remains the published LCMM
quantity and is not silently relabelled as a direct ALSDA kernel.

## Origin and response representation

The radial mesh includes `r=0`, where LR-04's physical volume weight is exactly
zero.  The solve excludes those origin columns and uses one unknown for every
site and positive-measure radial point.  The origin value of `U_LCMM` and
`B_eff` is the zero extension, and it cannot affect a response contraction.
No `1/W(0)`, epsilon floor, radial cutoff, or origin regularization is used.

The result remains a spherical angularly diagonal radial interaction.  The
current baseline requires one certified circular response channel and uses the
`L=0`/`L=0` static block.  If a spherical source generates a nonzero active
`L>0` response, the full LR-04 metric residual is reported and the result is
`BLOCKED`; the implementation does not hide that failure by promoting a
general response-space fit to the published local approximation.

## Stable solve and blocker policy

`solve_lr_goldstone_equation` uses the complex LAPACK SVD least-squares kernel
(`ZGELSS`) with `RCOND=-1`, meaning LAPACK's machine-precision rank rule.  This
is a numerical rank diagnostic, not a user-selected smoothing cutoff.  The
result reports:

* all singular values;
* numerical rank and the equation dimension;
* the 2-norm condition number for a full-rank equation;
* the Euclidean equation residual and relative residual; and
* the LR-04 metric residual, relative residual, and rigid-vector overlap of
  the generated static response.

A rank-deficient equation is marked `BLOCKED`, even though LAPACK returns its
minimum-norm diagnostic vector.  LCMM does not select a null-space interaction,
so that vector is not promoted to production status.  A failed equation or
full-response Goldstone residual is handled the same way.  No Tikhonov
regularization, empirical rescaling, pseudoinverse cutoff, or direct-ALSDA
fallback is present.

## Public API

The request form keeps the three source contracts explicit:

```fortran
type(lr_goldstone_sumrule_request) :: request
type(lr_goldstone_sumrule_result) :: result

request%response_space => space
request%static_susceptibility = chi_static_canonical
request%magnetization = pauli_m
request%magnetization_label = lr_gsr_magnetization_pauli
call evaluate_lr_goldstone_sumrule(request, result)
```

The explicit overload accepts
`(space, static_susceptibility, magnetization, result, magnetization_label)`.
The result carries `u_lcmm`, the pointwise `effective_field`, the canonical
`canonical_interaction`, `gamma`, singular values, rank/conditioning, equation
and metric residuals, generated response, overlap, and `blocked/status`
metadata.

## Static Goldstone identity

After solving `Gamma U=m_00`, the implementation builds
`K_eff=4*pi*U_LCMM` through the LR-04 complex local-operator helper and applies
it to `m_00`.  It then applies the supplied static canonical `chiKS` and
evaluates

\[
 r=\chi^{KS}(0)K_{\rm eff}m_{00}-m_{00},
 \qquad
 \|r\|_W=\sqrt{r^\dagger W r}.
\]

The full residual is used for the status decision.  The rigid overlap is the
LR-04 metric overlap with `m_00`; multiplying both vectors by the common
`-i theta^+` phase gives the equivalent LR-03 rigid-rotation vector.

## Test evidence

`UnitLrGoldstoneSumrule` is a standalone CTest target.  Its independent
fixtures cover:

1. an analytic complex `Gamma/U/m` solve;
2. a one-site radial local interaction;
3. multisite unequal accepted moments;
4. an exact manufactured nonlocal rigid-rotation identity;
5. a deliberate nonspherical response to verify the LR-04 metric residual and
   `BLOCKED` status;
6. an explicit singular-matrix blocker; and
7. a direct `B_eff/(4*pi*m_z)` comparison diagnostic that is printed but never
   used to construct `Gamma` or the solution.

The focused run reported:

```text
analytic Gamma/U/m solve max error =   2.9894E-16
Analytic solve condition/residual/rank =   2.1492E+00   7.9770E-16 3
Rank-deficient solve rank = 1, status = BLOCKED: SVD numerical rank is below equation dimension
one-site radial U max error =   1.6686E-16
one-site metric Goldstone residual error =   1.3732E-16
One-site metric residual/condition =   4.3672E-19   1.5862E+00
multisite unequal-moment U max error =   1.1114E-16
Multisite metric Goldstone residual =   1.6296E-18
manufactured nonlocal rigid U max error =   5.0179E-16
Manufactured identity residual/overlap =   3.9356E-18   1.9816E-20
metric-aware full-response residual error =   0.0000E+00
Metric-aware residual/identity status =   3.3956E-03 BLOCKED: full LR-04 Goldstone identity residual exceeds tolerance
Direct ALSDA comparison diagnostic (not used in solve) =   3.9197E-01
UnitLrGoldstoneSumrule: PASS (analytic, radial, multisite, rigid, metric, independence)
```

Commands:

```text
cmake --build build --target UnitLrGoldstoneSumrule -j2
ctest --test-dir build --output-on-failure -R '^UnitLrGoldstoneSumrule$'
```

The focused CTest result was 1/1 passed.  The evidence is algebraic and
finite-fixture numerical evidence only; a dynamical Dyson susceptibility and
material validation remain downstream tasks.

## Scope boundary

This task adds only the independent static LCMM interaction route.  It does
not implement a dynamical Dyson solve, loss matrix, Goldstone eigenvalue
correction, mode extraction, site-only reduction, noncollinear response, or a
new ground-state/XC functional.

## TDVK-06 independent compact GSR adapter

TDVK-06 adds an independent compact implementation for the accepted Fe
static state.  It does not call the direct ALSDA kernel.  The accepted Pauli
magnetization is projected to `m00=sqrt(4*pi)*m_pauli`, while each unknown
`U_LCMM` is represented in the retained L=0 product span.  Its field basis is
constructed by explicit point-space reconstruction and projection, then the
raw compact equation

```text
Gamma * u = m00,
Gamma(:,j) = chiKS_compact(0,eta) * [4*pi*U_j*m00]
```

is solved with LAPACK `ZGELSS` and `RCOND=-1`.  The GSR columns use the
compact-representable source `m00_compact_point=R*c`, so each column is
`f_j=P*K_j*R*c=Kc_j*c`, exactly the action used by the final compact operator.
The complete compact response residual is evaluated as
`chiKS_compact * K_eff * m00 - m00`, and is compared with the actual linear
system residual.  Rank, singular values, condition, coefficient norm, action
consistency, equation residual comparison, and status are all written.  Rank
deficiency, action inconsistency, or a residual above `1e-9` is reported as
`BLOCKED`; the returned least-squares vector is never promoted through
regularization or a null-space choice.

The point/product mapping is the certified TDVK-06 contract
`c=U^H sqrt(W)x`, `x=inv(sqrt(W))Uc`, and `Kc=U^H K U`.  The unit oracle and
the strict 232-mode Fe runtime inventory certify representation use only; they
do not certify a compact Dyson denominator or any dynamical/magnon claim.
