# LR-04 response-space metric and operator algebra

## Status and scope

This document locks the finite-dimensional algebra used by the direct radial
response space. It is an algebraic implementation contract only. It contains
no transition-vector evaluator, `chiKS`, XC kernel, Dyson equation, Goldstone
correction, or material validation.

This implementation was started on branch `fable_v4` at exact HEAD
`590da79c149d58cd9408ee3c2eb9ca7419d50246`. The supplied post-LR-03 prompt pack
was present as untracked task input; production changes are limited to the
module, its unit test/build registration, and this contract.

The implementation is [`source/lr_response_space.f90`](../source/lr_response_space.f90),
registered in the `rslmto` library. It reuses the LR-02R super-index and radial
mesh helpers from `response_basis_mapping_mod`.

## Coordinates and metric

The response coordinate is the established LR-02R order

\[
 I=(a,L,M,i,\mu),
 \qquad
 a=1\ldots N_s,\quad L=0\ldots L_{\rm max},\quad
 i=1\ldots N_r,\quad \mu=1\ldots N_\mu,
\]

with `(L,M)` in the live complex-harmonic order and site/radial/channel
flattening supplied by `response_basis_mapping_mod`. A response vector stores
pointwise harmonic coefficients; it does not store legacy weighted `RHO`.

For the LR-01 logarithmic mesh,

\[
 r_i=B[\exp(A(i-1))-1],
 \qquad
 \frac{dr}{di}=A(r_i+B).
\]

The radial metric used for a physical pointwise field is

\[
 W_i=s_i\,r_i^2\,A(r_i+B),
\]

where

\[
 s_1=s_{N_r}=\frac13,\qquad
 s_i=\begin{cases}4/3,&i\text{ even},\\2/3,&i\text{ odd},\end{cases}
\]

for interior points. The full super-index metric is diagonal,

\[
 W_{IJ}=\delta_{IJ}W_{i(I)}.
\]

There is no extra `4\pi` in this response metric: the normalized spherical
harmonics have already handled the angular projection. A spherical physical
scalar `n(r)` is represented in the `L=0` harmonic coefficient as
`sqrt(4*pi)*n(r)` when that conversion is needed. The legacy quantity
`RHO=4*pi*r**2*n` uses a different storage convention and must not be used as
the response vector or substituted for `W`.

For complex response vectors,

\[
 \langle x,y\rangle_W=x^\dagger W y
 =\sum_I x_I^*W_{i(I)}y_I,
 \qquad
 \|x\|_W=\sqrt{\operatorname{Re}\langle x,x\rangle_W}.
\]

The vector units inherit LR-03: number-density response fields are in
electrons bohr`^-3` (with the corresponding harmonic coefficient units),
radial `W` has bohr`^3`, and a contraction has the physical volume measure.
An operator has the units needed to map its input field to its output field;
the canonical right weighting itself carries the radial volume factor.

## Canonical stored representation

The production matrix is the right-weighted representation

\[
 \boxed{B_{IJ}=A_{IJ}W_J},
\]

where `A` is a raw pointwise kernel. Its action is therefore

\[
 y_I=\sum_J B_{IJ}x_J.
\]

This choice follows from the physical contraction

\[
 y_a^{LM\mu}(r_i)=
 \sum_{bL'M'\nu j}A_{ab}^{LM\mu,L'M'\nu}(r_i,r_j)
 W_jx_b^{L'M'\nu}(r_j).
\]

It keeps pointwise fields as the vectors consumed by later vertices, makes the
identity and local pointwise actions exact finite matrices, and makes operator
composition unambiguous. Since the vector metric remains explicit, no later
consumer may insert an additional `W` into multiplication of canonical matrices.

The raw/canonical conversions implemented by
`response_raw_to_canonical` and `response_canonical_to_raw` are

\[
 B_{IJ}=A_{IJ}W_J,
 \qquad
 A_{IJ}=B_{IJ}/W_J\quad(W_J>0).
\]

The origin columns have `W_J=0`. A finite quadrature kernel uses the exact
convention `B(:,J)=0` for those columns, and the inverse conversion returns a
zero raw origin column. A canonical matrix with a nonzero origin column is
rejected by the inverse conversion because it has no finite raw pointwise
kernel under this quadrature. Point-local operators at the origin are
canonical-native and are not forced through a singular raw delta matrix.

## Origin treatment

The production mesh includes `r_1=0`, so `W_1=0` exactly from the physical
volume factor. No epsilon, `1/W_1`, or guessed delta normalization is used.

The full stored vector retains the origin value because local fields are
pointwise objects. The metric is a seminorm on that full array: all origin
entries are a null subspace. The positive-measure active response dimension is
stored as `response_space_layout%active_dimension` and contains all
site/harmonic/channel entries at radial points `i=2...N_r`.

For a quadrature integral operator, an origin input cannot affect an active
output. The adjoint helper checks this condition. Origin output values are
allowed and remain available for pointwise evaluation, but they do not enter
metric contractions or the trace. This is the mathematically clean active-grid
interpretation of the singular origin weight.

## Operator rules

For canonical matrices `B` and `C`,

\[
 (B\circ C)x=B(Cx),
 \qquad
 B\circ C=BC.
\]

This ordinary multiplication is valid because each operand already contains
the metric on its source column. Reconstructing raw kernels and multiplying
them with ordinary dense multiplication would be wrong; raw composition is

\[
 (A\circ C)_{IK}=\sum_J A_{IJ}W_JC_{JK}.
\]

The identity is `I_{IJ}=delta_{IJ}` in canonical storage, including its exact
pointwise action at the origin. For a local radial scalar `K_a(r_i)`,

\[
 (Kx)_{aLMi\mu}=K_a(r_i)x_{aLMi\mu},
 \qquad
 B^K_{IJ}=\delta_{IJ}K_a(r_i).
\]

Thus a local operator is diagonal in site, `L`, `M`, radial coordinate, and
channel. The rank-2 API accepts `K_a(r_i)` shared by all channels; the rank-3
API permits a separate scalar for each channel. No angular quadrature is
introduced. Both APIs have real and complex overloads; the latter preserves a
finite-broadening complex static interaction without changing the pointwise
operator rule.

The metric adjoint is fixed by

\[
 \langle x,By\rangle_W=\langle B^\ddagger x,y\rangle_W.
\]

For positive-measure indices,

\[
 (B^\ddagger)_{IJ}=\frac{W_J}{W_I}B_{JI}^*.
\]

The null-origin block uses the ordinary conjugate transpose and active/null
mixing is zero. This selects one definite representative of the adjoint on
the metric null space while preserving the defining contraction identity.

The trace diagnostic is the positive-measure trace

\[
 \operatorname{tr}_W(B)=\sum_{I:W_I>0}B_{II};
\]

origin point values are excluded. Rigid-rotation vectors are not constructed
by LR-04; `response_rigid_vector_norm` and
`response_rigid_vector_overlap` apply exactly the same metric to a supplied
LR-03 rigid vector.

## Angular and channel structure

The angular basis remains the LR-02R complex basis and the response cutoff is
the existing product-space cutoff. For a spherical ground-state local scalar,
analytic harmonic projection makes the operator diagonal in site, `(L,M)`, and
radial point, with channel diagonal when a channel-resolved scalar is supplied.
There is no response-space angular quadrature and no site-only reduction.

## Test evidence

`UnitLrResponseSpace` is a standalone Fortran unit test registered by CMake. It
checks:

- analytic radial fields contracted with the explicit logarithmic volume metric;
- exact zero origin weight and active-grid dimension;
- identity action and positive-measure identity trace;
- pointwise local multiplication, including at the origin;
- multisite separation and spherical diagonal structure;
- raw-to-canonical-to-raw round trip on the active source grid;
- composition against an independently written nested-loop quadrature oracle;
- the metric adjoint relation;
- rigid-vector norm and overlap helpers.

The composition oracle does not call the production matrix-composition helper.
The test is algebraic/numerical evidence only; it is not a susceptibility,
Goldstone, converged-material, or literature-agreement claim.

The focused verification completed with:

    cmake --build build --target UnitLrResponseSpace -j2
    ctest --test-dir build --output-on-failure -R '^UnitLrResponseSpace$'

Result: `UnitLrResponseSpace` passed (1/1 test, 0.02 s in the local build).
The full configured CTest unit label also passed: 102/102 tests.

## Unsupported representations and scope boundary

The module does not implement the raw unweighted matrix as the production
representation, a symmetric `W^{1/2}AW^{1/2}` representation, or any
origin-regularized alternative. The symmetric representation is intentionally
not used because the production origin weight vanishes and no arbitrary
regularization is justified. It also does not implement transitions,
`chiKS`, XC physics, Dyson algebra, or any site-only fallback.
