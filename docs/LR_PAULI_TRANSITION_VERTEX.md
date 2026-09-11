# LR-05 Pauli transition-vector evaluator

## Status and claim level

The production one-electron Pauli/no-SOC transition-vector primitive is
implemented in [`source/lr_pauli_transition_vertex.f90`](../source/lr_pauli_transition_vertex.f90)
and registered in the `rslmto` library. It emits the complete direct response
vector in the LR-04 coordinate order

\[
 I=(a,L,M,i,\mu),
 \qquad a=1\ldots N_s,\quad L=0\ldots L_{\rm max},
 \quad i=1\ldots N_r,\quad \mu=1\ldots N_\mu.
\]

This establishes the LR03/LR04 algebraic capability that production Pauli
transition vectors exist in the certified response space. It does not establish
susceptibility correctness, an XC kernel, Goldstone compliance, a magnon
spectrum, or material/literature agreement.

## Public API

The module is `lr_pauli_transition_vertex_mod`.

`pauli_endpoint_state` contains one eigenvalue and one complete production
coefficient vector. The vector is packed exactly as the reciprocal eigensystem:

```text
(site 1, spin-up orbitals, spin-down orbitals,
 site 2, spin-up orbitals, spin-down orbitals, ...)
```

Within one spin block the live orbital order is
`(s),(p,-1:1),(d,-2:2),...`. The `initialize` type-bound procedure copies the
caller-owned coefficients; the evaluator never rephases them.

`pauli_vertex_capabilities` carries the representation tuple checked at the
public evaluator boundary:

```text
reciprocal_mode       = ham_only
hamiltonian_order     = second
orthogonal            = true
collinear             = true
has_soc               = false
has_extra_operator    = false
```

The radial basis array must contain one complete
`lmto_radial_basis` per response site, with a common certified direct mesh and
`sp` or `spd` angular basis. The evaluator calls each basis object's existing
`require_supported` method itself. An unsupported tuple therefore cannot be
silently accepted by a caller that supplies compatible-shaped arrays.

For one response channel, call:

```fortran
call evaluate_pauli_transition_vertex(space, radial_bases, left_state, right_state, &
   gamma, capabilities, transition_vector)
```

Here `gamma` is an explicit complex `2 x 2` Pauli-space matrix and the response
layout has `nchannel=1`. For several response channels, the rank-3 overload
accepts `gamma(2,2,nchannel)` and emits all channels in the same vector. The
module also supplies convention-owned constructors
`pauli_charge_matrix`, `pauli_sigma_x_matrix`, `pauli_sigma_y_matrix`,
`pauli_sigma_z_matrix`, `pauli_sigma_plus_matrix`, and
`pauli_sigma_minus_matrix`. Circular channels are selected by these matrices,
not by a magic integer.

The evaluator has no occupation, frequency, denominator, susceptibility, XC,
Dyson, or empirical normalization arguments.

## Exact evaluated equation

For the certified large-component Pauli projection,

\[
 \Psi^P_{n\mathbf k,a}(\mathbf r)
 =\frac1r\sum_{lm\sigma}
 c^{n\mathbf k}_{alm\sigma}U^{n\mathbf k}_{al\sigma}(r)
 Y^{\rm code}_{lm}(\hat r)|\sigma\rangle,
\]

with first-order production augmentation

\[
 U^{n\mathbf k}_{al\sigma}(r)
 =G_{al\sigma}(r)+
 (\epsilon_{n\mathbf k}-E^{\rm work}_{\nu,al\sigma})
 \dot G_{al\sigma}(r).
\]

The emitted transition coefficient is

\[
\boxed{
 T^{P,\mu}_{nm;aLM}(r_i)
 =\frac1{r_i^2}
 \sum_{\substack{lm,l'm'\\\sigma\sigma'}}
 (c^{n\mathbf k}_{alm\sigma})^*c^{m,\mathbf k+\mathbf q}_{al'm'\sigma'}
 U^{n\mathbf k}_{al\sigma}(r_i)U^{m,\mathbf k+\mathbf q}_{al'\sigma'}(r_i)
 \Gamma^\mu_{\sigma\sigma'}
 {\cal G}^{LM}_{lm,l'm'} .}
\]

The production `response_gaunt` utility supplies

\[
 {\cal G}^{LM}_{lm,l'm'}
 =\int d\Omega\,(Y^{\rm code}_{LM})^*
 (Y^{\rm code}_{lm})^*Y^{\rm code}_{l'm'}.
\]

All response harmonics through the configured complete product space are
emitted. For `sp`, `Lmax=2`; for `spd`, `Lmax=4`. If a larger response cutoff
is supplied, the additional product channels are emitted as zero.

## Radial, angular, and origin normalization

The evaluator reuses `phi_large`, `phidot_large`, and `enu_work` from the
production `lmto_radial_basis`. It never reads `phi_small`, `phidot_small`,
`gfac`, or `phiddot_*`. Consequently this is the documented Pauli projection,
not the exact scalar-relativistic bilinear.

The transition vector is pointwise. It does not absorb the LR-04 radial metric

\[
 W_i=s_i r_i^2 A(r_i+B),
\]

and it does not contain an additional `4*pi`. The metric enters only in later
response contractions. The `1/r_i^2` in the vertex is the physical wavefunction
factor and cancels the `r_i^2` in a later volume contraction as required by the
LR02R mapping.

The production mesh contains `r_1=0`. No epsilon radius is introduced. Regular
LMTO numerators have `U_l=O(r^{l+1})`, so the origin limit is zero when
`l+l' > 0`. The only possibly nonzero limit is `s-s`; it is obtained by the
same two-positive-point regular extrapolation used by the accepted radial
origin-density contract. The origin remains a pointwise vector entry while its
LR-04 metric weight is exactly zero.

## Endpoint and reciprocal gauge convention

The caller supplies the arbitrary-`k` endpoint eigenvectors returned by the
production reciprocal service, with the second state evaluated at the exact
`k+q` endpoint after reciprocal folding. The vertex consumes those coefficients
verbatim and performs no manual rephasing.

The live Fourier convention is `exp(+i 2*pi*k.R)`. For a reciprocal shift
`G`, the endpoint site gauge is

\[
 D_a(G)=\exp[-i2\pi G\cdot\tau_a],
 \qquad c_a(k+G)=D_a(G)c_a(k).
\]

The complete local transition vector consequently transforms by the same
second-endpoint site phase when only the right endpoint is folded:

\[
 T_{aLMi}^{\rm folded}=D_a(G)T_{aLMi}^{\rm raw}.
\]

This is an endpoint/eigenvector contract. It is not a susceptibility phase,
and the evaluator does not introduce a separate finite-`q` phase.

For a Hermitian explicit `Gamma` and equal radial endpoints, the evaluated
vectors obey

\[
 T_{mn;aLM}(r_i)=(-1)^M
 \left[T_{nm;a,L,-M}(r_i)\right]^*.
\]

Individual eigenvectors in an exactly degenerate subspace have no independent
physical meaning. The primitive may evaluate a selected representative;
quantities summed over the degenerate subspace must be invariant under its
unitary rotations. No degeneracy rotation or gauge fixing is hidden in this
API.

## Tests and evidence

`UnitLrPauliTransitionVertex` is a focused standalone Fortran test registered
with CTest. It covers:

- independent Gauss--Legendre/uniform-phi sphere quadrature against the
  production Gaunt evaluator for `s-s`, off-diagonal `p-p`, `p-d`, off-diagonal
  `d-d`, charge, `sigma-z`, and `sigma-plus`;
- the Hermitian endpoint relation;
- a two-site nontrivial-basis-position folded-endpoint test on the complete
  response vector;
- occupied diagonal charge and `sigma-z` closure against independent
  large-component LMTO radial products on production RSEQSR/PHDFSR radial
  arrays;
- hard rejection of SOC, generalized overlap, noncollinear, `spdf`, and
  additive-operator capability requests.

The focused commands are:

```text
cmake --build build --target UnitLrPauliTransitionVertex -j2
ctest --test-dir build --output-on-failure -R '^UnitLrPauliTransitionVertex$'
ctest --test-dir build --output-on-failure -R '^UnitLrPauliTransitionVertexReject'
```

These are algebraic and independent numerical-oracle claims. They do not
promote the scalar-relativistic-to-Pauli discrepancy documented in
`docs/LR_SR_PAULI_NUMERICAL_CLOSURE.md` into exact-SR response correctness,
and they do not claim Fe/Ni convergence, Goldstone correctness, or literature
agreement.

## Scope boundary

This commit adds only the one-electron Pauli transition vertex. It does not add
occupation differences, frequency denominators, `chiKS`, a reciprocal-GF bubble,
an XC kernel, Ward enforcement, a Dyson solver, a Goldstone correction, or a
site-only fallback.
