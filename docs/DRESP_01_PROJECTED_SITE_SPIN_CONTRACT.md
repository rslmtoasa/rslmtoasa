# DRESP-01 — Projected Site-Spin Operator and Matching Moment Contract

Status: implemented and certified on `fable_v4` at the DRESP-00 starting
revision.  This document is the boundary contract for the projected site-spin
operator.  It does not authorize a susceptibility, exchange kernel, or Dyson
implementation.

## Scope

DRESP-01 owns one-electron projected site-spin vertices and the matching
ground-state moment audit for the accepted baseline:

- orthogonal, `ham_only`, second-order/HOH LMTO augmentation;
- collinear two-spin Pauli blocks, with the accepted large-component endpoint
  convention;
- complete `spd` orbital basis (`l=0,1,2`) and the complete product response
  cutoff (`response_lmax = 4`);
- one explicit response channel and the existing logarithmic radial mesh.

The implementation is `lr_projected_site_spin_mod` in
`source/lr_projected_site_spin.f90`.

## Representation and selectors

The response/product representation is unchanged.  The product basis keeps
all Gaunt-allowed response harmonics `(L,M)`.  The two DRESP selectors are
one-electron orbital selectors:

| selector | retained orbital content |
| --- | --- |
| `d` | `l=2` only |
| `spd` | `l=0,1,2` |

Neither selector truncates response `L`.  `spdf` is unsupported and fails
closed at initialization; `projected_selector_supported('spdf')` returns
false.

The selector is applied to candidate records before the candidate-to-compact
map is contracted.  This preserves the existing retained product modes and
does not introduce a dense response vertex.

## Site integration

For a product block, the live product map reconstructs the harmonic-projected
radial coefficient as

\[
x_{LM}(r_i)=\frac{1}{\sqrt{W_i}}\sum_a U_{ia}z_{LM,a}.
\]

The full-sphere scalar integral is carried only by `L=0,M=0` because the code
harmonic is normalized with

\[
Y_{00}=\frac{1}{\sqrt{4\pi}},\qquad
\int Y_{00}\,d\Omega=4\pi Y_{00}.
\]

The implementation derives the scalar factor as
`4*pi*real(response_harmonic(0,0,0,0))`.  For each site the compact
functional is therefore

\[
p_a^* = (4\pi Y_{00})\sum_i \sqrt{W_i}\,U_{ia},
\qquad I_{\mathrm{site}}=p^H z,
\]

with the code’s complex conjugation convention stored in `p`.  Functional
entries for every `L>0`, `M`, and every other site are exactly zero.

The origin has zero response-space volume weight.  No epsilon-radius or
origin phase is introduced.

## Operator contract

The exact DRESP-02-facing API is intentionally small:

```fortran
type(projected_site_spin_contract) :: contract
call contract%initialize(space, radial_bases, 'd') ! or 'spd'

call contract%selected_transition_coordinates(product_plus, left, right, z)
call contract%site_integration_functional(product_plus, p)
call contract%transition_amplitudes(product_plus, left, right, t_plus)

call contract%direct_operator_matrix(radial_bases, energy, projected_operator_plus, v_plus)
call contract%direct_transition_amplitudes(radial_bases, left, right, &
     projected_operator_plus, t_direct)
```

`product_plus` and `product_minus` are existing
`lmto_product_response_basis` objects initialized with the same `space` and
radial basis.  The `left` and `right` arguments are already prepared
`pauli_endpoint_state` objects.  Thus their energy, k-point, folding, and
gauge provenance is consumed, not reconstructed by DRESP-01.

The direct matrix is an independent oracle for the same accepted large-
component endpoint product.  It is energy-parameterized and is not a static
response operator.  The supported operator kinds are
`projected_operator_plus`, `projected_operator_minus`, and
`projected_operator_z`.  The minus matrix is the adjoint of plus in the
accepted real radial baseline; z is Hermitian.

## q=0 and finite q

At `q=0`, the endpoint pair uses the same accepted k-space state and the
transition amplitude is the site integral above.  At finite q, DRESP-01
accepts the already prepared right endpoint from
`lr_q_endpoint_from_reciprocal` (or an equivalent certified endpoint
service).  It does not add a site phase, fold a k-point, or rephase an
eigenvector.  Any multi-site `q`/tau phase is part of the endpoint/basis
contract upstream.

The unit test uses a mesh-compatible nonzero endpoint gauge to prove that the
same no-extra-phase rule holds for finite q.

## Matching projected moments

The contract provides two independent valence-only `sigma_z` moment routes:

```fortran
call contract%moment_from_density(eigenvalues, eigenvectors, k_weights, ef, temperature, &
     ground_states, moment_density)
call contract%moment_from_operator(eigenvalues, eigenvectors, k_weights, ef, temperature, &
     ground_states, moment_operator)
```

The eigensystem arrays are the fields of an accepted immutable reciprocal
snapshot: `eigenvalues(nband,nk)`,
`eigenvectors(nbasis,nband,nk)`, and normalized or raw nonnegative
`k_weights(nk)`.  DRESP normalizes the weight sum once.  Occupations use the
same Ry/K Fermi-Dirac convention as the reciprocal state (`kT=max(T*kB,
1e-10 Ry)`).

The density route first accumulates the selected radial profiles

\[
\rho_{l\sigma}(r_i)=\sum_{k n}w_k f_{nk}
 |c_{nk,lm\sigma}|^2 U_{l\sigma}(r_i,E_{nk})^2,
\]

then applies the accepted logarithmic radial integral to
`rho_up-rho_down`.  The operator route builds the selected site-local
`sigma_z` matrix from independent radial norms and contracts it with each
one-particle state.  Agreement is the normalization and selector closure
required before a dynamical consumer may use the vertices.

The code moment convention is positive up-minus-down spin number.  With
`g=2`, the reported scalar has the same numerical value in `mu_B`; no extra
factor is inserted into the Pauli response operator.

Frozen core is not silently added.  The persisted `radial_ground_state`
contract keeps core Pauli density only as total per-spin radial arrays, with
no l resolution.  DRESP-01 therefore uses
`valence-only; frozen core excluded` for both `d` and `spd` projected moments.
`core_spin_number` may report the available total core spin as context, but it
is never merged into a selected orbital moment.

## Certification

The focused unit `UnitLrProjectedSiteSpin` covers:

- `d` versus `spd` one-electron selection;
- exact `L=0,M=0` site functional and zero cross-site leakage;
- q=0 and finite-q prepared endpoints;
- direct radial/operator versus compact product contraction;
- plus/minus adjoint and z Hermitian matrix relations;
- independent density-profile and operator-expectation moments;
- fail-closed `spdf` behavior.

On the finite analytic spd radial fixture, the largest product/direct
transition residual is `2.2e-19`; the density/operator moment residual is
`3.5e-18`.  The bcc-Fe material gate
`tests/validation/val_dresp01_projected_site_spin.py` uses the accepted
artifacts under `example/susceptibility/bccFe` and reports:

| quantity | value (`mu_B`) |
| --- | ---: |
| d projected moment | 2.096590 |
| spd projected moment | 2.003489 |
| accepted total moment | 2.003488859 |

The spd/total residual is `1.4e-7` in that artifact, below the `2e-5`
material-gate tolerance set by six-decimal band-moment output rounding.

No projected susceptibility or dynamical response is certified by DRESP-01.
