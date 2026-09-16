# DRESP-03T — Local-Torque and Force-Theorem Hessian Repair

Status: **PASS for the complete orthogonal-H fixture; contact term required.**
The native LKAG implementation remains the read-only Rung-0 reference.  No
file under `source/exchange.f90` was changed.

## Implemented representation

`source/lmto_magnetic_tangent.f90` now exposes the first and mixed derivatives
of the live `lmto_bond_value` algebra.  The derivatives include the onsite
`c1` term and all bond dot/cross terms.  `source/lr_kl_hessian.f90` assembles
the complete site-major k-space operator

\[
 H(k)=B(k)-Q(k)B(k)+E_\nu,
 \qquad Q(k)=\sum_R ee(R)\,obar_j e^{i2\pi k\cdot R},
\]

when HOH is enabled.  Its torque and mixed derivative include the Fourier
phase and the full product rule for `Q*B`:

\[
 T_i=B_i-Q_iB-QB_i,
\]

\[
 C_{ij}=B_{ij}-Q_{ij}B-Q_iB_j-Q_jB_i-QB_{ij}+E_{\nu,ij}.
\]

For distinct sites, the onsite `c1*m` and `e_nu` terms have zero mixed
derivative, while the mixed bond dot/cross terms remain nonzero.

## Hessian convention

For `G=(E-H)^{-1}` and the grand-potential convention used by the force
theorem, the exact second derivative is

\[
 \Omega_{ij}=-\frac{1}{\pi}\operatorname{Im}\int^{E_F}dE\,
 \operatorname{Tr}\left[T_iGT_jG+C_{ij}G\right].
\]

The implementation reports the two pieces independently.  The corresponding
gapped finite-spectrum form is the ordinary eigenvector-response term plus
`sum_(n occupied) <n|C_ij|n>`.  This is the explicit decision required by
DRESP-03T: the response is not torque-torque alone; the mixed contact term is
part of a representation-complete Hessian whenever `C_ij` is nonzero.

## Tests

`UnitDresp03tLmtoHessian` uses a two-site, nontrivial complex bond fixture at
several k points.  Its channels are generated from one underlying set of
`c`, `enu`, `srdel`, and `qpar` inputs by the same `predls` transformation as
production, then converted to the live `cx/wx/obx/e_nu` channels.  It checks:

* full complex-matrix first derivatives against central differences;
* `C_ij` against four independent Hamiltonian evaluations;
* direct gapped grand-potential finite differences against TT-only and
  TT-plus-contact spectral Hessians; and
* nonzero offsite torque and nonzero contact response.

The current run gives first-derivative error `3.2e-12`, mixed-derivative error
`3.4e-10`, and complete grand-potential Hessian error `6.9e-9`.

The production-input adapter `lmto_fixture_from_hamiltonian` is read-only and
uses the same directed lattice blocks and fractional neighbor vectors as the
reciprocal Fourier path.  The old
`build_ham_only_exchange_vertex` name is retained only as a compatibility
wrapper around the explicitly generic `build_supplied_ham_exchange_vertex`;
it is not a physical LMTO torque constructor.

## Remaining gate

This rung certifies the orthogonal-H torque/Hessian algebra and its independent
energy oracle.  A native canonical LKAG pair-by-pair Fe gate still has to be
run on the same production potential, mesh, Fermi level, orbital subset, and
approximation order.  Any disagreement belongs first to the new finite-H
bridge; the native `source/exchange.f90` route must remain unchanged.
