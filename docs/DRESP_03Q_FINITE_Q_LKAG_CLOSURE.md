# DRESP-03Q — Finite-q Torque-Hessian and Native LKAG Closure

Status: **BLOCKED — COMMON-STATE LKAG q ORACLE NOT CERTIFIED**.

The finite-q operator algebra and the live reciprocal-H adapter are implemented
and independently tested.  The old Fe comparison is retained below as a
historical two-shell diagnostic only: its `jij.out` source requested
`njij = 2`, so its Fourier curve is a truncated model and is not a
representation-equivalence gate.  The new same-state full-range comparison is
owned by the `exchange_q` workflow, but its q-space native-LKAG oracle still
requires independent validation and metallic convergence evidence.  No
rescaling or fitted tolerance is used.  DRESP-04 remains forbidden until this
status changes.

`source/exchange.f90` was not modified.

## Scope and implementation

The new finite-q API is in `source/lr_kl_hessian.f90`:

- `assemble_lmto_finite_q_torque` and its batched form return the complete
  `k -> k+q` derivative, including every directed bond and both endpoint
  phases;
- `assemble_lmto_finite_q_mixed_derivative` returns the `q,-q` contact
  derivative, including direct mixed dot/cross terms, same-site curvature,
  phase products, and the complete HOH product rule;
- `force_theorem_finite_q_hessian_from_eigenbasis` and its BZ batch wrapper
  return complex sublattice matrices for the TT-only, contact-only, and
  complete grand-potential Hessians; and
- `lmto_fixture_adapter_residual` audits the fixture against the production
  reciprocal assembler without writing production Hamiltonian state.

For a texture convention

\[
 \theta_{aR}=\theta_a(q)\exp(+i2\pi q\cdot(R+\tau_a)),
\]

the first vertex is the `k -> k+q` block.  Its endpoint phases are retained
separately for the source and target of each directed bond.  In the production
Bloch basis the common absolute endpoint phase cancels: the source factor is
one and the target factor is `exp(+i2*pi*q.D)`, with
`D=R+tau_target-tau_source`, while the underlying real-space texture remains
`exp(+i2*pi*q.(R+tau))`.  The reverse vertex is evaluated at `k+q` with `-q`.
A reciprocal-vector shift is not silently folded inside the vertex routine; it
remains the corresponding Bloch gauge phase.

The fixture has an explicit `include_enu` convention.  Synthetic historical
fixtures retain the `E_nu` term, while a fixture extracted from production
uses the reciprocal `ham_only` convention: first order is `H=B`, and the
second-order HOH path is `H=B-QB+E_nu`.  This distinction is required because
the reciprocal first-order production assembler deliberately omits `E_nu`.

## Independent finite-q oracle

`tests/unit/test_dresp03q_finite_q.f90` uses a one-orbital period-two
supercell.  Its finite differences construct the rotated Hamiltonian directly
from `lmto_bond_value`; no torque, mixed derivative, or Hessian helper is used
inside the oracle.  It checks complete versus TT-only, contact necessity, q=0
reduction to DRESP-03T, q reversal, reciprocal-period gauge covariance, and
the finite-q Hessian reversal relation.  A separate two-sublattice fixture
checks the finite-q vertex adjoint relation and full complex Hessian
Hermiticity.

Representative result from the current build:

```text
complete-vs-FD              1.10e-08
TT-only error               5.00e-01
contact                     5.00e-01
q=0 reduction               0
q/-q Hessian covariance     0
reciprocal-period covariance 0
```

The nonzero contact term is required for the complete finite-q curvature.

## Production reciprocal-H adapter

`tests/unit/test_dresp03q_production_adapter.f90` builds bcc Fe through the
live control/lattice/charge/Hamiltonian/reciprocal stack.  The gate is clean:
scalar-relativistic, no SOC, no Hubbard, no CCOR, no constraining field,
orthogonal `ham_only`, full `s,p,d` orbital space, and the production HOH
ordering is checked.  Five off-mesh k points are compared.

```text
max |H_fixture(k)-H_reciprocal(k)|       5.55e-17
real-space adapter residual              0
```

The spherical-orbital conversion is applied to the copied live bond channels
because production `hamiltonian_build` converts each of the four spin blocks
before storing `ee`.  This is a representation correction, not a numerical
fit.

## Hessian derivation and prefactors

Use the resolvent convention \(G(z)=(z-H)^{-1}\).  For a fixed reference
potential and \(\Omega=\operatorname{Tr} f(H-\mu)\), differentiating the
trace-log gives the finite-q mixed coefficient

\[
 \Omega_{ab}(q,-q)=
 -\frac{1}{\pi}\operatorname{Im}\int^{E_F}dE\,
 \operatorname{Tr}\left[
 T_a(q)G_{k+q}T_b(-q)G_k + C_{ab}(q,-q)G_k
 \right].
\]

The first term is the two-vertex response.  The second is the direct contact
term from the second derivative of the Hamiltonian.  In a gapped finite
spectrum, the same result is evaluated independently as

\[
 \sum_{n\,\mathrm{occ}}\langle n|C_{ab}|n\rangle
 +\sum_{n\,\mathrm{occ},m}
 \frac{
 \langle n|T_a(q)|m\rangle\langle m|T_b(-q)|n\rangle
 +\langle n|T_b(q)|m\rangle\langle m|T_a(-q)|n\rangle
 }{\epsilon_n-\epsilon_{m,k+q}}.
\]

This fixes the sign and the absence of an extra factor of two directly from
the resolvent identity.  Degeneracies wholly on one side of the zero-
temperature occupation edge are skipped as a basis-invariant occupied/empty
subspace contribution; a Fermi-edge degeneracy is rejected as non-analytic.

For the native scalar exchange convention, derive the Fourier comparison from
the classical pair Hamiltonian

\[
 H_{\mathrm{spin}}=-\frac12\sum_{R}J(R)\,e_0\cdot e_R,
 \qquad
 J(q)=\sum_R J(R)\exp(iq\cdot R).
\]

The mixed transverse coefficient is therefore `J(0)-J(q)` per primitive
cell.  No moment divisor, gyromagnetic factor, or TDDFT kernel belongs in this
comparison.

## Native LKAG Fourier audit

The native kernel is read directly from `source/exchange.f90`:

- `dGdG_Jnc` forms the charge-minus-three-Pauli-component product;
- `imtrace9` supplies the integrated scalar kernel;
- `simpson_f` integrates to the native Fermi level; and
- `calculate_exchange` writes `J_internal*1000/(4*pi)` in mRy.

The checked input contains only the two bcc shell representatives

```text
R = (-1/2,-1/2,-1/2)   J = 0.739154 mRy
R = ( 0,  0, -1)       J = 0.467405 mRy
```

`tests/unit/test_dresp03q_native_lkag.f90` expands these over the exact eight
first-shell and six second-shell vectors.  This matters for the selected
finite-q line: the second-shell axis vectors do not all have the same phase.
The finite-q API uses primitive reciprocal coordinates, whereas `jij.out`
prints conventional Cartesian `R` in units of `a`; for bcc the conversion is

```text
q_cart = (q2+q3, q1+q3, q1+q2),
phase  = exp(i*2*pi*q_cart dot R).
```

Expanding those representatives over their cubic orbits gives a complete
two-shell Fourier transform, not the infinite/full native LKAG interaction.
The required 8^3 full mesh was exercised at q values 1/8, 1/4, 3/8, and 1/2
along the bcc primitive reciprocal direction.  The historical diagnostic
reported:

```text
q=1/8   finite-H  4.396331e-03   native J(0)-J(q)  1.413572e-03
q=1/4   finite-H -1.103059e-02   native J(0)-J(q)  4.826236e-03
q=3/8   finite-H  5.889109e-03   native J(0)-J(q)  8.238900e-03
q=1/2   finite-H  7.573012e-03   native J(0)-J(q)  9.652472e-03
```

An independent reciprocal Lehmann diagnostic with `kspace_ham_order='first'`
and `green_eta=10^-3` also changes the native two-shell values to approximately
`0.326186` and `0.451905` mRy.  This is numerical-state evidence about the
historical truncated diagnostic, not proof of a native/finite-H representation
mismatch.  It must not be used as the DRESP-03Q closure gate.

## Same-state full-range LKAG-q comparison

`post_processing='exchange_q'` consumes the accepted reciprocal SCF snapshot.
The finite-H fixture is copied before `predls()` is used for the native vertex;
the production-H adapter residual is checked before and after that mutation.
The q-space reference evaluates the exact native `dGdG_Jnc` contraction in
reciprocal space, with the native spin ordering, `imtrace9` imaginary trace,
native `simpson_f` energy mesh, and the native `1/(4*pi)` output conversion.
It is explicitly called an LKAG q-space oracle, not a replacement for
`source/exchange.f90`.

The production gate is not closed until all of the following are recorded:

1. real-space/Fourier identity on an independently evaluated finite/gapped
   fixture;
2. native-pair spot checks against unchanged `exchange.f90`;
3. energy-mesh and regulator convergence for metallic Fe;
4. dense same-state curves for finite-H and LKAG-q with no fitted scale.

The allowed final verdicts are `PASS — NATIVE LKAG / FINITE-H q BRIDGE
CLOSED`, `PASS — CONTACT TERM REQUIRED FOR CLOSURE`, `BLOCKED — COMMON-STATE
LKAG q ORACLE NOT CERTIFIED`, `BLOCKED — METALLIC INTEGRATION NOT CONVERGED`,
and `BLOCKED — TRUE NATIVE/FH REPRESENTATION MISMATCH`.  Only the last one
authorizes reopening the deeper representation mapping.

## User-facing exchange_q workflow

Use an accepted bulk k-space SCF state:

```fortran
&calculation
  pre_processing  = 'bravais'
  post_processing = 'exchange_q'
/
&self
  use_kspace = .true.
/
&exchange_q
  q_coordinates = 'direct'
  q_file         = 'qpath.dat'
  n_q_points     = 0
  output_file    = 'exchange_q.dat'
  write_components = .true.
  native_crosscheck = .false.
  native_green_eta = 1.0e-3
  native_energy_points = 0
  rotation_axis = 1.0, 0.0, 0.0
/
```

For short paths, replace `q_file` with `n_q_points` and columns
`q_list(:,1)`, `q_list(:,2)`, etc.  A q file follows the frozen-magnon syntax:

```text
NQ
q1 q2 q3
...
```

`direct` values are reciprocal-lattice coordinates.  `cartesian` values are
Cartesian units of `2*pi/alat`, matching the established frozen-magnon
convention.  The electronic k mesh and magnetic q path are independent; exact
`k+q` endpoint diagonalization is used, so q values are not silently rounded.
When the optional native cross-check is enabled, `native_green_eta` controls
the retarded regulator and positive `native_energy_points` rebuilds the native
Simpson energy mesh for a reproducible convergence sweep.

The primary output is `DeltaJ(q)=J(Gamma)-J(q)`, with TT, contact, total,
`mRy`, and nonzero-q `DeltaJ/q^2` columns.  If `native_crosscheck=.true.`,
the native LKAG-q columns and their residual are appended.  No absolute
finite-H `J(q)` is manufactured by adding a native Gamma constant.

## Re-opening condition

Re-open DRESP-03Q only after a common-state native/finite-H comparison has an
independently controlled metallic integration, identical Hamiltonian order and
orbital representation, and agreement for the full required q mesh.  Until
then, DRESP-04 must not be started.
