# RSGF-00 — Native-GF endpoint augmentation

## Scope and claim

This document certifies the endpoint seam between the native coefficient-space
Green function and the Pauli/no-SOC spatial one-electron representation.  It
does not certify a susceptibility energy integral, a GF bubble, Dyson response,
an XC kernel, or a material driver.

The implementation is `source/lr_gf_endpoint_augmentation.f90`, in
`lr_gf_endpoint_augmentation_mod`.  It is response-owned and has no dependency
on `green_mod`.  This preserves the distinction between physical augmentation
and `green%auxiliary_gij`: the latter performs a screened/auxiliary
representation transformation and remains unchanged.

The PASS claim is only:

> Native coefficient-space GF blocks can be promoted to the same certified
> Pauli/no-SOC augmented one-electron spatial representation used by the
> reciprocal response, within the supported orthogonal second-order LMTO
> baseline.

## Supported capability

The public adapter accepts only this tuple and fails closed otherwise:

| capability | supported value |
|---|---|
| reciprocal representation | `ham_only` |
| Hamiltonian order | `second` / HOH |
| overlap | orthogonal |
| magnetic state | collinear |
| SOC | absent |
| Hubbard/additive operators | absent |
| orbital basis | `sp` or `spd` |
| endpoint radial response | Pauli large component |

The radial basis and `E_nu_work` values come from
`lmto_radial_augmentation_mod`.  Angular factors use
`response_angular_basis_mod`; no second spherical-harmonic convention is
introduced.  The factored result stores the large-component radial maps and
can evaluate a 2-by-2 Pauli spin Green function at two radial/angular points.

## Formal derivation

For the orthogonal LMTO basis, let

\[
  H=E_\nu+h^\gamma, \qquad G(z)=[zI-H]^{-1},
\]

where `E_nu` is the live working linearization energy (`enu_work`), not a raw
atomic energy selected independently of the production Hamiltonian.  The
physical Pauli endpoint map is

\[
  \mathcal A=\Phi+\dot\Phi h^\gamma.
\]

Therefore the spatial Green function between endpoints `a` and `b` is

\[
  G^P_{ab}(z)=\mathcal A_a G_{ab}(z)\mathcal A_b^\dagger.
\]

Set `D(z)=zI-E_nu`.  Since `h^gamma=H-E_nu`, the resolvent identities are

\[
  h^\gamma G=DG-I,\qquad Gh^\gamma=GD-I,
\]

and

\[
  h^\gamma G h^\gamma=DGD-D-h^\gamma.
\]

For a site/orbital block these are

\[
 (hG)_{ab}=D_aG_{ab}-\delta_{ab}I,
\]

\[
 (Gh)_{ab}=G_{ab}D_b-\delta_{ab}I,
\]

\[
 (hGh)_{ab}=D_aG_{ab}D_b-\delta_{ab}D_a-h^\gamma_{ab}.
\]

The production adapter stores all four coefficient-space terms explicitly and
evaluates

\[
\begin{aligned}
G^P_{ab}={}&\Phi_aG_{ab}\Phi_b^\dagger
 +\dot\Phi_a(hG)_{ab}\Phi_b^\dagger\\
 &+\Phi_a(Gh)_{ab}\dot\Phi_b^\dagger
 +\dot\Phi_a(hGh)_{ab}\dot\Phi_b^\dagger.
\end{aligned}
\]

The two `-I` terms are present only for an onsite block.  The `-h^gamma_ab`
term is present for both onsite and offsite blocks; for an offsite block it is
the only nonzero contact-like term in the double-sided identity.  In
particular, the adapter never replaces either `h^gamma` by a bare
`z-E_nu` shortcut.

## Provenance of `h^gamma_ab`

The API accepts a directly available effective block, or a
`lr_gf_effective_hamiltonian_provider`.  The provider contract is deliberately
an action contract: it applies the production `H_eff` action to a localized
identity seed at source site `b` and returns the destination block at site `a`.

For `a != b`, this returned block is `h^gamma_ab`.  For `a=b`, the adapter
subtracts the diagonal `E_nu_work` block derived from the endpoint radial
basis.  This is the required `H_eff`-action provenance and prevents an HOH
calculation from silently using a duplicated raw first-hop formula.  The
included dense provider is a small fixture/provider implementation; a native
recursion caller can wrap its already-completed effective-Hamiltonian action
behind the same abstract interface.

The reverse pair API accepts `G_ba` and `h^gamma_ba` explicitly when both
directions are needed.  The forward block does not infer a reverse block or
silently conjugate a non-Hermitian native GF.

## Radial and angular representation

The stored radial arrays are the certified large-component LMTO numerators
`phi_large` and `phidot_large`, repeated in the production orbital/spin
ordering.  The point evaluator forms

\[
  \frac{U_l(r)}{r}Y^{\rm code}_{lm}(\Omega),
\]

with `response_harmonic`, where `Y^code` is the existing response angular
convention.  The output is factored rather than a materialized four-coordinate
array: coefficient branches and radial maps are retained, and
`evaluate_point` returns the same pointwise Pauli spinor object needed by a
later radial/angular bubble.  No scalar-relativistic lower-component angular
response is claimed here; this is explicitly the Pauli/no-SOC projection.

## Tests and controls

`tests/unit/test_lr_gf_endpoint_augmentation.f90` constructs a finite,
two-site orthogonal LMTO-like model with nontrivial diagonal `E_nu`, onsite
terms, nonzero offsite hopping, and representative `phi`/`phidot` radial
arrays.  It checks four complex energies and both onsite and offsite endpoint
blocks.

The independent checks are:

1. Dense reference: form the full matrix `H=E_nu+h^gamma`, invert `zI-H`,
   construct the global endpoint maps explicitly as `A=Phi+Phidot*hgamma`,
   and compare `A G A^dagger` with the factored adapter.
2. Effective-action provenance: obtain `hgamma_ab` through the dense
   `H_eff` action provider and repeat the offsite comparison.
3. Spectral reference: diagonalize `H`, form `Psi_n^P=A c_n`, and compare the
   Lehmann sum `sum_n Psi_n^P Psi_n^{P dagger}/(z-epsilon_n)`.
4. Negative controls deliberately omit, one at a time, onsite `-I` contacts,
   offsite `-hgamma_ab` in `hGh`, and all `phidot` augmentation.  Each is
   required to differ from the dense reference by a nonzero floor.

The accepted block terminator or native recursion accuracy is not part of this
unit oracle.  A future native-GF integration test should report coefficient-GF
error and endpoint-augmentation error separately.

## Claim level

This is a foundational endpoint implementation and numerical algebra
certification for the stated finite-model and capability tuple.  It does not
promote `green%auxiliary_gij` into a radial-response routine, and it does not
validate the GF bubble, TD-DFT response, or any unsupported SOC, overlap,
noncollinear, Hubbard, additive-operator, or `spdf` extension.
