# DRESP-03 — Katsnelson–Lichtenstein Static Bridge to LKAG

Status: **audit complete; overall result BLOCKED — FINITE-HAMILTONIAN /
LMTO-LKAG REPRESENTATION MISMATCH**.

The finite orthogonal `ham_only` bridge is closed as an explicit operator
vertex.  The scalar compression is conditional and is rejected when the
supplied operator is orbital dependent.  A common-state numerical Fe
comparison to the mature native LKAG path is not certified because the two
paths currently expose different exchange vertices.  No interaction kernel,
Dyson dressing, Mills/Jülich U, Goldstone correction, rescaling, or magnon
interpretation is introduced here.

## 1. Native LKAG formula/code map

The mature combined exchange path is
`source/exchange.f90:1143-1326`, called here the Rung-0 reference.  For a
requested directed pair `(i,j)`, `calculate_exchange` evaluates, at every
energy point,

\[
 K_J(E) = \operatorname{ImTr}\left[
 \Delta_i(E)G^0_{ij}(E)\Delta_j(E)G^0_{ji}(E)
 - \sum_{a=x,y,z}\Delta_i(E)G^a_{ij}(E)\Delta_j(E)G^a_{ji}(E)
 \right].
\]

The exact matrix products are in `dGdG_Jnc`,
`source/exchange.f90:640-666`.  The forward blocks are
`Ginmag/Gix/Giy/Giz`; the directed reverse blocks are
`Gjnmag/Gjx/Gjy/Gjz`.  The trace is `imtrace9` at
`source/exchange.f90:1204-1210`.  Simpson integration to `en%fermi` is at
`source/exchange.f90:1227-1229`, followed by the writer convention

\[
 J_{ij}^{\rm file}[\mathrm{mRy}] =
 \frac{10^3}{4\pi}J_{ij}^{\rm internal}[\mathrm{Ry}].
\]

The seven-column `jij.out` record is `type_i type_j R_x R_y R_z J distance`.
It contains one selected pair representative, not a shell multiplicity.  The
repository's independent Fourier helper therefore requires explicit shell
multiplicities; it uses the established bcc convention
`q` in Cartesian `2*pi/alat` coordinates and `R` in `alat` units.

The canonical path has no separate on-site/additive branch: it evaluates the
configured pair list with the same `dGdG_Jnc` kernel.  The distinct auxiliary
route, `calculate_jij_auxgreen` at `source/exchange.f90:188-352`, does have an
explicit `i == j` branch (`source/exchange.f90:309-316`) containing the local
term and its `-1` factor.  It forms `DeltaP=P_up-P_down`, uses auxiliary
Green functions, and is not a production caller in the current workflow.

The x/y/z subtraction above is the combined isotropic exchange kernel.  DMI
and anisotropic terms are accumulated separately in the same production
routine; the two-index SO/FO decomposition is not substituted for the
combined observable.  There is no Goldstone or other response correction in
this native LKAG calculation.

## 2. Exchange-vertex taxonomy

The live source contains several objects that must not be conflated.

| object | source contract | energy dependence | status in DRESP-03 |
|---|---|---:|---|
| `symbolic_atom%d_matrix(E)` | canonical `dGdG_Jnc` vertex | yes | native LKAG reference |
| `DeltaP(E)=P_up(E)-P_down(E)` | auxiliary/path-operator route | yes | separate auxiliary candidate |
| `H_up-H_down` | `exchange.f90:1368-1371` `ee` blocks | no | finite bridge |
| DRESP-01 `V_i^+` | site-integrated Pauli/radial vertex | left/right endpoint energies | projected vertex, not LKAG `Delta` |
| radial `B_xc(r)` | physical local field primitive | radial/local | no certified coefficient-to-LKAG map exposed |

The canonical `d_matrix` implementation is
`source/symbolic_atom.f90:313-341`.  For each `l,m` it sets

\[
 d_l(E)=\frac{c_d w_u^2-c_u w_d^2+(w_d^2-w_u^2)E}{w_u w_d},
\]

with the Madelung-shifted `c` values and spin-dependent `dele` values from
the symbolic atom.  It is diagonal in the nine `s,p,d` orbital channels, but
it is not a fixed scalar site parameter.  The auxiliary `P` object is built
as `(E-(c+vmad))/dele^2` in `source/symbolic_atom.f90:477-506`.

The finite bridge therefore declares its vertex provenance explicitly as
`finite ham_only local H_up-H_down spin-flip vertex`.  The declaration is not
an assertion that this matrix equals `d_matrix(E)`.

## 3. Finite-Hamiltonian K/L derivation

The finite bridge uses the same immutable collinear, orthogonal, second-order
`ham_only` electronic snapshots used by DRESP-02.  For a site-local supplied
coefficient-space operator \(V_i^+\), define

\[
 X_i^{nm}(k,q)=
 \langle n,k\vert V_i^+\vert m,k+q\rangle.
\]

The DRESP-02 static circular susceptibility convention is

\[
 \chi^{0,+}_{ij}(q,0)=
 \frac{2}{\sum_k w_k}\sum_{k,n,m}w_k
 \frac{f_{nk}-f_{m,k+q}}
 {\epsilon_{nk}-\epsilon_{m,k+q}+i\eta}
 X_i^{\sigma,+}X_j^{\sigma,+*},
\]

where \(X^{\sigma,+}\) is the unit spin-flip transition amplitude.  The
finite source-convention K/L contraction implemented in
`source/lr_kl_static_bridge.f90` is

\[
 J^{H}_{ij}(q;\eta)=
 -\frac14\frac{1}{\sum_k w_k}\sum_{k,n,m}w_k
 \frac{f_{nk}-f_{m,k+q}}
 {\epsilon_{nk}-\epsilon_{m,k+q}+i\eta}
 X_i^{nm}X_j^{nm*}.
\]

The `-1/4` is the locked `kl_static_prefactor`; it is not a fit.  The sign
comes from reducing the source `J=Im Tr[...]` real-axis LKAG convention to
the spectral denominator.  Because DRESP-02 carries the factor two in its
unit-spin-flip \(\chi^0\), the exact scalar special case is

\[
 J^{H}_{ij}=-\frac18\,\Delta_i\,
 \chi^{0,+}_{ij}\,\Delta_j,
\]

implemented only as the algebraic helper
`contract_scalar_site_chi0`.  It accepts a predeclared scalar and never
infers one from the response.  The bridge result intentionally remains in
the internal source normalization; it does not apply the native
`10^3/(4*pi)` display conversion.

The `eta` in this finite spectral expression is a static pole regulator.  It
is distinct from DRESP-02's real-axis `integration_eta`; the latter is not
reopened or folded into this bridge.

## 4. Mapping to native LMTO LKAG

There is no certified identity in the current tree between the finite
operator \(H_\uparrow-H_\downarrow\) and canonical LMTO
`symbolic_atom%d_matrix(E)`.  The finite operator is a fixed coefficient
Hamiltonian difference.  The canonical native vertex is an energy-dependent
potential-function difference, and the auxiliary route instead uses
`DeltaP(E)` together with endpoint/path-operator transformations.

This matters in two ways:

1. A fixed \(H_\uparrow-H_\downarrow\) cannot reproduce the canonical
   `d_matrix(E)` over the LKAG energy integration without a derived
   energy-dependent representation transformation.
2. DRESP-01 site vertices carry separate left/right endpoint energies through
   radial overlaps and four affine endpoint components.  A one-energy
   `DeltaP(E)` or a fixed Hamiltonian block is not, by itself, the required
   two-energy coefficient-to-site transformation.

The auxiliary code's `sqrt(dele)` endpoint scaling and its special local
branch do not close this gap.  The existing architectural audit records the
same distinction: the auxiliary route is reusable only under a narrow
representation contract, while canonical `d_matrix(E)` is the Rung-0 LKAG
object.  Consequently the finite bridge is rigorous as a finite-Hamiltonian
bridge, but its equality to the mature native LMTO result is **blocked**.

## 5. Scalar-vs-operator vertex analysis

For a scalar compression to be exact in the selected site-blocked basis, the
supplied plus-channel operator must equal

\[
 V_i^+ = \Delta_i I_{\mathrm{selected\ orbitals}}
\]

with no selected-orbital dependence or selected-orbital off-diagonal terms.
The checker `assess_scalar_exchange_vertex` compares against a supplied
scalar and reports the residual; it does not construct a best-fit scalar.

The finite fixture gives:

| supplied operator | maximum absolute residual | result |
|---|---:|---|
| exact site scalar splitting | `0.0000e+00` | accepted |
| deliberately orbital-dependent/off-diagonal operator | `1.5000e-01` | rejected |

The operator-valued contraction remains valid for the second row.  This is
the required behavior: the scalar form is a theorem only after the operator
audit, not a parameterization chosen to improve agreement.

## 6. Finite fixture

`tests/unit/test_lr_kl_static_bridge.f90` constructs an independent two-site
finite Hamiltonian with two orbitals per site, eight coefficient/spin basis
states, real Hermitian intersite hopping, and a known local spin splitting.
The endpoint is a controlled site-phase version at
`q=(0.13,-0.07,0)`.  The fixture uses one k point, explicit Fermi
occupations, and `eta=0.025` for the primary contraction.

The expected value is rebuilt in a separate band-pair loop in the test; it
does not call the bridge evaluator or its helper.  The checks cover:

- plus spin-flip vertex placement and conjugated site outer-product order;
- site off-diagonal response and positive/negative-q phase ordering;
- the scalar `-1/8` compression against an independently accumulated unit
  spin-flip susceptibility;
- rejection of a non-scalar exchange operator while retaining the direct
  operator result; and
- finiteness through the static regulator ladder.

Observed primary residuals are:

| check | result |
|---|---:|
| scalar audit | `0.0000e+00` |
| operator audit residual | `1.5000e-01` |
| direct operator vs independent oracle | `6.0591e-17` |
| independent scalar compression residual | below `2e-13` test tolerance |

No expected number is derived from the implementation helper, and no
prefactor, sign, moment, on-site shift, or Goldstone correction is adjusted.

## 7. Product-space oracle

DRESP-02 already certifies the one-site complete `spd` LMTO product space as
232 compact coordinates.  The product component vertex tensor is built by
`source/lr_lmto_product_response_basis.f90:454-510`; its four components are
the left/right affine endpoint powers.  DRESP-01 contracts those components
to a site functional in
`source/lr_projected_site_spin.f90:308-338`.  The existing
`UnitLrProjectedReciprocalChi0` test independently checks the compact
Lehmann/GF oracle and the projected site matrix.

That 232-dimensional Pauli product basis is not silently treated as an
exchange-splitting basis.  A new product-space K/L contraction using it would
require the missing identification of the native `d_matrix(E)` or
`DeltaP(E)` with the four endpoint-dependent DRESP-01 operator components.
That is the representation question under audit.  For the finite fixture,
the direct coefficient-space band-pair sum is the appropriate full
vertex-weighted oracle, and it agrees at `6.0591e-17`.  The reduced scalar
site relation is independently checked only in the exact scalar fixture.

Thus the product-space evidence establishes that DRESP-02's 232 oracle is
healthy, but it does not manufacture a 232-dimensional native-LKAG vertex.

## 8. Static-limit convergence

The finite fixture was evaluated with the same operator and endpoint state at

`eta = 0.08, 0.04, 0.02, 0.01, 0.005, 0.0`.

The maximum matrix changes between successive values were:

| step | max absolute change |
|---|---:|
| `0.08 -> 0.04` | `1.9828e-02` |
| `0.04 -> 0.02` | `9.9792e-03` |
| `0.02 -> 0.01` | `4.9978e-03` |
| `0.01 -> 0.005` | `2.4999e-03` |
| `0.005 -> 0.0` | `2.5002e-03` |

All matrices, including the zero-regulator result, are finite for this
non-degenerate fixture.  The smooth ladder establishes the finite fixture's
static contraction and exposes the expected regulator dependence; it is not
a claim of material convergence for the native Fe contour.  DRESP-02's
`integration_eta` remains a separate real-axis quadrature control.

## 9. Fe \(J(R)\)

The mature ordinary bcc-Fe post-processing reference is
`tests/postproc/references/Example_exchange_bccFe/ref.json`, generated from
the production `post_processing='exchange'` path with `nsp=2`, block
recursion, no HOH, and the two configured pair records.  Its rounded rows
are:

| representative bond distance (`alat`) | native LKAG `J(R)` (`mRy`) |
|---:|---:|
| `0.866025` | `0.718730` |
| `1.000000` | `0.485399` |

The separate fixed bcc-Fe validation artifact
`results/validation/VAL-18_bccFe/jij_fixed/jij.out` records another accepted
two-shell state (`0.739154`, `0.467405` mRy).  It is not merged with the
post-processing reference here.  The DRESP-02C accepted projected response
state is also a separate 4x4x4, 64-k-point, 300 K `ham_only` second-order
snapshot documented in `docs/DRESP_02C_MATERIAL_CLOSURE.md:161-180`.

These provenance differences are precisely why no Fe K/L equality is claimed
from mixing rounded `J(R)` rows with a different response snapshot.

## 10. Fe \(J(q)\)

For the ordinary two-row reference, the established one-sublattice
centrosymmetric transform uses shell multiplicities 8 and 6 and

\[
 J(q)=\sum_R m_R J_R\cos(2\pi q\cdot R),
\]

with representative vectors `(-0.5,-0.5,-0.5)` and `(0,0,-1)`.  The
resulting two-shell Fourier table below is a Rung-0 transform only; it is not
a range-converged material dispersion.

| q (`2*pi/alat`) | mature LKAG `J(q)` (`mRy`) | `J(0)-J(q)` (`mRy`) |
|---|---:|---:|
| `(0,0,0)` | `8.662234` | `0.000000` |
| `(0,0,0.125)` | `7.371533` | `1.290701` |
| `(0,0,0.25)` | `4.065751` | `4.596483` |
| `(0,0,0.5)` | `-2.912394` | `11.574628` |
| `(0.125,0,0)` | `8.224553` | `0.437681` |
| `(0.25,0,0)` | `6.978145` | `1.684089` |

The calculation is reproducible with
`tools/analyze_gbt_wp11.py`; the helper's convention and the limitation of a
single two-shell range are documented in
`docs/dev/GBT_REAUDIT_WP11_LKAG_JQ.md`.  The table supplies the required
Gamma, small mesh-compatible, and additional q dependence for the mature
reference.  A same-state K/L table at those q points was not produced because
the native vertex mapping is blocked in Section 4.

## 11. K/L comparison

The finite-Hamiltonian comparison is closed at the operator level:

| comparison | outcome |
|---|---|
| explicit operator contraction vs independent finite band-pair oracle | PASS, `6.0591e-17` max residual |
| exact scalar operator vs `-1/8 Delta chi0 Delta` | PASS, below `2e-13` |
| orbital-dependent operator forced into scalar form | FAIL as intended; residual `0.15` |
| native `d_matrix(E)` LKAG vs finite `H_up-H_down` bridge | not certified; representation mismatch |
| Fe `spd` same-state `J(R)`, `J(q)`, K/L comparison | BLOCKED; no common vertex/state contract |

No fitting of `Delta`, prefactor, moment, sign, on-site term, ALSDA, or
Goldstone correction was used.  The Fe material gate therefore remains open,
not numerically “close.”

## 12. Projection dependence

`spd` is the primary projected-response selector because the mature canonical
LKAG trace is a full nine-channel `s,p,d` orbital trace.  DRESP-02 supports
both `spd` and `d`; its accepted Fe material campaign records both projected
responses, with `spd` as the complete valence comparison and `d` as a
diagnostic.

The `d` response must not be compared to the all-orbital native LKAG result as
if it were a d-only LKAG calculation.  The live canonical exchange routine
uses fixed 9x9 matrices and `d_matrix(E)` fills the supported `s,p,d`
channels.  No rigorously constructed d-only native LKAG reference is
available in this milestone.  A future masked native vertex would need its
own operator, trace, and provenance audit.

## 13. Final bridge verdict

1. **Finite-Hamiltonian bridge:** PASS, but only with the explicit
   operator-valued exchange vertex in the certified collinear orthogonal
   `ham_only` basis.
2. **Scalar site compression:** PASS only for an audited scalar operator;
   it is not generally valid and is rejected for the non-scalar fixture.
3. **Native LMTO/LKAG bridge:** BLOCKED.  Canonical LKAG uses the
   energy-dependent `symbolic_atom%d_matrix(E)` and its directed real-axis GF
   path; the finite service uses a fixed `H_up-H_down` coefficient operator.
   The repository has no proven two-energy representation map preserving the
   DRESP-01 affine endpoint structure and native on-site conventions.
4. **Fe material gate:** not passed and not hidden by tuning.  The mature
   `spd` `J(R)`/`J(q)` references are recorded, but no unsupported K/L equality
   is reported.
5. **Milestone boundary:** DRESP-04 has not been started.  The next required
   work is the missing native LMTO exchange-vertex mapping and a common-state
   Fe `spd` gate; only after a PASS outcome should the later interaction rungs
   be considered.

Implementation added for this audit:

- `source/lr_kl_static_bridge.f90` — explicit finite operator bridge and
  scalar-audit helper;
- `tests/unit/test_lr_kl_static_bridge.f90` — independent finite fixture,
  operator oracle, q/site ordering, and static ladder;
- CMake registration for `UnitLrKlStaticBridge`.

Verification on the current build:

```text
UnitLrProjectedReciprocalChi0: PASS
UnitLrKlStaticBridge: PASS
```
