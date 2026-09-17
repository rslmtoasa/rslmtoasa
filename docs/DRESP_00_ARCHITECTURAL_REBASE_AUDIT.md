# DRESP-00 — Architectural Rebase and Dynamical-Response Inventory

**Status: ARCHITECTURAL REBASE COMPLETE — IMPLEMENTATION BLOCKED AT THE
RUNG_1_KL PROJECTED SITE-SPIN CONTRACT.**

This is an audit of the live tree. No source physics, response convention,
production route, or historical evidence was modified.

## 1. Live repository provenance

| item | inspected value |
|---|---|
| branch | `fable_v4` |
| HEAD | `e9fb086624c2f2f343903c4f7880372ab76b11ba` |
| audit date | 2026-09-16 |
| target document | `docs/DRESP_00_ARCHITECTURAL_REBASE_AUDIT.md` |
| source changes | none |
| commit | not created; documentation-only audit remains uncommitted |

The pre-audit worktree already contained unrelated untracked campaign files,
including the post-LR03 prompt pack, Fe response outputs, state artifacts, and
TDVK smoke outputs. They were preserved. The after-audit worktree contains the
same pre-existing untracked files plus this document.

**Evidence: MATERIAL EXECUTION; SOURCE TRACE.** Branch, HEAD, and worktree
state were read from Git immediately before the audit document was created.

## 2. Executive verdict

The rebase produces the following architecture:

```text
accepted one-electron state
  -> coefficient-space reciprocal or native-RS GF
  -> explicit LMTO endpoint/radial Pauli primitives
  -> full point response or compact product response
  -> [missing certified site-integrated S_i^+ vertex]
  -> projected site chi0 and the first KL dynamic rung
  -> Mills/Stoner and Jülich/GSR branches
  -> later full Halle/BES ALSDA/GCR
  -> future longitudinal/charge response
```

The reciprocal Lehmann route, reciprocal real-axis GF route, and native
real-space GF route are mathematically distinct backends with explicit
representation contracts. The first two already provide the endpoint-energy
information needed for a dynamical bubble. The native route additionally
provides explicit endpoint branches and local contact terms. These are usable
as narrowly scoped infrastructure and oracles.

The missing production object is not a generic Green function. It is the
site-domain matrix representation of

\[
 \hat S_i^+=\int_{\Omega_i}d^3r\,\hat\psi^\dagger(\mathbf r)\sigma^+
 \hat\psi(\mathbf r),
\]

with a shared site moment, radial/core policy, orbital selection, and
normalization contract. Existing pointwise transition vertices and occupied
Pauli radial profiles provide most primitives for `d` and `spd`, but no live
public API certifies this site-integrated object. `spdf` is outside the live
radial and vertex capability guards.

The existing 232-dimensional product response is valid as a finite weighted
LMTO transition span and as a compact full-space diagnostic. It is not a
requirement for the first site-resolved projected route, and it is not proof
of a site-projected or material dynamical response. The existing GSR/Lounis
machinery is an independent sum-rule construction, not an ALSDA rescaling.

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; UNRESOLVED.**

## 3. Convention ledger

The live convention authority is `docs/LR_TDDFT_CONVENTIONS.md`, supplemented
by `docs/LR-GF-01_GF_CONTRACT_EVIDENCE.md` and the source guards.

| item | live convention |
|---|---|
| coefficient spinor order | site-major; within each site spin-up orbital block followed by spin-down; orbitals `(s),(p,-1:1),(d,-2:2),...` |
| circular operators | `sigma+ = (sigma_x+i sigma_y)/2`, `sigma- = (sigma_x-i sigma_y)/2`; matrix elements `(1,2)=1` and `(2,1)=1` respectively |
| number-spin | `s_z = n_up-n_down`; response `m`/`m_z` is Pauli number-spin density, not a signed physical magnetic moment |
| physical moment | `M_e = -(g mu_B/2)s`; reported code moment is `M^code=+mu_B s` |
| XC Pauli coefficient | `B_xc^sigma=(V_xc,up-V_xc,down)/2`, in Ry |
| LCMM field naming | `B_eff=V_down-V_up=-2 B_xc^sigma` for the XC-only mapped field; it is an energy-valued coefficient, not tesla |
| direct ALSDA | `K_xc^(P<-SR)=B_xc^(sigma,SR)/s^P`; the accepted route explicitly documents the SR-to-Pauli approximation |
| response order | measurement operator first; `chi_plus` uses the plus spin-flip block and `chi_minus` the reverse block |
| Fourier phase | positive live phase `exp(+i 2*pi*q.(R+tau_right-tau_left))`; reciprocal GF pair phase is `exp(+i k.(R_i-R_j))` |
| literature momentum | `q_BES=-q_RS` when the cell-origin convention is shared |
| retarded sign | `z=E+i eta`, `eta>0`; advanced uses `E-i eta` |
| response representation | point matrices are canonical LR-04 `B=chi_raw*W`; compact product matrices are weighted-orthonormal coordinates |
| loss | `L=-(chi-chi^dagger_W)/(2*i*pi)`; LR-04 metric adjoint in point space, ordinary dagger in orthonormal compact space |
| energy units | electronic energy, frequency, and broadening are Ry; output eV conversion is display-only |
| radial units/measure | radius in bohr; physical radial volume is `4*pi*r^2 dr`; LR-04 stores the explicit Simpson/volume metric |
| BZ normalization | reciprocal response uses normalized k weights; representative source factors are `2*k_weight/sum(k_weights)` in the transverse bubble |
| origin | `r=0` is an exact zero-measure response row; non-`s-s` products vanish and `s-s` uses the existing two-positive-point extension |
| core policy | reciprocal eigenpairs are valence-only; accepted frozen core is projected separately through the same accepted large-component Pauli rule and then added explicitly |

No convention repair is made here. The names `Delta`, `dele`, `B_eff`,
`B_xc^sigma`, `U_LCMM`, and `K_xc` refer to different source objects and
must not be collapsed.

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; LITERATURE-MAPPED.**

## 4. One-electron and Green-function inventory

### 4.1 Reciprocal Hamiltonian and accepted state

`source/reciprocal_fourier.f90:113-184,268-413,415-629` assembles the
site-packed reciprocal Hamiltonian with fractional/direct k coordinates and
phase `exp(+i 2*pi*k_frac.R_frac)`. In second-order mode the live effective
Hamiltonian is assembled from first-order hopping, the `hoh` subtraction, and
onsite `enim`/spin-orbit terms; optional SOC and extra operators are tracked in
state provenance. `source/hamiltonian_build.f90:1211-1326` builds the spin blocks
from `H0`, `Hz`, and `Hx +/- i Hy`.

The arbitrary-k service folds exactly into the canonical half-open fractional
zone, solves a complete eigensystem, and does not select the nearest normal
mesh point. `lr_snapshot_from_reciprocal` and `lr_q_endpoint_from_reciprocal`
make the accepted state and exact folded `k+q` endpoint explicit.

**Classification: COMMON_INFRASTRUCTURE. Reuse: REUSE_AS_CERTIFIED.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; INDEPENDENT NUMERICAL EVIDENCE.**

### 4.2 Hamiltonian resolvent

The reciprocal response service uses the complete orthogonal eigensystem to
form

\[
 G_c^{(p)}(k,z)=\sum_n\frac{\epsilon_{nk}^{p}|n k\rangle\langle n k|}
 {z-\epsilon_{nk}},\qquad p=0,1,2,
\]

in `source/lr_gf_susceptibility.f90:167-193`. The `p=0` object is an
orthogonal coefficient-space resolvent; `p=1,2` moments are required by the
affine LMTO endpoint vertex. Independent reciprocal direct inversion is
available in the Dyson/eigenpair validation path and is compared with the
Lehmann resolvent.

This is a finite-Hamiltonian coefficient propagator. It is not, by its name or
array shape, a full spatial `G(r,r')`, a screened path operator, or a physical
radial GF.

**Classification: COMMON_INFRASTRUCTURE. Reuse: REUSE_AS_CERTIFIED.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; INDEPENDENT NUMERICAL EVIDENCE.**

### 4.3 Auxiliary and screened/path-operator objects

`source/green.f90:336-473` exposes two historically named operations:

* `auxiliary_gij` applies endpoint diagonal `sqrt(Delta)` factors,
  `g_aux=D_i G_in D_j`;
* `transform_auxiliary_gij` applies the screened representation transform with
  endpoint ratios and an onsite additive term.

The live `green%gij/gji` arrays populated by reciprocal and native fillers are
coefficient-space blocks. Therefore `auxiliary_gij` is not automatically a
physical LMTO GF. The source-level `p_matrix` and screening operations exist,
but the dynamical response routes do not establish that their input is the
physical path operator `g=(P-S)^{-1}` in the sense required by the new
projected response.

This is the principal representation-name hazard in the tree. The old
auxiliary-GF exchange route is valid only under its narrow matched
orthogonal/channel contract; it must not be reused as a physical response GF
without a new representation proof.

**Classification: HYBRID_AMBIGUOUS. Reuse: REUSE_WITH_NARROW_CONTRACT.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; UNRESOLVED.**

### 4.4 Physical LMTO endpoint augmentation

`source/lmto_radial_augmentation.f90` stores accepted large/small radial
functions, derivatives, `gfac`, and energy references. The live response
contract uses the large-component Pauli projection and reconstructs

\[
 U_{l\sigma}(r,E)=\phi_{l\sigma}(r)+
 (E-\epsilon_{\nu,l\sigma}^{work})\dot\phi_{l\sigma}(r).
\]

The reciprocal point and product vertex modules use this endpoint dependence
through four affine products. The native route uses
`source/lr_gf_endpoint_augmentation.f90`, which stores `G`, `hG`, `Gh`, and
`hGh` branches and includes the algebraic contact terms from
`hG=DG-I`, `Gh=GD-I`, and `hGh=DGD-D-h^gamma`.

This is a physical/radial endpoint adapter for the supported Pauli/no-SOC
response, not a proof that `green%auxiliary_gij` has become the same object.
A complete general screened/path-operator-to-physical-GF transformation,
including all `dot P` semantics for every route, is not established.

**Classification: RUNG_4_HALLE. Reuse: REUSE_FOR_PROJECTED_CONTRACTION.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; INDEPENDENT NUMERICAL EVIDENCE.**

### 4.5 Real-space GF

The native provider contract in `source/lr_rs_gf_susceptibility.f90:35-145`
returns directed coefficient blocks `gij`, `gji`, and directed `hgamma` blocks
at arbitrary complex `z`. The reverse block is separately supplied, never
inferred by conjugation or transpose. `source/tddft_native_rsgf_provider.f90`
connects block recursion and Chebyshev to this seam; a dense inverse provider
is available as a finite exact oracle.

Block recursion and Chebyshev are approximations to the same coefficient-space
object, with independent recursion/polynomial controls. The native response
then applies endpoint augmentation, forms the bubble, applies the real-space
phase, and converts to canonical LR-04 coordinates.

**Classification: HYBRID_AMBIGUOUS. Reuse: REUSE_AS_ORACLE.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; NUMERICAL CONVERGENCE.**

## 5. Two-energy GF semantics

The live GF routes do not treat every `G` variable as interchangeable.

| route | first endpoint | second endpoint | transformation terms | local terms | conclusion |
|---|---|---|---|---|---|
| reciprocal `lr_gf_susceptibility` | spectral `A_L(E)` and advanced `G_L(E-omega)` | spectral `A_R(E)` and retarded `G_R(E+omega)` | four affine vertex components using `p,q=0,1`; resolvent moments through `p=2` | included through affine endpoint algebra, not `green%auxiliary_gij` | usable real-axis GF bubble |
| compact `lr_product_gf_susceptibility` | same two Kubo terms, matrix-backed or scalar | same | same component vertex tensor; factorized backend commutes band sums through energy integral | same as reciprocal product vertex | independent compact GF oracle |
| native `lr_rs_gf_susceptibility` | augmented directed reverse block at `E` and advanced `E-omega` | augmented directed block at `E+omega` | explicit `G`, `hG`, `Gh`, `hGh` branches at every endpoint and energy | `-I`, `-D`, and `-hgamma` contacts retained by endpoint adapter | semantically aligned full point GF candidate |
| canonical LKAG exchange | one real energy `E` | reverse intersite block at same `E` | no dynamical two-energy transformation | intersite route has no radial physical-GF contact algebra; onsite auxiliary branch is separate | static exchange only |
| auxiliary LKAG | one real energy `E` | reverse auxiliary block at `E` | `sqrt(Delta)` endpoint scaling and optional screening transform | onsite screened additive term exists in transform; no response radial augmentation | narrow auxiliary exchange route |

For reciprocal and product GF routes, the first Kubo term is explicitly
`f(E) Tr[A_L(E) V_I G_R^R(E+omega) V_J^dagger]` and the second is
`f(E) Tr[A_R(E) V_J^dagger G_L^A(E-omega) V_I]`. For the native route, the
same structure is assembled pairwise after applying endpoint branches to both
directions and multiplying the real-space phase.

Consequently, cross terms from the physical LMTO endpoint map are present in
the response GF routes. They are represented by the four affine branches, not
by an auxiliary-GF call. The current tree does not establish a corresponding
two-energy physical transformation for generic `green%auxiliary_gij`, so that
object is excluded from the projected dynamic route until re-derived.

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; INDEPENDENT NUMERICAL EVIDENCE.**

## 6. LKAG / Rung 0 inventory

The mature combined exchange path is `source/exchange.f90:1143-1326`.
At each energy it forms the canonical Pauli-decomposed matrix kernel

\[
 K_J(E)=\operatorname{ImTr}\left[
 \Delta_iG^0_{ij}\Delta_jG^0_{ji}
 -\sum_{a=x,y,z}\Delta_iG^a_{ij}\Delta_jG^a_{ji}\right]
\]

through `dGdG_Jnc`. Here `Delta_i` is the source `symbolic_atom%d_matrix`,
not literally `H_down-H_up`. `source/symbolic_atom.f90:313-341` constructs it
as an energy-dependent diagonal LMTO potential-function difference from
`c`, `dele`, `vmad`, and `E`. The source trace and output scale are therefore
the exchange module's contract, not a silently substituted textbook symbol.

The canonical path consumes `Ginmag/Gjnmag` and `Gix/Giy/Giz` with their
separately directed reverse blocks. It uses Simpson integration to the Fermi
energy and the source display scale `10^3/(4*pi)`. The same low-level kernels
also form DMI and anisotropic tensors. The two-index path is a perturbative
decomposition into charge/spin density/current pieces; its SO/FO results are
not the same observable as the combined full result.

`calculate_jij_auxgreen` at `source/exchange.f90:188-352` is different. It
forms `DeltaP=P_up-P_down`, calls `green%auxiliary_gij` for both directions,
and evaluates channel-resolved products. Its source comment says it works
better with `hoh` because it supposes an orthogonal representation. It has no
production caller in the current exchange workflow.

Reusable pre-integration objects are directed coefficient/Pauli GF blocks,
`d_matrix(E)` for canonical LKAG, `P_up/down(E)` and `DeltaP(E)` for the
auxiliary route, pair geometry, trace conventions, and accepted quadrature.
These support a future Katsnelson–Lichtenstein bridge only after its dynamic
spin vertex and representation are defined; they do not define projected
`S_i^+` response.

**Classification: RUNG_0_LKAG. Reuse: REUSE_AS_CERTIFIED for the narrow
canonical exchange contract; REUSE_AS_ORACLE for comparisons.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; INDEPENDENT NUMERICAL EVIDENCE;
MATERIAL EXECUTION.** Existing exchange validation is scoped to its recorded
GF producer, representation, mesh, and observable.

## 7. Projected-spin-operator status

There is no live certified matrix API whose public mathematical object is the
site-integrated `V_i^+`. `evaluate_pauli_transition_vertex` constructs a
complete pointwise LR-04 transition vector with the correct Pauli matrix, LMTO
radial products, Gaunt map, endpoint energy dependence, and site-local
support. `lmto_product_response_basis%component_vertex_tensor` stores the
corresponding four endpoint-energy components and is a valid primitive for a
later linear contraction.

The required site operator is not `I` in a local orbital block. It requires the
site-domain radial/angular integral, selected orbital content, and the same
Pauli/core/moment normalization used by the response.

| orbital request | live primitive coverage | projected-vertex status |
|---|---|---|
| `d` | d orbitals are present inside the supported `spd` coefficient/radial basis; no dedicated d-only site-domain API or mask is exposed | `PARTIAL / RE-DERIVATION_REQUIRED` |
| `spd` | complete supported large-component radial snapshot, Pauli vertex, Gaunt products, and accepted radial mesh exist | `PARTIAL / RE-DERIVATION_REQUIRED`; primitives are sufficient, site-integrated contract is not certified |
| `spdf` | live radial augmentation, product basis, and response vertex guards support only `lmax<=2` | `BLOCKED_BY_MISSING_GROUND_STATE_DATA` |

The occupied Pauli projection in `source/pauli_ground_state_projection.f90`
does provide radial `m_z(site,r)` and explicit valence/core components. Thus
raw supported-state moment data exist, but a shared site-integrated
response/moment API is still absent.

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; UNRESOLVED.**

## 8. Bare-susceptibility inventory

| route | mathematical object and indices | representation/normalization | validation and status |
|---|---|---|---|
| `lr_ks_susceptibility` | complete reciprocal band-pair Lehmann `chiKS(I,J,q,omega)` with `omega+E_L-E_R+i eta`, occupation difference, exact folded `k+q` endpoint | full site × `(L,M)` × radial × channel; canonical `B=chi_raw*W`; `chi_plus`/`chi_minus` | exact finite-state reference for declared Pauli/no-SOC baseline; not site-projected |
| `lr_gf_susceptibility` | reciprocal real-axis Kubo bubble with spectral, retarded, and advanced resolvents | same full point LR-04 space; Simpson integral, separate integration and physical broadening | independent GF route against Lehmann only when endpoint, mesh, eta, and state contracts match |
| `lr_rs_gf_susceptibility` | native pairwise real-space GF bubble followed by phase assembly | same full point LR-04 space; explicit directed pair set and site gauge | endpoint augmentation/full-route seams exist; provider convergence and projected site equivalence remain open |
| `lr_product_ks_susceptibility` | exact band-pair response after compact product transition coordinates | weighted-orthonormal retained product coordinates; no point response allocation | product-span and compact response evidence; not site-indexed |
| `lr_product_gf_susceptibility` | same real-axis GF bubble in compact product coordinates; scalar, optimized, and factorized contractions | compact 232 for accepted one-site `spd`, with GF controls | independent compact GF oracle; no KXC/Dyson implied by module |
| legacy `green`/exchange consumers | one-energy intersite exchange or auxiliary quantities, not `chiKS` | coefficient/Pauli or auxiliary exchange arrays | mature static exchange evidence only |

The exact point and product routes are not the same claim as

\[
 \bar\chi_{0,ij}^{+-}(q,\omega)=
 \langle V_i^+\,\chi_0(q,\omega)\,V_j^+\rangle
\]

until the site-integrated `V_i^+` and site moment are fixed. Once that linear
functional is certified, contracting an exact point or retained product
transition vector is a representation change. Choosing a finite orbital/site
domain or dropping radial/angular content without that derivation would
instead be a physical projection or approximation. No site-only fallback is
introduced by this audit.

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; INDEPENDENT NUMERICAL EVIDENCE;
UNRESOLVED.**

## 9. Kernel taxonomy

### 9.1 LKAG exchange vertex — Rung 0

The canonical exchange vertex is `symbolic_atom%d_matrix(E)` and the
auxiliary variant is `DeltaP(E)`. They are energy-dependent static-exchange
objects. Neither is the dynamic site spin vertex.

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT.**

### 9.2 Mills/Stoner mapping — Rung 2

No genuine live implementation of a projected Stoner mapping
`U_i^split ~ Delta_i/M_i` was found. `NOT IMPLEMENTED`. The direct ALSDA ratio
`B_xc^sigma/m` is not silently relabelled as Mills, and LKAG `Delta` is not
silently divided by a moment to make a Stoner interaction.

**Classification: RUNG_2_MILLS. Reuse: BLOCKED.**

**Evidence: SOURCE TRACE; UNRESOLVED.**

### 9.3 Jülich / LCMM sum-rule interaction — Rung 3

`source/lr_goldstone_sumrule.f90` solves an explicitly independent static
equation. In the active positive-measure spherical sector it forms

\[
 \Gamma_{ab}=4\pi\,\chi_{0,ab}(0)\,m_{00,b},\qquad
 \Gamma U_{LCMM}=m_{00},
\]

then constructs a canonical local interaction `K_eff=4*pi*U_LCMM`. The
unknown is site × positive radial point; the full input susceptibility may
couple the full response space, but the solved field is an L=0 radial field.
The origin is excluded and rank deficiency blocks promotion.

`lr_compact_static_interaction.f90` is the corresponding compact product-space
static diagnostic. It solves a retained L=0 local ansatz and checks its action
against compact magnetization. It is not a scalar rescaling of direct ALSDA.

**Classification: RUNG_3_JULICH. Reuse: REUSE_AS_ORACLE.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; INDEPENDENT NUMERICAL EVIDENCE.**

### 9.4 Full ALSDA — Rung 4

`source/lr_alsda_kernel.f90` forms the explicit mixed contract

\[
 K_{xc}^{P\leftarrow SR}(r)=
 \frac{(V_{xc,up}^{SR}-V_{xc,down}^{SR})/2}{s^P(r)}
\]

on the accepted radial mesh, validates XC provenance, rejects zero magnetization
on positive-measure points, and maps the scalar field to canonical LR-04
operator form. It has no GSR, Dyson, site reduction, or empirical correction
inside the module. The compact driver may project this operator as
`U^H K_point U`, but that is not proof that the projected dynamic route has
been built.

**Classification: RUNG_4_HALLE. Reuse: DEFER_TO_HALLE for the full production
claim; REUSE_WITH_NARROW_CONTRACT for diagnostics.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; MATERIAL EXECUTION.**

### 9.5 Goldstone/BES correction

`source/lr_goldstone_correction.f90` is a separate optional eigenvalue repair:
it forms `D=I-chiKS(0)Kxc`, identifies one rigid mode under explicit isolation,
overlap, eta, k-mesh, provenance, and no-SOC/field gates, sets only that
eigenvalue to zero, and reconstructs a corrected kernel. The production driver
rejects `goldstone_correction=true` unless separate validation evidence is
supplied. It remains separate from both ALSDA and Jülich/GSR.

**Classification: RUNG_4_HALLE. Reuse: DEFER_TO_HALLE.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; UNRESOLVED for production use.**

### 9.6 Longitudinal/charge — Rung 5

No live coupled charge/longitudinal TDDFT kernel exposing the required second
derivatives and charge-spin response blocks was found. The production driver
accepts only `chi_plus` or `chi_minus`; transverse ALSDA is not a longitudinal
kernel.

**Classification: RUNG_5_LONGITUDINAL. Reuse: BLOCKED.**

**Evidence: SOURCE TRACE; UNRESOLVED.**

## 10. Dyson and observable inventory

`source/tddft_dyson.f90` consumes already-canonical operators and implements

\[
 D=I-\chi_{KS}K,\qquad D\chi=\chi_{KS},
\]

with LAPACK `zgesv`, no explicit inverse, per-frequency residuals, singular
values, condition numbers, minimum eigenvalue magnitudes, and numerical-pole
flags. Full point space uses LR-04 canonical composition; compact orthonormal
space uses ordinary matrix multiplication.

Loss is computed as the metric-adjoint anti-Hermitian spectral product in point
space and ordinary-dagger form in compact orthonormal space. The module does
not find modes, fit linewidths, track modes across q, or apply pole
regularization.

The production driver has two architectural uses: the ordinary sweep can run
full point Lehmann, reciprocal GF, native-RSGF, direct ALSDA, and GSR paths;
`backend=compact_dyson` runs a dense retained product-space Dyson problem,
currently with a complete rank-stable 232-mode `spd` basis for accepted Fe.
The latter is a real existing dynamical route, but its final index is product
coordinate rather than the new initial site/sublattice index.

No current TDDFT code was found that infers a dispersion from a loss-trace
peak alone. The separate `frozen_magnon` post-processing path computes an
energy-surface observable and its own branch matrix; it is not dynamic loss
analysis and must not validate `chiKS` or projected Dyson modes. Existing
TDDFT reports defer mode extraction, linewidths, stiffness, and material
spectrum claims.

**Classification: COMMON_INFRASTRUCTURE. Reuse: REUSE_AS_CERTIFIED for the
declared canonical/compact contracts; REUSE_AS_ORACLE under the rebase.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; INDEPENDENT NUMERICAL EVIDENCE.**

## 11. Reclassification of the 232-dimensional product response

The live product basis is constructed from accepted radial endpoint products,
Gaunt triangle/parity rules, and direct weighted SVD. For one-site `spd`,
`response_lmax=4`, the exact unpruned/retained count is

\[
 1(12)+3(16)+5(16)+7(8)+9(4)=232.
\]

| question | DRESP-00 answer |
|---|---|
| exact within declared LMTO transition space? | **Yes, as a finite span claim.** Candidate construction and synthetic point-to-product closure are certified; accepted-Fe rank stability is recorded. This is not exactness of the full continuum response. |
| retains radial information for full ALSDA? | **It retains radial product information for the retained span and supports `U^H K_point U`.** It does not prove local ALSDA closure for every physical field or full material convergence. |
| site-spin vertices contract exactly through it? | **Potentially as a linear contraction once site-domain `V_i^+` is derived.** The current basis does not expose or certify that functional. |
| useful GF oracle? | **Yes.** `lr_product_gf_susceptibility` is an independent real-axis GF implementation with scalar/optimized/factorized options and comparison to product Lehmann. It is not yet a site-projected oracle. |
| necessary for projected Jülich? | **No.** A projected site `chi0_ij` and site `U_i` can be represented directly in a small site matrix. The 232 Dyson problem is optional. |
| natural later Halle/BES role? | **Yes.** Its full site/angular/radial product semantics are naturally useful for later full-spatial ALSDA/BES, subject to closure evidence. |
| previous evidence that survives? | **Product count, weighted-SVD/rank evidence, synthetic span closure, endpoint-affine identities, product GF algebra, and compact Dyson residual machinery survive within their original objects.** |
| evidence that does not transfer? | **No automatic transfer to projected site response, Mills/Stoner, Jülich site response, Goldstone compliance, material magnons, or longitudinal response.** |

**Classification: HYBRID_AMBIGUOUS. Reuse: REUSE_FOR_PROJECTED_CONTRACTION
and REUSE_AS_ORACLE; DEFER_TO_HALLE for full-space claims.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; INDEPENDENT NUMERICAL EVIDENCE;
NUMERICAL CONVERGENCE; UNRESOLVED.**

## 12. Real-space route classification

Previous native-RSGF work can supply the same *coefficient-space input*
definition used by the reciprocal route. The native provider returns directed
blocks, the endpoint adapter applies the same affine radial/Pauli map, the
bubble uses the same two-energy ordering, and the response backend applies the
same canonical phase/metric conversion.

The remaining semantic gap is the same as in reciprocal space: no certified
site-integrated `V_i^+` and no demonstrated equality between projected native
and projected reciprocal responses. Additional native conditions remain:
complete pair coverage, serial replicated workset in the current adapter,
explicit translation/site-gauge provenance, and recursion/Chebyshev
convergence. The current driver initializes native site positions to zero, so
it does not establish a general nontrivial sublattice gauge.

The route is a credible future DRESP projected-response backend, not a reason
to restart native integration in DRESP-00.

**Classification: HYBRID_AMBIGUOUS. Reuse: REUSE_WITH_NARROW_CONTRACT and
REUSE_AS_ORACLE.**

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; NUMERICAL CONVERGENCE;
UNRESOLVED.**

## 13. Complete component ledger

| Component | File/routine | Mathematical object | Representation | Basis/metric | Spin/channel | Energy arguments | Units/normalization | Existing evidence | New ladder classification | Reuse verdict | Open issue |
|---|---|---|---|---|---|---|---|---|---|---|---|
| response layout | `lr_response_space.f90` | finite coordinate/metric algebra | LR-04 point | site, angular, radial, channel; right weighted | one explicit circular channel | static coordinates | radial volume metric; origin null row | source and metric tests | `COMMON_INFRASTRUCTURE` | `REUSE_AS_CERTIFIED` | no site projection API |
| index/gauge map | `response_basis_mapping.f90` | super-index and phase map | point/real-space | canonical site-`L,M,r,channel` | plus/minus labels external | q and translations | positive Fourier phase | source/docs | `COMMON_INFRASTRUCTURE` | `REUSE_AS_CERTIFIED` | site operator must use same gauge |
| reciprocal assembler | `reciprocal_fourier.f90` | `H(k)` and exact arbitrary-k endpoint | coefficient space | site-major orthogonal baseline | spin blocks up/down; optional tracked SOC | k and k+q | Ry; fractional k | dense inverse/arbitrary-k tests | `COMMON_INFRASTRUCTURE` | `REUSE_AS_CERTIFIED` | generalized-overlap response unsupported |
| coefficient resolvent | `lr_gf_susceptibility.f90`, `reciprocal` | `G_c^(p)(k,z)` | coefficient-space Hamiltonian resolvent | orthogonal | spinor blocks; Pauli later | E, E±omega; p=0..2 | Ry and `+i eta` | independent inverse and causal tests | `COMMON_INFRASTRUCTURE` | `REUSE_AS_CERTIFIED` | not physical spatial GF alone |
| auxiliary endpoint scale | `green.f90:auxiliary_gij` | `D_i G D_j`, `D=sqrt(Delta)` | auxiliary candidate | coefficient endpoint blocks | channel-resolved | one energy | source LMTO scaling | representation oracle | `HYBRID_AMBIGUOUS` | `REUSE_WITH_NARROW_CONTRACT` | historical physical-GF name unsafe |
| screened transform | `green.f90:transform_auxiliary_gij` | screened/path representation transform | auxiliary/screened | diagonal spherical P; onsite additive term | up/down channels | one energy | representation-dependent | onsite/offsite inverse tests | `RUNG_0_LKAG` | `REUSE_WITH_NARROW_CONTRACT` | no dynamic physical proof |
| radial snapshot | `lmto_radial_augmentation.f90` | `phi`, `phidot`, `enu_work`, radial channels | physical endpoint primitive | accepted large component and mesh | collinear no SOC | endpoint E | bohr; large-component Pauli | radial/augmentation evidence | `COMMON_INFRASTRUCTURE` | `REUSE_WITH_NARROW_CONTRACT` | lmax limited to spd |
| point Pauli vertex | `lr_pauli_transition_vertex.f90` | full point `T_I=<L|sigma±|R>` | LR-04 point | Gaunt plus explicit radial products | `chi_plus`/`chi_minus` | left/right endpoint E | unweighted point value; metric later | point vertex tests | `RUNG_4_HALLE` | `REUSE_FOR_PROJECTED_CONTRACTION` | no site-integrated public vertex |
| product basis | `lr_lmto_product_response_basis.f90` | finite LMTO product span | compact weighted-orthonormal | retained radial product modes | one channel per basis | affine endpoint components | `U=sqrt(W)V` | 232 count, SVD, closure | `HYBRID_AMBIGUOUS` | `REUSE_FOR_PROJECTED_CONTRACTION` | full material response not proven |
| Lehmann chiKS | `lr_ks_susceptibility.f90` | band-pair bare susceptibility | full point | canonical LR-04 | plus/minus | `omega+E_L-E_R+i eta` | normalized k weights | exact finite-state route | `RUNG_4_HALLE` | `REUSE_AS_ORACLE` | no site final index |
| reciprocal GF chiKS | `lr_gf_susceptibility.f90` | real-axis Kubo bubble | full point | canonical LR-04 | plus/minus | E, E±omega, retarded/advanced | Simpson; two eta controls | GF/Lehmann crosscheck scope | `RUNG_1_KL` | `REUSE_AS_ORACLE` | projected vertex missing |
| compact Lehmann | `lr_ks_susceptibility.f90` product path | compact bare response | product coordinates | weighted-orthonormal 232 baseline | plus/minus | same band-pair denominator | compact matrix | product evidence | `HYBRID_AMBIGUOUS` | `REUSE_AS_ORACLE` | not site indexed |
| compact GF | `lr_product_gf_susceptibility.f90` | compact real-axis Kubo bubble | product coordinates | weighted-orthonormal | plus/minus | E, E±omega | Simpson; factorized/scalar/optimized | independent contraction oracle | `RUNG_1_KL` | `REUSE_AS_ORACLE` | projection still unimplemented |
| native provider | `tddft_native_rsgf_provider.f90` | directed coefficient `gij/gji` | native RS coefficient GF | site blocks; serial pair set | spinor coefficient blocks | arbitrary complex z | recursion/Chebyshev controls | provider integration evidence | `COMMON_INFRASTRUCTURE` | `REUSE_WITH_NARROW_CONTRACT` | pair/gauge/convergence closure |
| native response | `lr_rs_gf_susceptibility.f90` | pairwise augmented GF bubble | full point LR-04 | radial/angular metric after phase | plus/minus | E, E±omega | Simpson; explicit contacts | endpoint/full-route evidence | `RUNG_1_KL` | `REUSE_AS_ORACLE` | no projected site equivalence |
| canonical LKAG | `exchange.f90:dGdG_Jnc`, `calculate_exchange` | static J/D/A intersite kernels | canonical Pauli GF arrays | 9x9 orbital trace | global spin frame | one E | Simpson to EF; source scale | exchange validation map/material runs | `RUNG_0_LKAG` | `REUSE_AS_CERTIFIED` | not dynamic site chi |
| auxiliary LKAG | `exchange.f90:calculate_jij_auxgreen` | `DeltaP g_aux DeltaP g_aux` | auxiliary orthogonal route | orbital/channel trace | collinear up/down | one E | Simpson to EF | source only; no production caller | `RUNG_0_LKAG` | `REUSE_WITH_NARROW_CONTRACT` | representation match required |
| direct ALSDA | `lr_alsda_kernel.f90` | `B_xc^sigma/s^P` local kernel | full point canonical | diagonal local operator with LR-04 metric | Pauli projected | static radial | Ry bohr^3 | provenance/low-m/origin guards | `RUNG_4_HALLE` | `DEFER_TO_HALLE` | distinct from projected/Jülich kernel |
| GSR/LCMM | `lr_goldstone_sumrule.f90` | `Gamma U=m`, `K_eff=4pi U` | full point with L=0 solve | site × positive radial unknowns | Pauli `m_00` | static omega=0 | Ry bohr^3 | SVD rank/residual tests | `RUNG_3_JULICH` | `REUSE_AS_ORACLE` | no site-indexed projected contract |
| compact GSR | `lr_compact_static_interaction.f90` | compact local sum-rule action | product coordinates | L=0 retained compact span | Pauli magnetization | static | compact operator | action/residual diagnostics | `RUNG_3_JULICH` | `REUSE_AS_ORACLE` | not necessary for projected route |
| BES/GCR | `lr_goldstone_correction.f90` | one-mode denominator correction | canonical full point | active metric/non-Hermitian classification | rigid transverse mode | static q=0 | Ry bohr^3 | gate logic/reconstruction tests | `RUNG_4_HALLE` | `DEFER_TO_HALLE` | production selection blocked |
| Dyson/loss | `tddft_dyson.f90` | `D=I-chi0 K`, enhanced chi, loss | full canonical or compact orthonormal | metric adjoint in point space | supplied channel | frequency sweep | residual/condition diagnostics | algebraic and driver tests | `COMMON_INFRASTRUCTURE` | `REUSE_AS_ORACLE` | no mode/linewidth pipeline |
| production adapter | `tddft_production_driver.f90` | state handoff and route orchestration | full point or compact | accepted state; 232 compact branch | transverse only | q, omega, eta | provenance headers | smoke/static/product evidence | `HYBRID_AMBIGUOUS` | `REUSE_WITH_NARROW_CONTRACT` | full-product assumptions conflict |
| Pauli moment | `pauli_ground_state_projection.f90` | `m_z(site,r)` plus valence/core | radial physical profile | accepted large-component radial measure | up/down | occupied E through EF | electrons bohr^-3; integrated Pauli number | SR/Pauli diagnostic/static use | `COMMON_INFRASTRUCTURE` | `REUSE_WITH_NARROW_CONTRACT` | share site integral with `V_i` |
| frozen magnon | `calculation.f90:frozen_magnon*` | energy-surface/q branch diagnostic | separate MFT/SCF route | sublattice cone-angle space | spin-spiral, not `chi±` | q and energy differences | separate output normalization | legacy documented path | `OBSOLETE_UNSUPPORTED` for DRESP chi claims | `OBSOLETE` as projected evidence | do not use as loss/mode oracle |
| longitudinal kernel | no live implementation | coupled charge/longitudinal second derivatives | absent | absent | only transverse driver guard | absent | absent | no evidence | `RUNG_5_LONGITUDINAL` | `BLOCKED` | new route required |

## 14. Architecture mismatch table

| historical/live assumption | conflict under new ladder | disposition |
|---|---|---|
| full site × angular × radial response precedes any production response | new initial final index is site/sublattice for projected `S_i^+` response | keep full point/product infrastructure; do not treat it as prerequisite |
| `backend=compact_dyson` is the production dynamical response | it solves a 232-coordinate product Dyson problem, not initial site-indexed projected Dyson | retain as oracle and later Halle-capable route |
| native provider says site-only fallback is forbidden | that guard concerns incomplete real-space pair worksets, not the physical site-domain operator | do not reinterpret it as prohibition on DRESP projection |
| `green%auxiliary_gij` comments call input physical | live fillers provide coefficient-space `gij/gji`; endpoint scaling is not radial Pauli augmentation | require explicit representation provenance |
| `Delta`, `dele`, `DeltaP`, and `B_xc^sigma` share exchange-splitting language | they are different LMTO, auxiliary, and XC objects with different energy/unit roles | preserve names and derive bridges separately |
| GSR/Lounis code sits next to ALSDA code | GSR solves independent `Gamma U=m`; ALSDA forms `B_xc^sigma/s^P` | classify GSR as RUNG_3 and ALSDA as later RUNG_4 |
| compact product closure proves physical response closure | closure proves retained transition span and representation algebra only | do not inherit projected, Goldstone, or material claims |
| historical LR-REP text says native response is blocked | live source now contains endpoint augmentation and native bubble; old text is dated | use current source for implementation and old docs as dated evidence |
| frozen-magnon output is called a dispersion | it is an energy-surface diagnostic, not a TDDFT loss pole or projected `chi0` | keep outside this response ladder |

**Evidence: SOURCE TRACE; UNRESOLVED.** No mismatch is repaired in DRESP-00.

## 15. Surviving previous certifications

The following evidence remains valid because the mathematical object is
unchanged:

| certification | surviving scope |
|---|---|
| reciprocal Hamiltonian/eigenpair and arbitrary-k tests | orthogonal reciprocal coefficient-space H and exact folded endpoint service |
| Lehmann versus direct reciprocal inverse | declared coefficient-space resolvent, not a physical spatial GF |
| block-recursion/Chebyshev evidence | native coefficient-GF approximation and convergence controls |
| radial augmentation and Pauli point-vertex tests | supported `sp`/`spd`, collinear, no-SOC large-component endpoint algebra |
| LR-04 response metric/index/phase tests | point canonical operator algebra and origin treatment |
| product candidate count and 232 rank stability | finite `spd`, `Lmax=4`, weighted product representation |
| synthetic product-span closure | point transition vectors on synthetic fixture, both circular channels |
| compact GF scalar/optimized/factorized comparisons | compact GF bubble implementation for compact product object |
| LKAG exchange validation | static J/D/A observables for exact GF, Delta, mesh, and output conventions |
| GSR SVD/residual oracles | independent full/compact static sum-rule equations as implemented |
| Dyson residual/loss tests | canonical point and compact orthonormal matrix operations |
| accepted Pauli ground-state projection | radial valence/core decomposition and SR-versus-Pauli diagnostic |

None automatically proves the new projected site susceptibility. Product
closure and Lehmann/GF equality must be re-run after applying the same
site-integrated vertex on both sides.

**Evidence: INDEPENDENT NUMERICAL EVIDENCE; NUMERICAL CONVERGENCE; SOURCE TRACE.**

## 16. Invalidated or deferred assumptions

The following historical assumptions are not carried into the new production
ladder:

* full 232-dimensional response is the first required production index;
* a site-indexed response is necessarily an impermissible approximation;
* direct ALSDA and GSR are interchangeable names for one kernel;
* any array called `gij`, `aux_gij`, or `G` is already a physical response GF;
* static LKAG `Delta` can be reused as Mills/Stoner `U` without a mapping;
* compact Dyson proves projected Jülich or full Halle/BES physics;
* historical material validation applies after changing the final response
  object from full point/product to site-projected;
* `spdf` support follows from `spd` support;
* frozen-magnon energy output validates dynamic loss poles.

Deferred items are the site-integrated projected vertex, a shared projected
moment, projected reciprocal GF/Lehmann equivalence, projected native-RS
equivalence, a genuine Mills mapping, a site-resolved Jülich interaction if
needed, full Halle/BES closure, and longitudinal response.

**Evidence: UNRESOLVED; SOURCE TRACE.**

## 17. Blockers and gaps

These are hard-stop findings required by the audit, not implementation
failures of DRESP-00:

1. **Projected site-spin vertex is not certified.** A pointwise Pauli vertex
   exists, but the live tree lacks the site-domain `V_i^+` matrix/functional
   with explicit orbital selection and normalization.
2. **Projected response/moment sharing is not certified.** The occupied Pauli
   radial profile and core/valence split exist, but no shared site integration
   contract guarantees that `V_i^+`, `M_i`, and later Mills/GSR equations use
   exactly the same domain, radial measure, and core policy.
3. **`spdf` is blocked by live capability/data guards.** The radial snapshot,
   Pauli vertex, and product basis are currently limited to `lmax<=2`.
4. **Auxiliary versus physical GF semantics remain unsafe outside the narrow
   exchange contract.** `auxiliary_gij` is endpoint scaling; a complete dynamic
   physical path-operator transformation is not established.
5. **The Jülich/GSR object is not a site-projected Mills interaction.** Its
   current solve is radial L=0, full/product-space dependent, and independent
   of ALSDA. A projected site version must be specified separately if needed.
6. **Native projected equivalence is not certified.** Pair coverage, provider
   convergence, translation/site gauge, and equality to projected reciprocal
   response still require a later route test.
7. **No longitudinal/charge kernel exists.** The current response driver is
   transverse-only.
8. **No new projected material validation exists.** Existing evidence belongs
   to full point/product, static exchange, or documented diagnostics and cannot
   be transferred to the new projected object.

The two-energy semantics themselves are *not* a blocker for reciprocal and
native GF routes: they are explicitly represented. The blocker is the
operator to which those two-energy GFs must be contracted.

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; UNRESOLVED.**

## 18. Lowest unresolved layer

The lowest unresolved layer is:

> **RUNG_1_KL: certified projected site-spin operator and shared site-moment
> contract, immediately before projected `chi0_ij^{+-}(q,omega)`.**

This is below Mills/Stoner, Jülich/GSR, ALSDA, BES/GCR, Dyson, or loss claims.
It is above the already certified coefficient Hamiltonian, coefficient
resolvent, radial endpoint primitives, and point/product transition maps.

For supported `d`/`spd`, this is primarily a missing public contract and
certification seam rather than evidence that raw primitives do not exist. For
`spdf`, it is additionally missing live radial/augmentation capability.

**Evidence: SOURCE TRACE; UNRESOLVED.**

## 19. Recommended DRESP-01 scope

DRESP-01 should implement and certify exactly the first projected reciprocal
KL step:

1. define a site-domain integrated `V_i^+`/`V_i^-` matrix from accepted LMTO
   large-component radial functions, Pauli matrices, orbital selection,
   Gaunt/angular normalization, and radial measure;
2. define the matching `M_i`/Pauli number-spin site integral from the existing
   valence-plus-frozen-core radial projection, with one explicit core policy;
3. evaluate site-indexed reciprocal Lehmann
   `bar_chi0_ij^{+-}(q,omega)` on accepted `d`/`spd` states and exact folded
   `k+q` endpoints;
4. compare the same projected object against the existing reciprocal real-axis
   GF route, holding state, q, eta, endpoint map, and vertex fixed;
5. include algebraic fixtures and a material smoke reporting site-indexed
   matrices, convention provenance, and projection residuals;
6. leave the 232 product basis available as an optional contraction/oracle,
   but do not make a 232-dimensional Dyson solve a DRESP-01 prerequisite;
7. exclude Mills/Stoner, GSR/Jülich interaction construction, ALSDA changes,
   BES/GCR correction, Dyson/loss changes, native-RS production integration,
   `spdf`, and longitudinal/charge response from DRESP-01.

The acceptance gate should be equality of reciprocal Lehmann and reciprocal GF
*after the same site projection*, plus a documented projected moment/core
agreement. Only after that gate should the site response define a Mills mapping
or projected Jülich interaction.

**Evidence: SOURCE TRACE; ALGEBRAIC CONTRACT; UNRESOLVED.** This is a scoped
recommendation from the live dependency graph, not an implementation performed
by DRESP-00.

## Final status

`ARCHITECTURAL REBASE COMPLETE — IMPLEMENTATION BLOCKED AT RUNG_1_KL
PROJECTED SITE-SPIN VERTEX / SITE-MOMENT CONTRACT.`

No source physics was modified. Only this audit document was generated; no
commit was created.
