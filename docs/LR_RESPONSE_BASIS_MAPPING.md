# LR-02/02R radial/angular response-basis mapping

## Verdict

**PARTIAL PASS** for the formally defined Pauli/no-SOC radial/angular response
mapping, with two deliberately separate capability statements:

| target | disposition |
| --- | --- |
| BES/Halle-style Pauli, no-SOC radial/angular response | FORMALLY FEASIBLE |
| exact response of the present scalar-relativistic RS-LMTO density operator | BLOCKED / not yet derived |
| fully relativistic `j,kappa`-resolved response | OUT OF CURRENT BASELINE SCOPE |

The complex angular convention, Gaunt/product algebra, derived response cutoff,
production logarithmic-mesh measures, canonical response super-index, and
endpoint phase algebra have small physics-owned utilities and independent
oracles. The formal Pauli transition vector
`T^P_mu_nm;aLMi(k,q)` is constructible from the certified large radial
component, LMTO coefficients, explicit Pauli operator, and Gaunt coefficients.
This is a response-space mapping statement, not a claim that it reproduces the
full scalar-relativistic `NEWRHO` density. [BASIS MAPPING]

The remaining required numerical closure is to quantify the scalar-relativistic
ground-state density/moment difference from the Pauli projection for the LR-01
Fe fixture, preferably with Ni as a second case. No arbitrary pass tolerance is
assigned in LR-02R. [DEFERRED]

No `chiKS`/`chi0`, Dyson equation, XC kernel, Ward/Goldstone correction, or
mode extraction was added.

## Audit identity and prerequisite evidence

The audit began on branch `fable_v4` at exact starting HEAD
`3d43fa10e42efdf47fe1e38e67e73735a1f439d7d`. The starting worktree was clean.
The current branch remains `fable_v4`; the worktree is dirty only because of
the focused LR-02 additions described here.

The required prerequisite evidence was present before this change:

| prerequisite | evidence | focused tests present |
| --- | --- | --- |
| LR-GF-01 | [`docs/LR-GF-01_GF_CONTRACT_EVIDENCE.md`](LR-GF-01_GF_CONTRACT_EVIDENCE.md) | `UnitGreenLifecycle` and LR-GF validation sources |
| LR-BASIS-00 | [`docs/LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md`](LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md) | `UnitLrBasisRadial`, `UnitLrBasisAngular`, `UnitLrBasisAugmentation` |
| LR-01 | [`docs/LR_RADIAL_GROUND_STATE_AUDIT.md`](LR_RADIAL_GROUND_STATE_AUDIT.md) | `UnitLrRadialGroundState`, `Val22LrRadialGroundState` |

The literature target IDs are retained in
[`docs/dev/RS_LMTO_TDDFT_cleanroom_Luna/00A_LITERATURE_CONTRACTS.md`](dev/RS_LMTO_TDDFT_cleanroom_Luna/00A_LITERATURE_CONTRACTS.md):
BES-05, LCMM-02, and EEB-02/03/04. The older LR-02 feasibility file is treated
as a blueprint, not as proof of the live implementation.

## Scope

The audit scope is the already-certified tuple: reciprocal `ham_only`,
orthogonal second-order/HOH LMTO, collinear ground state, no SOC, no Hubbard or
other additive operator, scalar-relativistic radial functions, ordinary
periodic crystal, and `sp`/`spd` orbital spaces. Noncollinearity, SOC,
generalized overlap, additive operators, and `spdf` remain deferred.
[NUMERICAL REPRESENTATION] [DEFERRED]

## Physical target and literature object

The common ASA/KKR-style target is a site-, response-angular-, radial-, and
density/spin-channel-resolved kernel

\[
 \chi^{\mu\nu}_{ab;\Lambda\Lambda'}(r,r';\mathbf q,\omega),
 \qquad \Lambda=(L,M),
\]

with the radial and angular coordinates retained before any site-level
reduction. [LITERATURE]

For two one-electron states, the local transition field in sphere `a` is

\[
 \rho^\mu_{nm,a}(\mathbf r;\mathbf k,\mathbf q)
 =\Psi^\dagger_{n\mathbf k}(\mathbf r)\,\Gamma^\mu\,
  \Psi_{m,\mathbf k+\mathbf q}(\mathbf r),
\]

where the operator is supplied explicitly; the algebraic density/spin choices
are `Gamma^0=I` and `Gamma^(x,y,z)=sigma_(x,y,z)`. [LITERATURE]

The direct-mesh response amplitude that would populate the target is

\[
 T^\mu_{nm;a\Lambda i}(\mathbf k,\mathbf q)
 =\int d\Omega\,Y_\Lambda^*(\hat r)\,
   \Psi^\dagger_{n\mathbf k}(r_i,\hat r)\Gamma^\mu
   \Psi_{m,\mathbf k+\mathbf q}(r_i,\hat r).
\]

This is the full physical object whose computability, not a susceptibility
spectrum, was the original LR-02 pass criterion. LR-02R now assesses the
explicit Pauli projection separately and retains exact-SR closure as a distinct
blocked capability. [LITERATURE] [BASIS MAPPING]

The revised formal target is the Pauli/no-SOC member of that family. The cited
BES/Halle-style construction uses ordinary charge and Pauli spin operators and
ordinary radial Schrödinger-like regular/irregular solutions multiplied by
spherical harmonics; relativistic effects and diamagnetic response are outside
that target theory. [LITERATURE]

This distinction is consistent with the scalar-relativistic reference
description: scalar relativity retains selected relativistic effects while
omitting spin-orbit splitting and treating minor-component charge effects as an
approximation. See the [NIST ScRLDA description](https://math.nist.gov/DFTdata/atomdata/node20.html)
and the [Halle/KRR response formulation](https://arxiv.org/pdf/2603.03220).
[LITERATURE]

## Reused one-electron contracts

LR-BASIS-00 remains the source of truth for radial provenance, LMTO
augmentation, coefficient normalization, spin-block layout, and the
no-extra-operator first-order reconstruction. LR-01 remains the source of
truth for the accepted radial mesh, spin-channel naming, weighted `RHO`,
physical `n_up/down`, XC arrays, and accepted-state lifetime. Those audits were
not duplicated here.

The live `RSEQSR` comment states that `G(:,1)` is the large-component radial
numerator `U_l(r)` in `psi=(U/r)Y_lm`; the routine normalizes the packed
large/small pair and stores it in `G(:,1:2)`. This is corroborated by the
normalization loop and by the matching `GINTSR`/`NEWRHO` radial expressions,
not taken from the comment alone. [BASIS MAPPING]

For nonzero radius, the production scalar-relativistic radial metric is

\[
 \mathrm{GFAC}_l(r)=1+\frac{l(l+1)}{[TMC(r)r]^2},
 \qquad
 \mathcal N_l(r)=\mathrm{GFAC}_l(r)G_l(r)^2+G_{s,l}(r)^2.
\]

`RSEQSR`, `GINTSR`, and `NEWRHO` use this same metric for spherical radial
normalization and density accumulation. [BASIS MAPPING]

`lmto_radial_augmentation` stores the large and small radial arrays and their
energy derivatives, and reconstructs the first-order radial amplitudes. Its
documented `reconstruct_state` output deliberately does not multiply by
`Y_lm`; the stored scalar `G_s` is not accompanied by a live definition of the
full lower-component spin-angular function, its local spinor phase, or the
cross-angular operator matrix elements.
[BASIS MAPPING] [DEFERRED]

The live `NEWRHO` path also stores `phi_amp=G(:,1)` for its B6 radial matrix
elements and comments that the packed small component is a relativistic
correction not needed there. That is useful provenance for the existing radial
feature, but it is not a definition of the missing pointwise lower-component
angular field. [BASIS MAPPING] [DEFERRED]

### Exact scalar-relativistic boundary

The production scalar-relativistic contract certifies an angular-integrated
spherical norm. It does not expose the full scalar-relativistic density
operator needed to evaluate, for arbitrary `Gamma`, the pointwise
minor-component terms in `Psi^dagger Gamma Psi` when `l`, `m`, `l'`, and `m'`
differ. Applying `GFAC_l` to each arbitrary angular product would be an
implementation choice, not a consequence of the current interface. Thus the
exact scalar-relativistic bilinear remains **BLOCKED**, while the Pauli
projection below is formally feasible. [BASIS MAPPING] [DEFERRED]

The missing exact-SR seam is precise: expose, or prove a reconstruction of, the
effective scalar-relativistic density/spin operator, including its
minor-component angular action and phase for all supported orbital channels.
The existing `full_scalar_relativistic_spin_angular` guard remains an exact-SR
guard and returns false for the current production path; it is not a rejection
of the separately defined Pauli projection. [DEFERRED]

## Complex angular convention and Gaunt mapping

`hcpx` documents the live complex order as `(0,0),(1,-1),(1,0),(1,1),...` and
contains the explicit `sp`/`spd` transformation. Reconciling those columns
with the standard Condon--Shortley functions gives

\[
 Y^{\mathrm{code}}_{lm}=\bigl(Y^{\mathrm{CS}}_{lm}\bigr)^*.
\]

The production-facing angular utility uses that convention and order.
[BASIS MAPPING]

For the product required by the transition field,

\[
 (Y^{\mathrm{code}}_{lm})^*Y^{\mathrm{code}}_{l'm'}
 =\sum_{LM}{\cal G}^{LM}_{lm,l'm'}Y^{\mathrm{code}}_{LM},
\]

with

\[
 {\cal G}^{LM}_{lm,l'm'}
 =(-1)^{m'}\sqrt{\frac{(2L+1)(2l+1)(2l'+1)}{4\pi}}
 \begin{pmatrix}L&l&l'\\0&0&0\end{pmatrix}
 \begin{pmatrix}L&l&l'\\M&m&-m'\end{pmatrix}.
\]

The phase is the transformed Condon--Shortley Gaunt phase, not an imported
real-harmonic convention. [BASIS MAPPING]

The nonzero channels obey

\[
 |l-l'|\le L\le l+l',
 \qquad M=m'-m,
 \qquad L+l+l'\ \text{even},
\]

along with `|m|<=l`, `|m'|<=l'`, and `|M|<=L`; the last parity condition is
the zero-`m` 3-j selection rule. [BASIS MAPPING]

The complete product of a one-electron basis truncated at `lmax` therefore
requires

\[
 L_{\max}^{\mathrm{response}}=2l_{\max},
\]

unless symmetry removes channels. Thus `sp` requires response `Lmax=2` and
the certified `spd` basis requires response `Lmax=4`. This is a product-space
consequence, not a convergence or performance choice. [BASIS MAPPING]

`source/response_angular_basis.f90` owns only this convention, the 3-j/Gaunt
evaluation, the product coefficients, and the cutoff. It does not construct a
transition vertex or response kernel.

## Independent angular and L0 oracles

`UnitLrResponseBasis` evaluates independently implemented complex harmonics on
a Gauss--Legendre-in-cos(theta) by uniform-phi sphere, then compares the
analytic coefficients with direct angular quadrature and reconstructs the
product pointwise. It covers `s x s`, `s x p`, diagonal and off-diagonal `p x
p`, diagonal and off-diagonal `p x d`, and diagonal and off-diagonal `d x d`.
[NUMERICAL REPRESENTATION]

The same test includes an algebraic charge/spin-z `L=0` normalization oracle:
for a normalized orbital harmonic, the projected value is

\[
 \int d\Omega\,(Y^{\mathrm{code}}_{00})^*|Y^{\mathrm{code}}_{lm}|^2
 =\frac{1}{\sqrt{4\pi}}.
\]

For explicit up/down weights `0.75` and `0.25`, the charge projection is the
above value and the `sigma_z` projection is one half of it. This is an
algebraic angular oracle, not the physical converged-state ground-state
oracle requested below. [BASIS MAPPING]

## Radial representation and measures

The LR-01 production mesh is retained as the reference direct discretization.
Its logarithmic coordinate is

\[
 r_i=B[\exp(A(i-1))-1],
 \qquad \frac{dr}{di}=A(r_i+B).
\]

The physical volume measure is `r_i^2 (dr/di) di dOmega`, whereas the legacy
weighted density is

\[
 \mathrm{RHO}_\sigma(r_i)=4\pi r_i^2n_\sigma(r_i).
\]

Consequently, a physical radial volume integral uses `r_i^2` and the mesh
Jacobian, while an electron-number integral of `RHO` uses the Jacobian but no
second `r_i^2` because that factor is already inside `RHO`. [NUMERICAL REPRESENTATION]

For an odd-point Simpson mesh, `w_1=w_N=1/3`, even interior points carry
`4/3`, and odd interior points carry `2/3`, so a direct radial integral is
`sum_i w_i A(r_i+B) f(r_i)`. A pointwise transition amplitude `T(r_i)` itself
does not absorb this radial measure; the measure enters a later radial
contraction of response fields. [NUMERICAL REPRESENTATION]

`source/response_basis_mapping.f90` exposes separate
`response_log_mesh_jacobian`, `response_volume_measure`, weighted-density
conversion, physical-density conversion, and Simpson integration helpers. The
origin conversion is guarded rather than evaluating a `0/0` expression. A
finite mesh and Simpson rule are numerical discretizations; they are not
claimed to be an exact continuum representation. [NUMERICAL REPRESENTATION]

No reduced radial basis `p_alpha(r)` is selected. A reduced basis would require
the explicitly measured projection

\[
 T_\alpha=\sum_i w_i\,p_\alpha^*(r_i)T(r_i)
\]

with the established physical radial measure and a convergence bound back to
the direct mesh. Since the Pauli `T` is only formally mapped and the exact-SR
operator remains blocked, no reduced basis is introduced here; such a
projection would add an unproved model-space choice and is deferred.
[DEFERRED]

## Formal Pauli transition-amplitude mapping

For the BES/Halle-style Pauli/no-SOC target, define the Pauli projection of the
certified scalar-relativistic LMTO state by retaining the large component that
`RSEQSR` identifies with `(U/r)Y_lm`:

\[
 \Psi^P_{n\mathbf k,a}(\mathbf r)=\frac1r
 \sum_{lm\sigma}c^{n\mathbf k}_{alm\sigma}
 U^{n\mathbf k}_{al\sigma}(r)Y^{\mathrm{code}}_{lm}(\hat r)|\sigma\rangle.
\]

For the certified first-order augmentation, the radial factor is

\[
 U^{n\mathbf k}_{al\sigma}(r)=U_{al\sigma}(r)+
 (\epsilon_{n\mathbf k}-E^{\mathrm{work}}_{\nu,al\sigma})
 \dot U_{al\sigma}(r).
\]

The right-hand sides are available from the LMTO coefficients and the stored
large radial augmentation arrays. This defines a Pauli-projected response
state, not the exact scalar-relativistic state. [BASIS MAPPING]

For an arbitrary explicitly supplied Pauli-space operator `Gamma^mu`, define

\[
 T^{P,\mu}_{nm;a\Lambda}(r;\mathbf k,\mathbf q)=
 \int d\Omega\,(Y^{\mathrm{code}}_\Lambda)^*(\hat r)
 \Psi^{P\dagger}_{n\mathbf k,a}(\mathbf r)\Gamma^\mu
 \Psi^P_{m,\mathbf k+\mathbf q,a}(\mathbf r).
\]

Substitution gives the evaluable formal mapping

\[
\begin{aligned}
 T^{P,\mu}_{nm;aLM}(r)=\frac1{r^2}
 \sum_{\substack{lm,l'm'\\\sigma\sigma'}}
 &(c^{n\mathbf k}_{alm\sigma})^*
 c^{m,\mathbf k+\mathbf q}_{al'm'\sigma'}
 U^{n\mathbf k}_{al\sigma}(r)
 U^{m,\mathbf k+\mathbf q}_{al'\sigma'}(r)\\
 &\times\Gamma^\mu_{\sigma\sigma'}
 \mathcal G^{LM}_{lm,l'm'}.
\end{aligned}
\]

At a direct radial point `r_i`, this is the formal response vector
`T^P_mu_nm;aLMi(k,q)`. Every factor is defined by the certified coefficient,
large-component radial, explicit Pauli-operator, and angular contracts; no
lower-component scalar-relativistic object is required for this target.
[BASIS MAPPING]

For the known large-component convention, the pointwise product carries
`U_l(r_i)U_l'(r_i)/r_i^2`; a later physical-volume contraction supplies
`r_i^2 dr`, leaving the mesh Jacobian. This is a controlled projection choice,
not permission to insert `GFAC` into arbitrary Gaunt products. [BASIS MAPPING]

No production `transition_density_radial` or `T` evaluator is added in LR-02R;
the formal mapping closes at the audit level. A future implementation must name
the resulting response as Pauli/no-SOC and must not present it as the exact
scalar-relativistic density response. [NUMERICAL REPRESENTATION] [DEFERRED]

## Canonical response super-index

The direct-mesh storage type is
`response_super_index(site,response_l,response_m,radial_point,channel)`.
The canonical order is site first, then response `L=0..Lmax` with
`M=-L..L`, then radial point `i=1..Nr`, then explicit channel
`mu=1..Nchannel`. Storage is complex because the response harmonics and
transition amplitudes are complex; no Hermitian reduction is silently applied.
[NUMERICAL REPRESENTATION]

Using `h(L,M)=L^2+L+M` as the zero-based harmonic slot, the one-based flat
index is

\[
 I=\bigl((((a-1)(L_{\max}+1)^2+h(L,M))N_r+(i-1))N_\mu+\mu\bigr).
\]

The inverse mapping first extracts `mu` and `i` by remainder/division, then
extracts `h` and recovers `L=floor(sqrt(h+1))` with `M=h-L^2-L` in the
implemented bounded search. Both directions are tested for all entries of a
two-site `Lmax=4`, seven-point, four-channel layout. [NUMERICAL REPRESENTATION]

## Finite-q Fourier phase and endpoint gauge

The reciprocal assembler folds fractional `k` into `[-1/2,1/2)` and uses

\[
 \exp(+i2\pi\mathbf k\cdot\mathbf d).
\]

The neighbor vectors are formed from endpoint coordinate differences. For an
off-site block this has the endpoint form
`d=R+tau_b-tau_a` in the primitive-cell interpretation. [BASIS MAPPING]

For a reciprocal-lattice shift `G`, the assembled matrix therefore obeys

\[
 H(\mathbf k+\mathbf G)=D(\mathbf G)H(\mathbf k)D(\mathbf G)^\dagger,
 \qquad D_a(\mathbf G)=\exp(-i2\pi\mathbf G\cdot\tau_a),
\]

and a compatible coefficient representative transforms as

\[
 c_a(\mathbf k+\mathbf G)=D_a(\mathbf G)c_a(\mathbf k).
\]

The sign follows from the live positive Fourier phase and endpoint displacement;
it is not guessed. [BASIS MAPPING]

`UnitLrResponseBasis` contains a two-site nontrivial-`tau` fixture and checks
the coefficient phase and matrix covariance. This closes the algebraic
endpoint-gauge relation. Combined with the formal Pauli amplitude, this is
enough to establish finite-q Pauli mapping feasibility at the contract level;
the production `T^P` evaluator and an end-to-end folded-eigenvector test remain
unimplemented. Exact scalar-relativistic finite-q response remains blocked with
the exact-SR operator. [NUMERICAL REPRESENTATION] [DEFERRED]

## Ground-state closure status

The exact-SR ground-state oracle remains the occupied diagonal sum

\[
 \sum_{\mathbf k n}^{occ}w_{\mathbf k}
 \Psi^\dagger_{n\mathbf k}(\mathbf r)\Gamma^\mu
 \Psi_{n\mathbf k}(\mathbf r),
\]

with its `L=0` charge component compared with accepted `RHO`/`n_up+n_down`
and its `sigma_z` component compared with `n_up-n_down` from the converged
bcc-Fe LR-01 fixture. [LITERATURE]

For the formal Pauli route, the corresponding occupied diagonal sum defines
`n^P_up(r)` and `n^P_down(r)` from the diagonal `T^P` amplitudes, with

\[
 m^P(r)=n^P_\uparrow(r)-n^P_\downarrow(r),
 \qquad
 \delta n_\sigma(r)=n^{SR}_\sigma(r)-n^P_\sigma(r),
 \qquad
 \delta m(r)=m^{SR}(r)-m^P(r).
\]

This comparison is the narrow numerical closure required before a production
`chiKS` implementation: report the pointwise differences, integrated moment
difference, and relative differences in the magnetically important radial
region for the accepted Fe fixture, preferably with Ni as a second case. No
arbitrary tolerance is prescribed before observing the scale. [BASIS MAPPING]

The current LR-01 snapshot stores accepted `n_up/down` and `B_xc^Pauli`, while
the new LR-02 work provides the angular algebra but does not yet accumulate
occupied Pauli-projected radial states from a converged reciprocal eigensystem.
Therefore this is **PENDING NUMERICAL CLOSURE**, not a synthetic pass and not a
claim that `n^P=n^SR`. [DEFERRED]

## Backend compatibility

The literature spectral structure is schematically

\[
 \chi^{0,\mu\nu}_{IJ}
 =\frac1{N_k}\sum_{\mathbf k,nm}
 \frac{f_{n\mathbf k}-f_{m,\mathbf k+\mathbf q}}
 {\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta}
 T^\mu_{nm;I}(T^\nu_{nm;J})^*.
\]

It is recorded only to show why all routes must eventually provide the same
`I=(a,L,M,i,mu)` transition vector; the sum is not implemented.
[LITERATURE]

| route | assessment | reason |
| --- | --- | --- |
| reciprocal eigenpair Lehmann | FORMALLY FEASIBLE for Pauli/no-SOC; numerically unvalidated | `epsilon,c`, large-component augmentation, explicit Pauli operator, and Gaunt mapping provide `T^P`; SR-versus-Pauli ground-state closure is pending |
| reciprocal coefficient-space GF | FORMALLY FEASIBLE for Pauli/no-SOC; numerically unvalidated | LR-GF-01 supplies coefficient-space Green-function contracts; the same Pauli augmentation/angular vertex would be applied on both sides |
| native real-space/block-recursion GF | DEPENDS ON LR-REP-00 | endpoint/orbital representation transformations and augmentation into the same radial/angular space are not certified |

No backend is promoted to READY, and no auxiliary/screened-GF equality is
assumed. The formal Pauli mapping does not authorize coding the bubble or
claiming exact scalar-relativistic response.

## Relation to LR-01 XC fields

LR-01's accepted `m(r)` and `B_xc^Pauli(r)` live on the same production radial
mesh as the proposed direct response index. A later local adiabatic kernel may
therefore act pointwise in that radial coordinate, rather than after reducing
the response to one number per site. For the Pauli projection, however, the
scalar-relativistic ground-state fields and Pauli-projected response are not
automatically the same functional representation; that consistency must be
tested by the SR-versus-Pauli closure above. [LITERATURE] [DEFERRED]

A later ratio such as

\[
 K_{xc}(r)\sim B_{xc}(r)/m(r)
\]

has a domain issue at nodes or near-zero `m(r)`. No regularization is invented
in LR-02; the ratio domain and policy belong to the later XC-kernel task.
[DEFERRED]

## Evidence and claim levels

| evidence | claim level | result |
| --- | --- | --- |
| `UnitLrResponseBasis` independent harmonic quadrature | independent numerical cross-check | passes for representative diagonal/off-diagonal `m,m'` products |
| cutoff, radial measure, weighted-density conversion, super-index, gauge, and capability guard | algebraic consistency | passes |
| L0 charge/spin-z oracle in the new unit test | algebraic consistency only | passes; not a converged-state validation |
| LR-01 accepted bcc-Fe occupied-state SR-versus-Pauli reconstruction | physical converged-state validation | pending; required to quantify the controlled projection error |
| finite-q Pauli transition-density covariance | formal representation mapping | formally feasible from the coefficient gauge and `T^P`; end-to-end production test is not implemented |

The new unit test was run with GNU Fortran 13.3, OpenMP enabled, MPI disabled,
and the repository's configured libXC. It reports `UnitLrResponseBasis: PASS`
for the bounded algebraic/independent oracles while explicitly retaining the
exact-SR capability guard. The formal Pauli mapping is documented but not
implemented as a production evaluator.

## Capability matrix and checklist

| capability | status |
| --- | --- |
| `sp`, `spd` complex harmonic convention and ordering | PASS for algebraic mapping |
| complete product cutoff `Lresponse=2*lmax` | PASS |
| direct production radial mesh and measures | PASS as numerical representation |
| canonical `(site,L,M,radial point,channel)` index | PASS |
| explicit arbitrary Pauli `2x2` operator slot | PASS at formal mapping level; not physically evaluated in code |
| scalar-relativistic nonspherical local bilinear | BLOCKED |
| Pauli/no-SOC `T^P_mu_nm;aLMi(k,q)` | FORMALLY FEASIBLE; production evaluator deferred |
| exact-SR `T^mu_nm;aLMi(k,q)` | BLOCKED |
| converged SR-versus-Pauli charge/spin `L=0` reconstruction | PENDING NUMERICAL CLOSURE |
| endpoint phase algebra | PASS as algebraic fixture |
| finite-q Pauli transition-density gauge | FORMALLY FEASIBLE; end-to-end test deferred |
| exact-SR finite-q transition-density gauge | BLOCKED |
| reciprocal eigenpair route | FORMALLY FEASIBLE for Pauli/no-SOC; unvalidated |
| reciprocal coefficient-GF route | FORMALLY FEASIBLE for Pauli/no-SOC; unvalidated |
| native real-space GF route | DEPENDS ON LR-REP-00 |
| noncollinear, SOC, generalized overlap, additive operators, `spdf` | DEFERRED |

Checklist status is therefore: prerequisites verified; previous audits reused;
angular convention/Gaunt/cutoff/mesh/index/gauge primitives added and tested;
formal Pauli/no-SOC transition mapping closed; SR-versus-Pauli ground-state
closure and an end-to-end production `T^P` evaluator remain pending. Exact
scalar-relativistic response remains blocked. No site-only substitute,
susceptibility, or XC kernel was added.

The formal mapping can be recorded as the LR-02R **PARTIAL PASS**. Before
enabling production `chiKS`, the next narrow task is to measure the SR-to-Pauli
projection error on accepted Fe data, preferably with Ni as a second case; an
exact-SR response would instead require the separate effective density/spin
operator seam described above.
