# TDDFT Ward-repair blueprint

## Status and scope

This document freezes the convention contract for the Ward-identity repair.
It applies only to the CPU, collinear, no-SOC transverse response with
`reciprocal_mode='ham_only'`. It does not change the current production
spectra: the existing site-scalar kernel remains a legacy comparison route
until the pair-potential construction and its LMTO representation oracle are
implemented and validated.

The current `P_site sigma` coefficient-space vertices are not physical
operators for a generalized eigenproblem `H c = epsilon S c` without a
derived metric representation (and, where relevant, a `delta S` term).
Accordingly, `generalized_overlap_proxy` and `generalized_overlap_kanpur` are
outside this repair slice and must be rejected by the future pair-potential
route rather than treated as supported.

## Frozen transverse convention

The response coordinates and vertices are unhalved:

\[
m^\pm=m_x\pm i m_y,\qquad
\sigma^\pm=\sigma_x\pm i\sigma_y.
\]

Thus `RESPONSE_PLUS` and `RESPONSE_MINUS` have spin-flip matrix element two.
The magnetic Hamiltonian is

\[
H_{\rm mag}=B_x\sigma_x+B_y\sigma_y+B_z\sigma_z
=\tfrac12(B^-\sigma^+ + B^+\sigma^-)+B_z\sigma_z.
\]

For a collinear ground state along `z`, LSDA/ALSDA rotational invariance gives

\[
\delta B^+(\mathbf r)=\frac{B_{\rm xc}(\mathbf r)}{m_z(\mathbf r)}
\delta m^+(\mathbf r).
\]

The corresponding circular transverse kernel is therefore

\[
\kappa^{-+}(\mathbf r)=\frac{B_{\rm xc}(\mathbf r)}{2m_z(\mathbf r)}.
\]

The factor one half belongs to the unhalved circular convention. It must not
be copied to Cartesian coordinates: the Cartesian transverse derivative is
`Bxc/mz`, twice the circular scalar.

## Magnetization-shape coordinate and pair potential

For site `a`, use the converged signed collinear magnetization shape

\[
g_a(\mathbf r)=\begin{cases}m_z(\mathbf r)/M_a,&\mathbf r\in\Omega_a,\\
0,&\text{otherwise},\end{cases}\qquad
M_a=\int_{\Omega_a}m_z(\mathbf r)\,d\mathbf r.
\]

An input coordinate has

\[
\delta m^+(\mathbf r)=\sum_a g_a(\mathbf r)\delta M_a^+.
\]

Its kernel-weighted, or pair-potential, operator is

\[
\hat Q_a^-=g_a(\mathbf r)\kappa^{-+}(\mathbf r)\sigma^-
=\frac{B_{\rm xc}(\mathbf r)}{2M_a}\sigma^-\,\mathbf1_{\Omega_a}.
\]

`Bxc(r)` must remain an orbital/radial operator until matrix elements in the
active LMTO representation are evaluated. Replacing it with
`K_site * P_site sigma_minus` is only equivalent in a deliberately uniform
two-level oracle, not in a general material.

## Direct self-enhancement operator

Let

\[
V^+_{a,nm}=\langle n\mathbf k|P_a\sigma^+|m\mathbf k+\mathbf q\rangle,
\qquad
W^-_{b,mn}=\langle m\mathbf k+\mathbf q|Q_b^-|n\mathbf k\rangle.
\]

The direct self-enhancement operator is

\[
\Xi_{ab}(\mathbf q,\omega)=\sum_{\mathbf k,n,m}w_{\mathbf k}
\frac{f_{n\mathbf k}-f_{m\mathbf k+\mathbf q}}
{\omega+\epsilon_{n\mathbf k}-\epsilon_{m\mathbf k+\mathbf q}+i\eta}
V^+_{a,nm}W^-_{b,mn}.
\]

The observable bare response retains its ordinary right vertex `V^-`; the
dressed response is evaluated in the correct ordering as

\[
\chi=(I-\Xi)^{-1}\chi_{\rm KS}.
\]

No separately truncated `K` is inferred by inverting `chi_KS`.

## Raw Ward identity

For a rigid infinitesimal rotation, `delta M_a^+=M_a theta`, and therefore

\[
\sum_a M_a Q_a^-=\tfrac12B_{\rm xc}(\mathbf r)\sigma^-.
\]

When the same self-consistent spin-dependent operator is represented in the
eigensolver and response, the uncorrected static result obeys

\[
\Xi(\mathbf0,0)\mathbf M=\mathbf M,\qquad
r_G=\frac{\|\Xi\mathbf M-\mathbf M\|_2}{\|\mathbf M\|_2}=0.
\]

This raw identity is a release gate. A later controlled correction may be
reported beside it, but can never substitute for it.

## WR-00 executable evidence

`UnitTddftWardConventions` uses analytic two-level and two-orbital oracles.
It pins the unhalved circular matrices, the half in `Q^-`, Cartesian/circular
Dyson equivalence, the uniform-kernel reduction `Xi=chi_KS*K`, and the
unequal-exchange negative control. For the latter, the exact weighted
two-orbital vertex satisfies the Ward identity to numerical precision while
the old site-averaged scalar kernel has a nonzero residual.

## WR-01 generic direct-Xi layer

`tddft_xi_mod` accepts a set of dense weighted operators in the same
coefficient representation as the eigenvectors. It evaluates
`<m,k+q|Q_b|n,k>` without assuming Hermiticity or applying a conjugation, then
forms direct `Xi` with the same finite-temperature occupations, retarded
denominators, band window, k weights, and scalar/batched controls as
`chi_KS`. This layer deliberately knows nothing about radial XC data or the
LMTO construction of `Q`; WR-02 supplies that representation.

`UnitTddftDirectXi` establishes the uniform scalar reduction only as an
oracle, the unequal-orbital divergence from the old scalar route, scalar versus
batched equality for complex `q != 0` spinors, and the right-vertex-order
negative control. It does not alter `calculation.f90` dispatch.

## WR-01b endpoint-resolved LMTO magnetic tangents

`ham0m_nc` now delegates ordinary bond assembly to the side-effect-free
`lmto_magnetic_tangent_mod`. The same coefficient construction supplies the
fixed-potential derivative before `chbar_nc` collapses it into `hmag`, `hxc`,
and `ee`. This is provenance for a later vertex provider, not a TDDFT dispatch
change.

For a directed bond `i -> j`, with `S = hhh`, the actual implemented orbital
matrices are

\[
A_{00}=w_{0i}S w_{0j},\quad A_{11}=w_{1i}S w_{1j},\quad
A_{10}=w_{1i}S w_{0j},\quad A_{01}=w_{0i}S w_{1j}.
\]

With `e_i=potential%mom_i`, `e_j=potential%mom_j`, and onsite `D_i=cx1_i`
(`cex1_i` in the HOH value path), the core returns

\[
H_0=A_{00}+A_{11}(e_i\!\cdot e_j)+C_{0i},
\]
\[
H_\alpha=A_{10}e_{i\alpha}+A_{01}e_{j\alpha}
+iA_{11}(e_i\!\times e_j)_\alpha+D_i e_{i\alpha}.
\]

For explicit endpoint tangents `d_i=delta_mom_i` and `d_j=delta_mom_j`, the
new narrow API `ham0m_nc_tangent` returns

\[
\delta H_0=A_{11}(d_i\!\cdot e_j+e_i\!\cdot d_j),
\]
\[
\delta H_\alpha=A_{10}d_{i\alpha}+A_{01}d_{j\alpha}
+iA_{11}(d_i\!\times e_j+e_i\!\times d_j)_\alpha+D_i d_{i\alpha}.
\]

`potential%mom` is a **unit direction**, not a moment magnitude: it is
normalised in `potential.f90` at construction and input boundaries, and the
Hamiltonian uses it only as the orientation multiplying fixed spin-dependent
LMTO potential parameters. The response population `M_a` therefore remains a
separate, signed-normalisation problem for WR-02; WR-01b does not assume that
`potential%mom` supplies it.

The cartesian-to-spinor helper uses the existing convention

\[
H=H_0 I+H_x\sigma_x+H_y\sigma_y+H_z\sigma_z,
\]

so its blocks are `(H0+Hz, Hx-iHy; Hx+iHy, H0-Hz)`. WR-02 will form
`Q_a^-=(D_{a,x}-iD_{a,y})/(2M_a)` only after it assigns endpoint phases in the
reciprocal convention and establishes signed `M_a` provenance.

`ham0m_nc_endpoint_tangents` is an on-demand provider. Its record carries
source/neighbor **site** identities, their types, the directed neighbor/Fourier
bond slot, displacement, onsite ownership, and operator generation. Thus two
sites of type `1` do not alias. It returns left and right tangents with an
explicit final Cartesian-component index, so a future response layer obtains
all three components without trying to split an assembled block.

| Feature | WR-01b status |
| --- | --- |
| Ordinary CPU `periodic_nc` first-order LMTO, fixed potential | supported |
| `explicit_texture` endpoint lookup | supported only when its normal bulk-reuse guard holds |
| SOC, constraining/FSM-field provenance | rejected pending a separated derivation |
| HOH/generalized overlap | rejected: `delta S`/HOH tangent absent |
| GBT bond gauge | rejected |
| two-centre CCOR | rejected |
| Hubbard U/V and local-axis transformed route | rejected |

The capability result is false with a specific reason for every rejected
active route; it returns zero storage rather than a partial tangent. The core
unit oracle checks value preservation, complex spinor mapping, left/right and
common-rotation central differences, dot/cross negative controls, reverse-bond
covariance, same-type distinct-site records, and each capability rejection.

### WR-02 handoff

Consume `ham0m_nc_endpoint_tangents` before `hxc`/`ee` assembly, attach the
phase of the *perturbed endpoint* in the reciprocal Fourier convention, then
form `D_x`, `D_y`, and `Q^-`. Validate the q-dependent operator against a
commensurate-supercell finite rotation. No radial reconstruction or scalar
`K_site P_site sigma_minus` fallback is permitted.
