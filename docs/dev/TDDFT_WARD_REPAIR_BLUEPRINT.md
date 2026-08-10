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

## WR-02 LMTO representation audit — blocked

The required active-basis pair potential cannot be derived faithfully from the
stored ground-state state on commit `d3f5b1f`. This is a hard blocker, not a
license to use `K_site P_site sigma_minus`.

The eventual `Q_a` is a full real-space LMTO block operator, not an onsite
orbital identity. `ham0m_nc` forms each directed bond from separate
source-endpoint, target-endpoint, and cross-product magnetic terms; the normal
`build_bulkham` route then sums those terms into `hmag`, packs the spin blocks
into `hxc`, and Fourier transforms `ee`. Thus tails/hopping and any
representation-generated matrix elements must remain in `Q_a`.

Three facts prevent an analytic local derivative today:

1. `hmag` is a scratch array reset by every `chbar_nc` call. The persistent
   `hxc`/`ee` blocks contain only the **sum** of the two endpoint contributions
   to a bond. A derivative with respect to rotation of site `a` needs the
   source and target terms separately; they cannot be reconstructed from their
   sum.
2. `xc_response_kernel_provider` keeps only integrated radial quantities
   (`radial_spin_population`, `bxc_spin_moment`, and scalar averages), rather
   than radial samples or active-basis matrix elements of `Bxc(r)`. It cannot
   independently supply the missing orbital matrix.
3. The current provider/TDDFT driver stores a magnitude-only site population
   (`mtot` and, in the response driver, `sqrt(mx^2+my^2+mz^2)`). It therefore
   cannot provide the signed `M_a` normalization required for reversed or
   antiferromagnetic sublattices. The constraining-field contribution is also
   merged into the radial potential without an immutable response provenance.

The smallest prerequisite is a dedicated Hamiltonian derivative API. It must
persist or return the analytic source-endpoint and target-endpoint derivatives
of `ham0m_nc` before they are summed, assemble their full real-space and
Fourier-transformed matrices in the active `ham_only` basis, retain signed
site moments and XC-versus-constraining-field provenance, and expose a
side-effect-free `+theta/-theta` rebuild for the independent oracle. Only then
can WR-02 compare analytic `Q_a^+=(Q_a^-)^dagger` against central rotations,
including a rigid all-site rotation. Generalized-overlap modes remain rejected
until a metric-correct operator and any `delta S` contribution are derived.
