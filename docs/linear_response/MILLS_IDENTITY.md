# Mills-1U — controlled RS-LMTO model projection

**Classification:** `MILLS-1U — CONTROLLED RS-LMTO MODEL PROJECTION`

This route defines a single local d-shell interaction from native accepted
RS-LMTO quantities. The scalar interaction is not extracted from a fit to the
complete Hamiltonian spin difference. The bcc-Fe build, unit, and material
gates passed; the accepted-state measurements are recorded below.

## Literature and model definition

The one-parameter Mills mean-field Hamiltonian is

\[
H_i^C=-\frac{U_i m_i^d}{2}\sum_{\mu\in d}
  (n_{i\mu\uparrow}-n_{i\mu\downarrow}),
\qquad
m_i^d=\sum_{\mu\in d}(n_{i\mu\uparrow}-n_{i\mu\downarrow}).
\]

The spin-up potential is lowered and the spin-down potential raised for
positive \(m_i^d\). Therefore

\[
\Delta_i^d=H_{i,d,\downarrow}^C-H_{i,d,\uparrow}^C=U_i m_i^d,
\qquad
U_i^{\rm Mills}=\frac{\Delta_i^d}{m_i^d}.
\]

All five local d orbitals receive the same scalar field. The interaction has
energy units; the implementation reports Ry. The coefficient-space moment is
reported in \(\mu_B\), numerically equal to the dimensionless occupation
difference for \(g=2\).

## RS-LMTO mapping

For each accepted response site, the live transformed potential is
`potential%center_band(3,spin)`, with spin index 1 up and 2 down. The route
uses

\[
\boxed{\Delta_i^d=C_{i,d,\downarrow}-C_{i,d,\uparrow}},
\qquad
\boxed{U_i^{\rm Mills}=\Delta_i^d/m_i^d}.
\]

These values are read from the final accepted potential paired with the
reciprocal eigenpairs consumed by the response. They are not reconstructed
from an earlier SCF iteration or replaced with raw atomic parameters. The
spherical-ASA identity is checked against all five orbital-resolved `cx`
entries produced by `symbolic_atom%build_pot`.

The moment is the accepted orthogonal coefficient-space quantity

\[
m_i^d=\mu_B\sum_{kn}w_k f_{nk}\sum_{\mu\in d}
  (|c^{nk}_{i\mu\uparrow}|^2-|c^{nk}_{i\mu\downarrow}|^2).
\]

It uses the accepted reciprocal eigenvectors, occupations, k weights, Fermi
level, and temperature. It selects only \(l=2\), keeps the full spd
Hamiltonian eigenstates, and contains no radial norm, radial overlap, or
fitted-moment factor. A vanishing or non-finite moment is rejected.

## Controlled approximation and residual

The accepted RS-LMTO Hamiltonian has spin dependence beyond the local band
center. Width/hopping parameters contribute to \(H_\downarrow(k)-H_\uparrow(k)\),
and second-order terms contain the spin-dependent \(-HOH\) correction. The
route retains the native onsite d-center splitting and measures

\[
D(k)=H_\downarrow(k)-H_\uparrow(k),\qquad
D_{\rm Mills}=\sum_i\Delta_i^dP_i^d,\qquad
R(k)=D(k)-D_{\rm Mills}.
\]

No part of \(R(k)\) is fitted into \(U\). The output reports Frobenius norms
for the retained local d term, local d-shell remainder, remaining d sector,
non-d sector, k-dependent spin difference, unit-cell site-offdiagonal spin
dependence, total residual, and the second-order \(-HOH\) spin contribution.
The local d-shell remainder is a subset of the remaining d-sector norm; the
site-offdiagonal norm overlaps the orbital-sector norms. The k-dependent norm
is the weighted RMS of \(D(k)-\langle D\rangle_k\), so it captures periodic
hopping in a one-site unit cell, where the unit-cell site-offdiagonal block is
empty. The first-order and second-order Hamiltonians are evaluated through the
existing reciprocal Fourier builders; the explicit onsite \(e_\nu\) difference
is removed to isolate the HOH part.

The relative residual is provenance for this controlled model projection. It
is not a criterion that changes \(U\). The accepted state is described as
close to or far from the one-U Mills model only after the material values are
measured and reported.

## Accepted bcc-Fe gate

The `LinearResponseMills1U` integration gate passed with build version
`lr-campaign-archive-53-ge31d-dirty`. It converged one 4x4x4 accepted SCF state
in 19 iterations (SCF residual \(9.23\times10^{-7}\)), then handed the same
accepted reciprocal eigenpairs directly to the response. The state artifact
reports 64 k points, unit total k weight, 300 K, Fermi level
\(-0.0622448565095445\) Ry, and accepted total moment
\(2.00000130932976\,\mu_B\). The response reads the final transformed
potential from that same run; it does not reload a prior potential or re-SCF.

| Accepted bcc-Fe quantity | Value |
|---|---:|
| \(C_{d,\uparrow}\) | \(-0.204493711938692\) Ry |
| \(C_{d,\downarrow}\) | \(-0.0624298456361092\) Ry |
| \(\Delta_d=C_{d,\downarrow}-C_{d,\uparrow}\) | \(0.142063866302583\) Ry |
| Coefficient-space \(m_d\) | \(2.00785485887215\,\mu_B\) |
| Physical \(U_d=\Delta_d/m_d\) | \(0.0707540516062910\) Ry (\(0.962657913\) eV) |
| Five-channel spherical-center residual | 0 Ry |

For the accepted second-order Hamiltonian, weighted Frobenius RMS values were
\(\|D\|=0.361567049\) Ry, \(\|D_{\rm Mills}\|=0.317664462\) Ry, and
\(\|R\|=0.151284631\) Ry, giving \(\|R\|/\|D\|=0.418414\) and
\(\|R\|/\|D_{\rm Mills}\|=0.476240\). The local d-shell remainder was
0.0417066790 Ry (equal to the remaining d-sector norm for this one-site
cell); the non-d sector was 0.145422118 Ry. The k-dependent part of \(D(k)\)
was 0.0967165943 Ry, while the unit-cell site-offdiagonal norm was zero for
the one-site basis. The isolated second-order \(-HOH\) spin contribution was
0.0464586287 Ry; the total second-order-minus-first-order spin difference was
0.124473857 Ry. These measurements show a substantial controlled-model
remainder; none was used to retune \(U_d\).

At \(q=0,\omega=0\), the raw moment-action-relative denominator decreased
across the requested positive eta ladder:

| Eta (Ry) | Minimum singular value | \(\|D(0,0)m_d\|/\|m_d\|\) |
|---:|---:|---:|
| 0.04 | 0.250109 | 0.250109 |
| 0.02 | 0.136336 | 0.136336 |
| 0.01 | 0.0828154 | 0.0828154 |
| 0.005 | 0.0622635 | 0.0622635 |
| 0.0025 | 0.0559302 | 0.0559302 |

This is the uncorrected finite-eta result; it is not zeroed or used to select
the interaction. At the smallest eta, q/-q covariance for the mesh-compatible
pair \(q=\pm0.25\) was \(3.06\times10^{-12}\). The noncommensurate
\(\pm0.03\) and \(\pm0.125\) pairs gave \(4.09\times10^2\) and
\(1.02\times10^3\) on this 4x4x4 quadrature: these finite-mesh sums are not
invariant under a q translation that does not map the sampled k mesh onto
itself. They are reported as mesh diagnostics, not covariance passes. For
example, at \(q=(0.25,0,0)\), \(\omega=0.10\) Ry, and \(\eta=0.01\) Ry, the
code-normalized bare response was \(-59.2418-18.8191i\) Ry\(^{-1}\), the
enhanced response was \(47.1079-11.4470i\) Ry\(^{-1}\), and the loss was
3.64371 Ry\(^{-1}\). Complete rows and accepted-state checksums are in
[`mills_1u_fe_4k.dat`](../../tests/integration/mills_1u_bccfe/mills_1u_fe_4k.dat)
and its `.state` sidecar.

## Bare transverse response

The local interaction vertex is

\[
S_i^+=\sum_{\mu\in d}c^\dagger_{i\mu\uparrow}c_{i\mu\downarrow},
\qquad
T_i^{nm}(k,q)=\sum_{\mu\in d}
 c^{nk*}_{i\mu\uparrow}c^{m,k+q}_{i\mu\downarrow}.
\]

The reverse channel exchanges up and down. There is no radial-overlap
augmentation. The site-contracted bare response is accumulated directly from
the accepted full spd eigenstates:

\[
\chi_{0,ij}^{+-}(q,\omega)=\sum_{knm}
\frac{2w_k(f_{nk}-f_{m,k+q})}{
 \omega+\epsilon_{nk}-\epsilon_{m,k+q}+i\eta}
T_i^{nm}(k,q)T_j^{nm*}(k,q).
\]

The factor 2 is the repository's circular-channel normalization. The vertex
itself is the coefficient-space matrix element of
\(\sigma^+=(\sigma_x+i\sigma_y)/2=S^+\); there is no factor two in \(T\).

## Dyson normalization and sign

With the retarded commutator convention
\(\chi_{AB}=-i\theta(t)\langle[A(t),B(0)]\rangle\), a one-site split level
with \(\Delta=H_\downarrow-H_\uparrow>0\) has
\(\chi_{0,\mathrm{Kubo}}^{+-}(0)=-m/\Delta\). The mean-field transverse
perturbation is \(H'=-U\,\delta\langle S^+\rangle S^-+\mathrm{h.c.}\), so
linear response gives

\[
\chi_{\mathrm{Kubo}}=\chi_{0,\mathrm{Kubo}}-
 \chi_{0,\mathrm{Kubo}}U\chi_{\mathrm{Kubo}},
\qquad
\chi_{\mathrm{Kubo}}=\frac{\chi_{0,\mathrm{Kubo}}}
 {1+U\chi_{0,\mathrm{Kubo}}}.
\]

The shared site solver forms \(I-\chi_{\rm code}K\). Since
\(\chi_{\rm code}=2\chi_{\mathrm{Kubo}}\), its input must be
\(K=-U/2\), giving

\[
I-\chi_{\rm code}K=I+U\chi_{0,\mathrm{Kubo}}.
\]

This sign and factor follow from the stated mean field, the retarded response,
\(m=n_\uparrow-n_\downarrow\), \(\Delta=H_\downarrow-H_\uparrow\), and the
existing channel and solver conventions. The one-site analytic oracle checks
the bubble and Dyson denominator directly and includes wrong-sign and
wrong-factor negative controls. No normalization is selected from a Goldstone
condition, a Jülich comparison, or an expected numerical interaction.

## Goldstone result

The raw \(q=0,\omega=0\) denominator and its action on the accepted d moment
are reported at the requested positive broadenings. There is no interaction
sum-rule solve, pole shift, denominator correction, or post-hoc zeroing. For
the exact one-site Mills mean field the zero-broadening identity follows from
\(\Delta=Um\). A finite residual for an accepted RS-LMTO state is an
acceptance result for the controlled projection and should be compared with
the measured \(R(k)\) residual.

## Supported scope

The route accepts scalar-relativistic, collinear, no-SOC, orthogonal
`ham_only`, second-order reciprocal states with the full spd basis. The local
transverse interaction is d-only and has one scalar U per site. SOC,
noncollinear states, generalized overlap, other basis sizes, radial response,
longitudinal response, orbital-dependent or full Coulomb interactions, and
Jülich sum-rule determination are not supported by this Mills-1U route.

The former fitted `stoner_fit` interaction and its site-vertex least-squares
scalarization are removed. The old keyword fails explicitly as unsupported.
