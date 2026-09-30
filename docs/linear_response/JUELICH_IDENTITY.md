# Jülich projected-site identity

## Status and scope

The certified reference implemented here is **Juelich-d**: the Lounis–Costa–
Muniz–Mills d-projected construction, represented in the scalar-relativistic
RS-LMTO basis. The KKR and LMTO radial functions are not identical. Their
mapping is recorded as a controlled approximation, so this route does not
claim exact KKR radial physics.

This note covers only the d reference and the algebraic generalization gate
for spd and spdf. It does not certify the separate Mills or spatial TDDFT
routes.

Primary source: S. Lounis, A. T. Costa, R. B. Muniz, and D. L. Mills,
*Phys. Rev. B* **83**, 035109 (2011),
[DOI](https://doi.org/10.1103/PhysRevB.83.035109),
[Jülich hosted full text](https://juser.fz-juelich.de/record/14113/files/PhysRevB.83.035109.pdf).
Equation numbers below refer to that paper.

## Published construction

The retarded transverse Kohn–Sham bubble is Eq. 2, with its static form in
Eq. 3. The spin-rotation sum rule is derived in Eqs. 6–10. Its angular
contraction gives the explicit \(4\pi\) in Eq. 12. With
\(B_{\mathrm{eff}}=V_\downarrow-V_\uparrow\), the paper defines the local
kernel in Eq. 13 and obtains the static sum rule \(\Gamma U=M_z\) in
Eqs. 14–17.

The projection is specified in Sec. V. Equations 36–39 replace the full,
energy-dependent KKR Green function with a spectral expansion whose radial
functions are fixed. Equation 40 projects the full Green function onto those
radial functions and divides by their norms. Equation 41 chooses
\(\phi^{iL}(r)=R^{id}(r;E_F)\): the d-character regular KKR solution at the
Fermi energy, normalized to \(\psi^{id}\). Equation 42 gives the projected
Green function. The band energies, occupations, and spectral denominators
remain energy dependent; the projector radial shape does not.

Equations 43–44 reduce the interacting response to site space after inserting
the d-projected susceptibility. Equation 45 gives a radial matrix element for
the local interaction and Eq. 46 its ALDA specialization. The sum-rule
alternative is Eq. 47:

\[
\Gamma_{ij}=\bar{\bar\chi}^{ij}_0(0)M_{z,j},\qquad
\sum_j\Gamma_{ij}U_j=M_{z,i}.
\]

The moment in Eq. 47 is the d-projected site moment from the Sec. V Green
function construction. Its spin-density convention is the spin difference
of Eq. 9; the code reports that moment in \(\mu_B\). The paper's interaction
therefore has units of energy per \(\mu_B\). The interaction used here comes
from Eq. 47, not from Eq. 46.

The enhanced response follows Eq. 43,
\(\chi=\chi_0+\chi_0U\chi\), solved as
\((I-\chi_0U)\chi=\chi_0\), with \(U\) diagonal in site space. The static
identity is the retarded \(\omega=0,\eta\to0^+\) limit. \(U\) is determined
from that limit. No post-construction Goldstone shift is applied.

## Juelich-d mapping to RS-LMTO

| Jülich quantity | Paper equation | RS-LMTO quantity | Mapping type | Residual uncertainty |
|---|---|---|---|---|
| Frozen d radial projector | Eqs. 36–41 | Spin-resolved \(\phi_{d}(E_F)=\phi_d+ (E_F-\epsilon_\nu)\dot\phi_d\), from the scalar-relativistic LMTO linearization | CONTROLLED APPROXIMATION | The KKR regular solution at \(E_F\) and this LMTO linearized partial wave are not identical. Their shape difference is the principal representation error. |
| Projector normalization | Eqs. 40–41 | Log-mesh Simpson integration of the large/small radial components with the RS-LMTO scalar-relativistic norm metric; each spin and site is normalized once at \(E_F\) | CONTROLLED APPROXIMATION | The normalization is exact for the selected LMTO radial metric and mesh. It approximates the KKR cell norm in the paper. |
| Angular normalization | Eqs. 12, 22, 40–44 | Unit-normalized complex \(Y_{2m}\), \(m=-2,\ldots,2\); the Eq. 22 angular contraction is the explicit \(4\pi\) multiplier | EXACT REPRESENTATION TRANSFORM | This relies on the stated unit-normalized harmonic convention. |
| Projected moment | Eq. 9 and Sec. V, Eqs. 36–42; used in Eq. 47 | \(\mu_B\sum_{nk}w_k f_{nk}\sum_m(|P^\uparrow_{nkm}|^2-|P^\downarrow_{nkm}|^2)\), using the same frozen normalized projectors as the response | CONTROLLED APPROXIMATION | The projector is LMTO-represented KKR \(R_d(E_F)\); occupied band states and Fermi occupations remain those of the accepted reciprocal state. |
| Projected bare susceptibility | Eqs. 2, 22, 40–44 | Spin-flip Lehmann sum of frozen-projector transition amplitudes, multiplied by \(4\pi\), with one spin-flip matrix element and retarded denominator | CONTROLLED APPROXIMATION | The spectral contraction is the representation of the paper's bubble; the d radial projector uses the LMTO approximation above. |
| Circular spin operator | Eqs. 2 and 22 | One \(+\) spin-flip transition amplitude \(\sum_m P^{\uparrow *}_{im}P^\downarrow_{jm}\), paired with its site-matrix conjugate | EXACT REPRESENTATION TRANSFORM | Only the collinear, spin-conserving, no-SOC domain is certified. |
| Factor of two | Eqs. 2 and 22 | No extra factor of two in the spin-flip Lehmann numerator. The angular contraction contributes \(4\pi\), not a second spin-flip copy. | EXACT REPRESENTATION TRANSFORM | The factor convention is fixed by Eq. 2 together with the Eq. 22 angular contraction. |
| Moment units and spin sign | Eqs. 9, 13, 47 | \(M_d=\mu_B(n^\uparrow_d-n^\downarrow_d)\); \(B_{\rm eff}=V_\downarrow-V_\uparrow\); \(U\) is in Ry/\(\mu_B\) | EXACT REPRESENTATION TRANSFORM | The resulting sign of \(U\) follows the paper's potential-difference convention and the computed \(\chi_0\); no sign is fitted. |
| \(\Gamma\) | Eq. 47 | Column \(j\) of the static site matrix is multiplied by projected moment \(M_{d,j}\) | EXACT | This is exact algebra conditional on the projected \(\chi_0\) and moment above. |
| Local \(U\) | Eq. 47 | Solve \(\Gamma U=M_d\) once using the extrapolated static real \(\chi_0\) | CONTROLLED APPROXIMATION | The physical definition is the retarded static limit. The finite eta ladder and extrapolation are numerical approximations and are gated by convergence and rank checks. |
| Enhanced susceptibility | Eq. 43 | Site-space Dyson solve \((I-\chi_0U)\chi=\chi_0\) | EXACT REPRESENTATION TRANSFORM | Exact for the paper's site-level ansatz and supplied projected inputs. |
| Goldstone/magnetic mode | Eqs. 10, 17, 43, 47 | Report the relative action of \(I-\chi_0(0)U\) on the projected moment vector | EXACT | The diagnostic is not used to alter \(U\) or \(\chi_0\). Numerical residual depends on the static-limit estimate. |

### Radial and angular definitions

For each site and spin, the LMTO projector stored by the implementation is

\[
\widetilde\phi_{d\sigma}(r;E_F)
=\phi_{d\sigma}(r)+
(E_F-\epsilon_{\nu,d\sigma})\dot\phi_{d\sigma}(r).
\]

The normalizing integral uses the code's logarithmic radial measure
\(dr=a(r+b)d\xi\) and scalar-relativistic norm
\[
\mathcal N_{d\sigma}=
\int dr\,[g_{\rm EF}(r)|\widetilde\phi^L_{d\sigma}(r)|^2+
|\widetilde\phi^S_{d\sigma}(r)|^2].
\]
The stored projector is \(\psi_{d\sigma}=\widetilde\phi_{d\sigma}/
\sqrt{\mathcal N_{d\sigma}}\). The five complex \(d\) harmonics each have
unit angular norm. The \(4\pi\) is applied once in \(\chi_0\) to represent
the paper's angular contraction; it is not included in the radial norm or
the projected moment.

The same \(\psi_{d\uparrow},\psi_{d\downarrow}\) define both occupied-state
projection amplitudes and spin-flip transition amplitudes. Band radial
states still use their energy-dependent LMTO reconstruction
\(\phi+(\epsilon_{nk}-\epsilon_\nu)\dot\phi\), as required to project each
energy eigenstate onto the fixed projector. This is not an energy-dependent
projector.

## Static limit and fail-closed behavior

The driver requires at least four finite, positive, strictly descending
broadenings. It records the complex sum-rule interaction \(U(\eta)\) at each
configured eta as a diagnostic. It fits the real part of \(\chi_0(0+i\eta)\)
against \(\eta^2\) using the final four samples, compares the final-three and
final-four intercepts, and reports fit residual and the smallest-eta
imaginary-to-real ratio. The physical \(U\) is solved once from the
extrapolated static \(\chi_0\), with \(\eta=0\) in the static solve. It does
not use \(\operatorname{Re}U(\eta_{\min})\) as the static definition.

The route fails closed if the static solve is rank deficient or ill
conditioned, if the real sum-rule residual is too large, or if the intercept,
fit, imaginary response, or smallest-eta \(U(\eta)\) checks do not converge.
Configured finite-eta \(U\) values are never fitted independently and
selected by agreement. No Goldstone correction is added.

## Independent checks

The independent unit oracle reconstructs the \(E_F\) large/small d radial
function and normalization by direct loops, then independently calculates
selected band projections, the projected moment, and the Lehmann transitions.
Negative controls alter normalization, projector energy, angular weight,
and the spin-flip factor. It checks the \(4\pi\) contraction separately.
This certifies the implementation of the stated LMTO projection equations;
it does not eliminate the KKR-to-LMTO radial approximation.

The bcc-Fe check uses a single k-space SCF acceptance and then passes that
live reciprocal eigensystem directly into the Juelich-d runner. There is no
SCF after the acceptance and no response-stage Hamiltonian rebuild. The
runner writes the accepted-state provenance artifact and reports the raw
projected static \(\chi_0\), projected moment, every \(U(\eta)\), extrapolated
static \(U\), static sum-rule residual, rank/condition, and denominator action
on the magnetic mode. The input requests only Gamma and \(\omega=0\); it does
not calculate a dispersion.

Material check record:

- accepted reciprocal artifact SHA-256: 12f4378f9f3adb93d4f261d53d50a53baaa216a1d4d4999511a4a41f3b07fb39
- exact TDDFT left-state artifact SHA-256: 45b4a7fe8bbecaa153cf69318a5a2dd87176bf00a44f76f604c2b74c097c1076
- direct accepted-cache handoff: true; k mesh: \(4\times4\times4\), 64 k points; \(E_F=-0.062244856509544519\) Ry; \(T=300\) K
- accepted total moment: \(2.0000013093297557\,\mu_B\); SCF residual: \(9.23\times10^{-7}\)
- raw normalization integrals \(\mathcal N_{d\uparrow}=1.0491047071\), \(\mathcal N_{d\downarrow}=1.0252035830\); these are the pre-normalization integrals, so the normalized projector has unit norm
- projected d moment: \(1.8663102573116879\,\mu_B\)
- extrapolated \(\chi_0^{dd}(0)=-179.02959011566327\) Ry\(^{-1}\); \(U_d(0)=-0.0055856688235388537\) Ry/\(\mu_B\)
- static sum-rule rank: 1; condition: 1; relative residual: 0; relative denominator action on \(M_d\): 0
- final-three versus final-four intercept estimate: \(5.71\times10^{-9}\); fit residual: \(8.80\times10^{-9}\); smallest-eta imaginary-to-real ratio: \(1.9561\times10^{-3}\)
- smallest-eta \(U\) relative difference from static \(U\): \(3.06\times10^{-7}\); smallest-eta \(|\operatorname{Im}U|/|\operatorname{Re}U|=1.9561\times10^{-3}\)

The diagnostic complex interactions were:

| \(\eta\) (Ry) | \(\operatorname{Re}U(\eta)\) (Ry/\(\mu_B\)) | \(\operatorname{Im}U(\eta)\) (Ry/\(\mu_B\)) | Rank / condition | Real-constrained residual |
|---:|---:|---:|---:|---:|
| 0.0400000 | -0.00561398401225 | 0.00139904380086 | 1 / 1 | 0.2418113029 |
| 0.0200000 | -0.00559279062707 | 0.00069932646532 | 1 / 1 | 0.1240745005 |
| 0.0100000 | -0.00558745200890 | 0.00034963915267 | 1 / 1 | 0.0624536211 |
| 0.0050000 | -0.00558611476709 | 0.00017481657887 | 1 / 1 | 0.0312795286 |
| 0.0025000 | -0.00558578029377 | 0.00008740791516 | 1 / 1 | 0.0156463754 |
| 0.0012500 | -0.00558569666525 | 0.00004370391081 | 1 / 1 | 0.0078240148 |
| 0.0006250 | -0.00558567575748 | 0.00002185194956 | 1 / 1 | 0.0039121108 |
| 0.0003125 | -0.00558567053049 | 0.00001092597405 | 1 / 1 | 0.0019560683 |

The older worktree artifact has SHA-256
7a862caec0e96be5d3d2a6c2eb8cb4a08250c0676cc131703ca23382ff88171e and is
not an input loader. The material check created one accepted reciprocal
state in a single driver invocation, then used its live cache for the
response; its artifact hash differs from that earlier file. There was no
SCF after acceptance or Hamiltonian rebuild between acceptance and response.

## Juelich-spd and Juelich-spdf: derivation only

The d method has five \(m\) channels in the projector but the paper's Eq. 43
Dyson problem is still in site space. Projector size and interacting matrix
dimension are separate concepts.

For an enlarged orbital set \(\mathcal S\), let \(X_{i\alpha,j\beta}\) denote
the orbital-resolved static bare response and \(M_{j\beta}\) the projected
orbital moments. Applying one site scalar \(U_j\) to every orbital, the full
orbital sum rule is

\[
\sum_{j\beta}X_{i\alpha,j\beta}M_{j\beta}U_j=M_{i\alpha}.
\]

Summing over the output orbital \(\alpha\) gives

\[
M_i^{\mathcal S}=\sum_j U_j
\sum_{\alpha\beta}X_{i\alpha,j\beta}M_{j\beta}.
\]

When \(M_j^{\mathcal S}=\sum_\beta M_{j\beta}\ne0\), define the
moment-weighted site response
\[
\chi_{\mathcal S}^{0,ij}
=\frac{\sum_{\alpha\beta}X_{i\alpha,j\beta}M_{j\beta}}
{M_j^{\mathcal S}}.
\]
Then the aggregate identity is \(\Gamma_{\mathcal S}U=M_{\mathcal S}\),
with \(\Gamma_{\mathcal S}^{ij}=\chi_{\mathcal S}^{0,ij}M_j^{\mathcal S}\).
This is a controlled site-scalar sum-rule generalization, not a literal
Jülich projection: it introduces an enlarged orbital projector and a
moment-weighted aggregation.

Cross-orbital terms mean the unweighted orbital sum
\(\sum_{\alpha\beta}X_{i\alpha,j\beta}\) cannot in general be factored from
\(\sum_\beta M_{j\beta}\). The moment-weighted definition above preserves
the aggregate static identity, but it does not close the full orbital Dyson
equation in site space. Without an additional separability assumption, the
interacting response dimension is site times the number of selected
orbitals: d has 5, spd has 9, and spdf has 16 orbitals per site. No such
larger production solve is implemented here.

| Route | Classification | Final response dimension |
|---|---|---|
| Juelich-d | LITERATURE REFERENCE | Site space after the paper's d projection |
| Juelich-spd | CONTROLLED SITE-SCALAR GENERALIZATION | Aggregate sum rule can be site space; full orbital Dyson response is site × 9 without additional factorization |
| Juelich-spdf | CONTROLLED SITE-SCALAR GENERALIZATION | Aggregate sum rule can be site space; full orbital Dyson response is site × 16 without additional factorization |

## Production boundary

Only projection d enters the certified Juelich-d driver. The implementation
rejects other projection selections for this route. The energy-dependent
projected diagnostic remains a separate legacy path and is not used for
Juelich-d. No spd or spdf production route is defined by this task.
