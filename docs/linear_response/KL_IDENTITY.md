# Katsnelson–Lichtenstein identity and rotation semantics

Reference: M. I. Katsnelson and A. I. Lichtenstein, *Magnetic susceptibility,
exchange interactions and spin-wave spectra in the local spin density
approximation*, J. Phys.: Condens. Matter **16**, 7439 (2004),
[cond-mat/0406488](https://arxiv.org/abs/cond-mat/0406488).
Equation numbers below refer to that paper. This contract closes the existing
adiabatic observable mapping; it introduces no KL formulation or backend.

## A. Full KL formulation

**FULL KL FINITE-FREQUENCY SUSCEPTIBILITY — NOT IMPLEMENTED**

The auxiliary Kohn–Sham transverse susceptibility is \(\chi_0^{+-}\)
(Eq. 12). In ADA/TDDFT the physical spin response solves

\[
\chi^{+-}=\chi_0^{+-}+\chi_0^{+-}I_{\rm xc}\chi^{+-},
\qquad I_{\rm xc}(r)=B_{\rm xc}(r)/m(r)
\]

(Eqs. 9–10). Products are spatial operator products in their written order;
\(m\) and \(B_{\rm xc}\) act as local multiplication operators.
The divergence/gradient transition operator \(\Lambda(\omega)\) of Eq. 25
is the remainder in the identity used in Eqs. 21–23:

\[
\chi_0^{+-}B_{\rm xc}=m-\omega\chi_0^{+-}+\Lambda.
\]

Consequently the full response can be written as

\[
\boxed{\chi^{+-}=m[\omega-(\chi_0^{+-})^{-1}\Lambda]^{-1}}
\quad\text{(Eq. 24)},
\]

\[
\boxed{\chi^{+-}=(m+\Lambda)[\omega-I_{\rm xc}\Lambda]^{-1}}
\quad\text{(Eq. 26)},
\qquad \Omega(\omega)=I_{\rm xc}\Lambda(\omega)\quad\text{(Eq. 27)}.
\]

The numerator \(m+\Lambda\) matters for physical response amplitudes.
The adiabatic spin-wave limit replaces \(\Lambda(\omega)\) by
\(\Lambda(0)\). This additional spectral approximation is distinct from
using an adiabatic XC functional in a finite-frequency TDDFT calculation.
The current rotation route evaluates a different object, an electronic
local-rotation effective action. This classification does not imply a defect
in that formulation. Full KL finite-frequency susceptibility is deferred.

## B. Adiabatic/MFT limit

**KL ADIABATIC MFT SPIN-WAVE LIMIT — CERTIFIED**

Use the old-MFT exchange observable, with the production symmetrization

\[
J(q)\equiv J_{\rm sym}(q)=\tfrac12[J_{\uparrow\downarrow}(q)
 +J_{\downarrow\uparrow}(q)],\qquad \Delta J(q)=J(0)-J(q).
\]

The native static service constructs \({\cal C}(q)=2\Delta J(q)\).
Independently the finite-H route computes the grand-potential Hessian

\[
\boxed{H_{\theta\theta}(q)=\text{torque--torque}+\text{contact}
 =2[J(0)-J(q)]}.
\]

This is the validated static identity in the accepted Hamiltonian/fixture
sector. Numerical cross-checks retain their existing finite-mesh and contour
tolerances; certification does not assert universal material convergence.
The direct spin-commutator Berry oracle gives

\[
\boxed{B_{\rm Berry}=M_{\rm band}/2},\qquad
\partial_\omega K_+=-B_{\rm Berry},\qquad
\partial_\omega K_-=+B_{\rm Berry}.
\]

Thus, for the positive-moment convention \(M=M_{\rm band}>0\),

\[
\boxed{\omega(q)=\frac{H_{\theta\theta}(q)}{B_{\rm Berry}}
 =\frac{2H_{\theta\theta}(q)}{M}
 =\frac{4[J(0)-J(q)]}{M}}\quad\text{(KL Eq. 42)}.
\]

This certifies the adiabatic MFT equation, not a finite-frequency amplitude
identity. Signed circular roots remain \(K_+(q,0)/B\) and
\(-K_-(q,0)/B\); the positive-frequency channel depends on the signs.
The production `adiabatic_finiteH_meV` and `adiabatic_Turek_meV_or_missing`
columns use \(|M_{\rm band}|\), with the existing small-moment denominator
floor, as reporting scales. They do not replace signed circular pole analysis.
The automatic scalar pole workflow is one-site only.

## C. Native RS-LMTO realization

**EXACT OBSERVABLE-LEVEL MAPPING TO KL ADIABATIC MFT**

The Turek/LKAG service uses LMTO potential-function splitting, path
operators/Green functions, and ordered up/down and down/up contractions.
It computes the same old-MFT exchange observable entering KL Eq. 42.
It is not a literal matrix-by-matrix implementation of KL Eqs. 38–40
(the Hamiltonian XC splitting, KS bubble, and exchange sandwich).

| Production symbol / procedure | Identity |
|---|---|
| `native_turek_static_reference` in `source/lmto_path_operator_contour.f90` | `jq_sym=0.5*(jq_ud+jq_du)`; `delta_j=J(Gamma)-J(q)`; `curvature=2*delta_j` |
| `compute_static_rotation_curvature` in `source/linear_response_rotation.f90` | Independent finite-H/H2 torque–torque plus contact Hessian |
| `global_spin_berry` / `prepare_rotation_response` | Direct spin commutator and independent residual against `0.5*magnetization` |
| `run_rotation_response_workflow` | `finite_mev=2*finite_h(1,1,iq)/abs(M_band)*ry_to_mev`; native diagnostic `4*native_delta(iq)/abs(M_band)*ry_to_mev` (existing denominator floor retained) |

Native Turek remains an optional independent static cross-check. It does not
set the branch, pole window, acceptance, Berry normalization, or Goldstone
treatment. No second exchange service or inverse-response solver is needed.

## D. Dynamic rotation formulation

**LOCAL-ROTATION EFFECTIVE-ACTION KERNEL**

For the accepted second-order LMTO Hamiltonian, define the local transverse
rotation vertex \(T_A(q)=\partial H/\partial\theta_A(q)\). Production evaluates

\[
\boxed{K^R_{AB}(q,\omega)=C_{AB}(q)+\Pi^R_{AB}(q,\omega)},
\]

\[
\Pi^R_{AB}=\sum_{knm}w_k
 \frac{f_{nk}-f_{m,k+q}}{\omega+\epsilon_{nk}-\epsilon_{m,k+q}+i\eta}
 T_A^{mn}T_B^{nm},\qquad
C_{AB}=\left\langle\frac{\partial^2 H}
 {\partial\theta_A\partial\theta_B}\right\rangle.
\]

The exact static symmetric reduction is the force-theorem Hessian.
The dynamic rotation kernel is not identified with the full KL transverse
susceptibility. Its inverse is a rotation-coordinate propagator; physical
spin-susceptibility amplitudes require a separate source-coupling derivation.
Neither `K = chi^{-1}` nor `K^{-1} = chi` is an established production identity.

Zeros of \(\det K\), or of scalar circular \(K_\pm\), diagnose a collective
local-rotation pole. The real-kernel crossing, minimum \(|K|\), peak in
\(-\operatorname{Im}K^{-1}\), eta/FWHM behavior, and causality sign are
supporting rotation-mode diagnostics. Pole equivalence to a physical spin
response is expected in the rigid-rotation sector, but is not a proven full
finite-frequency response identity. Spectral amplitudes are not certified
TDDFT spin-susceptibility weights.

Public columns are `minus_Im_Kinv_plus_invRy`,
`minus_Im_Kinv_minus_invRy`, `rotation_spectral_peak_meV`, and
`max_minus_Im_Kinv_invRy`. The spectral peak and FWHM refer to the signed
inverse rotation kernel. The shared input names `chi_plus` / `chi_minus`
select circular channels for rotation; they do not assert a susceptibility
amplitude identity.

## Independent evidence retained

`UnitLrRotationProductionAdapter` retains accepted-Hamiltonian fixture
reconstruction, direct Berry oracle, signed slopes, negative controls,
Goldstone identity, and q/-q/-omega covariance. `UnitLrRotationStaticLkag`
retains the independent native static comparison against the historical
truncated two-shell input reference; that diagnostic alone is not an exact
full-spd equality gate.
`UnitLrRotationSecondOrder`, `UnitLrRotationTorqueHessian`, and the finite-q
and static-contour tests retain finite-H versus direct spectral-oracle and
contact-term validation. `UnitLrRotationPoleCapability` retains the pole
workflow boundary. The archived numeric rotation baseline is preserved;
only schema labels change. No tolerances or numerical formulas change.
