# Independent TD-DFT Ward derivation for RS-LMTO-ASA

Status: TDWARD-01 derivation complete. No production sign or factor was
changed as part of this task.

## Determination

The implementation fixes the following convention:

* \(H_{\rm spin}=H_0\sigma_0+\mathbf H\mathbin{\cdot}\boldsymbol{\sigma}\).
* \(B^E_{\rm xc}=(V_{\uparrow}-V_{\downarrow})/2\), with no leading minus.
* \(O_\pm=\sigma_x\pm i\sigma_y\) are measured without a half; a circular
  Hamiltonian source carries \(1/2\).
* In the co-rotating channel,
  \[
  \chi^{\rm co}_{\rm KS}(0,0)B^{\rm circ,src}_{\rm xc}=m_G.
  \]
* The direct pair object obeys
  \[
  \Xi^{\rm co}(0,0)m_G=m_G,\qquad \Xi=\chi_{\rm KS}Q.
  \]

An independent ground-state XC field exists: converged
`XCPOT_hybrid` output is accumulated by `VXC0SP` in
`xc_response_radial_projection%bxc_spin_moment`. TDWARD-02 exposes that
radial provenance and constructs the signed source from it; the fallback
`k_perp_circular*m_G` path is retained only for explicitly non-production
debug/unit callers and is labelled as derived. The radial site projection
and orbital LMTO pair tangent must not be assumed equivalent.

## 1. Ground-state Hamiltonian convention

The basis is site-major and spin-major:

```text
(site 1 orbitals up, ..., site N orbitals up,
 site 1 orbitals down, ..., site N orbitals down)
```

The mapping assembled in
[hamiltonian_build.f90](../source/hamiltonian_build.f90) is

\[
H_{\rm spin}=H_0\sigma_0+H_x\sigma_x+H_y\sigma_y+H_z\sigma_z,
\]
\[
H_{\uparrow\uparrow}=H_0+H_z,\quad
H_{\downarrow\downarrow}=H_0-H_z,\quad
H_{\uparrow\downarrow}=H_x-iH_y,\quad
H_{\downarrow\uparrow}=H_x+iH_y.
\]

The same mapping is used for `hxc`, but `hxc` is an assembled
magnetic Hamiltonian block, not an independently evaluated XC kernel.
FSM/constraining fields are inserted separately by
`add_constraining_field_hmag`.

In `symbolic_atom%predls`, ordinary LMTO parameters are split as

\[
c_0=\frac{c_{\uparrow}+c_{\downarrow}}2,\qquad
c_1=\frac{c_{\uparrow}-c_{\downarrow}}2.
\]

These are the code fields `cx0/cx1`; the same split is carried by
`wx0/wx1` and `cex0/cex1`. `ham0m_nc` consumes the
ordinary pair for the standard Hamiltonian and the shifted pair for its
explicit HOH path.

The magnetic tangent service in
[lmto_magnetic_tangent.f90](../source/lmto_magnetic_tangent.f90) uses
\(\mathbf H_i=\mathbf H_{i,\rm bond}+c_{1,i}\mathbf e_i\), with
`potential%mom` as the unit orientation. Thus `c1` is an
up-minus-down half-splitting coefficient, not the full difference and not,
by itself, a local XC field.

## 2. Independent \(B^E_{\rm xc}\) definition and provenance

The direct functional call is

```text
XCPOT_hybrid(rho_down, rho_up, rho_total, ..., v_down, v_up, exc)
```

as shown in [xc.f90](../source/xc.f90). Therefore

\[
V_{\rm xc}=v^0_{\rm xc}\sigma_0+B^E_{\rm xc}\sigma_z,\qquad
B^E_{\rm xc}=\frac{V_{\uparrow}-V_{\downarrow}}2.
\]

`evaluate_ground_state_xc_sample` in
[xc_response_kernel.f90](../source/xc_response_kernel.f90) stores this
half-difference in `sample%bxc_energy`. The production
`VXC0SP` route calls the same functional, adds FSM afterward, and
passes direct `VXC2,VXC1` output to
`xc_response_radial_projection%accumulate`. That routine forms

\[
m_z(r)=\rho_{\uparrow}(r)-\rho_{\downarrow}(r),\qquad
N_{B,i}=\int dr\,w_i(r)m_z(r)B^E_{{\rm xc},z}(r),
\]

stored as `bxc_spin_moment`. With positive response amplitude
\(M_i\), finalization gives

\[
B^E_{\parallel,i}=\frac{N_{B,i}}{M_i},\qquad
K^{\rm circ}_i=\frac{B^E_{\parallel,i}}{2M_i}
               =\frac{N_{B,i}}{2M_i^2}.
\]

For signed collinear \(m_{G,i}=s_iM_i\), \(s_i=\pm1\), the lab coefficient is
\(B^E_{z,i}=s_iB^E_{\parallel,i}\). The radial product is reversal-even, so
the finalized provider field is a local-aligned scalar; the Ward source must
restore \(s_i\). The legacy `Bxc_tot` is `Bxc_up-Bxc_dw`, the full
\(V_{\uparrow}-V_{\downarrow}=2B^E_{\rm xc}\), and is not this source.

## 3. Circular operators and source normalization

\[
O_+=\sigma_x+i\sigma_y=
\begin{pmatrix}0&2\\0&0\end{pmatrix}=2\sigma_+,\qquad
O_-=\sigma_x-i\sigma_y=
\begin{pmatrix}0&0\\2&0\end{pmatrix}=2\sigma_-.
\]

For \(B_\pm=B_x\pm iB_y\),

\[
\delta H=\delta B_x\sigma_x+\delta B_y\sigma_y
 =\frac12(\delta B_+O_-+\delta B_-O_+),
\]
\[
\frac{\partial H}{\partial B_+}=\frac{O_-}{2},\qquad
\frac{\partial H}{\partial B_-}=\frac{O_+}{2}.
\]

`response_operator(RESPONSE_PLUS/MINUS)` returns unhalved \(O_\pm\),
while `external_source_operator` returns the derivative with \(1/2\).
The constants `tddft_circular_operator_factor=2` and
`tddft_circular_source_factor=0.5` make this explicit. Hence the
circular kernel is \(B^E_\parallel/(2M)\), and the Cartesian kernel is twice it.

## 4. Static \(\chi_{\rm KS}\) and the Ward sign

The eigenpair engine contracts

\[
\chi^{AB}_{ij}(\mathbf q,\omega)=
\sum_{\mathbf k,n,m}w_{\mathbf k}
\frac{f_{n\mathbf k}-f_{m,\mathbf k+\mathbf q}}
{\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta}
\langle n\mathbf k|A_i|m,\mathbf k+\mathbf q\rangle
\langle m,\mathbf k+\mathbf q|B_j|n\mathbf k\rangle.
\]

The static route in [tddft_chi0.f90](../source/tddft_chi0.f90) uses
\((f_n-f_m)/(\epsilon_n-\epsilon_m)\), with the Fermi derivative at
degeneracy and no dynamic \(\eta\). For the one-site oracle
\(H=e_0+b\sigma_z\), \(b<0\), with \(M=f_\uparrow-f_\downarrow>0\),

\[
\chi^{+-}_{\rm KS}(0,0)
=\frac{4M}{\epsilon_\uparrow-\epsilon_\downarrow}
=\frac{2M}{b}.
\]

The independent circular source is \(b/2\), so

\[
\chi^{+-}_{\rm KS}(0,0)\frac b2=M.
\]

Equivalently \(K^{\rm circ}=b/(2M)\) gives
\(\Xi=\chi^{+-}_{\rm KS}K^{\rm circ}=1\) and \(\Xi M=M\). This fixes the
Ward sign to plus. In Cartesian variables \(\chi^{xx}=M/b\),
\(K^x=b/M\), and there is no circular half.

## 5. Signed Goldstone vector and reversal

The production path obtains

\[
m_{G,i}=s_iM_i,\qquad M_i>0,\qquad s_i=\operatorname{sign}(m_{z,i}),
\]

from `compute_kspace_spin_moments_spinor`. `calculation.f90`
stores signed `site_moments(3,:)` as `signed_mz` and separately
stores `moment_amplitudes=sqrt(sum(site_moments**2))` for pair
normalization. `signed_site_populations` returns the signed vector.

The independent lab-frame circular source is

\[
B^{\rm circ,src}_{{\rm xc},i}=\frac{s_iB^E_{\parallel,i}}2
=K^{\rm circ}_im_{G,i},
\]

and the exact site-resolved statement is

\[
\boxed{\chi^{\rm co}_{\rm KS}(0,0)B^{\rm circ,src}_{\rm xc}=m_G.}
\]

The equality to `K*m` is only the local projected-kernel
consistency relation; TDWARD-02 must form the left-hand source from direct
XC data first.

For global reversal, \(U=\sigma_x\) gives

\[
U\sigma_xU=\sigma_x,\quad U\sigma_yU=-\sigma_y,\quad
U\sigma_zU=-\sigma_z,\quad UO_+U=O_-,\quad UO_-U=O_+.
\]

All \(M_i\) remain positive, every \(s_i\) and lab \(B^E_z\) changes sign,
and the co-rotating channel switches from `+-` to `-+`.
Both sides of the Ward identity change sign, so the eigenvalue remains \(+1\),
not \(-1\). For a two-sublattice \(+z/-z\) state,
\(m_G=(+M_A,-M_B)\) and
\(B^{\rm circ,src}=(+B^E_{\parallel,A}/2,-B^E_{\parallel,B}/2)\); the full
site matrix must act on this signed vector.

## 6. Pair-potential \(\Xi\)

The LMTO pair route differentiates the assembled `ham_only` Hamiltonian
with respect to orientation. Its independent endpoint tangents are
\(D_{x,i}=\partial H/\partial e_{x,i}\) and
\(D_{y,i}=\partial H/\partial e_{y,i}\). The
[lmto_pair_potential.f90](../source/lmto_pair_potential.f90) service forms

\[
Q^-_i=\frac{D_{x,i}-iD_{y,i}}{2M_i},\qquad
Q^+_i=\frac{D_{x,i}+iD_{y,i}}{2M_i}.
\]

The reverse endpoint is assembled separately for \(Q^+\). The reciprocal
builder retains the same site-major `ham_only` representation as the
endpoint eigenvectors and uses \(k\) and \(k+q\) on the two endpoints.
The static builder in [tddft_xi.f90](../source/tddft_xi.f90) contracts

\[
\Xi^{A,Q}_{ij}(0,0)=
\sum_{\mathbf k,n,m}w_{\mathbf k}
\frac{f_n-f_m}{\epsilon_n-\epsilon_m}
\langle n\mathbf k|A_i|m\mathbf k\rangle
\langle m\mathbf k|Q_j|n\mathbf k\rangle.
\]

Rigid spin-rotation covariance gives
\(\Xi^{\rm co}(0,0)m_G=m_G\). Under reversal \(M_i\) stays positive,
`Q^-` and `Q^+` exchange circular sectors, and \(m_G\)
changes sign; the corresponding Xi eigenvalue remains \(+1\).

The legacy site-kernel construction is the separate
`circular_transverse_kernel`/`cartesian_transverse_kernel` path in
`xc_response_kernel.f90`; it returns \(K^{\rm circ}\) and exactly twice
that value in Cartesian components. `build_site_projected_k_perp` copies
the circular scalar into the Goldstone layer. `evaluate_goldstone` then
forms `xi_raw=construct_transverse_xi(chi_static,k_perp)` and, for its
production Ward report, forms
`bxc=signed_magnetization/moment_amplitude*bxc_spin_moment/(2*moment_amplitude)`
from the independent radial numerator.
`evaluate_raw_xi_diagnostics` instead evaluates the raw pair action
`Xi*m_G-m_G`. Finally, `evaluate_static_ward_identity` evaluates
`chi*bxc-m_G`, with `r_B` enabled only for the independent source. Xi-only
records carry `derived_identity_residual`; a debug `K*m_G` source is never
reported as independent by provenance label alone.

## 7. Provenance map and narrow TDWARD-02 plumbing

| quantity | exact source | use |
| --- | --- | --- |
| signed \(m_G\) | `compute_kspace_spin_moments_spinor`, `signed_mz`, provider `signed_spin_population` | signed Ward vector |
| positive \(M_i\) | `moment_amplitudes` | pair \(Q\) denominator |
| independent pointwise field | `evaluate_ground_state_xc_sample` and direct `XCPOT_hybrid` output | sign/factor oracle |
| independent production field | `VXC0SP` → `xc_response_radial_projection%accumulate` → `bxc_spin_moment` | form \(s_iN_{B,i}/(2M_i)\) |
| static \(\chi_{\rm KS}\) | `build_static_chi_ks_from_eigenpairs`, `tddft_static_divided_difference` | real static route |
| pair \(\Xi\) | `build_static_direct_xi_from_operator_source`, `build_lmto_pair_potential_at_kpoint` | independent tangent oracle |
| Ward report | `evaluate_static_ward_identity` | explicit source with provenance |

Expose `bxc_spin_moment` together with positive \(M_i\) and signed
orientation, then construct \(B^{\rm circ,src}_i=s_iN_{B,i}/(2M_i)\). Keep
direct `B_xc` and derived `K_xc*m` provenance separate.
Do not use `hamiltonian%hxc`, `potential%cx1`, or legacy
`Bxc_tot` as an independent XC source.

## 8. Falsification tests

1. **One-site \(+z\):** genuine LMTO state, direct field, static `+-`
   channel; require \(\chi B^{\rm circ,src}-m_G\), \(\Xi m_G-m_G\), and
   Goldstone eigenvalue \(+1\). Full difference, missing half, and minus
   sign are negative controls.
2. **One-site \(-z\):** reverse actual LMTO orientation and XC channels;
   require positive \(M\), signed source/vector reversal, channel switch to
   `-+`, and eigenvalue \(+1\).
3. **Two-site \(+z/-z\):** use \(m_G=(+M_A,-M_B)\), test the full site matrix
   and acoustic vector; independent absolute values must fail.
4. **Static direct versus pair Xi:** compare only at \(q=0,\omega=0\) with
   identical eigenpairs, occupations, gauge, `ham_only` basis,
   channel, and tangent source. A radial scalar versus a full orbital tangent
   is not automatically equivalent.

Focused checks on the current tree passed:

```text
UnitTddftWardConventions     RESULT: PASS
UnitTddftWard                RESULT: PASS
UnitTddftGoldstone           RESULT: PASS
UnitTddftDirectXi             RESULT: PASS
UnitLmtoPairPotential         max error 8.8558E-10; +z/-z eigenvalues 1/1
UnitLmtoMagneticTangents      max error 4.8566E-10
```

The pair fixture also reports two-sublattice signed action
`(+z,-z) = (1,-1)` and a pre-repair reverse-sign negative control.
`UnitTddftGoldstone` additionally verifies independent radial `r_B`,
phase-sensitive anti-Goldstone rejection, output labels, and branch policy;
`UnitTddftChiKS` exercises the static/dynamic eta ladder at finite `+q`.
These checks do not prove equivalence of the radial provider and full-LMTO
source; they verify that any mismatch remains observable rather than being
hidden by a tautological diagnostic.

## Risks and non-claims

* `site%bxc_energy` after radial finalization is a projected,
  local-aligned energy coefficient; `site%k_perp_circular*m_G` is derived.
* `hxc`, `cx1`, and pair tangents are effective LMTO
  Hamiltonian quantities and may contain transformed/hopping contributions.
* Finite-\(\eta\) dynamic response is continuity evidence only; it never
  supplies the static Ward field or a repaired kernel.
