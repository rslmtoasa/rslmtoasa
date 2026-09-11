# LR-03 — TDDFT transverse-response conventions and normalization

## Status and authority

**PASS — conventions locked for the certified no-SOC, collinear, orthogonal
RS-LMTO baseline, with one explicit and quantified SR→Pauli approximation.**

This file is the single convention contract for later transverse-response
work. A response consumer must use the symbols, operator order, factors,
measures, phases, and units defined here. It may not independently choose a
sign, a circular-channel label, a factor of two, a `4π`, or a unit conversion.

This is an authoritative derivation and convention audit. It adds no
production `chi^KS`, no XC kernel, no Dyson solve, no Goldstone correction, no
mode extraction, and no material-spectrum fit. The independent algebra check
described in §15 is not physical validation.

The certified scope is:

* scalar-relativistic, collinear, no-SOC, two-spin-channel RS-LMTO;
* reciprocal `ham_only`, orthogonal second-order/HOH or first-order LMTO;
* site-major coefficient storage with the up block followed by the down block;
* complex `sp`/`spd` angular basis and a direct logarithmic radial mesh;
* the Pauli/no-SOC response projection, explicitly distinguished from the
  exact scalar-relativistic density operator.

SOC, noncollinearity, generalized-overlap response, additive operators, and
`spdf` complex-harmonic response are not silently included. [CODE CONVENTION]
[DEFERRED]

## Evidence reused

This audit consumes the earlier contracts rather than redoing them.

| contract | authoritative evidence consumed | use here |
|---|---|---|
| LR-GF-01 | [`LR-GF-01_GF_CONTRACT_EVIDENCE.md`](LR-GF-01_GF_CONTRACT_EVIDENCE.md) | orthogonal Lehmann ordering, `+iη`, Ry energy, reciprocal `k` convention |
| LR-BASIS-00 | [`LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md`](LR-BASIS-00_RADIAL_AUGMENTATION_EVIDENCE.md) | LMTO state reconstruction, spin-block order, radial functions and units |
| LR-01 | [`LR_RADIAL_GROUND_STATE_AUDIT.md`](LR_RADIAL_GROUND_STATE_AUDIT.md) | `RHO`, physical densities, Vxc arrays, XC provenance, magnetic-moment reporting |
| LR-02R | [`LR_RESPONSE_BASIS_MAPPING.md`](LR_RESPONSE_BASIS_MAPPING.md) | complex harmonics, Gaunt normalization, response cutoff, super-index, radial measure, Fourier/site gauge |
| LR-02N | [`LR_SR_PAULI_NUMERICAL_CLOSURE.md`](LR_SR_PAULI_NUMERICAL_CLOSURE.md) | measured distinction between `m^SR` and `m^P` |

The literature contracts are the retained IDs `BES-01…06` and `LCMM-01…05` in
[`00A_LITERATURE_CONTRACTS.md`](dev/RS_LMTO_TDDFT_cleanroom_Luna/00A_LITERATURE_CONTRACTS.md).
The primary sources are [Buczek–Ernst–Sandratskii, Phys. Rev. B **84**,
174418 (2011)](https://www-old.mpi-halle.mpg.de/mpi/publi/pdf/10475_11.pdf)
and [Lounis–Costa–Muniz–Mills, Phys. Rev. B **83**, 035109
(2011)](https://juser.fz-juelich.de/record/14113/files/PhysRevB.83.035109.pdf).
[LITERATURE]

## 1. Canonical notation and claims discipline

The following labels are used throughout:

* **[LITERATURE]** — a definition or equation imported from a cited paper.
* **[CODE CONVENTION]** — a fact fixed by the certified RS-LMTO source/audit.
* **[DERIVED MAPPING]** — algebra translating a literature object to the
  canonical RS object.
* **[ALGEBRAIC TEST]** — a material-independent exact or numerical identity.
* **[NUMERICAL INPUT FROM LR-02N]** — a measured result from the SR→Pauli
  closure, not a TDDFT validation.
* **[DEFERRED]** — intentionally not implemented or not certified here.

No algebraic test below is a spectrum, a Goldstone test, or evidence that
TDDFT is physically correct. The Fe numbers from LR-02N quantify only the
chosen ground-state projection error. [CODE CONVENTION] [NUMERICAL INPUT FROM
LR-02N]

## 2. Fundamental density operators

The published Pauli density operators are taken literally:

\[
 \hat n^\mu(\mathbf r,t)=
 \hat\psi^\dagger(\mathbf r,t)\,\Gamma^\mu\,\hat\psi(\mathbf r,t),
 \qquad \mu\in\{0,x,y,z\}.
\]

The actual local RS-LMTO coefficient spinor is, for every orbital/site block,

\[
 \hat\psi=
 \begin{pmatrix}\hat\psi_\uparrow\\[2pt]\hat\psi_\downarrow\end{pmatrix},
 \qquad
 \text{coefficient order}=(\text{all orbitals, up};\ \text{all orbitals, down}).
\]

Equivalently, for a site with `norb` orbitals, entries `1…norb` are up and
entries `norb+1…2*norb` are down. This is the ordering used by
`reciprocal_spin_density` and `spin_off=norb`. [CODE CONVENTION]

In that order, the four matrices are exactly

\[
 \Gamma^0=I=
 \begin{pmatrix}1&0\\0&1\end{pmatrix},
 \quad
 \Gamma^x=\sigma_x=
 \begin{pmatrix}0&1\\1&0\end{pmatrix},
\]
\[
 \Gamma^y=\sigma_y=
 \begin{pmatrix}0&-i\\ i&0\end{pmatrix},
 \quad
 \Gamma^z=\sigma_z=
 \begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\]

The canonical densities are

\[
 n\equiv n^0,
 \qquad s_x\equiv n^x,
 \qquad s_y\equiv n^y,
 \qquad s_z\equiv n^z.
\]

For the collinear z baseline,

\[
 n(r)=n_\uparrow(r)+n_\downarrow(r),
 \qquad
 s_z(r)=n_\uparrow(r)-n_\downarrow(r).
\]

Here `s` is a **number-spin density**. It is not a magnetic moment until the
conversion in §3 is applied. [LITERATURE] [CODE CONVENTION]

## 3. Number-spin density versus magnetic-moment density

The electron spin operator is

\[
 \hat{\mathbf S}=\frac{\hbar}{2}\boldsymbol\sigma,
\]

and the electron magnetic moment is opposite to its spin:

\[
 \hat{\boldsymbol\mu}_e
 =-\frac{g\mu_B}{\hbar}\hat{\mathbf S}
 =-\frac{g\mu_B}{2}\boldsymbol\sigma.
\]

Therefore the physical electron magnetic-moment density is

\[
 \boxed{\mathbf M_e(\mathbf r)=
 -\frac{g\mu_B}{2}\,\mathbf s(\mathbf r)}.
 \]

In particular,

\[
 M_{e,z}(r)=-\frac{g\mu_B}{2}s_z(r).
\]

The BES paper uses the standard Pauli operators for its fundamental density
operators and couples a physical field through `-g μB B/2`; its published
quantity called `m` is consequently the Pauli number-spin density `s`, not
`M_e`. BES uses `g=2`. [LITERATURE] [DERIVED MAPPING]

RS-LMTO's accepted ground-state report uses a positive up-minus-down reporting
convention:

\[
 M_z^{\mathrm{code}}(r)=+\mu_B s_z(r),
 \qquad
 M_z^{\mathrm{code}}=\mu_B(N_\uparrow-N_\downarrow).
\]

This is a reported code moment, not the signed physical electron moment. For
`g=2`, `M_z^{code}=-M_{e,z}`. No extra electron minus sign is inserted into
the existing report. [CODE CONVENTION]

### Canonical density table

| symbol | meaning | canonical unit | sign/conversion |
|---|---|---|---|
| `n=n^0` | Pauli number density | electrons bohr`^-3` | `n_up+n_down` |
| `s_i=n^i` | Pauli number-spin density | electrons bohr`^-3` | `s_z=n_up-n_down` |
| `M_e` | physical electron magnetic-moment density | `μB` bohr`^-3` | `M_e=-(g μB/2)s` |
| `M^code` | RS reported positive up-minus-down moment | `μB` bohr`^-3` | `M^code=+μB s` |
| `m` in BES Eq. (22) | BES Pauli number-spin density | electrons bohr`^-3` | `m=s`, not a physical signed moment |
| `m_z` in LCMM sum rule | LCMM up-minus-down number density | electrons bohr`^-3` | `m_z=s_z` in the mapped contract |

The symbol `m` is not used in this document without a superscript or a
literature qualifier. [DERIVED MAPPING]

## 4. XC potential and field dictionary

Start only from the accepted LR-01 arrays `V_xc,up(r)` and
`V_xc,down(r)`, both in Ry. Define

\[
 V_{xc}^0(r)=\frac{V_{xc,\uparrow}(r)+V_{xc,\downarrow}(r)}{2},
 \qquad
 B_{xc}^{\sigma}(r)=
 \frac{V_{xc,\uparrow}(r)-V_{xc,\downarrow}(r)}{2}.
\]

Also retain, separately,

\[
 \Delta V_{xc}(r)=V_{xc,\uparrow}(r)-V_{xc,\downarrow}(r)
 =2B_{xc}^{\sigma}(r).
\]

The canonical spin-space matrix is

\[
 V_{xc}(r)=V_{xc}^0(r)I+B_{xc}^{\sigma}(r)\sigma_z
 =\begin{pmatrix}V_{xc,\uparrow}&0\\0&V_{xc,\downarrow}\end{pmatrix}.
\]

`B_xc^σ` is an **energy-valued Pauli coefficient**. It is not, by its name,
the physical magnetic field in tesla. The certified up/down ordering fixes the
sign; no opposite sign is permitted. [CODE CONVENTION] [DERIVED MAPPING]

If a physical field is introduced through the electron Zeeman coupling,

\[
 H_B= -\frac{g\mu_B}{2}\,\mathbf B\cdot\boldsymbol\sigma,
\]

then the field associated with the XC Pauli coefficient is

\[
 \boxed{B_{xc}^{\mathrm{phys}}
 =-\frac{2B_{xc}^{\sigma}}{g\mu_B}
 =-\frac{\Delta V_{xc}}{g\mu_B}.}
\]

For `g=2`, `B_xc^phys=-B_xc^σ/μB`. The minus sign is the electron sign; it
does not disappear when Ry atomic units are used. [DERIVED MAPPING]

### Names used by the literature and code

| published or code name | canonical meaning |
|---|---|
| `B_xc^σ` / `bxc_pauli` | `ΔV_xc/2`, energy-valued Pauli coefficient, Ry |
| `ΔV_xc` | `V_xc,up−V_xc,down`, energy splitting, Ry |
| BES `B_xc` | physical XC magnetic field in the coupling `-gμB B_xc/2`; map with `B_xc^σ=−gμB B_xc/2` |
| LCMM `B_eff` | energy-valued `V_down−V_up`; map with `B_eff=−ΔV_xc=−2B_xc^σ` for XC only |
| exchange splitting, if defined as `V_up−V_down` | `+ΔV_xc=+2B_xc^σ` |
| exchange splitting, if defined as `V_down−V_up` | `−ΔV_xc=−2B_xc^σ` |
| `B_fsm` | separate scalar constraining field in the RS Hamiltonian, not XC |

With LR-01's insertion `V_up += B_fsm`, `V_down -= B_fsm`, the constraining
contribution to the down-minus-up splitting is `-2B_fsm`. Thus, ignoring other
spin-dependent terms,

\[
 B_{\mathrm{eff}}^{\mathrm{total}}
 =V_{KS,\downarrow}-V_{KS,\uparrow}
 =-\Delta V_{xc}-2B_{fsm}.
\]

`B_fsm` must never be silently included in `B_xc^σ`. [CODE CONVENTION]

## 5. Circular operators and factors of two

Define the halved ladder matrices

\[
 \sigma^\pm=\frac{\sigma_x\pm i\sigma_y}{2},
\]

so that, in the certified `(up,down)` order,

\[
 \sigma^+=\begin{pmatrix}0&1\\0&0\end{pmatrix},
 \qquad
 \sigma^-=\begin{pmatrix}0&0\\1&0\end{pmatrix}.
\]

They act as `σ^+|down⟩=|up⟩` and `σ^-|up⟩=|down⟩`. The unhalved matrices are

\[
 \Sigma^\pm\equiv\sigma_x\pm i\sigma_y=2\sigma^\pm
 =\begin{cases}
 \begin{pmatrix}0&2\\0&0\end{pmatrix},&+\\[8pt]
 \begin{pmatrix}0&0\\2&0\end{pmatrix},&-.
 \end{cases}
\]

The density combinations are deliberately unhalved:

\[
 s^\pm=s_x\pm is_y
 =\hat\psi^\dagger\Sigma^\pm\hat\psi
 =2\hat\psi^\dagger\sigma^\pm\hat\psi.
\]

For the physical electron moment and the positive code moment,

\[
 M_e^\pm=-\frac{g\mu_B}{2}s^\pm,
 \qquad
 (M^{\mathrm{code}})^\pm=+\mu_Bs^\pm.
\]

There is no `plus_minus`, `minus_plus`, or unqualified `operator_factor` in
the canonical contract. A transition vertex must state whether it uses
`σ^±` or `Σ^±`. [DERIVED MAPPING]

## 6. Circular susceptibility ordering

The canonical retarded response is defined with the **measurement operator
first** and the **source operator second**:

\[
 \chi^{ij}_{AB}(t-t')=-i\theta(t-t')
 \langle[\hat A_i(t),\hat B_j(t')]\rangle.
\]

For the Cartesian Pauli density response,

\[
 \chi^{ij}=-i\theta\langle[\hat n^i,\hat n^j]\rangle.
\]

BES defines

\[
 m^\pm=m_x\pm im_y,
 \qquad B^\pm=B_x\pm iB_y,
 \qquad \chi^\pm=\chi^{xx}\mp i\chi^{xy},
\]

with `m^±=χ^± B^±`. Here the source `B^±` means the corresponding circular
component of the driving coefficient in the density-coupled Hamiltonian; if it
means a physical magnetic field, the source coefficient is first multiplied by
`-gμB/2`. [LITERATURE]

To derive the channel order, use

\[
 \chi^+=\chi_{\sigma_x,\,\sigma_x-i\sigma_y}
 =\frac12\chi_{\Sigma^+,\Sigma^-},
\]

after the collinear axial symmetry removes the non-spin-flip contribution.
Likewise `χ^-=(1/2)χ_{Σ^-,Σ^+}`. Therefore:

| channel | measurement operator | source operator | positive-frequency spin-flip term |
|---|---|---|---|
| `χ^+` | `s^+`, equivalently `σ^+` in the vertex | `s^-` / `σ^-` | occupied up `→` empty down |
| `χ^-` | `s^-`, equivalently `σ^-` in the vertex | `s^+` / `σ^+` | occupied down `→` empty up |

The last column follows from the retarded denominator, not from the sign of a
material magnon. With halved vertices `T^±` built from `σ^±`, the published
circular factor is `2`. With unhalved vertices `\widetilde T^±` built from
`Σ^±=2σ^±`, the same expression has factor `1/2`:

\[
 2T^\pm(T^\pm)^*
 =\frac12\widetilde T^\pm(\widetilde T^\pm)^*.
\]

This is the complete factor-of-two mapping. [LITERATURE] [DERIVED MAPPING]

## 7. Retarded KS susceptibility

For a density-coupled source field `b^ν_J` in

\[
 \hat H' =\sum_{J,\nu}\int d^3r\,
 \hat n^\nu_J(\mathbf r)b^\nu_J(\mathbf r,t),
\]

the canonical response is `δn^μ=χ^{μν}δb^ν`. In the approved direct response
super-space

\[
 I=(a,L,M,i,\mu),
\]

define the Pauli transition amplitude with the measurement operator first:

\[
 T^\mu_{nm;I}(\mathbf k,\mathbf q)
 =\int d\Omega\,Y_{LM}^{\mathrm{code}*}(\hat r)
 \Psi^{P\dagger}_{n\mathbf k,a}(r_i,\hat r)
 \Gamma^\mu
 \Psi^P_{m,\mathbf k+\mathbf q,a}(r_i,\hat r).
\]

The Pauli superscript is mandatory for the current baseline. The exact
scalar-relativistic vertex is not certified. [CODE CONVENTION] [DEFERRED]

The retarded KS spectral expression is

\[
 \boxed{
 \chi^{\mu\nu}_{KS,IJ}(\mathbf q,\omega)
 =\frac1{N_k}\sum_{\mathbf k,nm}
 \frac{f_{n\mathbf k}-f_{m,\mathbf k+\mathbf q}}
 {\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta}
 T^\mu_{nm;I}(\mathbf k,\mathbf q)
 \bigl[T^\nu_{nm;J}(\mathbf k,\mathbf q)\bigr]^* .}
\]

The contract fixes all of the following:

* occupations are `f_nk−f_m,k+q`;
* the first state is `(n,k)` and the second is `(m,k+q)`;
* the measurement vertex is the unstarred first factor at `I`;
* the source vertex is the conjugated second factor at `J`;
* the denominator is `ω+ε_n−ε_m+iη`;
* `η>0` is retarded;
* the endpoint is exactly `k+q`, folded only by the reciprocal-lattice rule in
  §8;
* no occupation, denominator, or conjugation change is allowed in a circular
  specialization.

Using halved ladder vertices,

\[
 \chi^+_{KS,IJ}=2\frac1{N_k}\sum_{\mathbf k,nm}
 \frac{f_{n\mathbf k}^{\uparrow}-f_{m,\mathbf k+\mathbf q}^{\downarrow}}
 {\omega+\epsilon_{n\mathbf k}^{\uparrow}-epsilon_{m,\mathbf k+\mathbf q}^{\downarrow}+i\eta}
 T^+_{nm;I}\,(T^+_{nm;J})^*,
\]

where `T^+` uses `σ^+` and the sum is the up-to-down transition. Similarly,

\[
 \chi^-_{KS,IJ}=2\frac1{N_k}\sum_{mathbf k,nm}
 \frac{f_{n\mathbf k}^{\downarrow}-f_{m,\mathbf k+\mathbf q}^{\uparrow}}
 {\omega+\epsilon_{n\mathbf k}^{\downarrow}-\epsilon_{m,\mathbf k+\mathbf q}^{\uparrow}+i\eta}
 T^-_{nm;I}\,(T^-_{nm;J})^*.
\]

These are definitions only; no susceptibility is implemented here. [LITERATURE]
[DERIVED MAPPING] [DEFERRED]

## 8. Fourier convention, finite-q endpoint, and folding

The canonical RS reciprocal coordinates are fractional/direct coordinates:

\[
 \mathbf k_f=(k_1,k_2,k_3),
 \qquad
 \mathbf R_f=(R_1,R_2,R_3),
\]

with the live phase

\[
 \exp\!\left(+i2\pi\mathbf k_f\cdot\mathbf R_f\right).
\]

In Cartesian coordinates, if `a_i` are direct primitive vectors and `b_i` are
reciprocal primitive vectors satisfying `a_i·b_j=2πδ_ij`, then

\[
 \mathbf k=\sum_i k_i\mathbf b_i,
 \qquad
 \mathbf R=\sum_i R_i\mathbf a_i,
 \qquad
 e^{+i\mathbf k\cdot\mathbf R}=e^{+i2\pi\mathbf k_f\cdot\mathbf R_f}.
\]

For an endpoint block from site `a` to site `b`, the certified displacement is

\[
 \mathbf d_f=\mathbf R_f+\boldsymbol\tau_{b,f}-\boldsymbol\tau_{a,f},
 \qquad e^{+i2\pi\mathbf k_f\cdot\mathbf d_f}.
\]

The canonical response momentum `q_RS` is the endpoint-shift momentum in this
positive-phase convention. BES Appendix A writes the lattice response phase as
`exp(-i q_BES·R)`, so the exact mapping is

\[
 \boxed{\mathbf q_{BES}=-\mathbf q_{RS}}
\]

when the same cell-origin convention is used. The site endpoint phase is not
discarded: it is represented by `τ_b−τ_a` and the associated site gauge.
[CODE CONVENTION] [LITERATURE] [DERIVED MAPPING]

An arbitrary fractional reciprocal point is folded canonically as

\[
 \operatorname{fold}(x)=x-\lfloor x+\tfrac12\rfloor,
 \qquad \operatorname{fold}(x)\in[-\tfrac12,\tfrac12).
\]

For a supplied pair `(k,q)`, form the physical endpoint `k+q`, then fold that
endpoint by an integer reciprocal vector `G` only. Do not replace it by the
nearest mesh point. The arbitrary-k service and the response endpoint use the
same rule. [CODE CONVENTION]

Under `k→k+G`, the site gauge is

\[
 D_a(G)=e^{-i2\pi\mathbf G_f\cdot\boldsymbol\tau_{a,f}},
 \qquad
 H(k+G)=D(G)H(k)D(G)^\dagger,
 \qquad
 c_a(k+G)=D_a(G)c_a(k).
\]

The response vertex must transform with this same endpoint representative on
the state at `k+q`; the physical `T`, and hence `χ`, is gauge invariant after
the bra/ket endpoint phases are combined. No independent response phase may
be introduced. [CODE CONVENTION] [DERIVED MAPPING]

## 9. Radial, angular, and response-index normalization

### Angular convention

The live complex harmonics are ordered

```text
sp :  (0,0), (1,-1), (1,0), (1,1)
spd:  (0,0), (1,-1), (1,0), (1,1), (2,-2), (2,-1), (2,0), (2,1), (2,2)
```

and obey

\[
 Y^{\mathrm{code}}_{lm}=Y^{\mathrm{CS}*}_{lm},
 \qquad
 \int d\Omega\,Y_{LM}^{\mathrm{code}*}Y_{L'M'}^{\mathrm{code}}
 =\delta_{LL'}\delta_{MM'}.
\]

The product coefficient used by the response mapping is

\[
 (Y^{\mathrm{code}}_{lm})^*Y^{\mathrm{code}}_{l'm'}
 =\sum_{LM}{\cal G}^{LM}_{lm,l'm'}Y^{\mathrm{code}}_{LM},
\]
\[
 {cal G}^{LM}_{lm,l'm'}
 =(-1)^{m'}\sqrt{\frac{(2L+1)(2l+1)(2l'+1)}{4\pi}}
 \begin{pmatrix}L&l&l'\\0&0&0\end{pmatrix}
 \begin{pmatrix}L&l&l'\\M&m&-m'\end{pmatrix}.
\]

The complete product cutoff is

\[
 L_{\max}^{response}=2l_{\max};
 \qquad sp:2,\quad spd:4.
\]

The canonical direct response super-index is

\[
 I=(a,L,M,i,\mu),
\]

with site first, then `L=0…Lmax`, `M=-L…L`, radial point, and explicit density
or spin channel. The one-based flat index is

\[
 I_{flat}=\left((((a-1)(L_{max}+1)^2+h(L,M))N_r+(i-1))N_\mu+\mu\right),
 \quad h(L,M)=L^2+L+M.
\]

No Hermitian reduction or real-harmonic relabeling is permitted. [CODE
CONVENTION]

### Radial functions and measures

The production logarithmic mesh is

\[
 r_i=B[e^{A(i-1)}-1],
 \qquad
 \frac{dr}{di}=A(r_i+B).
\]

The physical spatial measure is

\[
 d^3r=r^2dr\,d\Omega,
 \qquad
 dr_i= A(r_i+B)\,di.
\]

The legacy stored radial density is instead

\[
 \mathrm{RHO}_\sigma(r_i)=4\pi r_i^2n_\sigma(r_i).
\]

Thus

\[
 N_\sigma=\int 4\pi r^2n_\sigma(r)dr
 =\int \mathrm{RHO}_\sigma(r)dr
 \simeq\sum_i w_iA(r_i+B)\mathrm{RHO}_\sigma(r_i),
\]

where the production composite-Simpson weights are `1/3,4/3,2/3,…,4/3,1/3`.
There is no second `r²` in the `RHO dr` contraction. Conversely, a physical
density contraction has `r² A(r+B)` and the angular integral. [CODE CONVENTION]

For a spherical quantity `q(r)`, its normalized `L=0` coefficient is

\[
 q_{00}(r)=\int d\Omega,Y_{00}^*q(r)=\sqrt{4\pi}\,q(r),
 \qquad Y_{00}=\frac1{\sqrt{4\pi}}.
\]

This gives the two equivalent charge forms

\[
 \int d\Omega,q(r)=4\pi q(r)=\sqrt{4\pi}\,q_{00}(r).
\]

The `4π` in `RHO` is a physical angular integral; the `Y00` coefficient is
normalized and contains `\sqrt{4π}`, not `4π`. These are not interchangeable
storage conventions. [DERIVED MAPPING]

For the Pauli-projected LMTO state,

\[
 \Psi^P_{n\mathbf k,a}(\mathbf r)=\frac1r
 \sum_{lm\sigma}c^{n\mathbf k}_{alm\sigma}
 U^{n\mathbf k}_{al\sigma}(r)Y^{\mathrm{code}}_{lm}(\hat r)|\sigma\rangle,
\]

with

\[
 U^{n\mathbf k}_{al\sigma}(r)=G_{al\sigma}(r)+
 (\epsilon_{n\mathbf k}-E^{work}_{\nu,al\sigma})\dot G_{al\sigma}(r).
\]

The direct radial transition value is therefore

\[
 T^{P,\mu}_{nm;aLM}(r)=\frac1{r^2}
 \sum_{lm,l'm',\sigma\sigma'}
 (c^{n\mathbf k}_{alm\sigma})^*c^{m,\mathbf k+\mathbf q}_{al'm'\sigma'}
 U^n_{al\sigma}(r)U^m_{al'\sigma'}(r)
 \Gamma^\mu_{\sigma\sigma'}{cal G}^{LM}_{lm,l'm'}.
\]

The `r²` in the physical volume measure cancels the explicit `1/r²` only in
a later volume contraction. It must not be removed from the pointwise direct
vertex. The scalar-relativistic `GFAC` and the packed lower component belong to
the certified spherical `NEWRHO` metric; they are not inserted into arbitrary
Pauli Gaunt products. [CODE CONVENTION] [DERIVED MAPPING]

### Direct values versus weighted objects

| object | canonical meaning | unit/contraction |
|---|---|---|
| `n`, `s`, `T(r_i)` | pointwise physical density/transition value | electrons bohr`^-3`; `T` has the same density unit |
| `RHO` | weighted legacy radial density | electrons bohr`^-1`; contract with `dr` |
| `dr/di` | log-mesh Jacobian | bohr |
| `r² dr` | radial part of physical volume | bohr`^3` |
| `χ(r,r')` | pointwise density response kernel to an energy source | electrons bohr`^-6` Ry`^-1`; both spatial measures are explicit |
| direct-mesh `χ_{IJ}` | same kernel sampled at two points | electrons bohr`^-6` Ry`^-1`; quadrature weights are not hidden |
| response contraction | an explicitly measured integral/sum | carries the declared quadrature measure; no bare array is a new physical unit |

No reduced radial basis is selected in LR-03. Any later contracted matrix must
publish its basis normalization and whether its quadrature weights are in the
matrix or in the contraction. [DEFERRED]

## 10. BES ALSDA mapping

BES defines the external coupling by

\[
 H_{ex}=\sum_i\int d^3x\,\hat n^i(x){\cal V}_i(x),
 \qquad
 {\cal V}_{x,y,z}=-\frac{g\mu_B}{2}B_{x,y,z},
\]

and uses Pauli density operators. Its transverse ALSDA relation is

\[
 K_{xc}^{BES}(x)=-\mu_B\frac{B_{xc}^{BES}(x)}{m(x)}.
\]

The `m(x)` in this equation is `s(x)`, while `B_xc^BES` is the physical field
whose energy coefficient is `-gμB B_xc^BES/2`. [LITERATURE]

Since

\[
 B_{xc}^{\sigma}=-\frac{g\mu_B}{2}B_{xc}^{BES},
\]

the exact algebraic RS mapping is

\[
 \boxed{
 K_{xc}^{RS}(r)=\frac{2}{g}\frac{B_{xc}^{\sigma}(r)}{s(r)}.}
\]

For BES/RS `g=2`,

\[
 \boxed{K_{xc}^{RS}(r)=\frac{B_{xc}^{\sigma}(r)}{s(r)}
 =\frac{\Delta V_{xc}(r)}{2s(r)}.}
\]

The sign is positive in terms of the up-minus-down RS Pauli coefficient because
the electron sign already occurred in `B_xc^σ=-μB B_xc^BES`. If one instead
uses the physical electron moment, the same relation is

\[
 K_{xc}^{RS}=-\frac{\mu_B B_{xc}^{\sigma}}{M_e}
 \quad (g=2),
\]

whereas with the positive code moment it is

\[
 K_{xc}^{RS}=+\frac{\mu_B B_{xc}^{\sigma}}{M^{code}}.
\]

No extra factor of `μB`, no extra factor of two, and no Ha-to-Ry conversion may
be applied after `B_xc^σ` is read from LR-01. [DERIVED MAPPING]

### SR versus Pauli baseline

LR-01's accepted XC arrays are generated from the accepted scalar-relativistic
radial densities. Define

\[
 s^{SR}=n^{SR}_\uparrow-n^{SR}_\downarrow,
 \qquad B_{xc}^{\sigma,SR}=\frac{\Delta V_{xc}^{SR}}2.
\]

This pair is **exactly matched to the ground-state XC functional**:

\[
 K_{xc}^{SR}=B_{xc}^{\sigma,SR}/s^{SR}\quad(g=2).
\]

The Pauli transition vertices and their occupied diagonal closure define

\[
 s^P=m^P=n^P_\uparrow-n^P_\downarrow.
\]

This is **exactly matched to the certified Pauli response basis**, but it is
not the density used to generate the stored LR-01 XC field. [NUMERICAL INPUT FROM
LR-02N]

The unique initial Pauli-response contract used by later work is therefore
explicitly the controlled mixed approximation

\[
 \boxed{K_{xc}^{P\leftarrow SR}(r)=
 \frac{B_{xc}^{\sigma,SR}(r)}{s^P(r)}\quad(g=2),}
\]

with Pauli vertices and Pauli rigid-rotation vector. It is not silently called
an exact functional derivative of the Pauli-projected density. The exact
functional pair and the response-basis pair remain visible in all reports.
At points where `s^P=0`, the formula is undefined; LR-03 invents no density
floor or regularization. A later kernel implementation must report `BLOCKED`
if its approved basis has nonzero weight on such a singular region. [DERIVED
MAPPING] [DEFERRED]

LR-02N gives the measured Fe scale for this approximation: the relative
volume-weighted magnetization L2 difference between `s^SR` and `s^P` is
`1.43312057e-2` over the full ASA sphere and `1.40811445e-2` over the
90%-absolute-magnetization prefix. The integrated spin numbers are
`2.1038582410` (SR) and `2.1487266935` (Pauli), a difference of
`-4.48684524e-2` electrons. These are projection discrepancies, not TDDFT or
Goldstone validation. [NUMERICAL INPUT FROM LR-02N]

An exact Pauli-functional contract would require evaluating the same XC
functional on a certified Pauli ground-state density, or exposing the exact
scalar-relativistic response density operator. Neither is silently inferred
here. [DEFERRED]

## 11. Lounis/Costa/Muniz/Mills static sum-rule mapping

The LCMM sum rule is written in terms of their effective splitting

\[
 B_{eff}^{LCMM}=V_\downarrow-V_\uparrow,
\]

and their radial magnetization `m_z`, which is the up-minus-down number density
in the Green-function identity. Their full spatial sum rule is

\[
 \sum_j\int d^3r'\,
 \chi^{ij}_0(\mathbf r,\mathbf r';0)B_{eff}^j(\mathbf r')
 =m_z^i(\mathbf r).
\]

For a spherical ASA field, their angular reduction is

\[
 \sum_j\int dr'\sum_{LL_1}
 \chi^{iLL_1;jL_1L}_0(r,r';0)B_{eff}^j(r')
 =4\pi m_z^i(r),
\]

and their local definition is

\[
 U^j(r)=\frac{B_{eff}^j(r)}{4\pi m_z^j(r)}.
\]

These equations are mapped conventions only. LR-03 does not construct `χ0`,
solve `ΓU=m`, or implement the sum-rule interaction. [LITERATURE] [DEFERRED]

For the canonical RS fields,

\[
 B_{eff}^{XC,RS}=V_{xc,\downarrow}-V_{xc,\uparrow}
 =-\Delta V_{xc}=-2B_{xc}^{\sigma},
\]

and, including the separately stored constraining field,

\[
 B_{eff}^{total,RS}=-\Delta V_{xc}-2B_{fsm}.
\]

No `μB` is inserted into this LCMM energy-valued splitting. If an external
physical field is included, its energy coefficient must first be formed with
`-gμB B/2`, and only then may it be added to a spin-channel potential.
[DERIVED MAPPING]

### Origin of the published `4π`

LCMM's `4π` is from integrating the spherical right-hand side over the output
solid angle:

\[
 \int d\Omega\,m_z(r)=4\pi m_z(r).
\]

In the normalized complex response basis,

\[
 m_{z,00}(r)=\sqrt{4\pi}\,m_z(r),
 \qquad
 B_{eff,00}(r)=\sqrt{4\pi}\,B_{eff}(r).
\]

Thus a normalized `L=0` coefficient equation has `m_{z,00}`, not `4πm_z`,
on its right-hand side. The LCMM radial equation has `4πm_z` because it uses
the angular-integrated radial convention. If `U` is imported into the
normalized RS coefficients, the literal relation is

\[
 U_{LCMM}=\frac{B_{eff}}{4\pi m_z}
 =\frac{B_{eff,00}/m_{z,00}}{4\pi}.
\]

This is the required normalization map. A later GSR implementation must pick
one representation and retain its measure; it must not add or remove `4π` by
numerical comparison. [DERIVED MAPPING]

## 12. Goldstone rigid-rotation vector

For a spherical collinear Pauli density matrix,

\[
 \rho_a^P(r)=\frac12\left[n_a^P(r)I+s_a^P(r)\sigma_z\right].
\]

For a global infinitesimal rotation vector
`\boldsymbol\theta=(\theta_x,\theta_y,0)`,

\[
 \delta\mathbf s_a^P(r)=\boldsymbol\theta\times
 [s_a^P(r)\hat{\mathbf z}],
\]
\[
 \delta s_{x,a}^P=\theta_y s_a^P,
 \qquad
 \delta s_{y,a}^P=-\theta_x s_a^P,
 \qquad
 \delta s_a^{P,+}=-i(\theta_x+i\theta_y)s_a^P.
\]

In the normalized response super-space, the Pauli rigid-rotation vector is

\[
 g^{P,+}_{a,L,M,i}
 =-i\,\theta^+\,\delta_{L0}\delta_{M0}
 \sqrt{4\pi}\,s_a^P(r_i),
 \qquad \theta^+=\theta_x+i\theta_y.
\]

The `−` channel is the conjugate rotation,

\[
 g^{P,-}_{a,L,M,i}=+i\,\theta^-\,\delta_{L0}\delta_{M0}
 \sqrt{4\pi}\,s_a^P(r_i).
\]

There is one common rotation angle per site for a rigid global rotation; the
site dependence is only the certified radial profile `s_a^P(r)`. For the SR
ground-state diagnostic, replace `s_a^P` by `s_a^SR`; the two vectors are not
silently identified. A physical-moment vector is obtained by multiplying by
`-gμB/2`, and a code-moment vector by `+μB`. [DERIVED MAPPING]

No zero eigenvalue is imposed, no matrix is modified, and no Goldstone
correction is applied in LR-03. [DEFERRED]

## 13. Symmetry and covariance identities

The following statements apply only to the certified no-SOC collinear baseline,
with Hermitian Pauli densities and the response ordering of §6. They are
identities or conditional covariances, not material-spectrum checks.

### Retarded/advanced and frequency covariance

For Hermitian density operators, in a fixed compatible endpoint gauge, the
spectral definition gives

\[
 \boxed{
 \chi_R^{\mu\nu}{}_{IJ}(\mathbf q,\omega)
 =\left[
 \chi_R^{\mu\nu}{}_{IJ}(-\mathbf q,-\omega^*)
 \right]^*.}
\]

Equivalently, at real `ω`, complex conjugation changes the retarded `+iη`
prescription to the advanced `−iη` prescription, with `q→−q`. The operator
order remains measurement first/source second; it is not silently exchanged.
This is the safe retarded/advanced statement; a bare `χ(q,ω)^*=χ(q,ω)` is
false in general. [DERIVED MAPPING]

In a collinear axially symmetric state this becomes

\[
 \boxed{
 \chi_R^+(\mathbf q,\omega)=
 \left[\chi_R^-(-\mathbf q,-\omega^*)\right]^*,
 \qquad
 \chi_R^-(\mathbf q,\omega)=
 \left[\chi_R^+(-\mathbf q,-\omega^*)\right]^*.}
\]

Thus the `+` and `−` channels exchange under the combined `q→−q`,
`ω→−ω`, complex-conjugation covariance. The exchange is not a statement that
the two channels have the same positive-frequency spectrum. [DERIVED MAPPING]

### `q→−q`

There is no standalone identity `χ^±(q,ω)=χ^±(−q,ω)` for an arbitrary
multi-site basis. It requires an additional spatial inversion or a certified
site-permuting symmetry, including its endpoint gauge. If a symmetry `P`
exists with `P H(k)P^{-1}=H(-k)` and maps response indices, then

\[
 \chi^\pm_{IJ}(\mathbf q,\omega)
 =\chi^\pm_{P(I)P(J)}(-\mathbf q,\omega).
\]

That conditional identity does not exchange circular channels. Without the
additional spatial symmetry, use only the retarded covariance above. [DERIVED
MAPPING] [DEFERRED]

### Global spin reversal

Use the explicit spin reversal rotation `S=iσ_y` (an irrelevant overall phase
is immaterial):

\[
 S\sigma_xS^\dagger=-\sigma_x,
 \quad S\sigma_yS^\dagger=+\sigma_y,
 \quad S\sigma_zS^\dagger=-\sigma_z,
\]
\[
 S\sigma^+S^\dagger=-\sigma^-.
\]

The reversed collinear state has `s_z→−s_z`,
`V_up↔V_down`, `ΔV_xc→−ΔV_xc`, and `B_xc^σ→−B_xc^σ`. The two minus signs
from the measurement/source circular operators cancel, so

\[
 \boxed{
 \chi^{+}_{[\mathrm{reversed}]}(\mathbf q,\omega)
 =\chi^{-}_{[\mathrm{original}]}(\mathbf q,\omega),}
\]

with the analogous `+↔−` relation for the other channel. The ALSDA ratio
`B_xc^σ/s_z` is invariant under this simultaneous reversal. [ALGEBRAIC TEST]
[DERIVED MAPPING]

## 14. Units

The canonical production energy unit is Ry. The response source `b` is an
energy-valued coefficient multiplying a Pauli density operator.

| quantity | canonical unit | locked conversion/meaning |
|---|---|---|
| electronic `ε`, `ω`, `η` | Ry | retarded denominator uses `+iη`, `η>0` |
| optional displayed energy | eV | `1 Ry = 13.605693122994 eV`; conversion is output-only |
| Hartree | Ha | `1 Ha=2 Ry=27.211386245988 eV` |
| radial `r`, `R`, `B` | bohr | not Å; lattice input units do not change this radial contract |
| `n`, `s` | electrons bohr`^-3` | number densities |
| `RHO` | electrons bohr`^-1` | `RHO=4πr²n`; `RHO dr` integrates to electrons |
| `M_e`, `M^code` | `μB` bohr`^-3` | `M_e=−gμB s/2`; `M^code=+μB s` |
| `V_xc,σ`, `V_xc^0`, `ΔV_xc`, `B_xc^σ`, `B_fsm` | Ry | multiplicative energy coefficients |
| LCMM `B_eff` | Ry | `V_down−V_up`, not a tesla field |
| physical `B` | tesla at an SI boundary | first form `b_B=−gμB B/2`, then convert `b_B` to Ry |
| `χ(r,r')` | electrons bohr`^-6` Ry`^-1` | pointwise density kernel; spatial measures are explicit |
| direct `T(r)` | electrons bohr`^-3` | unintegrated transition density |
| `K_xc` | Ry bohr`^3` | `δb/δs`; for `g=2`, `B_xc^σ/s` |
| normalized `L=0` coefficient | same density unit as the field | coefficient is `√(4π)` times a spherical radial value |

The BES paper states that Rydberg atomic units are used. In ordinary Hartree
atomic units, `e=\hbar=m_e=4πε_0=1` and `μB=1/2` in energy per atomic magnetic
field. If only the energy unit is changed from Ha to Ry while the atomic
magnetic-field unit is retained, the numerical value of `μB` is `1` in
Ry per atomic field. These numerical simplifications do **not** remove the
electron minus sign, the explicit `g/2`, `ΔV=2B^σ`, ladder-operator factors,
or the `4π`/`Y00` distinction. With tesla input, retain physical constants and
convert to Ry explicitly. [LITERATURE] [DERIVED MAPPING]

## 15. Independent algebraic falsification

The standalone script [`tests/validation/val24_lr03_conventions.py`](../tests/validation/val24_lr03_conventions.py)
uses only Python's standard library and explicitly constructs all matrices. It
does not import response production code, ground-state data, Fe/Ni data, or a
susceptibility implementation. It checks:

1. `[σx,σy]=2iσz`;
2. the exact halved circular matrices in the certified `(up,down)` order;
3. `Σ^±=2σ^±` and the spin-reversal map `σ^+→−σ^-`, which exchanges the two
   circular channels;
4. reconstruction of arbitrary unequal `V_up,V_down` from `V^0 I+B^σσ_z`.

The expected result is `PASS` for algebraic conventions only. It is not a
physical response validation. [ALGEBRAIC TEST]

## 16. PASS boundary and deferred implementation

The following are locked with no hidden choice:

* density definition and actual RS spinor order;
* Pauli matrices and halved/unhalved circular operators;
* number-spin versus physical electron moment and positive code moment;
* `V_xc^0`, `B_xc^σ`, `ΔV_xc`, BES `B_xc`, and LCMM `B_eff`;
* BES circular source/measurement ordering and published factor `2`;
* retarded occupation difference, endpoint, denominator, `+iη`, and conjugation;
* RS positive Fourier phase, BES momentum sign map, folding, and site gauge;
* harmonic/Gaunt normalization, `L_max=2l_max`, radial measure, `RHO`, and
  `Y00` normalization;
* BES ALSDA mapping and Ry/Ha/μB dimensions;
* LCMM `4π` origin and normalized-basis conversion;
* Pauli rigid-rotation vector and covariance identities;
* the explicit SR→Pauli approximation `K^{P←SR}=B_xc^{σ,SR}/s^P`.

The only approximation in the certified baseline is therefore named,
quantified, and visible. It must not be replaced by whichever choice gives a
smaller later Goldstone residual. An exact SR response would require the
missing scalar-relativistic nonspherical density operator; an exact Pauli
functional response would require XC evaluation on the Pauli-projected
density. Both are deferred. [NUMERICAL INPUT FROM LR-02N] [DEFERRED]

No production response code is authorized by this document. Later code must
cite this file and consume its canonical objects rather than duplicate old
response conventions. [DEFERRED]

## Checklist

- [x] LR-02N numerical evidence consumed.
- [x] Density operators fixed.
- [x] RS spinor ordering and Pauli matrices fixed.
- [x] Number-spin and magnetic-moment densities separated.
- [x] Circular operators and every factor of two derived.
- [x] Circular susceptibility ordering fixed.
- [x] Retarded denominator, occupation difference, endpoint, and conjugation fixed.
- [x] XC scalar/splitting/Pauli-field dictionary fixed.
- [x] BES ALSDA mapping derived but not implemented.
- [x] Lounis sum-rule normalization mapped but not implemented.
- [x] `m^SR` versus `m^P` kept explicit; mixed approximation quantified.
- [x] Radial `4π`, `r²`, `dr/di`, and `Y00` factors mapped.
- [x] Fourier convention and finite-q folding/site gauge imported.
- [x] Units fixed.
- [x] Rigid-rotation vector defined without enforcing a zero mode.
- [x] Covariance identities derived with assumptions stated.
- [x] Independent algebraic falsification test added.
- [x] No sign/factor chosen from Fe/Ni behavior.
- [x] No susceptibility, kernel, Dyson, or Goldstone implementation.
- [x] PASS declared with an explicit controlled SR→Pauli approximation.
