# DRESP-02F — GF/Lehmann Formulation Equivalence Audit

## Status and verdict

**Verdict:** `PASS — CURRENT GF FORMULATION PROVEN EQUIVALENT`

The reciprocal real-axis GF/Kubo expression is analytically equivalent to the
finite-spectrum Lehmann response under the repository convention.  The
equivalence is exact in the controlled limit

\[
 \eta_{\rm int}\to0^+,
 \qquad E_{\min}\to-\infty,
 \qquad E_{\max}\to+\infty,
 \qquad h/\eta_{\rm int}\to0,
\]

at fixed positive response broadening \(\eta_{\rm response}\).  At the finite
`integration_eta` used by the production GF route, equality to the Lehmann
answer is not expected: the spectral delta functions have deliberately been
replaced by Lorentzians.  For the present Kubo form this adds

\[
 \eta_{\rm int}
\]

to the response pole width, giving \(\eta_{\rm response}+\eta_{\rm int}\)
for an isolated transition when the Fermi factor is locally constant.  The
finite-window and Simpson errors are additional, separate errors.

No source correction was made.  In particular, no sign, spin-ordering,
occupation, conjugation, q-endpoint, or LMTO-GF representation defect was
identified by this audit.

This audit concerns the bare Kohn-Sham response only.  It does not use or
evaluate an XC kernel, Goldstone condition, ALSDA, GSR, Mills (U), Dyson
poles, or an interacting susceptibility.

## 1. Repository contract

For a fixed k point, write

\[
 L=(\mathbf k,\uparrow),\qquad R=(\mathbf k+\mathbf q,\downarrow),
\]

with

\[
 H_L|nL\rangle=\epsilon_{nL}|nL\rangle,
 \qquad
 H_R|mR\rangle=\epsilon_{mR}|mR\rangle.
\]

The plus-channel vertex is

\[
 T_i^{nm}(\mathbf k,\mathbf q)
 =\langle nL|V_i^+|mR\rangle .
\]

The source encodes \(\sigma_+=(\sigma_x+i\sigma_y)/2\) as the spin matrix
with element `(up,down)=1`; the minus channel has `(down,up)=1`.
The corresponding source definitions are
`source/lr_pauli_transition_vertex.f90:102-116`.

The common normalization used by both backends is

\[
 C_k=\frac{2w_k}{\sum_{k'}w_{k'}}.
\]

The factor 2 is a repository response normalization and is not part of the
GF-versus-Lehmann question.  The GF result is converted to the same LR-04
canonical right-weighted representation after accumulation.

The relevant implementation paths are:

| object | implementation |
|---|---|
| Lehmann response | `source/lr_ks_susceptibility.f90:371-446` |
| real-axis GF response | `source/lr_gf_susceptibility.f90:53-160` |
| GF Kubo contractions | `source/lr_gf_susceptibility.f90:195-246` |
| compact independent GF oracle | `source/lr_product_gf_susceptibility.f90` |
| projected DRESP-02 GF route | `source/lr_projected_reciprocal_chi0.f90:240-470` |

The exact folded q endpoint, occupation provenance, EF, temperature, energy
zero, dimensions, and baseline representation are checked by
`validate_gf_inputs` (`source/lr_gf_susceptibility.f90:360-432`).

## 2. Retarded correlator and Lehmann derivation

Use the repository convention

\[
 \chi^{+-}_{ij}(t)
 =-i\theta(t)\langle[S_i^+(t),S_j^-(0)]\rangle,
 \qquad S_j^-=(S_j^+)^\dagger .
\]

For a general retarded correlator, with eigenstates \(|a\rangle\), energies
\(E_a\), and statistical weights \(p_a\),

\[
 \chi^R_{AB}(\omega)
 =\sum_{ab}
 \frac{p_a-p_b}{\omega+E_a-E_b+i0^+}
 \langle a|A|b\rangle\langle b|B|a\rangle .
\]

The sign follows directly from

\[
 -i\int_0^\infty dt\,
 e^{i(\omega+E_a-E_b)t-0^+t}
 =\frac{1}{\omega+E_a-E_b+i0^+}.
\]

For a one-body transverse operator, a particle-hole transition takes an
occupied \(nL\) state to an \(mR\) state.  Wick contraction gives

\[
 p_a-p_b\longrightarrow f_{nL}-f_{mR},
\]

and

\[
 \langle a|S_i^+|b\rangle
 \langle b|S_j^-|a\rangle
 =T_i^{nm}T_j^{nm*}.
\]

Therefore

\[
 \boxed{
 \chi^{+-}_{0,ij}(\mathbf q,\omega)=
 \sum_{\mathbf k,n,m} C_k
 \frac{f_{nL}-f_{mR}}
 {\omega+\epsilon_{nL}-\epsilon_{mR}+i\eta_{\rm response}}
 T_i^{nm}T_j^{nm*}}
 \]

which is the repository expression.  In particular:

* the denominator is `frequency + left energy - right energy + i eta`;
* the numerator is left/up occupation minus right/down occupation;
* the left endpoint is \(\mathbf k,\uparrow\);
* the right endpoint is the exact folded \(\mathbf k+\mathbf q,\downarrow\);
* the response matrix element is \(T_iT_j^*\), not (T_iT_j) or
  (T_i^*T_j).

The source accumulation at `lr_ks_susceptibility.f90:408-420` maps to these
factors one-for-one.

## 3. GF representation derived from the same spectrum

Define

\[
 G_{L}^{R/A}(E)=
 \sum_n\frac{|nL\rangle\langle nL|}
 {E-\epsilon_{nL}\pm i0^+},
\]

and the normalized spectral function

\[
 A_L(E)=\frac{i}{2\pi}\left[G_L^R(E)-G_L^A(E)\right]
 =\sum_n|nL\rangle\langle nL|\delta(E-\epsilon_{nL}).
\]

For a positive-frequency transverse response, the exact real-axis expression
with the source's operator ordering is

\[
\begin{aligned}
 \chi^{(1)}_{ij}(\omega)
 &=\int_{-\infty}^{\infty}dE\,f(E)\,
 \operatorname{Tr}\left[A_L(E)V_iG_R^R(E+\omega)V_j^\dagger\right],\\
 \chi^{(2)}_{ij}(\omega)
 &=\int_{-\infty}^{\infty}dE\,f(E)\,
 \operatorname{Tr}\left[A_R(E)V_j^\dagger G_L^A(E-\omega)V_i\right],\\
 \chi_{ij}(\omega)&=\chi^{(1)}_{ij}(\omega)+\chi^{(2)}_{ij}(\omega).
\end{aligned}
\]

The k sum and common factor \(C_k\) are understood.  Insert the spectral
representation into the first term:

\[
 \chi^{(1)}_{ij}
 =\sum_{nm}
 \frac{f(\epsilon_{nL})T_i^{nm}T_j^{nm*}}
 {\epsilon_{nL}+\omega-\epsilon_{mR}+i\eta_{\rm response}}.
\]

The second term gives

\[
 \chi^{(2)}_{ij}
 =\sum_{nm}
 \frac{f(\epsilon_{mR})T_i^{nm}T_j^{nm*}}
 {\epsilon_{mR}-\omega-\epsilon_{nL}-i\eta_{\rm response}}.
\]

Since

\[
 \epsilon_{mR}-\omega-\epsilon_{nL}-i\eta_{m response}
 =-left(\omega+\epsilon_{nL}-\epsilon_{mR}+i\eta_{m response}\right),
\]

the sum is

\[
 \chi^{(1)}_{ij}+\chi^{(2)}_{ij}
 =\sum_{nm}
 \frac{f(\epsilon_{nL})-f(\epsilon_{mR})}
 {\omega+\epsilon_{nL}-\epsilon_{mR}+i\eta_{\rm response}}
 T_i^{nm}T_j^{nm*}.
\]

This is exactly the Lehmann result.  The proof is finite-spectrum and does
not use a material approximation.

The source implements precisely these two terms:

* first term: `left_a(E)`, `right_gr(E+omega)`, and
  `sum(temporary*conjg(vertex_j))`;
* second term: `right_a(E)`, `left_ga(E-omega)`, and
  `conjg(transpose(vertex_j))` followed by `transpose(vertex_i)`.

The `conjg(transpose(...))` is the matrix Hermitian adjoint.  The final
elementwise sums are trace contractions, so the source multiplication order
is algebraically correct for complex vertices.

The result is exact at finite temperature as well as at zero temperature,
provided the real-energy integral is over the full line, the spectral delta
functions are exact, and the same Fermi function is used in both terms.

## 4. Finite-temperature audit

The exact finite-temperature identity is not a zero-temperature shortcut:

\[
 \int dE\,f_T(E)\delta(E-\epsilon_n)=f_T(\epsilon_n).
\]

Consequently the two terms reconstruct

\[
 f_T(\epsilon_{nL})-f_T(\epsilon_{mR})
\]

for arbitrary finite (T).  The source evaluates the same scalar

\[
 f_T(E)=\frac{1}{e^{(E-E_F)/(k_BT)}+1}
\]

for both Kubo terms at the unshifted integration variable; see
`lr_gf_susceptibility.f90:123-139`.  The Lehmann route stores the corresponding
endpoint occupations explicitly and forms their difference at
`lr_ks_susceptibility.f90:408-417`.

There is an equivalent shifted representation, but the Fermi factor must be
shifted with the variable.  For example, in the second term let

\[
 x=E-\omega.
\]

Then

\[
 \int dE,f(E)A_R(E)V_j^\dagger G_L^A(E-\omega)V_i
 =\int dx,f(x+\omega)A_R(x+\omega)V_j^\dagger G_L^A(x)V_i.
\]

Thus `f(E)` in both source terms is correct in the source's unshifted form;
using `f(E-omega)` or `f(E+omega)` without the corresponding variable change
would be wrong.  There is no missing finite-temperature contour term in this
representation.  A contour derivation may distribute terms differently, but
after the contour/variable transformation it must reduce to the same two
real-axis terms.

At `T=0`, the source uses the numerically stabilized Fermi helper with a
`kT` floor (`lr_ks_susceptibility.f90:155-167`).  Away from a level exactly at
EF this has the expected step-function limit.  At 300 K it is the ordinary
finite-temperature Fermi function in the material run.

## 5. Broadening taxonomy and exact single-transition analysis

There are three different widths:

1. `eta_response` is the physical retarded pole regularization in
   \(G_R^R(E+\omega)\) and \(G_L^A(E-\omega)\).
2. `integration_eta` is the numerical width used in
   \(A(E)=i(G^R-G^A)/(2\pi)\).
3. A direct product of two broadened spectral functions has a convolution
   width; for identical Lorentzians it is (2\eta_{m int}).

The current Kubo implementation does **not** contain two spectral functions in
one term.  It contains one broadened (A) and one resolvent whose width is
`eta_response`.  That distinction matters.

Let

\[
 \delta_\gamma(x)=\frac{1}{\pi}\frac{\gamma}{x^2+\gamma^2},
 \qquad \gamma=\eta_{\rm int},
\]

and consider one transition with

\[
 \Delta=\epsilon_m-\epsilon_n,
 \qquad x=\omega-\Delta.
\]

The source's scalar bubble is

\[
\begin{aligned}
 \chi_\gamma(\omega)=|T|^2\int dE\,f(E)\bigg[&
 \frac{\delta_\gamma(E-\epsilon_n)}{E+\omega-\epsilon_m+i\eta_r}\\
 &+\frac{\delta_\gamma(E-\epsilon_m)}{E-\omega-\epsilon_n-i\eta_r}\bigg],
\end{aligned}
\]

where \(\eta_r=\eta_{\rm response}\).

If (f) is locally constant over the Lorentzian support, residue calculus
gives

\[
 \int_{-\infty}^{\infty}dE\,
 \frac{\delta_\gamma(E-a)}{E-b+i\eta_r}
 =\frac{1}{a-b+i(\gamma+\eta_r)},
\]

and the advanced counterpart gives the corresponding negative term.  Thus

\[
 \boxed{
 \chi_\gamma(\omega)=
 \frac{(f_n-f_m)|T|^2}
 {\omega-\Delta+i(\eta_r+\gamma)}}
 \]

in that locally constant-Fermi limit.  At resonance,

\[
 -\operatorname{Im}\chi_\gamma(\Delta)
 =\frac{(f_n-f_m)|T|^2}{\eta_r+\gamma},
\]

whereas the Lehmann value is the same expression with \(\eta_r\) only.  The
leading fixed-η response error is therefore (O(\gamma/\eta_r)), not an
unknown normalization factor.

For a varying Fermi function, define the broadened occupation

\[
 F_\gamma(a)=\int dE\,f(E)\delta_\gamma(E-a).
\]

The integrated spectral weight of the current Kubo form is

\[
 -\frac{1}{\pi}\int_{-\infty}^{\infty}d\omega\,
 \operatorname{Im}\chi_\gamma(\omega)
 =|T|^2\left[F_\gamma(\epsilon_n)-F_\gamma(\epsilon_m)\right],
\]

not exactly \(|T|^2(f_n-f_m)\) at nonzero γ.  At zero temperature,

\[
 F_\gamma(a)=\frac12+\frac1\pi
 \arctan\frac{E_F-a}{\gamma},
\]

so the finite-γ integrated-weight error is analytically predictable.  This
is the occupation-smearing component of the residual, distinct from pole
width broadening.

For comparison, a direct double-spectral construction

\[
 \int dE\,A_{L,\gamma}(E)A_{R,\gamma}(E+\omega)
\]

has a Cauchy convolution width (2\gamma).  That is not the width of the
current one-(A\)-one-(G\) Kubo terms.  If the dynamic resolvent were also
given the integration width, the current formula would similarly produce
(2\gamma); it is not: the source uses `integration_eta` only in `left_a`/
`right_a` and uses `request%eta` in the shifted retarded/advanced resolvents.

### Isolated numerical transition

The following independent scalar evaluation uses

\[
 \epsilon_n=-0.10,\quad \epsilon_m=0.20,\quad E_F=0,
 \quad k_BT=0.02,\quad \eta_r=0.01,
 \quad |T|^2=1.
\]

The exact Lehmann resonance is

\[
 f_n=0.9933071491,\quad f_m=0.0000453979,
 \quad \chi_L(\Delta)=-i,99.32617512.
\]

The finite-γ current integral was evaluated directly over the real line:

| γ | GF at ω=Δ | Lehmann relative error |
|---:|---:|---:|
| 0.0400 | −1.09762 − 19.76009i | 0.80113 |
| 0.0200 | −0.69008 − 33.02949i | 0.66750 |
| 0.0100 | −0.41279 − 49.60853i | 0.50057 |
| 0.0050 | −0.23666 − 66.18304i | 0.33369 |
| 0.0025 | −0.12990 − 79.44094i | 0.20021 |
| 0.0010 | −0.05567 − 90.28762i | 0.09100 |
| 0.0005 | −0.02856 − 94.59172i | 0.04767 |
| 0.0001 | −0.00584 − 98.34179i | 0.00991 |

The peak stays at \(\omega=\Delta\) for this isolated transition.  In the
locally constant-Fermi analytic limit its half-width at half-maximum is

\[
 \Gamma_{1/2}=\eta_r+\gamma,
\]

and its integrated weight is unchanged.  With the finite-temperature Fermi
factor above, the numerical integrated weights are 0.94039 for
\(\gamma=0.01\) and 0.97986 for \(\gamma=0.0025\), approaching the exact
\(f_n-f_m=0.99326\) as γ decreases.  Both the width and weight trends are
therefore the expected spectral-regularization law.

## 6. Audit of both Kubo terms

For one (n,m) pair define

\[
 D_{nm}(\omega)=\omega+\epsilon_{nL}-\epsilon_{mR}+i\eta_r.
\]

The term-by-term spectral insertion is:

| source term | spectral contribution | role |
|---|---|---|
| \(\chi^{(1)}\) | \(f_{nL}T_iT_j^*/D_{nm}\) | occupied-left/up contribution and retarded right/down propagation |
| \(\chi^{(2)}\) | \(-f_{mR}T_iT_j^*/D_{nm}\) | right/down occupation counterpart and advanced left/up propagation |
| sum | \((f_{nL}-f_{mR})T_iT_j^*/D_{nm}\) | Lehmann transition |

The second term is sometimes described as the advanced or negative-frequency
partner.  Algebraically, in the positive-frequency representation used here,
it supplies the (-f_{mR}) part of the same Lehmann denominator.  The complete
sum over (n,m) still contains the reverse transitions and their negative
frequency structure.

The source labels and evaluates the terms independently in
`accumulate_gf_bubble`:

* `left_a(E) V_i right_gr(E+omega) V_j^dagger` is χ(^{(1)});
* `right_a(E) V_j^dagger left_ga(E-omega) V_i` is χ(^{(2)}).

The optimized projected backend retains these two arrays independently as
`kubo_term_one` and `kubo_term_two`.  The existing DRESP-02P Fe rerun found an
internal term-sum residual of (2.07\times10^{-12}) in the fine diagnostic
sample.  Thus the material residual is not caused by omitting or double
counting one Kubo term.

## 7. Fermi placement and real-energy window

The exact derivation is over ((-∞,\infty)).  The source uses

\[
 E_{\min}=\min(\epsilon_L,\epsilon_R)-\texttt{energy\_margin},
 \qquad
 E_{\max}=\max(\epsilon_L,\epsilon_R)+\texttt{energy\_margin},
\]

with composite Simpson quadrature.  This is a numerical truncation, not part
of the exact GF identity.

At finite `integration_eta`, Lorentzian tails extend beyond both endpoints;
the omitted contribution is therefore nonzero even when every pole lies
inside the nominal window.  A window increase can reduce this error, but it
does not remove the \(\eta_{\rm int}\)-induced width and occupation changes.

The DRESP-02R 300 K diagnostic gives the relevant separation at the material
sample: the fine-window Fermi-weighted spectral residual was (7.90\times10^{-2}),
while the window ladder changed the GF/Lehmann residual only modestly.  The
controlled `eta_int` ladder was the dominant monotonic trend.

## 8. Zero-temperature and limit-order audits

The zero-temperature unit fixture in
`tests/unit/test_lr_gf_susceptibility.f90` uses exact 0/1 occupations and a
finite spectrum.  It reports:

| check | result |
|---|---:|
| GF coarse full-response difference, `integration_eta=0.010` | (2.1443\times10^{-9}) absolute, 2.4633% relative |
| GF fine full-response difference, `integration_eta=0.001` | (3.6257\times10^{-9}) absolute, 4.1652% relative |
| isolated finite-spectrum analytic matrix element | (1.1277\times10^{-9}) absolute |
| static q=0 maximum difference | (9.2484\times10^{-10}) |

The small absolute values in that unit fixture do not imply zero relative
error.  The isolated analytic element validates the sign and denominator; the
remaining finite-frequency relative residual is consistent with the finite
spectral width and quadrature.

The safe hierarchy is:

1. Hold \(\eta_r>0\) fixed, decrease \(\eta_{\rm int}\), and resolve it with
   (h/\eta_{\rm int}\ll1).
2. Increase the energy window until Lorentzian-tail and Fermi-weighted
   moment errors are negligible.
3. Only after that take \(\eta_r\to0^+\), if a distributional real-axis
   response is required.

Taking \(\eta_r\to0\) at fixed nonzero \(\eta_{\rm int}\) leaves the
spectral broadening in the response.  Taking both to zero with a pole exactly
on the real axis is not a controlled numerical operation.

## 9. Direct spectral-function oracle

The compact product GF route provides an independent intermediate oracle.  It
forms scalar spectral factors explicitly,

\[
 A_{n,\gamma}(E)=\frac{i}{2\pi}\left[
 \frac{1}{E-\epsilon_n+i\gamma}
 -\frac{1}{E-\epsilon_n-i\gamma}\right],
\]

and contracts them with independently formed transition amplitudes.  It does
not call the Lehmann denominator accumulator.  Its scalar, optimized, and
factorized contraction backends agree at approximately machine precision in
`UnitLrProductGfSusceptibility`.

The mixed-eigenvector audit additionally verifies that the component vertex
tensor and eigenbasis transition amplitudes agree to (3.21\times10^{-16})
relative error for the plus channel and (1.44\times10^{-16}) for the minus
channel.  This separates the GF energy algebra from the LMTO affine-vertex
construction.

## 10. Complex-arithmetic, q, and representation audit

The source and fixtures establish the following:

* `build_weighted_resolvent` forms
  \(U\,\mathrm{diag}(\epsilon_n^p/(z-\epsilon_n))U^\dagger\), with the
  required `conjg(eigenvector)` on the right endpoint;
* the plus vertex maps right/down to left/up, and the minus vertex reverses
  this ordering;
* the response outer product is (T_iT_j^*\);
* the advanced endpoint is (G_L^A(E-\omega)), not a retarded or
  (E+\omega) endpoint;
* the exact folded q endpoint is checked rather than replaced by a nearest
  mesh point;
* the mixed fixture has nonzero imaginary off-diagonal Hamiltonian elements
  from (2.5\times10^{-3}) to (1.58\times10^{-2}), so a real/diagonal
  accidental pass is excluded;
* the orbital-rotation oracle gives a Hamiltonian eigenpair residual of
  (6.66\times10^{-16}), while GF and Lehmann responses are each invariant
  under the common unitary rotation to approximately (10^{-15}) relative
  error.

No transpose/Hermitian-transpose, imaginary-sign, spin, or q-order defect was
found.

## 11. Finite-model numerical closure

The primary finite fixture is
`tests/unit/test_lr_product_gf_quadrature_audit.f90` in `mixed` mode.  It has

* an 8-dimensional one-electron basis and 8 nondegenerate levels per endpoint;
* two k points;
* dense Hermitian Hamiltonians with genuinely complex off-diagonal entries;
* complex eigenvectors from numerical diagonalization;
* finite temperature and nontrivial LMTO radial/site/operator vertices;
* plus and minus channels and the q/−q endpoint path.

Its controlled ladder holds the Simpson resolution at

\[
 h/\eta_{\rm int}=0.435981420
\]

and uses \(\eta_{\rm response}=0.04\).  The GF-to-Lehmann Frobenius residuals
were:

| \(\eta_{\rm int}\) | \(\omega=0\) relative residual | \(\omega=0.17\) relative residual |
|---:|---:|---:|
| 0.0100 | 6.2557e−2 | 1.6778e−1 |
| 0.0050 | 3.3211e−2 | 9.3138e−2 |
| 0.0025 | 1.7115e−2 | 4.9181e−2 |

The observed approximately linear decrease is the predicted fixed-η
spectral-broadening law.  The three contraction backends agree independently;
the unitary-rotation oracle agrees with both production paths.  This is the
required finite complex Hamiltonian evidence that the Kubo representation,
rather than only a real diagonal fixture, is closing toward Lehmann.

## 12. Fe return and material interpretation

The warranted Fe rerun was already completed by DRESP-02P using the optimized
projected GF backend and the frozen accepted 4x4x4 state at 300 K.  No source
change occurred between that rerun and this formulation audit, so a duplicate
material run would not add evidence.

The controlled `eta_int` ladder was:

| selector | \(\eta_{\rm int}=0.008\) | 0.004 | 0.002 |
|---|---:|---:|---:|
| `d` | 1.35830 | 0.652224 | 0.320268 |
| `spd` | 3.43760 | 2.03628 | 1.17844 |

These are GF−Lehmann Frobenius differences with the mesh selected to keep
`h/eta_int` near 0.40.  The residual decreases as the derived law predicts;
it has not yet been extrapolated to zero at the material resolution.  The
window ladder changed the fine `d` difference only from 0.320680 to 0.313493
and the `spd` difference from 1.17120 to 1.16375 over margins 0.30 to 1.00 Ry.
The material term-sum residual was (2.07\times10^{-12}), so internal Kubo
term cancellation is not the issue.

Because the generic finite-Hamiltonian formulation closes analytically and the
Fe residual follows the same nonzero-`eta_int` trend, this audit does not
escalate to an LMTO representation mismatch.  If a future Fe ladder remains
nonzero after controlled \(\eta_{\rm int}\), window, and quadrature
extrapolations, the next audit is the already anticipated comparison between
the response GF and the physical Hamiltonian resolvent; that condition has
not been reached by the present evidence.

## 13. Literature mapping

### Katsnelson–Lichtenstein (2004)

Katsnelson and Lichtenstein start from the dynamical Kohn-Sham spin
susceptibility and distinguish the bare Kohn-Sham response from the later
exchange/kernel construction.  Their paper gives the spectral finite-state
structure of the transverse response and uses the orthogonal LMTO Hamiltonian
resolvent for the one-electron Green function.  This is the conceptual bridge
used here for the Lehmann denominator and the distinction between bare response
and later interaction physics.

The present audit uses only the bare \(\chi_0^{+-}\) part.  Their subsequent
frequency-dependent exchange and spin-wave discussion is outside DRESP-02F.

Primary source: [Katsnelson and Lichtenstein, JPCM 16, 7439 (2004), arXiv:cond-mat/0406488](https://arxiv.org/abs/cond-mat/0406488).

### Lounis–Costa–Muniz–Mills (2011)

Equation (2) of the Jülich/KKR-GF paper writes the Kohn-Sham transverse
response as a real-energy integral containing

\[
 -\frac1\pi\int dz\,f(z)\left[
 G^R_\downarrow(z+\omega)\operatorname{Im}G^R_\uparrow(z)
 +\operatorname{Im}G^R_\downarrow(z)G^A_\uparrow(z-\omega)
 \right].
\]

Using

\[
 -\frac1\pi\operatorname{Im}G^R
 =\frac{i}{2\pi}(G^R-G^A)=A
\]

maps it directly to the two source Kubo terms.  Their Eq. (26) further splits
the real-axis contour treatment into analytic and nonanalytic pieces; that
split is a contour-evaluation organization, not a different finite-spectrum
response.  Their later Eqs. (43) onward introduce a Dyson equation and an
effective (U), which are explicitly excluded here.  Their Eqs. (36)-(39)
also make clear that a Green function represented through energy-independent
wave functions and spectral poles must still be contracted with the correct
physical vertex; that representation question is separate from the generic
finite-Hamiltonian proof above.

Primary source: [Lounis, Costa, Muniz, and Mills, PRB 83, 035109 (2011), full PDF](https://juser.fz-juelich.de/record/14113/files/PhysRevB.83.035109.pdf).

The finite-temperature exactness statement in this document follows from the
finite-spectrum spectral representation and is not attributed to a
zero-temperature upper-limit form in the literature.

## 14. Root-cause verdict

The observed GF/Lehmann nonclosure at the production controls is explained by
three numerical facts, in descending order of importance:

1. `integration_eta` broadens the spectral delta in each one-(A)-one-(G)
   Kubo term, changing an isolated pole from width
   \(\eta_{\rm response}\) to
   \(\eta_{\rm response}+\eta_{\rm int}\);
2. the same Lorentzian broadening replaces exact endpoint occupations in
   integrated spectral weight by (F_{\eta_{\rm int}}(\epsilon));
3. the finite energy window and Simpson mesh omit Lorentzian tails and add
   quadrature error.

The following possible causes were audited and rejected:

* wrong retarded denominator sign;
* wrong (+i0^+) prescription;
* wrong up/down ordering;
* missing second Kubo term;
* shifted or duplicated Fermi factor;
* q versus −q endpoint order;
* transpose versus Hermitian transpose;
* real-only fixture cancellation;
* optimized-GF contraction mismatch;
* generic finite-Hamiltonian versus GF representation mismatch.

Therefore the exact requested classification is:

\[
\boxed{\texttt{PASS — CURRENT GF FORMULATION PROVEN EQUIVALENT}}
\]

The correct operational interpretation is that a finite-`integration_eta`
GF result is a controlled broadened approximation to the Lehmann response,
not an algebraically identical value at the same `eta_response`.  No
interacting TDDFT step is authorized by this result.
