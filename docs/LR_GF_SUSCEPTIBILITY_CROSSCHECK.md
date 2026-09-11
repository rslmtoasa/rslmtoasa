# LR-GF-02 — reciprocal-GF cross-check of the bare transverse susceptibility

## Disposition

**PASS — independent reciprocal-GF real-axis bubble agrees with the LR-06
spectral/Lehmann susceptibility on the certified finite response fixture.**

This is a representation-level numerical cross-check. It does not validate
the ALSDA kernel, a Dyson enhancement, Goldstone behavior, magnons, or a
converged Fe/Ni material result.

## Scope and provenance

The implementation is limited to the established baseline:

- reciprocal `ham_only`;
- orthogonal second-order/HOH eigensystems;
- collinear, no SOC, no additive response operator;
- complete finite band sets and explicit Fermi occupations;
- the LR-05 Pauli/no-SOC radial/angular vertex;
- the LR-04 canonical right-weighted response matrix.

The left state is at `k`; the endpoint state is the exact folded `k+q`
eigensystem. The GF route uses the left k-point weights and the same endpoint
ordering as LR-06. It does not call the LR-06 band-pair susceptibility
accumulator.

## Retarded/advanced derivation

For one circular sector, let (V_I(E_L,E_R)) be the coefficient-space
operator whose matrix element is the LR-05 transition vector,

\[
 T_{nm;I}=c_{n\mathbf k}^{\dagger}
 V_I(\epsilon_{n\mathbf k},\epsilon_{m,\mathbf k+\mathbf q})
 c_{m,\mathbf k+\mathbf q}.
\]

The LR-05 LMTO large-component reconstruction is affine in each endpoint
energy. Therefore the response vertex is represented exactly as

\[
 V_I(E_L,E_R)=\sum_{p,q=0}^{1}E_L^pE_R^q V_I^{pq}.
\]

The implementation evaluates the associated energy moments of the resolvent,
so this energy dependence is not frozen at a band energy or silently dropped.

For endpoint (a\in\{L,R\}), define the orthogonal reciprocal resolvents

\[
 G_a^R(E)=\sum_n\frac{|c_{na}\rangle\langle c_{na}|}
 {E-\epsilon_{na}+i\eta},\qquad
 G_a^A(E)=\sum_n\frac{|c_{na}\rangle\langle c_{na}|}
 {E-\epsilon_{na}-i\eta},
\]

and the spectral discontinuity

\[
 A_a(E)=\frac{i}{2\pi}\left[G_a^R(E)-G_a^A(E)\right].
\]

With the LR-03 measurement-first/source-second ordering, the retarded bubble
used by the implementation is

\[
\begin{aligned}
 \chi^R_{IJ}(\mathbf q,\omega)=
 \frac{2}{N_k}\sum_{\mathbf k}\int dE\,f(E)\,\big[&
 \operatorname{Tr}\{A_L(E)V_I G_R^R(E+\omega)V_J^\dagger\}\\
 &+\operatorname{Tr}\{A_R(E)V_J^\dagger G_L^A(E-\omega)V_I\}\big].
\end{aligned}
\]

The placements are therefore:

| object | placement |
|---|---|
| Fermi function | (f(E)), multiplying both Kubo terms |
| spectral/imaginary part | (A=i(G^R-G^A)/(2\pi)) at the integration energy |
| retarded GF | right endpoint, (G_R^R(E+\omega)), with (+i\eta) |
| advanced GF | left endpoint, (G_L^A(E-\omega)), with (-i\eta) |
| endpoint ordering | left (k), right exact folded (k+q) |
| circular vertex | (sigma^+) for `chi_plus`, (sigma^-) for `chi_minus` |
| LR-04 metric | applied once after the raw GF bubble as (B=\chi_{raw}W) |

Expanding the first term on the left pole gives

\[
 \frac{f_{n\mathbf k}}
 {\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta}.
\]

For the second term, the right pole is at (E=\epsilon_m), and

\[
 \frac{1}
 {\epsilon_m-\omega-\epsilon_n-i\eta}
 =-\frac{1}
 {\omega+\epsilon_n-\epsilon_m+i\eta},
\]

which supplies the (-f_m) contribution. The sum is consequently

\[
 \chi^R_{IJ}=\frac{2}{N_k}\sum_{\mathbf k,nm}
 \frac{f_{n\mathbf k}-f_{m,\mathbf k+\mathbf q}}
 {\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta}
 T_{nm;I}T_{nm;J}^{*},
\]

which is exactly the LR-03/LR-06 convention, including the circular factor
of two and the retarded (+i\eta) sign.

## Numerical integration

The first implementation uses a fixed real-axis composite Simpson rule. The
energy interval is the complete finite spectrum of both endpoints plus an
explicit `energy_margin`; the number of points is the odd request field
`integration_points`.

The spectral discontinuity has its own explicit `integration_eta`, separate
from the physical response broadening `eta`. This makes the two numerical
errors distinguishable:

- reducing `integration_eta` controls the finite-width approximation to the
  spectral discontinuity;
- increasing `integration_points` controls Simpson discretization;
- `eta` controls the requested retarded response itself.

The default is correctness-oriented rather than production-performance
oriented. No DOS/GF mesh, contour convention, or hidden broadening is reused.

## Independent implementation boundary

The backend builds (G^{(p)}(z)=\sum_n\epsilon_n^p|n\rangle\langle n|/(z-
\epsilon_n)) for (p=0,1,2), contracts the four affine vertex components,
and performs the energy integral. It shares only the certified eigensystem,
response-space, radial augmentation, angular Gaunt, and Pauli conventions.
There is no call to `evaluate_lr_ks_susceptibility` or to its nested band-pair
accumulator.

## Tests and evidence

Focused command:

```text
cmake --build build --target UnitLrGfSusceptibility -j2
ctest --test-dir build --output-on-failure -R '^UnitLrGfSusceptibility$'
```

Observed output from the certified fixture:

```text
GF coarse full-response difference (abs, rel) =   2.1443E-09   2.4633E-02
GF fine full-response difference (abs, rel) =     3.6257E-09   4.1652E-02
Finite-spectrum analytic GF error =                1.1277E-09
Static q=0 full-response GF difference =           9.2484E-10
GF off-mesh-q full-response difference (abs, rel) = 1.7422E-08 2.0014E-01
GF alternate eta full-response difference (abs, rel) = 1.7456E-08 2.0017E-01
UnitLrGfSusceptibility: PASS (GF bubble, convergence, full response, static limit)
```

The comparison norm is over every response-space matrix element and frequency
in the one-site `sp` direct-radial fixture; it is not a site-projected scalar.
The off-mesh case uses (q_x=0.23) and an exact folded endpoint.
The relative norm is reported for scale context, while the acceptance oracle
uses the absolute full-response norm because the direct radial metric contains
small-volume entries.

The finite-spectrum analytic oracle is the occupied-up to empty-down two-level
transition with its LR-04 radial metric. The same test varies the energy
integration resolution and `integration_eta`, compares the complete reciprocal
fixture against LR-06, and checks the independent `q=0, omega=0` limit.

LR-03 supplies retarded/advanced covariance and the static response identity,
but no separate finite-cutoff frequency-integral sum rule for this circular
matrix. No unprovided sum rule is enforced here.

## Error budget and claim boundary

| source | shared or independent? | disposition |
|---|---|---|
| energy integration / spectral-discontinuity width | independent GF route | varied in the focused test and exposed in the request |
| response `eta` | common requested physics parameter | held equal between GF and LR-06 |
| k mesh and finite electronic basis | shared | not independently validated by this cross-check |
| SR→Pauli approximation | shared LR-05 vertex | inherited, not revalidated here |
| radial/angular response basis | shared certified mapping | full response matrix is compared |

The PASS claim is therefore an independent reciprocal-GF representation check
of the bare KS susceptibility. It is not ALSDA validation, Dyson validation,
Goldstone validation, converged-material validation, or literature agreement.

## Completion checklist

- [x] LR-03 Kubo ordering and circular convention derived before coding;
- [x] retarded and advanced GF placements documented;
- [x] real-axis integration strategy documented;
- [x] GF backend independent of the LR-06 spectral accumulator;
- [x] finite-spectrum analytic oracle passes;
- [x] integration resolution/broadening varied;
- [x] full response-space reciprocal fixture compared;
- [x] off-mesh q with exact folded endpoint compared;
- [x] q=0, omega=0 static limit compared;
- [x] error budget separates independent and shared errors;
- [x] no Kxc, Dyson, Goldstone, or site-only oracle added.

## Commit

`linear-response: cross-check KS susceptibility with reciprocal GF`
