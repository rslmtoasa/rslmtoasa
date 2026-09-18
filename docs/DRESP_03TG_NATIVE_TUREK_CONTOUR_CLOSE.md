# DRESP-03TG-CLOSE — Native Turek contour and (J_{ij}/J(q)) closure

Status: **PASS — R1 native algebra, material convergence, and unit gates
pass.**

The fixed-complex-energy representation audit is extended here to a native
finite-temperature contour evaluator.  Every native resolvent is formed by a
direct complex solve of

\[
g^\alpha(k,z)=[P^\alpha(z)-S^\alpha(k)]^{-1};
\]

no finite-(H) resolvent, eigenbasis, fitted scale, moment divisor, or
physical/true-Green-function substitution is used.

The historical representation audit remains in
[`DRESP_03TG_NATIVE_TUREK_GF.md`](DRESP_03TG_NATIVE_TUREK_GF.md).  Its frozen
`source/exchange.f90` invariant is unchanged.

## Native contour

`source/lr_lmto_turek_contour.f90` provides:

- an explicitly parameterized counter-clockwise ellipse with
  `dz/(2*pi*i)` weights;
- pole-subtracted `f_reg` weighting with explicit enclosed Matsubara-pole
  residues;
- direct native (P-S) solves at every contour node, k point, and q point;
- structure-constant caching across contour nodes; and
- native absolute (J_{ij}(q)) and
  \(\Delta J(q)=J(\Gamma)-J(q)\) production values.

The certified finite-temperature prescription is

\[
 I=\sum_n w_n f_{reg}(z_n)R(z_n)+k_BT\sum_{p\in C}R(p),
 \qquad
 f_{reg}(z)=f(z)+k_BT\sum_{p\in C}\frac{1}{z-p},
\]

where the contour is counter-clockwise and each weight contains
`dz/(2*pi*i)`.  The scalar oracle evaluates this identity at several real
energies and verifies the Fermi-pole residue sign and normalization.  The
conversion from the occupied contour to the real-axis LKAG convention is
implemented as

\[
J_{ij}(q)=-\frac14\operatorname{Re}
 \oint_C\frac{dz}{2\pi i}\,f(z)
 \operatorname{Tr}[\Delta P_i g^\uparrow_{ij}(k,z)
 \Delta P_j g^\downarrow_{ji}(k+q,z)].
\]

The real-space closure is explicit and normalized:

\[
J_{ij}(R)=\frac1{N_q}\sum_q J_{ij}(q)e^{+i2\pi q\cdot R},
\qquad
J_{ij}(q)=\sum_R J_{ij}(R)e^{-i2\pi q\cdot R}.
\]

The former `native_exchange_jij_contour` entry point remains an explicitly
labelled DFT helper.  It is not the independent closure.  The R1 route first
solves the native path operator on the complete electronic mesh and forms

\[
g^\uparrow_{ij}(R,z)=\sum_k w_k e^{-i2\pi k\cdot R}g^\uparrow_{ij}(k,z),
\quad
g^\downarrow_{ji}(-R,z)=\sum_k w_k e^{+i2\pi k\cdot R}g^\downarrow_{ji}(k,z),
\]

then contour-integrates the ordered `ud` and `du` products directly.  Direct
reciprocal `J(q)` is compared with an independently evaluated real-space
`J(R)` on complete compatible meshes; no fitted scale or moment divisor is
used.

The `exchange_q` namelist accepts `native_turek=.true.`.  The older
`native_crosscheck` switch remains a compatibility alias for enabling the same
native contour route.  `native_contour_points`, `native_contour_margin`,
`native_contour_height_fraction`, and
`native_contour_account_fermi_poles` and the optional even
`native_contour_target_fermi_poles` control the native contour independently
of the finite-(H) contour.

## Hard-gate evidence

`UnitDresp03tgNativeContour` uses the live bcc-Fe spd fixture and checks:

1. contour construction and vanishing total contour weight;
2. the scalar pole oracle, including raw-Fermi and `f_reg` forms;
3. finite native `J(q)` at Gamma, +q, and -q;
4. the algebraic DFT helper identity;
5. independent native ordered-pair `J(R)` versus direct `J(q)` on a complete
   2³ mesh;
6. `ud(q)=du(-q)`, one-site inversion, and compatible-mesh `ud=du`;
7. native spectral bounds from the exact native coefficient problem; and
8. the analytic complex-amplitude Heisenberg factor
   `K(q)=2*(J_sym(Gamma)-J_sym(q))`.

Representative output:

```text
contour weight residual       4.654752e-17
contour nodes/poles           1024/512     (complete 2^3 live mesh)
native J(q)                  fixture-dependent; see test log
synthetic Fourier residual    1.110223e-16
independent pair/q residual   <8e-16
native max spectral ellipse  4.942282e-1
```

The targeted regression set also passes `UnitDresp03tgNativeFixedZ`,
`UnitDresp03qProductionAdapter`, and `UnitExchangeQ`.

## bcc-Fe production validation

`example/exchange_q/bccFe/input_dresp03tg_close_12.nml` enables the native
route on the accepted 12x12x12 bcc-Fe k mesh, a four-point commensurate q path,
and a fixed 256-pole target.  The production output writes ordered `J_ud`,
`J_du`, visible `J_sym`, raw `DeltaJ_sym`, `DeltaJ_sym/q^2`, the expected
curvature `2*DeltaJ_sym`, and `finiteH_minus_native_curvature`.  The finite-H
columns are retained as diagnostics; they are not substituted for the native
path-operator result.  Native spectral bounds come from the exact native
coefficient problem and the report records the maximum ellipse value.

The previous same-state values are retained only as a column-format example;
they are not the R1 convergence record:

| q | native (J(q)) (Ry) | native \\(\Delta J\\) (Ry) | finite-H \\(\Delta J\\) (Ry) | finite-H minus native (Ry) |
|---:|---:|---:|---:|---:|
| 0 | 1.938217541e-1 | 0 | -4.70e-19 | -4.70e-19 |
| 1/12 | 1.926296575e-1 | 1.192096588e-3 | 1.620624491e-3 | 4.285279025e-4 |
| 1/6 | 1.891041813e-1 | 4.717572754e-3 | 1.073458918e-3 | -3.644113837e-3 |
| 1/4 | 1.876184770e-1 | 6.203277102e-3 | 7.760220650e-3 | 1.556943548e-3 |

No fitted scale or q-dependent correction is applied.  The final R1 table
below reports the 64/96/128 contour comparison and the 8³/12³/16³ electronic
mesh comparison at common small-q points.  Native and finite-H remain
separate representations and observables.

### R1 convergence record

The production deck uses target 256 Fermi poles.  All reported native
spectral poles are strictly inside the ellipse.  The contour sweep is at a
fixed 12³ electronic mesh; the k-mesh sweep uses 64 contour points.  The
finite-H comparison is retained as a diagnostic and is not used to define the
native result.

| contour points | native max ellipse | max `DeltaJ` (Ry) | `DeltaJ/q^2` at q=1/12, 1/6, 1/4 (Ry A²) | max `|finiteH-native curvature|` (Ry) |
|---:|---:|---:|---:|---:|
| 64  | 0.937909 | 5.941679182e-3 | 1.605920433e-2, 1.695239825e-2, 9.856798571e-3 | 8.010014321e-3 |
| 96  | 0.937909 | 5.943569665e-3 | 1.606316202e-2, 1.695600614e-2, 9.859934740e-3 | 8.011947506e-3 |
| 128 | 0.937909 | 5.940479554e-3 | 1.605711301e-2, 1.695031017e-2, 9.854808479e-3 | 8.008895478e-3 |

The maximum change in native `DeltaJ` over the contour sweep is
`3.09e-6 Ry`; the small-q stiffness diagnostics vary by less than
`3.7e-6 Ry A²`.  The absolute `Jsym` values have a nearly q-independent
quadrature offset, so convergence is assessed on `DeltaJ` and the reported
curvature rather than on the absolute contact trace.

| electronic k mesh | native max ellipse | `DeltaJ` at q=1/12, 1/6, 1/4 (Ry) | max `|finiteH-native curvature|` (Ry) |
|---:|---:|---:|---:|
| 8³  | 0.955924 | 4.322092089e-4, 3.255325653e-3, 8.178703727e-3 | 1.846960583e-2 |
| 12³ | 0.937909 | 1.075610000e-3, 4.541736619e-3, 5.941679182e-3 | 8.010014321e-3 |
| 16³ | 0.988038 | 1.830421048e-3, 6.450457874e-3, 7.554837997e-3 | 7.240543352e-3 |

The k-mesh rows are accepted-state material validation runs; their Fermi
levels are independently converged for each mesh.  They therefore test
electronic-mesh stability, while the contour table tests the integration
resolution at a common 12³ state.

## Verification

```text
cmake --build build -j2
ctest --test-dir build --output-on-failure -R \
  'UnitDresp03tgNativeContour|UnitDresp03tgNativeFixedZ|UnitDresp03qProductionAdapter|UnitExchangeQ'
```

The independent unit gates and the production convergence runs pass on branch
`fable_v4`.  Generated production outputs remain outside the repository.
