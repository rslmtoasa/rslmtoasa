# DRESP-03TG-CLOSE — Native Turek contour and (J_{ij}/J(q)) closure

Status: **CLOSED for the native path-operator and Fourier closure scope.**

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
- stable complex Fermi weighting with explicit enclosed Matsubara-pole
  residues;
- direct native (P-S) solves at every contour node, k point, and q point;
- structure-constant caching across contour nodes; and
- native absolute (J_{ij}(q)) and
  \(\Delta J(q)=J(\Gamma)-J(q)\) production values.

The conversion from the occupied contour to the real-axis LKAG convention is
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

The `exchange_q` namelist accepts `native_turek=.true.`.  The older
`native_crosscheck` switch remains a compatibility alias for enabling the same
native contour route.  `native_contour_points`, `native_contour_margin`,
`native_contour_height_fraction`, and
`native_contour_account_fermi_poles` control the native contour independently
of the finite-(H) contour.

## Hard-gate evidence

`UnitDresp03tgNativeContour` uses the live bcc-Fe spd fixture and checks:

1. contour construction and vanishing total contour weight;
2. finite native (J(q)) at (Gamma), (+q), and (-q);
3. the algebraic (J(q)\leftrightarrow J_{ij}) Fourier identity; and
4. the same (J(q)\to J_{ij}\to J(q)) closure applied to a complete live
   (2^3) native q mesh.

Representative output:

```text
contour weight residual       4.654752e-17
contour nodes/poles           256/1520     (complete 2^3 live mesh)
native J(q)                  -2.245755e-1  -1.820462e-1  -1.820462e-1 Ry
synthetic Fourier residual    1.110223e-16
native Jij/J(q) residual     5.551115e-17
```

The targeted regression set also passes `UnitDresp03tgNativeFixedZ`,
`UnitDresp03qProductionAdapter`, and `UnitExchangeQ`.

## bcc-Fe production validation

`example/exchange_q/bccFe/input_dresp03tg_close_12.nml` enables the native
route on the accepted 12x12x12 bcc-Fe k mesh and a four-point commensurate q
path.  The production output is
`exchange_q_dresp03tg_close_12.dat`; it writes native absolute (J(q)), native
\(\Delta J(q)\), the finite-(H) value, and their direct residual in Ry.

Representative same-state values from the deck are:

| q | native (J(q)) (Ry) | native \\(\Delta J\\) (Ry) | finite-H \\(\Delta J\\) (Ry) | finite-H minus native (Ry) |
|---:|---:|---:|---:|---:|
| 0 | 1.938217541e-1 | 0 | -4.70e-19 | -4.70e-19 |
| 1/12 | 1.926296575e-1 | 1.192096588e-3 | 1.620624491e-3 | 4.285279025e-4 |
| 1/6 | 1.891041813e-1 | 4.717572754e-3 | 1.073458918e-3 | -3.644113837e-3 |
| 1/4 | 1.876184770e-1 | 6.203277102e-3 | 7.760220650e-3 | 1.556943548e-3 |

No fitted scale or q-dependent correction is applied.  The finite-H spectral
and finite-H contour columns are also both written by the deck; their largest
reported same-state residual over these four points is (8.84\times10^{-6})
Ry.  The native values are therefore a material validation record and a
direct comparison surface, not a hidden replacement of the certified
finite-H observable.

## Verification

```text
cmake --build build -j2
ctest --test-dir build --output-on-failure -R \
  'UnitDresp03tgNativeContour|UnitDresp03tgNativeFixedZ|UnitDresp03qProductionAdapter|UnitExchangeQ'
```

All four targeted tests pass on branch `fable_v4`.
