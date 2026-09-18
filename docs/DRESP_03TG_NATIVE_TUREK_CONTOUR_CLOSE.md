# DRESP-03TG-FINAL — Native Turek contour and exchange closure

Status: **CORE PASS.** The native pole algebra, direct (P-S) pole oracle,
scalar contour oracle, independent ordered-pair Fourier closure, and fixed-z
regressions pass. The bounded Fe study reports **contour convergence PASS**,
**Fe electronic k-mesh convergence NOT CONVERGED THROUGH (24^3)**, and
**small-q Fe window NOT RESOLVED AT (24^3)**. The latter two are production
material limitations, not native-representation blockers.

The native production object is always the direct path operator

\[
g^\alpha(k,z)=[P^\alpha(z)-S^\alpha(k)]^{-1},
\]

with no finite-(H) replacement, fitted scale, moment divisor, or
true-Green-function substitution. `source/exchange.f90` remains frozen.

## Native pole repair

The pole constructor now exposes the coefficient ordering explicitly:

\[
Q=QI-\alpha,\qquad M=I-SQ,
\]
\[
A=MD,\qquad B=MDC+S,
\]

and solves (Az=Bz) through the existing LAPACK path. Thus the generalized
problem is

\[
[(I-SQ)D]z=(I-SQ)DC+S,
\]

with diagonal factors applied by columns in the stated order. The direct
production (P^\alpha(z)-S^\alpha(k)) SVD oracle is independent of the pole
construction and reports:

```text
old generalized-root max relative P-S residual = 3.140490e-01
corrected generalized-root max relative P-S residual = 3.987466e-15
direct P-S sigma_min max / relative max       = 6.471473e-14 / 3.987466e-15
representative sigma_min / relative / ||P-S|| = 8.881784e-16 / 2.404530e-17 / 3.693772e+01
```

The corrected roots therefore pass the required (10^{-9}) relative pole
residual gate, while the negative control demonstrates that the old
((I-S)D) construction would have been caught.

## Ordered-pair and non-self-inverse-q closure

`native_exchange_pairs_contour` constructs all four Fourier blocks
independently:

```text
gup_R, gup_minus_R, gdown_R, gdown_minus_R
```

and uses

\[
J_{ud}(R)\sim\Delta P_i g^\uparrow_{ij}(R)\Delta P_j
g^\downarrow_{ji}(-R),
\]
\[
J_{du}(R)\sim\Delta P_i g^\downarrow_{ij}(R)\Delta P_j
g^\uparrow_{ji}(-R).
\]

The hard regression uses the cyclic mesh (k=0,1/3,2/3),
(q=0,+1/3,-1/3), and (R=0,1,2), so (+q\ne-q\pmod G). The old
2³ Fourier utility check remains, but is labelled only as a DFT helper.

```text
DFT helper roundtrip                         = 3.469447e-18
independent ud pair/q residual              = 1.695382e-17
independent du pair/q residual              = 2.036634e-17
ordered maximum                             = 2.036634e-17
Jud(q)-Jdu(-q) covariance residual           = 2.480342e-17
Jud(q)-Jdu(q) one-site symmetry residual     = 3.760107e-17
```

These are below the (10^{-10}) closure target and do not rely on the
self-inverse (q=1/2) case.

## Scalar contour and finite-q oracles

The raw-Fermi and pole-subtracted scalar contour forms remain independently
checked, including explicit enclosed Fermi-pole residues:

```text
scalar contour oracle coarse / fine         = 9.785529e-09 / 4.551914e-15
scalar LKAG -1/4 oracle coarse / fine       = 2.446382e-09 / 9.436896e-16
contour weight residual                     = 1.321838e-17
contour nodes / poles                       = 192 / 96
curvature-factor oracle                     = 8.326673e-17
```

The independently derived complex-amplitude normalization remains

\[
\boxed{K(q)=2[J_{\rm sym}(0)-J_{\rm sym}(q)]}.
\]

Production output keeps `J_ud`, `J_du`, `J_sym`, `DeltaJ_sym`,
`DeltaJ_sym/q^2`, `native_expected_curvature=2*DeltaJ_sym`, and the
finite-(H) diagnostic separate.

The unit fixture has native maximum ellipse value
`8.227106e-1`, giving safety margin `1.772894e-1`. All production poles in
the following tables are likewise strictly inside the certified ellipse.

## Bounded bcc-Fe campaign

The contour sweep uses a fixed accepted (12^3) electronic state and 64,
96, and 128 contour nodes. Values below use

\[
\Delta J=J_{\rm sym}(0)-J_{\rm sym}(q).
\]

| contour | Fermi level (Ry) | max ellipse | `J_sym(0)` (Ry) | `DeltaJ` at (1/12,1/6,1/4) (Ry) |
|---:|---:|---:|---:|---:|
| 64  | -6.7656145763e-2 | 0.822621179 | -2.307630796e-1 | 1.075075406e-3, 4.539708756e-3, 5.937549878e-3 |
| 96  | -6.7656145914e-2 | 0.822621178 | -2.284329747e-1 | 1.074019238e-3, 4.535753383e-3, 5.929604706e-3 |
| 128 | -6.7656146771e-2 | 0.822621178 | -2.324144088e-1 | 1.074721047e-3, 4.538388585e-3, 5.934918833e-3 |

The maximum 96-to-128 change in `DeltaJ` is approximately
`5.314e-6 Ry`, so:

```text
CONTOUR CONVERGENCE = PASS
```

The fixed 64-node electronic-mesh campaign is:

| k mesh | Fermi level (Ry) | native poles | max ellipse | margin | `J_sym(0)` (Ry) |
|---:|---:|---:|---:|---:|---:|
| 12³ | -6.7656145763e-2 | 124416 | 0.822621179 | 0.177378821 | -2.307630796e-1 |
| 16³ | -6.8155280645e-2 | 294912 | 0.822631567 | 0.177368433 | -2.4454212877e-1 |
| 20³ | -6.9564231299e-2 | 576000 | 0.822624757 | 0.177375243 | -2.4017927379e-1 |
| 24³ | -6.9063005866e-2 | 995328 | 0.822637221 | 0.177362779 | -2.4112868568e-1 |

| k mesh | `DeltaJ` at (1/12,1/6,1/4) (Ry) | `2*DeltaJ` at (1/12,1/6,1/4) (Ry) | `DeltaJ/q²` at (1/12,1/6,1/4) (Ry Å²) |
|---:|---:|---:|---:|
| 12³ | 1.075075406e-3, 4.539708756e-3, 5.937549878e-3 | 2.150150812e-3, 9.079417513e-3, 1.187509975e-2 | 1.605122266e-2, 1.694482909e-2, 9.849948366e-3 |
| 16³ | 9.220898669e-4, 3.247599257e-3, 7.592281908e-3 | 1.844179734e-3, 6.495198514e-3, 1.518456382e-2 | 1.376709920e-2, 1.212192617e-2, 1.259502426e-2 |
| 20³ | 7.335189719e-4, 3.528718873e-3, 6.611110092e-3 | 1.467037944e-3, 7.057437746e-3, 1.322222018e-2 | 1.095167490e-2, 1.317122781e-2, 1.096733407e-2 |
| 24³ | 9.246643758e-4, 3.514857729e-3, 7.000061928e-3 | 1.849328752e-3, 7.029715458e-3, 1.400012386e-2 | 1.380553745e-2, 1.311948997e-2, 1.161257589e-2 |

For the two finest meshes, the relative (20^3\to24^3) change is about
26% at (q=1/12) and 0.39% at (q=1/6). The first required q point fails
the bounded criterion, so the campaign terminates without inventing another
implementation task:

```text
MATERIAL Fe k-MESH CONVERGENCE = NOT CONVERGED THROUGH 24^3
```

On the finest completed (24^3) mesh, the small-q diagnostic is:

| \(\xi\) | `DeltaJ` (Ry) | `DeltaJ/|q|²` (Ry Å²) |
|---:|---:|---:|
| 1/48 | 2.362088203e-5 | 5.642678232e-3 |
| 2/48 | 2.453698281e-4 | 1.465380278e-2 |
| 3/48 | 4.798173970e-4 | 1.273569518e-2 |
| 4/48 | 9.246643768e-4 | 1.380553747e-2 |

The final three coefficients span about 14%, above the approximate 10%
diagnostic window. No stiffness is declared:

```text
SMALL-q MATERIAL LIMIT = NOT RESOLVED AT AVAILABLE 24^3 MESH
```

The finite-(H) `finiteH_minus_native_curvature` column is a production-H1
diagnostic only. It is not an acceptance criterion for the exact native
Turek representation. A numerically equivalent historical shell comparison
would require matching the old exchange route's state, energy integration,
and normalization; it was not practical within this bounded campaign:

```text
HISTORICAL SHELL CROSSCHECK = DEFERRED NUMERICALLY
```

## Final decision

```text
CORE DRESP-03TG             = PASS
CONTOUR CONVERGENCE         = PASS
Fe ELECTRONIC k-MESH        = NOT CONVERGED THROUGH 24^3
SMALL-q Fe WINDOW           = NOT RESOLVED AT AVAILABLE MESH
HISTORICAL exchange.f90     = DEFERRED NUMERICALLY
```

The DRESP-03TG implementation and physics closure is complete. Remaining Fe
mesh refinement is a production convergence study, not an implementation
blocker. DRESP-04 remains outside this task.

## Verification

```text
cmake --build build -j2
ctest --test-dir build --output-on-failure -R \
  'UnitDresp03tgNativeContour|UnitDresp03tgNativeFixedZ|UnitDresp03qProductionAdapter|UnitDresp03|UnitExchangeQ|UnitLrKlStaticBridge'
```

The protected historical source has SHA-256:

```text
6e9ae7da6af46367aaf18a4fef37bdd4a9a0d6f69fe8e08b589e5c5f278a70e6
```
