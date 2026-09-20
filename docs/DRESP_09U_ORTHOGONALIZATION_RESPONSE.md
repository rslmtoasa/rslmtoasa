# DRESP-09U — production LMTO representation tangent

Status: **PASS-B for the Hamiltonian representation chain; density coefficient
map remains open**

Primary classification: **`FIELD_CLOSED_DENSITY_REPRESENTATION_OPEN`**

This milestone is restricted to a global `L=0` rigid spin rotation of the
accepted collinear, scalar-relativistic, orthogonal `ham_only` second-order
state on the 4x4x4 (64-point) Fe mesh. It does not add arbitrary `L>0` fields, ALSDA Ward/Dyson/spectra,
BES/Halle, or a Goldstone fit.

## Production map

The live code does not contain a generic Löwdin matrix `X=S**(-1/2)`. The
production representation is the following structured map:

```text
C, E_nu, SRDEL, QPAR
        | analytic predls tangent
        v
predls: center_band, shifted_band, width_band, obar
        |
        v
build_pot: cx / wx / cex / obx
        |
        | spin average/vector decomposition after the transform
        v
cx0/cx1, wx0/wx1, cex0/cex1, obx0/obx1
        |
        +--> chbar_nc --------------------> ee
        +--> build_obarm ------------------> obarm
        +--> build_enim -------------------> enim
        v
h(k), o(k), E_nu(k)
        |
        v
H2(k) = E_nu(k) + h(k) - h(k)o(k)h(k)
```

The source trace is `symbolic_atom.f90:predls` and
`symbolic_atom.f90:build_pot`, followed by `hamiltonian_build.f90` in
`build_bulkham`, `build_obarm`, and `build_enim`. `predls` uses
`WOW=wsm/ws_r`, `I=l+1`, and the canonical or explicitly supplied `QM(l)`.

For each spin channel:

\[
 d= C-E_\nu,\qquad A=QI-Q_M,\qquad
 \Delta=\mathrm{SRDEL}\,\mathrm{WOW}^{1/2-I},\qquad
 QI=\mathrm{QPAR}\,\mathrm{WOW}^{1-2I},
\]

\[
 X=1-\frac{Ad}{\Delta^2},\qquad
 Y=\frac{A}{dA-\Delta^2}.
\]

The analytic tangent is:

\[
 \delta d=\delta C-\delta E_\nu,\quad
 \delta\Delta=\mathrm{WOW}^{1/2-I}\delta\mathrm{SRDEL},\quad
 \delta QI=\mathrm{WOW}^{1-2I}\delta\mathrm{QPAR},
\]

\[
 \delta X=-\frac{\delta A\,d+A\,\delta d}{\Delta^2}
              +\frac{2Ad}{\Delta^3}\delta\Delta,
\]

\[
 N=dA-\Delta^2,\qquad
 \delta Y=\frac{\delta A}{N}-\frac{A(\delta d\,A+d\,\delta A-2\Delta\,\delta\Delta)}{N^2}.
\]

Hence:

\[
\begin{aligned}
 \delta C^{TB} &= \delta d\,X+d\,\delta X+\delta E_\nu,\\
 \delta C_{ex}^{TB} &= \delta d\,X+d\,\delta X,\\
 \delta W^{TB} &= \delta\Delta\,X+\Delta\,\delta X,\\
 \delta o^{TB} &= \delta Y.
\end{aligned}
\]

The implementation is the pure service
`source/lr_lmto_representation_tangent.f90`; it reproduces the live algebra
without modifying `predls`. Central differences appear only in
`tests/unit/test_lr_lmto_representation_tangent.f90`.

## Local and laboratory spin frames

`predls` and `build_pot` operate on local radial eigenchannels. Under the
rigid rotation used here, the eigenchannel parameters are invariant:

```text
delta C = delta E_nu = delta SRDEL = delta QPAR = 0
```

The laboratory spin representation changes through the moment direction in
the Pauli lift. For a positive `y` rotation,
`delta moment = [moment_z, 0, -moment_x]`. For every transformed scalar/vector
pair `a0,a1`, the production matrix is

\[
 a=a_0 I+a_1\,\mathbf m\cdot\boldsymbol\sigma,
\]

and, with channel derivatives retained for the general tangent,

\[
 \delta a=\delta a_0I+
 (\delta a_1\mathbf m+a_1\delta\mathbf m)\cdot\boldsymbol\sigma.
\]

The spin decomposition is therefore performed after `predls`, exactly as in
`build_pot`. The service `lmto_spin_parameter_tangent` applies the same
cartesian-to-spherical `hcpx` seam as `build_obarm` and `build_enim`.

The `enim` source semantics are not guessed: the live builder uses
`cx-cex` for each spin channel, which equals the live `center_band -
shifted_band` channel and therefore `E_nu + VMAD` in the rigid local-frame
case. `obarm` uses `obx` directly.

## Bond tangent

`chbar_nc` consumes the transformed `wx` and either `cex` (HOH) or `cx`
(first order). The DRESP-09U service differentiates the actual four-channel
bond polynomial, including both endpoint moments, the dot/cross terms, and
the onsite `c0/c1` term. It then applies the production `hhmag -> spinor`
packing and `hcpx` conversion. No native tangent is subtracted or used as a
construction input.

The material bridge compares this independently assembled `ee`, `obarm`, and
`enim` tangent with the DRESP-08 native path and then differentiates the live
second-order product:

\[
 \delta H^{(2)}=\delta E_\nu+\delta h
 -\delta h\,o\,h-h\,\delta o\,h-h\,o\,\delta h.
\]

## Residual ledger

| arrow | independent check | DRESP-09U result |
|---|---|---:|
| orthogonal parameters → `predls` | analytic tangent vs central difference | **closed**; `theta^2` convergence |
| transformed channels → spin matrix | finite rotated spin blocks | **closed** |
| transformed channels → `ee` | onsite, left endpoint, right endpoint, combined | **closed** |
| `obx` → `obarm` | live builder tangent | **closed** |
| `cx-cex` → `enim` | live builder tangent | **closed** |
| `ee` → `h(k)` | every accepted k point | **closed** |
| `h,o,E_nu` → `H2(k)` | independent product rule vs native tangent | **closed** |
| old DRESP-09T endpoint Gram diagnostic | `S K=Gamma` in radial endpoint space | not the production map |
| production density coefficients | explicit coefficient/operator `X` | **open** |

The old DRESP-09T metric residual `1.7646` is not repaired by adding a
second connection. It is an invalid diagnostic for this final structured
orthogonal/TB map: the live code exposes transformed potential parameters and
the `obar` operator, not an orthogonalizing coefficient matrix. Accordingly,
the DRESP-09U bridge does not publish a fictitious `X` or `delta X`.

The DRESP-09T density correction is retired from interpretation. A common
field/density coefficient map must be established before a corrected density
residual can be reported.

## Negative controls and scope

The unit fixture includes a nonmagnetic control with equal up/down channels;
all transverse transformed-parameter and bond tangents vanish. It also uses a
diagonal commuting structure-constant block as an independent closed-form
oracle. The accepted material bridge remains global `L=0`; arbitrary `L>0`
ALSDA work is deliberately not started.

## Verification

```text
cmake --build build -j2
ctest --test-dir build -R UnitLrLmtoRepresentationTangent --output-on-failure
```

The optional material gate uses
`tests/integration/tddft_driver_smoke/input_dresp09u_fe.nml` with backend
`representation_tangent` and is checked by
`tests/validation/dresp09u_fe_artifact.py`. DRESP-08, DRESP-09R, DRESP-09S,
and DRESP-09T remain independent paths.

Primary classification: **`FIELD_CLOSED_DENSITY_REPRESENTATION_OPEN`**
